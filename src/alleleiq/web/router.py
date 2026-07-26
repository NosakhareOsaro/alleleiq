"""Routes for the minimal server-rendered web interface.

Every route here delegates the actual analysis to
:func:`alleleiq.pipeline.build_evidence_report` -- the same function
:mod:`alleleiq.cli` and :mod:`alleleiq.api.routers.variants` call -- so
this interface is a presentation layer, not a second implementation of
variant analysis.
"""

from __future__ import annotations

import asyncio
import tempfile
import uuid
from collections.abc import AsyncIterator
from pathlib import Path

from fastapi import APIRouter, Depends, File, Form, HTTPException, Request, UploadFile, status
from fastapi.responses import FileResponse, HTMLResponse, StreamingResponse
from fastapi.templating import Jinja2Templates

from alleleiq.api.dependencies import get_settings_dependency
from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.clients.vep import VEPClient
from alleleiq.config import Settings
from alleleiq.exceptions import AlleleIQError, ConfigurationError, VariantParseError
from alleleiq.logging import get_logger
from alleleiq.models.report import EvidenceReport
from alleleiq.models.variant import GenomeBuild
from alleleiq.parsers.vcf_parser import parse_vcf
from alleleiq.pipeline import ProgressCallback, build_evidence_report
from alleleiq.reports import ReportOrchestrator

logger = get_logger(__name__)
router = APIRouter(tags=["web"])

_TEMPLATES_DIR = Path(__file__).parent / "templates"
templates = Jinja2Templates(directory=str(_TEMPLATES_DIR))

# Report filenames must never contain a path separator, and are always
# exactly `{variant_id}.json` / `{variant_id}.html` (see
# ReportOrchestrator.generate) -- '..' is rejected explicitly as
# defence in depth even though a filename without '/' or '\\' cannot
# otherwise escape its parent directory via Path resolution.
_ALLOWED_REPORT_SUFFIXES = (".json", ".html")

# In-memory, per-process progress queues keyed by job_id, one per
# in-flight analysis. Like RateLimitMiddleware's sliding-window counter
# (see alleleiq.api.middleware), this is a known, documented
# single-process scope limitation, not an oversight: a multi-worker
# deployment would need a shared pub/sub layer (e.g. Redis) for the SSE
# progress feed to reach a client whose GET .../progress request lands
# on a different worker than the POST /analyse that is driving it.
_PROGRESS_SENTINEL = "__complete__"
_progress_queues: dict[str, asyncio.Queue[str]] = {}


def _get_or_create_progress_queue(job_id: str) -> asyncio.Queue[str]:
    """Return job_id's progress queue, creating it if this is the first touch.

    Whichever of the SSE GET or the analysing POST reaches a given
    job_id first creates the queue; the other one just finds it
    already there, so the two requests can arrive in either order.
    """
    if job_id not in _progress_queues:
        _progress_queues[job_id] = asyncio.Queue()
    return _progress_queues[job_id]


def _make_progress_notifier(
    queue: asyncio.Queue[str], variant_index: int, total_variants: int
) -> ProgressCallback:
    """Build an :data:`~alleleiq.pipeline.ProgressCallback` for one variant in a batch.

    A plain closure over a loop variable (``queue``) would be rebound
    on every iteration of the caller's loop by the time the callback
    actually runs -- a real bug class, not just a lint nitpick -- so
    ``queue``, ``variant_index``, and ``total_variants`` are captured
    here as this factory's own arguments instead, each iteration
    getting its own independent closure.
    """

    def _notify_progress(stage: str) -> None:
        queue.put_nowait(f"variant {variant_index}/{total_variants}: {stage}")

    return _notify_progress


def _parse_optional_job_id(raw_job_id: str) -> str | None:
    """Validate a client-supplied job_id form field, if one was given.

    Returns None (silently -- an absent or malformed job_id just means
    "no progress tracking for this request", not an error) unless it's
    present and a syntactically valid UUID.
    """
    if not raw_job_id:
        return None
    try:
        uuid.UUID(raw_job_id)
    except ValueError:
        return None
    return raw_job_id


def _clinical_error_message(exc: Exception) -> str:
    """Translate an exception into a clinical-grade message, never a raw traceback.

    Args:
        exc: The exception raised while handling an upload.

    Returns:
        A clear, actionable, non-alarming message safe to show a
        reviewing scientist -- specific enough to act on for the
        anticipated :class:`~alleleiq.exceptions.AlleleIQError`
        subclasses, and a generic fallback for anything else (so an
        unanticipated internal error still never reflects raw
        exception text or a stack trace back to the browser).
    """
    if isinstance(exc, VariantParseError):
        return (
            f"This VCF could not be parsed: {exc}. Please check the file and try again. "
            "AlleleIQ requires well-formed VCF v4.x and supports GRCh38 only."
        )
    if isinstance(exc, ConfigurationError):
        return f"The server is misconfigured and cannot complete this analysis: {exc}"
    if isinstance(exc, AlleleIQError):
        return f"This analysis could not be completed: {exc}"
    logger.error("unexpected error handling web upload", extra={"error": str(exc)})
    return (
        "An unexpected error occurred while analysing this file. "
        "Please try again; if the problem persists, contact your system administrator."
    )


async def _connectivity_status(settings: Settings) -> dict[str, bool]:
    """Check gnomAD/VEP/ClinVar reachability, the same checks `/api/v1/health/ready` uses."""
    async with (
        GnomADClient(settings) as gnomad_client,
        VEPClient(settings) as vep_client,
        ClinVarClient(settings) as clinvar_client,
    ):
        return {
            "gnomad_ok": await gnomad_client.health_check(),
            "vep_ok": await vep_client.health_check(),
            "clinvar_ok": await clinvar_client.health_check(),
        }


@router.get("/", response_class=HTMLResponse, summary="Upload a VCF for analysis")
async def index(
    request: Request, settings: Settings = Depends(get_settings_dependency)
) -> HTMLResponse:
    """Render the upload form with live evidence-source connectivity status."""
    status_flags = await _connectivity_status(settings)
    return templates.TemplateResponse(
        request, "index.html", {**status_flags}, status_code=status.HTTP_200_OK
    )


async def _save_upload_to_tempfile(upload: UploadFile, max_bytes: int) -> Path:
    """Persist an uploaded VCF to a temporary file, enforcing a size limit.

    Mirrors :func:`alleleiq.api.routers.variants._save_upload_to_tempfile`,
    kept separate rather than shared because the two callers need to
    turn an oversized upload into different response shapes (a JSON
    413 for the REST API, a rendered ``error.html`` here).
    """
    contents = await upload.read()
    if len(contents) > max_bytes:
        raise ValueError(f"upload exceeds the {max_bytes}-byte limit")

    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as handle:
        handle.write(contents)
        tmp_path = Path(handle.name)
    return tmp_path


def _report_urls(job_id: str, report: EvidenceReport) -> dict[str, str]:
    base_name = report.variant.variant_id
    return {
        "json_url": f"/reports/{job_id}/{base_name}.json",
        "html_url": f"/reports/{job_id}/{base_name}.html",
    }


@router.post("/analyse", response_class=HTMLResponse, summary="Analyse an uploaded VCF")
async def web_analyse(
    request: Request,
    vcf_file: UploadFile = File(..., description="A VCF file containing one or more variants."),
    genome_build: str = Form("GRCh38"),
    include_clinvar: bool = Form(True),
    job_id: str = Form(
        "",
        description=(
            "Optional client-generated UUID (see index.html) used to correlate this "
            "request with a GET /analyse/{job_id}/progress SSE subscription. Analysis "
            "proceeds identically whether or not one is supplied."
        ),
    ),
    settings: Settings = Depends(get_settings_dependency),
) -> HTMLResponse:
    """Parse and analyse an uploaded VCF, returning a rendered result page or fragment.

    Renders a full page (extending ``base.html``) for a plain form
    submission, or just the result fragment when the request carries
    HTMX's ``HX-Request`` header -- both are backed by the same
    ``_result_content.html`` partial, so the two response shapes can
    never show different results for the same analysis.

    If a valid ``job_id`` was supplied, real progress events -- sourced
    from :mod:`alleleiq.pipeline`'s ``on_progress`` hook, not a
    simulated/fake progress bar -- are pushed to that job's queue as
    each variant moves through evidence-gathering and scoring, for any
    concurrently-open ``GET /analyse/{job_id}/progress`` SSE connection
    to pick up. The queue's completion sentinel is always sent from a
    ``finally`` block, so a subscribed SSE connection is guaranteed to
    terminate even if analysis fails outright (e.g. the VCF cannot be
    parsed) rather than hanging indefinitely.
    """
    is_htmx = request.headers.get("HX-Request") == "true"
    result_template = "_result_content.html" if is_htmx else "result.html"
    error_template = "_error_content.html" if is_htmx else "error.html"
    resolved_job_id = _parse_optional_job_id(job_id)

    try:
        try:
            tmp_path = await _save_upload_to_tempfile(vcf_file, settings.api_max_upload_bytes)
        except ValueError as exc:
            return templates.TemplateResponse(
                request,
                error_template,
                {"message": _clinical_error_message(exc)},
                status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            )

        try:
            variants, _stats = parse_vcf(tmp_path, expected_build=GenomeBuild(genome_build))
        except (VariantParseError, ValueError) as exc:
            return templates.TemplateResponse(
                request,
                error_template,
                {"message": _clinical_error_message(exc)},
                status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            )
        finally:
            tmp_path.unlink(missing_ok=True)

        response_job_id = resolved_job_id or str(uuid.uuid4())
        job_dir = Path(settings.web_reports_dir) / response_job_id
        orchestrator = ReportOrchestrator(job_dir)

        results: list[dict[str, object]] = []
        errors: list[str] = []
        total = len(variants)
        for index_, variant in enumerate(variants, start=1):
            on_progress = None
            if resolved_job_id is not None:
                on_progress = _make_progress_notifier(
                    _get_or_create_progress_queue(resolved_job_id), index_, total
                )
            try:
                report = await build_evidence_report(
                    variant,
                    settings=settings,
                    include_clinvar=include_clinvar,
                    on_progress=on_progress,
                )
            except AlleleIQError as exc:
                logger.error(
                    "web variant analysis failed",
                    extra={"variant_id": variant.variant_id, "error": str(exc)},
                )
                errors.append(f"{variant.variant_id}: {exc}")
                continue
            orchestrator.generate(report, formats=("json", "html"))
            results.append({"report": report, **_report_urls(response_job_id, report)})

        return templates.TemplateResponse(
            request,
            result_template,
            {"results": results, "errors": errors, "job_id": response_job_id},
            status_code=status.HTTP_200_OK,
        )
    finally:
        if resolved_job_id is not None:
            _get_or_create_progress_queue(resolved_job_id).put_nowait(_PROGRESS_SENTINEL)


@router.get(
    "/analyse/{job_id}/progress",
    summary="Server-Sent Events stream of real per-variant pipeline progress",
)
async def analyse_progress(job_id: str) -> StreamingResponse:
    """Stream progress events for one analysis job as Server-Sent Events.

    Each event's ``data:`` payload is a stage description sourced
    directly from :mod:`alleleiq.pipeline`'s ``on_progress`` callback
    (e.g. ``"variant 1/2: querying gnomAD"``) -- not a fabricated or
    time-based fake progress bar. The stream ends (after one final
    ``event: complete`` message) once the corresponding
    ``POST /analyse`` request finishes, whether it succeeded or failed.
    """
    try:
        uuid.UUID(job_id)
    except ValueError as exc:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="not found") from exc

    queue = _get_or_create_progress_queue(job_id)

    async def event_stream() -> AsyncIterator[str]:
        try:
            while True:
                message = await queue.get()
                if message == _PROGRESS_SENTINEL:
                    yield "event: complete\ndata: complete\n\n"
                    return
                yield f"data: {message}\n\n"
        finally:
            _progress_queues.pop(job_id, None)

    return StreamingResponse(event_stream(), media_type="text/event-stream")


@router.get("/reports/{job_id}/{filename}", summary="Download a generated report file")
async def download_report(
    job_id: str, filename: str, settings: Settings = Depends(get_settings_dependency)
) -> FileResponse:
    """Serve one report file from a completed analysis job.

    Security: ``job_id`` must parse as a UUID (rejecting anything that
    could be interpreted as a path segment), and ``filename`` is
    rejected outright if it contains a path separator, is exactly
    ``.`` / ``..``, or does not end in an allowed report suffix --
    only after those checks does the path get resolved and verified to
    still live inside the job's own directory, as defence in depth
    against path traversal.
    """
    try:
        uuid.UUID(job_id)
    except ValueError as exc:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="not found") from exc

    if (
        "/" in filename
        or "\\" in filename
        or filename in (".", "..")
        or not filename.endswith(_ALLOWED_REPORT_SUFFIXES)
    ):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="not found")

    job_dir = (Path(settings.web_reports_dir) / job_id).resolve()
    target = (job_dir / filename).resolve()
    if not target.is_relative_to(job_dir) or not target.is_file():
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="not found")

    media_type = "application/json" if filename.endswith(".json") else "text/html"
    return FileResponse(target, media_type=media_type, filename=filename)
