"""Routes for the minimal server-rendered web interface.

Every route here delegates the actual analysis to
:func:`alleleiq.pipeline.build_evidence_report` -- the same function
:mod:`alleleiq.cli` and :mod:`alleleiq.api.routers.variants` call -- so
this interface is a presentation layer, not a second implementation of
variant analysis.
"""

from __future__ import annotations

import tempfile
import uuid
from pathlib import Path

from fastapi import APIRouter, Depends, File, Form, HTTPException, Request, UploadFile, status
from fastapi.responses import FileResponse, HTMLResponse
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
from alleleiq.pipeline import build_evidence_report
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
    settings: Settings = Depends(get_settings_dependency),
) -> HTMLResponse:
    """Parse and analyse an uploaded VCF, returning a rendered result page or fragment.

    Renders a full page (extending ``base.html``) for a plain form
    submission, or just the result fragment when the request carries
    HTMX's ``HX-Request`` header -- both are backed by the same
    ``_result_content.html`` partial, so the two response shapes can
    never show different results for the same analysis.
    """
    is_htmx = request.headers.get("HX-Request") == "true"
    result_template = "_result_content.html" if is_htmx else "result.html"
    error_template = "_error_content.html" if is_htmx else "error.html"

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

    job_id = str(uuid.uuid4())
    job_dir = Path(settings.web_reports_dir) / job_id
    orchestrator = ReportOrchestrator(job_dir)

    results: list[dict[str, object]] = []
    errors: list[str] = []
    for variant in variants:
        try:
            report = await build_evidence_report(
                variant, settings=settings, include_clinvar=include_clinvar
            )
        except AlleleIQError as exc:
            logger.error(
                "web variant analysis failed",
                extra={"variant_id": variant.variant_id, "error": str(exc)},
            )
            errors.append(f"{variant.variant_id}: {exc}")
            continue
        orchestrator.generate(report, formats=("json", "html"))
        results.append({"report": report, **_report_urls(job_id, report)})

    return templates.TemplateResponse(
        request,
        result_template,
        {"results": results, "errors": errors, "job_id": job_id},
        status_code=status.HTTP_200_OK,
    )


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
