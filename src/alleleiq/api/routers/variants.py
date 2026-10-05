"""Variant analysis endpoints."""

from __future__ import annotations

import tempfile
from pathlib import Path

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status

from alleleiq.api.dependencies import get_settings_dependency
from alleleiq.api.schemas import BatchAnalyseResponse, ErrorResponse
from alleleiq.config import Settings
from alleleiq.exceptions import AlleleIQError, VariantParseError
from alleleiq.logging import get_logger
from alleleiq.models.report import EvidenceReport
from alleleiq.models.variant import GenomeBuild, Variant
from alleleiq.parsers.vcf_parser import ParsingStatistics, parse_vcf
from alleleiq.pipeline import build_evidence_report

logger = get_logger(__name__)
router = APIRouter(prefix="/variants", tags=["variants"])


async def _save_upload_to_tempfile(upload: UploadFile, max_bytes: int) -> Path:
    """Persist an uploaded VCF to a temporary file, enforcing a size limit.

    Raises:
        HTTPException: 413 if the upload exceeds ``max_bytes``.
    """
    contents = await upload.read()
    if len(contents) > max_bytes:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail=f"upload exceeds the {max_bytes}-byte limit",
        )

    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as handle:
        handle.write(contents)
        tmp_path = Path(handle.name)
    return tmp_path


def _parse_uploaded_vcf(path: Path, genome_build: str) -> tuple[list[Variant], ParsingStatistics]:
    try:
        return parse_vcf(path, expected_build=GenomeBuild(genome_build))
    except VariantParseError as exc:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT, detail=str(exc)
        ) from exc
    except ValueError as exc:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            detail=f"unsupported genome_build {genome_build!r}: {exc}",
        ) from exc
    finally:
        path.unlink(missing_ok=True)


@router.post(
    "/analyse",
    response_model=EvidenceReport,
    summary="Analyse a single-variant VCF upload",
    responses={
        status.HTTP_422_UNPROCESSABLE_CONTENT: {"model": ErrorResponse},
        status.HTTP_502_BAD_GATEWAY: {"model": ErrorResponse},
    },
)
async def analyse_variant(
    vcf_file: UploadFile = File(..., description="A VCF file containing exactly one variant."),
    genome_build: str = Form("GRCh38", description="Genome build; GRCh38 only."),
    include_clinvar: bool = Form(
        True,
        description="Whether to query ClinVar. PS1, PM5, and PM1's density signal need it.",
    ),
    settings: Settings = Depends(get_settings_dependency),
) -> EvidenceReport:
    """Analyse a single variant and return its complete evidence report.

    Returns 422 for a malformed VCF or a VCF containing anything other
    than exactly one variant (use ``/variants/analyse/batch`` for
    multi-variant files), and 502 if an upstream evidence source
    (gnomAD, VEP, or ClinVar) could not be reached.
    """
    tmp_path = await _save_upload_to_tempfile(vcf_file, settings.api_max_upload_bytes)
    variants, _stats = _parse_uploaded_vcf(tmp_path, genome_build)

    if len(variants) != 1:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            detail=(
                f"expected exactly one variant, found {len(variants)}; "
                "use /variants/analyse/batch for multi-variant VCFs"
            ),
        )

    try:
        return await build_evidence_report(
            variants[0], settings=settings, include_clinvar=include_clinvar
        )
    except AlleleIQError as exc:
        raise HTTPException(status_code=status.HTTP_502_BAD_GATEWAY, detail=str(exc)) from exc


@router.post(
    "/analyse/batch",
    response_model=BatchAnalyseResponse,
    summary="Analyse a multi-variant VCF upload",
    responses={status.HTTP_422_UNPROCESSABLE_CONTENT: {"model": ErrorResponse}},
)
async def analyse_batch(
    vcf_file: UploadFile = File(..., description="A VCF file containing one or more variants."),
    genome_build: str = Form("GRCh38", description="Genome build; GRCh38 only."),
    include_clinvar: bool = Form(
        True,
        description="Whether to query ClinVar. PS1, PM5, and PM1's density signal need it.",
    ),
    settings: Settings = Depends(get_settings_dependency),
) -> BatchAnalyseResponse:
    """Analyse every variant in a VCF and return one report per variant.

    A per-variant evidence-source failure does not abort the batch: it
    is recorded in ``errors`` and analysis continues with the remaining
    variants, so one unreachable upstream query cannot silently drop an
    entire batch's results.
    """
    tmp_path = await _save_upload_to_tempfile(vcf_file, settings.api_max_upload_bytes)
    variants, _stats = _parse_uploaded_vcf(tmp_path, genome_build)

    reports: list[EvidenceReport] = []
    errors: list[str] = []
    for variant in variants:
        try:
            reports.append(
                await build_evidence_report(
                    variant, settings=settings, include_clinvar=include_clinvar
                )
            )
        except AlleleIQError as exc:
            logger.error(
                "batch variant analysis failed",
                extra={"variant_id": variant.variant_id, "error": str(exc)},
            )
            errors.append(f"{variant.variant_id}: {exc}")

    return BatchAnalyseResponse(
        total_variants=len(variants),
        successful=len(reports),
        reports=reports,
        errors=errors,
    )
