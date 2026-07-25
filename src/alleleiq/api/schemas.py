"""Request and response schemas for the AlleleIQ REST API.

:class:`~alleleiq.models.report.EvidenceReport` itself is used directly
as the response model for the single-variant endpoint — it is already a
complete, versioned Pydantic model, and duplicating its fields into a
separate API-layer schema would be exactly the kind of drift-prone
duplication the project avoids elsewhere (see
:mod:`alleleiq.reports.json_report`, which renders from the same
report rather than a parallel representation). The schemas here cover
what genuinely differs at the API boundary: batch responses and health
signalling.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field

from alleleiq.models.report import EvidenceReport


class BatchAnalyseResponse(BaseModel):
    """Response for a multi-variant VCF analysis request."""

    model_config = ConfigDict(frozen=True)

    total_variants: int = Field(description="Number of variants parsed from the uploaded VCF.")
    successful: int = Field(description="Number of variants successfully analysed.")
    reports: list[EvidenceReport] = Field(
        description="One EvidenceReport per successfully analysed variant."
    )
    errors: list[str] = Field(
        default_factory=list,
        description="Human-readable errors for variants that failed analysis, if any.",
    )


class HealthResponse(BaseModel):
    """Liveness check response."""

    model_config = ConfigDict(frozen=True)

    status: str = Field(description="Always 'healthy' if the process is responding at all.")
    version: str


class ReadinessResponse(BaseModel):
    """Readiness check response: upstream evidence-source connectivity."""

    model_config = ConfigDict(frozen=True)

    status: str = Field(
        description="'ready' if all upstream services are reachable, else 'degraded'."
    )
    gnomad: bool
    vep: bool
    clinvar: bool


class ErrorResponse(BaseModel):
    """Standard error response shape, matching FastAPI's default HTTPException body."""

    model_config = ConfigDict(frozen=True)

    detail: str
