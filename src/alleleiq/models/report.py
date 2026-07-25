"""Evidence report models, including the mandatory clinical disclaimer.

``CLINICAL_DISCLAIMER`` is the primary patient-safety mechanism at
prototype stage: it is typed as a ``Literal`` on
:class:`EvidenceReport.clinical_disclaimer`, which means no code path —
however it constructs a report — can produce a report object with a
missing, blank, or altered disclaimer without failing validation.
"""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Final, Literal
from uuid import UUID, uuid4

from pydantic import BaseModel, ConfigDict, Field

from alleleiq.models.evidence import BayesianScore, CriterionEvaluation
from alleleiq.models.population import (
    AncestryGroup,
    CoverageGap,
    OverallCoverageAssessment,
    PopulationFrequency,
)
from alleleiq.models.variant import Variant

CLINICAL_DISCLAIMER: Final = (
    "PRELIMINARY EVIDENCE SUMMARY — FOR CLINICAL SCIENTIST REVIEW ONLY. "
    "This output does not constitute a variant classification. "
    "Classification must be performed by a qualified Clinical Genetic "
    "Scientist following ACGS 2024 guidelines with full consideration "
    "of clinical phenotype, family history, and all available evidence."
)


class ReportMetadata(BaseModel):
    """Provenance metadata for a single evidence report."""

    model_config = ConfigDict(frozen=True)

    report_id: UUID = Field(default_factory=uuid4)
    generated_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    alleleiq_version: str
    database_versions: dict[str, str] = Field(
        default_factory=dict,
        description="e.g. {'gnomad': 'gnomad_r4', 'clinvar': '2026-06-01', 'vep': '112'}.",
    )
    query_timestamps: dict[str, datetime] = Field(
        default_factory=dict,
        description="When each upstream database was queried, keyed by database name.",
    )


class DatabaseQueryRecord(BaseModel):
    """A single logged query to an external evidence source."""

    model_config = ConfigDict(frozen=True)

    database: str = Field(description="e.g. 'gnomAD', 'ClinVar', 'VEP'.")
    query_description: str
    parameters: dict[str, str] = Field(default_factory=dict)
    timestamp: datetime
    success: bool
    error_message: str | None = None


class AuditTrail(BaseModel):
    """The complete, ordered record of external queries made for a report."""

    model_config = ConfigDict(frozen=True)

    queries: list[DatabaseQueryRecord] = Field(default_factory=list)


class CoverageSummary(BaseModel):
    """Report-level roll-up of per-ancestry coverage gaps."""

    model_config = ConfigDict(frozen=True)

    populations_with_sufficient_data: list[AncestryGroup] = Field(default_factory=list)
    populations_with_insufficient_data: list[AncestryGroup] = Field(default_factory=list)
    populations_with_no_data: list[AncestryGroup] = Field(default_factory=list)
    overall_coverage_assessment: OverallCoverageAssessment
    clinical_advisory: str = Field(
        description="Prose summary of interpretation implications for the reviewer."
    )


class EvidenceReport(BaseModel):
    """The complete structured output of one variant's evidence aggregation.

    This is the canonical in-memory representation from which both the
    JSON and HTML report renderers (:mod:`alleleiq.reports`) derive their
    output, so the two formats never drift out of sync.
    """

    model_config = ConfigDict(frozen=True)

    variant: Variant
    metadata: ReportMetadata
    population_frequencies: dict[AncestryGroup, PopulationFrequency] = Field(default_factory=dict)
    coverage_gaps: dict[AncestryGroup, CoverageGap] = Field(default_factory=dict)
    coverage_summary: CoverageSummary
    criterion_evaluations: list[CriterionEvaluation] = Field(default_factory=list)
    bayesian_score: BayesianScore
    audit_trail: AuditTrail
    # NOTE: mypy cannot resolve `Literal[CLINICAL_DISCLAIMER]` (Literal parameters must be
    # literal expressions, not variable references — see python/mypy#8917), so the text is
    # duplicated here. tests/unit/test_report_models.py asserts these never drift apart.
    clinical_disclaimer: Literal[
        "PRELIMINARY EVIDENCE SUMMARY — FOR CLINICAL SCIENTIST REVIEW ONLY. "
        "This output does not constitute a variant classification. "
        "Classification must be performed by a qualified Clinical Genetic "
        "Scientist following ACGS 2024 guidelines with full consideration "
        "of clinical phenotype, family history, and all available evidence."
    ] = CLINICAL_DISCLAIMER
