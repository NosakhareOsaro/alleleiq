"""Public Pydantic v2 data models shared across AlleleIQ.

No raw ``dict`` is passed across module boundaries in AlleleIQ's public
interfaces; every cross-module data exchange uses a model defined here.
"""

from alleleiq.models.evidence import (
    ACMGCode,
    BayesianScore,
    ClassificationConfidence,
    CriterionEvaluation,
    EvidenceStrength,
    PreliminaryClassification,
)
from alleleiq.models.population import (
    AncestryGroup,
    CoverageFlagLevel,
    CoverageGap,
    OverallCoverageAssessment,
    PopulationFrequency,
)
from alleleiq.models.report import (
    CLINICAL_DISCLAIMER,
    AuditTrail,
    CoverageSummary,
    DatabaseQueryRecord,
    EvidenceReport,
    ReportMetadata,
)
from alleleiq.models.variant import GenomeBuild, GenomicCoordinates, Variant, VariantType

__all__ = [
    "CLINICAL_DISCLAIMER",
    "ACMGCode",
    "AncestryGroup",
    "AuditTrail",
    "BayesianScore",
    "ClassificationConfidence",
    "CoverageFlagLevel",
    "CoverageGap",
    "CoverageSummary",
    "CriterionEvaluation",
    "DatabaseQueryRecord",
    "EvidenceReport",
    "EvidenceStrength",
    "GenomeBuild",
    "GenomicCoordinates",
    "OverallCoverageAssessment",
    "PopulationFrequency",
    "PreliminaryClassification",
    "ReportMetadata",
    "Variant",
    "VariantType",
]
