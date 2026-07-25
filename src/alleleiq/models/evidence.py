"""ACMG/AMP evidence criterion and Bayesian scoring models.

Modelled per Richards et al. (2015), the ACGS 2024 UK best-practice
guidelines, and the Tavtigian et al. (2020) Bayesian point-scoring
framework. A criterion's *code* (e.g. ``PM2``) and its *applied strength*
(e.g. ``SUPPORTING``) are kept as separate fields rather than encoding
strength modifiers as distinct enum members (``PM2_SUPPORTING``), since
ACGS 2024 routinely adjusts the strength of a criterion's application
without changing which criterion is being applied.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field


class EvidenceStrength(StrEnum):
    """Strength at which a criterion was applied.

    Point weights for each strength (per Tavtigian et al. 2020) are
    defined in :mod:`alleleiq.scoring.bayesian`, not here, keeping this
    module free of scoring policy.
    """

    VERY_STRONG = "VERY_STRONG"
    STRONG = "STRONG"
    MODERATE = "MODERATE"
    SUPPORTING = "SUPPORTING"
    BENIGN_SUPPORTING = "BENIGN_SUPPORTING"
    BENIGN_MODERATE = "BENIGN_MODERATE"
    BENIGN_STRONG = "BENIGN_STRONG"
    BENIGN_STANDALONE = "BENIGN_STANDALONE"


class ACMGCode(StrEnum):
    """ACMG/AMP evidence criterion codes (Richards et al. 2015)."""

    PVS1 = "PVS1"
    PS1 = "PS1"
    PS2 = "PS2"
    PS3 = "PS3"
    PS4 = "PS4"
    PM1 = "PM1"
    PM2 = "PM2"
    PM3 = "PM3"
    PM4 = "PM4"
    PM5 = "PM5"
    PM6 = "PM6"
    PP1 = "PP1"
    PP2 = "PP2"
    PP3 = "PP3"
    PP4 = "PP4"
    PP5 = "PP5"
    BA1 = "BA1"
    BS1 = "BS1"
    BS2 = "BS2"
    BS3 = "BS3"
    BS4 = "BS4"
    BP1 = "BP1"
    BP2 = "BP2"
    BP3 = "BP3"
    BP4 = "BP4"
    BP5 = "BP5"
    BP6 = "BP6"
    BP7 = "BP7"


class PreliminaryClassification(StrEnum):
    """Five-tier ACMG/AMP classification band.

    Named "preliminary" throughout AlleleIQ to reinforce that this is an
    evidence summary, not a clinical classification.
    """

    PATHOGENIC = "PATHOGENIC"
    LIKELY_PATHOGENIC = "LIKELY_PATHOGENIC"
    UNCERTAIN_SIGNIFICANCE = "UNCERTAIN_SIGNIFICANCE"
    LIKELY_BENIGN = "LIKELY_BENIGN"
    BENIGN = "BENIGN"


class ClassificationConfidence(StrEnum):
    """Qualitative confidence band for a preliminary classification.

    Reflects the completeness of population coverage and the presence of
    conflicting evidence, not the statistical strength of any single
    criterion. See :mod:`alleleiq.scoring.bayesian` for the precise rule
    used to derive this value.
    """

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"


class CriterionEvaluation(BaseModel):
    """The outcome of evaluating a single ACMG/AMP criterion for one variant."""

    model_config = ConfigDict(frozen=True)

    code: ACMGCode
    strength: EvidenceStrength
    met: bool = Field(description="Whether the criterion's conditions were satisfied.")
    evidence_basis: str = Field(
        description="Human-readable explanation of why the criterion was, or was not, met."
    )
    data_source: str = Field(
        description="Evidence source, e.g. 'gnomAD v4', 'Ensembl VEP / REVEL', 'ClinVar'."
    )
    data_source_version: str | None = Field(default=None)
    requires_clinical_input: bool = Field(
        default=False,
        description="True when automated evaluation was inconclusive and needs clinical review.",
    )
    notes: str | None = Field(default=None)


class BayesianScore(BaseModel):
    """Aggregate Tavtigian (2020) evidence-point score for a variant.

    Point contributions and classification thresholds are computed by
    :mod:`alleleiq.scoring.bayesian`; this model is a pure data container
    for the result.
    """

    model_config = ConfigDict(frozen=True)

    evidence_points: int
    criteria_applied: list[ACMGCode] = Field(
        default_factory=list, description="Codes that contributed points."
    )
    criteria_not_met: list[ACMGCode] = Field(
        default_factory=list, description="Codes evaluated but not satisfied."
    )
    criteria_not_assessed: list[ACMGCode] = Field(
        default_factory=list,
        description="Codes requiring clinical input the automated pipeline lacks.",
    )
    preliminary_classification: PreliminaryClassification
    classification_confidence: ClassificationConfidence
    evidence_points_confidence_interval: tuple[int, int] | None = Field(
        default=None, description="95% confidence interval on evidence_points, if determinable."
    )
    conflict_flags: list[str] = Field(
        default_factory=list,
        description="Human-readable notes on conflicting evidence, e.g. BA1 alongside PVS1.",
    )
