"""Bayesian evidence-point accumulation (Tavtigian et al. 2020), as adopted by ACGS 2024.

Point weights are keyed on :class:`~alleleiq.models.evidence.EvidenceStrength`,
not on the ACMG code itself — the whole point of decoupling code from
strength in :mod:`alleleiq.models.evidence` is that ACGS 2024 routinely
adjusts a criterion's applied strength (e.g. PM2 downgraded to
Supporting, PVS1 downgraded to Strong or Supporting) without changing
which criterion fired, and the point total must follow the applied
strength, not a fixed per-code weight.
"""

from __future__ import annotations

from alleleiq.models.evidence import (
    ACMGCode,
    BayesianScore,
    ClassificationConfidence,
    CriterionEvaluation,
    EvidenceStrength,
    PreliminaryClassification,
)
from alleleiq.models.population import OverallCoverageAssessment

# Tavtigian et al. (2020) evidence-point weights, as adopted by ACGS 2024.
# EvidenceStrength.BENIGN_MODERATE has no corresponding standard ACMG/AMP
# criterion code (Richards et al. 2015 defines no "moderate" benign tier);
# it is included here at -2, mirroring PM's +2, purely for symmetry should
# any future criterion ever be evaluated at that strength.
_STRENGTH_POINTS: dict[EvidenceStrength, int] = {
    EvidenceStrength.VERY_STRONG: 8,
    EvidenceStrength.STRONG: 4,
    EvidenceStrength.MODERATE: 2,
    EvidenceStrength.SUPPORTING: 1,
    EvidenceStrength.BENIGN_SUPPORTING: -1,
    EvidenceStrength.BENIGN_MODERATE: -2,
    EvidenceStrength.BENIGN_STRONG: -4,
    EvidenceStrength.BENIGN_STANDALONE: -8,
}

# Tavtigian (2020) / ACGS 2024 classification bands.
_PATHOGENIC_THRESHOLD = 10
_LIKELY_PATHOGENIC_THRESHOLD = 6
_UNCERTAIN_SIGNIFICANCE_THRESHOLD = 0
_LIKELY_BENIGN_THRESHOLD = -6

_PATHOGENIC_STRENGTHS = frozenset(
    {
        EvidenceStrength.VERY_STRONG,
        EvidenceStrength.STRONG,
        EvidenceStrength.MODERATE,
        EvidenceStrength.SUPPORTING,
    }
)
_BENIGN_STRENGTHS = frozenset(
    {
        EvidenceStrength.BENIGN_SUPPORTING,
        EvidenceStrength.BENIGN_MODERATE,
        EvidenceStrength.BENIGN_STRONG,
        EvidenceStrength.BENIGN_STANDALONE,
    }
)


def accumulate_evidence(
    evaluations: list[CriterionEvaluation],
    overall_coverage_assessment: OverallCoverageAssessment = OverallCoverageAssessment.COMPLETE,
) -> BayesianScore:
    """Aggregate criterion evaluations into a preliminary Bayesian evidence score.

    Args:
        evaluations: Every criterion evaluated for this variant, whether
            met, not met, or not assessable.
        overall_coverage_assessment: The variant's population coverage
            completeness (see :mod:`alleleiq.scoring.coverage`), factored
            into the classification confidence band.

    Returns:
        A :class:`BayesianScore` summarising the point total, which
        criteria contributed, which were evaluated but not met, which
        could not be assessed, any internal conflicts, and a preliminary
        classification with an explicit uncertainty band.
    """
    applied = [e for e in evaluations if e.met]
    not_met = [e for e in evaluations if not e.met and not e.requires_clinical_input]
    not_assessed = [e for e in evaluations if not e.met and e.requires_clinical_input]

    evidence_points = sum(_STRENGTH_POINTS[e.strength] for e in applied)
    classification = _classify(evidence_points)
    conflict_flags = _detect_conflicts(applied)
    confidence = _assess_confidence(
        conflict_flags=conflict_flags,
        not_assessed_count=len(not_assessed),
        applied_count=len(applied),
        overall_coverage_assessment=overall_coverage_assessment,
    )
    confidence_interval = _evidence_points_confidence_interval(evidence_points, not_assessed)

    return BayesianScore(
        evidence_points=evidence_points,
        criteria_applied=[e.code for e in applied],
        criteria_not_met=[e.code for e in not_met],
        criteria_not_assessed=[e.code for e in not_assessed],
        preliminary_classification=classification,
        classification_confidence=confidence,
        evidence_points_confidence_interval=confidence_interval,
        conflict_flags=conflict_flags,
    )


def _classify(points: int) -> PreliminaryClassification:
    if points >= _PATHOGENIC_THRESHOLD:
        return PreliminaryClassification.PATHOGENIC
    if points >= _LIKELY_PATHOGENIC_THRESHOLD:
        return PreliminaryClassification.LIKELY_PATHOGENIC
    if points >= _UNCERTAIN_SIGNIFICANCE_THRESHOLD:
        return PreliminaryClassification.UNCERTAIN_SIGNIFICANCE
    if points >= _LIKELY_BENIGN_THRESHOLD:
        return PreliminaryClassification.LIKELY_BENIGN
    return PreliminaryClassification.BENIGN


def _detect_conflicts(applied: list[CriterionEvaluation]) -> list[str]:
    """Flag internally contradictory evidence, e.g. BA1 alongside pathogenic criteria.

    A variant should not simultaneously present standalone-benign and
    pathogenic evidence; when it does, that is not a scoring edge case
    to quietly resolve — it is a signal the automated evaluation itself
    may be wrong (e.g. a misannotated consequence) and must be flagged
    for clinical review rather than averaged away.
    """
    pathogenic = [e for e in applied if e.strength in _PATHOGENIC_STRENGTHS]
    benign = [e for e in applied if e.strength in _BENIGN_STRENGTHS]

    if not pathogenic or not benign:
        return []

    flags = []
    ba1_applied = any(e.code == ACMGCode.BA1 for e in benign)
    if ba1_applied:
        pathogenic_codes = ", ".join(sorted(e.code.value for e in pathogenic))
        flags.append(
            f"BA1 (benign standalone) applied alongside pathogenic evidence "
            f"({pathogenic_codes}); this combination should not occur and warrants "
            "review of the underlying annotations."
        )
    else:
        pathogenic_codes = ", ".join(sorted(e.code.value for e in pathogenic))
        benign_codes = ", ".join(sorted(e.code.value for e in benign))
        flags.append(
            f"Conflicting evidence: pathogenic criteria ({pathogenic_codes}) applied "
            f"alongside benign criteria ({benign_codes})."
        )
    return flags


def _assess_confidence(
    *,
    conflict_flags: list[str],
    not_assessed_count: int,
    applied_count: int,
    overall_coverage_assessment: OverallCoverageAssessment,
) -> ClassificationConfidence:
    if conflict_flags or overall_coverage_assessment == OverallCoverageAssessment.SEVERELY_LIMITED:
        return ClassificationConfidence.LOW
    if not_assessed_count > applied_count:
        return ClassificationConfidence.LOW
    if overall_coverage_assessment == OverallCoverageAssessment.PARTIAL or not_assessed_count > 0:
        return ClassificationConfidence.MODERATE
    return ClassificationConfidence.HIGH


def _evidence_points_confidence_interval(
    evidence_points: int, not_assessed: list[CriterionEvaluation]
) -> tuple[int, int]:
    """Best-/worst-case point range if every unassessed criterion resolved.

    This is not a statistical confidence interval; it is an explicit
    uncertainty band answering "how far could the score move if every
    criterion currently marked requires_clinical_input turned out, on
    review, to be met at its evaluated strength?".
    """
    optimistic_gain = sum(max(_STRENGTH_POINTS[e.strength], 0) for e in not_assessed)
    pessimistic_loss = sum(min(_STRENGTH_POINTS[e.strength], 0) for e in not_assessed)
    return evidence_points + pessimistic_loss, evidence_points + optimistic_gain
