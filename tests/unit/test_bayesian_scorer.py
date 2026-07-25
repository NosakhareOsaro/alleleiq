"""Unit tests for alleleiq.scoring.bayesian.accumulate_evidence."""

from __future__ import annotations

from alleleiq.models.evidence import (
    ACMGCode,
    ClassificationConfidence,
    CriterionEvaluation,
    EvidenceStrength,
    PreliminaryClassification,
)
from alleleiq.models.population import OverallCoverageAssessment
from alleleiq.scoring.bayesian import accumulate_evidence


def _evaluation(
    code: ACMGCode,
    strength: EvidenceStrength,
    met: bool = True,
    requires_clinical_input: bool = False,
) -> CriterionEvaluation:
    return CriterionEvaluation(
        code=code,
        strength=strength,
        met=met,
        evidence_basis="test",
        data_source="test",
        requires_clinical_input=requires_clinical_input,
    )


class TestKnownPathogenicVariant:
    def test_pvs1_plus_pm2_plus_pp3_reaches_pathogenic(self) -> None:
        evaluations = [
            _evaluation(ACMGCode.PVS1, EvidenceStrength.VERY_STRONG),
            _evaluation(ACMGCode.PM2, EvidenceStrength.SUPPORTING),
            _evaluation(ACMGCode.PP3, EvidenceStrength.SUPPORTING),
        ]

        score = accumulate_evidence(evaluations)

        assert score.evidence_points == 10
        assert score.preliminary_classification == PreliminaryClassification.PATHOGENIC
        assert score.conflict_flags == []

    def test_strong_plus_two_moderate_reaches_likely_pathogenic(self) -> None:
        evaluations = [
            _evaluation(ACMGCode.PS1, EvidenceStrength.STRONG),
            _evaluation(ACMGCode.PM1, EvidenceStrength.MODERATE),
            _evaluation(ACMGCode.PM2, EvidenceStrength.MODERATE),
        ]

        score = accumulate_evidence(evaluations)

        assert score.evidence_points == 8
        assert score.preliminary_classification == PreliminaryClassification.LIKELY_PATHOGENIC


class TestKnownBenignVariant:
    def test_ba1_alone_reaches_benign(self) -> None:
        evaluations = [_evaluation(ACMGCode.BA1, EvidenceStrength.BENIGN_STANDALONE)]

        score = accumulate_evidence(evaluations)

        assert score.evidence_points == -8
        assert score.preliminary_classification == PreliminaryClassification.BENIGN

    def test_two_benign_strong_reaches_likely_benign_boundary(self) -> None:
        evaluations = [
            _evaluation(ACMGCode.BS1, EvidenceStrength.BENIGN_STRONG),
            _evaluation(ACMGCode.BP4, EvidenceStrength.BENIGN_SUPPORTING),
        ]

        score = accumulate_evidence(evaluations)

        assert score.evidence_points == -5
        assert score.preliminary_classification == PreliminaryClassification.LIKELY_BENIGN


class TestVUSWithConflictingEvidence:
    def test_ba1_and_pvs1_both_met_flags_conflict_and_lands_at_vus(self) -> None:
        evaluations = [
            _evaluation(ACMGCode.BA1, EvidenceStrength.BENIGN_STANDALONE),
            _evaluation(ACMGCode.PVS1, EvidenceStrength.VERY_STRONG),
        ]

        score = accumulate_evidence(evaluations)

        assert score.evidence_points == 0
        assert score.preliminary_classification == PreliminaryClassification.UNCERTAIN_SIGNIFICANCE
        assert len(score.conflict_flags) == 1
        assert "BA1" in score.conflict_flags[0]
        assert score.classification_confidence == ClassificationConfidence.LOW

    def test_non_ba1_pathogenic_benign_conflict_still_flagged(self) -> None:
        evaluations = [
            _evaluation(ACMGCode.BS3, EvidenceStrength.BENIGN_STRONG),
            _evaluation(ACMGCode.PS3, EvidenceStrength.STRONG),
        ]

        score = accumulate_evidence(evaluations)

        assert len(score.conflict_flags) == 1
        assert "BA1" not in score.conflict_flags[0]


class TestNotAssessedAndConfidence:
    def test_unassessed_criteria_populate_not_assessed_list(self) -> None:
        evaluations = [
            _evaluation(
                ACMGCode.PM2,
                EvidenceStrength.SUPPORTING,
                met=False,
                requires_clinical_input=True,
            ),
        ]

        score = accumulate_evidence(evaluations)

        assert score.criteria_not_assessed == [ACMGCode.PM2]
        assert score.criteria_applied == []
        assert score.classification_confidence == ClassificationConfidence.LOW

    def test_not_met_without_clinical_input_is_confidently_not_met(self) -> None:
        evaluations = [
            _evaluation(ACMGCode.BA1, EvidenceStrength.BENIGN_STANDALONE, met=False),
        ]

        score = accumulate_evidence(evaluations)

        assert score.criteria_not_met == [ACMGCode.BA1]
        assert score.criteria_not_assessed == []

    def test_severely_limited_coverage_forces_low_confidence(self) -> None:
        evaluations = [_evaluation(ACMGCode.PVS1, EvidenceStrength.VERY_STRONG)]

        score = accumulate_evidence(
            evaluations, overall_coverage_assessment=OverallCoverageAssessment.SEVERELY_LIMITED
        )

        assert score.classification_confidence == ClassificationConfidence.LOW

    def test_complete_coverage_no_conflicts_no_gaps_is_high_confidence(self) -> None:
        evaluations = [_evaluation(ACMGCode.PVS1, EvidenceStrength.VERY_STRONG)]

        score = accumulate_evidence(
            evaluations, overall_coverage_assessment=OverallCoverageAssessment.COMPLETE
        )

        assert score.classification_confidence == ClassificationConfidence.HIGH

    def test_confidence_interval_reflects_unassessed_upside_and_downside(self) -> None:
        evaluations = [
            _evaluation(
                ACMGCode.PM2,
                EvidenceStrength.SUPPORTING,
                met=False,
                requires_clinical_input=True,
            ),
            _evaluation(
                ACMGCode.BP4,
                EvidenceStrength.BENIGN_SUPPORTING,
                met=False,
                requires_clinical_input=True,
            ),
        ]

        score = accumulate_evidence(evaluations)

        assert score.evidence_points_confidence_interval == (-1, 1)
