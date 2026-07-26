"""Unit tests for alleleiq.validation: gold-standard evaluation and outcome model."""

from __future__ import annotations

import pytest

import alleleiq.validation as validation_module
from alleleiq.models.evidence import (
    ACMGCode,
    BayesianScore,
    ClassificationConfidence,
    PreliminaryClassification,
)
from alleleiq.models.population import OverallCoverageAssessment
from alleleiq.models.report import AuditTrail, CoverageSummary, EvidenceReport, ReportMetadata
from alleleiq.validation import (
    GoldStandardVariant,
    Outcome,
    _direction,
    evaluate_variant,
    render_results,
)


def _report_with_classification(classification: PreliminaryClassification) -> EvidenceReport:
    from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType

    coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")
    return EvidenceReport(
        variant=Variant(coordinates=coords, variant_type=VariantType.SNV, gene_symbol="TESTGENE"),
        metadata=ReportMetadata(alleleiq_version="0.1.0-alpha"),
        coverage_summary=CoverageSummary(
            overall_coverage_assessment=OverallCoverageAssessment.COMPLETE,
            clinical_advisory="ok",
        ),
        bayesian_score=BayesianScore(
            evidence_points=0,
            criteria_applied=[ACMGCode.PP3],
            preliminary_classification=classification,
            classification_confidence=ClassificationConfidence.MODERATE,
        ),
        audit_trail=AuditTrail(),
    )


_GOLD = GoldStandardVariant(
    label="TESTGENE test variant",
    chrom="17",
    pos=43094692,
    ref="G",
    alt="A",
    gene_symbol="TESTGENE",
    clinvar_accession="VCV999999999",
    expected_classification=PreliminaryClassification.PATHOGENIC,
)


class TestDirection:
    @pytest.mark.parametrize(
        ("classification", "expected"),
        [
            (PreliminaryClassification.PATHOGENIC, "PATHOGENIC"),
            (PreliminaryClassification.LIKELY_PATHOGENIC, "PATHOGENIC"),
            (PreliminaryClassification.UNCERTAIN_SIGNIFICANCE, "UNCERTAIN"),
            (PreliminaryClassification.LIKELY_BENIGN, "BENIGN"),
            (PreliminaryClassification.BENIGN, "BENIGN"),
        ],
    )
    def test_direction_mapping(
        self, classification: PreliminaryClassification, expected: str
    ) -> None:
        assert _direction(classification) == expected


class TestEvaluateVariant:
    async def test_concordant_when_directions_match(self, monkeypatch: pytest.MonkeyPatch) -> None:
        async def fake_build(*args: object, **kwargs: object) -> EvidenceReport:
            return _report_with_classification(PreliminaryClassification.LIKELY_PATHOGENIC)

        monkeypatch.setattr(validation_module, "build_evidence_report", fake_build)

        result = await evaluate_variant(_GOLD)

        assert result.outcome == Outcome.CONCORDANT
        assert result.explanation is None

    async def test_conservative_when_actual_is_uncertain(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        async def fake_build(*args: object, **kwargs: object) -> EvidenceReport:
            return _report_with_classification(PreliminaryClassification.UNCERTAIN_SIGNIFICANCE)

        monkeypatch.setattr(validation_module, "build_evidence_report", fake_build)

        result = await evaluate_variant(_GOLD)

        assert result.outcome == Outcome.CONSERVATIVE
        assert result.explanation is not None

    async def test_dangerous_when_direction_is_opposite(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        async def fake_build(*args: object, **kwargs: object) -> EvidenceReport:
            return _report_with_classification(PreliminaryClassification.BENIGN)

        monkeypatch.setattr(validation_module, "build_evidence_report", fake_build)

        result = await evaluate_variant(_GOLD)

        assert result.outcome == Outcome.DANGEROUS
        assert "expected PATHOGENIC" in (result.explanation or "")


class TestRenderResults:
    def test_render_does_not_raise_with_mixed_outcomes(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        from alleleiq.validation import ValidationResult

        results = [
            ValidationResult(
                gold_standard=_GOLD,
                actual_classification=PreliminaryClassification.PATHOGENIC,
                evidence_points=10,
                criteria_applied=["PVS1"],
                outcome=Outcome.CONCORDANT,
                explanation=None,
            ),
            ValidationResult(
                gold_standard=_GOLD,
                actual_classification=PreliminaryClassification.UNCERTAIN_SIGNIFICANCE,
                evidence_points=0,
                criteria_applied=[],
                outcome=Outcome.CONSERVATIVE,
                explanation="expected PATHOGENIC but got VUS",
            ),
            ValidationResult(
                gold_standard=_GOLD,
                actual_classification=PreliminaryClassification.BENIGN,
                evidence_points=-8,
                criteria_applied=["BA1"],
                outcome=Outcome.DANGEROUS,
                explanation="expected PATHOGENIC but got BENIGN",
            ),
        ]

        render_results(results)

        captured = capsys.readouterr()
        assert "Dangerous discordances" in captured.out
