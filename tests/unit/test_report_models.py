"""Unit tests for alleleiq.models: variant, population, evidence, report."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from alleleiq.models import (
    CLINICAL_DISCLAIMER,
    AncestryGroup,
    AuditTrail,
    BayesianScore,
    ClassificationConfidence,
    CoverageSummary,
    EvidenceReport,
    GenomicCoordinates,
    PopulationFrequency,
    PreliminaryClassification,
    ReportMetadata,
    Variant,
    VariantType,
)


class TestGenomicCoordinates:
    @pytest.mark.parametrize(
        ("raw_chrom", "expected"),
        [
            ("17", "17"),
            ("chr17", "17"),
            ("Chr17", "17"),
            ("X", "X"),
            ("chrX", "X"),
            ("chrM", "MT"),
            ("MT", "MT"),
            ("chrMT", "MT"),
        ],
    )
    def test_chromosome_normalisation(self, raw_chrom: str, expected: str) -> None:
        coords = GenomicCoordinates(chrom=raw_chrom, pos=100, ref="A", alt="G")
        assert coords.chrom == expected

    def test_rejects_unrecognised_chromosome(self) -> None:
        with pytest.raises(ValidationError):
            GenomicCoordinates(chrom="99", pos=100, ref="A", alt="G")

    def test_alleles_are_uppercased(self) -> None:
        coords = GenomicCoordinates(chrom="1", pos=100, ref="a", alt="g")
        assert coords.ref == "A"
        assert coords.alt == "G"

    def test_rejects_invalid_allele_characters(self) -> None:
        with pytest.raises(ValidationError):
            GenomicCoordinates(chrom="1", pos=100, ref="A", alt="Z")

    def test_accepts_symbolic_structural_allele(self) -> None:
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="<DEL>")
        assert coords.alt == "<DEL>"

    def test_variant_id_format(self) -> None:
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")
        assert coords.variant_id == "17-43094692-G-A"

    def test_coordinates_are_immutable(self) -> None:
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="G")
        with pytest.raises(ValidationError):
            coords.pos = 200  # type: ignore[misc]


class TestVariant:
    def test_snv_with_nucleotide_allele_is_valid(
        self, brca1_coordinates: GenomicCoordinates
    ) -> None:
        variant = Variant(coordinates=brca1_coordinates, variant_type=VariantType.SNV)
        assert variant.variant_id == "17-43094692-G-A"

    def test_snv_type_with_symbolic_allele_is_rejected(self) -> None:
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="<DEL>")
        with pytest.raises(ValidationError):
            Variant(coordinates=coords, variant_type=VariantType.SNV)

    def test_structural_type_with_nucleotide_allele_is_rejected(
        self, brca1_coordinates: GenomicCoordinates
    ) -> None:
        with pytest.raises(ValidationError):
            Variant(coordinates=brca1_coordinates, variant_type=VariantType.STRUCTURAL)

    def test_cnv_with_symbolic_allele_is_valid(self) -> None:
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="<DUP>")
        variant = Variant(coordinates=coords, variant_type=VariantType.CNV)
        assert variant.variant_type == VariantType.CNV


class TestPopulationFrequency:
    def test_zero_allele_number_forces_frequency_to_none(self) -> None:
        pf = PopulationFrequency(ancestry=AncestryGroup.SAS, allele_count=0, allele_number=0)
        assert pf.allele_frequency is None
        assert pf.has_data is False

    def test_frequency_is_computed_when_not_supplied(self) -> None:
        pf = PopulationFrequency(ancestry=AncestryGroup.NFE, allele_count=1, allele_number=100_000)
        assert pf.allele_frequency == pytest.approx(0.00001)
        assert pf.has_data is True

    def test_allele_count_exceeding_allele_number_is_rejected(self) -> None:
        with pytest.raises(ValidationError):
            PopulationFrequency(ancestry=AncestryGroup.AFR, allele_count=10, allele_number=5)


class TestEvidenceReport:
    def test_clinical_disclaimer_matches_module_constant(
        self,
        brca1_variant: Variant,
        report_metadata: ReportMetadata,
        empty_audit_trail: AuditTrail,
        complete_coverage_summary: CoverageSummary,
        vus_bayesian_score: BayesianScore,
    ) -> None:
        """Guards against the field-level Literal and CLINICAL_DISCLAIMER drifting apart.

        mypy cannot verify `Literal[CLINICAL_DISCLAIMER] == CLINICAL_DISCLAIMER` statically
        (see the NOTE in alleleiq/models/report.py), so this test verifies it at runtime.
        """
        report = EvidenceReport(
            variant=brca1_variant,
            metadata=report_metadata,
            coverage_summary=complete_coverage_summary,
            bayesian_score=vus_bayesian_score,
            audit_trail=empty_audit_trail,
        )
        assert report.clinical_disclaimer == CLINICAL_DISCLAIMER

    def test_clinical_disclaimer_cannot_be_overridden(
        self,
        brca1_variant: Variant,
        report_metadata: ReportMetadata,
        empty_audit_trail: AuditTrail,
        complete_coverage_summary: CoverageSummary,
        vus_bayesian_score: BayesianScore,
    ) -> None:
        with pytest.raises(ValidationError):
            EvidenceReport(
                variant=brca1_variant,
                metadata=report_metadata,
                coverage_summary=complete_coverage_summary,
                bayesian_score=vus_bayesian_score,
                audit_trail=empty_audit_trail,
                clinical_disclaimer="This variant is definitely pathogenic.",
            )

    def test_report_serialises_to_json(
        self,
        brca1_variant: Variant,
        report_metadata: ReportMetadata,
        empty_audit_trail: AuditTrail,
        complete_coverage_summary: CoverageSummary,
        vus_bayesian_score: BayesianScore,
    ) -> None:
        report = EvidenceReport(
            variant=brca1_variant,
            metadata=report_metadata,
            coverage_summary=complete_coverage_summary,
            bayesian_score=vus_bayesian_score,
            audit_trail=empty_audit_trail,
        )
        payload = report.model_dump_json()
        assert "CLINICAL SCIENTIST REVIEW ONLY" in payload
        assert "BRCA1" in payload


class TestBayesianScore:
    def test_preliminary_classification_and_confidence_round_trip(self) -> None:
        score = BayesianScore(
            evidence_points=10,
            preliminary_classification=PreliminaryClassification.PATHOGENIC,
            classification_confidence=ClassificationConfidence.HIGH,
        )
        assert score.preliminary_classification == PreliminaryClassification.PATHOGENIC
        assert score.classification_confidence == ClassificationConfidence.HIGH

    def test_conflict_flags_default_to_empty(self) -> None:
        score = BayesianScore(
            evidence_points=0,
            preliminary_classification=PreliminaryClassification.UNCERTAIN_SIGNIFICANCE,
            classification_confidence=ClassificationConfidence.LOW,
        )
        assert score.conflict_flags == []
