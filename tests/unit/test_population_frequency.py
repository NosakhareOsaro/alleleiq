"""Unit tests for alleleiq.scoring.criteria.population_frequency (BA1, BS1, PM2)."""

from __future__ import annotations

from alleleiq.models.evidence import ACMGCode
from alleleiq.models.population import (
    AncestryGroup,
    CoverageFlagLevel,
    CoverageGap,
    PopulationFrequency,
)
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.scoring.criteria.base import CriterionContext
from alleleiq.scoring.criteria.population_frequency import BA1Criterion, BS1Criterion, PM2Criterion


def _variant() -> Variant:
    coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")
    return Variant(coordinates=coords, variant_type=VariantType.SNV, gene_symbol="BRCA1")


def _pf(ancestry: AncestryGroup, ac: int, an: int) -> PopulationFrequency:
    return PopulationFrequency(ancestry=ancestry, allele_count=ac, allele_number=an)


def _gap(ancestry: AncestryGroup, an: int = 0, threshold: int = 10_000) -> CoverageGap:
    flag = CoverageFlagLevel.NO_DATA if an == 0 else CoverageFlagLevel.WARNING
    return CoverageGap(
        ancestry=ancestry,
        allele_number=an,
        threshold=threshold,
        flag_level=flag,
        clinical_guidance_message="see supplementary databases",
    )


class TestBA1:
    def test_high_frequency_in_all_populations_triggers_ba1(self) -> None:
        freqs = {a: _pf(a, 5000, 20_000) for a in AncestryGroup}
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs)

        result = BA1Criterion().evaluate(ctx)

        assert result.code == ACMGCode.BA1
        assert result.met is True
        assert result.requires_clinical_input is False

    def test_high_frequency_in_one_population_only_triggers_ba1_with_note(self) -> None:
        freqs = {
            AncestryGroup.NFE: _pf(AncestryGroup.NFE, 8000, 100_000),
            AncestryGroup.AFR: _pf(AncestryGroup.AFR, 5, 40_000),
        }
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs)

        result = BA1Criterion().evaluate(ctx)

        assert result.met is True
        assert "Non-Finnish European" in result.evidence_basis

    def test_high_frequency_only_in_coverage_gapped_population_does_not_trigger_ba1(self) -> None:
        freqs = {
            AncestryGroup.SAS: _pf(AncestryGroup.SAS, 400, 2_000),  # 20% AF but AN below threshold
            AncestryGroup.NFE: _pf(AncestryGroup.NFE, 5, 100_000),
        }
        gaps = {AncestryGroup.SAS: _gap(AncestryGroup.SAS, an=2_000)}
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs, coverage_gaps=gaps)

        result = BA1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False

    def test_all_populations_gapped_requires_clinical_input(self) -> None:
        freqs = {AncestryGroup.SAS: _pf(AncestryGroup.SAS, 0, 0)}
        gaps = {AncestryGroup.SAS: _gap(AncestryGroup.SAS, an=0)}
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs, coverage_gaps=gaps)

        result = BA1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True


class TestPM2:
    def test_absent_in_all_well_covered_populations_triggers_pm2_supporting(self) -> None:
        freqs = {
            AncestryGroup.NFE: _pf(AncestryGroup.NFE, 0, 68_034),
            AncestryGroup.AFR: _pf(AncestryGroup.AFR, 0, 41_438),
        }
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs)

        result = PM2Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength.value == "SUPPORTING"
        assert result.requires_clinical_input is False

    def test_absent_in_some_coverage_gap_in_others_still_applies_with_warning(self) -> None:
        freqs = {
            AncestryGroup.NFE: _pf(AncestryGroup.NFE, 0, 68_034),
            AncestryGroup.SAS: _pf(AncestryGroup.SAS, 0, 0),
        }
        gaps = {AncestryGroup.SAS: _gap(AncestryGroup.SAS, an=0)}
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs, coverage_gaps=gaps)

        result = PM2Criterion().evaluate(ctx)

        assert result.met is True
        assert result.notes is not None
        assert "South Asian" in result.notes

    def test_low_but_nonzero_frequency_does_not_trigger_pm2(self) -> None:
        freqs = {AncestryGroup.NFE: _pf(AncestryGroup.NFE, 1, 100_000)}
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs)

        result = PM2Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True

    def test_gnomad_data_entirely_absent_requires_clinical_input(self) -> None:
        ctx = CriterionContext(variant=_variant(), population_frequencies={}, coverage_gaps={})

        result = PM2Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True


class TestBS1:
    def test_no_inheritance_mode_requires_clinical_input(self) -> None:
        freqs = {AncestryGroup.NFE: _pf(AncestryGroup.NFE, 5000, 100_000)}
        ctx = CriterionContext(variant=_variant(), population_frequencies=freqs)

        result = BS1Criterion().evaluate(ctx)

        assert result.requires_clinical_input is True

    def test_dominant_threshold_exceeded(self) -> None:
        freqs = {
            AncestryGroup.NFE: _pf(AncestryGroup.NFE, 2000, 100_000)
        }  # AF=0.02 > 0.01 dominant
        ctx = CriterionContext(
            variant=_variant(), population_frequencies=freqs, inheritance_mode="dominant"
        )

        result = BS1Criterion().evaluate(ctx)

        assert result.met is True

    def test_recessive_threshold_not_exceeded_by_dominant_level_frequency(self) -> None:
        freqs = {
            AncestryGroup.NFE: _pf(AncestryGroup.NFE, 2000, 100_000)
        }  # AF=0.02 < 0.03 recessive
        ctx = CriterionContext(
            variant=_variant(), population_frequencies=freqs, inheritance_mode="recessive"
        )

        result = BS1Criterion().evaluate(ctx)

        assert result.met is False
