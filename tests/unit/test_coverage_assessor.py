"""Unit tests for alleleiq.scoring.coverage.assess_coverage."""

from __future__ import annotations

from alleleiq.models.population import (
    AncestryGroup,
    CoverageFlagLevel,
    OverallCoverageAssessment,
    PopulationFrequency,
)
from alleleiq.scoring.coverage import assess_coverage


def _well_covered() -> dict[AncestryGroup, PopulationFrequency]:
    return {
        a: PopulationFrequency(ancestry=a, allele_count=5, allele_number=50_000)
        for a in AncestryGroup
    }


class TestCompleteCoverage:
    def test_all_populations_above_threshold_yields_complete(self) -> None:
        gaps, summary = assess_coverage(_well_covered())

        assert gaps == {}
        assert summary.overall_coverage_assessment == OverallCoverageAssessment.COMPLETE
        assert len(summary.populations_with_sufficient_data) == len(AncestryGroup)


class TestZeroAlleleNumber:
    def test_explicit_zero_allele_number_yields_no_data_gap(self) -> None:
        freqs = _well_covered()
        freqs[AncestryGroup.SAS] = PopulationFrequency(
            ancestry=AncestryGroup.SAS, allele_count=0, allele_number=0
        )

        gaps, summary = assess_coverage(freqs)

        assert gaps[AncestryGroup.SAS].flag_level == CoverageFlagLevel.NO_DATA
        assert AncestryGroup.SAS in summary.populations_with_no_data

    def test_ancestry_missing_from_mapping_entirely_is_also_no_data(self) -> None:
        freqs = _well_covered()
        del freqs[AncestryGroup.OTH]

        gaps, summary = assess_coverage(freqs)

        assert gaps[AncestryGroup.OTH].flag_level == CoverageFlagLevel.NO_DATA
        assert AncestryGroup.OTH in summary.populations_with_no_data


class TestPartialCoverage:
    def test_mix_of_sufficient_and_insufficient_yields_partial(self) -> None:
        freqs = _well_covered()
        freqs[AncestryGroup.SAS] = PopulationFrequency(
            ancestry=AncestryGroup.SAS, allele_count=0, allele_number=0
        )

        _, summary = assess_coverage(freqs)

        assert summary.overall_coverage_assessment == OverallCoverageAssessment.PARTIAL
        assert AncestryGroup.NFE in summary.populations_with_sufficient_data

    def test_warning_vs_critical_flag_thresholds(self) -> None:
        freqs = _well_covered()
        # AFR threshold is 10_000; 9_000 is a WARNING (above 50% of threshold).
        freqs[AncestryGroup.AFR] = PopulationFrequency(
            ancestry=AncestryGroup.AFR, allele_count=1, allele_number=9_000
        )
        # AMR threshold is 5_000; 1_000 is CRITICAL (below 50% of threshold).
        freqs[AncestryGroup.AMR] = PopulationFrequency(
            ancestry=AncestryGroup.AMR, allele_count=1, allele_number=1_000
        )

        gaps, _ = assess_coverage(freqs)

        assert gaps[AncestryGroup.AFR].flag_level == CoverageFlagLevel.WARNING
        assert gaps[AncestryGroup.AMR].flag_level == CoverageFlagLevel.CRITICAL

    def test_guidance_message_names_a_supplementary_resource(self) -> None:
        freqs = _well_covered()
        freqs[AncestryGroup.AFR] = PopulationFrequency(
            ancestry=AncestryGroup.AFR, allele_count=0, allele_number=0
        )

        gaps, _ = assess_coverage(freqs)

        assert "H3Africa" in gaps[AncestryGroup.AFR].clinical_guidance_message


class TestSeverelyLimitedCoverage:
    def test_all_populations_below_threshold_yields_severely_limited(self) -> None:
        gaps, summary = assess_coverage({})

        assert len(gaps) == len(AncestryGroup)
        assert summary.overall_coverage_assessment == OverallCoverageAssessment.SEVERELY_LIMITED
        assert summary.populations_with_sufficient_data == []
