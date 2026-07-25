"""Population coverage gap assessment: AlleleIQ's clinical-equity core.

Coverage gaps are never suppressed, hidden, or averaged away. Every
ancestry group below its configured minimum allele-number threshold is
surfaced as an explicit :class:`~alleleiq.models.population.CoverageGap`
with a severity flag and a specific, actionable clinical guidance
message — not folded into a single "data quality: mixed" note that a
reviewer could skim past.
"""

from __future__ import annotations

from alleleiq.config import Settings, get_settings
from alleleiq.models.population import (
    AncestryGroup,
    CoverageFlagLevel,
    CoverageGap,
    OverallCoverageAssessment,
    PopulationFrequency,
)
from alleleiq.models.report import CoverageSummary

# If this fraction or more of all ancestry groups have a coverage gap
# (WARNING, CRITICAL, or NO_DATA), the overall assessment is downgraded
# from PARTIAL to SEVERELY_LIMITED rather than treating widespread gaps
# as merely "partial" coverage.
_SEVERELY_LIMITED_FRACTION = 0.75

_SUPPLEMENTARY_RESOURCE_HINTS: dict[AncestryGroup, str] = {
    AncestryGroup.AFR: (
        "Consider consulting H3Africa or other African-ancestry-specific variant "
        "resources, which remain substantially under-represented in gnomAD."
    ),
    AncestryGroup.AMR: (
        "Consider consulting Latino/Admixed American-specific cohort studies "
        "(e.g. regional biobanks), which are under-represented in gnomAD relative "
        "to their global population share."
    ),
    AncestryGroup.SAS: (
        "Consider consulting GenomeAsia 100K or other South Asian-specific " "reference resources."
    ),
    AncestryGroup.EAS: (
        "Consider consulting ChinaMAP, ToMMo, or other East Asian-specific "
        "reference resources for corroborating frequency data."
    ),
    AncestryGroup.ASJ: (
        "Ashkenazi Jewish reference data is comparatively strong in gnomAD, but "
        "consider founder-population-specific registries for rare disease genes."
    ),
    AncestryGroup.FIN: (
        "Finnish reference data is comparatively strong in gnomAD given the "
        "well-studied Finnish founder population; corroborate with FinnGen if "
        "available for this gene."
    ),
    AncestryGroup.NFE: (
        "Non-Finnish European reference data is typically the most complete in "
        "gnomAD; a gap here may indicate an unusually rare or novel variant."
    ),
    AncestryGroup.OTH: (
        "The 'Other/Unassigned' group is inherently heterogeneous; consider "
        "whether a more specific ancestry assignment and targeted database "
        "consultation is possible for this patient."
    ),
}


def _guidance_message(ancestry: AncestryGroup, flag_level: CoverageFlagLevel) -> str:
    severity_text = {
        CoverageFlagLevel.NO_DATA: (
            f"No {ancestry.display_name} allele data is available in gnomAD v4 for this "
            "variant. Absence-based criteria (e.g. PM2) cannot be applied for this population."
        ),
        CoverageFlagLevel.CRITICAL: (
            f"{ancestry.display_name} allele number is far below the minimum threshold for "
            "reliable frequency estimation; any observed frequency in this population should "
            "be treated with substantial caution."
        ),
        CoverageFlagLevel.WARNING: (
            f"{ancestry.display_name} allele number is below the minimum threshold for fully "
            "reliable frequency estimation; treat any observed frequency in this population "
            "with some caution."
        ),
    }[flag_level]
    return f"{severity_text} {_SUPPLEMENTARY_RESOURCE_HINTS[ancestry]}"


def assess_coverage(
    population_frequencies: dict[AncestryGroup, PopulationFrequency],
    settings: Settings | None = None,
) -> tuple[dict[AncestryGroup, CoverageGap], CoverageSummary]:
    """Assess per-ancestry gnomAD v4 coverage adequacy for one variant.

    Args:
        population_frequencies: Per-ancestry frequency data returned by
            :class:`~alleleiq.clients.gnomad.GnomADClient`. An ancestry
            group absent from this mapping is treated identically to
            one present with ``allele_number=0``.
        settings: Configuration providing per-ancestry thresholds;
            defaults to :func:`alleleiq.config.get_settings`.

    Returns:
        A tuple of (per-ancestry coverage gaps, report-level summary).
        Every :class:`~alleleiq.models.population.AncestryGroup` is
        accounted for in exactly one of the summary's three population
        lists.
    """
    settings = settings or get_settings()
    thresholds = settings.coverage_thresholds_by_population_id
    critical_fraction = settings.coverage_critical_fraction

    gaps: dict[AncestryGroup, CoverageGap] = {}
    sufficient: list[AncestryGroup] = []
    insufficient: list[AncestryGroup] = []
    no_data: list[AncestryGroup] = []

    for ancestry in AncestryGroup:
        threshold = thresholds[ancestry.value]
        frequency = population_frequencies.get(ancestry)
        allele_number = frequency.allele_number if frequency is not None else 0

        if allele_number == 0:
            flag_level = CoverageFlagLevel.NO_DATA
            no_data.append(ancestry)
        elif allele_number < threshold:
            flag_level = (
                CoverageFlagLevel.CRITICAL
                if allele_number < threshold * critical_fraction
                else CoverageFlagLevel.WARNING
            )
            insufficient.append(ancestry)
        else:
            sufficient.append(ancestry)
            continue

        gaps[ancestry] = CoverageGap(
            ancestry=ancestry,
            allele_number=allele_number,
            threshold=threshold,
            flag_level=flag_level,
            clinical_guidance_message=_guidance_message(ancestry, flag_level),
        )

    overall = _overall_assessment(sufficient, insufficient, no_data)
    advisory = _clinical_advisory(sufficient, insufficient, no_data, overall)

    summary = CoverageSummary(
        populations_with_sufficient_data=sufficient,
        populations_with_insufficient_data=insufficient,
        populations_with_no_data=no_data,
        overall_coverage_assessment=overall,
        clinical_advisory=advisory,
    )
    return gaps, summary


def _overall_assessment(
    sufficient: list[AncestryGroup],
    insufficient: list[AncestryGroup],
    no_data: list[AncestryGroup],
) -> OverallCoverageAssessment:
    total = len(AncestryGroup)
    problematic = len(insufficient) + len(no_data)
    if problematic == 0:
        return OverallCoverageAssessment.COMPLETE
    if problematic / total >= _SEVERELY_LIMITED_FRACTION:
        return OverallCoverageAssessment.SEVERELY_LIMITED
    return OverallCoverageAssessment.PARTIAL


def _clinical_advisory(
    sufficient: list[AncestryGroup],
    insufficient: list[AncestryGroup],
    no_data: list[AncestryGroup],
    overall: OverallCoverageAssessment,
) -> str:
    if overall == OverallCoverageAssessment.COMPLETE:
        return (
            "All ancestry groups met the minimum gnomAD v4 allele number threshold. "
            "Population-frequency-derived criteria (BA1, BS1, PM2) can be applied with "
            "confidence across all groups."
        )

    gapped = insufficient + no_data
    gapped_names = ", ".join(sorted(a.display_name for a in gapped))

    if overall == OverallCoverageAssessment.SEVERELY_LIMITED:
        return (
            f"Reference data coverage is severely limited: {len(gapped)} of "
            f"{len(AncestryGroup)} ancestry groups ({gapped_names}) fall below the "
            "minimum allele number threshold. Population-frequency-derived criteria "
            "should be applied with significant caution and supplemented with "
            "population-specific resources and clinical judgement wherever possible."
        )

    sufficient_names = (
        ", ".join(sorted(a.display_name for a in sufficient)) if sufficient else "none"
    )
    return (
        f"Reference data coverage is partial: {gapped_names} fall below the minimum "
        f"allele number threshold, while {sufficient_names} are adequately covered. "
        "Population-frequency-derived criteria should be applied only against the "
        "adequately covered populations listed above."
    )
