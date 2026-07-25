"""Population-frequency ACMG/AMP criteria: BA1, BS1, PM2.

This module is AlleleIQ's clinical-equity core. Every criterion here
evaluates each ancestry group *independently* and only draws on
populations with adequate gnomAD v4 coverage (i.e. absent from
:attr:`~alleleiq.scoring.criteria.base.CriterionContext.coverage_gaps`).
A population's frequency is never averaged with another population's
before a criterion is applied — doing so is exactly the step that hides
ancestry-specific reference-data gaps, which is the problem AlleleIQ
exists to surface rather than paper over.
"""

from __future__ import annotations

from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength
from alleleiq.models.population import AncestryGroup, PopulationFrequency

from .base import CriterionContext, CriterionEvaluator

_DATA_SOURCE = "gnomAD v4"


def _reliable_frequencies(
    context: CriterionContext,
) -> dict[AncestryGroup, PopulationFrequency]:
    """Populations with adequate coverage and an observed allele frequency.

    A population counts as "reliable" only if it has no coverage gap
    flagged *and* was actually returned by gnomAD (defensive: if the
    coverage assessor has not run, an ancestry missing from both
    ``population_frequencies`` and ``coverage_gaps`` is treated as
    unreliable rather than silently ignored).
    """
    return {
        ancestry: pf
        for ancestry, pf in context.population_frequencies.items()
        if ancestry not in context.coverage_gaps and pf.has_data
    }


def _excluded_ancestries(context: CriterionContext) -> list[AncestryGroup]:
    """Every ancestry group excluded from consideration due to a coverage gap."""
    return [
        ancestry
        for ancestry in AncestryGroup
        if ancestry in context.coverage_gaps or ancestry not in context.population_frequencies
    ]


def _coverage_note(excluded: list[AncestryGroup]) -> str | None:
    if not excluded:
        return None
    names = ", ".join(sorted(a.display_name for a in excluded))
    return f"Excluded from consideration due to insufficient coverage: {names}."


class BA1Criterion(CriterionEvaluator):
    """BA1 (Benign Standalone): allele frequency too high to be pathogenic.

    ACGS 2024: applied per-population against ``ba1_frequency_threshold``
    (default 5%). A population whose only high-frequency observation
    comes from an under-covered cohort does not trigger BA1 — a small
    denominator makes the frequency estimate unreliable, and treating it
    as reliable would be exactly the kind of silent, ancestry-blind
    shortcut this criterion must not take.
    """

    code = ACMGCode.BA1
    default_strength = EvidenceStrength.BENIGN_STANDALONE

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        threshold = context.settings.ba1_frequency_threshold
        reliable = _reliable_frequencies(context)
        excluded = _excluded_ancestries(context)

        exceeding = [
            (ancestry, pf.allele_frequency)
            for ancestry, pf in reliable.items()
            if pf.allele_frequency is not None and pf.allele_frequency > threshold
        ]

        if exceeding:
            ancestry, frequency = max(exceeding, key=lambda pair: pair[1] or 0.0)
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=True,
                evidence_basis=(
                    f"Allele frequency {frequency:.4f} in {ancestry.display_name} exceeds "
                    f"the BA1 threshold of {threshold:.2f}."
                ),
                data_source=_DATA_SOURCE,
                notes=_coverage_note(excluded),
            )

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=False,
            evidence_basis=(
                "No adequately covered population showed an allele frequency above "
                f"the BA1 threshold of {threshold:.2f}."
            ),
            data_source=_DATA_SOURCE,
            requires_clinical_input=not reliable,
            notes=_coverage_note(excluded),
        )


class BS1Criterion(CriterionEvaluator):
    """BS1 (Benign Strong): allele frequency higher than expected for the disorder.

    Requires disease-specific context (inheritance mode) to select the
    appropriate threshold; without it, this criterion cannot be assessed
    and is marked ``requires_clinical_input``.
    """

    code = ACMGCode.BS1
    default_strength = EvidenceStrength.BENIGN_STRONG

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        if context.inheritance_mode is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "BS1 requires a disease inheritance-mode/prevalence context to select "
                    "the appropriate frequency threshold, which was not provided."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes=(
                    "Supply inheritance_mode ('dominant' or 'recessive') to enable "
                    "automated BS1 assessment."
                ),
            )

        threshold = (
            context.settings.bs1_threshold_dominant
            if context.inheritance_mode == "dominant"
            else context.settings.bs1_threshold_recessive
        )
        reliable = _reliable_frequencies(context)
        excluded = _excluded_ancestries(context)

        exceeding = [
            (ancestry, pf.allele_frequency)
            for ancestry, pf in reliable.items()
            if pf.allele_frequency is not None and pf.allele_frequency > threshold
        ]

        if exceeding:
            ancestry, frequency = max(exceeding, key=lambda pair: pair[1] or 0.0)
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=True,
                evidence_basis=(
                    f"Allele frequency {frequency:.4f} in {ancestry.display_name} exceeds "
                    f"the BS1 {context.inheritance_mode} threshold of {threshold:.2f}."
                ),
                data_source=_DATA_SOURCE,
                notes=_coverage_note(excluded),
            )

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=False,
            evidence_basis=(
                f"No adequately covered population exceeded the BS1 {context.inheritance_mode} "
                f"threshold of {threshold:.2f}."
            ),
            data_source=_DATA_SOURCE,
            requires_clinical_input=not reliable,
            notes=_coverage_note(excluded),
        )


class PM2Criterion(CriterionEvaluator):
    """PM2_Supporting (ACGS 2024): absent from adequately covered gnomAD v4 populations.

    ACGS 2024 downgrades PM2 from Moderate to Supporting strength by
    default (absent strong gene-disease evidence justifying Moderate,
    which AlleleIQ does not have access to). PM2 is only assessed as
    "met" against populations with adequate coverage; a variant absent
    only because a population was never sampled is a coverage gap, not
    evidence of rarity.

    Known simplification: this evaluator assesses true absence
    (allele count zero in every adequately covered population). ACGS
    2024's "or at extremely low frequency" nuance is not automated here
    (no clinically validated universal rarity threshold was specified),
    so a variant observed at very low but non-zero frequency is reported
    as PM2 not met with ``requires_clinical_input=True``, directing the
    clinical scientist to assess whether the observed frequency is low
    enough to still support PM2_Supporting.
    """

    code = ACMGCode.PM2
    default_strength = EvidenceStrength.SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        reliable = _reliable_frequencies(context)
        excluded = _excluded_ancestries(context)

        if not reliable:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No population had sufficient gnomAD v4 coverage to assess absence or rarity."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes=_coverage_note(excluded),
            )

        observed = {ancestry: pf for ancestry, pf in reliable.items() if pf.allele_count > 0}
        if observed:
            names = ", ".join(sorted(a.display_name for a in observed))
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"Variant observed in adequately covered population(s): {names}; "
                    "true absence does not hold."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes=(
                    "Automated evaluation only assesses true absence, not the 'extremely low "
                    "frequency' variant of PM2_Supporting; clinical review of observed "
                    "frequency is required. " + (_coverage_note(excluded) or "")
                ).strip(),
            )

        names = ", ".join(sorted(a.display_name for a in reliable))
        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=True,
            evidence_basis=(
                f"Variant absent from all adequately covered gnomAD v4 populations: {names}."
            ),
            data_source=_DATA_SOURCE,
            notes=_coverage_note(excluded),
        )
