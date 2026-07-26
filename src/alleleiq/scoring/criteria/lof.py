"""PVS1: loss-of-function ACMG/AMP criterion.

Implements a simplified version of the ACGS 2024 PVS1 decision tree.
Gene-disease mechanism strength is resolved from ClinGen's bundled
Gene-Disease Validity and Dosage Sensitivity Map snapshots (see
:mod:`alleleiq.clients.clingen`) as the primary evidence source, tiered
as follows:

- **Definitive/Strong** gene-disease validity *and* a dosage-map
  haploinsufficiency score of 3 -> full-strength PVS1 (Very Strong,
  subject to the existing last/penultimate-exon downgrade to Strong).
- **Moderate** gene-disease validity *and* a haploinsufficiency score of
  2 -> a fixed Strong-strength PVS1 (below full strength, but not
  contingent on exon position — the classification ceiling, not
  transcript-truncation timing, is the limiting factor here).
- Any other combination (Limited/Disputed/Refuted/No Known Disease
  Relationship validity, a haploinsufficiency score below the tier
  thresholds, or the gene simply absent from ClinGen's curation) ->
  Supporting-strength PVS1, flagged ``requires_clinical_input``.

When ClinGen has no curation at all for a gene (as opposed to a curation
that fails the tier thresholds above), PVS1 falls back to the original
v0.1.0-alpha heuristic — a small bundled list of well-established
loss-of-function-mechanism genes, or gnomAD's own pLI/LOEUF constraint
metrics — before dropping to Supporting strength. This keeps a
best-effort answer available for genes ClinGen has not yet curated,
without letting that fallback override an explicit ClinGen classification
that says otherwise.

Three aspects of the full decision tree remain explicitly out of scope
and are documented as known limitations rather than silently
approximated:

1. **Transcript-specific expression regions** (a real PVS1 exclusion —
   a LoF variant in an exon poorly expressed in the disease-relevant
   tissue should not get full-strength PVS1) are not checked, since
   AlleleIQ has no data source for tissue-specific transcript usage.
2. The ACGS decision tree's finer-grained NMD-escape and reading-frame
   sub-rules beyond the last/penultimate-exon check are not implemented.
3. Structural exon-deletion LoF (CNV/STRUCTURAL variant types) is not
   assessed; PVS1 here only evaluates SNV/INDEL consequence terms.
"""

from __future__ import annotations

from dataclasses import dataclass

from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength

from .base import CriterionContext, CriterionEvaluator

_DATA_SOURCE = "Ensembl VEP + ClinGen gene-disease validity/dosage sensitivity"

_FULL_STRENGTH_CLASSIFICATIONS = frozenset({"Definitive", "Strong"})
_FULL_STRENGTH_HI_SCORE = 3
_MODERATE_CLASSIFICATION = "Moderate"
_MODERATE_HI_SCORE = 2
_DISPUTED_OR_REFUTED = frozenset({"Disputed", "Refuted"})

_LOF_CONSEQUENCES = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "start_lost",
    }
)
_NON_CANONICAL_SPLICE_CONSEQUENCES = frozenset(
    {
        "splice_region_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
    }
)

# Minimal bundled fallback of genes with a ClinGen-definitive haploinsufficiency /
# loss-of-function disease mechanism. Illustrative only, pending a live ClinGen
# Gene-Disease Validity API integration — do not treat as clinically exhaustive.
_ESTABLISHED_LOF_MECHANISM_GENES = frozenset(
    {
        "BRCA1",
        "BRCA2",
        "TP53",
        "MLH1",
        "MSH2",
        "MSH6",
        "PMS2",
        "APC",
        "MUTYH",
        "PTEN",
        "STK11",
        "CDH1",
        "SMAD4",
        "BMPR1A",
        "VHL",
        "RB1",
        "NF1",
        "NF2",
        "TSC1",
        "TSC2",
        "MEN1",
        "RET",
        "SDHB",
        "SDHD",
        "LDLR",
    }
)


@dataclass(frozen=True)
class _MechanismTier:
    """Resolved haploinsufficiency-mechanism strength for one gene.

    ``strength`` is the ceiling PVS1 strength this tier permits — for
    :attr:`~alleleiq.models.evidence.EvidenceStrength.VERY_STRONG` it may
    still be downgraded to Strong by the last/penultimate-exon check;
    for a ClinGen Moderate-tier result (``fixed_strength=True``) the
    caller applies ``strength`` (Strong) directly, without the exon check.
    """

    strength: EvidenceStrength
    reason: str
    fixed_strength: bool = False
    note: str | None = None


def _parse_exon(exon_field: str | None) -> tuple[int, int] | None:
    """Parse VEP's ``exon`` field, e.g. ``'18/23'`` -> ``(18, 23)``."""
    if not exon_field or "/" not in exon_field:
        return None
    numerator, _, denominator = exon_field.partition("/")
    try:
        return int(numerator), int(denominator)
    except ValueError:
        return None


class PVS1Criterion(CriterionEvaluator):
    """PVS1: null variant in a gene where loss-of-function is a known disease mechanism."""

    code = ACMGCode.PVS1
    default_strength = EvidenceStrength.SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        annotation = context.vep_annotation
        if annotation is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No VEP annotation available to assess loss-of-function consequence."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
            )

        consequence = annotation.most_severe_consequence

        if consequence in _NON_CANONICAL_SPLICE_CONSEQUENCES:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"Consequence '{consequence}' is a non-canonical splice-region variant "
                    "(not the canonical ±1/±2 splice donor/acceptor dinucleotide); PVS1 does "
                    "not apply."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
                notes="Assess splicing impact via the SpliceAI-based PP3 criterion instead.",
            )

        if consequence not in _LOF_CONSEQUENCES:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"Consequence '{consequence}' is not a loss-of-function type; "
                    "PVS1 not applicable."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        gene_symbol = context.variant.gene_symbol or annotation.gene_symbol
        tier = self._resolve_mechanism_tier(gene_symbol, context)

        if tier.strength is EvidenceStrength.SUPPORTING:
            return CriterionEvaluation(
                code=self.code,
                strength=EvidenceStrength.SUPPORTING,
                met=True,
                evidence_basis=(
                    f"Consequence '{consequence}' is loss-of-function, but "
                    f"{gene_symbol or 'this gene'}'s haploinsufficiency mechanism is not "
                    f"established at full or moderate tier ({tier.reason}); applying PVS1 "
                    "at Supporting strength only."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes=tier.note,
            )

        if tier.strength is EvidenceStrength.STRONG and tier.fixed_strength:
            # ClinGen Moderate validity + HI score 2: a fixed Strong-strength
            # outcome, not subject to the exon-position downgrade below — the
            # classification ceiling is what limits this, not NMD timing.
            return CriterionEvaluation(
                code=self.code,
                strength=EvidenceStrength.STRONG,
                met=True,
                evidence_basis=(
                    f"Consequence '{consequence}' in {gene_symbol}: ClinGen classifies this "
                    f"gene-disease relationship as Moderate with a haploinsufficiency score "
                    "of 2; PVS1 applied at Strong strength."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        exon_position = _parse_exon(annotation.exon)
        if exon_position is not None:
            exon_number, total_exons = exon_position
            if exon_number >= total_exons - 1:
                return CriterionEvaluation(
                    code=self.code,
                    strength=EvidenceStrength.STRONG,
                    met=True,
                    evidence_basis=(
                        f"Consequence '{consequence}' at exon {exon_number}/{total_exons} "
                        "(last or penultimate exon) in a gene with an established "
                        f"loss-of-function mechanism ({tier.reason}); nonsense-mediated decay "
                        "may be incomplete this close to the transcript's end, so PVS1 is "
                        "applied at Strong rather than Very Strong strength."
                    ),
                    data_source=_DATA_SOURCE,
                )

        return CriterionEvaluation(
            code=self.code,
            strength=EvidenceStrength.VERY_STRONG,
            met=True,
            evidence_basis=(
                f"Consequence '{consequence}' in {gene_symbol}, a gene with an established "
                f"loss-of-function disease mechanism ({tier.reason}), outside the "
                "last/penultimate exon; full-strength PVS1 applies."
            ),
            data_source=_DATA_SOURCE,
            notes=(
                "Transcript-specific expression regions are not checked (no data source "
                "integrated); confirm the affected exon is broadly expressed in "
                "disease-relevant tissue."
            ),
        )

    @classmethod
    def _resolve_mechanism_tier(
        cls, gene_symbol: str | None, context: CriterionContext
    ) -> _MechanismTier:
        classification = context.clingen_classification
        if classification is not None:
            hi_score = context.clingen_haploinsufficiency_score or 0
            if (
                classification in _FULL_STRENGTH_CLASSIFICATIONS
                and hi_score == _FULL_STRENGTH_HI_SCORE
            ):
                return _MechanismTier(
                    strength=EvidenceStrength.VERY_STRONG,
                    reason=f"ClinGen classification '{classification}', haploinsufficiency score 3",
                )
            if classification == _MODERATE_CLASSIFICATION and hi_score == _MODERATE_HI_SCORE:
                return _MechanismTier(
                    strength=EvidenceStrength.STRONG,
                    reason="ClinGen classification 'Moderate', haploinsufficiency score 2",
                    fixed_strength=True,
                )
            note = (
                "Confirm gene-disease mechanism via ClinGen Gene-Disease Validity curation "
                "before relying on this criterion."
            )
            if classification in _DISPUTED_OR_REFUTED:
                note = (
                    f"ClinGen classifies this gene-disease relationship as '{classification}' — "
                    "treat any PVS1 application with particular caution and confirm manually."
                )
            return _MechanismTier(
                strength=EvidenceStrength.SUPPORTING,
                reason=(
                    f"ClinGen classification '{classification}', haploinsufficiency score "
                    f"{hi_score}, below the Definitive/Strong+3 or Moderate+2 tier thresholds"
                ),
                note=note,
            )

        # No ClinGen curation was consulted at all (gene absent from ClinGen's
        # curation, or the caller did not supply ClinGen fields) — fall back
        # to the v0.1.0-alpha heuristic rather than defaulting straight to
        # Supporting, since a well-known LoF-mechanism gene ClinGen simply
        # has not yet curated should not lose all automated signal.
        if cls._has_established_lof_mechanism(gene_symbol, context):
            return _MechanismTier(
                strength=EvidenceStrength.VERY_STRONG,
                reason="bundled fallback gene list / gnomAD constraint (no ClinGen curation)",
            )
        return _MechanismTier(
            strength=EvidenceStrength.SUPPORTING,
            reason="not established from bundled ClinGen data or gnomAD constraint metrics",
            note=(
                "Confirm gene-disease mechanism via ClinGen Gene-Disease Validity "
                "curation before relying on this criterion."
            ),
        )

    @staticmethod
    def _has_established_lof_mechanism(gene_symbol: str | None, context: CriterionContext) -> bool:
        if gene_symbol and gene_symbol.upper() in _ESTABLISHED_LOF_MECHANISM_GENES:
            return True
        constraint = context.gene_constraint
        if constraint is None:
            return False
        settings = context.settings
        pli_supports = (
            constraint.pli is not None
            and constraint.pli > settings.pli_haploinsufficiency_threshold
        )
        loeuf_supports = (
            constraint.loeuf is not None
            and constraint.loeuf < settings.loeuf_haploinsufficiency_threshold
        )
        return pli_supports or loeuf_supports
