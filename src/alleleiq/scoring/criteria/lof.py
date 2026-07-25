"""PVS1: loss-of-function ACMG/AMP criterion.

Implements a simplified version of the ACGS 2024 PVS1 decision tree.
Two aspects of the full decision tree are explicitly out of scope for
v0.1.0-alpha and documented as known limitations rather than silently
approximated:

1. **Gene-disease validity** comes from a small bundled fallback list of
   well-established loss-of-function-mechanism genes, not a live query
   against the ClinGen Gene-Disease Validity API (no such integration
   exists yet). This list is illustrative, not clinically exhaustive.
2. **Transcript-specific expression regions** (a real PVS1 exclusion —
   a LoF variant in an exon poorly expressed in the disease-relevant
   tissue should not get full-strength PVS1) are not checked, since
   AlleleIQ has no data source for tissue-specific transcript usage.
3. The spec's "PVS1_moderate: uncertainty about LoF mechanism" and
   "PVS1_supporting: LoF in gene where mechanism uncertain" tiers
   describe the same underlying situation without a clear
   differentiator; this implementation collapses both into a single
   Supporting-strength outcome, flagged ``requires_clinical_input``.
4. Structural exon-deletion LoF (CNV/STRUCTURAL variant types) is not
   assessed; PVS1 here only evaluates SNV/INDEL consequence terms.
"""

from __future__ import annotations

from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength

from .base import CriterionContext, CriterionEvaluator

_DATA_SOURCE = "Ensembl VEP + gnomAD constraint + bundled ClinGen fallback"

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
        established = self._has_established_lof_mechanism(gene_symbol, context)

        if not established:
            return CriterionEvaluation(
                code=self.code,
                strength=EvidenceStrength.SUPPORTING,
                met=True,
                evidence_basis=(
                    f"Consequence '{consequence}' is loss-of-function, but "
                    f"{gene_symbol or 'this gene'}'s haploinsufficiency mechanism is not "
                    "established from bundled ClinGen data or gnomAD constraint metrics; "
                    "applying PVS1 at Supporting strength only."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes=(
                    "Confirm gene-disease mechanism via ClinGen Gene-Disease Validity "
                    "curation before relying on this criterion."
                ),
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
                        "loss-of-function mechanism; nonsense-mediated decay may be "
                        "incomplete this close to the transcript's end, so PVS1 is applied "
                        "at Strong rather than Very Strong strength."
                    ),
                    data_source=_DATA_SOURCE,
                )

        return CriterionEvaluation(
            code=self.code,
            strength=EvidenceStrength.VERY_STRONG,
            met=True,
            evidence_basis=(
                f"Consequence '{consequence}' in {gene_symbol}, a gene with an established "
                "loss-of-function disease mechanism, outside the last/penultimate exon; "
                "full-strength PVS1 applies."
            ),
            data_source=_DATA_SOURCE,
            notes=(
                "Transcript-specific expression regions are not checked (no data source "
                "integrated); confirm the affected exon is broadly expressed in "
                "disease-relevant tissue."
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
