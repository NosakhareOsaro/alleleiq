"""Unit tests for alleleiq.scoring.criteria.lof.PVS1Criterion, covering the
full simplified ACGS 2024 PVS1 decision tree.
"""

from __future__ import annotations

from alleleiq.clients.gnomad import GeneConstraint
from alleleiq.clients.vep import VEPAnnotation
from alleleiq.models.evidence import EvidenceStrength
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.scoring.criteria.base import CriterionContext
from alleleiq.scoring.criteria.lof import PVS1Criterion

_COORDS = GenomicCoordinates(chrom="17", pos=43093157, ref="C", alt="A")


def _variant(gene_symbol: str = "BRCA1") -> Variant:
    return Variant(coordinates=_COORDS, variant_type=VariantType.SNV, gene_symbol=gene_symbol)


def _annotation(**overrides: object) -> VEPAnnotation:
    defaults: dict[str, object] = {
        "variant_id": _COORDS.variant_id,
        "most_severe_consequence": "stop_gained",
        "gene_symbol": "BRCA1",
        "exon": "5/23",
    }
    defaults.update(overrides)
    return VEPAnnotation(**defaults)  # type: ignore[arg-type]


class TestEstablishedGeneCanonicalPosition:
    def test_stop_gained_in_established_gene_not_last_exon_triggers_full_pvs1(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(exon="5/23"))

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.VERY_STRONG
        assert result.requires_clinical_input is False


class TestLastExon:
    def test_stop_gained_in_last_exon_downgrades_to_strong(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(exon="23/23"))

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.STRONG

    def test_stop_gained_in_penultimate_exon_downgrades_to_strong(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(exon="22/23"))

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.STRONG


class TestUnknownHaploinsufficiency:
    def test_frameshift_in_unestablished_gene_downgrades_to_supporting_and_flags(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="UNKNOWNGENE"),
            vep_annotation=_annotation(
                most_severe_consequence="frameshift_variant",
                gene_symbol="UNKNOWNGENE",
                exon="5/10",
            ),
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.requires_clinical_input is True

    def test_gnomad_constraint_can_establish_mechanism_without_bundled_list(self) -> None:
        constraint = GeneConstraint(gene_symbol="OTHERGENE", pli=0.98, loeuf=0.20)
        ctx = CriterionContext(
            variant=_variant(gene_symbol="OTHERGENE"),
            vep_annotation=_annotation(gene_symbol="OTHERGENE", exon="3/10"),
            gene_constraint=constraint,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.VERY_STRONG
        assert result.requires_clinical_input is False

    def test_low_pli_and_high_loeuf_do_not_establish_mechanism(self) -> None:
        constraint = GeneConstraint(gene_symbol="OTHERGENE", pli=0.01, loeuf=1.5)
        ctx = CriterionContext(
            variant=_variant(gene_symbol="OTHERGENE"),
            vep_annotation=_annotation(gene_symbol="OTHERGENE"),
            gene_constraint=constraint,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.requires_clinical_input is True


class TestNotApplicable:
    def test_missense_variant_is_not_applicable(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(most_severe_consequence="missense_variant"),
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False

    def test_no_vep_annotation_requires_clinical_input(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=None)

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True


class TestClinGenMechanismTiers:
    """ClinGen classification + haploinsufficiency score determine PVS1 strength.

    See alleleiq.scoring.criteria.lof's module docstring for the tier
    thresholds: Definitive/Strong + HI 3 -> full strength (subject to the
    exon downgrade); Moderate + HI 2 -> fixed Strong; anything else ->
    Supporting + requires_clinical_input.
    """

    def test_definitive_classification_with_hi_3_grants_full_strength(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Definitive",
            clingen_haploinsufficiency_score=3,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.VERY_STRONG
        assert result.requires_clinical_input is False
        assert "Definitive" in result.evidence_basis

    def test_strong_classification_with_hi_3_grants_full_strength(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Strong",
            clingen_haploinsufficiency_score=3,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.VERY_STRONG
        assert result.requires_clinical_input is False

    def test_definitive_classification_with_hi_3_in_last_exon_downgrades_to_strong(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="23/23"),
            clingen_classification="Definitive",
            clingen_haploinsufficiency_score=3,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.STRONG

    def test_moderate_classification_with_hi_2_grants_fixed_strong(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Moderate",
            clingen_haploinsufficiency_score=2,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.STRONG
        assert result.requires_clinical_input is False

    def test_moderate_classification_with_hi_2_is_not_downgraded_by_last_exon(self) -> None:
        """The Moderate+HI2 tier is a fixed Strong outcome, unlike the
        exon-position downgrade that applies to the Definitive/Strong+HI3 tier."""
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="23/23"),
            clingen_classification="Moderate",
            clingen_haploinsufficiency_score=2,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.STRONG

    def test_limited_classification_drops_to_supporting_and_flags(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Limited",
            clingen_haploinsufficiency_score=1,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.requires_clinical_input is True

    def test_definitive_classification_but_hi_score_below_3_drops_to_supporting(self) -> None:
        """Definitive/Strong gene-disease validity alone is not sufficient — the
        dosage map must independently confirm haploinsufficiency (HI score 3);
        a Definitive classification could describe a gain-of-function mechanism."""
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Definitive",
            clingen_haploinsufficiency_score=1,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.requires_clinical_input is True

    def test_disputed_classification_flags_particular_caution(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Disputed",
            clingen_haploinsufficiency_score=0,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.requires_clinical_input is True
        assert result.notes is not None
        assert "Disputed" in result.notes

    def test_refuted_classification_flags_particular_caution(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Refuted",
            clingen_haploinsufficiency_score=0,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.notes is not None
        assert "Refuted" in result.notes

    def test_definitive_recessive_gene_grants_full_strength_flagged_for_review(self) -> None:
        """Dosage code 30 (autosomal recessive) normalises to HI score 0, but
        biallelic loss of function is the expected mechanism for such genes."""
        ctx = CriterionContext(
            variant=_variant(gene_symbol="CFTR"),
            vep_annotation=_annotation(gene_symbol="CFTR", exon="5/27"),
            clingen_classification="Definitive",
            clingen_haploinsufficiency_score=0,
            clingen_recessive_gene=True,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.VERY_STRONG
        assert result.requires_clinical_input is True
        assert result.notes is not None and "recessive" in result.notes

    def test_recessive_gene_still_downgraded_in_last_exon(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="CFTR"),
            vep_annotation=_annotation(gene_symbol="CFTR", exon="27/27"),
            clingen_classification="Definitive",
            clingen_haploinsufficiency_score=0,
            clingen_recessive_gene=True,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.STRONG
        assert result.requires_clinical_input is True

    def test_limited_recessive_gene_stays_supporting(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="SOMEGENE"),
            vep_annotation=_annotation(gene_symbol="SOMEGENE", exon="5/23"),
            clingen_classification="Limited",
            clingen_haploinsufficiency_score=0,
            clingen_recessive_gene=True,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.SUPPORTING

    def test_no_clingen_data_falls_back_to_bundled_list(self) -> None:
        """clingen_classification=None (the default) means ClinGen was not
        consulted or had no curation for this gene; PVS1 falls back to the
        original bundled-list/gnomAD-constraint heuristic rather than
        dropping straight to Supporting."""
        ctx = CriterionContext(
            variant=_variant(gene_symbol="BRCA1"),
            vep_annotation=_annotation(gene_symbol="BRCA1", exon="5/23"),
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.VERY_STRONG
        assert result.requires_clinical_input is False
        assert "no ClinGen curation" in result.evidence_basis

    def test_no_clingen_data_and_unestablished_gene_drops_to_supporting(self) -> None:
        ctx = CriterionContext(
            variant=_variant(gene_symbol="UNKNOWNGENE"),
            vep_annotation=_annotation(gene_symbol="UNKNOWNGENE", exon="5/23"),
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.requires_clinical_input is True

    def test_explicit_clingen_classification_overrides_bundled_list(self) -> None:
        """A gene on the bundled fallback list (e.g. BRCA1) with an explicit
        ClinGen result that fails the tier thresholds must NOT fall back to
        the bundled list's automatic full-strength grant — an explicit
        ClinGen classification always takes precedence over the fallback."""
        ctx = CriterionContext(
            variant=_variant(gene_symbol="BRCA1"),
            vep_annotation=_annotation(gene_symbol="BRCA1", exon="5/23"),
            clingen_classification="Limited",
            clingen_haploinsufficiency_score=0,
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.strength == EvidenceStrength.SUPPORTING
        assert result.requires_clinical_input is True


class TestSpliceSitePositions:
    def test_canonical_splice_donor_is_evaluated_as_lof(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="splice_donor_variant", exon="10/23"
            ),
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.VERY_STRONG

    def test_canonical_splice_acceptor_is_evaluated_as_lof(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="splice_acceptor_variant", exon="10/23"
            ),
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is True

    def test_non_canonical_splice_region_variant_is_not_applicable(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(most_severe_consequence="splice_region_variant"),
        )

        result = PVS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False
        assert result.notes is not None
        assert "SpliceAI" in result.notes
