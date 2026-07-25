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
