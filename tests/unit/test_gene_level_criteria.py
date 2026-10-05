"""Unit tests for alleleiq.scoring.criteria.gene_level (PP2/BP1), using small
synthetic fixtures under tests/fixtures/gene_annotations/ rather than the real
bundled 50-gene curated file (see data/gene_annotations/README.md).
"""

from __future__ import annotations

from pathlib import Path

from alleleiq.clients.vep import VEPAnnotation
from alleleiq.config import Settings
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.scoring.criteria.base import CriterionContext
from alleleiq.scoring.criteria.gene_level import (
    BP1Criterion,
    PP2Criterion,
    load_gene_annotations,
)

_FIXTURES = Path(__file__).parent.parent / "fixtures" / "gene_annotations"
_COORDS = GenomicCoordinates(chrom="1", pos=100, ref="C", alt="T")


def _settings(path: Path) -> Settings:
    return Settings(gene_annotations_path=str(path))


def _variant(gene_symbol: str) -> Variant:
    return Variant(coordinates=_COORDS, variant_type=VariantType.SNV, gene_symbol=gene_symbol)


def _annotation(gene_symbol: str, consequence: str) -> VEPAnnotation:
    return VEPAnnotation(
        variant_id=_COORDS.variant_id,
        most_severe_consequence=consequence,
        gene_symbol=gene_symbol,
    )


def _context(gene_symbol: str, consequence: str, path: Path) -> CriterionContext:
    return CriterionContext(
        variant=_variant(gene_symbol),
        vep_annotation=_annotation(gene_symbol, consequence),
        settings=_settings(path),
    )


class TestPP2Criterion:
    def test_missense_in_pp2_applicable_gene_is_applied(self) -> None:
        ctx = _context("PP2GENE", "missense_variant", _FIXTURES / "missense_constraint.json")

        result = PP2Criterion().evaluate(ctx)

        assert result.met is True
        assert result.requires_clinical_input is False
        assert "PP2GENE" in result.evidence_basis

    def test_missense_in_bp1_only_gene_is_not_applied(self) -> None:
        ctx = _context("BP1GENE", "missense_variant", _FIXTURES / "missense_constraint.json")

        result = PP2Criterion().evaluate(ctx)

        assert result.met is False

    def test_gene_absent_from_curated_list_is_not_assessed(self) -> None:
        ctx = _context("UNKNOWNGENE", "missense_variant", _FIXTURES / "missense_constraint.json")

        result = PP2Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False
        assert result.notes is not None and "expert curation" in result.notes

    def test_non_missense_consequence_is_not_applicable_even_in_pp2_gene(self) -> None:
        ctx = _context("PP2GENE", "stop_gained", _FIXTURES / "missense_constraint.json")

        result = PP2Criterion().evaluate(ctx)

        assert result.met is False
        assert "missense" in result.evidence_basis.lower()

    def test_missing_annotation_file_yields_not_assessed_for_every_gene(self) -> None:
        ctx = _context("PP2GENE", "missense_variant", _FIXTURES / "does_not_exist.json")

        result = PP2Criterion().evaluate(ctx)

        assert result.met is False

    def test_no_gene_symbol_available_is_not_assessed(self) -> None:
        ctx = CriterionContext(
            variant=_variant(""),
            vep_annotation=None,
            settings=_settings(_FIXTURES / "missense_constraint.json"),
        )

        result = PP2Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False


class TestBP1Criterion:
    def test_missense_in_bp1_applicable_gene_is_applied(self) -> None:
        ctx = _context("BP1GENE", "missense_variant", _FIXTURES / "missense_constraint.json")

        result = BP1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.requires_clinical_input is False
        assert "BP1GENE" in result.evidence_basis

    def test_missense_in_pp2_only_gene_is_not_applied(self) -> None:
        ctx = _context("PP2GENE", "missense_variant", _FIXTURES / "missense_constraint.json")

        result = BP1Criterion().evaluate(ctx)

        assert result.met is False

    def test_gene_absent_from_curated_list_is_not_assessed(self) -> None:
        ctx = _context("UNKNOWNGENE", "missense_variant", _FIXTURES / "missense_constraint.json")

        result = BP1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False

    def test_non_missense_consequence_is_not_applicable_even_in_bp1_gene(self) -> None:
        ctx = _context("BP1GENE", "frameshift_variant", _FIXTURES / "missense_constraint.json")

        result = BP1Criterion().evaluate(ctx)

        assert result.met is False
        assert "missense" in result.evidence_basis.lower()


def _domain_context(
    hgvs_p: str | None, spliceai: float | None = None, gene: str = "DOMAINGENE"
) -> CriterionContext:
    return CriterionContext(
        variant=_variant(gene),
        vep_annotation=VEPAnnotation(
            variant_id=_COORDS.variant_id,
            most_severe_consequence="missense_variant",
            gene_symbol=gene,
            hgvs_p=hgvs_p,
            spliceai_ds_ag=spliceai,
        ),
        settings=_settings(_FIXTURES / "missense_constraint.json"),
    )


class TestBP1Withholding:
    def test_missense_inside_key_domain_is_withheld_for_review(self) -> None:
        result = BP1Criterion().evaluate(_domain_context("ENSP1.1:p.Arg150Gly"))

        assert result.met is False
        assert result.requires_clinical_input is True
        assert "codon 150" in result.evidence_basis

    def test_missense_outside_key_domains_applies(self) -> None:
        result = BP1Criterion().evaluate(_domain_context("ENSP1.1:p.Arg350Gly"))

        assert result.met is True

    def test_missing_protein_position_in_domain_gene_is_withheld(self) -> None:
        result = BP1Criterion().evaluate(_domain_context(None))

        assert result.met is False
        assert result.requires_clinical_input is True

    def test_predicted_splicing_effect_withholds_bp1(self) -> None:
        result = BP1Criterion().evaluate(
            _domain_context("ENSP1.1:p.Arg350Gly", spliceai=0.6, gene="BP1GENE")
        )

        assert result.met is False
        assert result.requires_clinical_input is True
        assert "splicing" in result.evidence_basis

    def test_bundled_brca1_annotation_withholds_bp1_in_the_ring_domain(self) -> None:
        ctx = CriterionContext(
            variant=_variant("BRCA1"),
            vep_annotation=VEPAnnotation(
                variant_id=_COORDS.variant_id,
                most_severe_consequence="missense_variant",
                gene_symbol="BRCA1",
                hgvs_p="ENSP00000350283.3:p.Cys61Gly",
            ),
        )

        result = BP1Criterion().evaluate(ctx)

        assert result.met is False
        assert "RING" in result.evidence_basis


class TestLoadGeneAnnotations:
    def test_loads_valid_entries(self) -> None:
        result = load_gene_annotations(str(_FIXTURES / "missense_constraint.json"))

        assert "PP2GENE" in result
        assert result["PP2GENE"].pp2_applicable is True

    def test_missing_file_returns_empty_mapping(self) -> None:
        result = load_gene_annotations(str(_FIXTURES / "does_not_exist.json"))

        assert result == {}

    def test_incoherent_both_true_entry_is_skipped_but_siblings_still_load(self) -> None:
        result = load_gene_annotations(str(_FIXTURES / "incoherent.json"))

        assert "BOTHGENE" not in result
        assert "GOODGENE" in result
