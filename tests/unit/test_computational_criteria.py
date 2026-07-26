"""Unit tests for alleleiq.scoring.criteria.computational (PP3, BP4)."""

from __future__ import annotations

from alleleiq.clients.vep import VEPAnnotation
from alleleiq.models.evidence import ACMGCode, EvidenceStrength
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.scoring.criteria.base import CriterionContext
from alleleiq.scoring.criteria.computational import BP4Criterion, PP3Criterion

_COORDS = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")


def _variant() -> Variant:
    return Variant(coordinates=_COORDS, variant_type=VariantType.SNV, gene_symbol="BRCA1")


def _annotation(**overrides: object) -> VEPAnnotation:
    defaults: dict[str, object] = {
        "variant_id": _COORDS.variant_id,
        "most_severe_consequence": "missense_variant",
    }
    defaults.update(overrides)
    return VEPAnnotation(**defaults)  # type: ignore[arg-type]


class TestPP3:
    def test_high_revel_triggers_moderate_strength(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(revel_score=0.95))

        result = PP3Criterion().evaluate(ctx)

        assert result.code == ACMGCode.PP3
        assert result.met is True
        assert result.strength == EvidenceStrength.MODERATE
        assert result.requires_clinical_input is False

    def test_revel_at_supporting_threshold_only(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(revel_score=0.75))

        result = PP3Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.SUPPORTING

    def test_revel_below_threshold_not_met(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(revel_score=0.50))

        result = PP3Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False

    def test_high_spliceai_triggers_moderate_strength(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="splice_donor_variant",
                revel_score=None,
                spliceai_ds_dl=0.6,
            ),
        )

        result = PP3Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.MODERATE

    def test_low_spliceai_supporting_only(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="splice_donor_variant",
                revel_score=None,
                spliceai_ds_ag=0.25,
            ),
        )

        result = PP3Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.SUPPORTING

    def test_missense_with_no_revel_requires_clinical_input(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="missense_variant", revel_score=None, spliceai_ds_ag=0.3
            ),
        )

        result = PP3Criterion().evaluate(ctx)

        # SpliceAI is present and evaluated, but REVEL was expected (missense) and absent.
        assert result.requires_clinical_input is True

    def test_non_missense_missing_revel_is_not_flagged(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="splice_donor_variant",
                revel_score=None,
                spliceai_ds_ag=0.3,
            ),
        )

        result = PP3Criterion().evaluate(ctx)

        assert result.requires_clinical_input is False

    def test_no_scores_at_all_requires_clinical_input(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(revel_score=None))

        result = PP3Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True

    def test_no_vep_annotation_requires_clinical_input(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=None)

        result = PP3Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True


class TestBP4:
    def test_low_revel_triggers_benign_supporting(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(revel_score=0.05))

        result = BP4Criterion().evaluate(ctx)

        assert result.code == ACMGCode.BP4
        assert result.met is True
        assert result.strength == EvidenceStrength.BENIGN_SUPPORTING

    def test_revel_above_bp4_threshold_not_met(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=_annotation(revel_score=0.5))

        result = BP4Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False

    def test_low_spliceai_triggers_benign_supporting(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="splice_donor_variant",
                revel_score=None,
                spliceai_ds_ag=0.01,
                spliceai_ds_al=0.0,
                spliceai_ds_dg=0.0,
                spliceai_ds_dl=0.0,
            ),
        )

        result = BP4Criterion().evaluate(ctx)

        assert result.met is True

    def test_missense_with_no_revel_requires_clinical_input(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="missense_variant", revel_score=None
            ),
        )

        result = BP4Criterion().evaluate(ctx)

        assert result.requires_clinical_input is True

    def test_no_scores_at_all_requires_clinical_input(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=None)

        result = BP4Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True
