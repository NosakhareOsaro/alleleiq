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
                most_severe_consequence="splice_region_variant",
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
                most_severe_consequence="splice_region_variant",
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
                most_severe_consequence="synonymous_variant",
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


class TestConsequenceAwareness:
    """PP3/BP4 use only the scores that fit the variant's consequence (ClinGen SVI)."""

    def test_canonical_splice_variant_does_not_double_count_with_pvs1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="splice_donor_variant", spliceai_ds_dl=0.91
            ),
        )

        result = PP3Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False
        assert "PVS1" in result.evidence_basis

    def test_stop_gained_never_gets_bp4_from_a_low_spliceai_score(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(most_severe_consequence="stop_gained", spliceai_ds_ag=0.01),
        )

        result = BP4Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False

    def test_missense_with_intermediate_revel_and_no_splice_signal_is_not_bp4(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(revel_score=0.50, spliceai_ds_ag=0.01),
        )

        result = BP4Criterion().evaluate(ctx)

        assert result.met is False

    def test_missense_with_benign_revel_but_splice_signal_is_not_bp4(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(revel_score=0.05, spliceai_ds_dg=0.35),
        )

        result = BP4Criterion().evaluate(ctx)

        assert result.met is False

    def test_missense_with_benign_revel_and_no_splice_signal_is_bp4(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(revel_score=0.05, spliceai_ds_ag=0.02),
        )

        result = BP4Criterion().evaluate(ctx)

        assert result.met is True

    def test_revel_is_ignored_for_non_missense_consequences(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(
                most_severe_consequence="synonymous_variant",
                revel_score=0.95,
                spliceai_ds_ag=0.01,
            ),
        )

        pp3 = PP3Criterion().evaluate(ctx)
        bp4 = BP4Criterion().evaluate(ctx)

        assert pp3.met is False
        assert "REVEL" not in pp3.evidence_basis
        assert bp4.met is True


class TestREVELClientIntegration:
    """context.revel_score (as populated by pipeline.build_evidence_report via
    REVELClient) takes precedence over vep_annotation.revel_score, since it may
    have been resolved from AlleleIQ's local REVEL index rather than VEP's
    incomplete REST coverage.
    """

    def test_revel_score_on_context_overrides_annotation_for_pp3(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(revel_score=0.50),
            revel_score=0.95,
        )

        result = PP3Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.MODERATE

    def test_revel_score_on_context_overrides_annotation_for_bp4(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(revel_score=0.50),
            revel_score=0.05,
        )

        result = BP4Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.BENIGN_SUPPORTING

    def test_falls_back_to_annotation_when_context_score_absent(self) -> None:
        ctx = CriterionContext(
            variant=_variant(), vep_annotation=_annotation(revel_score=0.95), revel_score=None
        )

        result = PP3Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.MODERATE

    def test_context_score_used_even_without_annotation(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=None, revel_score=0.95)

        result = PP3Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.MODERATE
        # No annotation means consequence is unknown, so the missense-specific
        # "expected but absent" flag cannot apply; met=True either way here.
        assert result.requires_clinical_input is False
