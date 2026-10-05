"""Unit tests for alleleiq.scoring.criteria.same_position (PS1, PM5)."""

from __future__ import annotations

from alleleiq.clients.clinvar import ClinVarPositionalMatch
from alleleiq.clients.vep import VEPAnnotation
from alleleiq.models.evidence import EvidenceStrength
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.scoring.criteria.base import CriterionContext
from alleleiq.scoring.criteria.same_position import (
    PM5Criterion,
    PS1Criterion,
    extract_protein_change,
)

_COORDS = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")


def _variant() -> Variant:
    return Variant(coordinates=_COORDS, variant_type=VariantType.SNV, gene_symbol="BRCA1")


def _annotation(**overrides: object) -> VEPAnnotation:
    defaults: dict[str, object] = {
        "variant_id": _COORDS.variant_id,
        "most_severe_consequence": "missense_variant",
        "gene_symbol": "BRCA1",
        "hgvs_p": "ENSP00000418960.2:p.Arg1699Gln",
    }
    defaults.update(overrides)
    return VEPAnnotation(**defaults)  # type: ignore[arg-type]


def _match(
    *,
    protein_change: str | None,
    clinical_significance: str = "Pathogenic",
    star_rating: int = 3,
    accession: str = "VCV000017661",
    location: tuple[int, str, str] | None = (43094693, "C", "T"),
) -> ClinVarPositionalMatch:
    """Build a match. ``location`` defaults to a different nucleotide change from
    ``_COORDS``; pass ``_COORDS`` values for the variant's own record, or None for
    a record without VCF-style location detail."""
    position, ref, alt = location if location is not None else (None, None, None)
    return ClinVarPositionalMatch(
        accession=accession,
        clinical_significance=clinical_significance,
        review_status="reviewed by expert panel" if star_rating >= 3 else "criteria provided",
        star_rating=star_rating,
        submitter_count=1,
        variation_name=(
            f"NM_007294.4(BRCA1):c.5096G>A (p.{protein_change})"
            if protein_change
            else "NM_007294.4(BRCA1):c.5096G>A"
        ),
        protein_change=protein_change,
        genomic_position=position,
        reference_allele=ref,
        alternate_allele=alt,
    )


class TestExtractProteinChange:
    def test_extracts_from_vep_transcript_prefixed_form(self) -> None:
        assert extract_protein_change("ENSP00000418960.2:p.Arg1699Gln") == "Arg1699Gln"

    def test_extracts_from_clinvar_parenthesised_form(self) -> None:
        assert extract_protein_change("NM_007294.4(BRCA1):c.5096G>A (p.Arg1699Gln)") == "Arg1699Gln"

    def test_none_when_no_protein_change_present(self) -> None:
        assert extract_protein_change(None) is None
        assert extract_protein_change("NM_007294.4:c.5096G>A") is None


class TestPS1SameAminoAcidChange:
    def test_same_change_three_star_pathogenic_applies_ps1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Gln", star_rating=3),
            ],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.STRONG
        assert result.requires_clinical_input is True

    def test_same_change_one_star_does_not_apply_ps1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Gln", star_rating=1),
            ],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False
        assert "star" in result.evidence_basis

    def test_no_clinvar_data_at_position_not_assessed(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False
        assert result.notes is not None and "Not assessed" in result.notes

    def test_no_protein_change_available_requires_clinical_input(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(hgvs_p=None),
            same_position_clinvar_matches=[_match(protein_change="Arg1699Gln")],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True

    def test_different_change_at_position_does_not_apply_ps1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Trp", star_rating=3),
            ],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is False

    def test_splice_signal_withholds_ps1_despite_matching_change(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(spliceai_ds_dl=0.85),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Gln", star_rating=3),
            ],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True
        assert "splic" in result.evidence_basis.lower()


class TestPS1ExcludesTheVariantItself:
    def test_variants_own_clinvar_record_does_not_satisfy_ps1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(
                    protein_change="Arg1699Gln",
                    star_rating=3,
                    location=(_COORDS.pos, _COORDS.ref, _COORDS.alt),
                ),
            ],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False
        assert "different variant" in result.evidence_basis

    def test_own_record_ignored_but_distinct_established_record_still_applies(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(
                    protein_change="Arg1699Gln",
                    accession="VCV-SELF",
                    location=(_COORDS.pos, _COORDS.ref, _COORDS.alt),
                ),
                _match(protein_change="Arg1699Gln", accession="VCV-OTHER"),
            ],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is True
        assert "VCV-OTHER" in result.evidence_basis

    def test_record_without_location_detail_is_withheld_for_review(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Gln", star_rating=3, location=None),
            ],
        )

        result = PS1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True
        assert "location detail" in result.evidence_basis


class TestPM5DifferentAminoAcidChangeSameCodon:
    def test_different_change_same_codon_two_star_applies_pm5(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Trp", star_rating=2),
            ],
        )

        result = PM5Criterion().evaluate(ctx)

        assert result.met is True
        assert result.strength == EvidenceStrength.MODERATE
        assert result.requires_clinical_input is True

    def test_same_change_is_ps1_territory_not_pm5(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Gln", star_rating=3),
            ],
        )

        result = PM5Criterion().evaluate(ctx)

        assert result.met is False

    def test_no_clinvar_data_at_position_not_assessed(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[],
        )

        result = PM5Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is False

    def test_different_codon_does_not_apply_pm5(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Val1838Met", star_rating=3),
            ],
        )

        result = PM5Criterion().evaluate(ctx)

        assert result.met is False

    def test_low_star_different_change_does_not_apply_pm5(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            same_position_clinvar_matches=[
                _match(protein_change="Arg1699Trp", star_rating=1),
            ],
        )

        result = PM5Criterion().evaluate(ctx)

        assert result.met is False
