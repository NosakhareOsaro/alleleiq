"""Unit tests for alleleiq.scoring.criteria.hotspot (PM1), using small synthetic
fixtures under tests/fixtures/hotspots/ rather than the real bundled 6-gene
curated file (see data/hotspots/README.md). All ClinVar/UniProt data here is
constructed directly as CriterionContext fields — no live network calls.
"""

from __future__ import annotations

from pathlib import Path

from alleleiq.clients.clinvar import ClinVarPositionalMatch
from alleleiq.clients.uniprot import ProteinRegion
from alleleiq.clients.vep import VEPAnnotation
from alleleiq.config import Settings
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.scoring.criteria.base import CriterionContext
from alleleiq.scoring.criteria.hotspot import PM1Criterion, load_curated_hotspots

_FIXTURES = Path(__file__).parent.parent / "fixtures" / "hotspots"
_COORDS = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")


def _variant() -> Variant:
    return Variant(coordinates=_COORDS, variant_type=VariantType.SNV, gene_symbol="HOTSPOTGENE")


def _annotation(hgvs_p: str | None = "ENSP00000000001.1:p.Arg175His") -> VEPAnnotation:
    return VEPAnnotation(
        variant_id=_COORDS.variant_id,
        most_severe_consequence="missense_variant",
        gene_symbol="HOTSPOTGENE",
        hgvs_p=hgvs_p,
    )


def _match(
    *,
    clinical_significance: str,
    star_rating: int = 3,
    accession: str = "VCV000000001",
    location: tuple[int, str, str] | None = (43094700, "C", "T"),
) -> ClinVarPositionalMatch:
    """Build a nearby match. ``location`` defaults to a different variant from
    ``_COORDS``; pass ``_COORDS`` values for the variant's own record."""
    position, ref, alt = location if location is not None else (None, None, None)
    return ClinVarPositionalMatch(
        accession=accession,
        clinical_significance=clinical_significance,
        review_status="reviewed by expert panel" if star_rating >= 3 else "criteria provided",
        star_rating=star_rating,
        submitter_count=1,
        variation_name="NM_000000.1(HOTSPOTGENE):c.524G>A (p.Arg175His)",
        protein_change="Arg175His",
        genomic_position=position,
        reference_allele=ref,
        alternate_allele=alt,
    )


def _region(start: int = 170, end: int = 180, source: str = "test") -> ProteinRegion:
    return ProteinRegion(
        start_aa=start, end_aa=end, description="Fixture domain region", source=source
    )


class TestPM1CuratedAndDomainHits:
    def test_curated_hotspot_region_applies_pm1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[_region(source="curated hotspot list")],
            pm1_nearby_clinvar_matches=[],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.requires_clinical_input is True
        assert "curated hotspot list" in result.evidence_basis

    def test_uniprot_domain_with_no_nearby_benign_applies_pm1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[_region(source="UniProt P00000 (HOTSPOTGENE)")],
            pm1_nearby_clinvar_matches=[_match(clinical_significance="Pathogenic", star_rating=3)],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is True
        assert result.requires_clinical_input is True
        assert "UniProt" in result.evidence_basis

    def test_domain_region_with_disqualifying_benign_variants_withholds_pm1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[_region()],
            pm1_nearby_clinvar_matches=[
                _match(clinical_significance="Benign", accession="VCV1"),
                _match(clinical_significance="Likely benign", accession="VCV2"),
            ],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True
        assert "Benign" in result.evidence_basis or "benign" in result.evidence_basis

    def test_single_benign_variant_does_not_disqualify_by_default(self) -> None:
        # Settings default pm1_max_allowed_benign_count=1: exactly one tolerated benign
        # record must not withhold PM1, only exceeding that count does.
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[_region()],
            pm1_nearby_clinvar_matches=[_match(clinical_significance="Benign", accession="VCV1")],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is True


class TestPM1ClinVarDensityHotspot:
    def test_three_established_pathogenic_hits_with_no_domain_applies_pm1(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[],
            pm1_nearby_clinvar_matches=[
                _match(clinical_significance="Pathogenic", star_rating=3, accession="VCV1"),
                _match(clinical_significance="Likely pathogenic", star_rating=2, accession="VCV2"),
                _match(clinical_significance="Pathogenic", star_rating=2, accession="VCV3"),
            ],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is True
        assert "established" in result.evidence_basis

    def test_two_established_pathogenic_hits_is_below_threshold(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[],
            pm1_nearby_clinvar_matches=[
                _match(clinical_significance="Pathogenic", star_rating=3, accession="VCV1"),
                _match(clinical_significance="Pathogenic", star_rating=2, accession="VCV2"),
            ],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True

    def test_low_star_pathogenic_hits_do_not_count_toward_density(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[],
            pm1_nearby_clinvar_matches=[
                _match(clinical_significance="Pathogenic", star_rating=1, accession="VCV1"),
                _match(clinical_significance="Pathogenic", star_rating=1, accession="VCV2"),
                _match(clinical_significance="Pathogenic", star_rating=1, accession="VCV3"),
            ],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False


class TestPM1ExcludesTheVariantItself:
    def test_variants_own_record_does_not_count_toward_density(self) -> None:
        own = (_COORDS.pos, _COORDS.ref, _COORDS.alt)
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[],
            pm1_nearby_clinvar_matches=[
                _match(clinical_significance="Pathogenic", accession="VCV-SELF", location=own),
                _match(clinical_significance="Pathogenic", accession="VCV1"),
                _match(clinical_significance="Pathogenic", accession="VCV2"),
            ],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False

    def test_same_protein_change_without_location_is_treated_as_the_variant(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[],
            pm1_nearby_clinvar_matches=[
                _match(clinical_significance="Pathogenic", accession="VCV-SELF", location=None),
                _match(clinical_significance="Pathogenic", accession="VCV1"),
                _match(clinical_significance="Pathogenic", accession="VCV2"),
            ],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False


class TestPM1NotAssessed:
    def test_no_domain_no_density_no_curated_is_not_assessed(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[],
            pm1_nearby_clinvar_matches=[],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True
        assert "not assessed" in result.evidence_basis.lower()

    def test_no_protein_change_available_is_not_assessed(self) -> None:
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(hgvs_p=None),
            pm1_domain_regions=[_region()],
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True

    def test_no_vep_annotation_at_all_is_not_assessed(self) -> None:
        ctx = CriterionContext(variant=_variant(), vep_annotation=None)

        result = PM1Criterion().evaluate(ctx)

        assert result.met is False
        assert result.requires_clinical_input is True


class TestLoadCuratedHotspots:
    def test_loads_valid_entries(self) -> None:
        result = load_curated_hotspots(str(_FIXTURES / "curated_hotspots.json"))

        assert "HOTSPOTGENE" in result
        assert result["HOTSPOTGENE"][0].start_aa == 170

    def test_missing_file_returns_empty_mapping(self) -> None:
        result = load_curated_hotspots(str(_FIXTURES / "does_not_exist.json"))

        assert result == {}

    def test_malformed_entry_is_skipped_but_siblings_still_load(self) -> None:
        result = load_curated_hotspots(str(_FIXTURES / "malformed_entry.json"))

        assert "BADGENE" not in result
        assert "GOODGENE" in result


class TestPM1SettingsBoundaries:
    def test_custom_min_pathogenic_hotspot_count_is_respected(self) -> None:
        settings = Settings(pm1_min_pathogenic_hotspot_count=2)
        ctx = CriterionContext(
            variant=_variant(),
            vep_annotation=_annotation(),
            pm1_domain_regions=[],
            pm1_nearby_clinvar_matches=[
                _match(clinical_significance="Pathogenic", star_rating=3, accession="VCV1"),
                _match(clinical_significance="Pathogenic", star_rating=2, accession="VCV2"),
            ],
            settings=settings,
        )

        result = PM1Criterion().evaluate(ctx)

        assert result.met is True
