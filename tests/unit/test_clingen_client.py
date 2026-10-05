"""Unit tests for alleleiq.clients.clingen.ClinGenClient, using small synthetic
fixtures under tests/fixtures/clingen/ rather than the real bundled ClinGen
snapshots (which are large and change over time as ClinGen curates more genes).
"""

from __future__ import annotations

from pathlib import Path

from alleleiq.clients.clingen import ClinGenClient
from alleleiq.config import Settings

_FIXTURES = Path(__file__).parent.parent / "fixtures" / "clingen"


def _client(**overrides: object) -> ClinGenClient:
    settings = Settings(
        clingen_gene_validity_path=str(_FIXTURES / "gene_disease_validity.csv"),
        clingen_dosage_path=str(_FIXTURES / "dosage_sensitivity.tsv"),
        **overrides,  # type: ignore[arg-type]
    )
    return ClinGenClient(settings)


class TestGeneDiseaseValidity:
    def test_single_entry_gene_returns_its_classification(self) -> None:
        result = _client().get_disease_validity("FAKEGENE1")

        assert result is not None
        assert result.classification == "Definitive"
        assert result.gene_symbol == "FAKEGENE1"
        assert result.gcep == "Fake Gene Curation Expert Panel"

    def test_lookup_is_case_insensitive(self) -> None:
        result = _client().get_disease_validity("fakegene1")

        assert result is not None
        assert result.gene_symbol == "FAKEGENE1"

    def test_multi_entry_gene_returns_strongest_classification(self) -> None:
        """MULTIGENE has a Limited row and a Definitive row; Definitive wins."""
        result = _client().get_disease_validity("MULTIGENE")

        assert result is not None
        assert result.classification == "Definitive"
        assert result.disease_label == "multi disease B"

    def test_gene_absent_from_file_returns_none(self) -> None:
        assert _client().get_disease_validity("NOTAGENE") is None

    def test_load_gene_validity_maps_every_gene_to_best_classification(self) -> None:
        mapping = _client().load_gene_validity()

        assert mapping["FAKEGENE1"] == "Definitive"
        assert mapping["MULTIGENE"] == "Definitive"
        assert mapping["DISPUTEDGENE"] == "Disputed"
        assert "NOTAGENE" not in mapping

    def test_missing_validity_file_returns_empty_and_does_not_raise(self, tmp_path: Path) -> None:
        settings = Settings(
            clingen_gene_validity_path=str(tmp_path / "missing.csv"),
            clingen_dosage_path=str(_FIXTURES / "dosage_sensitivity.tsv"),
        )
        client = ClinGenClient(settings)

        assert client.load_gene_validity() == {}
        assert client.get_disease_validity("FAKEGENE1") is None


class TestHaploinsufficiencyScore:
    def test_score_3_returned_directly(self) -> None:
        assert _client().get_haploinsufficiency_score("FAKEGENE1") == 3

    def test_score_2_returned_directly(self) -> None:
        assert _client().get_haploinsufficiency_score("MULTIGENE") == 2

    def test_special_code_30_normalises_to_zero(self) -> None:
        """30 = 'gene associated with autosomal recessive phenotype', not
        evidence for haploinsufficiency."""
        assert _client().get_haploinsufficiency_score("SPECIALCODE") == 0

    def test_special_code_40_normalises_to_zero(self) -> None:
        """40 = 'dosage sensitivity unlikely'."""
        assert _client().get_haploinsufficiency_score("UNLIKELYGENE") == 0

    def test_gene_absent_from_dosage_map_returns_zero(self) -> None:
        assert _client().get_haploinsufficiency_score("NOTAGENE") == 0

    def test_dosage_sensitivity_preserves_raw_code(self) -> None:
        dosage = _client().get_dosage_sensitivity("SPECIALCODE")

        assert dosage is not None
        assert dosage.haploinsufficiency_score == 0
        assert dosage.haploinsufficiency_raw_code == 30
        assert "autosomal recessive" in dosage.haploinsufficiency_description

    def test_missing_dosage_file_returns_zero_and_does_not_raise(self, tmp_path: Path) -> None:
        settings = Settings(
            clingen_gene_validity_path=str(_FIXTURES / "gene_disease_validity.csv"),
            clingen_dosage_path=str(tmp_path / "missing.tsv"),
        )
        client = ClinGenClient(settings)

        assert client.get_haploinsufficiency_score("FAKEGENE1") == 0
        assert client.get_dosage_sensitivity("FAKEGENE1") is None


class TestSnapshotVersions:
    def test_reads_creation_dates_from_both_files(self) -> None:
        assert _client().snapshot_versions() == {
            "clingen_gene_validity": "2026-01-01",
            "clingen_dosage": "2026-01-01",
        }

    def test_missing_files_are_omitted(self, tmp_path: Path) -> None:
        settings = Settings(
            clingen_gene_validity_path=str(tmp_path / "missing.csv"),
            clingen_dosage_path=str(tmp_path / "missing.tsv"),
        )
        assert ClinGenClient(settings).snapshot_versions() == {}
