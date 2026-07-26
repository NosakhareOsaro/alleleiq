"""Unit tests for alleleiq.clients.revel.REVELClient.

Local lookups use small, real bgzipped/tabix-indexed fixtures (see
tests/fixtures/revel/) — building a real tabix index is cheap and
exercises the actual pysam query path rather than mocking it away.
Remote (VEP REST) fallback is exercised via a stubbed VEPClient, exactly
as test_vep_client.py mocks httpx transport for VEPClient itself.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from alleleiq.clients import revel as revel_module
from alleleiq.clients.revel import REVELClient
from alleleiq.clients.vep import VEPAnnotation
from alleleiq.config import Settings
from alleleiq.models.variant import GenomicCoordinates

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures" / "revel"
BARE_CHROM_INDEX = FIXTURES_DIR / "revel_fixture.tsv.gz"
CHR_PREFIXED_INDEX = FIXTURES_DIR / "revel_fixture_chr_prefixed.tsv.gz"


def _settings(revel_data_path: str | None) -> Settings:
    return Settings(revel_data_path=revel_data_path)


def _annotation(revel_score: float | None) -> VEPAnnotation:
    return VEPAnnotation(variant_id="stub", revel_score=revel_score)


class _StubVEPClient:
    """Replaces alleleiq.clients.revel.VEPClient for fallback-path tests."""

    def __init__(self, annotation: VEPAnnotation | None) -> None:
        self._annotation = annotation

    def __call__(self, _settings: Settings) -> _StubVEPClient:
        return self

    async def __aenter__(self) -> _StubVEPClient:
        return self

    async def __aexit__(self, *_args: object) -> None:
        return None

    async def fetch_annotation(self, _coordinates: GenomicCoordinates) -> VEPAnnotation | None:
        return self._annotation


class TestLocalLookup:
    async def test_exact_match_found(self) -> None:
        client = REVELClient(_settings(str(BARE_CHROM_INDEX)))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        result = await client.get_score(coords)

        assert result.score == pytest.approx(0.912)
        assert result.source == "local"
        client.close()

    async def test_multiallelic_matches_exact_alt_not_just_position(self) -> None:
        client = REVELClient(_settings(str(BARE_CHROM_INDEX)))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="T")

        result = await client.get_score(coords)

        assert result.score == pytest.approx(0.204)
        assert result.source == "local"
        client.close()

    async def test_wrong_alt_at_indexed_position_falls_back(self) -> None:
        client = REVELClient(_settings(str(BARE_CHROM_INDEX)))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="C")

        result = await client.get_score(coords, vep_annotation=_annotation(0.4))

        assert result.score == pytest.approx(0.4)
        assert result.source == "vep_rest"
        client.close()

    async def test_bare_chrom_query_matches_chr_prefixed_index(self) -> None:
        client = REVELClient(_settings(str(CHR_PREFIXED_INDEX)))
        coords = GenomicCoordinates(chrom="20", pos=200000, ref="T", alt="C")

        result = await client.get_score(coords)

        assert result.score == pytest.approx(0.85)
        assert result.source == "local"
        client.close()

    async def test_tabix_handle_reused_across_repeated_lookups(self) -> None:
        client = REVELClient(_settings(str(BARE_CHROM_INDEX)))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        first = await client.get_score(coords)
        second = await client.get_score(coords)

        assert first == second
        assert first.source == "local"
        client.close()

    async def test_position_absent_from_index_falls_back(self) -> None:
        client = REVELClient(_settings(str(BARE_CHROM_INDEX)))
        coords = GenomicCoordinates(chrom="9", pos=1, ref="A", alt="G")

        result = await client.get_score(coords, vep_annotation=_annotation(0.6))

        assert result.score == pytest.approx(0.6)
        assert result.source == "vep_rest"
        client.close()

    async def test_corrupt_index_file_degrades_to_fallback_without_raising(
        self, tmp_path: Path
    ) -> None:
        bogus = tmp_path / "not_actually_bgzipped.tsv.gz"
        bogus.write_text("this is not a valid bgzip/tabix file")
        client = REVELClient(_settings(str(bogus)))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        result = await client.get_score(coords, vep_annotation=_annotation(0.3))

        assert result.score == pytest.approx(0.3)
        assert result.source == "vep_rest"
        client.close()


class TestRemoteFallback:
    async def test_missing_local_file_uses_provided_vep_annotation(self) -> None:
        client = REVELClient(_settings("./data/does_not_exist.tsv.gz"))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        result = await client.get_score(coords, vep_annotation=_annotation(0.512))

        assert result.score == pytest.approx(0.512)
        assert result.source == "vep_rest"
        client.close()

    async def test_no_local_file_and_no_annotation_queries_vep_rest(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(revel_module, "VEPClient", _StubVEPClient(_annotation(0.77)))
        client = REVELClient(_settings(None))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        result = await client.get_score(coords)

        assert result.score == pytest.approx(0.77)
        assert result.source == "vep_rest"
        client.close()

    async def test_no_data_anywhere_returns_unavailable_without_raising(self) -> None:
        client = REVELClient(_settings(None))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        result = await client.get_score(coords, vep_annotation=_annotation(None))

        assert result.score is None
        assert result.source == "unavailable"
        client.close()

    async def test_local_file_configured_but_absent_and_no_data_available(self) -> None:
        client = REVELClient(_settings(str(BARE_CHROM_INDEX)))
        coords = GenomicCoordinates(chrom="9", pos=1, ref="A", alt="G")

        result = await client.get_score(coords, vep_annotation=_annotation(None))

        assert result.score is None
        assert result.source == "unavailable"
        client.close()


class TestNoConfiguredPath:
    async def test_none_data_path_skips_local_lookup_entirely(self) -> None:
        client = REVELClient(_settings(None))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        result = await client.get_score(coords, vep_annotation=_annotation(0.9))

        assert result.score == pytest.approx(0.9)
        assert result.source == "vep_rest"
        client.close()
