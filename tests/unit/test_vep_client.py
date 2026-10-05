"""Unit tests for alleleiq.clients.vep.VEPClient.

All tests use httpx.MockTransport — no real network calls are made.
"""

from __future__ import annotations

from collections.abc import Awaitable, Callable
from typing import Any

import httpx
import pytest

from alleleiq.clients.vep import VEPClient
from alleleiq.exceptions import VEPError
from alleleiq.models.variant import GenomicCoordinates

RequestHandler = Callable[[httpx.Request], Awaitable[httpx.Response]]


def _client_with_transport(handler: RequestHandler) -> VEPClient:
    client = VEPClient()
    client._client = httpx.AsyncClient(
        base_url="https://rest.ensembl.org", transport=httpx.MockTransport(handler)
    )
    return client


def _json_handler(payload: Any, status_code: int = 200) -> RequestHandler:
    async def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(status_code, json=payload)

    return handler


class TestMissenseVariant:
    async def test_parses_consequence_and_hgvs(
        self, vep_missense_response: list[dict[str, Any]]
    ) -> None:
        client = _client_with_transport(_json_handler(vep_missense_response))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        annotation = await client.fetch_annotation(coords)

        assert annotation is not None
        assert annotation.most_severe_consequence == "missense_variant"
        assert annotation.gene_symbol == "BRCA1"
        assert annotation.hgvs_p == "ENSP00000350283.3:p.Arg1699Gln"
        assert annotation.is_canonical is True
        assert annotation.existing_variant_id == "rs80357382"

    async def test_revel_score_parsed_when_available(
        self, vep_missense_response: list[dict[str, Any]]
    ) -> None:
        client = _client_with_transport(_json_handler(vep_missense_response))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        annotation = await client.fetch_annotation(coords)

        assert annotation is not None
        assert annotation.revel_score == pytest.approx(0.912)


class TestStopGainedVariant:
    async def test_classified_as_stop_gained(
        self, vep_stop_gained_response: list[dict[str, Any]]
    ) -> None:
        client = _client_with_transport(_json_handler(vep_stop_gained_response))
        coords = GenomicCoordinates(chrom="17", pos=43093157, ref="C", alt="A")

        annotation = await client.fetch_annotation(coords)

        assert annotation is not None
        assert annotation.most_severe_consequence == "stop_gained"
        assert annotation.hgvs_p == "ENSP00000350283.3:p.Tyr1737Ter"

    async def test_revel_score_absent_when_not_provided(
        self, vep_stop_gained_response: list[dict[str, Any]]
    ) -> None:
        client = _client_with_transport(_json_handler(vep_stop_gained_response))
        coords = GenomicCoordinates(chrom="17", pos=43093157, ref="C", alt="A")

        annotation = await client.fetch_annotation(coords)

        assert annotation is not None
        assert annotation.revel_score is None


class TestSpliceSiteVariant:
    async def test_spliceai_scores_parsed(
        self, vep_splice_site_response: list[dict[str, Any]]
    ) -> None:
        client = _client_with_transport(_json_handler(vep_splice_site_response))
        coords = GenomicCoordinates(chrom="17", pos=43091032, ref="T", alt="A")

        annotation = await client.fetch_annotation(coords)

        assert annotation is not None
        assert annotation.most_severe_consequence == "splice_donor_variant"
        assert annotation.spliceai_ds_dl == pytest.approx(0.91)


class TestNonCodingVariant:
    async def test_variant_with_no_transcript_consequences_degrades_gracefully(
        self, vep_non_coding_response: list[dict[str, Any]]
    ) -> None:
        client = _client_with_transport(_json_handler(vep_non_coding_response))
        coords = GenomicCoordinates(chrom="17", pos=43125270, ref="G", alt="A")

        annotation = await client.fetch_annotation(coords)

        assert annotation is not None
        assert annotation.most_severe_consequence == "intergenic_variant"
        assert annotation.gene_symbol is None
        assert annotation.transcript_id is None
        assert annotation.revel_score is None


class TestVariantOutsideAnnotatedRegions:
    async def test_empty_response_yields_no_annotation(self) -> None:
        client = _client_with_transport(_json_handler([]))
        coords = GenomicCoordinates(chrom="1", pos=1, ref="A", alt="G")

        annotation = await client.fetch_annotation(coords)

        assert annotation is None

    async def test_empty_coordinates_list_short_circuits(self) -> None:
        client = _client_with_transport(_json_handler([]))
        result = await client.fetch_annotations([])
        assert result == {}


class TestErrorHandling:
    async def test_vep_error_response_raises_vep_error(self) -> None:
        client = _client_with_transport(_json_handler({"error": "server overloaded"}))
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        with pytest.raises(VEPError, match="VEP error"):
            await client.fetch_annotation(coords)

    async def test_server_downtime_raises_vep_error(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("connection refused")

        client = _client_with_transport(handler)
        client._max_retries = 0
        coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")

        with pytest.raises(VEPError, match="VEP request failed"):
            await client.fetch_annotation(coords)

    async def test_health_check_false_on_downtime(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("connection refused")

        client = _client_with_transport(handler)
        client._max_retries = 0
        assert await client.health_check() is False


class TestFetchRelease:
    async def test_parses_and_caches_release(self) -> None:
        VEPClient._release = None
        calls = 0

        async def handler(request: httpx.Request) -> httpx.Response:
            nonlocal calls
            calls += 1
            assert request.url.path == "/info/software"
            assert request.headers["Content-Type"] == "application/json"
            return httpx.Response(200, json={"release": 112})

        try:
            client = _client_with_transport(handler)
            assert await client.fetch_release() == "112"
            assert await client.fetch_release() == "112"
            assert calls == 1
        finally:
            VEPClient._release = None

    async def test_failure_returns_none_without_raising(self) -> None:
        VEPClient._release = None

        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("connection refused")

        client = _client_with_transport(handler)
        client._max_retries = 0
        assert await client.fetch_release() is None
