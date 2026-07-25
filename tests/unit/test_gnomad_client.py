"""Unit tests for alleleiq.clients.gnomad.GnomADClient.

All tests use httpx.MockTransport — no real network calls are made.
"""

from __future__ import annotations

from collections.abc import Awaitable, Callable
from typing import Any

import httpx
import pytest

from alleleiq.clients.gnomad import GnomADClient
from alleleiq.exceptions import GnomADError
from alleleiq.models.population import AncestryGroup
from alleleiq.models.variant import GenomicCoordinates

RequestHandler = Callable[[httpx.Request], Awaitable[httpx.Response]]


def _client_with_transport(handler: RequestHandler) -> GnomADClient:
    client = GnomADClient()
    client._client = httpx.AsyncClient(
        base_url="https://gnomad.broadinstitute.org/api",
        transport=httpx.MockTransport(handler),
    )
    return client


def _json_handler(payload: dict[str, Any], status_code: int = 200) -> RequestHandler:
    async def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(status_code, json=payload)

    return handler


class TestFetchPopulationFrequencies:
    async def test_merges_genome_and_exome_per_population(
        self, gnomad_brca2_response: dict[str, Any]
    ) -> None:
        client = _client_with_transport(_json_handler(gnomad_brca2_response))
        coords = GenomicCoordinates(chrom="13", pos=32340301, ref="T", alt="C")

        result = await client.fetch_population_frequencies(coords)

        afr = result[AncestryGroup.AFR]
        assert afr.allele_count == 8 + 10
        assert afr.allele_number == 41438 + 41398

    async def test_never_aggregates_across_populations(
        self, gnomad_brca2_response: dict[str, Any]
    ) -> None:
        """Each ancestry group's frequency must be independently derivable
        from the response — never summed into one overall value."""
        client = _client_with_transport(_json_handler(gnomad_brca2_response))
        coords = GenomicCoordinates(chrom="13", pos=32340301, ref="T", alt="C")

        result = await client.fetch_population_frequencies(coords)

        assert len(result) == len(AncestryGroup)
        assert (
            result[AncestryGroup.NFE].allele_frequency != result[AncestryGroup.EAS].allele_frequency
        )

    async def test_population_present_only_in_exome_still_has_data(
        self, gnomad_brca2_response: dict[str, Any]
    ) -> None:
        """SAS has an=0 in the genome cohort but real data in exome; the
        merged result must reflect exome coverage, not "no data"."""
        client = _client_with_transport(_json_handler(gnomad_brca2_response))
        coords = GenomicCoordinates(chrom="13", pos=32340301, ref="T", alt="C")

        result = await client.fetch_population_frequencies(coords)

        assert result[AncestryGroup.SAS].has_data is True
        assert result[AncestryGroup.SAS].allele_number == 15308

    async def test_variant_not_found_returns_empty_mapping(
        self, gnomad_no_data_response: dict[str, Any]
    ) -> None:
        client = _client_with_transport(_json_handler(gnomad_no_data_response))
        coords = GenomicCoordinates(chrom="1", pos=123456, ref="A", alt="G")

        result = await client.fetch_population_frequencies(coords)

        assert result == {}

    async def test_mitochondrial_variant_raises_gnomad_error(self) -> None:
        client = _client_with_transport(_json_handler({"data": {"variant": None}}))
        coords = GenomicCoordinates(chrom="MT", pos=100, ref="A", alt="G")

        with pytest.raises(GnomADError, match="mitochondrial"):
            await client.fetch_population_frequencies(coords)

    async def test_network_timeout_raised_as_gnomad_error(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectTimeout("connection timed out")

        client = _client_with_transport(handler)
        client._max_retries = 0
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="G")

        with pytest.raises(GnomADError, match="gnomAD request failed"):
            await client.fetch_population_frequencies(coords)

    async def test_graphql_error_response_raises_gnomad_error(self) -> None:
        client = _client_with_transport(
            _json_handler({"errors": [{"message": "Variant ID is not valid"}]})
        )
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="G")

        with pytest.raises(GnomADError, match="GraphQL error"):
            await client.fetch_population_frequencies(coords)

    async def test_malformed_json_response_raises_gnomad_error(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, content=b"not json")

        client = _client_with_transport(handler)
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="G")

        with pytest.raises(GnomADError, match="non-JSON"):
            await client.fetch_population_frequencies(coords)

    async def test_rate_limit_response_raised_as_gnomad_error(self) -> None:
        client = _client_with_transport(_json_handler({"error": "rate limited"}, status_code=429))
        client._max_retries = 0
        coords = GenomicCoordinates(chrom="1", pos=100, ref="A", alt="G")

        with pytest.raises(GnomADError, match="gnomAD request failed"):
            await client.fetch_population_frequencies(coords)


class TestFetchGeneConstraint:
    async def test_parses_pli_and_loeuf(self) -> None:
        payload = {
            "data": {
                "gene": {
                    "gene_id": "ENSG00000012048",
                    "symbol": "BRCA1",
                    "gnomad_constraint": {"pli": 0.02, "oe_lof": 0.85, "oe_lof_upper": 1.12},
                }
            }
        }
        client = _client_with_transport(_json_handler(payload))

        constraint = await client.fetch_gene_constraint("BRCA1")

        assert constraint is not None
        assert constraint.pli == 0.02
        assert constraint.loeuf == 1.12

    async def test_gene_not_found_returns_none(self) -> None:
        client = _client_with_transport(_json_handler({"data": {"gene": None}}))

        constraint = await client.fetch_gene_constraint("NOTAREALGENE")

        assert constraint is None


class TestHealthCheck:
    async def test_healthy_endpoint_returns_true(self) -> None:
        client = _client_with_transport(_json_handler({"data": {"__typename": "Query"}}))
        assert await client.health_check() is True

    async def test_unreachable_endpoint_returns_false(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("connection refused")

        client = _client_with_transport(handler)
        client._max_retries = 0
        assert await client.health_check() is False
