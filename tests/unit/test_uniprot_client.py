"""Unit tests for alleleiq.clients.uniprot.UniProtClient.

All tests use httpx.MockTransport — no real network calls are made.
"""

from __future__ import annotations

from collections.abc import Awaitable, Callable
from typing import Any

import httpx
import pytest

from alleleiq.clients.uniprot import UniProtClient
from alleleiq.exceptions import UniProtError

RequestHandler = Callable[[httpx.Request], Awaitable[httpx.Response]]


def _client_with_transport(handler: RequestHandler) -> UniProtClient:
    client = UniProtClient()
    client._client = httpx.AsyncClient(
        base_url="https://rest.uniprot.org", transport=httpx.MockTransport(handler)
    )
    return client


def _json_handler(payload: Any, status_code: int = 200) -> RequestHandler:
    async def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(status_code, json=payload)

    return handler


_BRCA1_RESPONSE = {
    "results": [
        {
            "primaryAccession": "P38398",
            "genes": [{"geneName": {"value": "BRCA1"}}],
            "features": [
                {
                    "type": "Domain",
                    "location": {"start": {"value": 1642}, "end": {"value": 1736}},
                    "description": "BRCT 1",
                },
                {
                    "type": "Domain",
                    "location": {"start": {"value": 1756}, "end": {"value": 1855}},
                    "description": "BRCT 2",
                },
                {
                    "type": "Region",
                    "location": {"start": {"value": 1000}, "end": {"value": 1010}},
                    "description": "Disordered - should be excluded, not a critical type",
                },
            ],
        }
    ]
}


class TestFetchDomains:
    async def test_parses_domain_regions_from_real_shape_response(self) -> None:
        UniProtClient.clear_cache()
        client = _client_with_transport(_json_handler(_BRCA1_RESPONSE))

        regions = await client.fetch_domains("BRCA1")

        assert len(regions) == 2
        assert regions[0].start_aa == 1642
        assert regions[0].end_aa == 1736
        assert "BRCT 1" in regions[0].description
        assert "P38398" in regions[0].source

    async def test_region_feature_type_is_excluded(self) -> None:
        UniProtClient.clear_cache()
        client = _client_with_transport(_json_handler(_BRCA1_RESPONSE))

        regions = await client.fetch_domains("BRCA1")

        assert all("Disordered" not in r.description for r in regions)

    async def test_no_reviewed_entry_returns_empty_list(self) -> None:
        UniProtClient.clear_cache()
        client = _client_with_transport(_json_handler({"results": []}))

        regions = await client.fetch_domains("NOTAREALGENE")

        assert regions == []

    async def test_result_is_cached_across_calls(self) -> None:
        UniProtClient.clear_cache()
        call_count = 0

        async def handler(request: httpx.Request) -> httpx.Response:
            nonlocal call_count
            call_count += 1
            return httpx.Response(200, json=_BRCA1_RESPONSE)

        client = _client_with_transport(handler)

        await client.fetch_domains("BRCA1")
        await client.fetch_domains("brca1")  # case-insensitive cache key

        assert call_count == 1

    async def test_transport_error_raises_uniprot_error(self) -> None:
        UniProtClient.clear_cache()

        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("simulated network failure")

        client = _client_with_transport(handler)

        with pytest.raises(UniProtError):
            await client.fetch_domains("BRCA1")

    async def test_non_json_response_raises_uniprot_error(self) -> None:
        UniProtClient.clear_cache()

        async def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, text="not json")

        client = _client_with_transport(handler)

        with pytest.raises(UniProtError):
            await client.fetch_domains("BRCA1")


class TestHealthCheck:
    async def test_returns_true_on_200(self) -> None:
        client = _client_with_transport(_json_handler({"results": []}))

        assert await client.health_check() is True

    async def test_returns_false_on_transport_error(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("simulated network failure")

        client = _client_with_transport(handler)

        assert await client.health_check() is False
