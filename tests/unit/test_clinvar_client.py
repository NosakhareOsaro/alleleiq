"""Unit tests for alleleiq.clients.clinvar.ClinVarClient.

All tests use httpx.MockTransport — no real network calls are made.
"""

from __future__ import annotations

from collections.abc import Awaitable, Callable
from datetime import date

import httpx
import pytest

from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.exceptions import ClinVarError

RequestHandler = Callable[[httpx.Request], Awaitable[httpx.Response]]


def _client_with_transport(handler: RequestHandler) -> ClinVarClient:
    client = ClinVarClient()
    client._client = httpx.AsyncClient(
        base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
        transport=httpx.MockTransport(handler),
    )
    return client


def _esearch_then_efetch(uid: str | None, xml_body: str) -> RequestHandler:
    async def handler(request: httpx.Request) -> httpx.Response:
        if "esearch" in str(request.url):
            idlist = [uid] if uid else []
            return httpx.Response(200, json={"esearchresult": {"idlist": idlist}})
        return httpx.Response(200, text=xml_body)

    return handler


class TestExpertPanelClassification:
    async def test_parses_three_star_expert_panel_record(
        self, clinvar_expert_panel_xml: str
    ) -> None:
        client = _client_with_transport(_esearch_then_efetch("17661", clinvar_expert_panel_xml))

        result = await client.fetch_classification(rsid="rs80357382")

        assert result is not None
        assert result.accession == "VCV000017661"
        assert result.clinical_significance == "Pathogenic"
        assert result.review_status == "reviewed by expert panel"
        assert result.star_rating == 3
        assert result.submitter_count == 3
        assert result.condition == "Hereditary breast and ovarian cancer syndrome"
        assert result.last_evaluated == date(2023, 6, 15)


class TestConflictingSubmissions:
    async def test_parses_conflicting_classification_as_one_star(
        self, clinvar_conflicting_xml: str
    ) -> None:
        client = _client_with_transport(_esearch_then_efetch("98765", clinvar_conflicting_xml))

        result = await client.fetch_classification(rsid="rs11571833")

        assert result is not None
        assert result.star_rating == 1
        assert "Conflicting" in result.clinical_significance
        assert result.submitter_count == 2

    async def test_one_star_result_falls_below_strong_evidence_threshold(
        self, clinvar_conflicting_xml: str
    ) -> None:
        from alleleiq.config import get_settings

        client = _client_with_transport(_esearch_then_efetch("98765", clinvar_conflicting_xml))
        result = await client.fetch_classification(rsid="rs11571833")

        assert result is not None
        assert result.star_rating < get_settings().clinvar_min_stars_strong_evidence


class TestNovelVariant:
    async def test_no_esearch_hits_returns_none(self) -> None:
        client = _client_with_transport(_esearch_then_efetch(None, ""))

        result = await client.fetch_classification(rsid="rs999999999")

        assert result is None


class TestSearchTermConstruction:
    async def test_missing_rsid_and_hgvs_raises_clinvar_error(self) -> None:
        client = _client_with_transport(_esearch_then_efetch(None, ""))

        with pytest.raises(ClinVarError, match="cannot search"):
            await client.fetch_classification()

    async def test_falls_back_to_gene_and_hgvs_when_no_rsid(
        self, clinvar_expert_panel_xml: str
    ) -> None:
        client = _client_with_transport(_esearch_then_efetch("17661", clinvar_expert_panel_xml))

        result = await client.fetch_classification(gene_symbol="BRCA1", hgvs_c="c.5096G>A")

        assert result is not None
        assert result.accession == "VCV000017661"


class TestErrorHandling:
    async def test_malformed_xml_raises_clinvar_error(self) -> None:
        client = _client_with_transport(_esearch_then_efetch("17661", "not xml <<<"))

        with pytest.raises(ClinVarError, match="malformed XML"):
            await client.fetch_classification(rsid="rs80357382")

    async def test_network_failure_raises_clinvar_error(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectTimeout("timed out")

        client = _client_with_transport(handler)
        client._max_retries = 0

        with pytest.raises(ClinVarError, match="ClinVar esearch failed"):
            await client.fetch_classification(rsid="rs80357382")


class TestHealthCheck:
    async def test_healthy_endpoint_returns_true(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = _client_with_transport(handler)
        assert await client.health_check() is True

    async def test_unreachable_endpoint_returns_false(self) -> None:
        async def handler(request: httpx.Request) -> httpx.Response:
            raise httpx.ConnectError("connection refused")

        client = _client_with_transport(handler)
        client._max_retries = 0
        assert await client.health_check() is False
