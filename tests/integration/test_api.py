"""Integration tests for the FastAPI REST API, via FastAPI's TestClient.

These tests exercise the full app (routing, middleware, pipeline
wiring) but mock the external evidence sources (gnomAD, VEP, ClinVar)
at the transport level, the same technique used in the Phase 2 client
unit tests. This keeps the checked-in, CI-run test suite fully offline
and deterministic; the clients themselves were separately verified
against the real live APIs during development (see the commit that
fixed the VEP/ClinVar field-shape bugs those live checks caught).
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import httpx
import pytest
from fastapi.testclient import TestClient

import alleleiq.pipeline as pipeline_module
from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.clients.vep import VEPClient
from alleleiq.config import get_settings


def _mocked_gnomad_client_factory(brca2_response: dict[str, Any]) -> type[GnomADClient]:
    class _FakeGnomADClient(GnomADClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FakeGnomADClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                body = request.content.decode()
                if "gene(" in body:
                    return httpx.Response(200, json={"data": {"gene": None}})
                return httpx.Response(200, json=brca2_response)

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    return _FakeGnomADClient


def _mocked_vep_client_factory(missense_response: list[dict[str, Any]]) -> type[VEPClient]:
    class _FakeVEPClient(VEPClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FakeVEPClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                return httpx.Response(200, json=missense_response)

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    return _FakeVEPClient


def _mocked_clinvar_client_factory() -> type[ClinVarClient]:
    class _FakeClinVarClient(ClinVarClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FakeClinVarClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                return httpx.Response(200, json={"esearchresult": {"idlist": []}})

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    return _FakeClinVarClient


@pytest.fixture
def mocked_evidence_sources(
    monkeypatch: pytest.MonkeyPatch,
    gnomad_brca2_response: dict[str, Any],
    vep_missense_response: list[dict[str, Any]],
) -> None:
    """Patch the evidence-source clients pipeline.py constructs, at the module level."""
    monkeypatch.setattr(
        pipeline_module, "GnomADClient", _mocked_gnomad_client_factory(gnomad_brca2_response)
    )
    monkeypatch.setattr(
        pipeline_module, "VEPClient", _mocked_vep_client_factory(vep_missense_response)
    )
    monkeypatch.setattr(pipeline_module, "ClinVarClient", _mocked_clinvar_client_factory())


@pytest.fixture
def client(mocked_evidence_sources: None) -> Iterator[TestClient]:
    get_settings.cache_clear()
    from alleleiq.api.main import create_app

    with TestClient(create_app()) as test_client:
        yield test_client
    get_settings.cache_clear()


def _single_snv_vcf_bytes() -> bytes:
    return (
        b"##fileformat=VCFv4.2\n"
        b"##reference=GRCh38\n"
        b"##contig=<ID=13,length=114364328>\n"
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        b"13\t32340301\t.\tT\tC\t100\tPASS\t.\n"
    )


def _multi_variant_vcf_bytes() -> bytes:
    return (
        b"##fileformat=VCFv4.2\n"
        b"##reference=GRCh38\n"
        b"##contig=<ID=13,length=114364328>\n"
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        b"13\t32340301\t.\tT\tC\t100\tPASS\t.\n"
        b"13\t32340350\t.\tG\tA\t100\tPASS\t.\n"
    )


def _grch37_vcf_bytes() -> bytes:
    return (
        b"##fileformat=VCFv4.2\n"
        b"##reference=GRCh37\n"
        b"##contig=<ID=13,length=115169878>\n"
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        b"13\t32340301\t.\tT\tC\t100\tPASS\t.\n"
    )


class TestAnalyseSingleVariant:
    def test_valid_single_variant_vcf_returns_evidence_report(self, client: TestClient) -> None:
        response = client.post(
            "/api/v1/variants/analyse",
            files={"vcf_file": ("variant.vcf", _single_snv_vcf_bytes(), "text/plain")},
        )

        assert response.status_code == 200
        body = response.json()
        assert "bayesian_score" in body
        assert "preliminary_classification" in body["bayesian_score"]
        assert "clinical_disclaimer" in body
        assert "CLINICAL SCIENTIST REVIEW ONLY" in body["clinical_disclaimer"]

    def test_malformed_vcf_returns_422(self, client: TestClient) -> None:
        response = client.post(
            "/api/v1/variants/analyse",
            files={"vcf_file": ("bad.vcf", b"this is not a vcf file", "text/plain")},
        )

        assert response.status_code == 422
        assert "detail" in response.json()

    def test_grch37_vcf_returns_422_with_clear_message(self, client: TestClient) -> None:
        response = client.post(
            "/api/v1/variants/analyse",
            files={"vcf_file": ("variant.vcf", _grch37_vcf_bytes(), "text/plain")},
        )

        assert response.status_code == 422
        assert "GRCh37" in response.json()["detail"]


class TestAnalyseBatch:
    def test_multi_variant_vcf_returns_one_report_per_variant(self, client: TestClient) -> None:
        response = client.post(
            "/api/v1/variants/analyse/batch",
            files={"vcf_file": ("variants.vcf", _multi_variant_vcf_bytes(), "text/plain")},
        )

        assert response.status_code == 200
        body = response.json()
        assert body["total_variants"] == 2
        assert body["successful"] == 2
        assert len(body["reports"]) == 2
        assert body["errors"] == []

    def test_malformed_batch_vcf_returns_422(self, client: TestClient) -> None:
        response = client.post(
            "/api/v1/variants/analyse/batch",
            files={"vcf_file": ("bad.vcf", b"not a vcf", "text/plain")},
        )

        assert response.status_code == 422


class TestHealthEndpoints:
    def test_health_returns_200_with_version(self, client: TestClient) -> None:
        response = client.get("/api/v1/health")

        assert response.status_code == 200
        body = response.json()
        assert body["status"] == "healthy"
        assert "version" in body

    def test_ready_returns_200_or_503_with_service_status(self, client: TestClient) -> None:
        response = client.get("/api/v1/health/ready")

        assert response.status_code in (200, 503)
        body = response.json()
        assert set(body.keys()) == {"status", "gnomad", "vep", "clinvar"}


class TestRateLimiting:
    def test_rate_limit_exceeded_returns_429(
        self, monkeypatch: pytest.MonkeyPatch, mocked_evidence_sources: None
    ) -> None:
        monkeypatch.setenv("ALLELEIQ_API_RATE_LIMIT_PER_HOUR", "3")
        get_settings.cache_clear()
        from alleleiq.api.main import create_app

        with TestClient(create_app()) as limited_client:
            statuses = [limited_client.get("/api/v1/health").status_code for _ in range(5)]

        assert statuses[:3] == [200, 200, 200]
        assert 429 in statuses[3:]

        monkeypatch.delenv("ALLELEIQ_API_RATE_LIMIT_PER_HOUR", raising=False)
        get_settings.cache_clear()
