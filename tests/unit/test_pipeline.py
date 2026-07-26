"""Unit tests for alleleiq.pipeline.build_evidence_report.

Mocks GnomADClient/VEPClient/ClinVarClient at the transport level (as in
the Phase 2 client tests) so this exercises the real orchestration logic
— evidence gathering, audit trail recording (including failures),
variant enrichment from VEP, and criterion evaluation — without a
network dependency.
"""

from __future__ import annotations

from typing import Any

import httpx
import pytest

import alleleiq.pipeline as pipeline_module
from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.clients.vep import VEPClient
from alleleiq.models.evidence import EvidenceStrength
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.pipeline import build_evidence_report


def _patch_gnomad(monkeypatch: pytest.MonkeyPatch, response: dict[str, Any]) -> None:
    class _FakeGnomADClient(GnomADClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FakeGnomADClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                body = request.content.decode()
                if "gene(" in body:
                    return httpx.Response(200, json={"data": {"gene": None}})
                return httpx.Response(200, json=response)

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    monkeypatch.setattr(pipeline_module, "GnomADClient", _FakeGnomADClient)


def _patch_gnomad_failing(monkeypatch: pytest.MonkeyPatch) -> None:
    class _FailingGnomADClient(GnomADClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FailingGnomADClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                raise httpx.ConnectError("simulated network failure")

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    monkeypatch.setattr(pipeline_module, "GnomADClient", _FailingGnomADClient)
    monkeypatch.setattr(GnomADClient, "clear_cache", classmethod(lambda cls: None))


def _patch_vep(monkeypatch: pytest.MonkeyPatch, response: list[dict[str, Any]]) -> None:
    class _FakeVEPClient(VEPClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FakeVEPClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                return httpx.Response(200, json=response)

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    monkeypatch.setattr(pipeline_module, "VEPClient", _FakeVEPClient)


def _patch_clinvar_not_found(monkeypatch: pytest.MonkeyPatch) -> None:
    class _FakeClinVarClient(ClinVarClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FakeClinVarClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                return httpx.Response(200, json={"esearchresult": {"idlist": []}})

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    monkeypatch.setattr(pipeline_module, "ClinVarClient", _FakeClinVarClient)


def _variant() -> Variant:
    # Coordinates match vep_missense_response's "input" field (17 43094692 . G A . . .):
    # VEPClient looks results up by the variant_id parsed from the response body, so the
    # queried coordinates and the fixture's declared input must agree.
    coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")
    return Variant(coordinates=coords, variant_type=VariantType.SNV)


def _stop_gained_variant() -> Variant:
    # Coordinates match vep_stop_gained_response's "input" field (17 43093157 . C A . . .).
    coords = GenomicCoordinates(chrom="17", pos=43093157, ref="C", alt="A")
    return Variant(coordinates=coords, variant_type=VariantType.SNV)


class TestBuildEvidenceReport:
    async def test_successful_run_populates_report_and_audit_trail(
        self,
        monkeypatch: pytest.MonkeyPatch,
        gnomad_brca2_response: dict[str, Any],
        vep_missense_response: list[dict[str, Any]],
    ) -> None:
        _patch_gnomad(monkeypatch, gnomad_brca2_response)
        _patch_vep(monkeypatch, vep_missense_response)
        _patch_clinvar_not_found(monkeypatch)

        report = await build_evidence_report(_variant())

        assert report.population_frequencies != {}
        assert report.variant.gene_symbol == "BRCA1"
        assert report.variant.consequence == "missense_variant"
        assert any(q.database == "gnomAD" and q.success for q in report.audit_trail.queries)
        assert any(q.database == "VEP" and q.success for q in report.audit_trail.queries)
        assert len(report.criterion_evaluations) == 8

    async def test_skips_clinvar_when_disabled(
        self,
        monkeypatch: pytest.MonkeyPatch,
        gnomad_brca2_response: dict[str, Any],
        vep_missense_response: list[dict[str, Any]],
    ) -> None:
        _patch_gnomad(monkeypatch, gnomad_brca2_response)
        _patch_vep(monkeypatch, vep_missense_response)

        report = await build_evidence_report(_variant(), include_clinvar=False)

        assert not any(q.database == "ClinVar" for q in report.audit_trail.queries)

    async def test_gnomad_failure_is_recorded_not_raised(
        self,
        monkeypatch: pytest.MonkeyPatch,
        vep_missense_response: list[dict[str, Any]],
    ) -> None:
        _patch_gnomad_failing(monkeypatch)
        _patch_vep(monkeypatch, vep_missense_response)
        _patch_clinvar_not_found(monkeypatch)

        report = await build_evidence_report(_variant())

        gnomad_queries = [q for q in report.audit_trail.queries if q.database == "gnomAD"]
        assert any(not q.success for q in gnomad_queries)
        assert report.population_frequencies == {}
        # Coverage assessment must still run against the (empty) frequencies —
        # a failed query must never look identical to "no criteria assessed at all".
        assert report.coverage_summary is not None


class TestClinGenIntegration:
    """ClinGenClient is consulted for every variant with a resolved gene symbol and
    recorded in the audit trail; its result flows into PVS1's mechanism tiering.

    Uses the repository's real bundled ClinGen snapshot (data/clingen/, the
    default settings path) rather than a mock, since BRCA1's real ClinGen
    classification (Definitive, haploinsufficiency score 3) is a stable,
    well-established curation unlikely to regress to a lower tier — this
    doubles as a smoke test that the bundled data still supports BRCA1 at
    full strength.
    """

    async def test_clingen_query_recorded_in_audit_trail(
        self,
        monkeypatch: pytest.MonkeyPatch,
        gnomad_brca2_response: dict[str, Any],
        vep_stop_gained_response: list[dict[str, Any]],
    ) -> None:
        _patch_gnomad(monkeypatch, gnomad_brca2_response)
        _patch_vep(monkeypatch, vep_stop_gained_response)
        _patch_clinvar_not_found(monkeypatch)

        report = await build_evidence_report(_stop_gained_variant())

        clingen_queries = [q for q in report.audit_trail.queries if q.database == "ClinGen"]
        assert len(clingen_queries) == 1
        assert clingen_queries[0].success is True
        assert clingen_queries[0].parameters["gene_symbol"] == "BRCA1"

    async def test_clingen_classification_flows_into_pvs1_strength(
        self,
        monkeypatch: pytest.MonkeyPatch,
        gnomad_brca2_response: dict[str, Any],
        vep_stop_gained_response: list[dict[str, Any]],
    ) -> None:
        _patch_gnomad(monkeypatch, gnomad_brca2_response)
        _patch_vep(monkeypatch, vep_stop_gained_response)
        _patch_clinvar_not_found(monkeypatch)

        report = await build_evidence_report(_stop_gained_variant())

        pvs1 = next(e for e in report.criterion_evaluations if e.code.value == "PVS1")
        assert pvs1.met is True
        assert pvs1.strength == EvidenceStrength.VERY_STRONG
        assert pvs1.requires_clinical_input is False
        assert "Definitive" in pvs1.evidence_basis
