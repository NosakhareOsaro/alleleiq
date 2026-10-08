"""Integration tests for the minimal server-rendered web interface.

Mirrors tests/integration/test_api.py's approach: exercises the full
app (routing, templates, pipeline wiring) via FastAPI's TestClient,
with external evidence sources mocked at the transport level so the
suite stays offline and deterministic.
"""

from __future__ import annotations

import asyncio
import os
import re
import time
import uuid
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import httpx
import pytest
from fastapi.testclient import TestClient

import alleleiq.pipeline as pipeline_module
from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.clients.uniprot import UniProtClient
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


def _mocked_uniprot_client_factory() -> type[UniProtClient]:
    class _FakeUniProtClient(UniProtClient):  # type: ignore[misc]
        async def __aenter__(self) -> _FakeUniProtClient:
            async def handler(request: httpx.Request) -> httpx.Response:
                return httpx.Response(200, json={"results": []})

            self._client = httpx.AsyncClient(
                base_url=self._base_url, transport=httpx.MockTransport(handler)
            )
            return self

    return _FakeUniProtClient


def _healthy_client_factory(base: type[Any]) -> type[Any]:
    """Wrap a client class so health_check() reports reachable without a real request."""

    class _HealthyClient(base):  # type: ignore[misc,valid-type]
        async def health_check(self) -> bool:
            return True

    return _HealthyClient


@pytest.fixture
def mocked_evidence_sources(
    monkeypatch: pytest.MonkeyPatch,
    gnomad_brca2_response: dict[str, Any],
    vep_missense_response: list[dict[str, Any]],
) -> None:
    """Patch the evidence-source clients pipeline.py and web/router.py construct."""
    gnomad_cls = _healthy_client_factory(_mocked_gnomad_client_factory(gnomad_brca2_response))
    vep_cls = _healthy_client_factory(_mocked_vep_client_factory(vep_missense_response))
    clinvar_cls = _healthy_client_factory(_mocked_clinvar_client_factory())
    uniprot_cls = _healthy_client_factory(_mocked_uniprot_client_factory())

    monkeypatch.setattr(pipeline_module, "GnomADClient", gnomad_cls)
    monkeypatch.setattr(pipeline_module, "VEPClient", vep_cls)
    monkeypatch.setattr(pipeline_module, "ClinVarClient", clinvar_cls)
    monkeypatch.setattr(pipeline_module, "UniProtClient", uniprot_cls)

    import alleleiq.web.router as web_router_module

    monkeypatch.setattr(web_router_module, "GnomADClient", gnomad_cls)
    monkeypatch.setattr(web_router_module, "VEPClient", vep_cls)
    monkeypatch.setattr(web_router_module, "ClinVarClient", clinvar_cls)
    monkeypatch.setattr(web_router_module, "UniProtClient", uniprot_cls)


@pytest.fixture
def client(
    mocked_evidence_sources: None, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> Iterator[TestClient]:
    monkeypatch.setenv("ALLELEIQ_WEB_REPORTS_DIR", str(tmp_path / "web-reports"))
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


class TestIndexPage:
    def test_index_renders_with_connectivity_status(self, client: TestClient) -> None:
        response = client.get("/")

        assert response.status_code == 200
        assert "CLINICAL DECISION SUPPORT ONLY" in response.text
        assert "reachable" in response.text
        assert "UniProt:" in response.text
        assert "Analyse" in response.text

    def test_index_page_has_no_javascript_requirement_on_core_form(
        self, client: TestClient
    ) -> None:
        response = client.get("/")

        assert 'method="post"' in response.text
        assert 'enctype="multipart/form-data"' in response.text
        assert 'action="/analyse"' in response.text


class TestWebAnalyse:
    def test_plain_post_returns_full_result_page(self, client: TestClient) -> None:
        response = client.post(
            "/analyse", files={"vcf_file": ("variant.vcf", _single_snv_vcf_bytes(), "text/plain")}
        )

        assert response.status_code == 200
        assert "<!doctype html>" in response.text.lower()
        assert "CLINICAL DECISION SUPPORT ONLY" in response.text
        assert "Preliminary classification" in response.text
        assert "report-frame" in response.text
        assert "Download JSON" in response.text
        assert "Download HTML report" in response.text

    def test_htmx_post_returns_fragment_without_full_page_chrome(self, client: TestClient) -> None:
        response = client.post(
            "/analyse",
            files={"vcf_file": ("variant.vcf", _single_snv_vcf_bytes(), "text/plain")},
            headers={"HX-Request": "true"},
        )

        assert response.status_code == 200
        assert "<!doctype html>" not in response.text.lower()
        assert re.search(r"<head[ >]", response.text.lower()) is None
        assert "Preliminary classification" in response.text

    def test_plain_and_htmx_responses_agree_on_the_result(self, client: TestClient) -> None:
        full = client.post(
            "/analyse", files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")}
        )
        fragment = client.post(
            "/analyse",
            files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")},
            headers={"HX-Request": "true"},
        )

        classification_pattern = re.compile(r"Preliminary classification:</strong>\s*(\w+)")
        full_match = classification_pattern.search(full.text)
        fragment_match = classification_pattern.search(fragment.text)
        assert full_match is not None
        assert fragment_match is not None
        assert full_match.group(1) == fragment_match.group(1)

    def test_multi_variant_vcf_renders_one_card_per_variant(self, client: TestClient) -> None:
        response = client.post(
            "/analyse",
            files={"vcf_file": ("variants.vcf", _multi_variant_vcf_bytes(), "text/plain")},
        )

        assert response.status_code == 200
        assert response.text.count("report-card") == 2
        assert "2 variants analysed" in response.text

    def test_duplicate_variant_is_analysed_once(self, client: TestClient) -> None:
        vcf = _single_snv_vcf_bytes() + b"13\t32340301\t.\tT\tC\t100\tPASS\t.\n"

        response = client.post("/analyse", files={"vcf_file": ("dup.vcf", vcf, "text/plain")})

        assert response.status_code == 200
        assert response.text.count("report-card") == 1

    def test_malformed_vcf_renders_clinical_grade_error(self, client: TestClient) -> None:
        response = client.post(
            "/analyse", files={"vcf_file": ("bad.vcf", b"this is not a vcf file", "text/plain")}
        )

        assert response.status_code == 422
        assert "could not be" in response.text.lower()
        assert "Traceback" not in response.text
        assert "error-box" in response.text

    def test_grch37_vcf_renders_clear_message_not_traceback(self, client: TestClient) -> None:
        grch37_vcf = (
            b"##fileformat=VCFv4.2\n"
            b"##reference=GRCh37\n"
            b"##contig=<ID=13,length=115169878>\n"
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            b"13\t32340301\t.\tT\tC\t100\tPASS\t.\n"
        )
        response = client.post("/analyse", files={"vcf_file": ("v.vcf", grch37_vcf, "text/plain")})

        assert response.status_code == 422
        assert "GRCh37" in response.text
        assert "Traceback" not in response.text

    def test_oversized_upload_renders_clinical_grade_error(
        self, client: TestClient, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("ALLELEIQ_API_MAX_UPLOAD_BYTES", "10")
        get_settings.cache_clear()
        response = client.post(
            "/analyse", files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")}
        )

        assert response.status_code == 413
        assert "error-box" in response.text
        monkeypatch.delenv("ALLELEIQ_API_MAX_UPLOAD_BYTES", raising=False)
        get_settings.cache_clear()


class TestClinVarToggle:
    """An unchecked checkbox is omitted from a form submission, so leaving the
    ClinVar box unchecked must actually disable ClinVar."""

    @pytest.fixture
    def recorded_flags(self, monkeypatch: pytest.MonkeyPatch) -> list[bool]:
        import alleleiq.web.router as web_router_module

        flags: list[bool] = []
        real_build = web_router_module.build_evidence_report

        async def spy(*args: Any, **kwargs: Any) -> Any:
            flags.append(kwargs["include_clinvar"])
            return await real_build(*args, **kwargs)

        monkeypatch.setattr(web_router_module, "build_evidence_report", spy)
        return flags

    def test_unchecked_box_disables_clinvar(
        self, client: TestClient, recorded_flags: list[bool]
    ) -> None:
        client.post(
            "/analyse", files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")}
        )

        assert recorded_flags == [False]

    def test_checked_box_enables_clinvar(
        self, client: TestClient, recorded_flags: list[bool]
    ) -> None:
        client.post(
            "/analyse",
            data={"include_clinvar": "true"},
            files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")},
        )

        assert recorded_flags == [True]

    def test_index_renders_the_box_checked(self, client: TestClient) -> None:
        response = client.get("/")

        assert 'name="include_clinvar" value="true" checked' in response.text


class TestDownloadReport:
    def test_download_json_and_html_reports(self, client: TestClient) -> None:
        analyse_response = client.post(
            "/analyse", files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")}
        )
        job_id_match = re.search(r"/reports/([0-9a-f-]{36})/", analyse_response.text)
        assert job_id_match is not None
        job_id = job_id_match.group(1)

        json_response = client.get(f"/reports/{job_id}/13-32340301-T-C.json")
        assert json_response.status_code == 200
        assert json_response.headers["content-type"].startswith("application/json")

        html_response = client.get(f"/reports/{job_id}/13-32340301-T-C.html")
        assert html_response.status_code == 200

    def test_invalid_job_id_returns_404(self, client: TestClient) -> None:
        response = client.get("/reports/not-a-uuid/13-32340301-T-C.json")
        assert response.status_code == 404

    def test_path_traversal_filename_returns_404(self, client: TestClient) -> None:
        job_id = str(uuid.uuid4())
        response = client.get(f"/reports/{job_id}/../../../../etc/passwd")
        assert response.status_code == 404

    def test_nonexistent_report_file_returns_404(self, client: TestClient) -> None:
        job_id = str(uuid.uuid4())
        response = client.get(f"/reports/{job_id}/does-not-exist.json")
        assert response.status_code == 404

    def test_disallowed_extension_returns_404(self, client: TestClient) -> None:
        job_id = str(uuid.uuid4())
        response = client.get(f"/reports/{job_id}/config.py")
        assert response.status_code == 404


class TestProgressSSE:
    def test_progress_stream_reflects_real_pipeline_stages(self, client: TestClient) -> None:
        job_id = str(uuid.uuid4())
        analyse_response = client.post(
            "/analyse",
            data={"job_id": job_id},
            files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")},
        )
        assert analyse_response.status_code == 200

        with client.stream("GET", f"/analyse/{job_id}/progress") as response:
            assert response.status_code == 200
            body = "".join(response.iter_text())

        assert "variant 1/1: querying gnomAD" in body
        assert "variant 1/1: querying VEP" in body
        assert "variant 1/1: scoring" in body
        assert body.rstrip().endswith("event: complete\ndata: complete")

    def test_progress_stream_orders_stages_per_variant_in_multi_variant_batch(
        self, client: TestClient
    ) -> None:
        job_id = str(uuid.uuid4())
        client.post(
            "/analyse",
            data={"job_id": job_id},
            files={"vcf_file": ("v.vcf", _multi_variant_vcf_bytes(), "text/plain")},
        )

        with client.stream("GET", f"/analyse/{job_id}/progress") as response:
            body = "".join(response.iter_text())

        first_variant_index = body.index("variant 1/2: querying gnomAD")
        second_variant_index = body.index("variant 2/2: querying gnomAD")
        complete_index = body.index("event: complete")
        assert first_variant_index < second_variant_index < complete_index

    def test_progress_stream_completes_even_when_analysis_fails_to_parse(
        self, client: TestClient
    ) -> None:
        job_id = str(uuid.uuid4())
        client.post(
            "/analyse",
            data={"job_id": job_id},
            files={"vcf_file": ("bad.vcf", b"not a vcf", "text/plain")},
        )

        with client.stream("GET", f"/analyse/{job_id}/progress") as response:
            body = "".join(response.iter_text())

        assert body.strip() == "event: complete\ndata: complete"

    def test_progress_stream_with_invalid_job_id_returns_404(self, client: TestClient) -> None:
        response = client.get("/analyse/not-a-uuid/progress")
        assert response.status_code == 404

    def test_analyse_without_job_id_needs_no_progress_subscription(
        self, client: TestClient
    ) -> None:
        response = client.post(
            "/analyse", files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")}
        )
        assert response.status_code == 200


class TestProgressQueueLifecycle:
    async def test_unsubscribed_queue_is_released_after_the_grace_period(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import alleleiq.web.router as web_router_module

        monkeypatch.setattr(web_router_module, "_UNSUBSCRIBED_QUEUE_TTL_SECONDS", 0.01)
        job_id = str(uuid.uuid4())
        web_router_module._get_or_create_progress_queue(job_id)

        web_router_module._release_progress_queue_later(job_id)
        await asyncio.sleep(0.05)

        assert job_id not in web_router_module._progress_queues

    async def test_idle_subscription_ends_and_releases_its_queue(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import alleleiq.web.router as web_router_module

        monkeypatch.setattr(web_router_module, "_PROGRESS_IDLE_TIMEOUT_SECONDS", 0.01)
        job_id = str(uuid.uuid4())

        response = await web_router_module.analyse_progress(job_id)
        events = [chunk async for chunk in response.body_iterator]

        assert events == []
        assert job_id not in web_router_module._progress_queues


class TestReportRetention:
    @staticmethod
    def _make_job_dir(root: Path, name: str, age_hours: float) -> Path:
        job_dir = root / name
        job_dir.mkdir(parents=True)
        (job_dir / "report.json").write_text("{}")
        stamp = time.time() - age_hours * 3600
        os.utime(job_dir, (stamp, stamp))
        return job_dir

    def test_expired_job_directories_are_removed(self, tmp_path: Path) -> None:
        from alleleiq.web.router import _purge_expired_job_dirs

        old_job = self._make_job_dir(tmp_path, str(uuid.uuid4()), age_hours=30)
        recent_job = self._make_job_dir(tmp_path, str(uuid.uuid4()), age_hours=1)
        unrelated = self._make_job_dir(tmp_path, "keep-me", age_hours=30)

        _purge_expired_job_dirs(tmp_path, retention_hours=24)

        assert not old_job.exists()
        assert recent_job.exists()
        assert unrelated.exists()

    def test_zero_retention_keeps_everything(self, tmp_path: Path) -> None:
        from alleleiq.web.router import _purge_expired_job_dirs

        old_job = self._make_job_dir(tmp_path, str(uuid.uuid4()), age_hours=1000)

        _purge_expired_job_dirs(tmp_path, retention_hours=0)

        assert old_job.exists()

    def test_new_analysis_purges_expired_jobs(self, client: TestClient, tmp_path: Path) -> None:
        old_job = self._make_job_dir(tmp_path / "web-reports", str(uuid.uuid4()), age_hours=48)

        response = client.post(
            "/analyse", files={"vcf_file": ("v.vcf", _single_snv_vcf_bytes(), "text/plain")}
        )

        assert response.status_code == 200
        assert not old_job.exists()
