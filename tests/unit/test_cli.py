"""Unit tests for alleleiq.cli, using Click's CliRunner.

Evidence-gathering (`build_evidence_report`) and validation
(`run_validation`) are monkeypatched at the point of use in
`alleleiq.cli` — this suite tests CLI behaviour (argument parsing,
output formatting, error handling, exit codes), not the pipeline itself
(covered elsewhere).
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

import alleleiq.cli as cli_module
from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.clients.vep import VEPClient
from alleleiq.exceptions import GnomADError
from alleleiq.models.evidence import (
    BayesianScore,
    ClassificationConfidence,
    PreliminaryClassification,
)
from alleleiq.models.population import OverallCoverageAssessment
from alleleiq.models.report import AuditTrail, CoverageSummary, EvidenceReport, ReportMetadata
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.validation import GoldStandardVariant, Outcome, ValidationResult


def _report(
    variant_id_coords: tuple[str, int, str, str] = ("17", 43094692, "G", "A"),
    classification: PreliminaryClassification = PreliminaryClassification.UNCERTAIN_SIGNIFICANCE,
) -> EvidenceReport:
    chrom, pos, ref, alt = variant_id_coords
    coords = GenomicCoordinates(chrom=chrom, pos=pos, ref=ref, alt=alt)
    return EvidenceReport(
        variant=Variant(coordinates=coords, variant_type=VariantType.SNV, gene_symbol="BRCA1"),
        metadata=ReportMetadata(alleleiq_version="0.1.0-alpha"),
        coverage_summary=CoverageSummary(
            overall_coverage_assessment=OverallCoverageAssessment.COMPLETE,
            clinical_advisory="ok",
        ),
        bayesian_score=BayesianScore(
            evidence_points=0,
            preliminary_classification=classification,
            classification_confidence=ClassificationConfidence.MODERATE,
        ),
        audit_trail=AuditTrail(),
    )


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


class TestAnalyse:
    def test_successful_analysis_writes_reports_and_prints_summary(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        single_snv_vcf_path: Path,
    ) -> None:
        async def fake_build(
            variant: Variant, *, settings: object, include_clinvar: bool
        ) -> EvidenceReport:
            return _report()

        monkeypatch.setattr(cli_module, "build_evidence_report", fake_build)
        monkeypatch.setattr(
            GnomADClient,
            "cache_stats",
            classmethod(lambda cls: {"hits": 0, "misses": 0, "size": 0}),
        )

        result = runner.invoke(
            cli_module.cli,
            ["analyse", str(single_snv_vcf_path), "--output-dir", str(tmp_path)],
        )

        assert result.exit_code == 0, result.output
        assert "AlleleIQ analysis summary" in result.output
        assert "PRELIMINARY EVIDENCE SUMMARY" in result.output
        assert (tmp_path / "17-43094692-G-A.json").exists()
        assert (tmp_path / "17-43094692-G-A.html").exists()

    def test_malformed_vcf_exits_nonzero_with_suggestion(
        self, runner: CliRunner, malformed_vcf_path: Path, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            cli_module.cli,
            ["analyse", str(malformed_vcf_path), "--output-dir", str(tmp_path)],
        )

        assert result.exit_code == 1
        assert "Error" in result.output
        assert "Suggestion" in result.output

    def test_evidence_source_failure_is_reported_but_does_not_crash(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        single_snv_vcf_path: Path,
    ) -> None:
        async def failing_build(
            variant: Variant, *, settings: object, include_clinvar: bool
        ) -> EvidenceReport:
            raise GnomADError("simulated upstream failure")

        monkeypatch.setattr(cli_module, "build_evidence_report", failing_build)
        monkeypatch.setattr(
            GnomADClient,
            "cache_stats",
            classmethod(lambda cls: {"hits": 0, "misses": 0, "size": 0}),
        )

        result = runner.invoke(
            cli_module.cli,
            ["analyse", str(single_snv_vcf_path), "--output-dir", str(tmp_path)],
        )

        assert result.exit_code == 0
        assert "Error analysing" in result.output
        assert "0 report(s) written" in result.output

    def test_json_only_format_skips_html(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        single_snv_vcf_path: Path,
    ) -> None:
        async def fake_build(
            variant: Variant, *, settings: object, include_clinvar: bool
        ) -> EvidenceReport:
            return _report()

        monkeypatch.setattr(cli_module, "build_evidence_report", fake_build)
        monkeypatch.setattr(
            GnomADClient,
            "cache_stats",
            classmethod(lambda cls: {"hits": 0, "misses": 0, "size": 0}),
        )

        result = runner.invoke(
            cli_module.cli,
            [
                "analyse",
                str(single_snv_vcf_path),
                "--output-dir",
                str(tmp_path),
                "--format",
                "json",
            ],
        )

        assert result.exit_code == 0
        assert (tmp_path / "17-43094692-G-A.json").exists()
        assert not (tmp_path / "17-43094692-G-A.html").exists()


class TestValidate:
    def test_all_concordant_exits_zero(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        gold = GoldStandardVariant(
            label="test",
            chrom="17",
            pos=1,
            ref="A",
            alt="G",
            gene_symbol="X",
            clinvar_accession="VCV1",
            expected_classification=PreliminaryClassification.PATHOGENIC,
        )
        results = [
            ValidationResult(
                gold_standard=gold,
                actual_classification=PreliminaryClassification.PATHOGENIC,
                evidence_points=10,
                criteria_applied=["PVS1"],
                outcome=Outcome.CONCORDANT,
                explanation=None,
            )
        ]

        async def fake_run_validation(*args: object, **kwargs: object) -> list[ValidationResult]:
            return results

        monkeypatch.setattr(cli_module, "run_validation", fake_run_validation)

        result = runner.invoke(cli_module.cli, ["validate"])

        assert result.exit_code == 0

    def test_dangerous_outcome_exits_nonzero(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        gold = GoldStandardVariant(
            label="test",
            chrom="17",
            pos=1,
            ref="A",
            alt="G",
            gene_symbol="X",
            clinvar_accession="VCV1",
            expected_classification=PreliminaryClassification.PATHOGENIC,
        )
        results = [
            ValidationResult(
                gold_standard=gold,
                actual_classification=PreliminaryClassification.BENIGN,
                evidence_points=-8,
                criteria_applied=["BA1"],
                outcome=Outcome.DANGEROUS,
                explanation="expected PATHOGENIC but got BENIGN",
            )
        ]

        async def fake_run_validation(*args: object, **kwargs: object) -> list[ValidationResult]:
            return results

        monkeypatch.setattr(cli_module, "run_validation", fake_run_validation)

        result = runner.invoke(cli_module.cli, ["validate"])

        assert result.exit_code == 1


class TestVersion:
    def test_prints_version_and_connectivity(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        async def healthy(self: object) -> bool:
            return True

        monkeypatch.setattr(GnomADClient, "health_check", healthy)
        monkeypatch.setattr(VEPClient, "health_check", healthy)
        monkeypatch.setattr(ClinVarClient, "health_check", healthy)

        result = runner.invoke(cli_module.cli, ["version"])

        assert result.exit_code == 0
        assert "AlleleIQ" in result.output
        assert "gnomAD dataset" in result.output
        assert "reachable" in result.output

    def test_reports_unreachable_services(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        async def unhealthy(self: object) -> bool:
            return False

        monkeypatch.setattr(GnomADClient, "health_check", unhealthy)
        monkeypatch.setattr(VEPClient, "health_check", unhealthy)
        monkeypatch.setattr(ClinVarClient, "health_check", unhealthy)

        result = runner.invoke(cli_module.cli, ["version"])

        assert result.exit_code == 0
        assert "unreachable" in result.output


class TestCLIGroup:
    def test_version_option(self, runner: CliRunner) -> None:
        result = runner.invoke(cli_module.cli, ["--version"])

        assert result.exit_code == 0
        assert "alleleiq" in result.output.lower()

    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(cli_module.cli, ["--help"])

        assert result.exit_code == 0
        assert "analyse" in result.output
        assert "validate" in result.output
        assert "version" in result.output
