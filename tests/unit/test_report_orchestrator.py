"""Unit tests for alleleiq.reports.ReportOrchestrator."""

from __future__ import annotations

from pathlib import Path

from alleleiq.models.evidence import BayesianScore
from alleleiq.models.report import AuditTrail, CoverageSummary, EvidenceReport, ReportMetadata
from alleleiq.models.variant import Variant
from alleleiq.reports import ReportOrchestrator


def _report(
    brca1_variant: Variant,
    report_metadata: ReportMetadata,
    complete_coverage_summary: CoverageSummary,
    vus_bayesian_score: BayesianScore,
    empty_audit_trail: AuditTrail,
) -> EvidenceReport:
    return EvidenceReport(
        variant=brca1_variant,
        metadata=report_metadata,
        coverage_summary=complete_coverage_summary,
        bayesian_score=vus_bayesian_score,
        audit_trail=empty_audit_trail,
    )


class TestGenerate:
    def test_generates_both_formats_by_default(
        self,
        tmp_path: Path,
        brca1_variant: Variant,
        report_metadata: ReportMetadata,
        complete_coverage_summary: CoverageSummary,
        vus_bayesian_score: BayesianScore,
        empty_audit_trail: AuditTrail,
    ) -> None:
        report = _report(
            brca1_variant,
            report_metadata,
            complete_coverage_summary,
            vus_bayesian_score,
            empty_audit_trail,
        )
        orchestrator = ReportOrchestrator(tmp_path)

        paths = orchestrator.generate(report)

        assert paths.json_path is not None
        assert paths.html_path is not None
        assert paths.json_path.exists()
        assert paths.html_path.exists()
        assert paths.json_path.stat().st_size > 0
        assert paths.html_path.stat().st_size > 0

    def test_filenames_use_variant_id(
        self,
        tmp_path: Path,
        brca1_variant: Variant,
        report_metadata: ReportMetadata,
        complete_coverage_summary: CoverageSummary,
        vus_bayesian_score: BayesianScore,
        empty_audit_trail: AuditTrail,
    ) -> None:
        report = _report(
            brca1_variant,
            report_metadata,
            complete_coverage_summary,
            vus_bayesian_score,
            empty_audit_trail,
        )
        orchestrator = ReportOrchestrator(tmp_path)

        paths = orchestrator.generate(report)

        assert paths.json_path is not None
        assert paths.json_path.name == f"{report.variant.variant_id}.json"

    def test_json_only_leaves_html_path_none(
        self,
        tmp_path: Path,
        brca1_variant: Variant,
        report_metadata: ReportMetadata,
        complete_coverage_summary: CoverageSummary,
        vus_bayesian_score: BayesianScore,
        empty_audit_trail: AuditTrail,
    ) -> None:
        report = _report(
            brca1_variant,
            report_metadata,
            complete_coverage_summary,
            vus_bayesian_score,
            empty_audit_trail,
        )
        orchestrator = ReportOrchestrator(tmp_path)

        paths = orchestrator.generate(report, formats=("json",))

        assert paths.json_path is not None
        assert paths.html_path is None

    def test_creates_output_directory_if_missing(
        self,
        tmp_path: Path,
        brca1_variant: Variant,
        report_metadata: ReportMetadata,
        complete_coverage_summary: CoverageSummary,
        vus_bayesian_score: BayesianScore,
        empty_audit_trail: AuditTrail,
    ) -> None:
        report = _report(
            brca1_variant,
            report_metadata,
            complete_coverage_summary,
            vus_bayesian_score,
            empty_audit_trail,
        )
        nested_dir = tmp_path / "does" / "not" / "exist" / "yet"
        orchestrator = ReportOrchestrator(nested_dir)

        paths = orchestrator.generate(report)

        assert nested_dir.is_dir()
        assert paths.json_path is not None and paths.json_path.exists()
