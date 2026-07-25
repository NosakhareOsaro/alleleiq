"""Tests for JSON and HTML evidence report generation, including a snapshot test."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from uuid import UUID

from alleleiq.models import (
    ACMGCode,
    AncestryGroup,
    AuditTrail,
    BayesianScore,
    ClassificationConfidence,
    CoverageFlagLevel,
    CoverageGap,
    CoverageSummary,
    CriterionEvaluation,
    EvidenceStrength,
    GenomicCoordinates,
    OverallCoverageAssessment,
    PopulationFrequency,
    PreliminaryClassification,
    ReportMetadata,
    Variant,
    VariantType,
)
from alleleiq.models.report import CLINICAL_DISCLAIMER, EvidenceReport
from alleleiq.reports.html_report import render_html_report
from alleleiq.reports.json_report import render_json_report


def _brca1_variant() -> Variant:
    coords = GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")
    return Variant(
        coordinates=coords,
        variant_type=VariantType.SNV,
        gene_symbol="BRCA1",
        transcript_id="ENST00000357654",
        hgvs_c="c.5096G>A",
        hgvs_p="p.Arg1699Gln",
        consequence="missense_variant",
    )


def _pathogenic_report() -> EvidenceReport:
    evaluations = [
        CriterionEvaluation(
            code=ACMGCode.PVS1,
            strength=EvidenceStrength.VERY_STRONG,
            met=True,
            evidence_basis="stop_gained in an established LoF-mechanism gene",
            data_source="Ensembl VEP + gnomAD constraint",
        ),
        CriterionEvaluation(
            code=ACMGCode.PM2,
            strength=EvidenceStrength.SUPPORTING,
            met=True,
            evidence_basis="Absent from adequately covered gnomAD v4 populations",
            data_source="gnomAD v4",
        ),
    ]
    return EvidenceReport(
        variant=_brca1_variant(),
        metadata=ReportMetadata(alleleiq_version="0.1.0-alpha"),
        population_frequencies={
            AncestryGroup.NFE: PopulationFrequency(
                ancestry=AncestryGroup.NFE, allele_count=0, allele_number=68_034
            )
        },
        coverage_summary=CoverageSummary(
            overall_coverage_assessment=OverallCoverageAssessment.COMPLETE,
            clinical_advisory="All ancestry groups adequately covered.",
        ),
        criterion_evaluations=evaluations,
        bayesian_score=BayesianScore(
            evidence_points=9,
            criteria_applied=[ACMGCode.PVS1, ACMGCode.PM2],
            preliminary_classification=PreliminaryClassification.LIKELY_PATHOGENIC,
            classification_confidence=ClassificationConfidence.HIGH,
            evidence_points_confidence_interval=(9, 9),
        ),
        audit_trail=AuditTrail(),
    )


def _benign_report() -> EvidenceReport:
    evaluations = [
        CriterionEvaluation(
            code=ACMGCode.BA1,
            strength=EvidenceStrength.BENIGN_STANDALONE,
            met=True,
            evidence_basis=(
                "Allele frequency 0.1200 in Non-Finnish European exceeds the BA1 threshold"
            ),
            data_source="gnomAD v4",
        ),
    ]
    return EvidenceReport(
        variant=_brca1_variant(),
        metadata=ReportMetadata(alleleiq_version="0.1.0-alpha"),
        population_frequencies={
            AncestryGroup.NFE: PopulationFrequency(
                ancestry=AncestryGroup.NFE, allele_count=8000, allele_number=68_034
            )
        },
        coverage_summary=CoverageSummary(
            overall_coverage_assessment=OverallCoverageAssessment.COMPLETE,
            clinical_advisory="All ancestry groups adequately covered.",
        ),
        criterion_evaluations=evaluations,
        bayesian_score=BayesianScore(
            evidence_points=-8,
            criteria_applied=[ACMGCode.BA1],
            preliminary_classification=PreliminaryClassification.BENIGN,
            classification_confidence=ClassificationConfidence.HIGH,
            evidence_points_confidence_interval=(-8, -8),
        ),
        audit_trail=AuditTrail(),
    )


def _report_with_coverage_gaps() -> EvidenceReport:
    return EvidenceReport(
        variant=_brca1_variant(),
        metadata=ReportMetadata(alleleiq_version="0.1.0-alpha"),
        population_frequencies={
            AncestryGroup.NFE: PopulationFrequency(
                ancestry=AncestryGroup.NFE, allele_count=1, allele_number=68_034
            ),
        },
        coverage_gaps={
            AncestryGroup.AFR: CoverageGap(
                ancestry=AncestryGroup.AFR,
                allele_number=8_000,
                threshold=10_000,
                flag_level=CoverageFlagLevel.WARNING,
                clinical_guidance_message="Consult H3Africa.",
            ),
            AncestryGroup.SAS: CoverageGap(
                ancestry=AncestryGroup.SAS,
                allele_number=0,
                threshold=10_000,
                flag_level=CoverageFlagLevel.NO_DATA,
                clinical_guidance_message="No SAS data; consult GenomeAsia 100K.",
            ),
        },
        coverage_summary=CoverageSummary(
            populations_with_insufficient_data=[AncestryGroup.AFR],
            populations_with_no_data=[AncestryGroup.SAS],
            overall_coverage_assessment=OverallCoverageAssessment.PARTIAL,
            clinical_advisory="AFR and SAS coverage is limited.",
        ),
        criterion_evaluations=[],
        bayesian_score=BayesianScore(
            evidence_points=0,
            preliminary_classification=PreliminaryClassification.UNCERTAIN_SIGNIFICANCE,
            classification_confidence=ClassificationConfidence.MODERATE,
        ),
        audit_trail=AuditTrail(),
    )


class TestJSONReportSchema:
    def test_pathogenic_variant_json_report_matches_schema(self) -> None:
        payload = json.loads(render_json_report(_pathogenic_report()))

        assert payload["schema_version"] == "1.0.0"
        report = payload["report"]
        assert report["bayesian_score"]["preliminary_classification"] == "LIKELY_PATHOGENIC"
        assert report["bayesian_score"]["evidence_points"] == 9
        assert len(report["evidence_evaluation"]["criteria_applied"]) == 2
        assert report["clinical_disclaimer"] == CLINICAL_DISCLAIMER

    def test_benign_variant_json_report_matches_schema(self) -> None:
        payload = json.loads(render_json_report(_benign_report()))

        report = payload["report"]
        assert report["bayesian_score"]["preliminary_classification"] == "BENIGN"
        assert report["bayesian_score"]["evidence_points"] == -8
        assert report["evidence_evaluation"]["criteria_applied"][0]["code"] == "BA1"


class TestHTMLReportSections:
    def test_all_required_sections_present(self) -> None:
        html = render_html_report(_pathogenic_report())

        assert 'class="disclaimer-box"' in html
        assert 'class="variant-summary"' in html
        assert 'class="classification-banner' in html
        assert "Population frequency" in html
        assert "ACMG/AMP criteria evaluation" in html
        assert "Evidence not assessed" in html
        assert 'class="audit-trail"' in html
        assert "<footer>" in html

    def test_coverage_gap_rendering_shows_warning_and_critical_rows(self) -> None:
        html = render_html_report(_report_with_coverage_gaps())

        assert "gap-warning" in html
        assert "gap-critical" in html
        assert "H3Africa" in html
        assert "GenomeAsia" in html


class TestClinicalDisclaimerPresence:
    def test_disclaimer_present_in_json_output(self) -> None:
        assert CLINICAL_DISCLAIMER in render_json_report(_pathogenic_report())

    def test_disclaimer_present_in_html_output(self) -> None:
        assert CLINICAL_DISCLAIMER in render_html_report(_pathogenic_report())


class TestHTMLSnapshot:
    def test_html_report_snapshot(self, snapshot: object) -> None:
        report = EvidenceReport(
            variant=_brca1_variant(),
            metadata=ReportMetadata(
                report_id=UUID("00000000-0000-0000-0000-000000000000"),
                generated_at=datetime(2026, 1, 1, tzinfo=UTC),
                alleleiq_version="0.1.0-alpha",
                database_versions={"gnomad": "gnomad_r4", "vep": "112", "clinvar": "2026-01-01"},
            ),
            population_frequencies={
                AncestryGroup.NFE: PopulationFrequency(
                    ancestry=AncestryGroup.NFE, allele_count=0, allele_number=68_034
                )
            },
            coverage_gaps={
                AncestryGroup.SAS: CoverageGap(
                    ancestry=AncestryGroup.SAS,
                    allele_number=0,
                    threshold=10_000,
                    flag_level=CoverageFlagLevel.NO_DATA,
                    clinical_guidance_message="No SAS data; consult GenomeAsia 100K.",
                )
            },
            coverage_summary=CoverageSummary(
                populations_with_sufficient_data=[AncestryGroup.NFE],
                populations_with_no_data=[AncestryGroup.SAS],
                overall_coverage_assessment=OverallCoverageAssessment.PARTIAL,
                clinical_advisory="SAS coverage is absent.",
            ),
            criterion_evaluations=[
                CriterionEvaluation(
                    code=ACMGCode.PVS1,
                    strength=EvidenceStrength.VERY_STRONG,
                    met=True,
                    evidence_basis="stop_gained in an established LoF-mechanism gene",
                    data_source="Ensembl VEP + gnomAD constraint",
                ),
                CriterionEvaluation(
                    code=ACMGCode.PM2,
                    strength=EvidenceStrength.SUPPORTING,
                    met=True,
                    evidence_basis="Absent from adequately covered gnomAD v4 populations",
                    data_source="gnomAD v4",
                ),
            ],
            bayesian_score=BayesianScore(
                evidence_points=9,
                criteria_applied=[ACMGCode.PVS1, ACMGCode.PM2],
                preliminary_classification=PreliminaryClassification.LIKELY_PATHOGENIC,
                classification_confidence=ClassificationConfidence.HIGH,
                evidence_points_confidence_interval=(9, 9),
            ),
            audit_trail=AuditTrail(),
        )

        snapshot.snapshot_dir = str(Path(__file__).parent / "snapshots")  # type: ignore[attr-defined]
        snapshot.assert_match(render_html_report(report), "brca1_likely_pathogenic_report.html")  # type: ignore[attr-defined]
