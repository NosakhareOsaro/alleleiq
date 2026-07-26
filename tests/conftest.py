"""Shared pytest fixtures for the AlleleIQ test suite."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from alleleiq.__version__ import VERSION
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.models import (
    AncestryGroup,
    AuditTrail,
    BayesianScore,
    ClassificationConfidence,
    CoverageSummary,
    GenomicCoordinates,
    OverallCoverageAssessment,
    PopulationFrequency,
    PreliminaryClassification,
    ReportMetadata,
    Variant,
    VariantType,
)

FIXTURES_DIR = Path(__file__).parent / "fixtures"
VCF_FIXTURES_DIR = FIXTURES_DIR / "vcf"
RESPONSE_FIXTURES_DIR = FIXTURES_DIR / "responses"


@pytest.fixture(autouse=True)
def _clear_gnomad_client_cache() -> None:
    """GnomADClient's response cache is a process-wide class-level singleton by
    design (see commit 40 — it needs to survive across the fresh client instances
    alleleiq.pipeline constructs per variant, including the fake subclasses the API
    integration tests monkeypatch in). Applied suite-wide, not just in
    test_gnomad_client.py, since a test in one file could otherwise silently reuse
    another test's cached mock response for the same variant coordinates.
    """
    GnomADClient.clear_cache()


@pytest.fixture
def vcf_fixtures_dir() -> Path:
    """Directory containing sample VCF files used across parser tests."""
    return VCF_FIXTURES_DIR


@pytest.fixture
def single_snv_vcf_path() -> Path:
    """A minimal, well-formed single-SNV VCF file (BRCA1, GRCh38)."""
    return VCF_FIXTURES_DIR / "single_snv.vcf"


@pytest.fixture
def multi_variant_vcf_path() -> Path:
    """A VCF containing a multi-allelic site and an indel, GRCh38."""
    return VCF_FIXTURES_DIR / "multi_variant.vcf"


@pytest.fixture
def malformed_vcf_path() -> Path:
    """A VCF with a non-numeric POS field, for error-path testing."""
    return VCF_FIXTURES_DIR / "malformed.vcf"


def _load_response_fixture(name: str) -> dict[str, Any]:
    with (RESPONSE_FIXTURES_DIR / name).open() as handle:
        result: dict[str, Any] = json.load(handle)
        return result


@pytest.fixture
def gnomad_brca2_response() -> dict[str, Any]:
    """A realistic gnomAD v4 GraphQL response with mixed population coverage."""
    return _load_response_fixture("gnomad_brca2.json")


@pytest.fixture
def gnomad_no_data_response() -> dict[str, Any]:
    """A gnomAD v4 GraphQL response for a variant absent from gnomAD entirely."""
    return _load_response_fixture("gnomad_no_data.json")


def _load_vep_fixture(name: str) -> list[dict[str, Any]]:
    with (RESPONSE_FIXTURES_DIR / name).open() as handle:
        result: list[dict[str, Any]] = json.load(handle)
        return result


@pytest.fixture
def vep_missense_response() -> list[dict[str, Any]]:
    """A realistic Ensembl VEP REST response for a BRCA1 missense variant."""
    return _load_vep_fixture("vep_missense.json")


@pytest.fixture
def vep_stop_gained_response() -> list[dict[str, Any]]:
    """A VEP response for a BRCA1 nonsense (stop_gained) variant."""
    return _load_vep_fixture("vep_stop_gained.json")


@pytest.fixture
def vep_splice_site_response() -> list[dict[str, Any]]:
    """A VEP response for a canonical splice donor variant, with SpliceAI scores."""
    return _load_vep_fixture("vep_splice_site.json")


@pytest.fixture
def vep_non_coding_response() -> list[dict[str, Any]]:
    """A VEP response for an intergenic (non-coding) variant with no transcript consequences."""
    return _load_vep_fixture("vep_non_coding.json")


@pytest.fixture
def clinvar_expert_panel_xml() -> str:
    """A ClinVar VCV XML record: 3-star expert panel Pathogenic classification."""
    return (RESPONSE_FIXTURES_DIR / "clinvar_expert_panel.xml").read_text()


@pytest.fixture
def clinvar_conflicting_xml() -> str:
    """A ClinVar VCV XML record: 1-star conflicting-classifications submission."""
    return (RESPONSE_FIXTURES_DIR / "clinvar_conflicting.xml").read_text()


@pytest.fixture
def brca1_coordinates() -> GenomicCoordinates:
    """Canonical BRCA1 missense variant coordinates (17-43094692-G-A, GRCh38)."""
    return GenomicCoordinates(chrom="17", pos=43094692, ref="G", alt="A")


@pytest.fixture
def brca1_variant(brca1_coordinates: GenomicCoordinates) -> Variant:
    """A fully annotated BRCA1 missense Variant."""
    return Variant(
        coordinates=brca1_coordinates,
        variant_type=VariantType.SNV,
        gene_symbol="BRCA1",
        transcript_id="ENST00000357654",
        hgvs_c="c.5096G>A",
        hgvs_p="p.Arg1699Gln",
        consequence="missense_variant",
    )


@pytest.fixture
def sufficient_population_frequency() -> PopulationFrequency:
    """A well-covered population with a rare observed allele."""
    return PopulationFrequency(ancestry=AncestryGroup.NFE, allele_count=3, allele_number=68034)


@pytest.fixture
def no_data_population_frequency() -> PopulationFrequency:
    """A population with zero allele number: absent data, not zero frequency."""
    return PopulationFrequency(ancestry=AncestryGroup.SAS, allele_count=0, allele_number=0)


@pytest.fixture
def report_metadata() -> ReportMetadata:
    """Minimal valid ReportMetadata for the currently installed AlleleIQ version."""
    return ReportMetadata(alleleiq_version=VERSION)


@pytest.fixture
def empty_audit_trail() -> AuditTrail:
    """An AuditTrail with no recorded queries, for tests that don't exercise it."""
    return AuditTrail()


@pytest.fixture
def complete_coverage_summary() -> CoverageSummary:
    """A CoverageSummary representing full population coverage."""
    return CoverageSummary(
        overall_coverage_assessment=OverallCoverageAssessment.COMPLETE,
        clinical_advisory="All ancestry groups met the minimum allele number threshold.",
    )


@pytest.fixture
def vus_bayesian_score() -> BayesianScore:
    """A BayesianScore representing an uninformative, low-evidence VUS."""
    return BayesianScore(
        evidence_points=0,
        preliminary_classification=PreliminaryClassification.UNCERTAIN_SIGNIFICANCE,
        classification_confidence=ClassificationConfidence.LOW,
    )
