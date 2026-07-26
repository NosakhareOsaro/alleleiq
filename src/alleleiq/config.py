"""Application configuration.

All configurable parameters are declared here as a single
:class:`pydantic_settings.BaseSettings` model, sourced from environment
variables (prefixed ``ALLELEIQ_``) or a ``.env`` file. See
``.env.example`` for the full list of settings with documentation of
each value's clinical or operational rationale.

Centralising configuration here — rather than scattering default
thresholds through the scoring modules — satisfies the "no magic
numbers" requirement: every clinically meaningful threshold is a named,
documented, overridable setting.
"""

from __future__ import annotations

from functools import lru_cache

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Runtime configuration for AlleleIQ.

    Instantiate via :func:`get_settings` rather than directly, so that
    the process uses a single cached instance.
    """

    model_config = SettingsConfigDict(
        env_prefix="ALLELEIQ_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # ── General ──────────────────────────────────────────────────────
    log_level: str = Field(default="INFO", description="DEBUG, INFO, WARNING, or ERROR.")
    genome_build: str = Field(
        default="GRCh38",
        description="Genome build accepted by the VCF parser. GRCh37 input is rejected.",
    )

    # ── gnomAD v4 GraphQL client ─────────────────────────────────────
    gnomad_api_url: str = Field(default="https://gnomad.broadinstitute.org/api")
    gnomad_dataset: str = Field(default="gnomad_r4")
    gnomad_timeout_seconds: float = Field(default=30.0, gt=0)
    gnomad_max_retries: int = Field(default=3, ge=0)
    gnomad_min_request_interval_seconds: float = Field(default=0.5, ge=0)
    gnomad_cache_ttl_seconds: int = Field(default=3600, ge=0)

    # ── Ensembl VEP REST client ──────────────────────────────────────
    vep_api_url: str = Field(default="https://rest.ensembl.org")
    vep_timeout_seconds: float = Field(default=30.0, gt=0)
    vep_max_retries: int = Field(default=3, ge=0)

    # ── Local REVEL tabix lookup ──────────────────────────────────────
    revel_data_path: str | None = Field(
        default=None,
        description=(
            "Path to a bgzipped, tabix-indexed REVEL scores TSV (see "
            "scripts/download_revel.sh). Optional: PP3/BP4 fall back to "
            "VEP REST's REVEL score, then to requires_clinical_input, "
            "when unset or the file is absent."
        ),
    )

    # ── Local ClinGen gene-disease validity + dosage sensitivity lookup ─
    clingen_gene_validity_path: str = Field(
        default="data/clingen/gene_disease_validity.csv",
        description=(
            "Path to the bundled ClinGen Gene-Disease Validity CSV snapshot "
            "(see scripts/update_clingen.sh). Missing file degrades PVS1's "
            "haploinsufficiency-mechanism strength to the bundled fallback "
            "gene list / gnomAD constraint path."
        ),
    )
    clingen_dosage_path: str = Field(
        default="data/clingen/dosage_sensitivity.tsv",
        description=(
            "Path to the bundled ClinGen Dosage Sensitivity Map TSV snapshot "
            "(see scripts/update_clingen.sh)."
        ),
    )

    # ── ClinVar (NCBI E-utilities) client ────────────────────────────
    clinvar_api_url: str = Field(default="https://eutils.ncbi.nlm.nih.gov/entrez/eutils")
    clinvar_api_key: str | None = Field(default=None)
    clinvar_timeout_seconds: float = Field(default=30.0, gt=0)
    clinvar_max_retries: int = Field(default=3, ge=0)
    clinvar_requests_per_second_no_key: float = Field(default=3.0, gt=0)
    clinvar_requests_per_second_with_key: float = Field(default=10.0, gt=0)
    clinvar_min_stars_strong_evidence: int = Field(default=2, ge=0, le=4)
    same_position_window_bp: int = Field(
        default=3,
        ge=0,
        description=(
            "Base pairs either side of a variant's genomic position searched in ClinVar "
            "for PS1/PM5's same-position comparison (see "
            "alleleiq.scoring.criteria.same_position). AlleleIQ does not perform "
            "strand-aware transcript-to-codon mapping, so this is a symmetric genomic "
            "window rather than a true codon boundary; 3bp either side of a single-base "
            "substitution covers its own codon in the common case while staying narrow "
            "enough to avoid pulling in neighbouring, unrelated codons."
        ),
    )

    # ── Population coverage thresholds (minimum gnomAD v4 allele number) ─
    coverage_threshold_afr: int = Field(default=10_000, ge=0)
    coverage_threshold_amr: int = Field(default=5_000, ge=0)
    coverage_threshold_asj: int = Field(default=3_000, ge=0)
    coverage_threshold_eas: int = Field(default=5_000, ge=0)
    coverage_threshold_fin: int = Field(default=3_000, ge=0)
    coverage_threshold_nfe: int = Field(default=20_000, ge=0)
    coverage_threshold_sas: int = Field(default=10_000, ge=0)
    coverage_threshold_oth: int = Field(default=1_000, ge=0)
    coverage_critical_fraction: float = Field(default=0.5, gt=0, lt=1)

    # ── ACMG/AMP criteria thresholds (ACGS 2024 aligned) ─────────────
    ba1_frequency_threshold: float = Field(default=0.05, gt=0, lt=1)
    bs1_threshold_dominant: float = Field(default=0.01, gt=0, lt=1)
    bs1_threshold_recessive: float = Field(default=0.03, gt=0, lt=1)
    revel_pp3_supporting_threshold: float = Field(default=0.70, ge=0, le=1)
    revel_pp3_moderate_threshold: float = Field(default=0.90, ge=0, le=1)
    revel_bp4_threshold: float = Field(default=0.15, ge=0, le=1)
    spliceai_pp3_supporting_threshold: float = Field(default=0.20, ge=0, le=1)
    spliceai_pp3_moderate_threshold: float = Field(default=0.50, ge=0, le=1)
    spliceai_bp4_threshold: float = Field(default=0.10, ge=0, le=1)
    pli_haploinsufficiency_threshold: float = Field(default=0.90, ge=0, le=1)
    loeuf_haploinsufficiency_threshold: float = Field(default=0.35, ge=0)

    # ── API server ────────────────────────────────────────────────────
    api_cors_allow_origins: list[str] = Field(default_factory=lambda: ["http://localhost:3000"])
    api_rate_limit_per_hour: int = Field(default=100, gt=0)
    api_max_upload_bytes: int = Field(default=52_428_800, gt=0)

    @field_validator("api_cors_allow_origins", mode="before")
    @classmethod
    def _split_csv(cls, value: object) -> object:
        """Allow ``ALLELEIQ_API_CORS_ALLOW_ORIGINS`` to be a comma-separated string."""
        if isinstance(value, str):
            return [origin.strip() for origin in value.split(",") if origin.strip()]
        return value

    @field_validator("log_level")
    @classmethod
    def _validate_log_level(cls, value: str) -> str:
        allowed = {"DEBUG", "INFO", "WARNING", "ERROR"}
        normalised = value.upper()
        if normalised not in allowed:
            raise ValueError(f"log_level must be one of {sorted(allowed)}, got {value!r}")
        return normalised

    @property
    def coverage_thresholds_by_population_id(self) -> dict[str, int]:
        """Map gnomAD v4 population IDs (lowercase) to minimum allele number.

        Returns:
            Dict such as ``{"afr": 10000, "amr": 5000, ...}``, keyed by
            the lowercase population identifiers used in gnomAD's API,
            not by AlleleIQ's :class:`~alleleiq.models.population.AncestryGroup`
            enum (kept independent here to avoid a config -> models
            import cycle).
        """
        return {
            "afr": self.coverage_threshold_afr,
            "amr": self.coverage_threshold_amr,
            "asj": self.coverage_threshold_asj,
            "eas": self.coverage_threshold_eas,
            "fin": self.coverage_threshold_fin,
            "nfe": self.coverage_threshold_nfe,
            "sas": self.coverage_threshold_sas,
            "oth": self.coverage_threshold_oth,
        }

    @property
    def clinvar_min_request_interval_seconds(self) -> float:
        """Minimum seconds between ClinVar requests, given NCBI's key-gated rate limits."""
        requests_per_second = (
            self.clinvar_requests_per_second_with_key
            if self.clinvar_api_key
            else self.clinvar_requests_per_second_no_key
        )
        return 1.0 / requests_per_second


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return the process-wide cached :class:`Settings` instance.

    Returns:
        A single shared ``Settings`` instance, constructed once per
        process from the environment and ``.env`` file.
    """
    return Settings()
