"""Population frequency and ancestry-coverage models.

These models are the clinical-equity core of AlleleIQ: every
frequency-derived value is scoped to a single ancestry group, and every
ancestry group with insufficient reference data is represented by an
explicit :class:`CoverageGap` rather than being silently folded into an
aggregate.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, model_validator


class AncestryGroup(StrEnum):
    """gnomAD v4 genetic ancestry groups.

    Values match gnomAD's lowercase population identifiers so they can be
    used directly as GraphQL response keys.
    """

    AFR = "afr"
    AMR = "amr"
    ASJ = "asj"
    EAS = "eas"
    FIN = "fin"
    NFE = "nfe"
    SAS = "sas"
    OTH = "oth"

    @property
    def display_name(self) -> str:
        """Human-readable population label for reports."""
        return {
            AncestryGroup.AFR: "African/African American",
            AncestryGroup.AMR: "Latino/Admixed American",
            AncestryGroup.ASJ: "Ashkenazi Jewish",
            AncestryGroup.EAS: "East Asian",
            AncestryGroup.FIN: "Finnish",
            AncestryGroup.NFE: "Non-Finnish European",
            AncestryGroup.SAS: "South Asian",
            AncestryGroup.OTH: "Other/Unassigned",
        }[self]


class CoverageFlagLevel(StrEnum):
    """Severity of a population's reference-data coverage gap."""

    WARNING = "WARNING"
    CRITICAL = "CRITICAL"
    NO_DATA = "NO_DATA"


class OverallCoverageAssessment(StrEnum):
    """Report-level summary of coverage across all ancestry groups."""

    COMPLETE = "COMPLETE"
    PARTIAL = "PARTIAL"
    SEVERELY_LIMITED = "SEVERELY_LIMITED"


class PopulationFrequency(BaseModel):
    """Allele frequency data for a single variant in a single ancestry group.

    An ``allele_number`` of zero means the population was not sampled at
    this position (no data), which is not the same as observing the
    variant at zero frequency in a well-sampled population. This model
    encodes that distinction explicitly: ``allele_frequency`` is forced
    to ``None`` whenever ``allele_number`` is zero.
    """

    model_config = ConfigDict(frozen=True)

    ancestry: AncestryGroup
    allele_count: int = Field(ge=0, description="Number of observed alternate alleles (AC).")
    allele_number: int = Field(ge=0, description="Total alleles genotyped in this population (AN).")
    allele_frequency: float | None = Field(
        default=None, ge=0, le=1, description="AC / AN, or None if allele_number is 0."
    )
    homozygote_count: int = Field(default=0, ge=0)
    confidence_interval_95: tuple[float, float] | None = Field(
        default=None, description="95% confidence interval on the allele frequency."
    )

    @model_validator(mode="after")
    def _reconcile_frequency_with_allele_number(self) -> PopulationFrequency:
        if self.allele_count > self.allele_number:
            raise ValueError(
                f"allele_count ({self.allele_count}) cannot exceed "
                f"allele_number ({self.allele_number})"
            )
        if self.allele_number == 0:
            if self.allele_frequency is not None:
                object.__setattr__(self, "allele_frequency", None)
            return self
        expected_af = self.allele_count / self.allele_number
        if self.allele_frequency is None:
            object.__setattr__(self, "allele_frequency", expected_af)
        return self

    @property
    def has_data(self) -> bool:
        """Whether this population was sampled at all (AN > 0)."""
        return self.allele_number > 0


class CoverageGap(BaseModel):
    """A flagged insufficiency of reference data for one ancestry group."""

    model_config = ConfigDict(frozen=True)

    ancestry: AncestryGroup
    allele_number: int = Field(ge=0)
    threshold: int = Field(
        gt=0, description="Configured minimum allele number for this population."
    )
    flag_level: CoverageFlagLevel
    clinical_guidance_message: str = Field(
        description="Specific, actionable guidance for the reviewing clinical scientist."
    )
