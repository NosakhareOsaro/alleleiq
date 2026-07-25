"""Variant identity and genomic coordinate models."""

from __future__ import annotations

import re
from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

_VALID_CHROMOSOMES = frozenset(
    [str(n) for n in range(1, 23)] + ["X", "Y", "MT"],
)
_NUCLEOTIDE_PATTERN = re.compile(r"^[ACGTN]+$")
_SYMBOLIC_ALLELE_PATTERN = re.compile(r"^<[A-Z:]+>$")


class GenomeBuild(StrEnum):
    """Supported reference genome builds.

    AlleleIQ operates exclusively on GRCh38. GRCh37 is modelled here so
    the parser can *detect* it and raise a clear
    :class:`~alleleiq.exceptions.VariantParseError` rather than silently
    misinterpreting coordinates against the wrong reference.
    """

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"


class VariantType(StrEnum):
    """Broad structural category of a genomic variant."""

    SNV = "SNV"
    INDEL = "INDEL"
    CNV = "CNV"
    STRUCTURAL = "STRUCTURAL"


class GenomicCoordinates(BaseModel):
    """A normalised, zero-ambiguity genomic position and allele change.

    Coordinates are expected to already be left-aligned and normalised
    (see :mod:`alleleiq.parsers.vcf_parser`) by the time they reach this
    model; this model validates structural well-formedness, not
    normalisation correctness.
    """

    model_config = ConfigDict(frozen=True)

    chrom: str = Field(description="Chromosome without a 'chr' prefix, e.g. '17', 'X', 'MT'.")
    pos: int = Field(gt=0, description="1-based genomic position (VCF convention).")
    ref: str = Field(min_length=1, description="Reference allele.")
    alt: str = Field(min_length=1, description="Alternate allele.")
    genome_build: GenomeBuild = Field(default=GenomeBuild.GRCH38)

    @field_validator("chrom", mode="before")
    @classmethod
    def _normalise_chrom(cls, value: str) -> str:
        normalised = str(value).strip()
        if normalised.lower().startswith("chr"):
            normalised = normalised[3:]
        normalised = (
            normalised.upper() if normalised.upper() in {"X", "Y", "MT", "M"} else normalised
        )
        if normalised == "M":
            normalised = "MT"
        if normalised not in _VALID_CHROMOSOMES:
            raise ValueError(f"unrecognised chromosome {value!r}")
        return normalised

    @field_validator("ref", "alt", mode="before")
    @classmethod
    def _uppercase_allele(cls, value: str) -> str:
        return str(value).strip().upper()

    @field_validator("ref", "alt")
    @classmethod
    def _validate_allele_characters(cls, value: str) -> str:
        if _NUCLEOTIDE_PATTERN.match(value) or _SYMBOLIC_ALLELE_PATTERN.match(value):
            return value
        raise ValueError(
            f"allele {value!r} is neither a valid nucleotide sequence (A/C/G/T/N) "
            "nor a symbolic structural allele (e.g. '<DEL>')"
        )

    @property
    def variant_id(self) -> str:
        """A stable, human-readable identifier: ``chrom-pos-ref-alt``."""
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


class Variant(BaseModel):
    """A single genomic variant with its functional annotation context.

    Functional annotation fields (``gene_symbol`` through ``consequence``)
    are populated by :mod:`alleleiq.clients.vep` and are optional here
    because a bare parsed variant, before annotation, is still a valid
    ``Variant``.
    """

    model_config = ConfigDict(frozen=True)

    coordinates: GenomicCoordinates
    variant_type: VariantType
    gene_symbol: str | None = Field(default=None)
    transcript_id: str | None = Field(
        default=None, description="Ensembl transcript ID, e.g. ENST00000357654."
    )
    hgvs_c: str | None = Field(default=None, description="HGVS coding-sequence notation.")
    hgvs_p: str | None = Field(default=None, description="HGVS protein notation.")
    consequence: str | None = Field(
        default=None,
        description="Most severe Sequence Ontology consequence term, e.g. 'missense_variant'.",
    )

    @model_validator(mode="after")
    def _validate_type_allele_consistency(self) -> Variant:
        is_symbolic = _SYMBOLIC_ALLELE_PATTERN.match(self.coordinates.alt) is not None
        if is_symbolic and self.variant_type not in (VariantType.CNV, VariantType.STRUCTURAL):
            raise ValueError(
                f"symbolic allele {self.coordinates.alt!r} requires variant_type "
                f"CNV or STRUCTURAL, got {self.variant_type}"
            )
        if not is_symbolic and self.variant_type in (VariantType.CNV, VariantType.STRUCTURAL):
            raise ValueError(
                f"variant_type {self.variant_type} requires a symbolic allele (e.g. '<DEL>'), "
                f"got {self.coordinates.alt!r}"
            )
        return self

    @property
    def variant_id(self) -> str:
        """Delegate to :attr:`GenomicCoordinates.variant_id`."""
        return self.coordinates.variant_id
