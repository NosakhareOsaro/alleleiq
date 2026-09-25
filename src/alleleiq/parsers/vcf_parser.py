"""VCF parsing, normalisation, and validation.

Known limitation: this module performs *parsimonious*
indel normalisation — trimming shared trailing then leading bases
between REF and ALT and adjusting POS accordingly — which correctly
minimises the representation implied by the VCF record itself. It does
not perform full reference-guided left-alignment across a repeat region
that extends beyond the given REF/ALT window, which requires the
reference FASTA (not bundled in this release). VCF records produced by
standard variant callers are typically already left-aligned for exactly
this reason, so this covers the common case; callers working with
non-normalised upstream VCFs should run ``bcftools norm -f <reference>``
first.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from cyvcf2 import VCF
from pydantic import ValidationError

from alleleiq.exceptions import VariantParseError
from alleleiq.logging import get_logger
from alleleiq.models.variant import GenomeBuild, GenomicCoordinates, Variant, VariantType

logger = get_logger(__name__)

# Reference contig lengths (bp) for a handful of chromosomes, used to fingerprint
# genome build when a VCF lacks an explicit ##reference= header line. Source: NCBI
# GRCh38.p14 / GRCh37.p13 assembly reports.
_GRCH38_CONTIG_LENGTHS: dict[str, int] = {
    "1": 248_956_422,
    "13": 114_364_328,
    "17": 83_257_441,
    "21": 46_709_983,
    "X": 156_040_895,
}
_GRCH37_CONTIG_LENGTHS: dict[str, int] = {
    "1": 249_250_621,
    "13": 115_169_878,
    "17": 81_195_210,
    "21": 48_129_895,
    "X": 155_270_560,
}

_SYMBOLIC_ALT_PATTERN = re.compile(r"^<([A-Z:]+)>$")
_STRUCTURAL_TAGS = frozenset({"INV", "BND", "TRA", "INS"})
_CNV_TAGS = frozenset({"DEL", "DUP", "CNV"})


@dataclass(frozen=True)
class ParsingStatistics:
    """Summary counts from a single VCF parse, for logging and CLI display."""

    total_records: int
    variants_emitted: int
    variants_filtered: int
    variants_normalised: int
    genome_build_verified: bool


def parse_vcf(
    path: str | Path, expected_build: GenomeBuild = GenomeBuild.GRCH38
) -> tuple[list[Variant], ParsingStatistics]:
    """Parse a VCF file into a list of validated, normalised :class:`Variant` objects.

    Args:
        path: Path to the VCF file (plain or bgzipped).
        expected_build: The genome build this VCF is required to match.
            AlleleIQ supports GRCh38 only.

    Returns:
        A tuple of (parsed variants, parsing statistics). Multi-allelic
        sites are split into one biallelic ``Variant`` per alternate
        allele.

    Raises:
        VariantParseError: If the file does not exist, cannot be opened,
            contains no records, declares or fingerprints as a different
            genome build than ``expected_build``, or a record cannot be
            parsed by htslib.
    """
    vcf_path = Path(path)
    if not vcf_path.is_file():
        raise VariantParseError(f"VCF file not found: {vcf_path}")

    try:
        vcf_reader = VCF(str(vcf_path))
    except Exception as exc:
        raise VariantParseError(f"failed to open VCF file {vcf_path}: {exc}") from exc

    try:
        genome_build_verified = _verify_genome_build(vcf_reader, expected_build, vcf_path)
        variants, total_records, variants_filtered, variants_normalised = _parse_records(
            vcf_reader, expected_build
        )
    except VariantParseError:
        raise
    except Exception as exc:
        raise VariantParseError(f"error while parsing VCF file {vcf_path}: {exc}") from exc
    finally:
        vcf_reader.close()

    if total_records == 0:
        raise VariantParseError(f"VCF file contains no variant records: {vcf_path}")

    stats = ParsingStatistics(
        total_records=total_records,
        variants_emitted=len(variants),
        variants_filtered=variants_filtered,
        variants_normalised=variants_normalised,
        genome_build_verified=genome_build_verified,
    )
    logger.info(
        "VCF parsing complete",
        extra={
            "total_records": stats.total_records,
            "variants_emitted": stats.variants_emitted,
            "variants_filtered": stats.variants_filtered,
            "variants_normalised": stats.variants_normalised,
            "genome_build_verified": stats.genome_build_verified,
        },
    )
    return variants, stats


def _parse_records(
    vcf_reader: VCF, expected_build: GenomeBuild
) -> tuple[list[Variant], int, int, int]:
    variants: list[Variant] = []
    total_records = 0
    variants_filtered = 0
    variants_normalised = 0

    for record in vcf_reader:
        total_records += 1
        if not record.REF or not record.ALT:
            variants_filtered += 1
            logger.warning(
                "skipping record with missing REF/ALT",
                extra={"chrom": record.CHROM, "pos": record.POS},
            )
            continue

        for alt in record.ALT:
            try:
                ref_norm, alt_norm, pos_norm, was_normalised = _normalise_allele(
                    record.REF, alt, record.POS
                )
                coordinates = GenomicCoordinates(
                    chrom=record.CHROM,
                    pos=pos_norm,
                    ref=ref_norm,
                    alt=alt_norm,
                    genome_build=expected_build,
                )
                variant = Variant(
                    coordinates=coordinates,
                    variant_type=_classify_variant_type(coordinates.ref, coordinates.alt),
                )
            except (ValidationError, ValueError) as exc:
                variants_filtered += 1
                logger.warning(
                    "skipping invalid variant record",
                    extra={
                        "chrom": record.CHROM,
                        "pos": record.POS,
                        "alt": alt,
                        "reason": str(exc),
                    },
                )
                continue

            variants.append(variant)
            if was_normalised:
                variants_normalised += 1

    return variants, total_records, variants_filtered, variants_normalised


def _normalise_allele(ref: str, alt: str, pos: int) -> tuple[str, str, int, bool]:
    """Trim shared trailing then leading bases between REF and ALT.

    Symbolic alleles (``<DEL>`` etc.) and single-base substitutions are
    returned unchanged. At least one base is always retained in both
    REF and ALT (the VCF anchor-base convention), so trimming never
    empties either field.

    Returns:
        ``(ref, alt, pos, was_normalised)``.
    """
    if _SYMBOLIC_ALT_PATTERN.match(alt) or (len(ref) == 1 and len(alt) == 1):
        return ref, alt, pos, False

    original = (ref, alt, pos)
    ref_chars, alt_chars = list(ref), list(alt)

    while len(ref_chars) > 1 and len(alt_chars) > 1 and ref_chars[-1] == alt_chars[-1]:
        ref_chars.pop()
        alt_chars.pop()

    while len(ref_chars) > 1 and len(alt_chars) > 1 and ref_chars[0] == alt_chars[0]:
        ref_chars.pop(0)
        alt_chars.pop(0)
        pos += 1

    normalised = ("".join(ref_chars), "".join(alt_chars), pos)
    return *normalised, normalised != original


def _classify_variant_type(ref: str, alt: str) -> VariantType:
    symbolic_match = _SYMBOLIC_ALT_PATTERN.match(alt)
    if symbolic_match:
        tag = symbolic_match.group(1).split(":")[0]
        if tag in _CNV_TAGS:
            return VariantType.CNV
        return VariantType.STRUCTURAL
    if len(ref) == 1 and len(alt) == 1:
        return VariantType.SNV
    # Equal-length multi-base substitutions (MNVs) have no dedicated category in
    # AlleleIQ's current VariantType taxonomy and are classified as INDEL; this is
    # a known simplification and rare in practice (most callers emit MNVs as
    # decomposed adjacent SNVs).
    return VariantType.INDEL


def _verify_genome_build(vcf_reader: VCF, expected_build: GenomeBuild, path: Path) -> bool:
    """Confirm the VCF's declared or fingerprinted genome build matches expectations.

    Returns:
        True if the build was positively confirmed via header or contig
        length fingerprinting; False if it could not be determined
        (parsing proceeds, but this is surfaced in
        :attr:`ParsingStatistics.genome_build_verified` for visibility).

    Raises:
        VariantParseError: If the VCF declares or fingerprints as a
            genome build other than ``expected_build``.
    """
    declared = _declared_build_from_header(vcf_reader.raw_header)
    if declared is not None:
        if declared != expected_build:
            raise VariantParseError(
                f"VCF file {path} declares genome build {declared.value}, but AlleleIQ "
                f"requires {expected_build.value}. Coordinate systems are not "
                "interchangeable between builds; re-run with a build-converted VCF."
            )
        return True

    try:
        seqlens = vcf_reader.seqlens
    except AttributeError:
        seqlens = None
    fingerprinted = _fingerprint_build_from_contigs(vcf_reader.seqnames, seqlens)
    if fingerprinted is not None:
        if fingerprinted != expected_build:
            raise VariantParseError(
                f"VCF file {path} contig lengths match genome build {fingerprinted.value}, "
                f"but AlleleIQ requires {expected_build.value}."
            )
        return True

    logger.warning(
        "could not verify VCF genome build from header or contig lengths; "
        f"proceeding assuming {expected_build.value}",
        extra={"path": str(path)},
    )
    return False


def _declared_build_from_header(raw_header: str) -> GenomeBuild | None:
    match = re.search(r"^##reference=(.+)$", raw_header, re.MULTILINE)
    if match is None:
        return None
    reference_value = match.group(1).strip().lower()
    if "grch38" in reference_value or "hg38" in reference_value:
        return GenomeBuild.GRCH38
    if "grch37" in reference_value or "hg19" in reference_value:
        return GenomeBuild.GRCH37
    return None


def _fingerprint_build_from_contigs(
    seqnames: list[str], seqlens: list[int] | None
) -> GenomeBuild | None:
    if not seqlens:
        return None
    lengths_by_name = dict(zip(seqnames, seqlens, strict=False))
    for chrom, length in lengths_by_name.items():
        normalised_chrom = chrom.removeprefix("chr")
        if _GRCH38_CONTIG_LENGTHS.get(normalised_chrom) == length:
            return GenomeBuild.GRCH38
        if _GRCH37_CONTIG_LENGTHS.get(normalised_chrom) == length:
            return GenomeBuild.GRCH37
    return None
