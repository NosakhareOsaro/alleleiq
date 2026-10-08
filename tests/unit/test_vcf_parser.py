"""Unit tests for alleleiq.parsers.vcf_parser.parse_vcf."""

from __future__ import annotations

from pathlib import Path

import pytest

from alleleiq.exceptions import VariantParseError
from alleleiq.models.variant import VariantType
from alleleiq.parsers.vcf_parser import parse_vcf

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##reference={reference}\n"
    "##contig=<ID=17,length={contig_length}>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _write_vcf(
    tmp_path: Path, records: str, reference: str = "GRCh38", contig_length: int = 83_257_441
) -> Path:
    vcf_path = tmp_path / "variants.vcf"
    vcf_path.write_text(
        _VCF_HEADER.format(reference=reference, contig_length=contig_length) + records
    )
    return vcf_path


class TestSingleSNV:
    def test_parses_single_snv(self, single_snv_vcf_path: Path) -> None:
        variants, stats = parse_vcf(single_snv_vcf_path)

        assert len(variants) == 1
        assert stats.total_records == 1
        assert stats.variants_emitted == 1
        assert stats.variants_filtered == 0
        assert stats.genome_build_verified is True

        variant = variants[0]
        assert variant.coordinates.chrom == "17"
        assert variant.coordinates.pos == 43094692
        assert variant.coordinates.ref == "G"
        assert variant.coordinates.alt == "A"
        assert variant.variant_type == VariantType.SNV


class TestMultiAllelicSite:
    def test_splits_multiallelic_site_into_biallelic_variants(
        self, multi_variant_vcf_path: Path
    ) -> None:
        variants, stats = parse_vcf(multi_variant_vcf_path)

        assert stats.total_records == 3
        assert stats.variants_emitted == 4

        ids = {v.variant_id for v in variants}
        assert "17-43094692-G-A" in ids
        assert "17-43094692-G-T" in ids

    def test_indel_present_and_classified(self, multi_variant_vcf_path: Path) -> None:
        variants, _ = parse_vcf(multi_variant_vcf_path)

        indel = next(v for v in variants if v.variant_type == VariantType.INDEL)
        assert indel.coordinates.ref == "AG"
        assert indel.coordinates.alt == "A"


class TestIndelNormalisation:
    def test_trims_shared_suffix_and_prefix_and_shifts_position(self, tmp_path: Path) -> None:
        # REF=CAT, ALT=CATT: shared suffix 'T' trims once (CA/CAT), then shared
        # prefix 'C' trims once more (pos +1), leaving REF=A, ALT=AT at pos+1.
        vcf_path = _write_vcf(tmp_path, "17\t100\t.\tCAT\tCATT\t100\tPASS\t.\n")

        variants, stats = parse_vcf(vcf_path)

        assert stats.variants_normalised == 1
        variant = variants[0]
        assert variant.coordinates.pos == 101
        assert variant.coordinates.ref == "A"
        assert variant.coordinates.alt == "AT"

    def test_already_minimal_indel_is_not_flagged_as_normalised(
        self, multi_variant_vcf_path: Path
    ) -> None:
        _, stats = parse_vcf(multi_variant_vcf_path)
        assert stats.variants_normalised == 0


class TestGenomeBuildValidation:
    def test_grch37_header_raises_variant_parse_error(self, tmp_path: Path) -> None:
        vcf_path = _write_vcf(
            tmp_path,
            "17\t43094692\t.\tG\tA\t100\tPASS\t.\n",
            reference="GRCh37",
            contig_length=81_195_210,
        )

        with pytest.raises(VariantParseError, match="GRCh37"):
            parse_vcf(vcf_path)

    def test_grch37_contig_length_fingerprint_raises_without_header(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "no_reference_header.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=81195210>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t43094692\t.\tG\tA\t100\tPASS\t.\n"
        )

        with pytest.raises(VariantParseError, match="GRCh37"):
            parse_vcf(vcf_path)

    def test_unverifiable_build_proceeds_with_warning_flag(self, tmp_path: Path) -> None:
        vcf_path = tmp_path / "no_build_info.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t43094692\t.\tG\tA\t100\tPASS\t.\n"
        )

        variants, stats = parse_vcf(vcf_path)

        assert len(variants) == 1
        assert stats.genome_build_verified is False


class TestMalformedInput:
    def test_malformed_record_raises_variant_parse_error(self, malformed_vcf_path: Path) -> None:
        with pytest.raises(VariantParseError):
            parse_vcf(malformed_vcf_path)

    def test_invalid_allele_characters_are_filtered_not_fatal(self, tmp_path: Path) -> None:
        vcf_path = _write_vcf(
            tmp_path,
            "17\t43094692\t.\tG\tA\t100\tPASS\t.\n17\t43094700\t.\tG\tZ\t100\tPASS\t.\n",
        )

        variants, stats = parse_vcf(vcf_path)

        assert stats.total_records == 2
        assert stats.variants_filtered == 1
        assert len(variants) == 1

    def test_empty_file_raises_variant_parse_error(self, tmp_path: Path) -> None:
        vcf_path = _write_vcf(tmp_path, "")

        with pytest.raises(VariantParseError, match="no variant records"):
            parse_vcf(vcf_path)

    def test_missing_file_raises_variant_parse_error(self, tmp_path: Path) -> None:
        with pytest.raises(VariantParseError, match="not found"):
            parse_vcf(tmp_path / "does_not_exist.vcf")


class TestDuplicateVariants:
    def test_repeated_variant_is_kept_once(self, tmp_path: Path) -> None:
        vcf_path = _write_vcf(
            tmp_path,
            "17\t100\t.\tA\tG\t100\tPASS\t.\n"
            "17\t200\t.\tC\tT\t100\tPASS\t.\n"
            "17\t100\t.\tA\tG\t100\tPASS\t.\n",
        )

        variants, stats = parse_vcf(vcf_path)

        assert [v.variant_id for v in variants] == ["17-100-A-G", "17-200-C-T"]
        assert stats.total_records == 3
        assert stats.variants_emitted == 2
        assert stats.duplicates_removed == 1

    def test_records_that_normalise_to_the_same_variant_are_deduplicated(
        self, tmp_path: Path
    ) -> None:
        vcf_path = _write_vcf(
            tmp_path,
            # CAT>CATT at 100 normalises to A>AT at 101, the same variant as the second record.
            "17\t100\t.\tCAT\tCATT\t100\tPASS\t.\n" "17\t101\t.\tA\tAT\t100\tPASS\t.\n",
        )

        variants, stats = parse_vcf(vcf_path)

        assert [v.variant_id for v in variants] == ["17-101-A-AT"]
        assert stats.duplicates_removed == 1

    def test_no_duplicates_reports_zero(self, single_snv_vcf_path: Path) -> None:
        _, stats = parse_vcf(single_snv_vcf_path)

        assert stats.duplicates_removed == 0
