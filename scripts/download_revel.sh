#!/usr/bin/env bash
#
# Downloads the REVEL v1.3 precomputed scores (all chromosomes, GRCh38
# coordinates), and builds a bgzipped, tabix-indexed lookup file for
# alleleiq.clients.revel.REVELClient.
#
# The upstream archive (https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip)
# is ~670MB zipped and its extracted CSV is several GB; the resulting
# bgzipped TSV is still well over GitHub's 100MB per-file limit, so this
# repository does NOT bundle the real data file — only this script and a
# small synthetic fixture for tests (tests/fixtures/revel/). Run this
# script once, locally or on your deployment host, and point
# ALLELEIQ_REVEL_DATA_PATH at its output.
#
# The source CSV's columns (verified against REVEL's documented public
# release format) are:
#   chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
# AlleleIQ is GRCh38-only (see src/alleleiq/parsers/vcf_parser.py), so
# this script keys the tabix index on grch38_pos and drops hg19_pos.
#
# Usage: ./scripts/download_revel.sh [output_dir]
#   output_dir defaults to ./data

set -euo pipefail

REVEL_URL="https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
OUTPUT_DIR="${1:-./data}"
OUTPUT_BASENAME="revel_v1.3_all_chromosomes"
ZIP_ENTRY_NAME="revel_with_transcript_ids"

for tool in curl unzip bgzip tabix awk sort; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "error: required tool '${tool}' not found on PATH" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"
FINAL_TSV_GZ="${OUTPUT_DIR}/${OUTPUT_BASENAME}.tsv.gz"

if [[ -f "${FINAL_TSV_GZ}" && -f "${FINAL_TSV_GZ}.tbi" ]]; then
    echo "already present: ${FINAL_TSV_GZ} (and .tbi) — remove it first to rebuild"
    exit 0
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "downloading REVEL v1.3 (this is a large download, ~670MB)..."
curl -fL --progress-bar -o "${WORKDIR}/revel.zip" "${REVEL_URL}"

echo "extracting..."
unzip -q "${WORKDIR}/revel.zip" -d "${WORKDIR}"

CSV_PATH="${WORKDIR}/${ZIP_ENTRY_NAME}"
if [[ ! -f "${CSV_PATH}" ]]; then
    # Upstream has changed its archive layout; fail loudly rather than
    # silently indexing the wrong (or no) file.
    echo "error: expected entry '${ZIP_ENTRY_NAME}' not found after extraction" >&2
    echo "contents of ${WORKDIR}:" >&2
    ls -la "${WORKDIR}" >&2
    exit 1
fi

echo "transforming to (chrom, pos, ref, alt, revel_score) TSV, GRCh38 coordinates, position-sorted..."
# Column layout in the source CSV: chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
# REVEL lists one row per (position, ref, alt, transcript); multiple
# transcripts agree on the same score for a given substitution, so
# de-duplicating on (chrom, pos, ref, alt) after dropping the transcript
# column is safe and keeps the index small.
awk -F',' 'NR > 1 && $3 != "." {print $1"\t"$3"\t"$4"\t"$5"\t"$8}' "${CSV_PATH}" \
    | sort -u -k1,1 -k2,2n \
    > "${WORKDIR}/${OUTPUT_BASENAME}.tsv"

echo "bgzipping and tabix-indexing..."
bgzip -c "${WORKDIR}/${OUTPUT_BASENAME}.tsv" > "${FINAL_TSV_GZ}"
tabix -s 1 -b 2 -e 2 "${FINAL_TSV_GZ}"

echo "done: ${FINAL_TSV_GZ} (and .tbi)"
echo "set ALLELEIQ_REVEL_DATA_PATH=${FINAL_TSV_GZ} to use it"
