#!/usr/bin/env bash
#
# Refreshes AlleleIQ's bundled ClinGen snapshots:
#   data/clingen/gene_disease_validity.csv   — Gene-Disease Validity curations
#   data/clingen/dosage_sensitivity.tsv      — Dosage Sensitivity Map
#
# Both files are small (~1.4MB combined) and are committed to this
# repository as static snapshots, unlike REVEL's multi-hundred-megabyte
# index (see scripts/download_revel.sh) — ClinGen has no live per-variant
# query API for either dataset, only downloadable full-table exports, so
# "querying ClinGen" in AlleleIQ means reading these bundled files (see
# alleleiq.clients.clingen.ClinGenClient), refreshed periodically by
# re-running this script.
#
# ClinGen curations change as expert panels review more genes and
# diseases; re-run this weekly (or before any clinically-relevant
# deployment) and review the diff before committing, since a
# classification change (e.g. a gene downgraded from Definitive to
# Disputed) can change which PVS1 strength tier past and future reports
# for that gene receive. Also check ClinGen's own release notes
# (https://clinicalgenome.org/curation-activities/) for classification
# changes that could affect the interpretation of prior AlleleIQ reports.
#
# Usage: ./scripts/update_clingen.sh [output_dir]
#   output_dir defaults to ./data/clingen

set -euo pipefail

GENE_VALIDITY_URL="https://search.clinicalgenome.org/kb/gene-validity/download"
DOSAGE_URL="https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
OUTPUT_DIR="${1:-./data/clingen}"

for tool in curl; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "error: required tool '${tool}' not found on PATH" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "downloading ClinGen Gene-Disease Validity curations..."
curl -fsSL -o "${WORKDIR}/gene_disease_validity.csv" "${GENE_VALIDITY_URL}"

# The download endpoint serves a genuine CSV directly, but it sits behind
# a web app that could plausibly start requiring a browser session in the
# future; fail loudly rather than silently bundling an HTML error page.
if ! head -c 64 "${WORKDIR}/gene_disease_validity.csv" | grep -q "CLINGEN GENE DISEASE VALIDITY"; then
    echo "error: unexpected content from ${GENE_VALIDITY_URL}" >&2
    echo "expected a CSV starting with 'CLINGEN GENE DISEASE VALIDITY CURATIONS'; got:" >&2
    head -c 200 "${WORKDIR}/gene_disease_validity.csv" >&2
    echo >&2
    exit 1
fi

echo "downloading ClinGen Dosage Sensitivity Map (GRCh38)..."
curl -fsSL -o "${WORKDIR}/dosage_sensitivity.tsv" "${DOSAGE_URL}"

if ! head -1 "${WORKDIR}/dosage_sensitivity.tsv" | grep -q "^#ClinGen Gene Curation Results"; then
    echo "error: unexpected content from ${DOSAGE_URL}" >&2
    head -c 200 "${WORKDIR}/dosage_sensitivity.tsv" >&2
    echo >&2
    exit 1
fi

mv "${WORKDIR}/gene_disease_validity.csv" "${OUTPUT_DIR}/gene_disease_validity.csv"
mv "${WORKDIR}/dosage_sensitivity.tsv" "${OUTPUT_DIR}/dosage_sensitivity.tsv"

echo "done: ${OUTPUT_DIR}/gene_disease_validity.csv, ${OUTPUT_DIR}/dosage_sensitivity.tsv"
echo "review \`git diff -- ${OUTPUT_DIR}\` for classification changes before committing"
