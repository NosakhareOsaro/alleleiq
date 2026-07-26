# Getting Started

## Prerequisites

- Python 3.11+ (3.12 recommended — AlleleIQ's own development uses 3.12
  for the broadest wheel compatibility with `cyvcf2`/`pysam`).
- [uv](https://docs.astral.sh/uv/) for dependency management.
- Docker (optional, for containerised deployment).

## Installation

```bash
git clone https://github.com/NosakhareOsaro/alleleiq.git
cd alleleiq
./scripts/setup_dev.sh
```

This pins Python 3.12, installs all dependencies (including dev tools),
installs pre-commit hooks, and copies `.env.example` to `.env` if one
doesn't already exist. Review `.env` before running against real data —
every threshold and endpoint is configurable there.

## Your first analysis

AlleleIQ ships with a small sample VCF you can use immediately:

```bash
uv run alleleiq analyse tests/fixtures/vcf/single_snv.vcf --output-dir ./reports
```

This will:

1. Parse and normalise the VCF (`tests/fixtures/vcf/single_snv.vcf`
   contains one BRCA1 variant on GRCh38).
2. Query gnomAD v4, Ensembl VEP, and ClinVar for evidence.
3. Evaluate BA1, BS1, PM2, PP3, BP4, and PVS1.
4. Assess population coverage per ancestry group.
5. Accumulate a Bayesian evidence score and preliminary classification.
6. Write `./reports/<chrom>-<pos>-<ref>-<alt>.json` and `.html`.

You'll see a Rich-formatted progress bar and a colour-coded summary
table in your terminal, ending with the mandatory preliminary-evidence
disclaimer.

## Checking connectivity

```bash
uv run alleleiq version
```

Prints the installed version, configured gnomAD/VEP/ClinVar endpoints,
and live connectivity status for each — useful for confirming your
network/proxy configuration before running a real analysis.

## Setting up local REVEL scores (recommended)

PP3 and BP4 (computational evidence of a deleterious or benign effect)
depend on a REVEL score for missense variants. Ensembl VEP's REST API
*can* serve one, but its `dbNSFP` plugin coverage is incomplete — many
positions come back with no REVEL score at all, and PP3/BP4 then
degrade to `requires_clinical_input` not because the evidence is
genuinely absent, but because AlleleIQ couldn't retrieve it.

Building a local, tabix-indexed copy of REVEL's own precomputed-scores
dataset closes that gap — it has complete genome-wide coverage, and a
lookup takes microseconds:

```bash
./scripts/download_revel.sh
```

This downloads REVEL v1.3 (GRCh38 coordinates) from the upstream Mount
Sinai host, transforms it into a position-sorted TSV, and
bgzip/tabix-indexes it into `data/revel_v1.3_all_chromosomes.tsv.gz`
(plus a `.tbi` index alongside it).

!!! note "Disk and time"
    The source archive is roughly 670MB, and the extracted/transformed
    working files are several GB before compression — the script needs
    that much free space in a temporary directory (cleaned up
    automatically), plus the final index size (a few hundred MB) in
    your chosen output directory. It can take several minutes depending
    on your connection. This step is entirely optional and only needs
    to be run once per deployment host.

Then point AlleleIQ at it:

```bash
# in .env
ALLELEIQ_REVEL_DATA_PATH=./data/revel_v1.3_all_chromosomes.tsv.gz
```

If you skip this step, `ALLELEIQ_REVEL_DATA_PATH` is unset by default
and PP3/BP4 fall back to whatever REVEL score (if any) VEP REST
returns, and finally to `requires_clinical_input` if neither source has
one — the same behaviour AlleleIQ had before this feature existed, so
skipping it is safe, just less complete.

## ClinGen gene-disease validity and dosage sensitivity

PVS1 (null variant in a gene where loss-of-function is a known disease
mechanism) needs two independent facts about a gene: how strongly it is
established to cause *any* disease (ClinGen's Gene-Disease Validity
classification — Definitive, Strong, Moderate, Limited, ...), and
whether losing one copy of it is itself pathogenic (ClinGen's Dosage
Sensitivity Map haploinsufficiency score, 0–3). Unlike REVEL, both
datasets are small enough (~1.4MB combined) to ship bundled in this
repository at `data/clingen/` — there is nothing to download before
first use.

ClinGen re-curates genes continuously as expert panels review more
evidence, so the bundled snapshot should be refreshed periodically:

```bash
./scripts/update_clingen.sh
```

This re-downloads both files from ClinGen's own distribution endpoints
and overwrites `data/clingen/gene_disease_validity.csv` and
`data/clingen/dosage_sensitivity.tsv` in place. A weekly cadence is a
reasonable default; always review `git diff -- data/clingen` before
committing a refresh, since a classification change (e.g. a gene
downgraded from Definitive to Disputed) changes which PVS1 strength
tier that gene receives going forward. Also watch ClinGen's own
[curation activity notes](https://clinicalgenome.org/curation-activities/)
for changes that could affect how a *previously generated* AlleleIQ
report for that gene should be re-read.

If a gene has no ClinGen curation at all, PVS1 falls back to a small
bundled fallback list plus gnomAD's own pLI/LOEUF constraint metrics —
see `docs/architecture/acmg-framework.md` for the full tiering logic.

## Running the API locally

```bash
make run-api
# or directly:
uv run uvicorn alleleiq.api.main:app --reload --port 8000
```

Then visit `http://localhost:8000/docs` for interactive OpenAPI
documentation, or see the [API Reference](api-reference.md) guide.

## Running against your own VCF

```bash
uv run alleleiq analyse path/to/your.vcf \
    --output-dir ./reports \
    --format both \
    --genome-build GRCh38
```

!!! warning "GRCh38 only"
    AlleleIQ requires GRCh38 coordinates and will refuse to run against a
    VCF that declares (or fingerprints as) GRCh37 — coordinate systems
    are not interchangeable, and silently misinterpreting them would be
    far worse than refusing to proceed. Convert with a tool like
    `CrossMap` or `picard LiftoverVcf` first if needed.

## Next steps

- [Interpreting Reports](interpreting-reports.md) — how to read the
  evidence report you just generated.
- [API Reference](api-reference.md) — the REST API in full.
- [ACMG Framework](../architecture/acmg-framework.md) — what's being
  evaluated, and what isn't.
