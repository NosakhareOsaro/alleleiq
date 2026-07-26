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
