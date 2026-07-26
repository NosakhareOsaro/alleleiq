# Population Stratification Methodology

This page documents `alleleiq.scoring.coverage`, the module responsible
for AlleleIQ's core clinical-equity claim: population coverage gaps are
never suppressed, hidden, or averaged away.

## The problem with a single frequency number

A variant reported as "allele frequency 0.0001" from gnomAD's aggregate
count could mean genuinely rare across all ancestries, or it could mean
absent from a population that was never adequately sampled in the first
place — two very different facts for a clinical interpretation. Averaging
across gnomAD's 8 ancestry groups before applying a criterion like PM2
(absence) or BA1 (frequency too high to be pathogenic) hides exactly the
distinction that matters.

AlleleIQ never performs that averaging. Every population-frequency
criterion (`scoring/criteria/population_frequency.py`) evaluates each
`AncestryGroup` independently, drawing only on populations *without* a
coverage gap.

## Coverage assessment (`scoring/coverage.py`)

For each of gnomAD v4's 8 ancestry groups, `assess_coverage()` compares
the observed allele number (AN) against a configurable minimum threshold:

| Ancestry | Default threshold (AN) |
|---|---|
| AFR (African/African American) | 10,000 |
| AMR (Latino/Admixed American) | 5,000 |
| ASJ (Ashkenazi Jewish) | 3,000 |
| EAS (East Asian) | 5,000 |
| FIN (Finnish) | 3,000 |
| NFE (Non-Finnish European) | 20,000 |
| SAS (South Asian) | 10,000 |
| OTH (Other/Unassigned) | 1,000 |

These are prototype defaults informed by gnomAD v4 cohort sizes
(`ALLELEIQ_COVERAGE_THRESHOLD_*` in `.env.example`), not a clinically
validated standard — revisiting them is part of the planned formal
validation study.

Each population below threshold gets a severity flag:

- **NO_DATA** — `allele_number == 0` (or the ancestry is absent from
  gnomAD's response entirely — treated identically).
- **CRITICAL** — below `threshold × coverage_critical_fraction` (default
  50%).
- **WARNING** — below threshold but above the critical fraction.

Every flagged population gets a specific, actionable
`clinical_guidance_message` naming a supplementary resource — H3Africa
for AFR, GenomeAsia 100K for SAS, ChinaMAP/ToMMo for EAS, and so on (see
`_SUPPLEMENTARY_RESOURCE_HINTS` in `coverage.py`).

## Report-level roll-up

`CoverageSummary` classifies the variant's *overall* coverage:

- **COMPLETE** — no gaps.
- **PARTIAL** — some gaps, but at least one population adequately
  covered.
- **SEVERELY_LIMITED** — ≥75% of ancestry groups have a gap.

`SEVERELY_LIMITED` coverage directly forces `ClassificationConfidence.LOW`
in the Bayesian accumulator (see [ACMG Framework](acmg-framework.md)) —
a variant that could only be assessed against one or two well-covered
populations should never present as a high-confidence classification,
regardless of how clean the criteria that *did* fire look.

## How this reaches the report

Both the JSON and HTML reports surface `population_frequencies`,
`coverage_gaps`, and `coverage_summary` as first-class sections — not
buried in an appendix. The HTML report highlights gapped rows in
amber (WARNING) or red (CRITICAL/NO_DATA) directly in the per-ancestry
frequency table, so a reviewing scientist sees the gap in the same
glance as the frequency itself.
