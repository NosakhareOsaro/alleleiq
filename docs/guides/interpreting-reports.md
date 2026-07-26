# Interpreting Reports

This guide walks through an AlleleIQ evidence report section by section,
for a reviewing Clinical Genetic Scientist reading either the HTML or
JSON output. Both formats carry identical content — HTML formats it for
human reading, JSON is for structured/LIMS consumption.

!!! danger "This is not a classification"
    Every field on this page describes *evidence*, not a verdict. The
    "preliminary classification" is a mechanically derived summary of
    the automated criteria only — see
    [ACMG Framework](../architecture/acmg-framework.md) for exactly
    which criteria that does and does not include.

## 1. Header and disclaimer

Report ID, generation timestamp, and AlleleIQ version, followed
immediately by the clinical disclaimer. This disclaimer is not
decorative — see [Safety Design](../clinical/safety-design.md) for how
it's enforced at the type level.

## 2. Variant summary

Gene, transcript, genomic coordinates, HGVS c./p. notation, and the most
severe consequence term — the identity of what's being evaluated.

## 3. Preliminary classification banner

The five-tier band (Pathogenic / Likely Pathogenic / Uncertain
Significance / Likely Benign / Benign), colour-coded, with three
additional numbers:

- **Evidence points** — the Tavtigian (2020) point total.
- **Confidence** (High/Moderate/Low) — reflects population coverage
  completeness and whether any criteria conflicted or couldn't be
  assessed, *not* the strength of the point total itself. A "Pathogenic"
  call at Low confidence is a real, meaningful combination: the criteria
  that fired were strong, but the surrounding evidence base has gaps
  worth knowing about before relying on it.
- **Uncertainty range** — if every criterion currently marked "requires
  clinical input" turned out, on review, to be met at its evaluated
  strength, how far could the score move in either direction? This is
  not a statistical confidence interval; it's a concrete "here's how
  much this could change" bound.

## 4. Conflict flags (if present)

Only appears when pathogenic-strength and benign-strength criteria were
both met — see [Safety Design §4](../clinical/safety-design.md). Treat
this as a prompt to check the underlying annotations before trusting
either side of the score.

## 5. Population frequency table

One row per gnomAD v4 ancestry group: allele count, allele number,
frequency, and 95% CI where available. Rows for under-covered
populations are highlighted amber (WARNING) or red (CRITICAL/NO_DATA)
with a specific guidance message. **Read this table before trusting a
BA1/BS1/PM2 result** — those criteria only draw on the rows *without* a
coverage flag; a gapped row means that population's data, whatever it
shows, wasn't used.

## 6. ACMG/AMP criteria evaluation

Two tables: criteria that were **applied** (contributed to the score,
with their evidence basis and data source) and criteria that were
**evaluated but not met** (checked, and confidently didn't apply).

## 7. Evidence not assessed — requires clinical input

Criteria that could not be conclusively automated — missing REVEL for a
missense variant, no inheritance-mode context for BS1, unestablished
gene-disease mechanism for PVS1, and so on — each with a specific note
on what's needed to resolve it. **This section is often the most
clinically important one on the page**: it's AlleleIQ telling you
exactly where its own evidence runs out.

## 8. Audit trail (collapsible)

Every upstream database query made while building this report —
database, query description, timestamp, and success/failure — plus the
evidence-source versions used (`gnomad_r4`, VEP release, etc.). Useful
for reproducibility and for understanding *why* a criterion came back
`requires_clinical_input` (e.g. a failed ClinVar query, visible here,
explains an empty ClinVar context elsewhere in the report).

## Common patterns to recognise

- **High confidence + few criteria applied**: complete population
  coverage and no conflicts, but the variant simply doesn't meet many
  automatable criteria — likely to end up VUS, and appropriately so.
- **Pathogenic-direction criteria + Low confidence**: strong automated
  signal, but check the coverage table and the "not assessed" section
  before treating the classification as settled — the gaps matter here.
- **Conservative VUS on a variant you expect to be pathogenic from other
  sources** (e.g. ClinVar): expected when the deciding evidence is
  functional, segregation, or hotspot-based — see
  [ACMG Framework](../architecture/acmg-framework.md)'s validation
  section for why this is the designed, safe behaviour of a 6-criterion
  automated subset, not a detection failure.
