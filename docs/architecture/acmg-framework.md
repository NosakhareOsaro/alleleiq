# ACMG/AMP Framework Implementation

AlleleIQ automates 8 of the ~28 ACMG/AMP criteria (Richards et al. 2015),
as adjusted by ACGS 2024, and scores them using the Tavtigian et al.
(2020) Bayesian point framework. This page explains what is implemented,
what is deliberately out of scope, and why.

## Implemented criteria

| Code | Module | Summary |
|---|---|---|
| **BA1** | `scoring/criteria/population_frequency.py` | Allele frequency > 5% (configurable) in *any adequately covered* gnomAD v4 population. Applied per-population, per ACGS 2024 — a high frequency observed only in an under-covered population does not trigger it. |
| **BS1** | `scoring/criteria/population_frequency.py` | Allele frequency above a disease-specific threshold (dominant/recessive). Requires `inheritance_mode` context; without it, marked `requires_clinical_input`. |
| **PM2** | `scoring/criteria/population_frequency.py` | Downgraded to **Supporting** strength by default per ACGS 2024. Assesses true absence (AC=0) in every adequately covered population; a variant observed at low-but-nonzero frequency is reported not-met with `requires_clinical_input=True` rather than an invented rarity cutoff. |
| **PP3 / BP4** | `scoring/criteria/computational.py` | REVEL ≥0.70/≥0.90 → Supporting/Moderate (PP3); ≤0.15 → Supporting (BP4). SpliceAI max delta ≥0.20/≥0.50 → Supporting/Moderate (PP3); <0.10 → Supporting (BP4). REVEL's absence is only flagged as a gap for missense variants (where it's expected); SpliceAI is considered for any consequence. |
| **PVS1** | `scoring/criteria/lof.py` | A simplified ACGS 2024 decision tree: LoF consequence type → gene-disease mechanism, tiered from ClinGen's bundled Gene-Disease Validity classification + Dosage Sensitivity Map haploinsufficiency score (Definitive/Strong + score 3 → full strength; Moderate + score 2 → fixed Strong; falls back to a small bundled gene list / gnomAD pLI-LOEUF when ClinGen has no curation for the gene) → exon position (last/penultimate exon downgrades Strong instead of Very Strong, reflecting possible NMD escape). |
| **PS1** | `scoring/criteria/same_position.py` | This variant's amino-acid change (from VEP's HGVSp) matches an established (≥`clinvar_min_stars_strong_evidence`-star) Pathogenic/Likely Pathogenic ClinVar record at the same genomic position. See "PS1 / PM5 limitations" below — this is comparative evidence, not a fully automated determination. |
| **PM5** | `scoring/criteria/same_position.py` | A *different* amino-acid change at the same codon has an established Pathogenic/Likely Pathogenic ClinVar classification. Same data source and limitations as PS1. |

## Why these eight

They are the criteria answerable from gnomAD, VEP, ClinVar, and gnomAD's
own constraint metrics alone — no clinical phenotype, family history,
segregation data, or functional assay results required. Every other
ACMG/AMP criterion (PS2–PS4, PM1, PM3, PM4, PM6, PP1, PP2, PP4, PP5,
BS2–BS4, BP1–BP3, BP5–BP7) needs evidence AlleleIQ has no data source
for, and is never silently assumed to be "not met" — it simply never
appears in `criteria_applied` or `criteria_not_met`, only in the
reviewing scientist's own workup.

## PS1 / PM5 limitations

PS1 and PM5 are automatable in the sense that ClinVar can be queried
programmatically, but the comparison itself has real boundaries that are
surfaced rather than hidden:

- **Position window, not codon boundary.** AlleleIQ does not perform
  strand-aware transcript-to-codon mapping. "Same position" is
  approximated by a symmetric genomic window
  (`Settings.same_position_window_bp`, default ±3bp) around this
  variant's VCF position, not the codon's exact three bases. This can
  occasionally miss a codon-spanning match, or — rarely — pull in an
  adjacent codon.
- **PS1's splice-divergence exclusion is best-effort.** ACGS guidance
  withholds PS1 when the new variant might create a splicing effect the
  matched pathogenic variant did not have. AlleleIQ has no per-ClinVar-hit
  splice annotation to compare against, so this is approximated using
  *this variant's own* SpliceAI signal (the same scores PP3/BP4 already
  consult) as a conservative proxy, not a true comparison of mechanism
  between the two variants.
- **Always `requires_clinical_input=True` when met.** Unlike BA1/PM2/PVS1,
  which can resolve to `requires_clinical_input=False` on a clean match,
  PS1 and PM5 never do — a same-position comparison always needs a human
  to confirm the two variants plausibly share a pathogenic mechanism
  before the match is relied upon.

## Scoring: `scoring/bayesian.py`

Point weights are keyed on **strength**, not code:

| Strength | Points |
|---|---|
| Very Strong | +8 |
| Strong | +4 |
| Moderate | +2 |
| Supporting | +1 |
| Benign Supporting | −1 |
| Benign Moderate | −2 *(no standard ACMG code uses this tier; included for symmetry)* |
| Benign Strong | −4 |
| Benign Standalone | −8 |

Classification bands: ≥10 Pathogenic · 6–9 Likely Pathogenic · 0–5 VUS ·
−1 to −6 Likely Benign · ≤−7 Benign.

`accumulate_evidence()` also:

- **Detects conflicts** — pathogenic-strength and benign-strength
  criteria both met is flagged (with a dedicated message when BA1 is
  involved, since standalone-benign alongside pathogenic evidence
  should essentially never occur and likely signals an annotation
  problem).
- **Derives confidence** from conflicts, unassessed-criteria volume, and
  population coverage completeness — not from the point total.
- **Computes an uncertainty band**: not a statistical confidence
  interval, but an explicit best-/worst-case range — "if every criterion
  currently marked `requires_clinical_input` resolved toward pathogenic
  (or benign), the score could reach X."

## Validation reality check

`scripts/validate_known_variants.py` runs this pipeline against 10 real
ClinVar expert-panel-reviewed variants. At the time of writing: 4/10
directly concordant, 6/10 conservative (AlleleIQ lands on VUS where
ClinVar has a definitive call), **0/10 dangerous** (no
Pathogenic↔Benign direction flips). The 6 conservative cases are not
failures — they are variants whose ClinVar classification depends on
PS3 (functional studies), PM1 (hotspot), or PP1 (segregation) evidence
AlleleIQ has no access to. See the module docstring in
`alleleiq.validation` for the full reasoning behind treating this as the
*expected and safe* outcome of a criteria subset this size.
