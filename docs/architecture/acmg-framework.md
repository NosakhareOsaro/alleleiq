# ACMG/AMP Framework Implementation

AlleleIQ automates 6 of the ~28 ACMG/AMP criteria (Richards et al. 2015),
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

## Why these six

They are the criteria answerable from gnomAD, VEP, and gnomAD's own
constraint metrics alone — no clinical phenotype, family history,
segregation data, or functional assay results required. Every other
ACMG/AMP criterion (PS1–PS4, PM1, PM3–PM6, PP1, PP2, PP4, PP5, BS2–BS4,
BP1–BP3, BP5–BP7) needs evidence AlleleIQ has no data source for, and is
never silently assumed to be "not met" — it simply never appears in
`criteria_applied` or `criteria_not_met`, only in the reviewing
scientist's own workup.

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
