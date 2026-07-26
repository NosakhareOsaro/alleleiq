# Data Models

Every cross-module data exchange in AlleleIQ uses a Pydantic v2 model —
no raw `dict` crosses a public interface. This page documents the model
layer (`alleleiq.models`), organised by module.

## `alleleiq.models.variant`

| Model | Purpose |
|---|---|
| `GenomeBuild` | `GRCH38` / `GRCH37` (the latter exists so the parser can *detect and reject* it). |
| `VariantType` | `SNV`, `INDEL`, `CNV`, `STRUCTURAL`. |
| `GenomicCoordinates` | `chrom`, `pos`, `ref`, `alt`, `genome_build`. Validates chromosome naming (`chr17` → `17`, `chrM`/`chrMT` → `MT`) and allele characters (nucleotide or symbolic `<DEL>`-style). Frozen. |
| `Variant` | `coordinates`, `variant_type`, plus optional `gene_symbol`, `transcript_id`, `hgvs_c`, `hgvs_p`, `consequence`. Cross-validates that symbolic alleles pair with `CNV`/`STRUCTURAL` types and nucleotide alleles with `SNV`/`INDEL`. Frozen. |

## `alleleiq.models.population`

| Model | Purpose |
|---|---|
| `AncestryGroup` | The 8 gnomAD v4 population groups (`AFR`, `AMR`, `ASJ`, `EAS`, `FIN`, `NFE`, `SAS`, `OTH`), values matching gnomAD's lowercase API identifiers. `.display_name` gives the human-readable label used in reports. |
| `CoverageFlagLevel` | `WARNING` / `CRITICAL` / `NO_DATA`. |
| `OverallCoverageAssessment` | `COMPLETE` / `PARTIAL` / `SEVERELY_LIMITED`, the report-level roll-up. |
| `PopulationFrequency` | `ancestry`, `allele_count`, `allele_number`, `allele_frequency`, `homozygote_count`, `confidence_interval_95`. A validator forces `allele_frequency` to `None` whenever `allele_number == 0` — this is the model-level encoding of "no data" vs. "zero frequency observed," which matters for every downstream criterion. |
| `CoverageGap` | `ancestry`, `allele_number`, `threshold`, `flag_level`, `clinical_guidance_message` (a specific, actionable note naming supplementary resources). |

## `alleleiq.models.evidence`

| Model | Purpose |
|---|---|
| `EvidenceStrength` | `VERY_STRONG` … `SUPPORTING` (pathogenic side), `BENIGN_SUPPORTING` … `BENIGN_STANDALONE` (benign side). |
| `ACMGCode` | The full Richards et al. (2015) code set (`PVS1`, `PS1`–`PS4`, `PM1`–`PM6`, `PP1`–`PP5`, `BA1`, `BS1`–`BS4`, `BP1`–`BP7`). |
| `CriterionEvaluation` | `code`, `strength`, `met`, `evidence_basis`, `data_source`, `data_source_version`, `requires_clinical_input`, `notes`. Code and strength are independent fields *by design* — see [ACMG Framework](acmg-framework.md). |
| `PreliminaryClassification` | The five-tier band: `PATHOGENIC`, `LIKELY_PATHOGENIC`, `UNCERTAIN_SIGNIFICANCE`, `LIKELY_BENIGN`, `BENIGN`. |
| `ClassificationConfidence` | `HIGH` / `MODERATE` / `LOW`, derived from coverage completeness, unassessed-criteria volume, and conflicts — not from the point total itself. |
| `BayesianScore` | `evidence_points`, `criteria_applied`, `criteria_not_met`, `criteria_not_assessed`, `preliminary_classification`, `classification_confidence`, `evidence_points_confidence_interval`, `conflict_flags`. |

## `alleleiq.models.report`

| Model | Purpose |
|---|---|
| `CLINICAL_DISCLAIMER` | A module-level constant — the exact disclaimer text, typed as a `Literal` on `EvidenceReport.clinical_disclaimer` so no code path can construct a report with a missing or altered disclaimer. |
| `ReportMetadata` | `report_id` (UUID), `generated_at`, `alleleiq_version`, `database_versions`, `query_timestamps`. |
| `DatabaseQueryRecord` | One logged query: `database`, `query_description`, `parameters`, `timestamp`, `success`, `error_message`. |
| `AuditTrail` | `queries: list[DatabaseQueryRecord]` — every upstream call made while building a report, success or failure. |
| `CoverageSummary` | Report-level roll-up: which ancestries were sufficient/insufficient/absent, the overall assessment, and a prose clinical advisory. |
| `EvidenceReport` | The complete output: `variant`, `metadata`, `population_frequencies`, `coverage_gaps`, `coverage_summary`, `criterion_evaluations`, `bayesian_score`, `audit_trail`, `clinical_disclaimer`. This is the single source of truth both the JSON and HTML renderers derive from. |

## Design choices worth knowing

- **Code and strength are decoupled** (`CriterionEvaluation.code` /
  `.strength`) rather than encoded as combinatorial enum members like
  `PM2_SUPPORTING`. ACGS 2024 routinely adjusts a criterion's applied
  strength (PM2 downgraded to Supporting, PVS1 downgraded to Strong or
  Supporting) without changing which criterion fired.
- **`requires_clinical_input` is not "met=False."** A criterion can be
  `met=True` and still `requires_clinical_input=True` (e.g. PVS1 applied
  at Supporting strength because gene mechanism is unestablished) — it
  contributed points at a conservative strength, but a scientist should
  confirm the mechanism before relying on it further.
- **Everything frozen.** Every model uses `ConfigDict(frozen=True)` (or
  is immutable by construction) — an `EvidenceReport`, once built, is
  never mutated in place.
