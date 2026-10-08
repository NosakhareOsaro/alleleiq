# Safety Design: Human-in-the-Loop Architecture

AlleleIQ's IVDR Class C classification (see
[Intended Use](intended-use.md)) is not a label bolted on after the
fact — it shapes concrete design decisions throughout the codebase.
This page documents the specific mechanisms that keep a qualified
Clinical Genetic Scientist in the loop for every result.

## 1. The disclaimer cannot be constructed away

`alleleiq.models.report.CLINICAL_DISCLAIMER` is a fixed constant, and
`EvidenceReport.clinical_disclaimer` is typed as a `Literal` matching it
exactly. No code path — CLI, API, or any future integration — can build
an `EvidenceReport` with a missing, blank, or edited disclaimer without
failing Pydantic validation. `tests/unit/test_report_models.py` asserts
this cannot be overridden. The disclaimer appears in both the JSON and
HTML report (prominently, red-bordered, in the HTML) and in the CLI's
terminal summary.

## 2. "Requires clinical input" is a first-class outcome, not an error

Every `CriterionEvaluator.evaluate()` is guaranteed never to raise (see
[ACMG Framework](../architecture/acmg-framework.md)). When a criterion
cannot be conclusively assessed — missing REVEL for a missense variant,
no inheritance-mode context for BS1, an unestablished gene-disease
mechanism for PVS1 — the result is `requires_clinical_input=True`, which
routes it into the report's dedicated "requires clinical input" section
(both JSON `evidence_evaluation.criteria_not_assessed` and a highlighted
HTML section), rather than being silently omitted or defaulted to "not
met."

## 3. Population coverage gaps are never averaged away

See [Population Model](../architecture/population-model.md) in full.
The short version: a criterion is only evaluated against ancestry groups
with adequate gnomAD coverage, and every gap is surfaced with a specific
guidance message — never folded into a single pooled frequency that
would hide which patients' data this variant's evidence gap most
affects.

## 4. Conflicting evidence is flagged, not resolved silently

`scoring/bayesian.py`'s `_detect_conflicts()` explicitly checks whether
pathogenic-strength and benign-strength criteria were both met for the
same variant — a situation that should essentially never occur and, when
it does, likely indicates an annotation problem worth a scientist's
attention rather than a number the point-sum can quietly average out.

## 5. Every automated failure is visible, never silent absence of evidence

`alleleiq.exceptions.DatabaseQueryError` and its subclasses
(`GnomADError`, `VEPError`, `ClinVarError`) are distinguished from a
legitimate "no data" result at every layer. In `pipeline.py`, a failed
query is recorded in the report's `audit_trail` with `success=False` and
the error message — it never gets silently treated as "the variant is
absent from this database," which would be indistinguishable from a
real, confidently-absent result and could bias a PM2 (absence) call.

## 6. Logging never carries patient-identifying genomic content

`alleleiq.logging`'s structured JSON logger is documented (and, in
review, enforced) to carry only process metadata — variant coordinates,
counts, timings, request IDs — at INFO level and above. Allele
frequencies, classifications, and criterion outcomes are logged only at
DEBUG, which is not expected to be enabled against real patient data in
any deployed environment.

## 7. Uncertainty is quantified, not just labelled

`BayesianScore.classification_confidence` (HIGH/MODERATE/LOW) and
`evidence_points_confidence_interval` (an explicit best-/worst-case
range if unassessed criteria resolved either direction) give the
reviewing scientist a concrete sense of how much automated confidence to
place in the preliminary classification — not just the classification
itself.

## What this does not do

None of the above makes AlleleIQ a classifier. It makes AlleleIQ an
honest evidence aggregator: one that tells a Clinical Genetic Scientist
exactly what it knows, what it doesn't, and where the gaps in its
knowledge specifically fall — so the human decision that follows is
better informed, not replaced.
