# Intended Use Statement

**Document ID:** ALQ-REG-001
**Product:** AlleleIQ
**Version:** 0.1.0-alpha

## Statement of intended use

AlleleIQ is a clinical decision support software tool intended to assist
qualified Clinical Genetic Scientists in the interpretation of germline
genomic variants identified through next-generation sequencing. AlleleIQ:

- Aggregates population frequency evidence (gnomAD v4), in silico
  predictive evidence (REVEL, SpliceAI via Ensembl VEP), and curated
  variant classification evidence (ClinVar) for a given variant.
- Applies a defined subset of ACMG/AMP evidence criteria — as codified
  by ACGS 2024 — using the Tavtigian et al. (2020) Bayesian point
  framework, producing a **preliminary** evidence-point total and
  five-tier classification band.
- Explicitly identifies and flags gaps in population reference-data
  coverage on a per-ancestry basis.
- Produces structured, auditable evidence reports (JSON and HTML) for
  review.

**AlleleIQ does not issue a clinical variant classification.** Every
report is explicitly labelled a preliminary evidence summary and
requires review, interpretation, and final classification by a
qualified Clinical Genetic Scientist, in full consideration of clinical
phenotype, family history, segregation data, and any other evidence
outside AlleleIQ's automated scope (see
[ACMG Framework](../architecture/acmg-framework.md) for exactly which
criteria are and are not automated).

## Intended users

Qualified Clinical Genetic Scientists and clinical genomics laboratory
personnel operating within an accredited clinical genomics laboratory
setting. AlleleIQ is not intended for direct use by patients or for
direct-to-consumer interpretation.

## Intended patient population

Individuals undergoing germline genomic sequencing for diagnostic or
predictive purposes, across all ancestries — with the explicit caveat
that interpretation confidence for patients whose ancestry is
under-represented in current reference databases may be lower, which
AlleleIQ is specifically designed to surface rather than obscure.

## Regulatory classification

AlleleIQ is classified as:

- **IVDR Class C**, under EU IVDR 2017/746 Annex VIII Rule 3(i)
  (software providing information used to detect a physiological or
  pathological state) and Rule 3(k) (software intended to provide
  information for detecting a predisposition to a medical condition).
- **UK MDR 2002, Part IV, Annex II, List B**.

This classification reflects that AlleleIQ provides decision-support
information that could influence a clinical decision, but does not
itself make an autonomous diagnostic determination. Every design
decision documented in [Safety Design](safety-design.md) reflects this
classification — the human clinical scientist always classifies;
AlleleIQ always supports.

## Limitations

See [ACMG Framework](../architecture/acmg-framework.md) for the specific
criteria out of scope, and the module docstrings in
`alleleiq.clients.vep`, `alleleiq.clients.clinvar`, and
`alleleiq.scoring.criteria.lof` for documented data-source and
methodology limitations (e.g. mitochondrial variants unsupported,
PVS1's transcript-specific-expression-region exclusion not checked, the
bundled ClinGen gene-disease validity list is illustrative rather than
exhaustive).
