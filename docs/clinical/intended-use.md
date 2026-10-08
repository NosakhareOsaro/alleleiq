# Intended Use Statement

**Document ID:** AIQ-REG-001 (summary)
**Product:** AlleleIQ
**Document version:** 0.3 (draft)
**Software version:** 0.2.0

This page summarises the Intended Use Statement AIQ-REG-001 v0.3, which
is maintained as a separate controlled document. It describes what the
v0.2.0 software actually implements. Copy-number and structural
variants are out of scope at this version.

## Statement of intended use

AlleleIQ is a clinical decision support software tool intended to assist
qualified Clinical Genetic Scientists in the interpretation of germline
genomic variants identified through next-generation sequencing. AlleleIQ:

- Aggregates population frequency evidence (gnomAD v4), in silico
  predictive evidence (REVEL, from a local index or via Ensembl VEP, and
  SpliceAI via Ensembl VEP), curated variant classification evidence
  (ClinVar), gene-level curation (bundled ClinGen Gene-Disease Validity
  and Dosage Sensitivity snapshots), and protein domain annotations
  (UniProt) for a given variant.
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
  (devices intended for human genetic testing) and Rule 3(k) (devices
  intended for the management of patients suffering from a
  life-threatening disease or condition). InnoScot Health confirmed this
  classification basis in correspondence dated 19 June 2026; Rule 2 is
  not applicable. Rule 3(h) (screening, diagnosis or staging of cancer)
  also applies to the hereditary cancer genes; MDCG 2020-16 lists
  hereditary cancer predisposition tests under it. It gives the same
  class.
- **Great Britain:** expected to be a **general IVD** under UK MDR 2002
  Part IV, because AlleleIQ is not on the closed Annex II List A or
  List B. General IVDs are self-declared and registered with the MHRA,
  with no UK Approved Body. This is awaiting confirmation from InnoScot
  Health.

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
PVS1's transcript-specific-expression-region exclusion not checked, and
PVS1 falling back to a small bundled gene list or gnomAD constraint
where ClinGen has no curation for the gene).
