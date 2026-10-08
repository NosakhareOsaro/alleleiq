# AlleleIQ

**Population-inclusive evidence aggregation for clinical genomic variant interpretation.**

!!! danger "AlleleIQ does not classify variants"
    AlleleIQ aggregates evidence, applies ACMG/AMP and ACGS 2024 criteria in
    a structured and auditable way, and presents the result — including
    where the evidence is thin or absent — for interpretation by a
    qualified Clinical Genetic Scientist. The human always classifies;
    AlleleIQ always supports.

## The clinical problem

Around 80% of participants in genomic studies are of European ancestry
(Popejoy & Fullerton, *Nature* 2016), and the reference data underpinning
clinical genomic interpretation — population allele frequencies,
functional annotations, curated variant classifications — carries the
same skew. As a consequence, patients from African, South Asian,
East Asian, and Latino backgrounds receive a Variant of Uncertain
Significance (VUS) result at 1.5–2x the rate of European patients. This is
not because their variants are inherently more ambiguous — it is because
the reference data needed to interpret them with confidence is
systematically thinner or absent.

We are not aware of a commercially available interpretation tool that
routinely reports per-ancestry evidence coverage or surfaces these gaps
explicitly to the reviewing scientist. AlleleIQ does both, treating population coverage
completeness as a first-class output alongside the evidence itself.

## What AlleleIQ does

1. **Aggregates evidence** from gnomAD v4 (per-ancestry allele
   frequencies), Ensembl VEP (consequence, SpliceAI), a local REVEL index
   (with VEP as fallback), ClinVar (curated classifications, weighted by
   review-status star rating), bundled ClinGen gene-disease validity and
   dosage-sensitivity snapshots, and UniProt protein domain annotations.
2. **Applies eleven automatable ACMG/AMP criteria** — PVS1, BA1, BS1, PM2,
   PP3, BP4, PS1, PM5, PP2, BP1, and PM1 — per the Tavtigian et al. (2020)
   Bayesian point framework, as adopted by ACGS 2024.
3. **Assesses population coverage** for every gnomAD v4 ancestry group
   independently, flagging under-covered populations rather than
   silently averaging them into an overall frequency.
4. **Produces a structured evidence report** (JSON and HTML) with an
   explicit preliminary classification, an uncertainty band, and a
   mandatory clinical disclaimer — for review, not for acceptance.

See [Architecture Overview](architecture/overview.md) for how these
pieces fit together, or [Getting Started](guides/getting-started.md) to
run your first analysis.

## Regulatory context

AlleleIQ is classified as **IVDR Class C** under EU IVDR 2017/746 Rule
3(i) and Rule 3(k) (Rule 3(h) also applies, with the same class). In
Great Britain it is expected to be a **general IVD** under UK MDR 2002
Part IV, pending confirmation. See
[Intended Use](clinical/intended-use.md) and
[Safety Design](clinical/safety-design.md) for the full regulatory and
human-in-the-loop safety rationale.
