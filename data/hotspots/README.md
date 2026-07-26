# Curated mutational hotspot regions

`curated_hotspots.json` maps a small set of gene symbols to amino-acid
position ranges with well-established recurrent pathogenic ("hotspot")
codons, used by PM1 (see `alleleiq/scoring/criteria/hotspot.py`) as one
of three independent signals that a variant falls in a critical
functional region (alongside live UniProt domain/active-site/binding-site
annotation and ClinVar variant-density detection).

**This is a curated illustrative snapshot, not a comprehensive or
authoritative automated feed.** Each entry reflects a defensible reading
of well-established genetics literature and genotype-phenotype
correlation at the time it was written, but:

- it was **not** independently re-verified against a live structural or
  curation database for every gene;
- coverage is intentionally narrow (6 genes) — a gene's absence from
  this file means "not yet curated here", never "this gene has no
  mutational hotspots";
- amino-acid ranges for broad domains (e.g. TP53's DNA-binding domain)
  are approximate; single-codon entries (e.g. RET M918, FGFR3 G380) are
  precise but still benefit from confirmation against current transcript
  annotation.

Consistent with AlleleIQ's broader stance (see
`docs/clinical/intended-use.md`), PM1 evaluations sourced from this file
are decision-support inputs for a qualified Clinical Genetic Scientist,
not a substitute for checking current ClinVar variant density, UniProt
annotation, or primary literature for the gene in question. Expert
reconfirmation before clinical use is expected, and is noted per-entry in
each region's `source` field. PM1 also always returns
`requires_clinical_input=True` regardless of outcome, for exactly this
reason.

## Extending the list

There is no automated script for this file (unlike
`data/gene_annotations/`'s `update_gene_annotations.py`) — adding a gene
here is a manual curation exercise requiring literature review, and the
file is small enough to edit directly. Follow the existing entries'
shape: `{"start_aa": int, "end_aa": int, "description": str, "source":
str}`, one list per gene symbol (uppercase).
