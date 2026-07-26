# Gene-level missense-mechanism annotations

`missense_constraint.json` maps a small set of gene symbols to whether a
missense variant is a recognised disease mechanism for that gene
(supporting `PP2`), whether the gene's disease mechanism is established
as loss-of-function/truncating specifically (supporting `BP1`), or
neither.

**This is a curated illustrative snapshot, not a comprehensive or
authoritative automated feed.** Each entry reflects a defensible reading
of well-established genetics literature and ClinGen-consistent gene
curation consensus at the time it was written, but:

- it was **not** independently re-verified against a live ClinGen
  Gene-Disease Validity or Dosage Sensitivity query for every gene;
- disease mechanism assignments can be revised as new evidence emerges
  (ClinGen gene curations are themselves periodically re-reviewed);
- coverage is intentionally narrow (50 genes) — a gene's absence from
  this file means "not yet curated here", never "missense is not a
  mechanism for this gene".

Consistent with AlleleIQ's broader stance (see
`docs/clinical/intended-use.md`), PP2 and BP1 evaluations sourced from
this file are decision-support inputs for a qualified Clinical Genetic
Scientist, not a substitute for checking the gene's current ClinGen
curation or primary literature. Expert reconfirmation before clinical
use is expected, and is noted per-gene in each entry's `source` field.

## Extending the list

Use `scripts/update_gene_annotations.py` to add or update entries. The
script does not derive mechanism assignments automatically — populating
it further is a manual curation exercise requiring literature or
ClinGen review, exactly like the initial 50 genes were curated.
