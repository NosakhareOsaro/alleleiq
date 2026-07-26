"""PS1 and PM5: same-position pathogenicity ACMG/AMP criteria.

Both criteria compare this variant's amino-acid change against other
*classified* variants ClinVar reports near the same genomic position
(:attr:`~alleleiq.scoring.criteria.base.CriterionContext.same_position_clinvar_matches`,
resolved by :meth:`alleleiq.clients.clinvar.ClinVarClient.search_same_position`):

- **PS1** (Pathogenic Strong): an established Pathogenic/Likely Pathogenic
  ClinVar variant produces the *same* amino-acid change as this variant.
- **PM5** (Pathogenic Moderate): an established Pathogenic/Likely
  Pathogenic ClinVar variant affects the *same codon* but produces a
  *different* amino-acid change.

Both automation boundaries are real and are surfaced rather than hidden:

1. **Position window, not codon boundary.** AlleleIQ does not perform
   strand-aware transcript-to-codon mapping, so "same position" is
   approximated by a symmetric genomic window
   (:attr:`~alleleiq.config.Settings.same_position_window_bp`) around this
   variant's VCF position, not the codon's exact three bases. This can
   occasionally miss a codon-spanning match or (rarely) include an
   adjacent codon.
2. **PS1's splice-divergence exclusion is best-effort.** ACMG/ACGS
   guidance withholds PS1 when the new variant might create a splice
   effect the matched pathogenic variant did not have. AlleleIQ has no
   per-ClinVar-hit splice annotation to compare against, so this is
   approximated using *this variant's own* SpliceAI signal (the same
   scores PP3/BP4 already use) as a proxy: a high SpliceAI delta on the
   variant under review withholds PS1 pending manual review, even when
   the amino-acid change matches.
3. **Every applied match still requires clinical sign-off.** PS1 and PM5
   are always returned with ``requires_clinical_input=True`` when met, per
   ACGS 2024 guidance that same-position comparisons need a human to
   confirm the two variants plausibly share a mechanism — this is never
   fully automatable, regardless of how confident the ClinVar match looks.
"""

from __future__ import annotations

import re

from alleleiq.clients.clinvar import ClinVarPositionalMatch
from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength

from .base import CriterionContext, CriterionEvaluator

_DATA_SOURCE = "Ensembl VEP (protein change) / ClinVar (same-position classified variants)"

# Matches both VEP's transcript-prefixed HGVSp ("ENSP00000269305.4:p.Arg175His")
# and ClinVar's parenthesised VariationName form ("... (p.Arg175His)").
_PROTEIN_CHANGE_PATTERN = re.compile(r"p\.([A-Za-z0-9_*=]+)")
_CODON_NUMBER_PATTERN = re.compile(r"(\d+)")

# Deliberately an exact-match set, not a substring check: ClinVar's free-text
# significance field includes compound values like "Conflicting interpretations
# of pathogenicity" that contain the substring "pathogenicity" without meaning
# the variant is pathogenic. Anything outside this set is treated as "not an
# established pathogenic match" rather than guessed at.
_PATHOGENIC_SIGNIFICANCE = frozenset(
    {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
)


def extract_protein_change(hgvs_p: str | None) -> str | None:
    """Extract the bare protein-change token (e.g. ``'Arg175His'``) from an HGVSp string.

    Args:
        hgvs_p: A transcript-prefixed HGVSp string (VEP) or a ClinVar
            ``VariationName``-style string containing a parenthesised
            protein change.

    Returns:
        The bare change token, or None if no ``p.`` change could be found.
    """
    if not hgvs_p:
        return None
    match = _PROTEIN_CHANGE_PATTERN.search(hgvs_p)
    return match.group(1) if match else None


def _codon_number(protein_change: str | None) -> int | None:
    if not protein_change:
        return None
    match = _CODON_NUMBER_PATTERN.search(protein_change)
    return int(match.group(1)) if match else None


def _is_established_pathogenic(match: ClinVarPositionalMatch, min_stars: int) -> bool:
    return (
        match.star_rating >= min_stars
        and match.clinical_significance.strip().lower() in _PATHOGENIC_SIGNIFICANCE
    )


def _is_any_pathogenic(match: ClinVarPositionalMatch) -> bool:
    return match.clinical_significance.strip().lower() in _PATHOGENIC_SIGNIFICANCE


def _splice_signal_present(context: CriterionContext) -> bool:
    """Best-effort proxy for "this variant might have a distinct splicing effect".

    See module docstring point 2: this checks the *querying* variant's own
    SpliceAI delta (already computed for PP3/BP4), not the matched
    ClinVar record's mechanism, which AlleleIQ has no way to inspect.
    """
    annotation = context.vep_annotation
    if annotation is None:
        return False
    scores = [
        annotation.spliceai_ds_ag,
        annotation.spliceai_ds_al,
        annotation.spliceai_ds_dg,
        annotation.spliceai_ds_dl,
    ]
    available = [score for score in scores if score is not None]
    if not available:
        return False
    return max(available) >= context.settings.spliceai_pp3_supporting_threshold


class PS1Criterion(CriterionEvaluator):
    """PS1 (Pathogenic Strong): same amino-acid change as an established pathogenic variant."""

    code = ACMGCode.PS1
    default_strength = EvidenceStrength.STRONG

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        annotation = context.vep_annotation
        protein_change = extract_protein_change(annotation.hgvs_p if annotation else None)

        if protein_change is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No protein-level (HGVSp) change available from VEP annotation; "
                    "PS1 cannot be assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
            )

        matches = context.same_position_clinvar_matches
        if not matches:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"No classified ClinVar variants found near this genomic position to "
                    f"compare against {protein_change}."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
                notes="Not assessed: absence of same-position ClinVar data is not evidence "
                "of benignity.",
            )

        min_stars = context.settings.clinvar_min_stars_strong_evidence
        same_change_hits = [m for m in matches if m.protein_change == protein_change]
        established = [m for m in same_change_hits if _is_established_pathogenic(m, min_stars)]

        if not established:
            weak_hits = [m for m in same_change_hits if _is_any_pathogenic(m)]
            if weak_hits:
                weakest = weak_hits[0]
                return CriterionEvaluation(
                    code=self.code,
                    strength=self.default_strength,
                    met=False,
                    evidence_basis=(
                        f"{protein_change} matches a {weakest.clinical_significance} ClinVar "
                        f"entry ({weakest.accession}), but its review status "
                        f"({weakest.star_rating} star(s)) falls below the {min_stars}-star "
                        "strong-evidence threshold; PS1 not applied."
                    ),
                    data_source=_DATA_SOURCE,
                    requires_clinical_input=False,
                    notes="A higher-confidence resubmission or additional evidence review "
                    "may change this outcome.",
                )
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"No established Pathogenic/Likely Pathogenic ClinVar record shares "
                    f"{protein_change}; PS1 not applicable."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        best = established[0]
        if _splice_signal_present(context):
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"{protein_change} matches an established {best.clinical_significance} "
                    f"ClinVar variant ({best.accession}, {best.star_rating} stars), but this "
                    "variant's own SpliceAI signal suggests a possible splicing effect the "
                    "matched variant may not share; PS1 withheld pending manual review."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes="Excluded by the splice-impact-divergence check (best-effort; see "
                "alleleiq.scoring.criteria.same_position module docs); confirm mechanism "
                "manually.",
            )

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=True,
            evidence_basis=(
                f"{protein_change} matches {best.accession}, classified "
                f"{best.clinical_significance} with {best.star_rating} stars in ClinVar — "
                "the same amino-acid change as an established pathogenic variant."
            ),
            data_source=_DATA_SOURCE,
            requires_clinical_input=True,
            notes="Automated same-amino-acid-change match; confirm no evidence of differing "
            "functional impact before relying on this.",
        )


class PM5Criterion(CriterionEvaluator):
    """PM5 (Pathogenic Moderate): different amino-acid change at an established pathogenic codon."""

    code = ACMGCode.PM5
    default_strength = EvidenceStrength.MODERATE

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        annotation = context.vep_annotation
        protein_change = extract_protein_change(annotation.hgvs_p if annotation else None)
        codon = _codon_number(protein_change)

        if protein_change is None or codon is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No protein-level (HGVSp) change available from VEP annotation; "
                    "PM5 cannot be assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
            )

        matches = context.same_position_clinvar_matches
        if not matches:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"No classified ClinVar variants found near this genomic position to "
                    f"compare against {protein_change}."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
                notes="Not assessed: absence of same-position ClinVar data is not evidence "
                "of benignity.",
            )

        min_stars = context.settings.clinvar_min_stars_strong_evidence
        different_change_hits = [
            m
            for m in matches
            if m.protein_change is not None
            and m.protein_change != protein_change
            and _codon_number(m.protein_change) == codon
        ]
        established = [m for m in different_change_hits if _is_established_pathogenic(m, min_stars)]

        if not established:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"No different amino-acid change at codon {codon} has an established "
                    f"({min_stars}+ star) Pathogenic/Likely Pathogenic ClinVar classification; "
                    "PM5 not applicable."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        best = established[0]
        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=True,
            evidence_basis=(
                f"Codon {codon} carries an established {best.clinical_significance} ClinVar "
                f"variant with a different amino-acid change ({best.protein_change}, "
                f"{best.accession}, {best.star_rating} stars); this variant instead produces "
                f"{protein_change}."
            ),
            data_source=_DATA_SOURCE,
            requires_clinical_input=True,
            notes=(
                "PM5 cannot confirm this variant is equally or more damaging than the "
                "established one without functional evidence; always requires clinical "
                "review before use."
            ),
        )
