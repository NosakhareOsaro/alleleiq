"""PM1: mutational-hotspot / critical-functional-domain ACMG/AMP criterion.

PM1 (Pathogenic Moderate) applies when a variant's amino-acid position
falls in a well-established mutational hotspot or functional domain
without benign variation. AlleleIQ checks three independent signals,
any one of which can make a position a PM1 candidate
(:attr:`~alleleiq.scoring.criteria.base.CriterionContext.pm1_domain_regions`,
resolved by merging :func:`load_curated_hotspots` with
:meth:`alleleiq.clients.uniprot.UniProtClient.fetch_domains`, and
:attr:`~alleleiq.scoring.criteria.base.CriterionContext.pm1_nearby_clinvar_matches`,
resolved by :meth:`alleleiq.clients.clinvar.ClinVarClient.search_same_position`
over a wider window than PS1/PM5 use):

1. A small bundled curated list of well-known hotspot codons/regions
   (:data:`~alleleiq.config.Settings.curated_hotspots_path`, see
   ``data/hotspots/README.md``).
2. A live UniProt-annotated domain, active site, or binding site
   (:mod:`alleleiq.clients.uniprot`).
3. ClinVar variant density: at least
   :attr:`~alleleiq.config.Settings.pm1_min_pathogenic_hotspot_count`
   independent, 2+-star Pathogenic/Likely Pathogenic records within
   :attr:`~alleleiq.config.Settings.pm1_hotspot_window_aa` amino acids.

Any candidate region is then checked for disqualifying benign variation:
more than :attr:`~alleleiq.config.Settings.pm1_max_allowed_benign_count`
Benign/Likely Benign ClinVar records nearby withholds PM1, per ACMG/AMP
guidance that recurrent benign variation in a region undermines treating
it as intolerant of variation.

**PM1 always returns ``requires_clinical_input=True``, regardless of
outcome** — including when it is not assessed at all. This is a
deliberate departure from PS1/PM5 and PP2/BP1's convention of
``requires_clinical_input=False`` for a pure "not assessed" result:
those criteria's absence of data is genuinely uninformative (no ClinVar
record nearby says nothing either way), whereas PM1's domain/hotspot
annotation coverage is *itself* incomplete for many genes in a way a
"not assessed" result cannot distinguish from "genuinely not a critical
region" — exactly the ambiguity the user's specification calls out
("domain annotations are not always complete for all genes"). Flagging
every PM1 result for review reflects that this uncertainty is different
in kind, not merely repeating the other criteria's caution.

Like PS1/PM5, "same position" is approximated by a symmetric genomic
window (not true codon-boundary mapping — see
:mod:`alleleiq.scoring.criteria.same_position`'s module docstring for
the same limitation), converted from an amino-acid window via an
approximate 3bp/codon factor.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from alleleiq.clients.clinvar import ClinVarPositionalMatch
from alleleiq.clients.uniprot import ProteinRegion
from alleleiq.logging import get_logger
from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength

from .base import CriterionContext, CriterionEvaluator
from .same_position import extract_protein_change

logger = get_logger(__name__)

_DATA_SOURCE = "Curated hotspot list (data/hotspots/) + UniProt domain annotation + ClinVar density"

_CODON_NUMBER_PATTERN = re.compile(r"(\d+)")

# Deliberately exact-match sets, not substring checks — see
# alleleiq.scoring.criteria.same_position for why ClinVar's free-text
# significance field makes substring matching unsafe here.
_PATHOGENIC_SIGNIFICANCE = frozenset(
    {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
)
_BENIGN_SIGNIFICANCE = frozenset({"benign", "likely benign", "benign/likely benign"})


@lru_cache(maxsize=8)
def load_curated_hotspots(path: str) -> dict[str, list[ProteinRegion]]:
    """Load and cache the curated mutational-hotspot region file.

    Args:
        path: Filesystem path to the hotspot JSON (typically
            ``settings.curated_hotspots_path``). Cached per distinct
            path for the process lifetime — the file is small and
            static, so re-parsing it per variant would be wasted work.

    Returns:
        ``{gene_symbol_upper: [ProteinRegion, ...]}``. Empty if the file
        is missing, unreadable, or malformed; callers treat an empty
        mapping identically to "gene not curated" for every gene.
    """
    file_path = Path(path)
    if not file_path.exists():
        logger.warning(
            "curated mutational hotspot file not found; PM1's curated-list "
            "signal will not be consulted for any gene",
            extra={"path": str(file_path)},
        )
        return {}

    try:
        raw = json.loads(file_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.error(
            "curated mutational hotspot file could not be parsed",
            extra={"path": str(file_path), "error": str(exc)},
        )
        return {}

    hotspots: dict[str, list[ProteinRegion]] = {}
    for gene_symbol, entries in raw.items():
        regions: list[ProteinRegion] = []
        for entry in entries:
            try:
                regions.append(ProteinRegion.model_validate(entry))
            except ValueError as exc:
                logger.error(
                    "skipping malformed curated hotspot entry",
                    extra={"gene_symbol": gene_symbol, "error": str(exc)},
                )
        if regions:
            hotspots[gene_symbol.upper()] = regions

    logger.info(
        "loaded curated mutational hotspot regions",
        extra={"path": str(file_path), "gene_count": len(hotspots)},
    )
    return hotspots


def _codon_number(protein_change: str | None) -> int | None:
    if not protein_change:
        return None
    match = _CODON_NUMBER_PATTERN.search(protein_change)
    return int(match.group(1)) if match else None


def _is_this_variant(
    match: ClinVarPositionalMatch, pos: int, ref: str, alt: str, protein_change: str | None
) -> bool:
    same = match.is_same_variant(pos, ref, alt)
    if same is None:
        return protein_change is not None and match.protein_change == protein_change
    return same


def _matching_region(codon: int, regions: list[ProteinRegion]) -> ProteinRegion | None:
    for region in regions:
        if region.start_aa <= codon <= region.end_aa:
            return region
    return None


def _is_established_pathogenic(match: ClinVarPositionalMatch, min_stars: int) -> bool:
    return (
        match.star_rating >= min_stars
        and match.clinical_significance.strip().lower() in _PATHOGENIC_SIGNIFICANCE
    )


def _is_benign(match: ClinVarPositionalMatch) -> bool:
    return match.clinical_significance.strip().lower() in _BENIGN_SIGNIFICANCE


class PM1Criterion(CriterionEvaluator):
    """PM1 (Pathogenic Moderate): variant lies in a mutational hotspot / critical domain."""

    code = ACMGCode.PM1
    default_strength = EvidenceStrength.MODERATE

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        annotation = context.vep_annotation
        protein_change = extract_protein_change(annotation.hgvs_p if annotation else None)
        codon = _codon_number(protein_change)

        if codon is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No protein-level (HGVSp) change available from VEP annotation; "
                    "PM1 cannot be assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
            )

        settings = context.settings
        coords = context.variant.coordinates
        # The variant's own ClinVar record falls inside its own search window;
        # counting it would let an existing classification manufacture its own
        # hotspot. A record without enough location detail to compare is
        # treated as the variant itself only if it carries the same protein change.
        matches = [
            match
            for match in context.pm1_nearby_clinvar_matches
            if not _is_this_variant(match, coords.pos, coords.ref, coords.alt, protein_change)
        ]
        min_stars = settings.clinvar_min_stars_strong_evidence
        benign_count = sum(1 for match in matches if _is_benign(match))
        pathogenic_count = sum(
            1 for match in matches if _is_established_pathogenic(match, min_stars)
        )

        matched_region = _matching_region(codon, context.pm1_domain_regions)
        density_hotspot = pathogenic_count >= settings.pm1_min_pathogenic_hotspot_count

        if matched_region is None and not density_hotspot:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"Codon {codon} is not within a curated hotspot region, a UniProt-"
                    "annotated domain/active-site/binding-site, or a ClinVar-derived "
                    f"variant-density hotspot (fewer than "
                    f"{settings.pm1_min_pathogenic_hotspot_count} established pathogenic "
                    "records nearby); PM1 not assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes=(
                    "Domain/hotspot annotation coverage is incomplete for many genes; "
                    "absence of a match here is not evidence against a critical region."
                ),
            )

        basis_parts = []
        if matched_region is not None:
            basis_parts.append(
                f"codon {codon} falls within {matched_region.description} "
                f"({matched_region.source})"
            )
        if density_hotspot:
            basis_parts.append(
                f"{pathogenic_count} established ({min_stars}+-star) Pathogenic/Likely "
                f"Pathogenic ClinVar records lie within {settings.pm1_hotspot_window_aa} "
                "amino acids"
            )
        region_basis = " and ".join(basis_parts)

        if benign_count > settings.pm1_max_allowed_benign_count:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"Codon {codon} is a PM1 candidate ({region_basis}), but {benign_count} "
                    "Benign/Likely Benign ClinVar record(s) in the same region are evidence "
                    "of tolerated benign variation there; PM1 not applied."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes="Disqualified by the benign-variation-in-region check (ACMG/AMP PM1 "
                "guidance).",
            )

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=True,
            evidence_basis=(
                f"Codon {codon} is a mutational hotspot / critical functional region "
                f"({region_basis}), with no disqualifying benign variation nearby."
            ),
            data_source=_DATA_SOURCE,
            requires_clinical_input=True,
            notes="Domain/hotspot annotation is not always complete for all genes; confirm "
            "region boundaries and mechanism manually before relying on this.",
        )
