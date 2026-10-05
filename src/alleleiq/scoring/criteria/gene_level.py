"""PP2 and BP1: gene-level missense-mechanism ACMG/AMP criteria.

Both criteria depend on a fact about the *gene*, not the variant itself:
whether missense variation is a recognised disease mechanism for it.

- **PP2** (Pathogenic Supporting): a missense variant in a gene with a
  low rate of benign missense variation, where missense variants are a
  common mechanism of disease.
- **BP1** (Benign Supporting): a missense variant in a gene where
  disease is established to be caused by loss-of-function/truncating
  variants specifically, making a missense finding less likely to be
  disease-causing by that gene's usual mechanism.

Both are resolved from a small curated local file
(:data:`~alleleiq.config.Settings.gene_annotations_path`, see
``data/gene_annotations/README.md``), loaded once per path and cached
for the process lifetime via :func:`load_gene_annotations` — the same
"local snapshot, not a live API" shape used by
:mod:`alleleiq.clients.clingen` and :mod:`alleleiq.clients.revel`. A
gene absent from the file is **not assessed** by either criterion; that
absence is a fact about curation coverage, never evidence of benignity
or pathogenicity, exactly as an absent ClinGen curation or REVEL score
degrades gracefully elsewhere in this codebase rather than being
silently treated as "no mechanism". A gene should never be curated
``pp2_applicable`` and ``bp1_applicable`` at once — that combination
would mean "missense is both a common and an implausible mechanism for
this gene", which is not a coherent claim the loader accepts.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

from alleleiq.logging import get_logger
from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength

from .base import CriterionContext, CriterionEvaluator
from .same_position import extract_protein_change

logger = get_logger(__name__)

_DATA_SOURCE = "Curated gene-level missense-mechanism annotations (data/gene_annotations/)"

# The ACMG/AMP definition of PP2/BP1 is specific to missense substitutions;
# in-frame indels and other "protein-altering" consequences change residue
# count or identity in ways the mechanism literature behind this curated
# list was not written to cover, so they are deliberately excluded rather
# than assumed to behave like a missense_variant.
_MISSENSE_CONSEQUENCE = "missense_variant"
_CODON_NUMBER_PATTERN = re.compile(r"(\d+)")


class ExcludedRegion(BaseModel):
    """An amino-acid range in which a gene-level criterion must not be applied."""

    model_config = ConfigDict(frozen=True)

    start_aa: int
    end_aa: int
    description: str


class GeneMechanismAnnotation(BaseModel):
    """One gene's curated missense-disease-mechanism annotation."""

    model_config = ConfigDict(frozen=True)

    pp2_applicable: bool
    bp1_applicable: bool
    mechanism: str
    source: str = Field(
        description="Citation-style rationale; see data/gene_annotations/README.md."
    )
    bp1_excluded_regions: list[ExcludedRegion] = Field(
        default_factory=list,
        description=(
            "Key functional domains where missense variants can be pathogenic even in a "
            "loss-of-function gene (e.g. ENIGMA's BRCA1/BRCA2 domains). BP1 is withheld "
            "for missense variants inside them."
        ),
    )


@lru_cache(maxsize=8)
def load_gene_annotations(path: str) -> dict[str, GeneMechanismAnnotation]:
    """Load and cache the curated gene-mechanism annotation file.

    Args:
        path: Filesystem path to the annotation JSON (typically
            ``settings.gene_annotations_path``). Cached per distinct
            path for the process lifetime — the file is small and
            static, so re-parsing it per variant would be wasted work.

    Returns:
        ``{gene_symbol_upper: GeneMechanismAnnotation}``. Empty if the
        file is missing, unreadable, or malformed; callers treat an
        empty mapping identically to "gene not curated" for every gene.
    """
    file_path = Path(path)
    if not file_path.exists():
        logger.warning(
            "gene-level missense-mechanism annotation file not found; "
            "PP2/BP1 will not be assessed for any gene",
            extra={"path": str(file_path)},
        )
        return {}

    try:
        raw = json.loads(file_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.error(
            "gene-level missense-mechanism annotation file could not be parsed",
            extra={"path": str(file_path), "error": str(exc)},
        )
        return {}

    annotations: dict[str, GeneMechanismAnnotation] = {}
    for gene_symbol, entry in raw.items():
        try:
            annotation = GeneMechanismAnnotation.model_validate(entry)
        except ValueError as exc:
            logger.error(
                "skipping malformed gene-level annotation entry",
                extra={"gene_symbol": gene_symbol, "error": str(exc)},
            )
            continue
        if annotation.pp2_applicable and annotation.bp1_applicable:
            logger.error(
                "skipping incoherent gene-level annotation entry: gene marked "
                "both pp2_applicable and bp1_applicable",
                extra={"gene_symbol": gene_symbol},
            )
            continue
        annotations[gene_symbol.upper()] = annotation

    logger.info(
        "loaded gene-level missense-mechanism annotations",
        extra={"path": str(file_path), "gene_count": len(annotations)},
    )
    return annotations


def _lookup(context: CriterionContext) -> tuple[str | None, GeneMechanismAnnotation | None]:
    annotation = context.vep_annotation
    gene_symbol = (annotation.gene_symbol if annotation else None) or context.variant.gene_symbol
    if not gene_symbol:
        return None, None
    by_gene = load_gene_annotations(context.settings.gene_annotations_path)
    return gene_symbol, by_gene.get(gene_symbol.upper())


def _is_missense(context: CriterionContext) -> bool:
    annotation = context.vep_annotation
    return annotation is not None and annotation.most_severe_consequence == _MISSENSE_CONSEQUENCE


def _codon(context: CriterionContext) -> int | None:
    annotation = context.vep_annotation
    protein_change = extract_protein_change(annotation.hgvs_p if annotation else None)
    match = _CODON_NUMBER_PATTERN.search(protein_change) if protein_change else None
    return int(match.group(1)) if match else None


def _bp1_withheld_reason(
    context: CriterionContext, gene_annotation: GeneMechanismAnnotation
) -> str | None:
    """Why BP1 should not be applied to this missense variant, or None if it may be.

    BP1's premise — missense is an implausible mechanism in this gene —
    fails for missense variants inside key functional domains, and for
    missense variants that are predicted to disrupt splicing (and so may
    act as loss of function).
    """
    annotation = context.vep_annotation
    if annotation is not None:
        splice_scores = [
            score
            for score in (
                annotation.spliceai_ds_ag,
                annotation.spliceai_ds_al,
                annotation.spliceai_ds_dg,
                annotation.spliceai_ds_dl,
            )
            if score is not None
        ]
        threshold = context.settings.spliceai_pp3_supporting_threshold
        if splice_scores and max(splice_scores) >= threshold:
            return (
                f"missense variant with a predicted splicing effect (max SpliceAI delta "
                f"{max(splice_scores):.2f}) may act as loss of function"
            )

    if not gene_annotation.bp1_excluded_regions:
        return None
    codon = _codon(context)
    if codon is None:
        return (
            "protein position unavailable, so it cannot be confirmed to lie outside the "
            "gene's key functional domains"
        )
    for region in gene_annotation.bp1_excluded_regions:
        if region.start_aa <= codon <= region.end_aa:
            return (
                f"codon {codon} lies within {region.description} (aa {region.start_aa}-"
                f"{region.end_aa}), where missense variants can be pathogenic"
            )
    return None


class PP2Criterion(CriterionEvaluator):
    """PP2 (Pathogenic Supporting): missense variant in a missense-mechanism gene."""

    code = ACMGCode.PP2
    default_strength = EvidenceStrength.SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        gene_symbol, gene_annotation = _lookup(context)

        if gene_symbol is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No gene symbol available from VEP annotation; PP2 cannot be assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        if not _is_missense(context):
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "PP2 applies only to missense_variant consequences; this variant's "
                    "consequence is not missense, so the gene's missense-mechanism status "
                    "is not applicable here."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        if gene_annotation is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"{gene_symbol} is not present in AlleleIQ's curated gene-level "
                    "missense-mechanism annotations; PP2 not assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
                notes=(
                    f"{gene_symbol} would benefit from expert curation of its missense-vs-"
                    "loss-of-function disease mechanism (see data/gene_annotations/README.md)."
                ),
            )

        if not gene_annotation.pp2_applicable:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"{gene_symbol}'s curated disease mechanism ({gene_annotation.mechanism}) "
                    "does not support PP2; missense is not established as a common disease "
                    "mechanism for this gene."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=True,
            evidence_basis=(
                f"{gene_symbol} has a curated missense disease mechanism "
                f"({gene_annotation.mechanism}); this variant is a missense_variant."
            ),
            data_source=_DATA_SOURCE,
            requires_clinical_input=False,
            notes=gene_annotation.source,
        )


class BP1Criterion(CriterionEvaluator):
    """BP1 (Benign Supporting): missense variant in a loss-of-function-only gene."""

    code = ACMGCode.BP1
    default_strength = EvidenceStrength.BENIGN_SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        gene_symbol, gene_annotation = _lookup(context)

        if gene_symbol is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No gene symbol available from VEP annotation; BP1 cannot be assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        if not _is_missense(context):
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "BP1 applies only to missense_variant consequences; this variant's "
                    "consequence is not missense, so the gene's loss-of-function-only status "
                    "is not applicable here."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        if gene_annotation is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"{gene_symbol} is not present in AlleleIQ's curated gene-level "
                    "missense-mechanism annotations; BP1 not assessed."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
                notes=(
                    f"{gene_symbol} would benefit from expert curation of its missense-vs-"
                    "loss-of-function disease mechanism (see data/gene_annotations/README.md)."
                ),
            )

        if not gene_annotation.bp1_applicable:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    f"{gene_symbol}'s curated disease mechanism ({gene_annotation.mechanism}) "
                    "does not support BP1; loss-of-function is not established as the "
                    "exclusive disease mechanism for this gene."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=False,
            )

        withheld = _bp1_withheld_reason(context, gene_annotation)
        if withheld is not None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=f"{gene_symbol}: {withheld}; BP1 withheld.",
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
                notes=gene_annotation.source,
            )

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=True,
            evidence_basis=(
                f"{gene_symbol} has a curated loss-of-function-only disease mechanism "
                f"({gene_annotation.mechanism}); this variant is a missense_variant, which "
                "is less likely to be disease-causing by this gene's established mechanism."
            ),
            data_source=_DATA_SOURCE,
            requires_clinical_input=False,
            notes=gene_annotation.source,
        )
