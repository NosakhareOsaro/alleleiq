"""Abstract base class and shared evidence context for ACMG/AMP criterion evaluators.

Each criterion (BA1, PM2, PVS1, ...) is implemented as a small,
self-contained, independently testable :class:`CriterionEvaluator`. The
public :meth:`CriterionEvaluator.evaluate` entry point is guaranteed
never to raise: any unexpected failure in a subclass's evaluation logic
is caught and converted into a :class:`~alleleiq.models.evidence.CriterionEvaluation`
with ``requires_clinical_input=True`` rather than aborting the whole
report. A criterion that *cannot* be assessed for lack of data is not a
system failure — it is a fact about the evidence that must be surfaced
to the reviewing scientist, exactly like a population coverage gap.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import ClassVar, Literal

from alleleiq.clients.clinvar import ClinVarPositionalMatch, ClinVarSubmission
from alleleiq.clients.gnomad import GeneConstraint
from alleleiq.clients.uniprot import ProteinRegion
from alleleiq.clients.vep import VEPAnnotation
from alleleiq.config import Settings, get_settings
from alleleiq.logging import get_logger
from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength
from alleleiq.models.population import AncestryGroup, CoverageGap, PopulationFrequency
from alleleiq.models.variant import Variant

logger = get_logger(__name__)


@dataclass(frozen=True)
class CriterionContext:
    """The bundle of evidence available when evaluating criteria for one variant.

    Built once per variant by the scoring orchestrator and passed to every
    :class:`CriterionEvaluator`; individual evaluators read only the
    fields relevant to their criterion.
    """

    variant: Variant
    population_frequencies: dict[AncestryGroup, PopulationFrequency] = field(default_factory=dict)
    coverage_gaps: dict[AncestryGroup, CoverageGap] = field(default_factory=dict)
    vep_annotation: VEPAnnotation | None = None
    gene_constraint: GeneConstraint | None = None
    clinvar_submission: ClinVarSubmission | None = None
    same_position_clinvar_matches: list[ClinVarPositionalMatch] = field(
        default_factory=list,
        metadata={
            "description": "Classified ClinVar records within settings.same_position_window_bp "
            "of this variant's genomic position, resolved by ClinVarClient.search_same_position. "
            "Consumed by PS1/PM5 (see alleleiq.scoring.criteria.same_position). Empty if not "
            "consulted (e.g. ClinVar disabled) or if ClinVar has no classified records nearby."
        },
    )
    revel_score: float | None = field(
        default=None,
        metadata={
            "description": "REVEL score resolved by REVELClient (local tabix index, "
            "preferred, or VEP REST fallback). Takes precedence over "
            "vep_annotation.revel_score for PP3/BP4 when set; falls back to "
            "vep_annotation.revel_score when REVELClient was not consulted."
        },
    )
    revel_source: str | None = None
    clingen_classification: str | None = field(
        default=None,
        metadata={
            "description": "Best ClinGen gene-disease-validity classification for this "
            "variant's gene (e.g. 'Definitive'), resolved by ClinGenClient. None if the "
            "gene has no ClinGen curation, or ClinGen was not consulted (e.g. a test "
            "constructing CriterionContext directly) — PVS1 falls back to its bundled "
            "gene list / gnomAD constraint check in either case."
        },
    )
    clingen_haploinsufficiency_score: int | None = field(
        default=None,
        metadata={
            "description": "ClinGen dosage-map haploinsufficiency score (0-3) for this "
            "variant's gene, resolved by ClinGenClient. None if not consulted; 0 if "
            "consulted but the gene is absent from ClinGen's dosage map or carries one "
            "of its two non-numeric special codes (see alleleiq.clients.clingen)."
        },
    )
    clingen_recessive_gene: bool = field(
        default=False,
        metadata={
            "description": "True when ClinGen's dosage map codes this gene 30 (gene "
            "associated with an autosomal recessive phenotype). Haploinsufficiency scores "
            "do not apply to such genes; PVS1 uses this to recognise biallelic loss of "
            "function as the established mechanism."
        },
    )
    pm1_domain_regions: list[ProteinRegion] = field(
        default_factory=list,
        metadata={
            "description": "Merged curated-hotspot + UniProt-annotated domain/active-site/"
            "binding-site regions for this variant's gene, resolved by "
            "alleleiq.pipeline._fetch_pm1_evidence. Consumed by PM1 (see "
            "alleleiq.scoring.criteria.hotspot). Empty if not consulted or if neither "
            "source has data for this gene."
        },
    )
    pm1_nearby_clinvar_matches: list[ClinVarPositionalMatch] = field(
        default_factory=list,
        metadata={
            "description": "Classified ClinVar records within settings.pm1_hotspot_window_aa "
            "(converted to an approximate genomic window) of this variant's genomic position, "
            "resolved by ClinVarClient.search_same_position with a wider window than PS1/PM5 "
            "use. Consumed by PM1's variant-density hotspot detection and disqualifying-"
            "benign-variation check."
        },
    )
    inheritance_mode: Literal["dominant", "recessive"] | None = field(
        default=None,
        metadata={
            "description": "Disease inheritance pattern, required for BS1's "
            "disease-specific frequency threshold. None if not supplied by the caller."
        },
    )
    settings: Settings = field(default_factory=get_settings)


class CriterionEvaluator(ABC):
    """Base class for a single ACMG/AMP criterion evaluator.

    Subclasses declare ``code`` (which ACMG code they evaluate) and
    ``default_strength`` (the strength used both when no adjustment
    applies and as the fallback strength if evaluation fails
    unexpectedly), then implement :meth:`_evaluate`.
    """

    code: ClassVar[ACMGCode]
    default_strength: ClassVar[EvidenceStrength]

    def evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        """Evaluate this criterion for the given variant. Never raises.

        Args:
            context: The evidence bundle for the variant under evaluation.

        Returns:
            A :class:`CriterionEvaluation`. If the subclass's evaluation
            logic raises, the result has ``met=False`` and
            ``requires_clinical_input=True`` rather than propagating the
            exception.
        """
        try:
            result = self._evaluate(context)
        except Exception as exc:
            logger.error(
                "criterion evaluator raised unexpectedly; deferring to clinical review",
                extra={
                    "code": self.code.value,
                    "variant_id": context.variant.variant_id,
                    "error": str(exc),
                },
            )
            result = CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "Automated evaluation encountered an unexpected error; requires manual review."
                ),
                data_source="AlleleIQ scoring engine",
                requires_clinical_input=True,
                notes=f"internal error: {exc}",
            )

        logger.debug(
            "criterion evaluated",
            extra={
                "code": result.code.value,
                "met": result.met,
                "strength": result.strength.value,
                "requires_clinical_input": result.requires_clinical_input,
                "variant_id": context.variant.variant_id,
            },
        )
        return result

    @abstractmethod
    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        """Subclass-specific evaluation logic.

        May raise; exceptions are caught and translated by :meth:`evaluate`.
        """
