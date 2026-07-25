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
from typing import ClassVar

from alleleiq.clients.clinvar import ClinVarSubmission
from alleleiq.clients.gnomad import GeneConstraint
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
