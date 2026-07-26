"""Computational prediction ACMG/AMP criteria: PP3, BP4.

REVEL is a missense-specific predictor: its absence for a missense
variant is a genuine evidentiary gap (flagged
``requires_clinical_input``), but its absence for, say, a splice-donor
variant simply means REVEL doesn't apply there — that is not a gap and
is not flagged. SpliceAI, by contrast, is broadly applicable regardless
of the VEP-reported most-severe consequence (a splicing effect is not
limited to variants VEP labels as splice-region), so it is always
considered when present. Only when *neither* score is available at all
does the criterion have no computational signal to evaluate, which is
when ``requires_clinical_input`` is set unconditionally.

The REVEL score itself is resolved by :class:`~alleleiq.clients.revel.REVELClient`
(local tabix index, preferred, falling back to VEP REST) rather than
read directly off the VEP annotation, since VEP's own REVEL coverage is
incomplete — see :func:`_resolve_revel_score`.
"""

from __future__ import annotations

from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength

from .base import CriterionContext, CriterionEvaluator

_DATA_SOURCE = "REVEL (local index / VEP REST) / Ensembl VEP (SpliceAI)"
_MISSENSE_CONSEQUENCE = "missense_variant"


def _resolve_revel_score(context: CriterionContext) -> float | None:
    """The REVEL score to use, preferring REVELClient's result over the raw VEP annotation.

    ``context.revel_score`` is populated by the pipeline via
    :class:`~alleleiq.clients.revel.REVELClient`, which already
    incorporates the VEP-REST fallback; ``vep_annotation.revel_score`` is
    used directly only when a caller builds a :class:`CriterionContext`
    without going through REVELClient (e.g. existing unit tests).
    """
    if context.revel_score is not None:
        return context.revel_score
    return context.vep_annotation.revel_score if context.vep_annotation else None


def _max_spliceai_delta(vep_annotation: object) -> float | None:
    if vep_annotation is None:
        return None
    scores = [
        getattr(vep_annotation, "spliceai_ds_ag", None),
        getattr(vep_annotation, "spliceai_ds_al", None),
        getattr(vep_annotation, "spliceai_ds_dg", None),
        getattr(vep_annotation, "spliceai_ds_dl", None),
    ]
    available = [s for s in scores if s is not None]
    return max(available) if available else None


class PP3Criterion(CriterionEvaluator):
    """PP3 (Pathogenic Supporting): computational evidence supports a deleterious effect.

    ACGS 2024 thresholds: REVEL >= 0.70 -> Supporting; REVEL >= 0.90 ->
    Moderate. SpliceAI max delta >= 0.20 -> Supporting; >= 0.50 -> Moderate.
    """

    code = ACMGCode.PP3
    default_strength = EvidenceStrength.SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        annotation = context.vep_annotation
        revel = _resolve_revel_score(context)
        spliceai_max = _max_spliceai_delta(annotation)
        consequence = annotation.most_severe_consequence if annotation else None
        is_missense = consequence == _MISSENSE_CONSEQUENCE

        if revel is None and spliceai_max is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No REVEL or SpliceAI scores available; manual in silico review required."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
            )

        settings = context.settings
        revel_supporting = revel is not None and revel >= settings.revel_pp3_supporting_threshold
        revel_moderate = revel is not None and revel >= settings.revel_pp3_moderate_threshold
        splice_supporting = (
            spliceai_max is not None and spliceai_max >= settings.spliceai_pp3_supporting_threshold
        )
        splice_moderate = (
            spliceai_max is not None and spliceai_max >= settings.spliceai_pp3_moderate_threshold
        )

        met = revel_supporting or splice_supporting
        strength = (
            EvidenceStrength.MODERATE
            if (revel_moderate or splice_moderate)
            else EvidenceStrength.SUPPORTING
        )

        return CriterionEvaluation(
            code=self.code,
            strength=strength,
            met=met,
            evidence_basis=self._evidence_basis(revel, spliceai_max, met),
            data_source=_DATA_SOURCE,
            requires_clinical_input=is_missense and revel is None,
        )

    @staticmethod
    def _evidence_basis(revel: float | None, spliceai_max: float | None, met: bool) -> str:
        parts = []
        if revel is not None:
            parts.append(f"REVEL={revel:.3f}")
        if spliceai_max is not None:
            parts.append(f"max SpliceAI delta={spliceai_max:.3f}")
        scores = ", ".join(parts)
        verdict = "meets" if met else "does not meet"
        return f"{scores}; {verdict} the PP3 threshold."


class BP4Criterion(CriterionEvaluator):
    """BP4 (Benign Supporting): computational evidence suggests no deleterious effect.

    ACGS 2024 thresholds: REVEL <= 0.15 -> Supporting. SpliceAI max delta
    < 0.10 -> Supporting. No Moderate/Strong benign tier is defined for
    this criterion.
    """

    code = ACMGCode.BP4
    default_strength = EvidenceStrength.BENIGN_SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        annotation = context.vep_annotation
        revel = _resolve_revel_score(context)
        spliceai_max = _max_spliceai_delta(annotation)
        consequence = annotation.most_severe_consequence if annotation else None
        is_missense = consequence == _MISSENSE_CONSEQUENCE

        if revel is None and spliceai_max is None:
            return CriterionEvaluation(
                code=self.code,
                strength=self.default_strength,
                met=False,
                evidence_basis=(
                    "No REVEL or SpliceAI scores available; manual in silico review required."
                ),
                data_source=_DATA_SOURCE,
                requires_clinical_input=True,
            )

        settings = context.settings
        revel_benign = revel is not None and revel <= settings.revel_bp4_threshold
        splice_benign = spliceai_max is not None and spliceai_max < settings.spliceai_bp4_threshold
        met = revel_benign or splice_benign

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=met,
            evidence_basis=self._evidence_basis(revel, spliceai_max, met),
            data_source=_DATA_SOURCE,
            requires_clinical_input=is_missense and revel is None,
        )

    @staticmethod
    def _evidence_basis(revel: float | None, spliceai_max: float | None, met: bool) -> str:
        parts = []
        if revel is not None:
            parts.append(f"REVEL={revel:.3f}")
        if spliceai_max is not None:
            parts.append(f"max SpliceAI delta={spliceai_max:.3f}")
        scores = ", ".join(parts)
        verdict = "meets" if met else "does not meet"
        return f"{scores}; {verdict} the BP4 threshold."
