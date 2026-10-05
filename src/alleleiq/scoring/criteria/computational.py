"""Computational prediction ACMG/AMP criteria: PP3, BP4.

Which scores count depends on the variant's consequence, following
ClinGen SVI guidance (Pejaver et al. 2022 for REVEL; Walker et al. 2023
for SpliceAI):

- **Null variants** (the consequences PVS1 evaluates — stop-gained,
  frameshift, canonical splice, start-lost): neither criterion applies.
  Their predicted impact is already captured by PVS1, and adding PP3
  (often from the same SpliceAI signal at a canonical splice site) would
  count one line of evidence twice; BP4 would contradict PVS1 on the
  strength of a splice predictor that says nothing about truncation.
- **Missense variants**: REVEL is the protein-impact predictor and
  SpliceAI the splicing one. PP3 applies if either is high. BP4 requires a
  benign REVEL *and* no splicing signal — a low SpliceAI score alone only
  says the variant does not disrupt splicing, not that the amino-acid
  change is tolerated. A missing REVEL score is a genuine gap
  (``requires_clinical_input``).
- **Other consequences** (synonymous, intronic, UTR, splice-region, ...):
  only SpliceAI is informative; REVEL does not apply.
- **Unknown consequence** (no VEP annotation): whatever scores exist are
  used, with BP4 still blocked by any splicing signal.

Only when no applicable score is available at all is
``requires_clinical_input`` set unconditionally.

The REVEL score itself is resolved by :class:`~alleleiq.clients.revel.REVELClient`
(local tabix index, preferred, falling back to VEP REST) rather than
read directly off the VEP annotation, since VEP's own REVEL coverage is
incomplete — see :func:`_resolve_revel_score`.
"""

from __future__ import annotations

from enum import Enum

from alleleiq.models.evidence import ACMGCode, CriterionEvaluation, EvidenceStrength

from .base import CriterionContext, CriterionEvaluator
from .lof import LOF_CONSEQUENCES

_DATA_SOURCE = "REVEL (local index / VEP REST) / Ensembl VEP (SpliceAI)"
_MISSENSE_CONSEQUENCE = "missense_variant"


class _ConsequenceClass(Enum):
    NULL = "null"
    MISSENSE = "missense"
    OTHER = "other"
    UNKNOWN = "unknown"


def _consequence_class(context: CriterionContext) -> _ConsequenceClass:
    annotation = context.vep_annotation
    consequence = annotation.most_severe_consequence if annotation else None
    if consequence is None:
        return _ConsequenceClass.UNKNOWN
    if consequence in LOF_CONSEQUENCES:
        return _ConsequenceClass.NULL
    if consequence == _MISSENSE_CONSEQUENCE:
        return _ConsequenceClass.MISSENSE
    return _ConsequenceClass.OTHER


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


def _applicable_scores(
    context: CriterionContext, consequence_class: _ConsequenceClass
) -> tuple[float | None, float | None]:
    """``(revel, spliceai_max)`` restricted to the scores this consequence can use."""
    revel = _resolve_revel_score(context)
    spliceai_max = _max_spliceai_delta(context.vep_annotation)
    if consequence_class is _ConsequenceClass.OTHER:
        revel = None
    return revel, spliceai_max


def _null_variant_result(code: ACMGCode, strength: EvidenceStrength) -> CriterionEvaluation:
    return CriterionEvaluation(
        code=code,
        strength=strength,
        met=False,
        evidence_basis=(
            "Null variant: its predicted impact is assessed under PVS1, so computational "
            f"evidence is not added under {code.value}."
        ),
        data_source=_DATA_SOURCE,
        requires_clinical_input=False,
        notes="ClinGen SVI guidance: do not combine PP3/BP4 with PVS1 for the same variant.",
    )


def _no_scores_result(code: ACMGCode, strength: EvidenceStrength) -> CriterionEvaluation:
    return CriterionEvaluation(
        code=code,
        strength=strength,
        met=False,
        evidence_basis="No applicable REVEL or SpliceAI scores available; manual in silico "
        "review required.",
        data_source=_DATA_SOURCE,
        requires_clinical_input=True,
    )


class PP3Criterion(CriterionEvaluator):
    """PP3 (Pathogenic Supporting): computational evidence supports a deleterious effect.

    ACGS 2024 thresholds: REVEL >= 0.70 -> Supporting; REVEL >= 0.90 ->
    Moderate. SpliceAI max delta >= 0.20 -> Supporting; >= 0.50 -> Moderate.
    Not applied to null variants (see module docstring).
    """

    code = ACMGCode.PP3
    default_strength = EvidenceStrength.SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        consequence_class = _consequence_class(context)
        if consequence_class is _ConsequenceClass.NULL:
            return _null_variant_result(self.code, self.default_strength)

        revel, spliceai_max = _applicable_scores(context, consequence_class)
        if revel is None and spliceai_max is None:
            return _no_scores_result(self.code, self.default_strength)

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
            evidence_basis=_evidence_basis(self.code, revel, spliceai_max, met),
            data_source=_DATA_SOURCE,
            requires_clinical_input=consequence_class is _ConsequenceClass.MISSENSE
            and revel is None,
        )


class BP4Criterion(CriterionEvaluator):
    """BP4 (Benign Supporting): computational evidence suggests no deleterious effect.

    ACGS 2024 thresholds: REVEL <= 0.15 -> Supporting. SpliceAI max delta
    < 0.10 -> Supporting. For missense variants both must agree (a benign
    REVEL and no splicing signal); for other non-null consequences only
    SpliceAI applies. Not applied to null variants. No Moderate/Strong
    benign tier is defined for this criterion.
    """

    code = ACMGCode.BP4
    default_strength = EvidenceStrength.BENIGN_SUPPORTING

    def _evaluate(self, context: CriterionContext) -> CriterionEvaluation:
        consequence_class = _consequence_class(context)
        if consequence_class is _ConsequenceClass.NULL:
            return _null_variant_result(self.code, self.default_strength)

        revel, spliceai_max = _applicable_scores(context, consequence_class)
        if revel is None and spliceai_max is None:
            return _no_scores_result(self.code, self.default_strength)

        settings = context.settings
        revel_benign = revel is not None and revel <= settings.revel_bp4_threshold
        splice_benign = spliceai_max is not None and spliceai_max < settings.spliceai_bp4_threshold
        splice_clear = spliceai_max is None or splice_benign

        if consequence_class is _ConsequenceClass.MISSENSE:
            # A low SpliceAI score only rules out a splicing effect; the amino-acid
            # change itself must also be predicted benign.
            met = revel_benign and splice_clear
        elif consequence_class is _ConsequenceClass.OTHER:
            met = splice_benign
        else:
            met = (revel_benign if revel is not None else splice_benign) and splice_clear

        return CriterionEvaluation(
            code=self.code,
            strength=self.default_strength,
            met=met,
            evidence_basis=_evidence_basis(self.code, revel, spliceai_max, met),
            data_source=_DATA_SOURCE,
            requires_clinical_input=consequence_class is _ConsequenceClass.MISSENSE
            and revel is None,
        )


def _evidence_basis(
    code: ACMGCode, revel: float | None, spliceai_max: float | None, met: bool
) -> str:
    parts = []
    if revel is not None:
        parts.append(f"REVEL={revel:.3f}")
    if spliceai_max is not None:
        parts.append(f"max SpliceAI delta={spliceai_max:.3f}")
    scores = ", ".join(parts)
    verdict = "meets" if met else "does not meet"
    return f"{scores}; {verdict} the {code.value} threshold."
