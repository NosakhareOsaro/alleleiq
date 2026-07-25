#!/usr/bin/env python3
"""Validate AlleleIQ's preliminary classifications against a curated set of
10 ClinVar gold-standard variants, each a "reviewed by expert panel"
(3-star) record — ClinVar's highest review tier below an official practice
guideline.

This is the embryo of AlleleIQ's formal clinical validation study, not a
substitute for one: 10 variants drawn from a handful of well-studied
cancer-predisposition genes cannot establish clinical-grade sensitivity or
specificity. Its purpose is narrower and concrete: catch obvious
regressions in the scoring pipeline, and give a reproducible concordance
number to report while that formal study is designed.

Coordinates and expected classifications were resolved live against
ClinVar's VCV XML (``efetch`` with ``rettype=vcv&is_variationid=true``),
reading the GRCh38 ``forDisplay`` ``SequenceLocation``'s
``positionVCF``/``referenceAlleleVCF``/``alternateAlleleVCF`` attributes —
not guessed from memory, which is worth stating plainly given how easy
that is to get subtly wrong (see the fix in the commit that preceded this
script). Re-verify periodically: ClinVar classifications are occasionally
revised.

Deliberately does not query ClinVar as part of scoring: no criterion
evaluator in :mod:`alleleiq.scoring.criteria` consumes a ClinVar
submission, so comparing AlleleIQ's output against ClinVar's expert-panel
consensus is not circular.

On outcome categories — read this before interpreting a run's output.
AlleleIQ currently automates 6 of the ~28 ACMG/AMP criteria (PVS1, BA1,
BS1, PM2, PP3, BP4). It has no PS3 (functional studies), PM1 (mutational
hotspot), PP1 (segregation), PS2 (de novo), or PS4 (case-control) data —
criteria that frequently supply the deciding points for a real expert
panel's Pathogenic/Benign call. A run of this script against
expert-panel-reviewed variants should therefore expect three outcome
categories, not two:

- CONCORDANT: AlleleIQ's direction matches ClinVar's.
- CONSERVATIVE: AlleleIQ lands on VUS where ClinVar has a definitive
  call. This is the *expected, safe* outcome of a criteria subset that
  cannot see the evidence that justified ClinVar's call — it is
  AlleleIQ correctly declining to overreach, not a detection failure.
- DANGEROUS: AlleleIQ's direction is the *opposite* of ClinVar's
  (Pathogenic called Benign, or vice versa). This is the only category
  that indicates an actual correctness bug, and is what this script's
  exit status gates on.

Usage:
    uv run python scripts/validate_known_variants.py

Exit status is non-zero only if any variant is DANGEROUS, so this can be
wired into CI as a regression gate without failing on expected
conservative behaviour.
"""

from __future__ import annotations

import asyncio
import sys
from dataclasses import dataclass

from rich.console import Console
from rich.table import Table

from alleleiq.clients.gnomad import GeneConstraint, GnomADClient
from alleleiq.clients.vep import VEPAnnotation, VEPClient
from alleleiq.exceptions import AlleleIQError
from alleleiq.logging import configure_logging, get_logger
from alleleiq.models.evidence import CriterionEvaluation, PreliminaryClassification
from alleleiq.models.population import AncestryGroup, CoverageGap, PopulationFrequency
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.scoring.bayesian import accumulate_evidence
from alleleiq.scoring.coverage import assess_coverage
from alleleiq.scoring.criteria.base import CriterionContext, CriterionEvaluator
from alleleiq.scoring.criteria.computational import BP4Criterion, PP3Criterion
from alleleiq.scoring.criteria.lof import PVS1Criterion
from alleleiq.scoring.criteria.population_frequency import BA1Criterion, BS1Criterion, PM2Criterion

logger = get_logger(__name__)
console = Console()

_EVALUATORS: tuple[CriterionEvaluator, ...] = (
    BA1Criterion(),
    BS1Criterion(),
    PM2Criterion(),
    PP3Criterion(),
    BP4Criterion(),
    PVS1Criterion(),
)


@dataclass(frozen=True)
class GoldStandardVariant:
    """A single ClinVar expert-panel-reviewed benchmark variant."""

    label: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene_symbol: str
    clinvar_accession: str
    expected_classification: PreliminaryClassification


GOLD_STANDARD_SET: tuple[GoldStandardVariant, ...] = (
    GoldStandardVariant(
        "MLH1 c.543C>T (p.Gly181=)",
        "3",
        37008903,
        "C",
        "T",
        "MLH1",
        "VCV000633499",
        PreliminaryClassification.PATHOGENIC,
    ),
    GoldStandardVariant(
        "MSH2 c.2083G>A (p.Val695Met)",
        "2",
        47476444,
        "G",
        "A",
        "MSH2",
        "VCV000633496",
        PreliminaryClassification.PATHOGENIC,
    ),
    GoldStandardVariant(
        "MSH6 c.1439T>A (p.Val480Glu)",
        "2",
        47799422,
        "T",
        "A",
        "MSH6",
        "VCV000633495",
        PreliminaryClassification.PATHOGENIC,
    ),
    GoldStandardVariant(
        "APC c.2546_2552del (p.Asp849fs)",
        "5",
        112838136,
        "AAAGATAG",
        "A",
        "APC",
        "VCV003901082",
        PreliminaryClassification.PATHOGENIC,
    ),
    GoldStandardVariant(
        "CFTR c.[4C>T;7A>T]",
        "7",
        117480098,
        "C",
        "T",
        "CFTR",
        "VCV001320029",
        PreliminaryClassification.PATHOGENIC,
    ),
    GoldStandardVariant(
        "VHL c.257C>T (p.Pro86Leu)",
        "3",
        10142104,
        "C",
        "T",
        "VHL",
        "VCV000182977",
        PreliminaryClassification.PATHOGENIC,
    ),
    GoldStandardVariant(
        "BRCA1 c.4081A>T (p.Met1361Leu)",
        "17",
        43091450,
        "T",
        "A",
        "BRCA1",
        "VCV000812551",
        PreliminaryClassification.BENIGN,
    ),
    GoldStandardVariant(
        "BRCA2 c.8754+75A>G",
        "13",
        32376866,
        "A",
        "G",
        "BRCA2",
        "VCV000812553",
        PreliminaryClassification.BENIGN,
    ),
    GoldStandardVariant(
        "PMS2 c.804-243G>A",
        "7",
        5995876,
        "C",
        "T",
        "PMS2",
        "VCV000091371",
        PreliminaryClassification.BENIGN,
    ),
    GoldStandardVariant(
        "VHL c.154G>A (p.Glu52Lys)",
        "3",
        10142001,
        "G",
        "A",
        "VHL",
        "VCV000161402",
        PreliminaryClassification.BENIGN,
    ),
)


class Outcome:
    """The three possible outcome categories for one gold-standard comparison.

    See the module docstring's "On outcome categories" section for the
    clinical-safety reasoning behind treating CONSERVATIVE as expected
    and non-failing, while DANGEROUS is the only category this script's
    exit status reacts to.
    """

    CONCORDANT = "CONCORDANT"
    CONSERVATIVE = "CONSERVATIVE"
    DANGEROUS = "DANGEROUS"


@dataclass(frozen=True)
class ValidationResult:
    """The outcome of running one gold-standard variant through the pipeline."""

    gold_standard: GoldStandardVariant
    actual_classification: PreliminaryClassification
    evidence_points: int
    criteria_applied: list[str]
    outcome: str
    explanation: str | None


def _direction(classification: PreliminaryClassification) -> str:
    """Collapse the five-tier classification into PATHOGENIC/BENIGN/UNCERTAIN.

    Concordance is assessed on direction, not exact tier: AlleleIQ calling
    a ClinVar-Pathogenic variant "Likely Pathogenic" is the correct
    real-world outcome for a conservative automated tool and should not
    register as a discordance.
    """
    if classification in (
        PreliminaryClassification.PATHOGENIC,
        PreliminaryClassification.LIKELY_PATHOGENIC,
    ):
        return "PATHOGENIC"
    if classification in (
        PreliminaryClassification.BENIGN,
        PreliminaryClassification.LIKELY_BENIGN,
    ):
        return "BENIGN"
    return "UNCERTAIN"


async def _fetch_evidence(
    coordinates: GenomicCoordinates,
) -> tuple[dict[AncestryGroup, PopulationFrequency], VEPAnnotation | None, GeneConstraint | None]:
    async with GnomADClient() as gnomad_client, VEPClient() as vep_client:
        try:
            population_frequencies = await gnomad_client.fetch_population_frequencies(coordinates)
        except AlleleIQError as exc:
            logger.warning(
                "gnomAD query failed during validation",
                extra={"variant_id": coordinates.variant_id, "error": str(exc)},
            )
            population_frequencies = {}

        try:
            vep_annotation = await vep_client.fetch_annotation(coordinates)
        except AlleleIQError as exc:
            logger.warning(
                "VEP query failed during validation",
                extra={"variant_id": coordinates.variant_id, "error": str(exc)},
            )
            vep_annotation = None

        gene_constraint = None
        if vep_annotation and vep_annotation.gene_symbol:
            try:
                gene_constraint = await gnomad_client.fetch_gene_constraint(
                    vep_annotation.gene_symbol
                )
            except AlleleIQError as exc:
                logger.warning(
                    "gnomAD gene constraint query failed during validation",
                    extra={"gene_symbol": vep_annotation.gene_symbol, "error": str(exc)},
                )

    return population_frequencies, vep_annotation, gene_constraint


async def evaluate_variant(gold: GoldStandardVariant) -> ValidationResult:
    """Run one gold-standard variant through the full evidence-and-scoring pipeline."""
    coordinates = GenomicCoordinates(chrom=gold.chrom, pos=gold.pos, ref=gold.ref, alt=gold.alt)
    population_frequencies, vep_annotation, gene_constraint = await _fetch_evidence(coordinates)

    variant = Variant(
        coordinates=coordinates,
        variant_type=(
            VariantType.SNV if len(gold.ref) == 1 and len(gold.alt) == 1 else VariantType.INDEL
        ),
        gene_symbol=(vep_annotation.gene_symbol if vep_annotation else gold.gene_symbol),
        transcript_id=vep_annotation.transcript_id if vep_annotation else None,
        hgvs_c=vep_annotation.hgvs_c if vep_annotation else None,
        hgvs_p=vep_annotation.hgvs_p if vep_annotation else None,
        consequence=vep_annotation.most_severe_consequence if vep_annotation else None,
    )

    coverage_gaps: dict[AncestryGroup, CoverageGap]
    coverage_gaps, coverage_summary = assess_coverage(population_frequencies)

    context = CriterionContext(
        variant=variant,
        population_frequencies=population_frequencies,
        coverage_gaps=coverage_gaps,
        vep_annotation=vep_annotation,
        gene_constraint=gene_constraint,
    )

    evaluations: list[CriterionEvaluation] = [
        evaluator.evaluate(context) for evaluator in _EVALUATORS
    ]
    score = accumulate_evidence(
        evaluations, overall_coverage_assessment=coverage_summary.overall_coverage_assessment
    )

    expected_direction = _direction(gold.expected_classification)
    actual_direction = _direction(score.preliminary_classification)

    if actual_direction == expected_direction:
        outcome = Outcome.CONCORDANT
    elif actual_direction == "UNCERTAIN":
        outcome = Outcome.CONSERVATIVE
    else:
        outcome = Outcome.DANGEROUS

    explanation = None
    if outcome != Outcome.CONCORDANT:
        applied = ", ".join(c.value for c in score.criteria_applied) or "none"
        explanation = (
            f"expected {gold.expected_classification.value} ({expected_direction}) but AlleleIQ "
            f"produced {score.preliminary_classification.value} ({actual_direction}) at "
            f"{score.evidence_points} evidence points; contributing criteria: {applied}."
        )

    return ValidationResult(
        gold_standard=gold,
        actual_classification=score.preliminary_classification,
        evidence_points=score.evidence_points,
        criteria_applied=[c.value for c in score.criteria_applied],
        outcome=outcome,
        explanation=explanation,
    )


async def run_validation() -> list[ValidationResult]:
    results = []
    for index, gold in enumerate(GOLD_STANDARD_SET):
        if index > 0:
            # A courtesy pause between variants for the shared public gnomAD/VEP
            # endpoints, on top of each client's own per-request throttling.
            await asyncio.sleep(1.5)
        console.print(f"[dim]Evaluating {gold.label}...[/dim]")
        results.append(await evaluate_variant(gold))
    return results


_OUTCOME_STYLE = {
    Outcome.CONCORDANT: "[green]concordant[/green]",
    Outcome.CONSERVATIVE: "[yellow]conservative (VUS)[/yellow]",
    Outcome.DANGEROUS: "[bold red]DANGEROUS[/bold red]",
}


def render_results(results: list[ValidationResult]) -> None:
    table = Table(title="AlleleIQ known-variant validation")
    table.add_column("Variant")
    table.add_column("ClinVar (expert panel)")
    table.add_column("AlleleIQ preliminary")
    table.add_column("Points", justify="right")
    table.add_column("Criteria applied")
    table.add_column("Outcome")

    for result in results:
        table.add_row(
            result.gold_standard.label,
            result.gold_standard.expected_classification.value,
            result.actual_classification.value,
            str(result.evidence_points),
            ", ".join(result.criteria_applied) or "—",
            _OUTCOME_STYLE[result.outcome],
        )
    console.print(table)

    counts = {
        outcome: sum(1 for r in results if r.outcome == outcome)
        for outcome in (Outcome.CONCORDANT, Outcome.CONSERVATIVE, Outcome.DANGEROUS)
    }
    console.print(
        f"\n{counts[Outcome.CONCORDANT]}/{len(results)} concordant, "
        f"{counts[Outcome.CONSERVATIVE]}/{len(results)} conservative (expected, safe), "
        f"{counts[Outcome.DANGEROUS]}/{len(results)} dangerous."
    )

    dangerous = [r for r in results if r.outcome == Outcome.DANGEROUS]
    if dangerous:
        console.print("\n[bold red]Dangerous discordances (investigate immediately):[/bold red]")
        for result in dangerous:
            console.print(
                f"  - {result.gold_standard.label} ({result.gold_standard.clinvar_accession}): "
                f"{result.explanation}"
            )
    else:
        console.print("\n[green]No dangerous (opposite-direction) discordances.[/green]")


def main() -> int:
    configure_logging("WARNING")
    results = asyncio.run(run_validation())
    render_results(results)
    return 0 if all(r.outcome != Outcome.DANGEROUS for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
