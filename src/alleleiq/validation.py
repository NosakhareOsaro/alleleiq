"""Validation against a curated set of ClinVar expert-panel-reviewed gold-standard variants.

Shared by ``scripts/validate_known_variants.py`` and the ``alleleiq
validate`` CLI command, so the two entry points run identical logic.

Coordinates and expected classifications were resolved live against
ClinVar's VCV XML (``efetch`` with ``rettype=vcv&is_variationid=true``),
reading the GRCh38 ``forDisplay`` ``SequenceLocation``'s
``positionVCF``/``referenceAlleleVCF``/``alternateAlleleVCF`` attributes —
not guessed from memory, which is worth stating plainly given how easy
that is to get subtly wrong (see the VEP/ClinVar client fix that preceded
this module). Re-verify periodically: ClinVar classifications are
occasionally revised.

Evaluation runs through :func:`alleleiq.pipeline.build_evidence_report`
with ``include_clinvar=False``: no criterion evaluator in
:mod:`alleleiq.scoring.criteria` consumes a ClinVar submission, so this
merely avoids a redundant query — the comparison against ClinVar's
expert-panel consensus is not circular either way.

On outcome categories — read this before interpreting a run's output.
AlleleIQ currently automates 6 of the ~28 ACMG/AMP criteria (PVS1, BA1,
BS1, PM2, PP3, BP4). It has no PS3 (functional studies), PM1 (mutational
hotspot), PP1 (segregation), PS2 (de novo), or PS4 (case-control) data —
criteria that frequently supply the deciding points for a real expert
panel's Pathogenic/Benign call. A run against expert-panel-reviewed
variants should therefore expect three outcome categories, not two:

- CONCORDANT: AlleleIQ's direction matches ClinVar's.
- CONSERVATIVE: AlleleIQ lands on VUS where ClinVar has a definitive
  call. This is the *expected, safe* outcome of a criteria subset that
  cannot see the evidence that justified ClinVar's call — it is
  AlleleIQ correctly declining to overreach, not a detection failure.
- DANGEROUS: AlleleIQ's direction is the *opposite* of ClinVar's
  (Pathogenic called Benign, or vice versa). This is the only category
  that indicates an actual correctness bug.
"""

from __future__ import annotations

import asyncio
from dataclasses import dataclass

from rich.console import Console
from rich.table import Table

from alleleiq.models.evidence import PreliminaryClassification
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.pipeline import build_evidence_report

console = Console()


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
    and non-failing, while DANGEROUS is the only category callers should
    treat as a regression signal.
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


async def evaluate_variant(gold: GoldStandardVariant) -> ValidationResult:
    """Run one gold-standard variant through the full evidence-and-scoring pipeline."""
    coordinates = GenomicCoordinates(chrom=gold.chrom, pos=gold.pos, ref=gold.ref, alt=gold.alt)
    variant = Variant(
        coordinates=coordinates,
        variant_type=(
            VariantType.SNV if len(gold.ref) == 1 and len(gold.alt) == 1 else VariantType.INDEL
        ),
        gene_symbol=gold.gene_symbol,
    )
    report = await build_evidence_report(variant, include_clinvar=False)
    score = report.bayesian_score

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


async def run_validation(*, pause_seconds: float = 1.5) -> list[ValidationResult]:
    """Evaluate the full gold-standard set, pausing briefly between variants.

    Args:
        pause_seconds: Courtesy delay between variants for the shared
            public gnomAD/VEP endpoints, on top of each client's own
            per-request throttling.
    """
    results = []
    for index, gold in enumerate(GOLD_STANDARD_SET):
        if index > 0:
            await asyncio.sleep(pause_seconds)
        console.print(f"[dim]Evaluating {gold.label}...[/dim]")
        results.append(await evaluate_variant(gold))
    return results


_OUTCOME_STYLE = {
    Outcome.CONCORDANT: "[green]concordant[/green]",
    Outcome.CONSERVATIVE: "[yellow]conservative (VUS)[/yellow]",
    Outcome.DANGEROUS: "[bold red]DANGEROUS[/bold red]",
}


def render_results(results: list[ValidationResult]) -> None:
    """Print a Rich summary table and concordance breakdown for a validation run."""
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
