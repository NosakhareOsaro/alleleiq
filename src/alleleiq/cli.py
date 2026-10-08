"""AlleleIQ command-line interface."""

from __future__ import annotations

import asyncio
from pathlib import Path
from typing import cast

import click
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

from alleleiq.__version__ import VERSION
from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.clients.vep import VEPClient
from alleleiq.config import get_settings
from alleleiq.exceptions import AlleleIQError, VariantParseError
from alleleiq.logging import configure_logging, get_logger
from alleleiq.models.evidence import PreliminaryClassification
from alleleiq.models.report import EvidenceReport
from alleleiq.models.variant import GenomeBuild
from alleleiq.parsers.vcf_parser import parse_vcf
from alleleiq.pipeline import build_evidence_report
from alleleiq.reports import ReportFormat, ReportOrchestrator
from alleleiq.validation import Outcome, render_results, run_validation

logger = get_logger(__name__)
console = Console(stderr=False)
error_console = Console(stderr=True)

_CLASSIFICATION_STYLE: dict[PreliminaryClassification, str] = {
    PreliminaryClassification.PATHOGENIC: "bold red",
    PreliminaryClassification.LIKELY_PATHOGENIC: "red",
    PreliminaryClassification.UNCERTAIN_SIGNIFICANCE: "yellow",
    PreliminaryClassification.LIKELY_BENIGN: "green",
    PreliminaryClassification.BENIGN: "bold green",
}


@click.group()
@click.version_option(version=VERSION, prog_name="alleleiq")
def cli() -> None:
    """AlleleIQ: population-inclusive evidence aggregation for clinical genomic
    variant interpretation.

    AlleleIQ aggregates evidence — it does not classify variants. Every
    report requires review by a qualified Clinical Genetic Scientist.
    """


@cli.command()
@click.argument("vcf_file", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=Path("./reports"),
    show_default=True,
    help="Directory to write reports.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["json", "html", "both"]),
    default="both",
    show_default=True,
    help="Output format.",
)
@click.option(
    "--genome-build",
    type=click.Choice(["GRCh38"]),
    default="GRCh38",
    show_default=True,
    help="Genome build (GRCh38 only in this release).",
)
@click.option("--verbose", is_flag=True, help="Enable DEBUG logging.")
@click.option(
    "--no-clinvar",
    is_flag=True,
    help="Skip ClinVar queries (faster; PS1, PM5, and PM1's density signal go unassessed).",
)
def analyse(
    vcf_file: Path,
    output_dir: Path,
    output_format: str,
    genome_build: str,
    verbose: bool,
    no_clinvar: bool,
) -> None:
    """Analyse variants in VCF_FILE and generate evidence reports."""
    configure_logging("DEBUG" if verbose else "INFO")

    try:
        variants, stats = parse_vcf(vcf_file, expected_build=GenomeBuild(genome_build))
    except VariantParseError as exc:
        error_console.print(f"[bold red]Error:[/bold red] {exc}")
        error_console.print(
            "[dim]Suggestion: verify the VCF is well-formed and declares (or fingerprints as) "
            "GRCh38. AlleleIQ does not support GRCh37 input.[/dim]"
        )
        raise SystemExit(1) from exc

    console.print(
        f"Parsed [bold]{stats.variants_emitted}[/bold] variant(s) from {vcf_file} "
        f"({stats.variants_filtered} filtered, {stats.variants_normalised} normalised)."
    )
    if stats.duplicates_removed:
        console.print(
            f"[yellow]Note:[/yellow] {stats.duplicates_removed} duplicate variant record(s) "
            "were skipped; each variant is analysed once."
        )
    if not stats.genome_build_verified:
        console.print(
            "[yellow]Warning:[/yellow] genome build could not be verified from the VCF header; "
            "proceeding assuming GRCh38."
        )

    formats: tuple[ReportFormat, ...] = (
        ("json", "html")
        if output_format == "both"
        # Safe: click.Choice(["json", "html", "both"]) guarantees this value at runtime.
        else (cast(ReportFormat, output_format),)
    )
    orchestrator = ReportOrchestrator(output_dir)
    settings = get_settings()

    reports: list[EvidenceReport] = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Analysing variants", total=len(variants))
        for variant in variants:
            try:
                report = asyncio.run(
                    build_evidence_report(
                        variant, settings=settings, include_clinvar=not no_clinvar
                    )
                )
                orchestrator.generate(report, formats=formats)
                reports.append(report)
            except AlleleIQError as exc:
                error_console.print(
                    f"[bold red]Error analysing {variant.variant_id}:[/bold red] {exc}"
                )
                logger.error(
                    "variant analysis failed",
                    extra={"variant_id": variant.variant_id, "error": str(exc)},
                )
            progress.advance(task)

    _render_analysis_summary(reports, output_dir)
    logger.info("gnomAD client cache statistics", extra=GnomADClient.cache_stats())


def _render_analysis_summary(reports: list[EvidenceReport], output_dir: Path) -> None:
    table = Table(title="AlleleIQ analysis summary")
    table.add_column("Variant")
    table.add_column("Gene")
    table.add_column("Classification")
    table.add_column("Points", justify="right")
    table.add_column("Coverage gaps", justify="right")

    total_coverage_gaps = 0
    for report in reports:
        style = _CLASSIFICATION_STYLE.get(report.bayesian_score.preliminary_classification, "")
        total_coverage_gaps += len(report.coverage_gaps)
        table.add_row(
            report.variant.variant_id,
            report.variant.gene_symbol or "—",
            f"[{style}]{report.bayesian_score.preliminary_classification.value}[/{style}]",
            str(report.bayesian_score.evidence_points),
            str(len(report.coverage_gaps)),
        )
    console.print(table)
    console.print(
        f"\n{len(reports)} report(s) written to {output_dir}. "
        f"{total_coverage_gaps} total population coverage gap(s) flagged across all variants."
    )
    console.print(
        "\n[bold]PRELIMINARY EVIDENCE SUMMARY — FOR CLINICAL SCIENTIST REVIEW ONLY.[/bold] "
        "These reports do not constitute variant classifications."
    )


@cli.command()
def validate() -> None:
    """Run AlleleIQ against the bundled ClinVar gold-standard set for QC.

    See alleleiq.validation for the full explanation of the outcome
    categories (concordant / conservative / dangerous) this prints.
    """
    configure_logging("WARNING")
    results = asyncio.run(run_validation())
    render_results(results)
    if any(r.outcome == Outcome.DANGEROUS for r in results):
        raise SystemExit(1)


@cli.command()
def version() -> None:
    """Print version, database versions, and API connectivity status."""
    settings = get_settings()
    console.print(f"AlleleIQ [bold]{VERSION}[/bold]")
    console.print(f"gnomAD dataset: {settings.gnomad_dataset}")
    console.print(f"gnomAD API: {settings.gnomad_api_url}")
    console.print(f"VEP API: {settings.vep_api_url}")
    console.print(f"ClinVar API: {settings.clinvar_api_url}")

    async def _check_connectivity() -> dict[str, bool]:
        async with (
            GnomADClient(settings) as gnomad_client,
            VEPClient(settings) as vep_client,
            ClinVarClient(settings) as clinvar_client,
        ):
            return {
                "gnomAD": await gnomad_client.health_check(),
                "VEP": await vep_client.health_check(),
                "ClinVar": await clinvar_client.health_check(),
            }

    console.print("\nConnectivity:")
    statuses = asyncio.run(_check_connectivity())
    for name, healthy in statuses.items():
        style = "green" if healthy else "red"
        label = "reachable" if healthy else "unreachable"
        console.print(f"  {name}: [{style}]{label}[/{style}]")


if __name__ == "__main__":
    cli()
