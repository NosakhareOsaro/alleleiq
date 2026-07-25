"""End-to-end variant evidence pipeline: gather evidence, score, build a report.

Shared by the CLI (:mod:`alleleiq.cli`) and the REST API
(:mod:`alleleiq.api.routers.variants`) so the two entry points can never
drift apart on how a variant is actually analysed. Every upstream query
is recorded in the returned report's audit trail, including failures —
a failed gnomAD or VEP query degrades the relevant criteria to
``requires_clinical_input`` rather than aborting the whole analysis or
silently treating missing evidence as absence of evidence.
"""

from __future__ import annotations

from datetime import UTC, datetime

from alleleiq.__version__ import VERSION
from alleleiq.clients.clinvar import ClinVarClient, ClinVarSubmission
from alleleiq.clients.gnomad import GeneConstraint, GnomADClient
from alleleiq.clients.vep import VEPAnnotation, VEPClient
from alleleiq.config import Settings, get_settings
from alleleiq.exceptions import AlleleIQError
from alleleiq.logging import get_logger
from alleleiq.models.evidence import CriterionEvaluation
from alleleiq.models.population import AncestryGroup, PopulationFrequency
from alleleiq.models.report import AuditTrail, DatabaseQueryRecord, EvidenceReport, ReportMetadata
from alleleiq.models.variant import Variant
from alleleiq.scoring.bayesian import accumulate_evidence
from alleleiq.scoring.coverage import assess_coverage
from alleleiq.scoring.criteria.base import CriterionContext, CriterionEvaluator
from alleleiq.scoring.criteria.computational import BP4Criterion, PP3Criterion
from alleleiq.scoring.criteria.lof import PVS1Criterion
from alleleiq.scoring.criteria.population_frequency import BA1Criterion, BS1Criterion, PM2Criterion

logger = get_logger(__name__)

_EVALUATORS: tuple[CriterionEvaluator, ...] = (
    BA1Criterion(),
    BS1Criterion(),
    PM2Criterion(),
    PP3Criterion(),
    BP4Criterion(),
    PVS1Criterion(),
)


async def _fetch_gnomad_and_vep(
    variant: Variant, settings: Settings, queries: list[DatabaseQueryRecord]
) -> tuple[dict[AncestryGroup, PopulationFrequency], VEPAnnotation | None, GeneConstraint | None]:
    async with GnomADClient(settings) as gnomad_client, VEPClient(settings) as vep_client:
        population_frequencies: dict[AncestryGroup, PopulationFrequency] = {}
        started_at = datetime.now(UTC)
        try:
            population_frequencies = await gnomad_client.fetch_population_frequencies(
                variant.coordinates
            )
            queries.append(
                DatabaseQueryRecord(
                    database="gnomAD",
                    query_description="fetch_population_frequencies",
                    parameters={"variant_id": variant.variant_id},
                    timestamp=started_at,
                    success=True,
                )
            )
        except AlleleIQError as exc:
            queries.append(
                DatabaseQueryRecord(
                    database="gnomAD",
                    query_description="fetch_population_frequencies",
                    parameters={"variant_id": variant.variant_id},
                    timestamp=started_at,
                    success=False,
                    error_message=str(exc),
                )
            )
            logger.warning(
                "gnomAD query failed", extra={"variant_id": variant.variant_id, "error": str(exc)}
            )

        vep_annotation: VEPAnnotation | None = None
        vep_started_at = datetime.now(UTC)
        try:
            vep_annotation = await vep_client.fetch_annotation(variant.coordinates)
            queries.append(
                DatabaseQueryRecord(
                    database="VEP",
                    query_description="fetch_annotation",
                    parameters={"variant_id": variant.variant_id},
                    timestamp=vep_started_at,
                    success=True,
                )
            )
        except AlleleIQError as exc:
            queries.append(
                DatabaseQueryRecord(
                    database="VEP",
                    query_description="fetch_annotation",
                    parameters={"variant_id": variant.variant_id},
                    timestamp=vep_started_at,
                    success=False,
                    error_message=str(exc),
                )
            )
            logger.warning(
                "VEP query failed", extra={"variant_id": variant.variant_id, "error": str(exc)}
            )

        gene_constraint: GeneConstraint | None = None
        gene_symbol = (
            vep_annotation.gene_symbol if vep_annotation else None
        ) or variant.gene_symbol
        if gene_symbol:
            constraint_started_at = datetime.now(UTC)
            try:
                gene_constraint = await gnomad_client.fetch_gene_constraint(gene_symbol)
                queries.append(
                    DatabaseQueryRecord(
                        database="gnomAD",
                        query_description="fetch_gene_constraint",
                        parameters={"gene_symbol": gene_symbol},
                        timestamp=constraint_started_at,
                        success=True,
                    )
                )
            except AlleleIQError as exc:
                queries.append(
                    DatabaseQueryRecord(
                        database="gnomAD",
                        query_description="fetch_gene_constraint",
                        parameters={"gene_symbol": gene_symbol},
                        timestamp=constraint_started_at,
                        success=False,
                        error_message=str(exc),
                    )
                )
                logger.warning(
                    "gnomAD gene constraint query failed",
                    extra={"gene_symbol": gene_symbol, "error": str(exc)},
                )

    return population_frequencies, vep_annotation, gene_constraint


async def _fetch_clinvar(
    variant: Variant,
    vep_annotation: VEPAnnotation | None,
    settings: Settings,
    queries: list[DatabaseQueryRecord],
) -> ClinVarSubmission | None:
    rsid = vep_annotation.existing_variant_id if vep_annotation else None
    gene_symbol = (vep_annotation.gene_symbol if vep_annotation else None) or variant.gene_symbol
    hgvs_c = vep_annotation.hgvs_c if vep_annotation else None
    if not rsid and not (gene_symbol and hgvs_c):
        return None

    started_at = datetime.now(UTC)
    try:
        async with ClinVarClient(settings) as clinvar_client:
            submission = await clinvar_client.fetch_classification(
                rsid=rsid,
                gene_symbol=None if rsid else gene_symbol,
                hgvs_c=None if rsid else hgvs_c,
            )
        queries.append(
            DatabaseQueryRecord(
                database="ClinVar",
                query_description="fetch_classification",
                parameters={"variant_id": variant.variant_id},
                timestamp=started_at,
                success=True,
            )
        )
        return submission
    except AlleleIQError as exc:
        queries.append(
            DatabaseQueryRecord(
                database="ClinVar",
                query_description="fetch_classification",
                parameters={"variant_id": variant.variant_id},
                timestamp=started_at,
                success=False,
                error_message=str(exc),
            )
        )
        logger.warning(
            "ClinVar query failed", extra={"variant_id": variant.variant_id, "error": str(exc)}
        )
        return None


def _merge_vep_into_variant(variant: Variant, vep_annotation: VEPAnnotation | None) -> Variant:
    """Fill in gene/transcript/HGVS/consequence from VEP if the parsed variant lacked them."""
    if vep_annotation is None:
        return variant
    return Variant(
        coordinates=variant.coordinates,
        variant_type=variant.variant_type,
        gene_symbol=vep_annotation.gene_symbol or variant.gene_symbol,
        transcript_id=vep_annotation.transcript_id or variant.transcript_id,
        hgvs_c=vep_annotation.hgvs_c or variant.hgvs_c,
        hgvs_p=vep_annotation.hgvs_p or variant.hgvs_p,
        consequence=vep_annotation.most_severe_consequence or variant.consequence,
    )


async def build_evidence_report(
    variant: Variant,
    *,
    settings: Settings | None = None,
    include_clinvar: bool = True,
) -> EvidenceReport:
    """Run one variant through the full evidence-gathering and scoring pipeline.

    Args:
        variant: The parsed variant to analyse.
        settings: Configuration to use; defaults to
            :func:`alleleiq.config.get_settings`.
        include_clinvar: Whether to query ClinVar. ClinVar evidence is
            surfaced in the report's audit trail but is not consumed by
            any scoring criterion (see :mod:`alleleiq.scoring.criteria`),
            so disabling it only affects report completeness, not the
            preliminary classification.

    Returns:
        A complete :class:`~alleleiq.models.report.EvidenceReport`.
    """
    settings = settings or get_settings()
    queries: list[DatabaseQueryRecord] = []

    population_frequencies, vep_annotation, gene_constraint = await _fetch_gnomad_and_vep(
        variant, settings, queries
    )

    clinvar_submission = None
    if include_clinvar:
        clinvar_submission = await _fetch_clinvar(variant, vep_annotation, settings, queries)

    variant = _merge_vep_into_variant(variant, vep_annotation)
    coverage_gaps, coverage_summary = assess_coverage(population_frequencies, settings)

    context = CriterionContext(
        variant=variant,
        population_frequencies=population_frequencies,
        coverage_gaps=coverage_gaps,
        vep_annotation=vep_annotation,
        gene_constraint=gene_constraint,
        clinvar_submission=clinvar_submission,
        settings=settings,
    )

    evaluations: list[CriterionEvaluation] = [
        evaluator.evaluate(context) for evaluator in _EVALUATORS
    ]
    bayesian_score = accumulate_evidence(
        evaluations, overall_coverage_assessment=coverage_summary.overall_coverage_assessment
    )

    metadata = ReportMetadata(
        alleleiq_version=VERSION,
        database_versions={"alleleiq": VERSION, "gnomad": settings.gnomad_dataset},
        query_timestamps={query.database: query.timestamp for query in queries},
    )

    return EvidenceReport(
        variant=variant,
        metadata=metadata,
        population_frequencies=population_frequencies,
        coverage_gaps=coverage_gaps,
        coverage_summary=coverage_summary,
        criterion_evaluations=evaluations,
        bayesian_score=bayesian_score,
        audit_trail=AuditTrail(queries=queries),
    )
