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

from collections.abc import Callable
from datetime import UTC, datetime

from alleleiq.__version__ import VERSION
from alleleiq.clients.clingen import ClinGenClient
from alleleiq.clients.clinvar import ClinVarClient, ClinVarPositionalMatch, ClinVarSubmission
from alleleiq.clients.gnomad import GeneConstraint, GnomADClient
from alleleiq.clients.revel import REVELClient
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
from alleleiq.scoring.criteria.gene_level import BP1Criterion, PP2Criterion
from alleleiq.scoring.criteria.lof import PVS1Criterion
from alleleiq.scoring.criteria.population_frequency import BA1Criterion, BS1Criterion, PM2Criterion
from alleleiq.scoring.criteria.same_position import PM5Criterion, PS1Criterion

logger = get_logger(__name__)

ProgressCallback = Callable[[str], None]
"""A synchronous, best-effort progress hook: called with a short stage
description (e.g. ``"querying gnomAD"``) at each evidence-gathering
step. Used by the web interface (see :mod:`alleleiq.web.router`) to
drive a Server-Sent Events progress feed; the CLI and JSON REST API
simply never supply one, which is a complete no-op (see
:func:`_notify`), so this has zero effect on either of those paths.
"""


def _notify(on_progress: ProgressCallback | None, stage: str) -> None:
    """Call ``on_progress(stage)`` if a callback was supplied; otherwise do nothing."""
    if on_progress is not None:
        on_progress(stage)


_EVALUATORS: tuple[CriterionEvaluator, ...] = (
    BA1Criterion(),
    BS1Criterion(),
    PM2Criterion(),
    PP3Criterion(),
    BP4Criterion(),
    PVS1Criterion(),
    PS1Criterion(),
    PM5Criterion(),
    PP2Criterion(),
    BP1Criterion(),
)


async def _fetch_gnomad_and_vep(
    variant: Variant,
    settings: Settings,
    queries: list[DatabaseQueryRecord],
    on_progress: ProgressCallback | None = None,
) -> tuple[dict[AncestryGroup, PopulationFrequency], VEPAnnotation | None, GeneConstraint | None]:
    async with GnomADClient(settings) as gnomad_client, VEPClient(settings) as vep_client:
        population_frequencies: dict[AncestryGroup, PopulationFrequency] = {}
        _notify(on_progress, "querying gnomAD")
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
        _notify(on_progress, "querying VEP")
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
    on_progress: ProgressCallback | None = None,
) -> ClinVarSubmission | None:
    rsid = vep_annotation.existing_variant_id if vep_annotation else None
    gene_symbol = (vep_annotation.gene_symbol if vep_annotation else None) or variant.gene_symbol
    hgvs_c = vep_annotation.hgvs_c if vep_annotation else None
    if not rsid and not (gene_symbol and hgvs_c):
        return None

    _notify(on_progress, "querying ClinVar")
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


async def _fetch_revel(
    variant: Variant,
    vep_annotation: VEPAnnotation | None,
    settings: Settings,
    queries: list[DatabaseQueryRecord],
    on_progress: ProgressCallback | None = None,
) -> tuple[float | None, str | None]:
    """Resolve a REVEL score via REVELClient, recording the outcome in the audit trail.

    ``vep_annotation`` is passed through so REVELClient's VEP fallback
    reuses the annotation already fetched by :func:`_fetch_gnomad_and_vep`
    instead of issuing a second VEP request.
    """
    _notify(on_progress, "querying REVEL")
    started_at = datetime.now(UTC)
    client = REVELClient(settings)
    try:
        result = await client.get_score(variant.coordinates, vep_annotation=vep_annotation)
        queries.append(
            DatabaseQueryRecord(
                database="REVEL",
                query_description=f"get_score (source={result.source})",
                parameters={"variant_id": variant.variant_id},
                timestamp=started_at,
                success=True,
            )
        )
        return result.score, result.source
    except AlleleIQError as exc:
        queries.append(
            DatabaseQueryRecord(
                database="REVEL",
                query_description="get_score",
                parameters={"variant_id": variant.variant_id},
                timestamp=started_at,
                success=False,
                error_message=str(exc),
            )
        )
        logger.warning(
            "REVEL lookup failed", extra={"variant_id": variant.variant_id, "error": str(exc)}
        )
        return None, None
    finally:
        client.close()


async def _fetch_same_position_matches(
    variant: Variant,
    vep_annotation: VEPAnnotation | None,
    settings: Settings,
    queries: list[DatabaseQueryRecord],
) -> list[ClinVarPositionalMatch]:
    """Resolve same-position ClinVar matches for PS1/PM5.

    Skipped (no network call) when VEP has no protein-level change for
    this variant, since PS1/PM5 cannot use the result either way in
    that case — see :mod:`alleleiq.scoring.criteria.same_position`.
    """
    if vep_annotation is None or not vep_annotation.hgvs_p:
        return []

    window = settings.same_position_window_bp
    start = max(1, variant.coordinates.pos - window)
    end = variant.coordinates.pos + window

    started_at = datetime.now(UTC)
    try:
        async with ClinVarClient(settings) as clinvar_client:
            matches = await clinvar_client.search_same_position(
                chrom=variant.coordinates.chrom, start=start, end=end
            )
        queries.append(
            DatabaseQueryRecord(
                database="ClinVar",
                query_description="search_same_position",
                parameters={"variant_id": variant.variant_id, "window_bp": str(window)},
                timestamp=started_at,
                success=True,
            )
        )
        return matches
    except AlleleIQError as exc:
        queries.append(
            DatabaseQueryRecord(
                database="ClinVar",
                query_description="search_same_position",
                parameters={"variant_id": variant.variant_id, "window_bp": str(window)},
                timestamp=started_at,
                success=False,
                error_message=str(exc),
            )
        )
        logger.warning(
            "ClinVar same-position search failed",
            extra={"variant_id": variant.variant_id, "error": str(exc)},
        )
        return []


def _resolve_clingen(
    gene_symbol: str | None,
    settings: Settings,
    queries: list[DatabaseQueryRecord],
    on_progress: ProgressCallback | None = None,
) -> tuple[str | None, int | None]:
    """Resolve ClinGen gene-disease validity classification and haploinsufficiency score.

    A local, synchronous file lookup (see :class:`~alleleiq.clients.clingen.ClinGenClient`)
    rather than a network query, and one that degrades gracefully on its
    own (missing bundled files log a warning and return None/0, they do
    not raise) — so, unlike the other ``_fetch_*`` helpers, there is no
    failure path to catch here. The lookup is still recorded in the
    audit trail like every other evidence source consulted, since it is
    a real bundled reference dataset a reviewer may want to see was used.
    """
    if not gene_symbol:
        return None, None

    _notify(on_progress, "querying ClinGen")
    started_at = datetime.now(UTC)
    client = ClinGenClient(settings)
    validity = client.get_disease_validity(gene_symbol)
    hi_score = client.get_haploinsufficiency_score(gene_symbol)
    queries.append(
        DatabaseQueryRecord(
            database="ClinGen",
            query_description="get_disease_validity + get_haploinsufficiency_score",
            parameters={"gene_symbol": gene_symbol},
            timestamp=started_at,
            success=True,
        )
    )
    classification = validity.classification if validity is not None else None
    return classification, hi_score


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
    on_progress: ProgressCallback | None = None,
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
        on_progress: Optional callback invoked with a short stage
            description (``"querying gnomAD"``, ``"scoring"``, ...) as
            each evidence-gathering step starts. A no-op when omitted;
            the CLI and JSON REST API never pass one. See
            :mod:`alleleiq.web.router` for the one caller that does.

    Returns:
        A complete :class:`~alleleiq.models.report.EvidenceReport`.
    """
    settings = settings or get_settings()
    queries: list[DatabaseQueryRecord] = []

    population_frequencies, vep_annotation, gene_constraint = await _fetch_gnomad_and_vep(
        variant, settings, queries, on_progress
    )

    clinvar_submission = None
    same_position_matches: list[ClinVarPositionalMatch] = []
    if include_clinvar:
        clinvar_submission = await _fetch_clinvar(
            variant, vep_annotation, settings, queries, on_progress
        )
        same_position_matches = await _fetch_same_position_matches(
            variant, vep_annotation, settings, queries
        )

    revel_score, revel_source = await _fetch_revel(
        variant, vep_annotation, settings, queries, on_progress
    )

    gene_symbol_for_clingen = (
        vep_annotation.gene_symbol if vep_annotation else None
    ) or variant.gene_symbol
    clingen_classification, clingen_hi_score = _resolve_clingen(
        gene_symbol_for_clingen, settings, queries, on_progress
    )

    variant = _merge_vep_into_variant(variant, vep_annotation)
    coverage_gaps, coverage_summary = assess_coverage(population_frequencies, settings)

    _notify(on_progress, "scoring")
    context = CriterionContext(
        variant=variant,
        population_frequencies=population_frequencies,
        coverage_gaps=coverage_gaps,
        vep_annotation=vep_annotation,
        gene_constraint=gene_constraint,
        clinvar_submission=clinvar_submission,
        same_position_clinvar_matches=same_position_matches,
        revel_score=revel_score,
        revel_source=revel_source,
        clingen_classification=clingen_classification,
        clingen_haploinsufficiency_score=clingen_hi_score,
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
