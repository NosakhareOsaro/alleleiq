"""Exception hierarchy for AlleleIQ.

Every error condition the system can encounter is represented by a
dedicated exception class rather than a generic ``Exception`` or a bare
string. This lets callers (CLI, API, scoring engine) catch failures at the
right granularity and lets the report generator distinguish "no evidence
found" (a legitimate, reportable state) from "the system could not
determine whether evidence exists" (a failure that must be surfaced, not
silently treated as absence of evidence).
"""

from __future__ import annotations


class AlleleIQError(Exception):
    """Base class for all exceptions raised intentionally by AlleleIQ.

    Catching this exception (rather than ``Exception``) at API and CLI
    boundaries guarantees that only anticipated, classified failure modes
    are handled as "expected" errors; anything else is an unhandled bug
    and should propagate as such.
    """


class ConfigurationError(AlleleIQError):
    """Raised when application configuration is missing or invalid.

    Examples include an unparsable ``.env`` value, a required setting left
    unset, or a configured threshold outside its valid range.
    """


class VariantParseError(AlleleIQError):
    """Raised when a VCF file or variant record cannot be parsed or
    normalised into a valid :class:`alleleiq.models.variant.Variant`.

    Covers malformed VCF syntax, unsupported genome builds, invalid
    nucleotide sequences, and structurally implausible coordinates.
    """


class DatabaseQueryError(AlleleIQError):
    """Base class for failures querying an external evidence source.

    Raised when a query to an upstream database fails in a way that
    prevents AlleleIQ from determining the answer — as distinct from a
    successful query that legitimately returns "no data" (e.g. a variant
    genuinely absent from gnomAD). The distinction matters clinically: a
    failed query must become an explicit "not assessed" state in the
    report, never a silent absence of evidence.
    """


class GnomADError(DatabaseQueryError):
    """Raised when the gnomAD v4 GraphQL client cannot complete a query.

    Covers network failures, GraphQL error responses, malformed payloads,
    and exhausted retry budgets.
    """


class ClinVarError(DatabaseQueryError):
    """Raised when the ClinVar E-utilities client cannot complete a query.

    Covers network failures, NCBI rate-limit exhaustion, and malformed
    XML/JSON responses from ``esearch``/``efetch``.
    """


class VEPError(DatabaseQueryError):
    """Raised when the Ensembl VEP REST client cannot complete a query.

    Covers network failures, VEP server downtime, and malformed responses.
    Consumers should generally degrade gracefully (mark the relevant
    criteria ``requires_clinical_input``) rather than aborting the entire
    analysis, since VEP is one of several evidence sources.
    """


class UniProtError(DatabaseQueryError):
    """Raised when the UniProt REST client cannot complete a query.

    Covers network failures and malformed responses. Consumers should
    degrade PM1 to ``requires_clinical_input`` rather than aborting the
    entire analysis, since UniProt domain annotation is one of several
    PM1 evidence sources (see :mod:`alleleiq.scoring.criteria.hotspot`).
    """


class ScoringError(AlleleIQError):
    """Base class for failures in the ACMG/AMP evidence scoring engine."""


class InsufficientEvidenceError(ScoringError):
    """Raised when a criterion evaluator cannot proceed for lack of the
    minimum evidence it requires to produce even a "not met" result.

    Criterion evaluators should generally prefer returning a
    ``CriterionEvaluation`` with ``requires_clinical_input=True`` over
    raising this exception; it exists for the small number of cases where
    even that degraded response is not constructible (e.g. no variant
    coordinates at all).
    """


class ReportGenerationError(AlleleIQError):
    """Raised when a JSON or HTML evidence report cannot be generated or
    written to disk, including template rendering failures and output
    validation failures (e.g. an unexpectedly empty output file).
    """
