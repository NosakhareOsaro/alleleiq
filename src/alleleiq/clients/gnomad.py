"""gnomAD v4 GraphQL client: per-ancestry population frequencies.

This is AlleleIQ's most clinically significant integration. Every value
returned here is kept scoped to a single :class:`AncestryGroup` for the
lifetime of the pipeline — genome and exome allele counts are summed per
population (gnomAD reports them as separate cohorts), but never summed
*across* populations into a single aggregate frequency. That aggregation
is exactly the step that hides ancestry-specific coverage gaps, and
AlleleIQ's entire clinical-equity value proposition depends on not taking
it.

Known limitation (v0.1.0-alpha): mitochondrial variants use a different
gnomAD schema (heteroplasmy-based, not diploid AC/AN) and are not yet
supported; :meth:`GnomADClient.fetch_population_frequencies` raises
:class:`~alleleiq.exceptions.GnomADError` for chrMT coordinates rather
than silently returning incorrect nuclear-style frequencies.
"""

from __future__ import annotations

from typing import Any, ClassVar

import httpx
from pydantic import BaseModel, ConfigDict

from alleleiq.config import Settings, get_settings
from alleleiq.exceptions import GnomADError
from alleleiq.logging import get_logger
from alleleiq.models.population import AncestryGroup, PopulationFrequency
from alleleiq.models.variant import GenomicCoordinates

from ._cache import TTLCache, is_cache_miss
from .base import BaseClient

logger = get_logger(__name__)

_VARIANT_QUERY = """
query AlleleIQVariantQuery($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    genome { ac an ac_hom populations { id ac an ac_hom } }
    exome { ac an ac_hom populations { id ac an ac_hom } }
  }
}
"""

_GENE_CONSTRAINT_QUERY = """
query AlleleIQGeneConstraintQuery($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
  gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
    gene_id
    symbol
    gnomad_constraint {
      pli
      oe_lof
      oe_lof_upper
    }
  }
}
"""


class GeneConstraint(BaseModel):
    """gnomAD gene-level constraint metrics used by the PVS1 criterion.

    ``loeuf`` is gnomAD's ``oe_lof_upper``: the upper bound of the 90%
    confidence interval on the observed/expected ratio of loss-of-function
    variants. Lower LOEUF and higher pLI both indicate stronger
    intolerance to loss-of-function variation (haploinsufficiency).
    """

    model_config = ConfigDict(frozen=True)

    gene_symbol: str
    pli: float | None = None
    oe_lof: float | None = None
    loeuf: float | None = None


class GnomADClient(BaseClient):
    """Async client for the public gnomAD v4 GraphQL API.

    Population frequency and gene constraint lookups are cached in a
    process-wide, class-level :class:`~alleleiq.clients._cache.TTLCache`
    (default TTL: ``ALLELEIQ_GNOMAD_CACHE_TTL_SECONDS``, 1 hour) keyed on
    variant/gene *and* dataset — a dataset change (e.g. ``gnomad_r4`` to
    a future ``gnomad_r5``) is a different cache key, so no explicit
    invalidation step is needed when the configured dataset changes.
    Shared across instances deliberately: :mod:`alleleiq.pipeline`
    constructs a fresh client per variant, so a per-instance cache would
    never help across a multi-variant VCF.
    """

    _cache: ClassVar[TTLCache] = TTLCache()

    def __init__(self, settings: Settings | None = None) -> None:
        """Initialise the client from application settings.

        Args:
            settings: Configuration to use; defaults to
                :func:`alleleiq.config.get_settings`.
        """
        settings = settings or get_settings()
        super().__init__(
            base_url=settings.gnomad_api_url,
            timeout_seconds=settings.gnomad_timeout_seconds,
            max_retries=settings.gnomad_max_retries,
            min_request_interval_seconds=settings.gnomad_min_request_interval_seconds,
        )
        self._dataset = settings.gnomad_dataset
        self._cache_ttl_seconds = settings.gnomad_cache_ttl_seconds

    @classmethod
    def cache_stats(cls) -> dict[str, int]:
        """Return the shared cache's hit/miss/size counters, for shutdown logging."""
        return {"hits": cls._cache.hits, "misses": cls._cache.misses, "size": cls._cache.size}

    @classmethod
    def clear_cache(cls) -> None:
        """Clear the shared cache. Primarily for test isolation between test cases."""
        cls._cache.clear()

    async def fetch_population_frequencies(
        self, coordinates: GenomicCoordinates
    ) -> dict[AncestryGroup, PopulationFrequency]:
        """Fetch per-ancestry allele frequencies for one biallelic variant.

        Args:
            coordinates: A single, already-normalised biallelic variant.
                Multi-allelic sites must be split by the caller (see
                :mod:`alleleiq.parsers.vcf_parser`) before reaching this
                client, since gnomAD variant IDs are inherently biallelic.

        Returns:
            A mapping from ancestry group to :class:`PopulationFrequency`.
            Ancestry groups entirely absent from gnomAD's response (no
            genome *and* no exome data) are omitted from the returned
            mapping; callers should treat a missing key the same as
            ``allele_number=0`` when assessing coverage.

        Raises:
            GnomADError: On a network failure, GraphQL error, malformed
                response, or an unsupported (mitochondrial) variant.
        """
        if coordinates.chrom == "MT":
            raise GnomADError(
                "gnomAD queries for mitochondrial variants are not supported in this "
                "release: mitochondrial allele frequencies use a heteroplasmy-based "
                "schema incompatible with the diploid AC/AN model used here."
            )

        variant_id = coordinates.variant_id
        cache_key = f"frequencies:{self._dataset}:{variant_id}"
        cached = self._cache.get(cache_key)
        if not is_cache_miss(cached):
            logger.debug("gnomAD population frequency cache hit", extra={"variant_id": variant_id})
            result: dict[AncestryGroup, PopulationFrequency] = cached
            return result

        body = await self._execute_query(
            _VARIANT_QUERY, {"variantId": variant_id, "datasetId": self._dataset}
        )

        variant_data = body.get("data", {}).get("variant")
        if variant_data is None:
            logger.info(
                "variant not found in gnomAD",
                extra={"dataset": self._dataset, "chrom": coordinates.chrom},
            )
            self._cache.set(cache_key, {}, self._cache_ttl_seconds)
            return {}

        frequencies = self._parse_population_frequencies(variant_data)
        self._cache.set(cache_key, frequencies, self._cache_ttl_seconds)
        return frequencies

    async def fetch_gene_constraint(self, gene_symbol: str) -> GeneConstraint | None:
        """Fetch gnomAD gene-level constraint metrics (pLI, LOEUF) for a gene.

        Args:
            gene_symbol: HGNC gene symbol, e.g. ``"BRCA1"``.

        Returns:
            A :class:`GeneConstraint`, or None if the gene is not found in
            gnomAD's constraint table.

        Raises:
            GnomADError: On a network failure, GraphQL error, or malformed response.
        """
        cache_key = f"constraint:{self._dataset}:{gene_symbol}"
        cached = self._cache.get(cache_key)
        if not is_cache_miss(cached):
            logger.debug("gnomAD gene constraint cache hit", extra={"gene_symbol": gene_symbol})
            constraint_result: GeneConstraint | None = cached
            return constraint_result

        body = await self._execute_query(
            _GENE_CONSTRAINT_QUERY,
            {"geneSymbol": gene_symbol, "referenceGenome": "GRCh38"},
        )
        gene_data = body.get("data", {}).get("gene")
        if gene_data is None:
            self._cache.set(cache_key, None, self._cache_ttl_seconds)
            return None

        constraint = gene_data.get("gnomad_constraint") or {}
        result = GeneConstraint(
            gene_symbol=gene_symbol,
            pli=constraint.get("pli"),
            oe_lof=constraint.get("oe_lof"),
            loeuf=constraint.get("oe_lof_upper"),
        )
        self._cache.set(cache_key, result, self._cache_ttl_seconds)
        return result

    async def health_check(self) -> bool:
        """Check gnomAD API reachability with a minimal introspection query."""
        try:
            response = await self._request(
                "POST", "", json={"query": "{ __typename }"}, raise_for_status=True
            )
        except httpx.HTTPError:
            return False
        return response.status_code == httpx.codes.OK

    async def _execute_query(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        """Execute a GraphQL query and return the decoded JSON body.

        Raises:
            GnomADError: On network failure, non-JSON response, or a
                GraphQL-level error in the response body.
        """
        try:
            response = await self._request(
                "POST", "", json={"query": query, "variables": variables}
            )
        except httpx.HTTPError as exc:
            raise GnomADError(f"gnomAD request failed: {exc}") from exc

        try:
            body: dict[str, Any] = response.json()
        except ValueError as exc:
            raise GnomADError("gnomAD returned a non-JSON response") from exc

        if body.get("errors"):
            raise GnomADError(f"gnomAD GraphQL error: {body['errors']}")

        return body

    @staticmethod
    def _parse_population_frequencies(
        variant_data: dict[str, Any],
    ) -> dict[AncestryGroup, PopulationFrequency]:
        """Merge gnomAD's separate genome/exome cohorts into one per-ancestry frequency.

        gnomAD reports whole-genome and exome sequencing cohorts as
        separate breakdowns of the same population; summing their AC/AN
        per population is the correct combined estimate, in contrast to
        summing *across* populations, which this function never does.
        """
        genome = variant_data.get("genome") or {}
        exome = variant_data.get("exome") or {}
        genome_by_id = {p["id"]: p for p in (genome.get("populations") or [])}
        exome_by_id = {p["id"]: p for p in (exome.get("populations") or [])}

        result: dict[AncestryGroup, PopulationFrequency] = {}
        for ancestry in AncestryGroup:
            genome_pop = genome_by_id.get(ancestry.value)
            exome_pop = exome_by_id.get(ancestry.value)
            if genome_pop is None and exome_pop is None:
                continue

            allele_count = (genome_pop["ac"] if genome_pop else 0) + (
                exome_pop["ac"] if exome_pop else 0
            )
            allele_number = (genome_pop["an"] if genome_pop else 0) + (
                exome_pop["an"] if exome_pop else 0
            )
            homozygote_count = (genome_pop["ac_hom"] if genome_pop else 0) + (
                exome_pop["ac_hom"] if exome_pop else 0
            )
            result[ancestry] = PopulationFrequency(
                ancestry=ancestry,
                allele_count=allele_count,
                allele_number=allele_number,
                homozygote_count=homozygote_count,
            )
        return result
