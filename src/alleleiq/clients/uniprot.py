"""UniProt REST client: annotated functional domains, active sites, and binding sites.

Used by PM1 (:mod:`alleleiq.scoring.criteria.hotspot`) to check whether a
variant's amino-acid position falls within a well-established functional
region. Queried by gene symbol directly against UniProt's reviewed
(Swiss-Prot) human entries — verified live (2026-07) against
``rest.uniprot.org``:
``/uniprotkb/search?query=gene:{SYMBOL}+AND+organism_id:9606+AND+reviewed:true``
returns a ``results`` array whose first (best) hit's ``features`` list
contains ``{"type": "Domain" | "Active site" | "Binding site", "location":
{"start": {"value": int}, "end": {"value": int}}, "description": str}``
entries, 1-based amino-acid positions.
"""

from __future__ import annotations

from typing import Any, ClassVar

import httpx
from pydantic import BaseModel, ConfigDict

from alleleiq.config import Settings, get_settings
from alleleiq.exceptions import UniProtError
from alleleiq.logging import get_logger

from ._cache import TTLCache, is_cache_miss
from .base import BaseClient

logger = get_logger(__name__)

# Domain annotations for a given gene do not change within a process
# lifetime; a generous fixed TTL avoids re-querying UniProt for every
# variant in a multi-variant VCF that happens to share a gene, the same
# rationale GnomADClient's constraint cache uses.
_DOMAIN_CACHE_TTL_SECONDS = 3600.0

_HUMAN_TAXON_ID = 9606

# The three feature categories PM1's spec names: annotated domains, active
# sites, and binding sites. Broader categories (e.g. "Region", "Repeat") are
# deliberately excluded — they range from well-established to speculative,
# and PM1 needs annotation confident enough to treat as "well-established
# functional domain", not merely "protein has been annotated at all".
_CRITICAL_FEATURE_TYPES = frozenset({"Domain", "Active site", "Binding site"})


class ProteinRegion(BaseModel):
    """A single amino-acid-position range with a description of what it is."""

    model_config = ConfigDict(frozen=True)

    start_aa: int
    end_aa: int
    description: str
    source: str


class UniProtClient(BaseClient):
    """Async client for the UniProt REST API (reviewed human entries only).

    Domain lookups are cached process-wide (see
    :class:`alleleiq.clients._cache.TTLCache`), the same pattern
    :class:`alleleiq.clients.gnomad.GnomADClient` uses for gene
    constraint, since :mod:`alleleiq.pipeline` constructs a fresh client
    per variant and a per-instance cache would never help across a
    multi-variant VCF.
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
            base_url=settings.uniprot_api_url,
            timeout_seconds=settings.uniprot_timeout_seconds,
            max_retries=settings.uniprot_max_retries,
        )

    @classmethod
    def cache_stats(cls) -> dict[str, int]:
        """Return process-wide cache hit/miss/size statistics, for startup/shutdown logging."""
        return {"hits": cls._cache.hits, "misses": cls._cache.misses, "size": cls._cache.size}

    @classmethod
    def clear_cache(cls) -> None:
        """Clear the process-wide cache. Primarily for test isolation."""
        cls._cache.clear()

    async def fetch_domains(self, gene_symbol: str) -> list[ProteinRegion]:
        """Fetch annotated domain/active-site/binding-site regions for a gene.

        Args:
            gene_symbol: HGNC gene symbol.

        Returns:
            Critical-region annotations for the gene's best-matching
            reviewed human UniProt entry. Empty if there is no reviewed
            entry, or it has no domain/active-site/binding-site
            annotations — a legitimate outcome, not an error; a gene
            with genuinely sparse UniProt curation is common and does
            not imply "has no functional domains".

        Raises:
            UniProtError: On a network failure or malformed response.
        """
        cache_key = gene_symbol.upper()
        cached = self._cache.get(cache_key)
        if not is_cache_miss(cached):
            return list(cached)

        try:
            response = await self._request(
                "GET",
                "/uniprotkb/search",
                params={
                    "query": f"gene:{gene_symbol} AND organism_id:{_HUMAN_TAXON_ID} "
                    "AND reviewed:true",
                    "fields": "accession,gene_names,ft_domain,ft_binding,ft_act_site",
                    "format": "json",
                    "size": "1",
                },
            )
        except httpx.HTTPError as exc:
            raise UniProtError(f"UniProt search failed for gene {gene_symbol!r}: {exc}") from exc

        try:
            payload = response.json()
        except ValueError as exc:
            raise UniProtError("UniProt search returned a non-JSON response") from exc

        results = payload.get("results")
        if not isinstance(results, list) or not results:
            logger.info("no reviewed UniProt entry found", extra={"gene_symbol": gene_symbol})
            self._cache.set(cache_key, [], _DOMAIN_CACHE_TTL_SECONDS)
            return []

        regions = self._parse_features(results[0], gene_symbol)
        self._cache.set(cache_key, regions, _DOMAIN_CACHE_TTL_SECONDS)
        return regions

    async def health_check(self) -> bool:
        """Check UniProt REST API reachability via a minimal search."""
        try:
            response = await self._request(
                "GET",
                "/uniprotkb/search",
                params={"query": f"gene:BRCA1 AND organism_id:{_HUMAN_TAXON_ID}", "size": "1"},
            )
        except httpx.HTTPError:
            return False
        return response.status_code == httpx.codes.OK

    @staticmethod
    def _parse_features(entry: dict[str, Any], gene_symbol: str) -> list[ProteinRegion]:
        accession = entry.get("primaryAccession", "unknown")
        regions: list[ProteinRegion] = []
        for feature in entry.get("features") or []:
            feature_type = feature.get("type")
            if feature_type not in _CRITICAL_FEATURE_TYPES:
                continue
            location = feature.get("location") or {}
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")
            if start is None or end is None:
                continue
            description = feature.get("description") or feature_type
            regions.append(
                ProteinRegion(
                    start_aa=int(start),
                    end_aa=int(end),
                    description=f"{feature_type}: {description}",
                    source=f"UniProt {accession} ({gene_symbol})",
                )
            )
        return regions
