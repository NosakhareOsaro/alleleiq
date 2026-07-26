"""Local (tabix-indexed) and remote (VEP REST) REVEL score lookup.

The Ensembl VEP REST API *can* serve REVEL scores (see
:mod:`alleleiq.clients.vep`), but its ``dbNSFP`` plugin coverage is
incomplete for many positions, so PP3/BP4 frequently fall back to
``requires_clinical_input`` for lack of any computational signal. A
locally indexed copy of the full REVEL precomputed-scores dataset (see
``scripts/download_revel.sh``) has complete genome-wide coverage and is
queried in microseconds via ``tabix``, so :class:`REVELClient` treats it
as the primary source and only falls back to VEP's REST-served score —
sourced from an already-fetched :class:`~alleleiq.clients.vep.VEPAnnotation`
where the caller has one, or via a fresh VEP request otherwise — when no
local file is configured or the position isn't found in it.

Which source produced a score is recorded on :class:`REVELScore` and
logged, both because it is a legitimate audit-trail fact (which database
was consulted) and because "local" and "vep_rest" carry different
completeness guarantees a reviewer may care about.
"""

from __future__ import annotations

import asyncio
from pathlib import Path
from typing import Literal

import pysam
from pydantic import BaseModel, ConfigDict

from alleleiq.clients.vep import VEPAnnotation, VEPClient
from alleleiq.config import Settings, get_settings
from alleleiq.logging import get_logger
from alleleiq.models.variant import GenomicCoordinates

logger = get_logger(__name__)

REVELSource = Literal["local", "vep_rest", "unavailable"]


class REVELScore(BaseModel):
    """The result of a REVEL score lookup, with its provenance."""

    model_config = ConfigDict(frozen=True)

    score: float | None
    source: REVELSource


class REVELClient:
    """Resolves a variant's REVEL score from a local tabix index, falling back to VEP REST.

    Unlike :class:`~alleleiq.clients.base.BaseClient` subclasses, this is
    not primarily an HTTP client: its primary data path
    (:class:`pysam.TabixFile`) is a local, synchronous file read, so it
    is not an async context manager. It only opens a network connection
    (via an internal :class:`~alleleiq.clients.vep.VEPClient`) on the
    fallback path, and only when the caller has not already supplied a
    :class:`~alleleiq.clients.vep.VEPAnnotation` from its own VEP query.
    """

    def __init__(self, settings: Settings | None = None) -> None:
        """Initialise the client from application settings.

        Args:
            settings: Configuration to use; defaults to
                :func:`alleleiq.config.get_settings`. ``settings.revel_data_path``,
                if set and pointing at an existing, tabix-indexed file,
                is used as the local lookup source.
        """
        settings = settings or get_settings()
        self._data_path = Path(settings.revel_data_path) if settings.revel_data_path else None
        self._settings = settings
        self._tabix: pysam.TabixFile | None = None
        self._tabix_unavailable = False

    def close(self) -> None:
        """Release the underlying tabix file handle, if one was opened."""
        if self._tabix is not None:
            self._tabix.close()
            self._tabix = None

    async def get_score(
        self,
        coordinates: GenomicCoordinates,
        *,
        vep_annotation: VEPAnnotation | None = None,
    ) -> REVELScore:
        """Resolve a REVEL score for one variant.

        Args:
            coordinates: The variant to look up.
            vep_annotation: A VEP annotation already fetched for this
                variant elsewhere in the pipeline (e.g. by
                :func:`alleleiq.pipeline.build_evidence_report`). When
                given, its ``revel_score`` is used as the fallback
                without issuing a second VEP request. When omitted and
                the local lookup misses, this method makes its own VEP
                REST call.

        Returns:
            A :class:`REVELScore`. ``score`` is None, with
            ``source="unavailable"``, if neither the local index nor VEP
            has a score for this variant — this is a legitimate result,
            not an error.
        """
        local_score = await asyncio.to_thread(self._lookup_local, coordinates)
        if local_score is not None:
            logger.debug(
                "REVEL score resolved from local index",
                extra={"variant_id": coordinates.variant_id, "source": "local"},
            )
            return REVELScore(score=local_score, source="local")

        remote_score = await self._lookup_remote(coordinates, vep_annotation)
        if remote_score is not None:
            logger.debug(
                "REVEL score resolved via VEP REST fallback",
                extra={"variant_id": coordinates.variant_id, "source": "vep_rest"},
            )
            return REVELScore(score=remote_score, source="vep_rest")

        logger.debug(
            "no REVEL score available from local index or VEP",
            extra={"variant_id": coordinates.variant_id, "source": "unavailable"},
        )
        return REVELScore(score=None, source="unavailable")

    async def _lookup_remote(
        self, coordinates: GenomicCoordinates, vep_annotation: VEPAnnotation | None
    ) -> float | None:
        if vep_annotation is not None:
            return vep_annotation.revel_score
        async with VEPClient(self._settings) as vep_client:
            annotation = await vep_client.fetch_annotation(coordinates)
        return annotation.revel_score if annotation else None

    def _lookup_local(self, coordinates: GenomicCoordinates) -> float | None:
        """Synchronous tabix lookup; run off the event loop via ``asyncio.to_thread``."""
        tabix = self._ensure_tabix()
        if tabix is None:
            return None

        for chrom_candidate in self._chrom_candidates(coordinates.chrom):
            try:
                rows = tabix.fetch(chrom_candidate, coordinates.pos - 1, coordinates.pos)
            except ValueError:
                # Contig not present under this naming candidate; try the next one.
                continue
            for row in rows:
                fields = row.split("\t")
                if len(fields) < 5:
                    continue
                _chrom, pos_str, ref, alt, score_str = fields[:5]
                if int(pos_str) != coordinates.pos:
                    continue
                if ref != coordinates.ref or alt != coordinates.alt:
                    continue
                try:
                    return float(score_str)
                except ValueError:
                    return None
        return None

    def _ensure_tabix(self) -> pysam.TabixFile | None:
        if self._tabix_unavailable:
            return None
        if self._tabix is not None:
            return self._tabix
        if self._data_path is None or not self._data_path.exists():
            self._tabix_unavailable = True
            return None

        try:
            self._tabix = pysam.TabixFile(str(self._data_path))
        except (OSError, ValueError) as exc:
            logger.warning(
                "REVEL local index could not be opened; falling back to VEP REST",
                extra={"path": str(self._data_path), "error": str(exc)},
            )
            self._tabix_unavailable = True
            return None
        return self._tabix

    @staticmethod
    def _chrom_candidates(chrom: str) -> tuple[str, ...]:
        """Try the contig name as-is, then with/without a 'chr' prefix.

        ``GenomicCoordinates.chrom`` is already normalised without a
        'chr' prefix (see :mod:`alleleiq.models.variant`), but the local
        REVEL index's contig naming depends on how
        ``scripts/download_revel.sh`` built it, so both forms are tried.
        """
        if chrom.lower().startswith("chr"):
            bare = chrom[3:]
            return (chrom, bare)
        return (chrom, f"chr{chrom}")
