"""Ensembl VEP REST client: consequence, transcript, REVEL, and SpliceAI annotation.

Field shapes verified directly against the live ``rest.ensembl.org``
API (2026-07): with ``dbNSFP=REVEL_score`` and ``SpliceAI=1`` set, the
public endpoint *does* serve both, at the transcript-consequence level
(not nested under a ``dbnsfp`` sub-object, despite that being the
plugin's name) — ``numbers=1`` is required to get the ``exon`` field
PVS1 depends on. ``revel_score`` comes back as a comma-separated string
of per-RefSeq-transcript-match scores with ``.`` placeholders for
transcripts lacking a score (e.g. ``".,.,0.229,.,0.229"``), not a
single float, so it must be parsed accordingly.
"""

from __future__ import annotations

from typing import Any, ClassVar

import httpx
from pydantic import BaseModel, ConfigDict, Field

from alleleiq.config import Settings, get_settings
from alleleiq.exceptions import VEPError
from alleleiq.models.variant import GenomicCoordinates

from .base import BaseClient

_ANNOTATION_PARAMS = {
    "canonical": "1",
    "mane": "1",
    "hgvs": "1",
    "numbers": "1",
    "dbNSFP": "REVEL_score",
    "SpliceAI": "1",
}


class VEPAnnotation(BaseModel):
    """Functional annotation for a single variant from Ensembl VEP.

    ``gene_symbol`` through ``spliceai_ds_dl`` are drawn from the
    *primary* transcript consequence: the canonical transcript if VEP
    flags one, else the MANE Select transcript, else the first reported
    consequence. gnomAD's pLI/LOEUF gene constraint is deliberately not
    part of this model — see
    :meth:`alleleiq.clients.gnomad.GnomADClient.fetch_gene_constraint` —
    since a VEP client has no business querying gnomAD; callers combine
    both when evaluating criteria that need both (e.g. PVS1).
    """

    model_config = ConfigDict(frozen=True)

    variant_id: str
    most_severe_consequence: str | None = None
    gene_symbol: str | None = None
    gene_id: str | None = None
    transcript_id: str | None = None
    is_canonical: bool = False
    is_mane_select: bool = False
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    exon: str | None = Field(
        default=None, description="Exon position as reported by VEP, e.g. '18/23'."
    )
    revel_score: float | None = None
    spliceai_ds_ag: float | None = None
    spliceai_ds_al: float | None = None
    spliceai_ds_dg: float | None = None
    spliceai_ds_dl: float | None = None
    existing_variant_id: str | None = Field(
        default=None, description="Colocated dbSNP rsID, if any, per VEP's cache."
    )


class VEPClient(BaseClient):
    """Async client for the Ensembl VEP REST API."""

    _release: ClassVar[str | None] = None

    def __init__(self, settings: Settings | None = None) -> None:
        """Initialise the client from application settings.

        Args:
            settings: Configuration to use; defaults to
                :func:`alleleiq.config.get_settings`.
        """
        settings = settings or get_settings()
        super().__init__(
            base_url=settings.vep_api_url,
            timeout_seconds=settings.vep_timeout_seconds,
            max_retries=settings.vep_max_retries,
        )

    async def fetch_annotation(self, coordinates: GenomicCoordinates) -> VEPAnnotation | None:
        """Fetch VEP annotation for a single variant.

        Args:
            coordinates: The variant to annotate.

        Returns:
            A :class:`VEPAnnotation`, or None if VEP returned no
            annotation for this variant (e.g. it falls outside any
            annotated region).

        Raises:
            VEPError: On a network failure, malformed response, or VEP
                server error.
        """
        results = await self.fetch_annotations([coordinates])
        return results.get(coordinates.variant_id)

    async def fetch_annotations(
        self, coordinates_list: list[GenomicCoordinates]
    ) -> dict[str, VEPAnnotation]:
        """Fetch VEP annotations for a batch of variants in one request.

        Args:
            coordinates_list: Variants to annotate.

        Returns:
            A mapping from ``variant_id`` (``chrom-pos-ref-alt``) to
            :class:`VEPAnnotation`. Variants VEP could not annotate are
            simply absent from the mapping.

        Raises:
            VEPError: On a network failure, malformed response, or VEP
                server error.
        """
        if not coordinates_list:
            return {}

        region_notations = [self._to_region_notation(c) for c in coordinates_list]
        try:
            response = await self._request(
                "POST",
                "/vep/human/region",
                params=_ANNOTATION_PARAMS,
                json={"variants": region_notations},
            )
        except httpx.HTTPError as exc:
            raise VEPError(f"VEP request failed: {exc}") from exc

        try:
            payload = response.json()
        except ValueError as exc:
            raise VEPError("VEP returned a non-JSON response") from exc

        if isinstance(payload, dict) and "error" in payload:
            raise VEPError(f"VEP error: {payload['error']}")
        if not isinstance(payload, list):
            raise VEPError(f"unexpected VEP response shape: {type(payload).__name__}")

        annotations: dict[str, VEPAnnotation] = {}
        for entry in payload:
            annotation = self._parse_entry(entry)
            annotations[annotation.variant_id] = annotation
        return annotations

    async def fetch_release(self) -> str | None:
        """Return the Ensembl release serving VEP annotations, e.g. ``"112"``.

        Cached for the process lifetime (a release does not change under a
        running analysis). Never raises: the release is provenance metadata
        for the report, so a failed lookup is recorded as unknown rather
        than aborting the analysis.
        """
        if VEPClient._release is not None:
            return VEPClient._release
        try:
            response = await self._request("GET", "/info/software")
            release = response.json().get("release")
        except (httpx.HTTPError, ValueError, AttributeError):
            return None
        if release is None:
            return None
        VEPClient._release = str(release)
        return VEPClient._release

    async def health_check(self) -> bool:
        """Check VEP REST API reachability via its ping endpoint."""
        try:
            response = await self._request("GET", "/info/ping")
        except httpx.HTTPError:
            return False
        return response.status_code == httpx.codes.OK

    @staticmethod
    def _to_region_notation(coordinates: GenomicCoordinates) -> str:
        """Render coordinates as VEP's ``CHROM POS ID REF ALT QUAL FILTER INFO`` input line."""
        return f"{coordinates.chrom} {coordinates.pos} . {coordinates.ref} {coordinates.alt} . . ."

    @staticmethod
    def _parse_revel_score(revel_raw: object) -> float | None:
        """Parse VEP's ``revel_score`` field.

        Verified against the live API: this arrives as a comma-separated
        string of per-RefSeq-transcript-match scores, using ``.`` for
        transcripts with no score, e.g. ``".,.,0.229,.,0.229"``. In
        practice, all non-placeholder values for a given variant agree
        (they are redundant mappings to the same underlying prediction),
        so the first parseable value is used.
        """
        if revel_raw is None:
            return None
        if isinstance(revel_raw, (int, float)):
            return float(revel_raw)
        if not isinstance(revel_raw, str):
            return None

        for candidate in revel_raw.split(","):
            candidate = candidate.strip()
            if not candidate or candidate == ".":
                continue
            try:
                return float(candidate)
            except ValueError:
                continue
        return None

    @classmethod
    def _parse_entry(cls, entry: dict[str, Any]) -> VEPAnnotation:
        variant_id = cls._variant_id_from_input(str(entry.get("input", "")))
        transcript_consequences: list[dict[str, Any]] = entry.get("transcript_consequences") or []
        primary = cls._select_primary_transcript(transcript_consequences)

        revel_score: float | None = None
        spliceai: dict[str, Any] = {}
        if primary is not None:
            revel_score = cls._parse_revel_score(primary.get("revel_score"))
            spliceai = primary.get("spliceai") or {}

        return VEPAnnotation(
            variant_id=variant_id,
            most_severe_consequence=entry.get("most_severe_consequence"),
            gene_symbol=primary.get("gene_symbol") if primary else None,
            gene_id=primary.get("gene_id") if primary else None,
            transcript_id=primary.get("transcript_id") if primary else None,
            is_canonical=bool(primary.get("canonical")) if primary else False,
            is_mane_select=bool(primary.get("mane_select")) if primary else False,
            hgvs_c=primary.get("hgvsc") if primary else None,
            hgvs_p=primary.get("hgvsp") if primary else None,
            exon=primary.get("exon") if primary else None,
            revel_score=revel_score,
            spliceai_ds_ag=spliceai.get("DS_AG"),
            spliceai_ds_al=spliceai.get("DS_AL"),
            spliceai_ds_dg=spliceai.get("DS_DG"),
            spliceai_ds_dl=spliceai.get("DS_DL"),
            existing_variant_id=cls._existing_variant_id(entry),
        )

    @staticmethod
    def _existing_variant_id(entry: dict[str, Any]) -> str | None:
        colocated = entry.get("colocated_variants") or []
        if colocated:
            variant_id = colocated[0].get("id")
            return str(variant_id) if variant_id is not None else None
        return None

    @staticmethod
    def _select_primary_transcript(
        consequences: list[dict[str, Any]],
    ) -> dict[str, Any] | None:
        """Pick canonical > MANE Select > first, matching clinical reporting convention."""
        if not consequences:
            return None
        for consequence in consequences:
            if consequence.get("canonical") == 1:
                return consequence
        for consequence in consequences:
            if consequence.get("mane_select"):
                return consequence
        return consequences[0]

    @staticmethod
    def _variant_id_from_input(input_line: str) -> str:
        parts = input_line.split()
        if len(parts) >= 5:
            chrom, pos, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            return f"{chrom}-{pos}-{ref}-{alt}"
        return input_line
