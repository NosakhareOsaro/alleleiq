"""ClinVar (NCBI E-utilities) client: curated variant classifications as evidence.

A ClinVar classification is evidence to weigh, not ground truth: a
single 0- or 1-star submission from one lab is weaker supporting
evidence than a 3-star expert-panel review, and AlleleIQ's scoring layer
(Phase 3) treats them accordingly via :attr:`ClinVarSubmission.star_rating`
and the configurable ``clinvar_min_stars_strong_evidence`` threshold.
This client only fetches and structures the data; it does not decide how
much weight a given star rating deserves.

Integration note: this client targets NCBI's documented VCV
(VariationArchive) ``efetch`` XML schema. NCBI's ClinVar XML schema has
changed across versions and this implementation should be re-verified
against live responses during integration testing; parsing is written
defensively (missing elements yield ``None``/``0`` rather than raising)
so schema drift degrades gracefully instead of crashing the pipeline.
"""

from __future__ import annotations

import re
from datetime import date
from typing import TypedDict
from xml.etree.ElementTree import Element

import httpx
from defusedxml import ElementTree
from defusedxml.ElementTree import ParseError
from pydantic import BaseModel, ConfigDict, Field

from alleleiq.config import Settings, get_settings
from alleleiq.exceptions import ClinVarError
from alleleiq.logging import get_logger

from .base import BaseClient

logger = get_logger(__name__)

# ClinVar's official review-status-to-star-rating mapping.
# https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
_REVIEW_STATUS_STAR_RATING: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
}


class _ClassificationFields(TypedDict):
    """Fields shared by :class:`ClinVarSubmission` and :class:`ClinVarPositionalMatch`."""

    accession: str
    clinical_significance: str
    review_status: str
    star_rating: int
    submitter_count: int
    condition: str | None
    last_evaluated: date | None


class ClinVarSubmission(BaseModel):
    """A ClinVar record's classification summary for one variant."""

    model_config = ConfigDict(frozen=True)

    accession: str = Field(description="VCV accession, e.g. 'VCV000017661'.")
    clinical_significance: str
    review_status: str
    star_rating: int
    submitter_count: int
    condition: str | None = None
    last_evaluated: date | None = None


# ClinVar's ``VariationName`` embeds the protein change parenthetically, e.g.
# "NM_007294.4(BRCA1):c.5096G>A (p.Arg1699Gln)"; this also matches VEP's
# transcript-prefixed HGVSp form ("ENSP00000269305.4:p.Arg175His").
_PROTEIN_CHANGE_PATTERN = re.compile(r"p\.([A-Za-z0-9_*=]+)")


class ClinVarPositionalMatch(ClinVarSubmission):
    """A classified ClinVar record found near a queried genomic position.

    Used by the PS1/PM5 (same-position pathogenicity) criteria — see
    :mod:`alleleiq.scoring.criteria.same_position` — to compare this
    variant's amino-acid change against other classified variants at or
    near the same position, something a single-variant lookup by rsID or
    HGVS (:meth:`ClinVarClient.fetch_classification`) cannot answer.
    """

    variation_name: str = Field(description="ClinVar's ``VariationName``, e.g. the HGVS summary.")
    protein_change: str | None = Field(
        default=None,
        description="Bare protein-change token extracted from variation_name, e.g. 'Arg175His'.",
    )
    genomic_position: int | None = Field(
        default=None, description="GRCh38 VCF-style position, if reported by ClinVar."
    )


class ClinVarClient(BaseClient):
    """Async client for NCBI's ClinVar E-utilities (esearch + efetch)."""

    def __init__(self, settings: Settings | None = None) -> None:
        """Initialise the client from application settings.

        Args:
            settings: Configuration to use; defaults to
                :func:`alleleiq.config.get_settings`.
        """
        settings = settings or get_settings()
        super().__init__(
            base_url=settings.clinvar_api_url,
            timeout_seconds=settings.clinvar_timeout_seconds,
            max_retries=settings.clinvar_max_retries,
            min_request_interval_seconds=settings.clinvar_min_request_interval_seconds,
        )
        self._api_key = settings.clinvar_api_key

    async def fetch_classification(
        self,
        *,
        rsid: str | None = None,
        gene_symbol: str | None = None,
        hgvs_c: str | None = None,
    ) -> ClinVarSubmission | None:
        """Search ClinVar for a variant and fetch its classification summary.

        Args:
            rsid: dbSNP rsID (e.g. from VEP's ``existing_variant_id``), if
                known. Preferred over HGVS search when available, since
                it is unambiguous.
            gene_symbol: HGNC gene symbol, used with ``hgvs_c`` when no
                rsID is available.
            hgvs_c: HGVS coding-sequence notation, used with
                ``gene_symbol`` when no rsID is available.

        Returns:
            A :class:`ClinVarSubmission`, or None if no matching ClinVar
            record was found (a novel or unsubmitted variant — a
            legitimate outcome, not an error).

        Raises:
            ClinVarError: On a network failure, malformed response, or
                if no search terms could be constructed from the given
                arguments.
        """
        term = self._build_search_term(rsid=rsid, gene_symbol=gene_symbol, hgvs_c=hgvs_c)
        uid = await self._esearch(term)
        if uid is None:
            logger.info("no ClinVar record found", extra={"had_rsid": rsid is not None})
            return None
        return await self._efetch(uid)

    async def search_same_position(
        self, *, chrom: str, start: int, end: int, max_hits: int = 20
    ) -> list[ClinVarPositionalMatch]:
        """Search ClinVar for classified variants overlapping a genomic position window.

        Used by the PS1/PM5 same-position criteria to find other
        classified variants near this one, so their amino-acid changes
        can be compared. Queries NCBI's documented ``chrpos38`` genomic
        position field (GRCh38-only, matching AlleleIQ's own scope) —
        not the informal "position[Genomic Location]" field name
        sometimes seen in casual descriptions of ClinVar search, which
        is not a real Entrez field.

        Args:
            chrom: Chromosome without a 'chr' prefix, e.g. '17'.
            start: Inclusive GRCh38 range start (1-based, VCF convention).
            end: Inclusive GRCh38 range end.
            max_hits: Maximum number of ClinVar records to fetch details
                for, bounding both the esearch page size and the
                subsequent batch efetch request.

        Returns:
            Classified matches in the window, each with clinical
            significance, star rating, and (when parseable) the protein
            change. Empty if nothing is classified in the window — a
            legitimate outcome, not an error.

        Raises:
            ClinVarError: On a network failure or malformed response.
        """
        term = f"{chrom}[chr] AND {start}:{end}[chrpos38]"
        try:
            response = await self._request(
                "GET",
                "/esearch.fcgi",
                params=self._with_auth(
                    {
                        "db": "clinvar",
                        "term": term,
                        "retmode": "json",
                        "retmax": str(max_hits),
                    }
                ),
            )
        except httpx.HTTPError as exc:
            raise ClinVarError(f"ClinVar positional esearch failed: {exc}") from exc

        try:
            body = response.json()
        except ValueError as exc:
            raise ClinVarError("ClinVar positional esearch returned a non-JSON response") from exc

        uids = body.get("esearchresult", {}).get("idlist", [])
        if not uids:
            return []
        return await self._efetch_many(uids)

    async def health_check(self) -> bool:
        """Check NCBI E-utilities reachability via a minimal esearch call."""
        try:
            response = await self._request(
                "GET", "/esearch.fcgi", params=self._with_auth({"db": "clinvar", "term": "1"})
            )
        except httpx.HTTPError:
            return False
        return response.status_code == httpx.codes.OK

    def _with_auth(self, params: dict[str, str]) -> dict[str, str]:
        if self._api_key:
            return {**params, "api_key": self._api_key}
        return params

    @staticmethod
    def _build_search_term(*, rsid: str | None, gene_symbol: str | None, hgvs_c: str | None) -> str:
        if rsid:
            return rsid
        if gene_symbol and hgvs_c:
            return f"{gene_symbol}[gene] AND {hgvs_c}"
        raise ClinVarError(
            "cannot search ClinVar without either an rsID or both gene_symbol and hgvs_c"
        )

    async def _esearch(self, term: str) -> str | None:
        try:
            response = await self._request(
                "GET",
                "/esearch.fcgi",
                params=self._with_auth({"db": "clinvar", "term": term, "retmode": "json"}),
            )
        except httpx.HTTPError as exc:
            raise ClinVarError(f"ClinVar esearch failed: {exc}") from exc

        try:
            body = response.json()
        except ValueError as exc:
            raise ClinVarError("ClinVar esearch returned a non-JSON response") from exc

        id_list = body.get("esearchresult", {}).get("idlist", [])
        return id_list[0] if id_list else None

    async def _efetch(self, uid: str) -> ClinVarSubmission | None:
        try:
            response = await self._request(
                "GET",
                "/efetch.fcgi",
                params=self._with_auth(
                    {"db": "clinvar", "id": uid, "rettype": "vcv", "is_variationid": "true"}
                ),
            )
        except httpx.HTTPError as exc:
            raise ClinVarError(f"ClinVar efetch failed: {exc}") from exc

        try:
            root = ElementTree.fromstring(response.text)
        except ParseError as exc:
            raise ClinVarError("ClinVar efetch returned malformed XML") from exc

        return self._parse_variation_archive(root)

    async def _efetch_many(self, uids: list[str]) -> list[ClinVarPositionalMatch]:
        """Batch-fetch classification details for multiple ClinVar UIDs in one request."""
        try:
            response = await self._request(
                "GET",
                "/efetch.fcgi",
                params=self._with_auth(
                    {
                        "db": "clinvar",
                        "id": ",".join(uids),
                        "rettype": "vcv",
                        "is_variationid": "true",
                    }
                ),
            )
        except httpx.HTTPError as exc:
            raise ClinVarError(f"ClinVar batch efetch failed: {exc}") from exc

        try:
            root = ElementTree.fromstring(response.text)
        except ParseError as exc:
            raise ClinVarError("ClinVar batch efetch returned malformed XML") from exc

        matches = []
        for archive in root.findall(".//VariationArchive"):
            match = self._parse_positional_archive(archive)
            if match is not None:
                matches.append(match)
        return matches

    @classmethod
    def _parse_variation_archive(cls, root: Element) -> ClinVarSubmission | None:
        archive = root if root.tag == "VariationArchive" else root.find(".//VariationArchive")
        if archive is None:
            return None
        fields = cls._extract_classification_fields(archive)
        if fields is None:
            return None
        return ClinVarSubmission(**fields)

    @classmethod
    def _parse_positional_archive(cls, archive: Element) -> ClinVarPositionalMatch | None:
        fields = cls._extract_classification_fields(archive)
        if fields is None:
            return None
        variation_name = archive.get("VariationName") or ""
        protein_match = _PROTEIN_CHANGE_PATTERN.search(variation_name)
        genomic_position = cls._extract_grch38_position(archive)
        return ClinVarPositionalMatch(
            **fields,
            variation_name=variation_name,
            protein_change=protein_match.group(1) if protein_match else None,
            genomic_position=genomic_position,
        )

    @staticmethod
    def _extract_classification_fields(archive: Element) -> _ClassificationFields | None:
        """Fields shared by both single-record and positional-search parsing."""
        classification = archive.find(".//Classifications/GermlineClassification")
        if classification is None:
            return None

        clinical_significance = classification.findtext("Description") or "not provided"
        review_status = classification.findtext("ReviewStatus") or "no assertion provided"
        star_rating = _REVIEW_STATUS_STAR_RATING.get(review_status.strip().lower(), 0)

        last_evaluated_raw = classification.get("DateLastEvaluated")
        last_evaluated = None
        if last_evaluated_raw:
            try:
                last_evaluated = date.fromisoformat(last_evaluated_raw)
            except ValueError:
                logger.warning(
                    "could not parse ClinVar DateLastEvaluated",
                    extra={"raw_value": last_evaluated_raw},
                )

        submitter_count = len(archive.findall(".//ClinicalAssertionList/ClinicalAssertion"))

        condition = None
        trait_name = archive.find(".//TraitSet/Trait/Name/ElementValue[@Type='Preferred']")
        if trait_name is not None and trait_name.text:
            condition = trait_name.text

        accession = archive.get("Accession") or f"VCV-UNKNOWN-{archive.get('VariationID', '')}"

        return _ClassificationFields(
            accession=accession,
            clinical_significance=clinical_significance,
            review_status=review_status,
            star_rating=star_rating,
            submitter_count=submitter_count,
            condition=condition,
            last_evaluated=last_evaluated,
        )

    @staticmethod
    def _extract_grch38_position(archive: Element) -> int | None:
        location = archive.find(".//Location/SequenceLocation[@Assembly='GRCh38']")
        if location is None:
            return None
        raw = location.get("positionVCF") or location.get("start")
        if raw is None:
            return None
        try:
            return int(raw)
        except ValueError:
            return None
