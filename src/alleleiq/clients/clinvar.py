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

from datetime import date
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

    @staticmethod
    def _parse_variation_archive(root: Element) -> ClinVarSubmission | None:
        archive = root if root.tag == "VariationArchive" else root.find(".//VariationArchive")
        if archive is None:
            return None

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

        return ClinVarSubmission(
            accession=accession,
            clinical_significance=clinical_significance,
            review_status=review_status,
            star_rating=star_rating,
            submitter_count=submitter_count,
            condition=condition,
            last_evaluated=last_evaluated,
        )
