"""Local ClinGen gene-disease validity and dosage-sensitivity lookup.

ClinGen's Gene-Disease Validity curations (which diseases a gene is
established to cause, and how strongly) and its Dosage Sensitivity Map
(whether loss of one copy of a gene is itself pathogenic) are both
static, expert-curated snapshots published as downloadable files, not a
live queryable API — so, like :mod:`alleleiq.clients.revel`,
:class:`ClinGenClient` reads bundled local files rather than making
network requests per variant. The two bundled files
(``data/clingen/gene_disease_validity.csv`` and
``data/clingen/dosage_sensitivity.tsv``) are real snapshots downloaded
from ClinGen's own distribution endpoints (see ``scripts/update_clingen.sh``
for the refresh procedure) and are small enough (~1.4MB combined) to
commit directly, unlike REVEL's multi-hundred-megabyte index.

This client is the evidence source for
:class:`alleleiq.scoring.criteria.lof.PVS1Criterion`'s
haploinsufficiency-mechanism strength tiering: a gene's ClinGen
gene-disease *classification* (DEFINITIVE, STRONG, MODERATE, LIMITED, ...)
combined with its dosage map *haploinsufficiency score* (0-3, plus two
non-numeric special codes — see :meth:`ClinGenClient.get_haploinsufficiency_score`)
together determine how strongly PVS1 should be applied, rather than the
single bundled 25-gene fallback list used before this integration
existed.
"""

from __future__ import annotations

import csv
from datetime import date, datetime
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

from alleleiq.config import Settings, get_settings
from alleleiq.logging import get_logger

logger = get_logger(__name__)

# ClinGen's own classification scale, strongest to weakest, used to pick a
# single "best" classification when a gene has multiple gene-disease rows
# (the curation file is one row per gene-disease pair, not one row per
# gene — e.g. TTN has 7 rows for 7 distinct diseases).
_CLASSIFICATION_PRECEDENCE: dict[str, int] = {
    "Definitive": 0,
    "Strong": 1,
    "Moderate": 2,
    "Limited": 3,
    "No Known Disease Relationship": 4,
    "Disputed": 5,
    "Refuted": 6,
}

# ClinGen's dosage map uses two non-numeric-scale codes alongside the 0-3
# evidence scale, both discovered by inspecting a live download of
# https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
# rather than assumed from documentation:
#   30 = "Gene associated with autosomal recessive phenotype" (dosage
#        sensitivity is not assessable the same way for AR genes; this is
#        not evidence *for* haploinsufficiency and must not be conflated
#        with score 3).
#   40 = "Dosage sensitivity unlikely" (a real assertion, distinct from
#        "0 = no evidence available yet").
# Both are normalised to 0 for get_haploinsufficiency_score()'s 0-3
# contract, since neither supports a haploinsufficiency claim, but the
# raw code and its description are preserved on DosageSensitivity for
# callers that need the distinction.
_HI_SPECIAL_CODES: dict[int, str] = {
    30: "Gene associated with autosomal recessive phenotype",
    40: "Dosage sensitivity unlikely",
}


class ClinGenValidity(BaseModel):
    """One ClinGen gene-disease validity curation (one gene-disease pair)."""

    model_config = ConfigDict(frozen=True)

    gene_symbol: str
    disease_label: str
    mondo_id: str | None = None
    mode_of_inheritance: str | None = None
    classification: str = Field(
        description="One of ClinGen's seven classifications, e.g. 'Definitive', 'Limited'."
    )
    classification_date: date | None = None
    gcep: str = Field(default="", description="The ClinGen Gene Curation Expert Panel name.")


class DosageSensitivity(BaseModel):
    """One ClinGen dosage-sensitivity map row for a gene."""

    model_config = ConfigDict(frozen=True)

    gene_symbol: str
    haploinsufficiency_score: int = Field(
        description="Normalised to 0-3; special codes 30/40 are mapped to 0 (see module docstring)."
    )
    haploinsufficiency_raw_code: int = Field(
        description="The unnormalised code as published (0-3, 30, or 40)."
    )
    haploinsufficiency_description: str = ""
    triplosensitivity_score: int = 0


def _resolve_path(configured: str) -> Path:
    return Path(configured)


class ClinGenClient:
    """Reads ClinGen's bundled gene-disease validity and dosage-sensitivity files.

    Unlike :class:`~alleleiq.clients.base.BaseClient` subclasses, this is
    a local-file reader, not an HTTP client — there is no live ClinGen
    query API for either dataset. Both files are parsed once, on first
    use, and cached in memory for the lifetime of the client instance;
    :mod:`alleleiq.pipeline` constructs one instance per variant, so
    within a single multi-variant VCF run the OS page cache (not this
    client) is what keeps repeat parses cheap. A process-wide cache was
    deliberately not added here, unlike :mod:`alleleiq.clients._cache`'s
    use in :class:`~alleleiq.clients.gnomad.GnomADClient`: these files
    are tiny (~1.4MB combined) and parse in low-single-digit
    milliseconds, so the added complexity of a shared cache was not
    justified.
    """

    def __init__(self, settings: Settings | None = None) -> None:
        """Initialise the client from application settings.

        Args:
            settings: Configuration to use; defaults to
                :func:`alleleiq.config.get_settings`. Reads
                ``settings.clingen_gene_validity_path`` and
                ``settings.clingen_dosage_path``.
        """
        settings = settings or get_settings()
        self._validity_path = _resolve_path(settings.clingen_gene_validity_path)
        self._dosage_path = _resolve_path(settings.clingen_dosage_path)
        self._validity_by_gene: dict[str, list[ClinGenValidity]] | None = None
        self._dosage_by_gene: dict[str, DosageSensitivity] | None = None

    def load_gene_validity(self) -> dict[str, str]:
        """Map every curated gene symbol to its single best classification.

        Returns:
            ``{gene_symbol: classification}``, e.g. ``{"BRCA1": "Definitive"}``.
            When a gene has multiple gene-disease curations, the
            strongest classification wins (see ``_CLASSIFICATION_PRECEDENCE``);
            e.g. a gene classified "Definitive" for one disease and
            "Limited" for another is reported as "Definitive" here. Use
            :meth:`get_disease_validity` for the full per-disease detail.
            Returns an empty dict if the bundled file is missing.
        """
        by_gene = self._ensure_validity_loaded()
        return {gene: self._best_entry(entries).classification for gene, entries in by_gene.items()}

    def get_disease_validity(self, gene_symbol: str) -> ClinGenValidity | None:
        """Return the single best-classified gene-disease curation for a gene.

        Args:
            gene_symbol: HGNC gene symbol, case-insensitive.

        Returns:
            The :class:`ClinGenValidity` entry with the strongest
            classification for this gene, or None if the gene has no
            ClinGen gene-disease-validity curation at all.
        """
        by_gene = self._ensure_validity_loaded()
        entries = by_gene.get(gene_symbol.upper())
        if not entries:
            return None
        return self._best_entry(entries)

    def get_haploinsufficiency_score(self, gene_symbol: str) -> int:
        """Return a gene's dosage-map haploinsufficiency score, normalised to 0-3.

        Args:
            gene_symbol: HGNC gene symbol, case-insensitive.

        Returns:
            0-3. Genes not present in the dosage map, and the two
            non-numeric special codes (30 = "autosomal recessive
            phenotype", 40 = "dosage sensitivity unlikely"), all
            normalise to 0 — none of these represent evidence *for*
            haploinsufficiency. Use :meth:`get_dosage_sensitivity` for
            the un-normalised code and its description.
        """
        dosage = self.get_dosage_sensitivity(gene_symbol)
        return dosage.haploinsufficiency_score if dosage is not None else 0

    def get_dosage_sensitivity(self, gene_symbol: str) -> DosageSensitivity | None:
        """Return the full dosage-map row for a gene, if curated.

        Args:
            gene_symbol: HGNC gene symbol, case-insensitive.

        Returns:
            A :class:`DosageSensitivity`, or None if the gene is absent
            from the bundled dosage map.
        """
        by_gene = self._ensure_dosage_loaded()
        return by_gene.get(gene_symbol.upper())

    @staticmethod
    def _best_entry(entries: list[ClinGenValidity]) -> ClinGenValidity:
        return min(
            entries,
            key=lambda e: _CLASSIFICATION_PRECEDENCE.get(
                e.classification, len(_CLASSIFICATION_PRECEDENCE)
            ),
        )

    def _ensure_validity_loaded(self) -> dict[str, list[ClinGenValidity]]:
        if self._validity_by_gene is not None:
            return self._validity_by_gene

        by_gene: dict[str, list[ClinGenValidity]] = {}
        if not self._validity_path.exists():
            logger.warning(
                "ClinGen gene-disease validity file not found; PVS1 will fall back to "
                "bundled gene list / gnomAD constraint",
                extra={"path": str(self._validity_path)},
            )
            self._validity_by_gene = by_gene
            return by_gene

        with self._validity_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle))

        header_index = next(
            (i for i, row in enumerate(rows) if row and row[0] == "GENE SYMBOL"), None
        )
        if header_index is None:
            logger.error(
                "ClinGen gene-disease validity file has an unrecognised format; "
                "expected a 'GENE SYMBOL' header row",
                extra={"path": str(self._validity_path)},
            )
            self._validity_by_gene = by_gene
            return by_gene

        for row in rows[header_index + 1 :]:
            if len(row) < 10 or row[0] in ("", "GENE SYMBOL") or set(row[0]) == {"+"}:
                continue
            gene_symbol, _hgnc_id, disease_label, mondo_id, moi, _sop, classification = row[:7]
            gcep = row[9] if len(row) > 9 else ""
            classification_date = _parse_clingen_datetime(row[8] if len(row) > 8 else None)
            entry = ClinGenValidity(
                gene_symbol=gene_symbol,
                disease_label=disease_label,
                mondo_id=mondo_id or None,
                mode_of_inheritance=moi or None,
                classification=classification,
                classification_date=classification_date,
                gcep=gcep,
            )
            by_gene.setdefault(gene_symbol.upper(), []).append(entry)

        logger.info(
            "loaded ClinGen gene-disease validity curations",
            extra={"path": str(self._validity_path), "gene_count": len(by_gene)},
        )
        self._validity_by_gene = by_gene
        return by_gene

    def _ensure_dosage_loaded(self) -> dict[str, DosageSensitivity]:
        if self._dosage_by_gene is not None:
            return self._dosage_by_gene

        by_gene: dict[str, DosageSensitivity] = {}
        if not self._dosage_path.exists():
            logger.warning(
                "ClinGen dosage sensitivity file not found; PVS1 haploinsufficiency "
                "scoring will treat all genes as score 0",
                extra={"path": str(self._dosage_path)},
            )
            self._dosage_by_gene = by_gene
            return by_gene

        with self._dosage_path.open(newline="", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 13:
                    continue
                gene_symbol = fields[0]
                raw_score = _parse_int(fields[4])
                description = fields[5]
                triplo_score = _parse_int(fields[12]) or 0
                if raw_score is None:
                    continue
                normalised_score = 0 if raw_score in _HI_SPECIAL_CODES else raw_score
                if raw_score in _HI_SPECIAL_CODES:
                    logger.debug(
                        "gene has a non-numeric ClinGen dosage code, normalised to 0",
                        extra={
                            "gene_symbol": gene_symbol,
                            "raw_code": raw_score,
                            "meaning": _HI_SPECIAL_CODES[raw_score],
                        },
                    )
                by_gene[gene_symbol.upper()] = DosageSensitivity(
                    gene_symbol=gene_symbol,
                    haploinsufficiency_score=normalised_score,
                    haploinsufficiency_raw_code=raw_score,
                    haploinsufficiency_description=description,
                    triplosensitivity_score=triplo_score,
                )

        logger.info(
            "loaded ClinGen dosage sensitivity map",
            extra={"path": str(self._dosage_path), "gene_count": len(by_gene)},
        )
        self._dosage_by_gene = by_gene
        return by_gene


def _parse_int(value: str) -> int | None:
    try:
        return int(value.strip())
    except (ValueError, AttributeError):
        return None


def _parse_clingen_datetime(value: str | None) -> date | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).date()
    except ValueError:
        return None
