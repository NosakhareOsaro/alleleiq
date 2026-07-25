"""Standalone, printable HTML evidence report.

Renders from the same :func:`~alleleiq.reports.json_report.report_to_dict`
structure as the JSON report, so the two formats are always derived
from one source of truth and cannot silently drift apart. The template
inlines all CSS and uses no external assets or JavaScript, so the
output file is fully standalone (openable offline, printable to PDF
directly from a browser).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, select_autoescape

from alleleiq.exceptions import ReportGenerationError
from alleleiq.models.population import AncestryGroup
from alleleiq.models.report import EvidenceReport
from alleleiq.reports.json_report import report_to_dict

_TEMPLATES_DIR = Path(__file__).parent / "templates"

_CLASSIFICATION_CSS_CLASS: dict[str, str] = {
    "PATHOGENIC": "classification-pathogenic",
    "LIKELY_PATHOGENIC": "classification-pathogenic",
    "UNCERTAIN_SIGNIFICANCE": "classification-vus",
    "LIKELY_BENIGN": "classification-benign",
    "BENIGN": "classification-benign",
}
_GAP_CSS_CLASS: dict[str, str] = {
    "WARNING": "gap-warning",
    "CRITICAL": "gap-critical",
    "NO_DATA": "gap-critical",
}


def _environment() -> Environment:
    return Environment(
        loader=FileSystemLoader(str(_TEMPLATES_DIR)),
        autoescape=select_autoescape(["html", "jinja2"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )


def _population_rows(report: EvidenceReport) -> list[dict[str, Any]]:
    """Build one row per ancestry group, merging frequency and gap data.

    Every :class:`AncestryGroup` appears exactly once, in canonical
    order, regardless of which are present in the underlying dicts.
    """
    rows = []
    for ancestry in AncestryGroup:
        frequency = report.population_frequencies.get(ancestry)
        gap = report.coverage_gaps.get(ancestry)
        rows.append(
            {
                "display_name": ancestry.display_name,
                "ancestry_code": ancestry.value.upper(),
                "allele_count": frequency.allele_count if frequency else None,
                "allele_number": frequency.allele_number if frequency else 0,
                "allele_frequency": frequency.allele_frequency if frequency else None,
                "confidence_interval_95": frequency.confidence_interval_95 if frequency else None,
                "gap_flag_level": gap.flag_level.value if gap else None,
                "gap_message": gap.clinical_guidance_message if gap else None,
                "gap_css_class": _GAP_CSS_CLASS.get(gap.flag_level.value, "") if gap else "",
            }
        )
    return rows


def render_html_report(report: EvidenceReport) -> str:
    """Render an :class:`EvidenceReport` as a standalone HTML document.

    Args:
        report: The complete in-memory evidence report.

    Returns:
        A complete, self-contained HTML document as a string.
    """
    payload = report_to_dict(report)["report"]
    classification = payload["bayesian_score"]["preliminary_classification"]

    template = _environment().get_template("report.html.jinja2")
    return template.render(
        report=payload,
        population_rows=_population_rows(report),
        classification_css_class=_CLASSIFICATION_CSS_CLASS.get(
            classification, "classification-vus"
        ),
    )


def write_html_report(report: EvidenceReport, output_path: str | Path) -> Path:
    """Render and write an HTML evidence report to disk.

    Args:
        report: The complete in-memory evidence report.
        output_path: Destination file path. Parent directories are
            created if needed.

    Returns:
        The path written to.

    Raises:
        ReportGenerationError: If the written file is unexpectedly empty.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_html_report(report))

    if path.stat().st_size == 0:
        raise ReportGenerationError(f"generated HTML report at {path} is unexpectedly empty")
    return path
