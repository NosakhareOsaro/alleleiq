"""Evidence report generation: structured JSON and human-readable HTML output."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from alleleiq.logging import get_logger
from alleleiq.models.report import EvidenceReport
from alleleiq.reports.html_report import render_html_report, write_html_report
from alleleiq.reports.json_report import SCHEMA_VERSION, render_json_report, write_json_report

logger = get_logger(__name__)

ReportFormat = Literal["json", "html"]


@dataclass(frozen=True)
class GeneratedReportPaths:
    """Filesystem paths for the reports produced by one orchestration run."""

    json_path: Path | None
    html_path: Path | None


class ReportOrchestrator:
    """Coordinates JSON and/or HTML evidence report generation for one variant.

    A thin coordination layer over :mod:`alleleiq.reports.json_report`
    and :mod:`alleleiq.reports.html_report`: it owns the output
    directory and file-naming convention, and logs generation
    statistics. Per-format output validation (non-empty files) is
    delegated to the individual ``write_*_report`` functions, which
    raise :class:`~alleleiq.exceptions.ReportGenerationError` on failure.
    """

    def __init__(self, output_dir: str | Path) -> None:
        """Initialise the orchestrator with a target output directory.

        Args:
            output_dir: Directory reports are written into. Created if
                it does not already exist.
        """
        self._output_dir = Path(output_dir)

    def generate(
        self,
        report: EvidenceReport,
        *,
        formats: tuple[ReportFormat, ...] = ("json", "html"),
    ) -> GeneratedReportPaths:
        """Generate the requested report formats for one variant.

        Args:
            report: The complete in-memory evidence report.
            formats: Which formats to generate. Defaults to both.

        Returns:
            The paths written, with the field for any format not
            requested left as None.
        """
        self._output_dir.mkdir(parents=True, exist_ok=True)
        base_name = report.variant.variant_id

        json_path = (
            write_json_report(report, self._output_dir / f"{base_name}.json")
            if "json" in formats
            else None
        )
        html_path = (
            write_html_report(report, self._output_dir / f"{base_name}.html")
            if "html" in formats
            else None
        )

        logger.info(
            "report generation complete",
            extra={
                "variant_id": report.variant.variant_id,
                "formats": list(formats),
                "output_dir": str(self._output_dir),
            },
        )
        return GeneratedReportPaths(json_path=json_path, html_path=html_path)


__all__ = [
    "SCHEMA_VERSION",
    "GeneratedReportPaths",
    "ReportFormat",
    "ReportOrchestrator",
    "render_html_report",
    "render_json_report",
    "write_html_report",
    "write_json_report",
]
