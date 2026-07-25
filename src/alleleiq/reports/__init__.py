"""Evidence report generation: structured JSON and human-readable HTML output."""

from alleleiq.reports.html_report import render_html_report, write_html_report
from alleleiq.reports.json_report import SCHEMA_VERSION, render_json_report, write_json_report

__all__ = [
    "SCHEMA_VERSION",
    "render_html_report",
    "render_json_report",
    "write_html_report",
    "write_json_report",
]
