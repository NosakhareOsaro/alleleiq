"""Evidence report generation: structured JSON and human-readable HTML output."""

from alleleiq.reports.json_report import SCHEMA_VERSION, render_json_report, write_json_report

__all__ = ["SCHEMA_VERSION", "render_json_report", "write_json_report"]
