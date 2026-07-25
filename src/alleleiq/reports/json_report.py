"""Structured, machine-readable JSON evidence report.

The JSON report is the canonical machine-readable output — intended for
future LIMS integration — and carries a ``schema_version`` independent
of :data:`alleleiq.__version__.VERSION`, so downstream consumers can
detect breaking schema changes across AlleleIQ releases without tying
their parsing logic to the software version string.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from alleleiq.exceptions import ReportGenerationError
from alleleiq.models.evidence import CriterionEvaluation
from alleleiq.models.report import EvidenceReport

SCHEMA_VERSION = "1.0.0"


def _partition_criteria(
    evaluations: list[CriterionEvaluation],
) -> tuple[list[CriterionEvaluation], list[CriterionEvaluation], list[CriterionEvaluation]]:
    applied = [c for c in evaluations if c.met]
    not_met = [c for c in evaluations if not c.met and not c.requires_clinical_input]
    not_assessed = [c for c in evaluations if not c.met and c.requires_clinical_input]
    return applied, not_met, not_assessed


def report_to_dict(report: EvidenceReport) -> dict[str, Any]:
    """Convert an :class:`EvidenceReport` into the versioned JSON schema.

    Args:
        report: The complete in-memory evidence report.

    Returns:
        A JSON-serialisable dict matching schema version
        :data:`SCHEMA_VERSION`, reorganising ``EvidenceReport``'s flat
        fields into the ``population_analysis`` / ``evidence_evaluation``
        / ``bayesian_score`` / ``audit_trail`` sections documented in
        ``docs/architecture/data-models.md``.
    """
    dumped = report.model_dump(mode="json")
    applied, not_met, not_assessed = _partition_criteria(report.criterion_evaluations)

    return {
        "schema_version": SCHEMA_VERSION,
        "report": {
            "metadata": dumped["metadata"],
            "variant": dumped["variant"],
            "population_analysis": {
                "frequencies": dumped["population_frequencies"],
                "coverage_gaps": dumped["coverage_gaps"],
                "coverage_summary": dumped["coverage_summary"],
            },
            "evidence_evaluation": {
                "criteria_applied": [c.model_dump(mode="json") for c in applied],
                "criteria_not_met": [c.model_dump(mode="json") for c in not_met],
                "criteria_not_assessed": [c.model_dump(mode="json") for c in not_assessed],
            },
            "bayesian_score": {
                "evidence_points": dumped["bayesian_score"]["evidence_points"],
                "preliminary_classification": dumped["bayesian_score"][
                    "preliminary_classification"
                ],
                "classification_confidence": dumped["bayesian_score"]["classification_confidence"],
                "contributing_criteria": dumped["bayesian_score"]["criteria_applied"],
                "evidence_points_confidence_interval": dumped["bayesian_score"][
                    "evidence_points_confidence_interval"
                ],
                "conflict_flags": dumped["bayesian_score"]["conflict_flags"],
            },
            "audit_trail": {
                "database_queries": dumped["audit_trail"]["queries"],
                "database_versions": dumped["metadata"]["database_versions"],
            },
            "clinical_disclaimer": dumped["clinical_disclaimer"],
        },
    }


def render_json_report(report: EvidenceReport, *, indent: int = 2) -> str:
    """Render an :class:`EvidenceReport` as a formatted JSON string.

    Args:
        report: The complete in-memory evidence report.
        indent: JSON indentation width.

    Returns:
        A formatted JSON string.
    """
    return json.dumps(report_to_dict(report), indent=indent, sort_keys=False)


def write_json_report(report: EvidenceReport, output_path: str | Path) -> Path:
    """Render and write a JSON evidence report to disk.

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
    path.write_text(render_json_report(report))

    if path.stat().st_size == 0:
        raise ReportGenerationError(f"generated JSON report at {path} is unexpectedly empty")
    return path
