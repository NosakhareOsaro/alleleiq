"""Minimal server-rendered web interface (upload a VCF, view/download evidence reports).

A thin presentational layer over :func:`alleleiq.pipeline.build_evidence_report`
-- the same entry point the CLI and JSON REST API use -- so this
interface can never compute a different result than either of those.
See :mod:`alleleiq.web.router` for the routes and
:mod:`alleleiq.web.templates` for the Jinja2 templates.
"""

from __future__ import annotations
