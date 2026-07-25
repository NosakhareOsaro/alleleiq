#!/usr/bin/env python3
"""CLI entry point for AlleleIQ's known-variant validation.

Thin wrapper around :mod:`alleleiq.validation`, which holds the
gold-standard variant set, evaluation logic, and outcome model — see
that module's docstring for the full explanation of what this validates
and how to read its output. This script exists so the check can be run
directly (``uv run python scripts/validate_known_variants.py``) or
wired into CI without going through the ``alleleiq`` CLI.

Exit status is non-zero only if any variant is DANGEROUS (an
opposite-direction classification versus ClinVar), so this can be used
as a CI regression gate without failing on expected conservative
behaviour.
"""

from __future__ import annotations

import asyncio
import sys

from alleleiq.logging import configure_logging
from alleleiq.validation import Outcome, render_results, run_validation


def main() -> int:
    configure_logging("WARNING")
    results = asyncio.run(run_validation())
    render_results(results)
    return 0 if all(r.outcome != Outcome.DANGEROUS for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
