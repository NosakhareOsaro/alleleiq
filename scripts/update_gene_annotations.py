#!/usr/bin/env python3
"""Structured editor for AlleleIQ's curated gene missense-mechanism list.

``data/gene_annotations/missense_constraint.json`` (see
``data/gene_annotations/README.md``) drives the PP2/BP1 criteria in
:mod:`alleleiq.scoring.criteria.gene_level`. Unlike
``scripts/update_clingen.sh`` or ``scripts/download_revel.sh``, there is
no automated feed for "which genes have an established missense-vs-
loss-of-function disease mechanism" — that judgement requires reading
ClinGen gene curations or the primary literature. This script does not
make that judgement for you; it only provides a safe, structured way to
add or update one entry at a time (validating shape and the
pp2/bp1-mutual-exclusivity rule) and re-serialise the file in a
consistent, sorted, diff-friendly format, so curation work goes through
one reviewable mechanism instead of ad hoc hand-edits.

Examples:
    Add or update one gene::

        uv run python scripts/update_gene_annotations.py set \\
            --gene TTN --mechanism loss_of_function_only \\
            --bp1 --source "Dilated cardiomyopathy; A-band truncating variants ..."

    Just re-sort/re-format the file after a manual edit::

        uv run python scripts/update_gene_annotations.py normalise
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_DEFAULT_PATH = Path("data/gene_annotations/missense_constraint.json")


def _load(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    data: dict[str, dict[str, object]] = json.loads(path.read_text(encoding="utf-8"))
    return data


def _save(path: Path, data: dict[str, dict[str, object]]) -> None:
    ordered = {gene: data[gene] for gene in sorted(data)}
    path.write_text(json.dumps(ordered, indent=2) + "\n", encoding="utf-8")


def _cmd_set(args: argparse.Namespace) -> int:
    path = Path(args.path)
    data = _load(path)

    pp2_applicable = bool(args.pp2)
    bp1_applicable = bool(args.bp1)
    if pp2_applicable and bp1_applicable:
        print(
            "error: a gene cannot be both --pp2 and --bp1 applicable "
            "(that would claim missense is both a common and an implausible "
            "mechanism for the same gene)",
            file=sys.stderr,
        )
        return 1
    if not pp2_applicable and not bp1_applicable:
        print(
            "error: specify --pp2 or --bp1 (a gene with neither should simply "
            "be omitted from this file, not added with both false)",
            file=sys.stderr,
        )
        return 1

    gene_symbol = args.gene.upper()
    entry: dict[str, object] = {
        "pp2_applicable": pp2_applicable,
        "bp1_applicable": bp1_applicable,
        "mechanism": args.mechanism,
        "source": args.source,
    }
    # Domain exclusions are curated by hand-editing the file; keep them.
    existing_regions = data.get(gene_symbol, {}).get("bp1_excluded_regions")
    if existing_regions:
        entry["bp1_excluded_regions"] = existing_regions
    data[gene_symbol] = entry
    _save(path, data)
    print(f"Updated {gene_symbol} in {path} ({len(data)} genes total).")
    return 0


def _cmd_remove(args: argparse.Namespace) -> int:
    path = Path(args.path)
    data = _load(path)
    gene_symbol = args.gene.upper()
    if gene_symbol not in data:
        print(f"error: {gene_symbol} not present in {path}", file=sys.stderr)
        return 1
    del data[gene_symbol]
    _save(path, data)
    print(f"Removed {gene_symbol} from {path} ({len(data)} genes remaining).")
    return 0


def _cmd_normalise(args: argparse.Namespace) -> int:
    path = Path(args.path)
    data = _load(path)
    _save(path, data)
    print(f"Normalised {path} ({len(data)} genes).")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--path",
        default=str(_DEFAULT_PATH),
        help=f"Annotation file path (default: {_DEFAULT_PATH})",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    set_parser = subparsers.add_parser("set", help="Add or update one gene's annotation.")
    set_parser.add_argument("--gene", required=True, help="HGNC gene symbol.")
    set_parser.add_argument("--pp2", action="store_true", help="Mark pp2_applicable=true.")
    set_parser.add_argument("--bp1", action="store_true", help="Mark bp1_applicable=true.")
    set_parser.add_argument(
        "--mechanism", required=True, help="Short mechanism label, e.g. loss_of_function_only."
    )
    set_parser.add_argument(
        "--source", required=True, help="Citation-style rationale for this classification."
    )
    set_parser.set_defaults(func=_cmd_set)

    remove_parser = subparsers.add_parser("remove", help="Remove one gene's annotation.")
    remove_parser.add_argument("--gene", required=True, help="HGNC gene symbol.")
    remove_parser.set_defaults(func=_cmd_remove)

    normalise_parser = subparsers.add_parser(
        "normalise", help="Re-sort and re-format the file without changing content."
    )
    normalise_parser.set_defaults(func=_cmd_normalise)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
