#!/usr/bin/env python3
"""Exploratory analysis: how AlleleIQ's evidence differs by gnomAD ancestry group.

Question: for rare missense variants in hereditary cancer genes, do the
computational criteria (PP3/BP4, from REVEL and SpliceAI) and the final
preliminary classification differ depending on which ancestry group the
variant is found in?

The design is fixed in this file before any results are seen:

- Genes: ten hereditary cancer predisposition genes (``GENES``).
- Variants: gnomAD v4 missense variants whose allele frequency within an
  ancestry group is between 0.1% and 1% (the band used by Pathak et al.,
  "Pervasive ancestry bias in variant effect predictors"), counting exome
  and genome data together, in a group with at least ``MIN_GROUP_AN``
  sequenced alleles. A variant can belong to more than one group.
- Groups: afr, amr, asj, eas, fin, nfe, sas. gnomAD's Middle Eastern
  group is excluded because AlleleIQ does not yet model it, and
  "remaining" is excluded because it is not an ancestry group.
- Every selected variant is run once through the full AlleleIQ pipeline,
  then results are tallied per group.

This is a descriptive, exploratory analysis of a small sample, not a
validation study. Results are reported with Wilson 95% intervals.

Usage::

    uv run python scripts/analyse_ancestry_disparity.py select
    uv run python scripts/analyse_ancestry_disparity.py run
    uv run python scripts/analyse_ancestry_disparity.py summarise
"""

from __future__ import annotations

import asyncio
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import httpx

from alleleiq.logging import configure_logging
from alleleiq.models.evidence import ACMGCode, PreliminaryClassification
from alleleiq.models.variant import GenomicCoordinates, Variant, VariantType
from alleleiq.pipeline import build_evidence_report

GENES = ("BRCA1", "BRCA2", "PALB2", "TP53", "PTEN", "MLH1", "MSH2", "MSH6", "ATM", "CHEK2")
GROUPS = ("afr", "amr", "asj", "eas", "fin", "nfe", "sas")
AF_LOW, AF_HIGH = 0.001, 0.01
MIN_GROUP_AN = 2000
# One at a time: in a first run with three in parallel, gnomAD silently failed for
# about a third of variants (the pipeline records source failures rather than raising).
CONCURRENCY = 1

OUT_DIR = Path("reports/ancestry-analysis")
SELECTED = OUT_DIR / "selected_variants.json"
RESULTS = OUT_DIR / "results.jsonl"
SUMMARY = OUT_DIR / "summary.md"

_GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
_GENE_QUERY = """
query($gene: String!) {
  gene(gene_symbol: $gene, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id consequence hgvsp
      exome { populations { id ac an } }
      genome { populations { id ac an } }
    }
  }
}
"""
_REVEL_PATTERN = re.compile(r"REVEL=([0-9.]+)")
_BENIGN_SIDE = {PreliminaryClassification.LIKELY_BENIGN, PreliminaryClassification.BENIGN}
_PATHOGENIC_SIDE = {
    PreliminaryClassification.LIKELY_PATHOGENIC,
    PreliminaryClassification.PATHOGENIC,
}


# ── select ──────────────────────────────────────────────────────────────
def _group_frequencies(variant: dict[str, Any]) -> dict[str, tuple[int, int]]:
    """Sum exome and genome allele counts/numbers per ancestry group."""
    totals: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for source in ("exome", "genome"):
        data = variant.get(source) or {}
        for pop in data.get("populations") or []:
            if pop["id"] in GROUPS:
                totals[pop["id"]][0] += pop["ac"]
                totals[pop["id"]][1] += pop["an"]
    return {group: (ac, an) for group, (ac, an) in totals.items()}


def select() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    selected: list[dict[str, Any]] = []
    with httpx.Client(timeout=180) as client:
        for gene in GENES:
            response = client.post(
                _GNOMAD_URL, json={"query": _GENE_QUERY, "variables": {"gene": gene}}
            )
            response.raise_for_status()
            variants = response.json()["data"]["gene"]["variants"]
            kept = 0
            for variant in variants:
                if variant["consequence"] != "missense_variant":
                    continue
                freqs = _group_frequencies(variant)
                in_band = {
                    group: ac / an
                    for group, (ac, an) in freqs.items()
                    if an >= MIN_GROUP_AN and AF_LOW <= ac / an <= AF_HIGH
                }
                if not in_band:
                    continue
                selected.append(
                    {
                        "variant_id": variant["variant_id"],
                        "gene": gene,
                        "hgvsp": variant["hgvsp"],
                        "groups": sorted(in_band),
                        "group_af": in_band,
                        "all_group_af": {
                            g: (ac / an if an else None) for g, (ac, an) in freqs.items()
                        },
                    }
                )
                kept += 1
            print(f"{gene}: {len(variants)} variants, {kept} selected")
    SELECTED.write_text(json.dumps(selected, indent=1))
    print(f"\n{len(selected)} variants selected -> {SELECTED}")


# ── run ─────────────────────────────────────────────────────────────────
def _variant_from_id(variant_id: str) -> Variant:
    chrom, pos, ref, alt = variant_id.split("-")
    return Variant(
        coordinates=GenomicCoordinates(chrom=chrom, pos=int(pos), ref=ref, alt=alt),
        variant_type=VariantType.SNV if len(ref) == len(alt) == 1 else VariantType.INDEL,
    )


async def _analyse(entry: dict[str, Any], semaphore: asyncio.Semaphore) -> dict[str, Any]:
    async with semaphore:
        try:
            report = await build_evidence_report(_variant_from_id(entry["variant_id"]))
        except Exception as exc:  # recorded, not fatal: the run continues
            return {**entry, "error": f"{type(exc).__name__}: {exc}"}

    failed = sorted({q.database for q in report.audit_trail.queries if not q.success})
    if failed:
        # A failed source would silently change the evidence, so retry it instead.
        return {**entry, "error": f"evidence source failed: {', '.join(failed)}"}

    criteria = {e.code.value: e for e in report.criterion_evaluations}
    revel = None
    for code in (ACMGCode.PP3, ACMGCode.BP4):
        evaluation = criteria.get(code.value)
        if evaluation and (match := _REVEL_PATTERN.search(evaluation.evidence_basis)):
            revel = float(match.group(1))
    return {
        **entry,
        "consequence": report.variant.consequence,
        "revel": revel,
        "met": sorted(code for code, e in criteria.items() if e.met),
        "met_strength": {code: e.strength.value for code, e in criteria.items() if e.met},
        "points": report.bayesian_score.evidence_points,
        "classification": report.bayesian_score.preliminary_classification.value,
        "coverage_gaps": sorted(g.value for g in report.coverage_gaps),
        "versions": report.metadata.database_versions,
    }


async def _run_all() -> None:
    selected = json.loads(SELECTED.read_text())
    done: set[str] = set()
    if RESULTS.exists():
        for line in RESULTS.read_text().splitlines():
            record = json.loads(line)
            if "error" not in record:
                done.add(record["variant_id"])
    pending = [e for e in selected if e["variant_id"] not in done]
    print(f"{len(done)} already done, {len(pending)} to run")

    semaphore = asyncio.Semaphore(CONCURRENCY)
    tasks = [asyncio.create_task(_analyse(e, semaphore)) for e in pending]
    with RESULTS.open("a") as handle:
        for count, task in enumerate(asyncio.as_completed(tasks), start=1):
            record = await task
            handle.write(json.dumps(record) + "\n")
            handle.flush()
            status = record.get("error") or record["classification"]
            print(f"[{count}/{len(pending)}] {record['variant_id']} {record['gene']}: {status}")


def run() -> None:
    configure_logging("ERROR")
    asyncio.run(_run_all())


# ── summarise ───────────────────────────────────────────────────────────
def _wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    z = 1.96
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, centre - half), min(1.0, centre + half))


def _pct(k: int, n: int) -> str:
    if n == 0:
        return "-"
    lo, hi = _wilson(k, n)
    return f"{100 * k / n:.0f}% ({100 * lo:.0f}-{100 * hi:.0f})"


def summarise() -> None:
    latest: dict[str, dict[str, Any]] = {}
    for line in RESULTS.read_text().splitlines():
        record = json.loads(line)
        latest[record["variant_id"]] = record
    records = [r for r in latest.values() if "error" not in r]
    errors = [r for r in latest.values() if "error" in r]

    rows = []
    for group in GROUPS:
        members = [r for r in records if group in r["groups"]]
        n = len(members)
        with_revel = [r for r in members if r["revel"] is not None]
        rows.append(
            {
                "group": group,
                "n": n,
                "genes": len({r["gene"] for r in members}),
                "revel_median": (
                    sorted(r["revel"] for r in with_revel)[len(with_revel) // 2]
                    if with_revel
                    else None
                ),
                "pp3": sum("PP3" in r["met"] for r in members),
                "bp4": sum("BP4" in r["met"] for r in members),
                "bp1": sum("BP1" in r["met"] for r in members),
                "benign_side": sum(
                    PreliminaryClassification(r["classification"]) in _BENIGN_SIDE for r in members
                ),
                "vus": sum(
                    r["classification"] == PreliminaryClassification.UNCERTAIN_SIGNIFICANCE.value
                    for r in members
                ),
                "pathogenic_side": sum(
                    PreliminaryClassification(r["classification"]) in _PATHOGENIC_SIDE
                    for r in members
                ),
                "gap": sum(group in r["coverage_gaps"] for r in members),
            }
        )

    lines = [
        "# Ancestry analysis: AlleleIQ evidence by gnomAD group",
        "",
        f"Genes: {', '.join(GENES)}. Missense variants at {AF_LOW:.1%}-{AF_HIGH:.0%} "
        f"allele frequency within a group (group AN ≥ {MIN_GROUP_AN}).",
        f"Variants analysed: {len(records)} (errors: {len(errors)}). "
        "A variant can count in more than one group.",
        "Percentages with Wilson 95% intervals.",
        "",
        "| Group | n | Genes | Median REVEL | PP3 met | BP4 met | BP1 met | LB/B | VUS | LP/P | "
        "Coverage gap flagged |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        n = row["n"]
        median = f"{row['revel_median']:.2f}" if row["revel_median"] is not None else "-"
        lines.append(
            f"| {row['group']} | {n} | {row['genes']} | {median} | {_pct(row['pp3'], n)} | "
            f"{_pct(row['bp4'], n)} | {_pct(row['bp1'], n)} | {_pct(row['benign_side'], n)} | "
            f"{_pct(row['vus'], n)} | "
            f"{_pct(row['pathogenic_side'], n)} | {_pct(row['gap'], n)} |"
        )
    lines += [
        "",
        "BP1 is a gene-level criterion (curated loss-of-function genes such as BRCA1/BRCA2), so",
        "classification differences between groups partly reflect each group's gene mix. PP3, BP4",
        "and REVEL are the predictor-driven columns.",
        "",
        "## Gene mix per group",
        "",
        "| Group | " + " | ".join(GENES) + " |",
        "|---|" + "---|" * len(GENES),
    ]
    for group in GROUPS:
        counts = [sum(1 for r in records if group in r["groups"] and r["gene"] == g) for g in GENES]
        lines.append(f"| {group} | " + " | ".join(str(c) for c in counts) + " |")
    if records:
        lines += ["", f"Evidence-source versions: `{json.dumps(records[0]['versions'])}`"]
    SUMMARY.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    commands = {"select": select, "run": run, "summarise": summarise}
    if len(sys.argv) != 2 or sys.argv[1] not in commands:
        sys.exit(f"usage: {sys.argv[0]} {{{'|'.join(commands)}}}")
    commands[sys.argv[1]]()
