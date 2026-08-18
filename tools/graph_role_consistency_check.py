#!/usr/bin/env python
"""Reaction-role consistency check for ExposoGraph 2.0's graph-data.json.

Cross-checks each ACTIVATES/DETOXIFIES edge's free-text `evidence` field
against its own edge `type`, flagging edges whose evidence describes the
opposite biological effect (e.g. an ACTIVATES edge whose evidence text is
plainly a clearance/detoxification narrative with no bioactivation override
language). This does not change any data -- it only surfaces candidates for
manual curation review.

Run from the repository root:

    python tools/graph_role_consistency_check.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

GRAPH_DATA_PATH = Path(__file__).resolve().parent.parent / "ExposoGraph" / "map" / "graph-data.json"

# Evidence-text keywords that plainly describe detoxification/clearance.
_DETOX_KEYWORDS = re.compile(r"\b(detoxif\w*|clearance|cleared|excretion|excreted)\b", re.IGNORECASE)
# Evidence-text keywords that plainly describe bioactivation/reactive-metabolite
# formation. Presence of these overrides a detox-keyword match, since several
# real phase II bioactivation edges (e.g. GSTT1->DCVG) legitimately describe
# both a protective baseline role AND an inverted bioactivation role in the
# same evidence string.
_BIOACTIVATION_KEYWORDS = re.compile(
    r"\b(bioactivat\w*|reactive\s+\w*|inverted\s+role|genotoxic)\b", re.IGNORECASE
)

CHECKS: list[dict[str, str]] = []


def record_check(name: str, status: str, detail: str) -> None:
    """Record one consistency-check outcome."""
    normalized = status.upper().strip()
    if normalized not in {"PASS", "WARN", "FAIL"}:
        raise ValueError(f"Unsupported check status: {status!r}")
    CHECKS.append({"name": name, "status": normalized, "detail": detail})


def _edge_label(edge: dict) -> str:
    return f"{edge.get('source')}->{edge.get('target')} ({edge.get('type')}, carcinogen={edge.get('carcinogen')})"


def check_activates_edges_for_detox_language(edges: list[dict]) -> None:
    """Flag ACTIVATES edges whose evidence text reads as pure detoxification."""
    flagged = []
    for edge in edges:
        if edge.get("type") != "ACTIVATES":
            continue
        evidence = edge.get("evidence", "") or ""
        if _DETOX_KEYWORDS.search(evidence) and not _BIOACTIVATION_KEYWORDS.search(evidence):
            flagged.append((edge, evidence))

    if not flagged:
        record_check(
            "ACTIVATES edges vs. detox-language evidence",
            "PASS",
            "No ACTIVATES edges found with unqualified detoxification/clearance language in evidence.",
        )
        return

    detail_lines = [f"{_edge_label(edge)}: {evidence!r}" for edge, evidence in flagged]
    record_check(
        "ACTIVATES edges vs. detox-language evidence",
        "WARN",
        f"{len(flagged)} ACTIVATES edge(s) have evidence text describing clearance/detoxification "
        f"with no bioactivation override language -- candidates for manual review: "
        + "; ".join(detail_lines),
    )


def check_detoxifies_edges_for_bioactivation_language(edges: list[dict]) -> None:
    """Flag DETOXIFIES edges whose evidence text reads as bioactivation."""
    flagged = []
    for edge in edges:
        if edge.get("type") != "DETOXIFIES":
            continue
        evidence = edge.get("evidence", "") or ""
        if _BIOACTIVATION_KEYWORDS.search(evidence):
            flagged.append((edge, evidence))

    if not flagged:
        record_check(
            "DETOXIFIES edges vs. bioactivation-language evidence",
            "PASS",
            "No DETOXIFIES edges found with bioactivation language in evidence.",
        )
        return

    detail_lines = [f"{_edge_label(edge)}: {evidence!r}" for edge, evidence in flagged]
    record_check(
        "DETOXIFIES edges vs. bioactivation-language evidence",
        "WARN",
        f"{len(flagged)} DETOXIFIES edge(s) have evidence text describing bioactivation -- "
        f"candidates for manual review: " + "; ".join(detail_lines),
    )


def check_same_enzyme_carcinogen_type_conflicts(edges: list[dict]) -> None:
    """Flag (enzyme, carcinogen) pairs with both ACTIVATES and DETOXIFIES edges
    whose evidence text does not clearly distinguish sequential pathway steps.
    This is a coarse heuristic intended to surface candidates, not a definitive
    contradiction detector.
    """
    by_pair: dict[tuple[str, str], list[dict]] = {}
    for edge in edges:
        if edge.get("type") not in {"ACTIVATES", "DETOXIFIES"}:
            continue
        key = (edge.get("source"), edge.get("carcinogen"))
        by_pair.setdefault(key, []).append(edge)

    flagged = []
    for (enzyme, carcinogen), pair_edges in by_pair.items():
        types_present = {e.get("type") for e in pair_edges}
        if types_present == {"ACTIVATES", "DETOXIFIES"}:
            flagged.append((enzyme, carcinogen, pair_edges))

    if not flagged:
        record_check(
            "Same enzyme+carcinogen ACTIVATES/DETOXIFIES co-occurrence",
            "PASS",
            "No enzyme+carcinogen pairs carry both an ACTIVATES and a DETOXIFIES edge.",
        )
        return

    detail_lines = []
    for enzyme, carcinogen, pair_edges in flagged:
        edge_descriptions = "; ".join(_edge_label(e) for e in pair_edges)
        detail_lines.append(f"{enzyme}+{carcinogen}: {edge_descriptions}")
    record_check(
        "Same enzyme+carcinogen ACTIVATES/DETOXIFIES co-occurrence",
        "WARN",
        f"{len(flagged)} enzyme+carcinogen pair(s) carry both edge types -- verify each pair represents "
        f"distinct sequential pathway steps rather than a mislabeled duplicate: " + "; ".join(detail_lines),
    )


def main() -> int:
    data = json.loads(GRAPH_DATA_PATH.read_text(encoding="utf-8"))
    edges = data.get("edges", [])
    if not isinstance(edges, list) or not edges:
        record_check("Graph data load", "FAIL", f"No edges found in {GRAPH_DATA_PATH}")
    else:
        record_check("Graph data load", "PASS", f"Loaded {len(edges)} edges from {GRAPH_DATA_PATH.name}")
        check_activates_edges_for_detox_language(edges)
        check_detoxifies_edges_for_bioactivation_language(edges)
        check_same_enzyme_carcinogen_type_conflicts(edges)

    max_name = max(len(item["name"]) for item in CHECKS)
    print("ExposoGraph 2.0 reaction-role consistency check")
    print("=" * 55)
    for item in CHECKS:
        print(f"{item['status']:<5} {item['name']:<{max_name}}  {item['detail']}")

    failed = [item for item in CHECKS if item["status"] == "FAIL"]
    warned = [item for item in CHECKS if item["status"] == "WARN"]
    print("-" * 55)
    print(f"Summary: {len(CHECKS) - len(failed) - len(warned)} PASS, {len(warned)} WARN, {len(failed)} FAIL")
    if failed:
        print("Reaction-role consistency check: FAIL")
        return 1
    print("Reaction-role consistency check: PASS (WARN items are curation candidates, not hard failures)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
