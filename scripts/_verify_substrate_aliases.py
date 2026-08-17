"""One-off verification script: check whether any NodeType.SUBSTRATE node in
graph-data.json is actually an alias (abbreviation/canonical-name variant) of
an existing non-Substrate node, rather than a genuinely new entity.

Existing nodes carry a `canonical_label` field (populated when
`match_status == "alias"`, e.g. "TCE" -> canonical_label "Trichloroethylene")
that a naive id/label comparison misses. This script normalizes and compares
against id, label, AND canonical_label on both sides to catch that case.

Read-only -- does not modify graph-data.json.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

GRAPH_PATH = Path(__file__).resolve().parent.parent / "ExposoGraph" / "map" / "graph-data.json"


def normalize(text: str | None) -> str:
    """Aggressively normalize for fuzzy comparison: lowercase, strip all
    non-alphanumeric characters (so "N-OH-PhIP", "n_oh_phip", "N-OH PhIP"
    all collapse to the same key)."""
    if not text:
        return ""
    return re.sub(r"[^a-z0-9]", "", text.lower())


def main() -> None:
    data = json.loads(GRAPH_PATH.read_text(encoding="utf-8"))
    nodes = data["nodes"]
    existing = [n for n in nodes if n["type"] != "Substrate"]
    substrates = [n for n in nodes if n["type"] == "Substrate"]

    index: dict[str, list[tuple[str, str, str | None, str]]] = {}
    for n in existing:
        for field in ("id", "label", "canonical_label"):
            key = normalize(n.get(field))
            if key:
                index.setdefault(key, []).append(
                    (n["id"], n["label"], n.get("canonical_label"), n["type"])
                )

    hits = []
    for s in substrates:
        keys = {normalize(s["id"]), normalize(s["label"])}
        found = [match for key in keys for match in index.get(key, [])]
        if found:
            hits.append((s["id"], s["label"], found))

    print(f"Substrate nodes checked: {len(substrates)}")
    print(f"Aliasing matches found (id/label/canonical_label): {len(hits)}")
    for sid, slabel, found in hits:
        print(f"  {sid} ({slabel}) -> {found}")
    if not hits:
        print("  (none)")


if __name__ == "__main__":
    main()
