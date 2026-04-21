#!/usr/bin/env python3
"""Build a Benzo[a]pyrene metabolism knowledge graph — no API key required.

This script demonstrates the ExposoGraph library API by constructing
a complete BaP metabolism graph from a pre-built JSON file, running
graph analysis, and exporting to multiple formats.

Usage::

    pip install ExposoGraph
    python examples/build_bap_graph.py
"""

from __future__ import annotations

import json
from pathlib import Path

from ExposoGraph import (
    GraphEngine,
    KnowledgeGraph,
    centrality,
    metabolism_chain,
    shortest_path,
    variant_impact_score,
)
from ExposoGraph.exporter import to_interactive_html, to_json

EXAMPLE_DIR = Path(__file__).parent
OUTPUT_DIR = EXAMPLE_DIR / "output"


def main() -> None:
    # ── Load the pre-built BaP graph ──────────────────────────────────
    with open(EXAMPLE_DIR / "bap_graph.json") as f:
        data = json.load(f)

    kg = KnowledgeGraph(**data)
    engine = GraphEngine()
    engine.load(kg)

    print(f"Loaded graph: {engine.node_count} nodes, {engine.edge_count} edges\n")

    # ── Graph analysis ────────────────────────────────────────────────
    # Shortest path: BaP activation to DNA adduct
    path = shortest_path(engine, "CYP1A1", "BPDE_dG")
    print(f"Shortest path CYP1A1 → BPDE-dG: {' → '.join(path or [])}")

    # Centrality: which nodes are most connected?
    scores = centrality(engine, method="degree")
    top5 = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:5]
    print("\nTop-5 nodes by degree centrality:")
    for node_id, score in top5:
        print(f"  {node_id:12s}  {score:.3f}")

    # Metabolism chain for BaP
    chain = metabolism_chain(engine, "BaP")
    print(f"\nBaP metabolism chain: {len(chain.node_ids)} nodes, {len(chain.edges)} edges")
    print(f"  Activation edges:    {len(chain.activation_edges)}")
    print(f"  Detoxification edges: {len(chain.detox_edges)}")
    print(f"  Adduct edges:        {len(chain.adduct_edges)}")
    print(f"  Repair edges:        {len(chain.repair_edges)}")

    # Variant impact for CYP1A1
    impact = variant_impact_score(engine, "CYP1A1")
    if impact:
        print("\nVariant impact score for CYP1A1:")
        print(f"  Activity score:        {impact.activity_score}")
        print(f"  Downstream adducts:    {impact.downstream_adduct_count}")
        print(f"  Impact score:          {impact.score:.2f}")

    # ── Export ────────────────────────────────────────────────────────
    OUTPUT_DIR.mkdir(exist_ok=True)

    # JSON
    json_path = to_json(engine, OUTPUT_DIR / "bap_graph_export.json")
    print(f"\nExported JSON → {json_path}")

    # Interactive HTML
    html_path = to_interactive_html(engine, OUTPUT_DIR / "bap_graph.html")
    print(f"Exported HTML → {html_path}")
    print("\nDone! Open the HTML file in a browser to explore the graph interactively.")

    print("\nDone!")


if __name__ == "__main__":
    main()
