#!/usr/bin/env python3
"""Demonstrate graph modes and visibility-aware export without any API key.

This example builds a tiny mixed-quality graph in memory:

- one canonically groundable enzyme
- one alias-matched carcinogen
- one unmatched exploratory carcinogen

It then compares exploratory vs strict merge behavior, writes filtered JSON
exports, and saves a validated-only revision into a local SQLite database.
"""

from __future__ import annotations

from pathlib import Path

from ExposoGraph import (
    Edge,
    EdgeType,
    GraphEngine,
    GraphMode,
    GraphRepository,
    GraphVisibility,
    KnowledgeGraph,
    Node,
    NodeType,
    to_json,
)

EXAMPLE_DIR = Path(__file__).parent
OUTPUT_DIR = EXAMPLE_DIR / "output"


def build_mixed_graph() -> KnowledgeGraph:
    return KnowledgeGraph(
        nodes=[
            Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
            Node(id="n2", label="BaP", type=NodeType.CARCINOGEN),
            Node(id="n3", label="Mystery Chemical", type=NodeType.CARCINOGEN),
        ],
        edges=[
            Edge(source="n1", target="n2", type=EdgeType.ACTIVATES),
            Edge(source="n1", target="n3", type=EdgeType.ACTIVATES),
        ],
    )


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    raw_graph = build_mixed_graph()

    exploratory_engine = GraphEngine()
    exploratory_warnings = exploratory_engine.merge(
        raw_graph,
        mode=GraphMode.EXPLORATORY,
    )
    print("Exploratory merge")
    print(f"  nodes: {exploratory_engine.node_count}")
    print(f"  edges: {exploratory_engine.edge_count}")
    print(f"  warnings: {len(exploratory_warnings)}")

    strict_engine = GraphEngine()
    strict_warnings = strict_engine.merge(
        raw_graph,
        mode=GraphMode.STRICT,
    )
    print("\nStrict merge")
    print(f"  nodes: {strict_engine.node_count}")
    print(f"  edges: {strict_engine.edge_count}")
    for warning in strict_warnings:
        print(f"  warning: {warning}")

    validated_json = to_json(
        exploratory_engine,
        OUTPUT_DIR / "validated_only_demo.json",
        visibility=GraphVisibility.VALIDATED_ONLY,
    )
    exploratory_json = to_json(
        exploratory_engine,
        OUTPUT_DIR / "exploratory_only_demo.json",
        visibility=GraphVisibility.EXPLORATORY_ONLY,
    )
    print(f"\nValidated-only JSON → {validated_json}")
    print(f"Exploratory-only JSON → {exploratory_json}")

    with GraphRepository(OUTPUT_DIR / "mode_visibility_demo.sqlite3") as repo:
        saved = repo.save_engine(
            graph_key="mode_visibility_demo",
            graph_name="Mode Visibility Demo",
            engine=exploratory_engine,
            visibility=GraphVisibility.VALIDATED_ONLY,
            note="Validated slice from exploratory graph",
        )
        print(
            "\nSaved revision "
            f"r{saved.revision_number} with visibility={saved.visibility.value}"
        )


if __name__ == "__main__":
    main()
