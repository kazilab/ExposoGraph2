"""Tests for Cytoscape adapter helpers."""

from pathlib import Path

from ExposoGraph.config import GraphVisibility
from ExposoGraph.cytoscape_adapter import (
    ViewerLayoutMode,
    build_cytoscape_bundle,
    compute_viewer_positions,
    load_cytoscape_bundle,
    load_viewer_positions,
    normalize_viewer_layout_mode,
    write_cytoscape_bundle,
    write_viewer_positions,
)
from ExposoGraph.engine import GraphEngine
from ExposoGraph.models import Edge, EdgeType, MatchStatus, Node, NodeType, RecordOrigin


def _sample_engine() -> GraphEngine:
    engine = GraphEngine()
    engine.add_node(
        Node(
            id="BaP",
            label="Benzo[a]pyrene",
            type=NodeType.CARCINOGEN,
            group="PAH",
            origin=RecordOrigin.SEEDED,
            match_status=MatchStatus.CANONICAL,
        )
    )
    engine.add_node(
        Node(
            id="CYP1A1",
            label="CYP1A1",
            type=NodeType.ENZYME,
            origin=RecordOrigin.LLM,
            match_status=MatchStatus.CANONICAL,
        )
    )
    engine.add_node(
        Node(
            id="BPDE",
            label="BPDE",
            type=NodeType.METABOLITE,
            origin=RecordOrigin.LLM,
            match_status=MatchStatus.UNMATCHED,
        )
    )
    engine.add_edge(
        Edge(
            source="BaP",
            target="CYP1A1",
            type=EdgeType.INDUCES,
            match_status=MatchStatus.CANONICAL,
        )
    )
    engine.add_edge(
        Edge(
            source="CYP1A1",
            target="BPDE",
            type=EdgeType.ACTIVATES,
            carcinogen="BaP",
            match_status=MatchStatus.UNMATCHED,
        )
    )
    return engine


def test_normalize_layout_mode_aliases():
    assert normalize_viewer_layout_mode("force") == ViewerLayoutMode.COSE
    assert normalize_viewer_layout_mode("saved") == ViewerLayoutMode.PRESET


def test_compute_positions_returns_all_nodes():
    engine = _sample_engine()
    positions = compute_viewer_positions(engine)

    assert set(positions) == {"BaP", "CYP1A1", "BPDE"}
    assert all("x" in coords and "y" in coords for coords in positions.values())


def test_build_bundle_contains_metadata_positions_and_elements():
    engine = _sample_engine()
    bundle = build_cytoscape_bundle(engine, layout_mode=ViewerLayoutMode.PRESET)

    assert bundle.metadata["node_count"] == 3
    assert bundle.metadata["edge_count"] == 2
    assert bundle.positions
    assert bundle.layout["name"] == "preset"
    assert len(bundle.elements) == 5
    assert any(
        element["data"].get("shape") == "diamond"
        for element in bundle.elements
        if "source" not in element["data"]
    )
    assert all(
        "position" in element
        for element in bundle.elements
        if "source" not in element["data"]
    )


def test_bundle_and_position_exports_roundtrip(tmp_path: Path):
    engine = _sample_engine()
    positions_path = write_viewer_positions(engine, tmp_path / "layout.json")
    bundle_path = write_cytoscape_bundle(
        engine,
        tmp_path / "bundle.json",
        visibility=GraphVisibility.ALL,
        layout_mode=ViewerLayoutMode.CIRCLE,
    )

    loaded_positions = load_viewer_positions(positions_path)
    loaded_bundle = load_cytoscape_bundle(bundle_path)

    assert set(loaded_positions) == {"BaP", "CYP1A1", "BPDE"}
    assert loaded_bundle.metadata["layout_mode"] == "circle"
    assert loaded_bundle.metadata["node_count"] == 3
    assert len(loaded_bundle.elements) == 5
