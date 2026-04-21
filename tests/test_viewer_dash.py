"""Tests for Dash Cytoscape viewer helpers."""

from collections.abc import Iterable

from ExposoGraph.cytoscape_adapter import build_cytoscape_bundle
from ExposoGraph.engine import GraphEngine
from ExposoGraph.models import Edge, EdgeType, MatchStatus, Node, NodeType, RecordOrigin
from ExposoGraph.viewer_dash import (
    _merge_bundle_positions,
    _toggle_filter_values,
    apply_viewer_filters,
    build_detail_payload,
    create_dash_viewer_app,
)


def _viewer_engine() -> GraphEngine:
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
            id="NNK",
            label="NNK",
            type=NodeType.CARCINOGEN,
            group="Nitrosamine",
            origin=RecordOrigin.SEEDED,
            match_status=MatchStatus.CANONICAL,
        )
    )
    engine.add_node(Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME))
    engine.add_node(Node(id="CYP2A6", label="CYP2A6", type=NodeType.ENZYME))
    engine.add_node(Node(id="BPDE", label="BPDE", type=NodeType.METABOLITE, detail="Reactive diol epoxide"))
    engine.add_node(Node(id="BPDE_dG", label="BPDE-dG", type=NodeType.DNA_ADDUCT))
    engine.add_edge(Edge(source="BaP", target="CYP1A1", type=EdgeType.INDUCES))
    engine.add_edge(Edge(source="CYP1A1", target="BPDE", type=EdgeType.ACTIVATES, carcinogen="BaP"))
    engine.add_edge(Edge(source="BPDE", target="BPDE_dG", type=EdgeType.FORMS_ADDUCT, carcinogen="BaP"))
    engine.add_edge(Edge(source="NNK", target="CYP2A6", type=EdgeType.INDUCES))
    return engine


def _component_ids(component: object) -> set[str]:
    ids: set[str] = set()
    stack = [component]
    while stack:
        current = stack.pop()
        if current is None:
            continue
        component_id = getattr(current, "id", None)
        if isinstance(component_id, str):
            ids.add(component_id)
        children = getattr(current, "children", None)
        if isinstance(children, Iterable) and not isinstance(children, (str, bytes)):
            stack.extend(children)
        elif children is not None:
            stack.append(children)
    return ids


def test_apply_viewer_filters_supports_group_and_search():
    engine = _viewer_engine()
    base_bundle = build_cytoscape_bundle(engine)
    state = apply_viewer_filters(base_bundle, carcinogen_group="PAH")

    assert state.visible_node_count == 4
    assert state.visible_edge_count == 3
    assert "Pathway" not in state.node_type_counts

    searched = apply_viewer_filters(base_bundle, search_query="bpde")
    node_classes = {
        element["data"]["id"]: set(str(element.get("classes", "")).split())
        for element in searched.elements
        if "source" not in element["data"]
    }

    assert "dimmed" not in node_classes["BPDE"]
    assert "dimmed" in node_classes["NNK"]


def test_apply_viewer_filters_marks_focus_neighborhood():
    base_bundle = build_cytoscape_bundle(_viewer_engine())
    state = apply_viewer_filters(base_bundle, focus_node_id="BPDE")

    node_classes = {
        element["data"]["id"]: set(str(element.get("classes", "")).split())
        for element in state.elements
        if "source" not in element["data"]
    }

    assert "selected" in node_classes["BPDE"]
    assert "connected" in node_classes["CYP1A1"]
    assert "connected" in node_classes["BPDE_dG"]
    assert "dimmed" in node_classes["NNK"]


def test_build_detail_payload_for_node_and_edge():
    bundle = build_cytoscape_bundle(_viewer_engine())

    node_payload = build_detail_payload(bundle, {"kind": "node", "id": "BPDE"})
    edge_state = apply_viewer_filters(bundle)
    edge_id = next(
        element["data"]["id"]
        for element in edge_state.elements
        if element["data"].get("source") == "BPDE"
    )
    edge_payload = build_detail_payload(bundle, {"kind": "edge", "id": edge_id})

    assert node_payload["title"] == "BPDE"
    assert any(section[0] == "Outgoing" for section in node_payload["connections"])
    assert edge_payload["summary"] == "FORMS_ADDUCT"
    assert ("Source", "BPDE") in edge_payload["fields"]


def test_toggle_filter_values_preserves_order_and_toggles():
    ordered = ["Carcinogen", "Enzyme", "Metabolite"]

    reduced = _toggle_filter_values(["Carcinogen", "Enzyme"], "Carcinogen", available_values=ordered)
    expanded = _toggle_filter_values(["Enzyme"], "Metabolite", available_values=ordered)

    assert reduced == ["Enzyme"]
    assert expanded == ["Enzyme", "Metabolite"]


def test_merge_bundle_positions_updates_node_positions():
    bundle = build_cytoscape_bundle(_viewer_engine())
    moved_elements = []
    for element in bundle.elements:
        copied = dict(element)
        if copied["data"].get("id") == "BPDE" and "source" not in copied["data"]:
            copied["position"] = {"x": 123.0, "y": 456.0}
        moved_elements.append(copied)

    merged = _merge_bundle_positions(bundle.to_dict(), moved_elements)

    assert merged is not None
    assert merged["positions"]["BPDE"] == {"x": 123.0, "y": 456.0}
    bpde_element = next(
        element for element in merged["elements"] if element["data"].get("id") == "BPDE" and "source" not in element["data"]
    )
    assert bpde_element["position"] == {"x": 123.0, "y": 456.0}


def test_create_dash_viewer_app_exposes_expected_component_ids():
    app = create_dash_viewer_app(_viewer_engine(), title="Advanced Viewer Test")
    ids = _component_ids(app.layout)

    assert "advanced-viewer-graph" in ids
    assert "viewer-hover-tooltip" in ids
    assert "viewer-search" in ids
    assert "viewer-node-types" in ids
    assert "viewer-edge-types" in ids
    assert "viewer-carcinogen-group" in ids


