from __future__ import annotations

import json

from ExposoGraph.graph_filtering import (
    GraphFilterCriteria,
    filter_graph_by_criteria,
    graph_filter_to_json_safe,
    heavy_metal_node_ids,
)
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, Node, NodeType


def _graph() -> KnowledgeGraph:
    return KnowledgeGraph(
        nodes=[
            Node(
                id="chromium_vi",
                label="Chromium(VI)",
                type=NodeType.CARCINOGEN,
                group="Heavy Metals",
                tissue_weights={"Lung": 0.8, "Liver": 0.2},
            ),
            Node(
                id="benzo_a_pyrene",
                label="Benzo[a]pyrene",
                type=NodeType.CARCINOGEN,
                group="PAH",
                tissue_weights={"Lung": 0.7},
            ),
            Node(
                id="GSTM1",
                label="GSTM1",
                type=NodeType.ENZYME,
                group="Heavy Metals",
                tissue_weights={"Lung": 0.75, "Liver": 0.1},
            ),
            Node(
                id="CYP1A1",
                label="CYP1A1",
                type=NodeType.ENZYME,
                group="PAH",
                tissue_weights={"Lung": 0.7},
            ),
            Node(id="lung", label="Lung", type=NodeType.TISSUE),
        ],
        edges=[
            Edge(
                source="chromium_vi",
                target="GSTM1",
                type=EdgeType.DETOXIFIES,
                carcinogen="chromium_vi",
                tissue="Lung",
            ),
            Edge(
                source="benzo_a_pyrene",
                target="CYP1A1",
                type=EdgeType.ACTIVATES,
                carcinogen="benzo_a_pyrene",
                tissue="Lung",
            ),
            Edge(source="GSTM1", target="lung", type=EdgeType.EXPRESSED_IN, tissue="Lung"),
        ],
    )


def test_issue14_filter_intersection_logic() -> None:
    filtered = filter_graph_by_criteria(
        _graph(),
        GraphFilterCriteria.from_values(
            carcinogen_classes=["Heavy Metals"],
            tissues=["Lung"],
            min_tissue_weight=0.7,
        ),
    )

    assert {node.id for node in filtered.nodes} == {"chromium_vi", "GSTM1"}
    assert [(edge.source, edge.target) for edge in filtered.edges] == [
        ("chromium_vi", "GSTM1")
    ]


def test_issue14_filter_by_carcinogen_class() -> None:
    filtered = filter_graph_by_criteria(_graph(), carcinogen_classes=["PAH"])

    assert {node.id for node in filtered.nodes} == {"benzo_a_pyrene", "CYP1A1"}
    assert {edge.carcinogen for edge in filtered.edges} == {"benzo_a_pyrene"}


def test_issue14_filter_by_tissue_weight_metadata() -> None:
    filtered = filter_graph_by_criteria(
        _graph(),
        tissues=["Lung"],
        min_tissue_weight=0.75,
    )

    assert {node.id for node in filtered.nodes} == {"chromium_vi", "GSTM1"}
    assert filtered.edges


def test_issue14_heavy_metal_nodes_edges_retained_when_present() -> None:
    graph = _graph()
    filtered = filter_graph_by_criteria(graph, carcinogen_classes=["metallograph"])

    assert heavy_metal_node_ids(graph) == {"chromium_vi"}
    assert "chromium_vi" in {node.id for node in filtered.nodes}
    assert ("chromium_vi", "GSTM1") in {
        (edge.source, edge.target) for edge in filtered.edges
    }


def test_issue14_graph_export_filtering_json_safe() -> None:
    filtered = filter_graph_by_criteria(
        _graph(),
        carcinogen_classes=["Heavy Metals"],
        tissues=["Lung"],
    )
    payload = graph_filter_to_json_safe(filtered)

    assert set(payload) == {"nodes", "edges"}
    assert payload["nodes"]
    assert payload["edges"]
    json.dumps(payload, sort_keys=True, allow_nan=False)

