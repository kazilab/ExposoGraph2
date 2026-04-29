"""Regression tests for the bundled reference graph builders."""

from collections import Counter

from ExposoGraph import (
    build_reference_architecture_summary,
    build_reference_engine,
    build_reference_graph,
)


def test_reference_graph_matches_current_bundled_counts():
    kg = build_reference_graph()
    node_types = Counter(node.type.value for node in kg.nodes)
    edge_types = Counter(edge.type.value for edge in kg.edges)

    assert len(kg.nodes) == 214
    assert len(kg.edges) == 321
    assert node_types == {
        "Carcinogen": 56,
        "Enzyme": 60,
        "Metabolite": 59,
        "DNA_Adduct": 27,
        "Pathway": 12,
    }
    assert edge_types == {
        "ACTIVATES": 104,
        "DETOXIFIES": 47,
        "FORMS_ADDUCT": 37,
        "REPAIRS": 33,
        "TRANSPORTS": 9,
        "PATHWAY": 88,
        "INDUCES": 2,
        "INHIBITS": 1,
    }


def test_reference_engine_loads_and_validates():
    engine = build_reference_engine()

    assert engine.node_count == 214
    assert engine.edge_count == 321
    assert engine.validate() == []


def test_reference_architecture_summary_matches_reference_graph():
    summary = build_reference_architecture_summary()

    assert summary.node_count == 214
    assert summary.edge_count == 321
    assert summary.node_type_counts == {
        "Carcinogen": 56,
        "Enzyme": 60,
        "Metabolite": 59,
        "DNA_Adduct": 27,
        "Pathway": 12,
    }
    assert summary.edge_type_counts == {
        "ACTIVATES": 104,
        "DETOXIFIES": 47,
        "FORMS_ADDUCT": 37,
        "REPAIRS": 33,
        "TRANSPORTS": 9,
        "PATHWAY": 88,
        "INDUCES": 2,
        "INHIBITS": 1,
    }
