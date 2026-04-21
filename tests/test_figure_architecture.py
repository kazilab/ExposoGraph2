"""Tests for the architecture figure helpers."""

from pathlib import Path

import matplotlib

from ExposoGraph.engine import GraphEngine
from ExposoGraph.figure_architecture import (
    build_architecture_figure_data,
    paper_architecture_overrides,
    save_architecture_figure,
)
from ExposoGraph.models import Edge, EdgeType, Node, NodeType

matplotlib.use("Agg")


def _figure_engine() -> GraphEngine:
    engine = GraphEngine()
    engine.add_node(Node(id="BaP", label="Benzo[a]pyrene", type=NodeType.CARCINOGEN, group="PAH"))
    engine.add_node(Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME, phase="I"))
    engine.add_node(Node(id="BPDE", label="BPDE", type=NodeType.METABOLITE))
    engine.add_node(Node(id="BPDE_dG", label="BPDE-dG", type=NodeType.DNA_ADDUCT))
    engine.add_node(Node(id="hsa05204", label="Chemical carcinogenesis", type=NodeType.PATHWAY))
    engine.add_edge(Edge(source="CYP1A1", target="BPDE", type=EdgeType.ACTIVATES))
    engine.add_edge(Edge(source="BPDE", target="BPDE_dG", type=EdgeType.FORMS_ADDUCT))
    engine.add_edge(Edge(source="CYP1A1", target="hsa05204", type=EdgeType.PATHWAY))
    return engine


def test_build_architecture_figure_data_uses_graph_counts():
    data = build_architecture_figure_data(_figure_engine())

    assert data.layer_counts["Carcinogens"] == 1
    assert data.layer_counts["Enzymes"] == 1
    assert data.layer_counts["Metabolites"] == 1
    assert data.layer_counts["DNA Adducts"] == 1
    assert data.layer_counts["KEGG Pathways"] == 1
    assert data.summary_lines[0] == "5 nodes  ·  3 edges"
    assert any(item.title == "PAH" for item in data.carcinogen_classes)
    assert any(
        item.title == "Phase I -- Activation" and item.count == 15
        for item in data.enzyme_categories
    )


def test_build_architecture_figure_data_supports_legacy_overrides():
    overrides = paper_architecture_overrides()
    data = build_architecture_figure_data(
        _figure_engine(),
        title=overrides["title"],
        subtitle=overrides["subtitle"],
        layer_count_overrides=overrides["layer_count_overrides"],
        enzyme_category_count_overrides=overrides["enzyme_category_count_overrides"],
        edge_count_overrides=overrides["edge_count_overrides"],
        summary_overrides=overrides["summary_overrides"],
    )

    assert data.title.startswith("CarcinoGenomic Knowledge Graph")
    assert data.layer_counts["Carcinogens"] == 20
    assert data.summary_lines[0] == "133 nodes  ·  166 edges"
    assert any(item.label == "ACTIVATES" and item.count == 57 for item in data.edge_legend)


def test_save_architecture_figure_writes_file(tmp_path: Path):
    data = build_architecture_figure_data(_figure_engine())
    output_path = save_architecture_figure(data, tmp_path / "architecture.png")

    assert output_path.exists()
    assert output_path.stat().st_size > 0
