"""Integration test: build BaP graph → analyze → export → reimport → verify.

This test exercises the full pipeline end-to-end without any mocks,
using the pre-built example graph JSON.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ExposoGraph import (
    GraphEngine,
    KnowledgeGraph,
    all_shortest_paths,
    centrality,
    metabolism_chain,
    pathway_subgraph,
    shortest_path,
    variant_impact_score,
)
from ExposoGraph.db_clients.iarc import IARCClassifier, IARCGroup
from ExposoGraph.exporter import (
    parse_graph_data_text,
    to_gexf,
    to_graph_data_js,
    to_interactive_html_string,
    to_json,
)
from ExposoGraph.models import EdgeType, NodeType

EXAMPLE_JSON = Path(__file__).parent.parent / "examples" / "bap_graph.json"


@pytest.fixture(scope="module")
def bap_kg() -> KnowledgeGraph:
    """Load the pre-built BaP example graph."""
    with open(EXAMPLE_JSON) as f:
        data = json.load(f)
    return KnowledgeGraph(**data)


@pytest.fixture(scope="module")
def engine(bap_kg: KnowledgeGraph) -> GraphEngine:
    eng = GraphEngine()
    eng.load(bap_kg)
    return eng


# ── Load & validate ──────────────────────────────────────────────────────


class TestLoadAndValidate:
    def test_node_count(self, engine: GraphEngine):
        assert engine.node_count == 20

    def test_edge_count(self, engine: GraphEngine):
        assert engine.edge_count == 20

    def test_all_node_types_present(self, engine: GraphEngine):
        kg = engine.to_knowledge_graph()
        node_types = {n.type for n in kg.nodes}
        assert NodeType.CARCINOGEN in node_types
        assert NodeType.ENZYME in node_types
        assert NodeType.METABOLITE in node_types
        assert NodeType.DNA_ADDUCT in node_types
        assert NodeType.PATHWAY in node_types
        assert NodeType.TISSUE in node_types

    def test_all_edge_types_present(self, engine: GraphEngine):
        kg = engine.to_knowledge_graph()
        edge_types = {e.type for e in kg.edges}
        assert EdgeType.ACTIVATES in edge_types
        assert EdgeType.DETOXIFIES in edge_types
        assert EdgeType.TRANSPORTS in edge_types
        assert EdgeType.FORMS_ADDUCT in edge_types
        assert EdgeType.REPAIRS in edge_types
        assert EdgeType.PATHWAY in edge_types
        assert EdgeType.EXPRESSED_IN in edge_types

    def test_validation_passes(self, engine: GraphEngine):
        errors = engine.validate()
        assert errors == []

    def test_key_nodes_present(self, engine: GraphEngine):
        for node_id in ["BaP", "CYP1A1", "CYP1B1", "GSTM1", "BPDE", "BPDE_dG", "XPC"]:
            assert engine.get_node(node_id) is not None, f"Missing node: {node_id}"


# ── Analysis ─────────────────────────────────────────────────────────────


class TestAnalysis:
    def test_shortest_path_activation_to_adduct(self, engine: GraphEngine):
        path = shortest_path(engine, "CYP1A1", "BPDE_dG")
        assert path is not None
        assert path[0] == "CYP1A1"
        assert path[-1] == "BPDE_dG"

    def test_all_shortest_paths(self, engine: GraphEngine):
        paths = all_shortest_paths(engine, "CYP1A1", "BPDE_dG")
        assert len(paths) >= 1

    def test_centrality_cyp1a1_high(self, engine: GraphEngine):
        scores = centrality(engine, method="degree")
        assert scores["CYP1A1"] > scores["Lung"]

    def test_betweenness_centrality(self, engine: GraphEngine):
        scores = centrality(engine, method="betweenness")
        assert len(scores) == engine.node_count

    def test_metabolism_chain(self, engine: GraphEngine):
        chain = metabolism_chain(engine, "BaP")
        assert len(chain.node_ids) >= 5
        assert len(chain.activation_edges) >= 2
        assert len(chain.detox_edges) >= 1
        assert len(chain.adduct_edges) >= 1
        assert len(chain.repair_edges) >= 1

    def test_pathway_subgraph(self, engine: GraphEngine):
        members = pathway_subgraph(engine, "hsa05204")
        assert "BaP" in members
        assert "CYP1A1" in members

    def test_variant_impact_cyp1a1(self, engine: GraphEngine):
        impact = variant_impact_score(engine, "CYP1A1")
        assert impact is not None
        assert impact.activity_score == 1.0
        assert impact.downstream_adduct_count >= 1
        assert impact.score >= 0

    def test_variant_impact_gstm1_null(self, engine: GraphEngine):
        impact = variant_impact_score(engine, "GSTM1")
        assert impact is not None
        assert impact.activity_score == 0.0


# ── Export & reimport ────────────────────────────────────────────────────


class TestExportReimport:
    def test_json_roundtrip(self, engine: GraphEngine, tmp_path):
        out = to_json(engine, tmp_path / "graph.json")
        data = json.loads(out.read_text())
        restored = KnowledgeGraph(**data)
        assert len(restored.nodes) == engine.node_count

    def test_json_roundtrip_via_engine(self, engine: GraphEngine):
        json_str = engine.to_json()
        data = json.loads(json_str)
        new_engine = GraphEngine()
        new_engine.load(KnowledgeGraph(**data))
        assert new_engine.node_count == engine.node_count
        assert new_engine.edge_count == engine.edge_count

    def test_graph_data_js_roundtrip(self, engine: GraphEngine, tmp_path):
        out = to_graph_data_js(engine, tmp_path / "graph-data.js")
        js_text = out.read_text()
        assert "GRAPH_DATA" in js_text
        restored = parse_graph_data_text(js_text)
        assert len(restored.nodes) == engine.node_count

    def test_gexf_export(self, engine: GraphEngine, tmp_path):
        out = to_gexf(engine, tmp_path / "graph.gexf")
        gexf_str = out.read_text()
        assert "<?xml" in gexf_str
        assert "BaP" in gexf_str
        assert "CYP1A1" in gexf_str

    def test_interactive_html_export(self, engine: GraphEngine):
        try:
            html = to_interactive_html_string(engine)
        except FileNotFoundError:
            pytest.skip("Viewer HTML template not found in this environment")
        assert "<html" in html.lower()
        assert "GRAPH_DATA" in html

    def test_reimported_graph_validates(self, engine: GraphEngine, tmp_path):
        out = to_json(engine, tmp_path / "graph.json")
        data = json.loads(out.read_text())
        restored = KnowledgeGraph(**data)
        new_engine = GraphEngine()
        new_engine.load(restored)
        assert new_engine.validate() == []


# ── IARC enrichment ──────────────────────────────────────────────────────


class TestIARCEnrichment:
    def test_bap_classified_group_1(self):
        clf = IARCClassifier()
        assert clf.classify("Benzo[a]pyrene") == IARCGroup.GROUP_1

    def test_enrichment_matches_graph(self, bap_kg: KnowledgeGraph):
        clf = IARCClassifier()
        bap_node = next(n for n in bap_kg.nodes if n.id == "BaP")
        iarc_entry = clf.get_entry("Benzo[a]pyrene")
        assert iarc_entry is not None
        assert bap_node.iarc == iarc_entry["group"]
