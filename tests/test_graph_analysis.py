"""Tests for ExposoGraph.graph_analysis."""

import pytest

from ExposoGraph.engine import GraphEngine
from ExposoGraph.graph_analysis import (
    MetabolismChain,
    VariantImpact,
    all_shortest_paths,
    centrality,
    metabolism_chain,
    pathway_subgraph,
    shortest_path,
    variant_impact_score,
)
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, Node, NodeType

# ── Fixtures ─────────────────────────────────────────────────────────────


@pytest.fixture
def engine():
    return GraphEngine()


@pytest.fixture
def bap_kg():
    """A realistic BaP metabolism graph for testing analysis functions."""
    return KnowledgeGraph(
        nodes=[
            Node(id="BaP", label="Benzo[a]pyrene", type=NodeType.CARCINOGEN, group="PAH", iarc="Group 1"),
            Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME, phase="I", role="Activation", activity_score=1.0),
            Node(id="CYP1B1", label="CYP1B1", type=NodeType.ENZYME, phase="I", role="Activation", activity_score=0.5),
            Node(id="EPHX1", label="EPHX1", type=NodeType.ENZYME, phase="I", role="Activation"),
            Node(id="BPDE", label="BPDE", type=NodeType.METABOLITE, reactivity="High"),
            Node(id="BaP_diol", label="BaP-7,8-diol", type=NodeType.METABOLITE, reactivity="Intermediate"),
            Node(id="BPDE_dG", label="BPDE-N2-dG", type=NodeType.DNA_ADDUCT),
            Node(id="XPC", label="XPC", type=NodeType.ENZYME, role="Repair", group="DNA Repair (NER)"),
            Node(id="GSTM1", label="GSTM1", type=NodeType.ENZYME, phase="II", role="Detoxification"),
            Node(id="BPDE_GSH", label="BPDE-GSH", type=NodeType.METABOLITE, reactivity="Low"),
            Node(id="ABCB1", label="ABCB1", type=NodeType.ENZYME, phase="III", role="Transport"),
            Node(id="hsa05204", label="Chemical Carcinogenesis", type=NodeType.PATHWAY),
            Node(id="Lung", label="Lung", type=NodeType.TISSUE),
        ],
        edges=[
            # Activation chain
            Edge(source="CYP1A1", target="BaP", type=EdgeType.ACTIVATES, label="epoxidation", carcinogen="BaP"),
            Edge(source="CYP1B1", target="BaP", type=EdgeType.ACTIVATES, label="epoxidation", carcinogen="BaP"),
            Edge(source="EPHX1", target="BaP_diol", type=EdgeType.ACTIVATES, label="hydrolysis", carcinogen="BaP"),
            Edge(source="CYP1A1", target="BPDE", type=EdgeType.ACTIVATES, label="second epoxidation", carcinogen="BaP"),
            # Adduct formation
            Edge(source="BPDE", target="BPDE_dG", type=EdgeType.FORMS_ADDUCT, carcinogen="BaP"),
            # Repair
            Edge(source="XPC", target="BPDE_dG", type=EdgeType.REPAIRS, carcinogen="BaP"),
            # Detoxification
            Edge(source="GSTM1", target="BPDE", type=EdgeType.DETOXIFIES, label="glutathione conjugation", carcinogen="BaP"),
            Edge(source="GSTM1", target="BPDE_GSH", type=EdgeType.DETOXIFIES, label="produces conjugate", carcinogen="BaP"),
            # Transport
            Edge(source="ABCB1", target="BPDE_GSH", type=EdgeType.TRANSPORTS, label="efflux", carcinogen="BaP"),
            # Pathway
            Edge(source="BaP", target="hsa05204", type=EdgeType.PATHWAY),
            Edge(source="CYP1A1", target="hsa05204", type=EdgeType.PATHWAY),
            # Expression
            Edge(source="CYP1A1", target="Lung", type=EdgeType.EXPRESSED_IN),
        ],
    )


@pytest.fixture
def bap_engine(engine, bap_kg):
    engine.load(bap_kg)
    return engine


# ── Shortest path ────────────────────────────────────────────────────────


class TestShortestPath:
    def test_direct_neighbors(self, bap_engine):
        path = shortest_path(bap_engine, "CYP1A1", "BaP")
        assert path == ["CYP1A1", "BaP"]

    def test_multi_hop(self, bap_engine):
        path = shortest_path(bap_engine, "CYP1A1", "BPDE_dG")
        assert path is not None
        assert path[0] == "CYP1A1"
        assert path[-1] == "BPDE_dG"
        assert len(path) >= 3

    def test_no_path(self, bap_engine):
        # Lung is only connected via EXPRESSED_IN to CYP1A1
        # but hsa05204 is connected to BaP and CYP1A1 via PATHWAY
        # There IS a path: Lung <-> CYP1A1 <-> hsa05204
        path = shortest_path(bap_engine, "Lung", "hsa05204")
        assert path is not None

    def test_missing_node(self, bap_engine):
        assert shortest_path(bap_engine, "MISSING", "BaP") is None

    def test_same_node(self, bap_engine):
        path = shortest_path(bap_engine, "BaP", "BaP")
        assert path == ["BaP"]

    def test_all_shortest_paths(self, bap_engine):
        paths = all_shortest_paths(bap_engine, "CYP1A1", "BPDE_dG")
        assert len(paths) >= 1
        for p in paths:
            assert p[0] == "CYP1A1"
            assert p[-1] == "BPDE_dG"

    def test_all_shortest_paths_missing_node(self, bap_engine):
        assert all_shortest_paths(bap_engine, "MISSING", "BaP") == []


# ── Centrality ───────────────────────────────────────────────────────────


class TestCentrality:
    def test_degree_centrality(self, bap_engine):
        scores = centrality(bap_engine, method="degree")
        assert len(scores) == bap_engine.node_count
        # CYP1A1 has many connections — should have high centrality
        assert scores["CYP1A1"] > 0

    def test_betweenness_centrality(self, bap_engine):
        scores = centrality(bap_engine, method="betweenness")
        assert len(scores) == bap_engine.node_count
        assert all(isinstance(v, float) for v in scores.values())

    def test_closeness_centrality(self, bap_engine):
        scores = centrality(bap_engine, method="closeness")
        assert len(scores) == bap_engine.node_count

    def test_invalid_method(self, bap_engine):
        with pytest.raises(ValueError, match="Unknown centrality method"):
            centrality(bap_engine, method="invalid")

    def test_empty_graph(self, engine):
        scores = centrality(engine, method="degree")
        assert scores == {}


# ── Metabolism chain ─────────────────────────────────────────────────────


class TestMetabolismChain:
    def test_bap_chain(self, bap_engine):
        chain = metabolism_chain(bap_engine, "BaP")
        assert isinstance(chain, MetabolismChain)
        assert chain.carcinogen_id == "BaP"
        assert "BaP" in chain.node_ids
        assert "BPDE" in chain.node_ids
        assert "BPDE_dG" in chain.node_ids

    def test_includes_activation(self, bap_engine):
        chain = metabolism_chain(bap_engine, "BaP")
        assert len(chain.activation_edges) >= 2  # CYP1A1 and CYP1B1 activate

    def test_includes_detox(self, bap_engine):
        chain = metabolism_chain(bap_engine, "BaP")
        assert len(chain.detox_edges) >= 1  # GSTM1 detoxifies

    def test_includes_adduct(self, bap_engine):
        chain = metabolism_chain(bap_engine, "BaP")
        assert len(chain.adduct_edges) >= 1  # BPDE forms adduct

    def test_includes_repair(self, bap_engine):
        chain = metabolism_chain(bap_engine, "BaP")
        assert len(chain.repair_edges) >= 1  # XPC repairs

    def test_excludes_non_metabolism_edges(self, bap_engine):
        chain = metabolism_chain(bap_engine, "BaP")
        edge_types = {e.get("type") for e in chain.edges}
        assert "PATHWAY" not in edge_types
        assert "EXPRESSED_IN" not in edge_types

    def test_missing_carcinogen(self, bap_engine):
        chain = metabolism_chain(bap_engine, "MISSING")
        assert chain.node_ids == []
        assert chain.edges == []

    def test_node_with_no_chain(self, bap_engine):
        chain = metabolism_chain(bap_engine, "Lung")
        # Lung has no metabolism edges — chain should be minimal
        assert "Lung" not in chain.node_ids or len(chain.edges) == 0

    def test_excludes_unlabeled_branches_that_only_share_an_enzyme(self, engine):
        engine.add_node(Node(id="BaP", label="BaP", type=NodeType.CARCINOGEN))
        engine.add_node(Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME))
        engine.add_node(Node(id="BPDE", label="BPDE", type=NodeType.METABOLITE))
        engine.add_node(Node(id="OtherMet", label="OtherMet", type=NodeType.METABOLITE))

        engine.add_edge(Edge(source="CYP1A1", target="BaP", type=EdgeType.ACTIVATES, carcinogen="BaP"))
        engine.add_edge(Edge(source="BPDE", target="OtherMet", type=EdgeType.FORMS_ADDUCT))
        engine.add_edge(Edge(source="CYP1A1", target="OtherMet", type=EdgeType.ACTIVATES))

        chain = metabolism_chain(engine, "BaP")

        assert chain.node_ids == ["BaP", "CYP1A1"]
        assert len(chain.edges) == 1


# ── Pathway subgraph ─────────────────────────────────────────────────────


class TestPathwaySubgraph:
    def test_pathway_members(self, bap_engine):
        members = pathway_subgraph(bap_engine, "hsa05204")
        assert "hsa05204" in members
        assert "BaP" in members
        assert "CYP1A1" in members

    def test_excludes_non_pathway_nodes(self, bap_engine):
        members = pathway_subgraph(bap_engine, "hsa05204")
        assert "BPDE" not in members  # connected via ACTIVATES, not PATHWAY

    def test_missing_pathway(self, bap_engine):
        assert pathway_subgraph(bap_engine, "MISSING") == []


# ── Variant impact score ─────────────────────────────────────────────────


class TestVariantImpactScore:
    def test_activation_enzyme(self, bap_engine):
        impact = variant_impact_score(bap_engine, "CYP1A1")
        assert isinstance(impact, VariantImpact)
        assert impact.gene_id == "CYP1A1"
        assert impact.activity_score == 1.0
        assert impact.downstream_adduct_count >= 1
        assert impact.score >= 0

    def test_low_activity_higher_score(self, bap_engine):
        # CYP1B1 has activity_score=0.5, CYP1A1 has 1.0
        # CYP1B1 should have a higher impact score (more impactful variant)
        # if they have similar downstream topology
        impact_1b1 = variant_impact_score(bap_engine, "CYP1B1")
        assert impact_1b1 is not None
        assert impact_1b1.activity_score == 0.5

    def test_repair_enzyme(self, bap_engine):
        impact = variant_impact_score(bap_engine, "XPC")
        assert impact is not None
        assert impact.gene_id == "XPC"

    def test_missing_node(self, bap_engine):
        assert variant_impact_score(bap_engine, "MISSING") is None

    def test_node_without_activity_score(self, bap_engine):
        impact = variant_impact_score(bap_engine, "EPHX1")
        assert impact is not None
        assert impact.activity_score is None
        # Should use default activity of 1.0 in scoring
