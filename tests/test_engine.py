"""Tests for ExposoGraph.engine."""

import json

import pytest
from pydantic import ValidationError

from ExposoGraph.config import GraphMode
from ExposoGraph.engine import GraphEngine
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, MatchStatus, Node, NodeType


@pytest.fixture
def engine():
    return GraphEngine()


@pytest.fixture
def sample_kg():
    return KnowledgeGraph(
        nodes=[
            Node(
                id="BaP",
                label="Benzo[a]pyrene",
                type=NodeType.CARCINOGEN,
                group="PAH",
                iarc="Group 1",
            ),
            Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME, phase="I", role="Activation"),
            Node(id="BPDE", label="BPDE", type=NodeType.METABOLITE, reactivity="High"),
        ],
        edges=[
            Edge(source="CYP1A1", target="BPDE", type=EdgeType.ACTIVATES, carcinogen="BaP"),
        ],
    )


class TestAddNode:
    def test_add_single(self, engine):
        node = Node(id="X", label="X", type=NodeType.ENZYME)
        engine.add_node(node)
        assert engine.node_count == 1

    def test_add_duplicate_overwrites(self, engine):
        engine.add_node(Node(id="X", label="X", type=NodeType.ENZYME, detail="v1"))
        engine.add_node(Node(id="X", label="X", type=NodeType.ENZYME, detail="v2"))
        assert engine.node_count == 1
        assert engine.get_node("X")["detail"] == "v2"


class TestAddEdge:
    def test_add_valid_edge(self, engine):
        engine.add_node(Node(id="A", label="A", type=NodeType.ENZYME))
        engine.add_node(Node(id="B", label="B", type=NodeType.METABOLITE))
        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES))
        assert engine.edge_count == 1

    def test_missing_source_raises(self, engine):
        engine.add_node(Node(id="B", label="B", type=NodeType.METABOLITE))
        with pytest.raises(ValueError, match="Missing source node"):
            engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES))

    def test_missing_target_raises(self, engine):
        engine.add_node(Node(id="A", label="A", type=NodeType.ENZYME))
        with pytest.raises(ValueError, match="Missing target node"):
            engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES))

    def test_missing_carcinogen_raises(self, engine):
        engine.add_node(Node(id="A", label="A", type=NodeType.ENZYME))
        engine.add_node(Node(id="B", label="B", type=NodeType.METABOLITE))
        with pytest.raises(ValueError, match="Missing carcinogen context"):
            engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES, carcinogen="BaP"))

    def test_parallel_edges_are_preserved(self, engine):
        engine.add_node(Node(id="A", label="A", type=NodeType.ENZYME))
        engine.add_node(Node(id="B", label="B", type=NodeType.METABOLITE))

        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES, pmid="1"))
        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES, pmid="2"))

        assert engine.edge_count == 2
        edge_pmids = {edge["pmid"] for edge in engine.to_dict()["edges"]}
        assert edge_pmids == {"1", "2"}


class TestLoadAndMerge:
    def test_load(self, engine, sample_kg):
        warnings = engine.load(sample_kg)
        assert engine.node_count == 3
        assert engine.edge_count == 1
        assert warnings == []

    def test_load_rejects_bad_edges(self, engine):
        with pytest.raises(ValidationError, match="MISSING"):
            KnowledgeGraph(
                nodes=[Node(id="A", label="A", type=NodeType.ENZYME)],
                edges=[Edge(source="A", target="MISSING", type=EdgeType.ACTIVATES)],
            )

    def test_merge_additive(self, engine, sample_kg):
        engine.load(sample_kg)
        extra = KnowledgeGraph(
            nodes=[
                Node(
                    id="XPC",
                    label="XPC",
                    type=NodeType.ENZYME,
                    role="Repair",
                    group="DNA Repair (NER)",
                ),
            ],
            edges=[],
        )
        engine.merge(extra)
        assert engine.node_count == 4

    def test_merge_exploratory_keeps_unmatched_content(self, engine):
        exploratory = KnowledgeGraph(
            nodes=[
                Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
                Node(id="n2", label="Unknown Chemical", type=NodeType.CARCINOGEN),
            ],
            edges=[Edge(source="n1", target="n2", type=EdgeType.ACTIVATES)],
        )

        warnings = engine.merge(exploratory, mode=GraphMode.EXPLORATORY)

        assert warnings == []
        assert engine.node_count == 2
        assert engine.edge_count == 1
        assert engine.get_node("n2")["match_status"] == MatchStatus.UNMATCHED.value

    def test_merge_strict_drops_unmatched_content(self, engine):
        exploratory = KnowledgeGraph(
            nodes=[
                Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
                Node(id="n2", label="Unknown Chemical", type=NodeType.CARCINOGEN),
            ],
            edges=[Edge(source="n1", target="n2", type=EdgeType.ACTIVATES)],
        )

        warnings = engine.merge(exploratory, mode=GraphMode.STRICT)

        assert engine.node_count == 1
        assert engine.edge_count == 0
        assert engine.get_node("n1") is not None
        assert engine.get_node("n2") is None
        assert any("non-canonical node" in warning for warning in warnings)
        assert any("non-canonical edge" in warning for warning in warnings)


class TestRemove:
    def test_remove_node(self, engine, sample_kg):
        engine.load(sample_kg)
        engine.remove_node("BaP")
        assert engine.node_count == 2

    def test_remove_nonexistent_node_noop(self, engine):
        engine.remove_node("DOES_NOT_EXIST")
        assert engine.node_count == 0

    def test_remove_edge(self, engine, sample_kg):
        engine.load(sample_kg)
        edge_key = next(iter(engine.G.edges(keys=True)))[2]
        engine.remove_edge("CYP1A1", "BPDE", edge_key)
        assert engine.edge_count == 0

    def test_remove_edge_with_unknown_key_noop(self, engine, sample_kg):
        engine.load(sample_kg)
        engine.remove_edge("CYP1A1", "BPDE", "CYP1A1-ACTIVATES-BPDE")
        assert engine.edge_count == 1


class TestQueries:
    def test_get_node(self, engine, sample_kg):
        engine.load(sample_kg)
        data = engine.get_node("BaP")
        assert data is not None
        assert data["type"] == "Carcinogen"

    def test_get_missing_node(self, engine):
        assert engine.get_node("X") is None

    def test_neighbors(self, engine, sample_kg):
        engine.load(sample_kg)
        nbrs = engine.neighbors("CYP1A1")
        assert "BPDE" in nbrs

    def test_nodes_by_type(self, engine, sample_kg):
        engine.load(sample_kg)
        enzymes = engine.nodes_by_type("Enzyme")
        assert len(enzymes) == 1
        assert enzymes[0]["id"] == "CYP1A1"


class TestSerialization:
    def test_to_dict(self, engine, sample_kg):
        engine.load(sample_kg)
        d = engine.to_dict()
        assert len(d["nodes"]) == 3
        assert len(d["edges"]) == 1

    def test_to_json(self, engine, sample_kg):
        engine.load(sample_kg)
        s = engine.to_json()
        data = json.loads(s)
        assert "nodes" in data
        assert "edges" in data

    def test_to_knowledge_graph(self, engine, sample_kg):
        engine.load(sample_kg)
        kg = engine.to_knowledge_graph()
        assert len(kg.nodes) == 3

    def test_clear(self, engine, sample_kg):
        engine.load(sample_kg)
        engine.clear()
        assert engine.node_count == 0
        assert engine.edge_count == 0


class TestValidation:
    def test_valid_graph(self, engine, sample_kg):
        engine.load(sample_kg)
        assert engine.validate() == []

    def test_dangling_carcinogen(self, engine):
        engine.add_node(Node(id="A", label="A", type=NodeType.ENZYME))
        engine.add_node(Node(id="B", label="B", type=NodeType.METABOLITE))
        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES))
        # Manually inject a bad carcinogen ref
        for _, _, _, data in engine.G.edges(keys=True, data=True):
            data["carcinogen"] = "GONE"
        errors = engine.validate()
        assert len(errors) == 1
        assert "GONE" in errors[0]
