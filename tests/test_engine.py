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

    def test_get_node_nested_key_dot_path(self, engine):
        engine.add_node(
            Node(
                id="CYP1A1",
                label="CYP1A1",
                type=NodeType.ENZYME,
                tissue_weights={"Liver": 1.0, "Lung": 0.25},
            )
        )
        assert engine.get_node("CYP1A1", "tissue_weights.Liver") == 1.0
        assert engine.get_node("CYP1A1", ("tissue_weights", "Lung")) == 0.25

    def test_get_node_nested_key_missing_segment_returns_default(self, engine):
        engine.add_node(
            Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME, tissue_weights={"Liver": 1.0})
        )
        assert engine.get_node("CYP1A1", "tissue_weights.Kidney") is None
        assert engine.get_node("CYP1A1", "tissue_weights.Kidney", default=0.0) == 0.0
        assert engine.get_node("CYP1A1", "does_not_exist") is None

    def test_get_node_with_key_on_missing_node_returns_default(self, engine):
        assert engine.get_node("X", "tissue_weights.Liver") is None
        assert engine.get_node("X", "tissue_weights.Liver", default="n/a") == "n/a"

    def test_get_edge_full_dict(self, engine, sample_kg):
        engine.load(sample_kg)
        data = engine.get_edge("CYP1A1", "BPDE")
        assert data is not None
        assert data["type"] == "ACTIVATES"

    def test_get_missing_edge(self, engine):
        assert engine.get_edge("X", "Y") is None

    def test_get_edge_nested_key(self, engine):
        engine.add_node(Node(id="NDMA", label="NDMA", type=NodeType.CARCINOGEN))
        engine.add_node(Node(id="CYP2E1", label="CYP2E1", type=NodeType.ENZYME))
        engine.add_edge(
            Edge(
                source="NDMA",
                target="CYP2E1",
                type=EdgeType.ACTIVATES,
                kinetics={"Km_uM": 12.5, "Ki": {"acetaldehyde": 3.1, "ethanol": 8.4}},
            )
        )
        assert engine.get_edge("NDMA", "CYP2E1", "kinetics.Km_uM") == 12.5
        assert engine.get_edge("NDMA", "CYP2E1", ("kinetics", "Ki", "ethanol")) == 8.4
        assert engine.get_edge("NDMA", "CYP2E1", "kinetics.Ki.missing") is None
        assert engine.get_edge("NDMA", "CYP2E1", "kinetics.Ki.missing", default=0.0) == 0.0

    def test_get_edge_with_key_on_missing_edge_returns_default(self, engine):
        assert engine.get_edge("X", "Y", "kinetics.Km_uM") is None
        assert engine.get_edge("X", "Y", "kinetics.Km_uM", default=-1) == -1

    def test_get_edge_keys_and_edge_key_disambiguation(self, engine):
        engine.add_node(Node(id="A", label="A", type=NodeType.ENZYME))
        engine.add_node(Node(id="B", label="B", type=NodeType.METABOLITE))
        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES, pmid="1"))
        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES, pmid="2"))

        keys = engine.get_edge_keys("A", "B")
        assert len(keys) == 2
        assert {engine.get_edge("A", "B", edge_key=k)["pmid"] for k in keys} == {"1", "2"}

    def test_get_edge_keys_empty_when_no_edge(self, engine):
        assert engine.get_edge_keys("X", "Y") == []

    def test_get_data_routes_to_node_when_target_omitted(self, engine, sample_kg):
        engine.load(sample_kg)
        assert engine.get_data("BaP") == engine.get_node("BaP")

    def test_get_data_routes_to_node_with_nested_key(self, engine):
        engine.add_node(
            Node(
                id="CYP1A1",
                label="CYP1A1",
                type=NodeType.ENZYME,
                tissue_weights={"Liver": 1.0, "Lung": 0.25},
            )
        )
        assert engine.get_data("CYP1A1", key="tissue_weights.Liver") == 1.0
        assert engine.get_data("CYP1A1", key=("tissue_weights", "Lung")) == 0.25
        assert engine.get_data("CYP1A1", key="tissue_weights.Kidney", default=0.0) == 0.0

    def test_get_data_routes_to_edge_when_target_given(self, engine, sample_kg):
        engine.load(sample_kg)
        assert engine.get_data("CYP1A1", "BPDE") == engine.get_edge("CYP1A1", "BPDE")

    def test_get_data_routes_to_edge_with_nested_key(self, engine):
        engine.add_node(Node(id="NDMA", label="NDMA", type=NodeType.CARCINOGEN))
        engine.add_node(Node(id="CYP2E1", label="CYP2E1", type=NodeType.ENZYME))
        engine.add_edge(
            Edge(
                source="NDMA",
                target="CYP2E1",
                type=EdgeType.ACTIVATES,
                kinetics={"Km_uM": 12.5, "Ki": {"acetaldehyde": 3.1, "ethanol": 8.4}},
            )
        )
        assert engine.get_data("NDMA", "CYP2E1", key="kinetics.Km_uM") == 12.5
        assert engine.get_data("NDMA", "CYP2E1", key=("kinetics", "Ki", "ethanol")) == 8.4
        assert engine.get_data("NDMA", "CYP2E1", key="kinetics.Ki.missing", default=0.0) == 0.0

    def test_get_data_edge_key_disambiguation(self, engine):
        engine.add_node(Node(id="A", label="A", type=NodeType.ENZYME))
        engine.add_node(Node(id="B", label="B", type=NodeType.METABOLITE))
        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES, pmid="1"))
        engine.add_edge(Edge(source="A", target="B", type=EdgeType.ACTIVATES, pmid="2"))

        keys = engine.get_edge_keys("A", "B")
        assert {engine.get_data("A", "B", edge_key=k)["pmid"] for k in keys} == {"1", "2"}

    def test_get_data_missing_node_and_edge_defaults(self, engine):
        assert engine.get_data("X") is None
        assert engine.get_data("X", key="tissue_weights.Liver", default="n/a") == "n/a"
        assert engine.get_data("X", "Y") is None
        assert engine.get_data("X", "Y", key="kinetics.Km_uM", default=-1) == -1

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

    def test_to_knowledge_graph_preserves_tissue_weights_raw(self, engine):
        engine.load_reference_graph()
        kg = engine.to_knowledge_graph()
        cyp1a1 = next(n for n in kg.nodes if n.id == "CYP1A1")
        assert cyp1a1.tissue_weights_raw is not None
        assert cyp1a1.tissue_weights_raw == engine.get_data("CYP1A1", key="tissue_weights_raw")

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


class TestLoadReferenceGraph:
    @pytest.fixture
    def graph_data_path(self, tmp_path):
        kg = KnowledgeGraph(
            nodes=[
                Node(
                    id="CYP1A1",
                    label="CYP1A1",
                    type=NodeType.ENZYME,
                    # Stale node-level tissue weights, as if baked in by an
                    # older export -- _apply_tissue_expression should overwrite this.
                    tissue_weights={"Liver": 0.1, "Lung": 0.9},
                ),
                Node(id="NOEXPR", label="NOEXPR", type=NodeType.ENZYME),
                Node(id="BPDE", label="BPDE", type=NodeType.METABOLITE),
            ],
            edges=[Edge(source="CYP1A1", target="BPDE", type=EdgeType.ACTIVATES)],
        )
        path = tmp_path / "graph-data.json"
        path.write_text(json.dumps(kg.model_dump(mode="json")))
        return path

    @pytest.fixture
    def tissue_expression_path(self, tmp_path):
        data = {
            "metadata": {"tissues": ["Liver", "Lung"]},
            "expression": {"CYP1A1": {"Liver": 40.0, "Lung": 10.0}},
        }
        path = tmp_path / "tissue_expression_data.json"
        path.write_text(json.dumps(data))
        return path

    def test_apply_tissue_expression_overwrites_stale_node_data(
        self, engine, graph_data_path, tissue_expression_path
    ):
        from ExposoGraph.exporter import parse_graph_artifact

        engine.load(parse_graph_artifact(graph_data_path))
        assert engine.get_data("CYP1A1", key="tissue_weights.Lung") == 0.9  # stale, pre-overwrite

        warnings = engine._apply_tissue_expression(tissue_expression_path)

        assert engine.get_data("CYP1A1", key="tissue_weights_raw") == {"Liver": 40.0, "Lung": 10.0}
        assert engine.get_data("CYP1A1", key="tissue_weights") == {"Liver": 1.0, "Lung": 0.25}
        assert "No tissue expression data for enzyme: NOEXPR" in warnings

    def test_apply_tissue_expression_clears_missing_enzyme_stale_data(self, engine):
        engine.add_node(
            Node(
                id="NOEXPR",
                label="NOEXPR",
                type=NodeType.ENZYME,
                tissue_weights={"Liver": 0.5},
            )
        )
        engine._apply_tissue_expression()  # uses the bundled tissue_expression_data.json
        assert engine.get_data("NOEXPR", key="tissue_weights") is None
        assert engine.get_data("NOEXPR", key="tissue_weights_raw") is None

    def test_load_reference_graph_with_custom_paths(
        self, engine, graph_data_path, tissue_expression_path
    ):
        warnings = engine.load_reference_graph(
            graph_data_path=graph_data_path,
            tissue_expression_path=tissue_expression_path,
        )
        assert engine.node_count == 3
        assert engine.edge_count == 1
        assert engine.get_data("CYP1A1", key="tissue_weights") == {"Liver": 1.0, "Lung": 0.25}
        # The edge itself carries no tissue-weight attributes -- they live on the node.
        edge = engine.get_data("CYP1A1", "BPDE")
        assert "tissue_weights" not in edge
        assert "tissue_weights_raw" not in edge
        assert "No tissue expression data for enzyme: NOEXPR" in warnings

    def test_load_reference_graph_defaults_to_bundled_files(self, engine):
        engine.load_reference_graph()
        assert engine.node_count == 231
        assert engine.edge_count == 335
        raw = engine.get_data("CYP1A1", key="tissue_weights_raw")
        normalized = engine.get_data("CYP1A1", key="tissue_weights")
        assert raw is not None and normalized is not None
        assert max(normalized.values()) == 1.0
        assert normalized["Liver"] == raw["Liver"] / max(raw.values())
