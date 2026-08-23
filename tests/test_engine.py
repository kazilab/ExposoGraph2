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
        # graph-data.json was consolidated from feature_extend_knowledge_graph
        # (277 nodes / 461 edges, replacing the prior 231-node/335-edge legacy
        # set) plus 49 NodeType.SUBSTRATE nodes sourced from
        # interaction_parameters.json's substrate keys that have no existing
        # Carcinogen counterpart (verified against id/label/canonical_label),
        # for 325 nodes total. 461 legacy edges plus 58 enzyme->substrate
        # edges added for the competitive_inhibition pairs that had no
        # existing qualifying edge, plus 11 more enzyme->substrate edges added
        # for the phase2_conjugation pairs that likewise had no existing
        # qualifying edge (topology baked into graph-data.json; kinetics
        # populated below), for 530 edges total.
        assert engine.node_count == 325
        assert engine.edge_count == 530
        raw = engine.get_data("CYP1A1", key="tissue_weights_raw")
        normalized = engine.get_data("CYP1A1", key="tissue_weights")
        assert raw is not None and normalized is not None
        assert max(normalized.values()) == 1.0
        assert normalized["Liver"] == raw["Liver"] / max(raw.values())
        # Interaction-parameter kinetics are applied dynamically, not baked
        # into graph-data.json -- see _apply_interaction_parameters. This
        # pair already has a qualifying edge, so this is pure enrichment.
        activated = engine.get_edge("CYP2E1", "Benzene_oxide")
        assert activated["kinetics"]["product"] == "benzene_oxide"
        assert activated["kinetics"]["product_carcinogenic"] is True
        # New topology edge added this commit (no prior qualifying edge
        # existed for this pair) whose kinetics also come from
        # interaction_parameters.json via the same overlay method.
        new_edge = engine.get_edge("CYP2E1", "Ethanol")
        assert new_edge["type"] == "ACTIVATES"
        assert new_edge["kinetics"]["product"] == "acetaldehyde"

    def test_load_reference_graph_substrate_nodes(self, engine):
        # Substrate identity nodes only -- no Km/Vmax/kinetics data is baked
        # into graph-data.json; that stays sourced from
        # interaction_parameters.json at instantiation time (analogous to
        # tissue expression weights being sourced from
        # tissue_expression_data.json rather than baked in). See
        # docs/design/kg_parameter_loading_scope.md, Addendum 3.
        engine.load_reference_graph()
        substrate_nodes = engine.nodes_by_type(NodeType.SUBSTRATE)
        substrate_ids = [node["id"] for node in substrate_nodes]
        assert len(substrate_nodes) == 49
        assert len(set(substrate_ids)) == 49
        # trichloroethylene aliases the existing TCE Carcinogen node
        # (canonical_label="Trichloroethylene") and must NOT get its own
        # Substrate node.
        assert "trichloroethylene" not in substrate_ids
        assert "caffeine" in substrate_ids
        assert "naphthalene" in substrate_ids
        for node in substrate_nodes:
            assert node["type"] == NodeType.SUBSTRATE.value
            # No interaction-parameters data (Km/Vmax/product/etc.) was
            # copied onto the node itself.
            assert "tissue_weights" not in node
            assert "Km" not in node
            assert "Vmax" not in node


class TestMapViewerFiltering:
    """Tests for the new map-viewer filtering helpers used by ui_map_viewer.py.

    Fixture graph shape (see ``mv_engine``):
      BaP  (Carcinogen, group="PAH")   -> CYP1A1 (Enzyme) -> BPDE (Metabolite) -> DNAAdduct1 (DNA_Adduct)
      Benz (Carcinogen, group="PAH")   -> CYP2E1 (Enzyme) -> BQ   (Metabolite)
      As   (Carcinogen, group="Metal") -> GSTM1  (Enzyme, no outgoing edge -- leaf)
      CaffeineSub (Substrate) <- CYP1A1 (kept in graph, hidden from viewer by default)
      CYP1A1.tissue_weights = {"Liver": 0.9, "Lung": 0.1}
      CYP2E1.tissue_weights = {"Liver": 0.05}
      GSTM1.tissue_weights = {} (empty -- no entry for any tissue)
    """

    @pytest.fixture
    def mv_engine(self):
        engine = GraphEngine()
        engine.add_node(Node(id="BaP", label="Benzo[a]pyrene", type=NodeType.CARCINOGEN, group="PAH"))
        engine.add_node(Node(id="Benz", label="Benzene", type=NodeType.CARCINOGEN, group="PAH"))
        engine.add_node(Node(id="As", label="Arsenic", type=NodeType.CARCINOGEN, group="Metal"))
        engine.add_node(
            Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME, tissue_weights={"Liver": 0.9, "Lung": 0.1})
        )
        engine.add_node(Node(id="CYP2E1", label="CYP2E1", type=NodeType.ENZYME, tissue_weights={"Liver": 0.05}))
        engine.add_node(Node(id="GSTM1", label="GSTM1", type=NodeType.ENZYME))
        engine.add_node(Node(id="BPDE", label="BPDE", type=NodeType.METABOLITE))
        engine.add_node(Node(id="BQ", label="Benzoquinone", type=NodeType.METABOLITE))
        engine.add_node(Node(id="DNAAdduct1", label="BPDE-dG adduct", type=NodeType.DNA_ADDUCT))
        engine.add_node(Node(id="CaffeineSub", label="Caffeine", type=NodeType.SUBSTRATE))
        engine.add_edge(Edge(source="BaP", target="CYP1A1", type=EdgeType.ACTIVATES, carcinogen="BaP"))
        engine.add_edge(Edge(source="CYP1A1", target="BPDE", type=EdgeType.ACTIVATES, carcinogen="BaP"))
        engine.add_edge(Edge(source="BPDE", target="DNAAdduct1", type=EdgeType.FORMS_ADDUCT, carcinogen="BaP"))
        engine.add_edge(Edge(source="Benz", target="CYP2E1", type=EdgeType.ACTIVATES, carcinogen="Benz"))
        engine.add_edge(Edge(source="CYP2E1", target="BQ", type=EdgeType.ACTIVATES, carcinogen="Benz"))
        engine.add_edge(Edge(source="As", target="GSTM1", type=EdgeType.DETOXIFIES, carcinogen="As"))
        engine.add_edge(Edge(source="CYP1A1", target="CaffeineSub", type=EdgeType.TRANSPORTS))
        return engine

    # ── carcinogen_group_paths_subgraph ──────────────────────────────────

    def test_group_paths_union_across_carcinogens_in_group(self, mv_engine):
        result = mv_engine.carcinogen_group_paths_subgraph(["PAH"])
        node_ids = {n["id"] for n in result["nodes"]}
        edge_pairs = {(e["source"], e["target"]) for e in result["edges"]}
        # Both BaP's and Benz's paths are unioned since both are in group "PAH".
        # CYP1A1 has two outgoing edges (to BPDE and to CaffeineSub), so BaP
        # has two maximal paths and both are included.
        assert node_ids == {"BaP", "CYP1A1", "BPDE", "DNAAdduct1", "CaffeineSub", "Benz", "CYP2E1", "BQ"}
        assert edge_pairs == {
            ("BaP", "CYP1A1"),
            ("CYP1A1", "BPDE"),
            ("BPDE", "DNAAdduct1"),
            ("CYP1A1", "CaffeineSub"),
            ("Benz", "CYP2E1"),
            ("CYP2E1", "BQ"),
        }
        # Metal-group carcinogen As is not pulled in.
        assert "As" not in node_ids and "GSTM1" not in node_ids

    def test_group_paths_leaf_carcinogen_keeps_its_own_node(self, mv_engine):
        # As -> GSTM1 is a real path, so As, GSTM1 pass; but GSTM1 has no
        # further outgoing edges (leaf), which must not drop either node.
        result = mv_engine.carcinogen_group_paths_subgraph(["Metal"])
        node_ids = {n["id"] for n in result["nodes"]}
        assert node_ids == {"As", "GSTM1"}

    def test_group_paths_unknown_group_is_empty_not_error(self, mv_engine):
        result = mv_engine.carcinogen_group_paths_subgraph(["NoSuchGroup"])
        assert result == {"nodes": [], "edges": []}

    def test_group_paths_empty_groups_returns_empty(self, mv_engine):
        assert mv_engine.carcinogen_group_paths_subgraph([]) == {"nodes": [], "edges": []}

    # ── subgraph_by_node_types ────────────────────────────────────────────

    def test_node_types_none_excludes_substrate_by_default(self, mv_engine):
        result = mv_engine.subgraph_by_node_types()
        node_ids = {n["id"] for n in result["nodes"]}
        assert "CaffeineSub" not in node_ids
        assert node_ids == {"BaP", "Benz", "As", "CYP1A1", "CYP2E1", "GSTM1", "BPDE", "BQ", "DNAAdduct1"}
        # The CYP1A1 -> CaffeineSub edge drops too since its target is excluded.
        assert ("CYP1A1", "CaffeineSub") not in {(e["source"], e["target"]) for e in result["edges"]}

    def test_node_types_restricts_to_selected_types_only(self, mv_engine):
        result = mv_engine.subgraph_by_node_types([NodeType.ENZYME.value])
        node_ids = {n["id"] for n in result["nodes"]}
        assert node_ids == {"CYP1A1", "CYP2E1", "GSTM1"}
        # No direct Enzyme-Enzyme edges exist, so the induced edge set is empty.
        assert result["edges"] == []

    def test_node_types_exclude_types_override_still_applies(self, mv_engine):
        # Even explicitly requesting Substrate nodes, exclude_types defaults
        # to hiding them unless the caller overrides exclude_types too.
        result = mv_engine.subgraph_by_node_types([NodeType.SUBSTRATE.value])
        assert result == {"nodes": [], "edges": []}

    def test_node_types_can_opt_in_to_substrate_via_exclude_types(self, mv_engine):
        result = mv_engine.subgraph_by_node_types([NodeType.SUBSTRATE.value], exclude_types=())
        assert {n["id"] for n in result["nodes"]} == {"CaffeineSub"}

    # ── map_viewer_subgraph ───────────────────────────────────────────────

    def test_map_viewer_no_filters_is_full_graph_minus_substrate(self, mv_engine):
        result = mv_engine.map_viewer_subgraph()
        node_ids = {n["id"] for n in result["nodes"]}
        assert "CaffeineSub" not in node_ids
        assert len(node_ids) == 9

    def test_map_viewer_intersects_type_and_carcinogen_axes(self, mv_engine):
        # Enzyme type-filter intersected with PAH carcinogen-group filter:
        # only the enzymes reachable from a PAH-group carcinogen survive.
        result = mv_engine.map_viewer_subgraph(
            node_types=[NodeType.ENZYME.value], carcinogen_groups=["PAH"]
        )
        node_ids = {n["id"] for n in result["nodes"]}
        assert node_ids == {"CYP1A1", "CYP2E1"}
        assert "GSTM1" not in node_ids  # Metal group, filtered out by carcinogen axis.

    def test_map_viewer_type_filter_alone_ignores_carcinogen_axis(self, mv_engine):
        result = mv_engine.map_viewer_subgraph(node_types=[NodeType.CARCINOGEN.value])
        assert {n["id"] for n in result["nodes"]} == {"BaP", "Benz", "As"}

    def test_map_viewer_carcinogen_filter_alone_ignores_type_axis(self, mv_engine):
        result = mv_engine.map_viewer_subgraph(carcinogen_groups=["Metal"])
        assert {n["id"] for n in result["nodes"]} == {"As", "GSTM1"}

    def test_map_viewer_substrate_always_excluded_even_with_filters(self, mv_engine):
        result = mv_engine.map_viewer_subgraph(node_types=[NodeType.SUBSTRATE.value, NodeType.ENZYME.value])
        node_ids = {n["id"] for n in result["nodes"]}
        assert "CaffeineSub" not in node_ids
        assert node_ids == {"CYP1A1", "CYP2E1", "GSTM1"}

    # ── dim_by_tissue_threshold ───────────────────────────────────────────

    def test_dim_by_tissue_threshold_marks_low_expression_enzymes(self, mv_engine):
        subgraph = mv_engine.map_viewer_subgraph()
        dimmed = mv_engine.dim_by_tissue_threshold(subgraph, "Liver", 0.5)
        dimmed_by_id = {n["id"]: n["_dimmed"] for n in dimmed["nodes"]}
        assert dimmed_by_id["CYP1A1"] is False  # 0.9 >= 0.5
        assert dimmed_by_id["CYP2E1"] is True  # 0.05 < 0.5
        assert dimmed_by_id["GSTM1"] is True  # no Liver entry at all
        # Non-enzyme nodes are never dimmed.
        assert dimmed_by_id["BaP"] is False
        assert dimmed_by_id["BPDE"] is False

    def test_dim_by_tissue_threshold_dims_edges_touching_dimmed_enzyme(self, mv_engine):
        subgraph = mv_engine.map_viewer_subgraph()
        dimmed = mv_engine.dim_by_tissue_threshold(subgraph, "Liver", 0.5)
        dimmed_edges = {(e["source"], e["target"]): e["_dimmed"] for e in dimmed["edges"]}
        assert dimmed_edges[("Benz", "CYP2E1")] is True  # touches dimmed CYP2E1
        assert dimmed_edges[("CYP2E1", "BQ")] is True
        assert dimmed_edges[("BaP", "CYP1A1")] is False  # touches non-dimmed CYP1A1
        assert dimmed_edges[("CYP1A1", "BPDE")] is False

    def test_dim_by_tissue_threshold_never_removes_nodes_or_edges(self, mv_engine):
        subgraph = mv_engine.map_viewer_subgraph()
        dimmed = mv_engine.dim_by_tissue_threshold(subgraph, "Liver", 0.5)
        assert len(dimmed["nodes"]) == len(subgraph["nodes"])
        assert len(dimmed["edges"]) == len(subgraph["edges"])

    def test_dim_by_tissue_threshold_is_non_mutating(self, mv_engine):
        subgraph = mv_engine.map_viewer_subgraph()
        mv_engine.dim_by_tissue_threshold(subgraph, "Liver", 0.5)
        assert all("_dimmed" not in n for n in subgraph["nodes"])
        assert all("_dimmed" not in e for e in subgraph["edges"])
