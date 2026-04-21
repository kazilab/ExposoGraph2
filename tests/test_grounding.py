"""Tests for ExposoGraph.grounding."""

from ExposoGraph.config import GraphMode
from ExposoGraph.grounding import (
    build_default_grounding_index,
    build_iarc_grounding_index,
    ground_knowledge_graph,
    ground_node,
    normalize_grounding_key,
    prepare_knowledge_graph,
)
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, MatchStatus, Node, NodeType


class TestNormalizeGroundingKey:
    def test_normalizes_common_bap_variants(self):
        assert normalize_grounding_key("Benzo[a]pyrene") == normalize_grounding_key(
            "Benzo(a)pyrene"
        )
        assert normalize_grounding_key("CYP1A1") == normalize_grounding_key("cyp1a1")


class TestGroundNode:
    def test_gene_symbol_matches_reference_panel(self):
        node = Node(id="user_gene", label="CYP1A1", type=NodeType.ENZYME)

        grounded = ground_node(node)

        assert grounded.match_status == MatchStatus.CANONICAL
        assert grounded.canonical_id == "CYP1A1"
        assert grounded.canonical_label == "CYP1A1"
        assert grounded.canonical_namespace == "reference_panel"
        assert grounded.tier == 1

    def test_gene_alias_matches_reference_panel(self):
        node = Node(id="repair_alias", label="XPD", type=NodeType.GENE)

        grounded = ground_node(node)

        assert grounded.match_status == MatchStatus.ALIAS
        assert grounded.canonical_id == "ERCC2"
        assert grounded.canonical_label == "ERCC2/XPD"
        assert grounded.canonical_namespace == "reference_panel"

    def test_chemical_alias_matches_iarc(self):
        node = Node(id="chem1", label="BaP", type=NodeType.CARCINOGEN)

        grounded = ground_node(node)

        assert grounded.match_status == MatchStatus.ALIAS
        assert grounded.canonical_namespace == "iarc"
        assert grounded.canonical_id == "50-32-8"
        assert grounded.canonical_label == "Benzo[a]pyrene"
        assert grounded.iarc == "Group 1"
        assert grounded.group == "PAH"

    def test_nnk_matches_new_iarc_nitrosamine_entry(self):
        node = Node(id="chem2", label="NNK", type=NodeType.CARCINOGEN)

        grounded = ground_node(node)

        assert grounded.match_status == MatchStatus.ALIAS
        assert grounded.canonical_namespace == "iarc"
        assert grounded.iarc == "Group 1"
        assert grounded.group == "Nitrosamine"

    def test_custom_node_is_left_untouched(self):
        node = Node(
            id="novel",
            label="Novel Exposure",
            type=NodeType.CARCINOGEN,
            match_status=MatchStatus.CUSTOM,
            custom_type="Exposure",
        )

        grounded = ground_node(node)

        assert grounded.match_status == MatchStatus.CUSTOM
        assert grounded.custom_type == "Exposure"
        assert grounded.canonical_id is None

    def test_unknown_node_becomes_unmatched(self):
        node = Node(id="mystery", label="Completely Unknown Molecule", type=NodeType.CARCINOGEN)

        grounded = ground_node(node)

        assert grounded.match_status == MatchStatus.UNMATCHED


class TestGroundKnowledgeGraph:
    def test_grounded_nodes_produce_canonical_edge(self):
        kg = KnowledgeGraph(
            nodes=[
                Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
                Node(id="n2", label="BaP", type=NodeType.CARCINOGEN),
            ],
            edges=[
                Edge(source="n1", target="n2", type=EdgeType.ACTIVATES),
            ],
        )

        grounded = ground_knowledge_graph(kg)

        assert grounded.nodes[0].match_status == MatchStatus.CANONICAL
        assert grounded.nodes[1].match_status == MatchStatus.ALIAS
        assert grounded.edges[0].match_status == MatchStatus.CANONICAL
        assert grounded.edges[0].canonical_predicate == "ACTIVATES"

    def test_edge_with_unmatched_endpoint_stays_unmatched(self):
        kg = KnowledgeGraph(
            nodes=[
                Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
                Node(id="n2", label="Unknown Chemical", type=NodeType.CARCINOGEN),
            ],
            edges=[
                Edge(source="n1", target="n2", type=EdgeType.ACTIVATES),
            ],
        )

        grounded = ground_knowledge_graph(kg)

        assert grounded.nodes[0].match_status == MatchStatus.CANONICAL
        assert grounded.nodes[1].match_status == MatchStatus.UNMATCHED
        assert grounded.edges[0].match_status == MatchStatus.UNMATCHED

    def test_custom_edge_becomes_custom_when_grounded(self):
        kg = KnowledgeGraph(
            nodes=[
                Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
                Node(id="n2", label="BaP", type=NodeType.CARCINOGEN),
            ],
            edges=[
                Edge(
                    source="n1",
                    target="n2",
                    type=EdgeType.CUSTOM,
                    custom_predicate="MODULATES",
                ),
            ],
        )

        grounded = ground_knowledge_graph(kg)

        assert grounded.edges[0].match_status == MatchStatus.CUSTOM
        assert grounded.edges[0].custom_predicate == "MODULATES"


class TestGroundingIndexes:
    def test_iarc_index_prefers_canonical_name_for_bap_alias(self):
        index = build_iarc_grounding_index()
        match = index[normalize_grounding_key("BaP")]

        assert match.match_status == MatchStatus.ALIAS
        assert match.canonical_label == "Benzo[a]pyrene"
        assert match.canonical_id == "50-32-8"

    def test_default_index_contains_panel_and_iarc_entries(self):
        index = build_default_grounding_index()

        assert normalize_grounding_key("CYP1A1") in index
        assert normalize_grounding_key("Benzo[a]pyrene") in index

    def test_default_index_merges_extra_reference_graphs(self):
        extra_graph = KnowledgeGraph(
            nodes=[
                Node(
                    id="hsa05204",
                    label="Chemical Carcinogenesis",
                    type=NodeType.PATHWAY,
                    match_status=MatchStatus.CANONICAL,
                )
            ],
            edges=[],
        )

        index = build_default_grounding_index(reference_graphs=[("kegg", extra_graph)])

        assert normalize_grounding_key("CYP1A1") in index
        assert normalize_grounding_key("hsa05204") in index


class TestPrepareKnowledgeGraph:
    def test_exploratory_mode_keeps_unmatched_content(self):
        kg = KnowledgeGraph(
            nodes=[
                Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
                Node(id="n2", label="Unknown Chemical", type=NodeType.CARCINOGEN),
            ],
            edges=[Edge(source="n1", target="n2", type=EdgeType.ACTIVATES)],
        )

        prepared, warnings = prepare_knowledge_graph(kg, mode=GraphMode.EXPLORATORY)

        assert len(prepared.nodes) == 2
        assert prepared.nodes[1].match_status == MatchStatus.UNMATCHED
        assert prepared.edges[0].match_status == MatchStatus.UNMATCHED
        assert warnings == []

    def test_strict_mode_drops_unmatched_content(self):
        kg = KnowledgeGraph(
            nodes=[
                Node(id="n1", label="CYP1A1", type=NodeType.ENZYME),
                Node(id="n2", label="Unknown Chemical", type=NodeType.CARCINOGEN),
            ],
            edges=[Edge(source="n1", target="n2", type=EdgeType.ACTIVATES)],
        )

        prepared, warnings = prepare_knowledge_graph(kg, mode=GraphMode.STRICT)

        assert [node.id for node in prepared.nodes] == ["n1"]
        assert prepared.edges == []
        assert warnings
        assert "Strict mode dropped 1 non-canonical node" in warnings[0]
        assert "Strict mode dropped 1 non-canonical edge" in warnings[1]
