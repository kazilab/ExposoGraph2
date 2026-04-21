"""Tests for ExposoGraph.graph_filters."""

from ExposoGraph.config import GraphVisibility
from ExposoGraph.engine import GraphEngine
from ExposoGraph.graph_filters import (
    filter_knowledge_graph,
    filtered_engine,
    graph_visibility_label,
)
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, MatchStatus, Node, NodeType


def _sample_graph() -> KnowledgeGraph:
    return KnowledgeGraph(
        nodes=[
            Node(
                id="CYP1A1",
                label="CYP1A1",
                type=NodeType.ENZYME,
                match_status=MatchStatus.CANONICAL,
            ),
            Node(
                id="BaP",
                label="BaP",
                type=NodeType.CARCINOGEN,
                match_status=MatchStatus.ALIAS,
                canonical_id="50-32-8",
                canonical_label="Benzo[a]pyrene",
            ),
            Node(
                id="Mystery",
                label="Mystery Chemical",
                type=NodeType.CARCINOGEN,
                match_status=MatchStatus.UNMATCHED,
            ),
            Node(
                id="Novel",
                label="Novel Exposure",
                type=NodeType.CARCINOGEN,
                match_status=MatchStatus.CUSTOM,
                custom_type="Exposure",
            ),
        ],
        edges=[
            Edge(
                source="CYP1A1",
                target="BaP",
                type=EdgeType.ACTIVATES,
                match_status=MatchStatus.CANONICAL,
            ),
            Edge(
                source="CYP1A1",
                target="Mystery",
                type=EdgeType.ACTIVATES,
                match_status=MatchStatus.UNMATCHED,
            ),
            Edge(
                source="Mystery",
                target="Novel",
                type=EdgeType.CUSTOM,
                match_status=MatchStatus.CUSTOM,
                custom_predicate="CO_OCCURS_WITH",
            ),
        ],
    )


class TestGraphVisibilityLabel:
    def test_returns_human_readable_labels(self):
        assert graph_visibility_label(GraphVisibility.ALL) == "All"
        assert graph_visibility_label(GraphVisibility.VALIDATED_ONLY) == "Validated Only"
        assert graph_visibility_label(GraphVisibility.EXPLORATORY_ONLY) == "Exploratory Only"


class TestFilterKnowledgeGraph:
    def test_all_visibility_keeps_entire_graph(self):
        graph = filter_knowledge_graph(_sample_graph(), GraphVisibility.ALL)

        assert len(graph.nodes) == 4
        assert len(graph.edges) == 3

    def test_validated_visibility_keeps_only_validated_content(self):
        graph = filter_knowledge_graph(_sample_graph(), GraphVisibility.VALIDATED_ONLY)

        assert {node.id for node in graph.nodes} == {"CYP1A1", "BaP"}
        assert len(graph.edges) == 1
        assert graph.edges[0].match_status == MatchStatus.CANONICAL

    def test_exploratory_visibility_keeps_only_provisional_content(self):
        graph = filter_knowledge_graph(_sample_graph(), GraphVisibility.EXPLORATORY_ONLY)

        assert {node.id for node in graph.nodes} == {"Mystery", "Novel"}
        assert len(graph.edges) == 1
        assert graph.edges[0].type == EdgeType.CUSTOM
        assert graph.edges[0].match_status == MatchStatus.CUSTOM

    def test_edges_with_hidden_carcinogen_context_are_dropped(self):
        graph = KnowledgeGraph(
            nodes=[
                Node(
                    id="A",
                    label="CYP1A1",
                    type=NodeType.ENZYME,
                    match_status=MatchStatus.CANONICAL,
                ),
                Node(
                    id="B",
                    label="BaP",
                    type=NodeType.CARCINOGEN,
                    match_status=MatchStatus.CANONICAL,
                ),
                Node(
                    id="C",
                    label="Unknown",
                    type=NodeType.CARCINOGEN,
                    match_status=MatchStatus.UNMATCHED,
                ),
            ],
            edges=[
                Edge(
                    source="A",
                    target="B",
                    type=EdgeType.ACTIVATES,
                    carcinogen="C",
                    match_status=MatchStatus.CANONICAL,
                )
            ],
        )

        filtered = filter_knowledge_graph(graph, GraphVisibility.VALIDATED_ONLY)

        assert len(filtered.nodes) == 2
        assert filtered.edges == []

    def test_returns_independent_model_copies(self):
        source = _sample_graph()

        filtered = filter_knowledge_graph(source, GraphVisibility.ALL)
        filtered.nodes[0].label = "Changed"
        filtered.edges[0].label = "Updated"

        assert source.nodes[0].label == "CYP1A1"
        assert source.edges[0].label is None


class TestFilteredEngine:
    def test_returns_engine_for_selected_visibility(self):
        engine = GraphEngine()
        for node in _sample_graph().nodes:
            engine.add_node(node)
        for edge in _sample_graph().edges:
            engine.add_edge(edge)

        filtered = filtered_engine(engine, GraphVisibility.VALIDATED_ONLY)

        assert filtered.node_count == 2
        assert filtered.edge_count == 1
