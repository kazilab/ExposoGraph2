"""Helpers for filtering graphs into validated or exploratory views."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .config import GraphVisibility, normalize_graph_visibility
from .models import Edge, KnowledgeGraph, MatchStatus, Node

_VALIDATED_MATCH_STATUSES = frozenset({MatchStatus.CANONICAL, MatchStatus.ALIAS})

if TYPE_CHECKING:
    from .engine import GraphEngine


def graph_visibility_label(visibility: GraphVisibility | str) -> str:
    """Return a human-readable label for a graph visibility setting."""
    normalized = (
        visibility
        if isinstance(visibility, GraphVisibility)
        else normalize_graph_visibility(visibility)
    )
    labels = {
        GraphVisibility.ALL: "All",
        GraphVisibility.VALIDATED_ONLY: "Validated Only",
        GraphVisibility.EXPLORATORY_ONLY: "Exploratory Only",
    }
    return labels[normalized]


def filter_knowledge_graph(
    graph: KnowledgeGraph,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> KnowledgeGraph:
    """Return a filtered graph view without mutating the source graph."""
    normalized = (
        visibility
        if isinstance(visibility, GraphVisibility)
        else normalize_graph_visibility(visibility)
    )

    if normalized == GraphVisibility.ALL:
        return KnowledgeGraph(
            nodes=[node.model_copy(deep=True) for node in graph.nodes],
            edges=[edge.model_copy(deep=True) for edge in graph.edges],
        )

    if normalized == GraphVisibility.VALIDATED_ONLY:

        def keep_node(node: Node) -> bool:
            return node.match_status in _VALIDATED_MATCH_STATUSES

        def keep_edge(edge: Edge) -> bool:
            return edge.match_status in _VALIDATED_MATCH_STATUSES

    else:

        def keep_node(node: Node) -> bool:
            return node.match_status not in _VALIDATED_MATCH_STATUSES

        def keep_edge(edge: Edge) -> bool:
            return edge.match_status not in _VALIDATED_MATCH_STATUSES

    nodes: list[Node] = [node.model_copy(deep=True) for node in graph.nodes if keep_node(node)]
    node_ids = {node.id for node in nodes}
    edges: list[Edge] = [
        edge.model_copy(deep=True)
        for edge in graph.edges
        if keep_edge(edge)
        and edge.source in node_ids
        and edge.target in node_ids
        and (edge.carcinogen is None or edge.carcinogen in node_ids)
    ]
    return KnowledgeGraph(nodes=nodes, edges=edges)


def filtered_engine(
    engine: GraphEngine,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> GraphEngine:
    """Return a temporary GraphEngine with only the selected visible subgraph."""
    from .engine import GraphEngine

    filtered_graph = filter_knowledge_graph(engine.to_knowledge_graph(), visibility)
    filtered = GraphEngine()
    for node in filtered_graph.nodes:
        filtered.add_node(node)
    for edge in filtered_graph.edges:
        filtered.add_edge(edge)
    return filtered
