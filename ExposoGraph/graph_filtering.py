"""Source-level graph filtering helpers for public graph/export paths."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from .models import Edge, KnowledgeGraph, Node, NodeType


def _norm(value: object) -> str:
    return str(value or "").strip().lower().replace("_", " ").replace("-", " ")


def _norm_set(values: Iterable[object] | None) -> frozenset[str]:
    return frozenset(_norm(value) for value in values or () if _norm(value))


@dataclass(frozen=True)
class GraphFilterCriteria:
    """Intersection criteria for local graph filtering."""

    carcinogen_classes: frozenset[str] = frozenset()
    tissues: frozenset[str] = frozenset()
    min_tissue_weight: float | None = None
    node_types: frozenset[str] = frozenset()

    @classmethod
    def from_values(
        cls,
        *,
        carcinogen_classes: Iterable[str] | None = None,
        tissues: Iterable[str] | None = None,
        min_tissue_weight: float | None = None,
        node_types: Iterable[str] | None = None,
    ) -> "GraphFilterCriteria":
        return cls(
            carcinogen_classes=_norm_set(carcinogen_classes),
            tissues=_norm_set(tissues),
            min_tissue_weight=min_tissue_weight,
            node_types=_norm_set(node_types),
        )


_HEAVY_METAL_ALIASES = frozenset(
    {
        "heavy metal",
        "heavy metals",
        "heavymetal",
        "metal",
        "metals",
        "metallograph",
        "metallographic",
    }
)


def _node_text_values(node: Node) -> set[str]:
    return {
        _norm(node.id),
        _norm(node.label),
        _norm(node.group),
        _norm(node.iarc),
        _norm(node.type.value),
        _norm(node.detail),
    }


def _matches_class(node: Node, classes: frozenset[str]) -> bool:
    if not classes:
        return True
    values = _node_text_values(node)
    if classes & values:
        return True
    text = " ".join(values)
    if classes & _HEAVY_METAL_ALIASES and any(alias in text for alias in _HEAVY_METAL_ALIASES):
        return True
    return any(criteria in text for criteria in classes)


def _matches_node_type(node: Node, node_types: frozenset[str]) -> bool:
    if not node_types:
        return True
    return _norm(node.type.value) in node_types or _norm(node.custom_type) in node_types


def _node_tissue_weight(node: Node, tissue: str) -> float | None:
    if not node.tissue_weights:
        return None
    for key, value in node.tissue_weights.items():
        if _norm(key) == tissue:
            try:
                return float(value)
            except (TypeError, ValueError):
                return None
    return None


def _matches_tissue(node: Node, criteria: GraphFilterCriteria) -> bool:
    if not criteria.tissues and criteria.min_tissue_weight is None:
        return True
    tissues = criteria.tissues or frozenset(_norm(key) for key in (node.tissue_weights or {}))
    if not tissues:
        return False
    node_tissue = _norm(node.tissue)
    for tissue in tissues:
        if node_tissue and tissue in node_tissue:
            if criteria.min_tissue_weight is None:
                return True
        weight = _node_tissue_weight(node, tissue)
        if weight is None:
            continue
        if criteria.min_tissue_weight is None or weight >= criteria.min_tissue_weight:
            return True
    return False


def _edge_matches(edge: Edge, kept_node_ids: set[str], criteria: GraphFilterCriteria) -> bool:
    if edge.source not in kept_node_ids or edge.target not in kept_node_ids:
        return False
    if edge.carcinogen is not None and edge.carcinogen not in kept_node_ids:
        return False
    if criteria.tissues and edge.tissue:
        return any(tissue in _norm(edge.tissue) for tissue in criteria.tissues)
    return True


def _ids_matching_classes(graph: KnowledgeGraph, classes: frozenset[str]) -> set[str] | None:
    if not classes:
        return None
    class_node_ids = {node.id for node in graph.nodes if _matches_class(node, classes)}
    allowed = set(class_node_ids)
    for edge in graph.edges:
        if (
            edge.carcinogen in class_node_ids
            or edge.source in class_node_ids
            or edge.target in class_node_ids
        ):
            allowed.update({edge.source, edge.target})
            if edge.carcinogen:
                allowed.add(edge.carcinogen)
    return allowed


def _ids_matching_tissues(graph: KnowledgeGraph, criteria: GraphFilterCriteria) -> set[str] | None:
    if not criteria.tissues and criteria.min_tissue_weight is None:
        return None
    allowed = {node.id for node in graph.nodes if _matches_tissue(node, criteria)}
    if criteria.tissues and criteria.min_tissue_weight is None:
        for edge in graph.edges:
            if edge.tissue and any(tissue in _norm(edge.tissue) for tissue in criteria.tissues):
                allowed.update({edge.source, edge.target})
                if edge.carcinogen:
                    allowed.add(edge.carcinogen)
    return allowed


def filter_graph_by_criteria(
    graph: KnowledgeGraph,
    criteria: GraphFilterCriteria | None = None,
    **criteria_values: Any,
) -> KnowledgeGraph:
    """Return a detached graph filtered by the intersection of all criteria."""
    resolved = criteria or GraphFilterCriteria.from_values(**criteria_values)
    class_allowed_ids = _ids_matching_classes(graph, resolved.carcinogen_classes)
    tissue_allowed_ids = _ids_matching_tissues(graph, resolved)
    nodes = [
        node.model_copy(deep=True)
        for node in graph.nodes
        if (class_allowed_ids is None or node.id in class_allowed_ids)
        and (tissue_allowed_ids is None or node.id in tissue_allowed_ids)
        and _matches_node_type(node, resolved.node_types)
    ]
    kept_node_ids = {node.id for node in nodes}
    edges = [
        edge.model_copy(deep=True)
        for edge in graph.edges
        if _edge_matches(edge, kept_node_ids, resolved)
    ]
    return KnowledgeGraph(nodes=nodes, edges=edges)


def graph_filter_to_json_safe(graph: KnowledgeGraph) -> dict[str, Any]:
    """Serialize a filtered graph to a JSON-safe dictionary."""
    return graph.model_dump(mode="json")


def heavy_metal_node_ids(graph: KnowledgeGraph) -> set[str]:
    """Return IDs for heavy-metal/metallograph nodes present in a graph."""
    return {
        node.id
        for node in graph.nodes
        if node.type is NodeType.CARCINOGEN and _matches_class(node, _HEAVY_METAL_ALIASES)
    }
