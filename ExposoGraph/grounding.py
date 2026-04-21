"""Grounding helpers for reconciling nodes and graphs to canonical references.

This module adds a lightweight normalization layer on top of the existing
knowledge-graph schema. It does not rewrite node IDs; instead it annotates
nodes and edges with canonical metadata so later merge/validation stages can
decide whether to keep, normalize, or reject exploratory content.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

from .config import GraphMode, normalize_graph_mode
from .db_clients.iarc import IARCClassifier
from .models import Edge, EdgeType, KnowledgeGraph, MatchStatus, Node
from .reference_data import build_full_panel

_SPLIT_RE = re.compile(r"[/|;,+]")
_VALIDATED_MATCH_STATUSES = frozenset({MatchStatus.CANONICAL, MatchStatus.ALIAS})


@dataclass(frozen=True)
class GroundingMatch:
    """Canonical match metadata for a grounded term."""

    canonical_id: str
    canonical_label: str
    canonical_namespace: str
    match_status: MatchStatus
    extra_fields: dict[str, str | int | float | None] = field(default_factory=dict)


def normalize_grounding_key(value: str) -> str:
    """Normalize a label or identifier for fuzzy exact-match grounding."""
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _candidate_terms(*values: str | None) -> set[str]:
    terms: set[str] = set()
    for value in values:
        if value is None:
            continue
        stripped = value.strip()
        if not stripped:
            continue
        terms.add(stripped)
        for part in _SPLIT_RE.split(stripped):
            part = part.strip()
            if part:
                terms.add(part)
    return terms


def _canonical_match_status(term: str, *, canonical_id: str, canonical_label: str) -> MatchStatus:
    normalized = normalize_grounding_key(term)
    if normalized in {
        normalize_grounding_key(canonical_id),
        normalize_grounding_key(canonical_label),
    }:
        return MatchStatus.CANONICAL
    return MatchStatus.ALIAS


def _match_priority(status: MatchStatus) -> int:
    return 0 if status == MatchStatus.CANONICAL else 1


def _choose_canonical_label(names: Iterable[str]) -> str:
    return max(set(names), key=lambda name: (len(name), name))


def build_graph_grounding_index(
    graph: KnowledgeGraph,
    *,
    namespace: str,
) -> dict[str, GroundingMatch]:
    """Build a normalized lookup index from a reference KnowledgeGraph."""
    index: dict[str, GroundingMatch] = {}
    for node in graph.nodes:
        canonical_id = node.canonical_id or node.id
        canonical_label = node.canonical_label or node.label
        extra_fields: dict[str, str | int | float | None] = {
            "tier": node.tier,
            "group": node.group,
            "iarc": node.iarc,
        }
        for term in _candidate_terms(
            node.id,
            node.label,
            node.canonical_id,
            node.canonical_label,
        ):
            normalized = normalize_grounding_key(term)
            if not normalized:
                continue
            match = GroundingMatch(
                canonical_id=canonical_id,
                canonical_label=canonical_label,
                canonical_namespace=namespace,
                match_status=_canonical_match_status(
                    term,
                    canonical_id=canonical_id,
                    canonical_label=canonical_label,
                ),
                extra_fields=extra_fields,
            )
            existing = index.get(normalized)
            if existing is None or _match_priority(match.match_status) < _match_priority(
                existing.match_status
            ):
                index[normalized] = match
    return index


def build_iarc_grounding_index(
    classifier: IARCClassifier | None = None,
) -> dict[str, GroundingMatch]:
    """Build a normalized lookup index from the bundled IARC classifier data."""
    classifier = classifier or IARCClassifier()
    grouped: dict[str, tuple[list[str], dict[str, str]]] = {}
    for name in classifier.all_chemicals:
        entry = classifier.get_entry(name)
        if entry is None:
            continue
        group_key = entry.get("cas") or name
        bucket = grouped.setdefault(group_key, ([], entry))
        bucket[0].append(name)

    index: dict[str, GroundingMatch] = {}
    for group_key, (names, entry) in grouped.items():
        canonical_label = _choose_canonical_label(names)
        canonical_id = entry.get("cas") or canonical_label
        extra_fields: dict[str, str | int | float | None] = {
            "group": entry.get("category"),
            "iarc": entry.get("group"),
        }
        for term in _candidate_terms(*names, entry.get("cas")):
            normalized = normalize_grounding_key(term)
            if not normalized:
                continue
            match = GroundingMatch(
                canonical_id=canonical_id,
                canonical_label=canonical_label,
                canonical_namespace="iarc",
                match_status=_canonical_match_status(
                    term,
                    canonical_id=canonical_id,
                    canonical_label=canonical_label,
                ),
                extra_fields=extra_fields,
            )
            existing = index.get(normalized)
            if existing is None or _match_priority(match.match_status) < _match_priority(
                existing.match_status
            ):
                index[normalized] = match
    return index


def build_default_grounding_index(
    *,
    reference_graphs: Sequence[tuple[str, KnowledgeGraph]] | None = None,
    classifier: IARCClassifier | None = None,
) -> dict[str, GroundingMatch]:
    """Build the default grounding index from reference panels plus IARC."""
    graphs: list[tuple[str, KnowledgeGraph]] = [("reference_panel", build_full_panel())]
    if reference_graphs:
        graphs.extend(reference_graphs)
    index: dict[str, GroundingMatch] = {}
    for namespace, graph in graphs:
        for term, match in build_graph_grounding_index(graph, namespace=namespace).items():
            existing = index.get(term)
            if existing is None or _match_priority(match.match_status) < _match_priority(
                existing.match_status
            ):
                index[term] = match
    for term, match in build_iarc_grounding_index(classifier).items():
        existing = index.get(term)
        if existing is None or _match_priority(match.match_status) < _match_priority(
            existing.match_status
        ):
            index[term] = match
    return index


def ground_node(
    node: Node,
    *,
    grounding_index: dict[str, GroundingMatch] | None = None,
    reference_graphs: Sequence[tuple[str, KnowledgeGraph]] | None = None,
    classifier: IARCClassifier | None = None,
) -> Node:
    """Return a grounded copy of *node* using built-in or supplied references."""
    if node.match_status in {
        MatchStatus.CANONICAL,
        MatchStatus.ALIAS,
        MatchStatus.CUSTOM,
    }:
        return node

    grounding_index = grounding_index or build_default_grounding_index(
        reference_graphs=reference_graphs,
        classifier=classifier,
    )

    for term in (node.label, node.id):
        normalized = normalize_grounding_key(term)
        if not normalized:
            continue
        match = grounding_index.get(normalized)
        if match is None:
            continue
        updates: dict[str, Any] = {
            "match_status": match.match_status,
            "canonical_id": match.canonical_id,
            "canonical_label": match.canonical_label,
            "canonical_namespace": match.canonical_namespace,
        }
        if node.tier is None and match.extra_fields.get("tier") is not None:
            updates["tier"] = match.extra_fields["tier"]
        if node.group is None and match.extra_fields.get("group"):
            updates["group"] = match.extra_fields["group"]
        if node.iarc is None and match.extra_fields.get("iarc"):
            updates["iarc"] = match.extra_fields["iarc"]
        return node.model_copy(update=updates)

    return node.model_copy(update={"match_status": MatchStatus.UNMATCHED})


def ground_knowledge_graph(
    graph: KnowledgeGraph,
    *,
    grounding_index: dict[str, GroundingMatch] | None = None,
    reference_graphs: Sequence[tuple[str, KnowledgeGraph]] | None = None,
    classifier: IARCClassifier | None = None,
) -> KnowledgeGraph:
    """Ground all nodes in *graph* and derive edge grounding status."""
    grounding_index = grounding_index or build_default_grounding_index(
        reference_graphs=reference_graphs,
        classifier=classifier,
    )

    grounded_nodes = [
        ground_node(
            node,
            grounding_index=grounding_index,
            reference_graphs=reference_graphs,
            classifier=classifier,
        )
        for node in graph.nodes
    ]
    node_by_id = {node.id: node for node in grounded_nodes}

    grounded_edges: list[Edge] = []
    for edge in graph.edges:
        if edge.match_status in {MatchStatus.CANONICAL, MatchStatus.ALIAS, MatchStatus.CUSTOM}:
            grounded_edges.append(edge)
            continue
        if edge.type == EdgeType.CUSTOM:
            grounded_edges.append(edge.model_copy(update={"match_status": MatchStatus.CUSTOM}))
            continue
        source_node = node_by_id.get(edge.source)
        source_status = (
            source_node.match_status if source_node is not None else MatchStatus.UNMATCHED
        )
        target_node = node_by_id.get(edge.target)
        target_status = (
            target_node.match_status if target_node is not None else MatchStatus.UNMATCHED
        )
        if (
            source_status in _VALIDATED_MATCH_STATUSES
            and target_status in _VALIDATED_MATCH_STATUSES
        ):
            grounded_edges.append(
                edge.model_copy(
                    update={
                        "match_status": MatchStatus.CANONICAL,
                        "canonical_predicate": edge.type.value,
                        "canonical_namespace": "schema",
                    }
                )
            )
        else:
            grounded_edges.append(edge.model_copy(update={"match_status": MatchStatus.UNMATCHED}))

    return KnowledgeGraph(nodes=grounded_nodes, edges=grounded_edges)


def _format_drop_preview(values: Sequence[str], *, limit: int = 5) -> str:
    preview = list(values[:limit])
    suffix = "" if len(values) <= limit else ", ..."
    return ", ".join(preview) + suffix


def _strict_graph(graph: KnowledgeGraph) -> tuple[KnowledgeGraph, list[str]]:
    kept_nodes = [node for node in graph.nodes if node.match_status in _VALIDATED_MATCH_STATUSES]
    kept_node_ids = {node.id for node in kept_nodes}
    dropped_nodes = [node.id for node in graph.nodes if node.id not in kept_node_ids]

    kept_edges: list[Edge] = []
    dropped_unvalidated_edges: list[str] = []
    dropped_context_edges: list[str] = []
    for edge in graph.edges:
        edge_label = f"{edge.source}-{edge.type.value}->{edge.target}"
        if edge.match_status not in _VALIDATED_MATCH_STATUSES:
            dropped_unvalidated_edges.append(edge_label)
            continue
        if edge.source not in kept_node_ids or edge.target not in kept_node_ids:
            dropped_unvalidated_edges.append(edge_label)
            continue
        if edge.carcinogen and edge.carcinogen not in kept_node_ids:
            dropped_context_edges.append(edge_label)
            continue
        kept_edges.append(edge)

    warnings: list[str] = []
    if dropped_nodes:
        warnings.append(
            "Strict mode dropped "
            f"{len(dropped_nodes)} non-canonical node(s): {_format_drop_preview(dropped_nodes)}"
        )
    if dropped_unvalidated_edges:
        warnings.append(
            "Strict mode dropped "
            f"{len(dropped_unvalidated_edges)} non-canonical edge(s): "
            f"{_format_drop_preview(dropped_unvalidated_edges)}"
        )
    if dropped_context_edges:
        warnings.append(
            "Strict mode dropped "
            f"{len(dropped_context_edges)} edge(s) with non-canonical carcinogen context: "
            f"{_format_drop_preview(dropped_context_edges)}"
        )

    return KnowledgeGraph(nodes=kept_nodes, edges=kept_edges), warnings


def prepare_knowledge_graph(
    graph: KnowledgeGraph,
    *,
    mode: GraphMode | str = GraphMode.EXPLORATORY,
    grounding_index: dict[str, GroundingMatch] | None = None,
    reference_graphs: Sequence[tuple[str, KnowledgeGraph]] | None = None,
    classifier: IARCClassifier | None = None,
) -> tuple[KnowledgeGraph, list[str]]:
    """Ground *graph* and optionally filter it for strict-mode workflows."""
    normalized_mode = mode if isinstance(mode, GraphMode) else normalize_graph_mode(mode)
    grounded = ground_knowledge_graph(
        graph,
        grounding_index=grounding_index,
        reference_graphs=reference_graphs,
        classifier=classifier,
    )
    if normalized_mode == GraphMode.EXPLORATORY:
        return grounded, []
    return _strict_graph(grounded)
