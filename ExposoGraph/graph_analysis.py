"""Graph analysis functions built on top of the NetworkX engine.

All functions are pure — they read from a :class:`GraphEngine` but never
mutate it.  This module provides domain-aware queries (metabolism chains,
variant impact) alongside standard graph-theory algorithms (shortest path,
centrality).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import networkx as nx

from .engine import GraphEngine

# ── Types ─────────────────────────────────────────────────────────────────

# Edge types that form the core activation → adduct → repair chain.
_METABOLISM_EDGE_TYPES = frozenset({
    "ACTIVATES",
    "DETOXIFIES",
    "TRANSPORTS",
    "FORMS_ADDUCT",
    "REPAIRS",
})


@dataclass
class MetabolismChain:
    """Result of :func:`metabolism_chain`."""

    carcinogen_id: str
    node_ids: list[str] = field(default_factory=list)
    edges: list[dict[str, Any]] = field(default_factory=list)

    @property
    def activation_edges(self) -> list[dict[str, Any]]:
        return [e for e in self.edges if e.get("type") == "ACTIVATES"]

    @property
    def detox_edges(self) -> list[dict[str, Any]]:
        return [e for e in self.edges if e.get("type") == "DETOXIFIES"]

    @property
    def adduct_edges(self) -> list[dict[str, Any]]:
        return [e for e in self.edges if e.get("type") == "FORMS_ADDUCT"]

    @property
    def repair_edges(self) -> list[dict[str, Any]]:
        return [e for e in self.edges if e.get("type") == "REPAIRS"]


@dataclass
class VariantImpact:
    """Result of :func:`variant_impact_score`."""

    gene_id: str
    activity_score: float | None
    downstream_adduct_count: int
    downstream_repair_count: int
    score: float


# ── Path & centrality ────────────────────────────────────────────────────


def shortest_path(
    engine: GraphEngine,
    source: str,
    target: str,
) -> list[str] | None:
    """Return the shortest path between *source* and *target*, or ``None``.

    Operates on the undirected view of the graph so that both incoming
    and outgoing edges are considered.
    """
    G = engine.G.to_undirected(as_view=True)
    try:
        return list(nx.shortest_path(G, source, target))
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return None


def all_shortest_paths(
    engine: GraphEngine,
    source: str,
    target: str,
) -> list[list[str]]:
    """Return *all* shortest paths between *source* and *target*."""
    G = engine.G.to_undirected(as_view=True)
    try:
        return [list(p) for p in nx.all_shortest_paths(G, source, target)]
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return []


def centrality(
    engine: GraphEngine,
    method: str = "degree",
) -> dict[str, float]:
    """Compute centrality scores for all nodes.

    *method* must be one of ``"degree"``, ``"betweenness"``, or
    ``"closeness"``.
    """
    G = engine.G
    if method == "degree":
        return dict(nx.degree_centrality(G))
    if method == "betweenness":
        return dict(nx.betweenness_centrality(G))
    if method == "closeness":
        return dict(nx.closeness_centrality(G))
    raise ValueError(
        f"Unknown centrality method: {method!r} (use degree, betweenness, or closeness)"
    )


# ── Domain-specific queries ───────────────────────────────────────────────


def metabolism_chain(
    engine: GraphEngine,
    carcinogen_id: str,
) -> MetabolismChain:
    """Extract the full metabolic chain for a carcinogen.

    Traverses edges of type ACTIVATES, DETOXIFIES, TRANSPORTS,
    FORMS_ADDUCT, and REPAIRS that are annotated with the given
    *carcinogen_id* (via the ``carcinogen`` edge attribute) **or** that
    directly connect to the carcinogen node.
    """
    chain = MetabolismChain(carcinogen_id=carcinogen_id)

    if carcinogen_id not in engine.G:
        return chain

    visited_nodes: set[str] = set()
    chain_edges: list[dict[str, Any]] = []

    for u, v, data in engine.G.edges(data=True):
        etype = data.get("type", "")
        if etype not in _METABOLISM_EDGE_TYPES:
            continue
        edge_carcinogen = data.get("carcinogen")
        if (
            edge_carcinogen == carcinogen_id
            or u == carcinogen_id
            or v == carcinogen_id
        ):
            chain_edges.append(dict(data, _source=u, _target=v))
            visited_nodes.update((u, v))

    # BFS to include transitive edges carrying the same carcinogen annotation,
    # so we follow e.g. metabolite → adduct → repair along the carcinogen's
    # pathway without absorbing unlabeled branches that merely share an enzyme.
    changed = True
    while changed:
        changed = False
        for u, v, data in engine.G.edges(data=True):
            etype = data.get("type", "")
            if etype not in _METABOLISM_EDGE_TYPES:
                continue
            if data.get("carcinogen") != carcinogen_id:
                continue
            edge_dict = dict(data, _source=u, _target=v)
            if edge_dict in chain_edges:
                continue
            if u in visited_nodes or v in visited_nodes:
                chain_edges.append(edge_dict)
                new_nodes = {u, v} - visited_nodes
                if new_nodes:
                    visited_nodes.update(new_nodes)
                    changed = True

    chain.node_ids = sorted(visited_nodes)
    chain.edges = [{k: v for k, v in e.items() if not k.startswith("_")} for e in chain_edges]
    return chain


def pathway_subgraph(
    engine: GraphEngine,
    pathway_id: str,
) -> list[str]:
    """Return node IDs connected to *pathway_id* via PATHWAY edges."""
    if pathway_id not in engine.G:
        return []

    connected: set[str] = {pathway_id}
    for u, v, data in engine.G.edges(data=True):
        if data.get("type") != "PATHWAY":
            continue
        if u == pathway_id:
            connected.add(v)
        elif v == pathway_id:
            connected.add(u)
    return sorted(connected)


def variant_impact_score(
    engine: GraphEngine,
    gene_id: str,
) -> VariantImpact | None:
    """Compute a variant impact score for a gene node.

    The score combines the node's ``activity_score`` with the number of
    downstream adduct-forming and repair paths reachable from the gene.
    A gene with a low activity score that sits upstream of many adducts
    receives a higher impact score (higher = more impactful variant).

    Returns ``None`` if *gene_id* is not in the graph.
    """
    node_data = engine.get_data(gene_id)
    if node_data is None:
        return None

    activity = node_data.get("activity_score")

    # Count downstream adduct and repair edges reachable via directed traversal
    reachable: set[str] = set()
    queue = [gene_id]
    while queue:
        current = queue.pop()
        for successor in engine.G.successors(current):
            if successor not in reachable:
                reachable.add(successor)
                queue.append(successor)

    adduct_count = 0
    repair_count = 0
    for u, v, data in engine.G.edges(data=True):
        if u not in reachable and u != gene_id:
            continue
        etype = data.get("type", "")
        if etype == "FORMS_ADDUCT":
            adduct_count += 1
        elif etype == "REPAIRS":
            repair_count += 1

    # Score: topology component (adducts - repairs) weighted by inverse
    # activity.  A normal metabolizer (activity=1.0) with 2 downstream
    # adducts and 1 repair → score = (2−1) × (2−1.0) = 1.0
    # A poor metabolizer (activity=0.0) → score = (2−1) × (2−0.0) = 2.0
    topology = max(adduct_count - repair_count, 0)
    activity_weight = 2.0 - (activity if activity is not None else 1.0)
    score = round(topology * activity_weight, 3)

    return VariantImpact(
        gene_id=gene_id,
        activity_score=activity,
        downstream_adduct_count=adduct_count,
        downstream_repair_count=repair_count,
        score=score,
    )
