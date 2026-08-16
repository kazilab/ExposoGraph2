"""NetworkX-backed graph engine for building and querying the knowledge graph."""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping, Sequence
from typing import Any

import networkx as nx

from .config import GraphMode
from .grounding import prepare_knowledge_graph
from .models import Edge, KnowledgeGraph, Node

logger = logging.getLogger(__name__)

_MISSING = object()
"""Sentinel distinguishing "path segment absent" from a stored ``None`` value."""


def _resolve_path(data: Mapping[str, Any], key: str | Sequence[str]) -> Any:
    """Drill into a node/edge attribute mapping along *key*.

    ``key`` may be:

    - a plain attribute name, e.g. ``"tissue_weights"``
    - a dot-delimited path into nested attributes, e.g.
      ``"tissue_weights.Liver"`` or ``"kinetics.Km_uM"``
    - an explicit sequence of path segments, e.g. ``("Ki", "NDMA")`` --
      useful when a segment name might itself contain a literal ``"."``

    Returns the sentinel :data:`_MISSING` if any segment along the path is
    absent, so callers can distinguish "not found" from a stored ``None``.
    """
    segments = key.split(".") if isinstance(key, str) else list(key)
    current: Any = data
    for segment in segments:
        if isinstance(current, Mapping) and segment in current:
            current = current[segment]
        else:
            return _MISSING
    return current


class GraphEngine:
    """Thin wrapper around a NetworkX MultiDiGraph that speaks our domain model."""

    def __init__(self) -> None:
        self.G: nx.MultiDiGraph = nx.MultiDiGraph()

    # ── Mutations ────────────────────────────────────────────────────────

    def add_node(self, node: Node) -> None:
        self.G.add_node(node.id, **node.model_dump(exclude_none=True, mode="json"))

    def add_edge(self, edge: Edge) -> None:
        if edge.source not in self.G:
            raise ValueError(f"Missing source node: {edge.source}")
        if edge.target not in self.G:
            raise ValueError(f"Missing target node: {edge.target}")
        if edge.carcinogen and edge.carcinogen not in self.G:
            raise ValueError(f"Missing carcinogen context node: {edge.carcinogen}")

        self.G.add_edge(
            edge.source,
            edge.target,
            **edge.model_dump(exclude_none=True, mode="json"),
        )

    def remove_node(self, node_id: str) -> None:
        if node_id in self.G:
            self.G.remove_node(node_id)

    def remove_edge(self, source: str, target: str, key: str | None = None) -> None:
        if key is not None:
            if self.G.has_edge(source, target, key):
                self.G.remove_edge(source, target, key)
            return
        if self.G.has_edge(source, target):
            self.G.remove_edge(source, target)

    # ── Bulk operations ──────────────────────────────────────────────────

    def _validated_reference_graph(self) -> KnowledgeGraph | None:
        if self.node_count == 0:
            return None
        current_graph = self.to_knowledge_graph()
        validated_graph, _warnings = prepare_knowledge_graph(
            current_graph,
            mode=GraphMode.STRICT,
        )
        if not validated_graph.nodes:
            return None
        return validated_graph

    def load(
        self,
        kg: KnowledgeGraph,
        *,
        mode: GraphMode | str = GraphMode.EXPLORATORY,
    ) -> list[str]:
        """Replace the current graph with *kg*.

        Clears all existing nodes and edges before loading.
        Returns a list of warning messages for any skipped edges.
        """
        self.clear()
        return self.merge(kg, mode=mode)

    def merge(
        self,
        kg: KnowledgeGraph,
        *,
        mode: GraphMode | str = GraphMode.EXPLORATORY,
    ) -> list[str]:
        """Additive merge — new nodes/edges are added, existing ones updated.

        Returns a list of warning messages for any skipped edges.
        """
        reference_graphs: list[tuple[str, KnowledgeGraph]] = []
        validated_graph = self._validated_reference_graph()
        if validated_graph is not None:
            reference_graphs.append(("current_graph", validated_graph))

        prepared_graph, warnings = prepare_knowledge_graph(
            kg,
            mode=mode,
            reference_graphs=reference_graphs or None,
        )
        for node in prepared_graph.nodes:
            self.add_node(node)
        for edge in prepared_graph.edges:
            try:
                self.add_edge(edge)
            except ValueError as exc:
                warnings.append(str(exc))
                logger.warning("Skipped edge during merge: %s", exc)
        return warnings

    def clear(self) -> None:
        self.G.clear()

    # ── Queries ──────────────────────────────────────────────────────────

    @property
    def node_count(self) -> int:
        return int(self.G.number_of_nodes())

    @property
    def edge_count(self) -> int:
        return int(self.G.number_of_edges())

    def get_node(
        self,
        node_id: str,
        key: str | Sequence[str] | None = None,
        *,
        default: Any = None,
    ) -> Any:
        """Return a node's attributes, or one (possibly nested) value.

        With no *key*, returns the node's full attribute dict, or ``None``
        if *node_id* isn't in the graph -- unchanged from prior behavior.

        With *key* set, drills into the node's attributes the same way
        :meth:`get_edge` drills into an edge's, e.g.::

            engine.get_node("CYP1A1", "tissue_weights.Liver")
            engine.get_node("CYP1A1", ("tissue_weights", "Liver"))

        Returns *default* (``None`` unless overridden) if the node is
        missing or any segment of *key* isn't present.
        """
        if node_id not in self.G:
            return None if key is None else default
        data = dict(self.G.nodes[node_id])
        if key is None:
            return data
        resolved = _resolve_path(data, key)
        return default if resolved is _MISSING else resolved

    def get_edge_keys(self, source: str, target: str) -> list[Any]:
        """Return the parallel-edge keys between *source* and *target*.

        Empty if the two nodes have no edge. Most edges in this graph are
        singular, in which case this returns a single-item list.
        """
        if not self.G.has_edge(source, target):
            return []
        return list(self.G[source][target].keys())

    def get_edge(
        self,
        source: str,
        target: str,
        key: str | Sequence[str] | None = None,
        *,
        edge_key: Any | None = None,
        default: Any = None,
    ) -> Any:
        """Return an edge's attributes, or one (possibly nested) value.

        Mirrors :meth:`get_node`'s syntax. With no *key*, returns the
        edge's full attribute dict, or ``None`` if *source*/*target*
        aren't connected. With *key* set, drills into that dict the same
        way -- e.g. an enzyme-substrate-specific inhibition constant::

            engine.get_edge("NDMA", "CYP2E1", "Ki.acetaldehyde")
            engine.get_edge("NDMA", "CYP2E1", ("Ki", "acetaldehyde"))

        If *source*/*target* have more than one parallel edge, pass
        *edge_key* to pick a specific one (see :meth:`get_edge_keys`);
        otherwise the first parallel edge is used. Returns *default*
        (``None`` unless overridden) if the edge is missing, *edge_key*
        doesn't match an existing parallel edge, or any segment of *key*
        isn't present.
        """
        if not self.G.has_edge(source, target):
            return None if key is None else default
        edge_view = self.G[source][target]
        if edge_key is not None:
            if edge_key not in edge_view:
                return None if key is None else default
            data = dict(edge_view[edge_key])
        else:
            data = dict(next(iter(edge_view.values())))
        if key is None:
            return data
        resolved = _resolve_path(data, key)
        return default if resolved is _MISSING else resolved

    def get_data(
        self,
        source: str,
        target: str | None = None,
        *,
        key: str | Sequence[str] | None = None,
        edge_key: Any | None = None,
        default: Any = None,
    ) -> Any:
        """Generic node/edge lookup that routes to :meth:`get_node` or
        :meth:`get_edge` based on whether *target* is given.

        - ``target`` omitted (``None``) -> routes to :meth:`get_node`,
          treating *source* as a node id::

              engine.get_data("CYP1A1")                            # full node dict
              engine.get_data("CYP1A1", key="tissue_weights.Liver")  # nested value

        - ``target`` given -> routes to :meth:`get_edge`, treating
          *source*/*target* as an edge's endpoints::

              engine.get_data("NDMA", "CYP2E1")  # full edge dict
              engine.get_data("NDMA", "CYP2E1", key="kinetics.Ki.ethanol")  # nested value
              # disambiguate parallel edges:
              engine.get_data("A", "B", key="pmid", edge_key=some_key)

        *key* and *edge_key* are keyword-only here (unlike on
        :meth:`get_node`/:meth:`get_edge`, where *key* is positional) so
        that the second positional argument is never ambiguous between
        "this is a nested key" and "this is the edge's target node".
        """
        if target is None:
            return self.get_node(source, key, default=default)
        return self.get_edge(source, target, key, edge_key=edge_key, default=default)

    def neighbors(self, node_id: str) -> list[str]:
        if node_id not in self.G:
            return []
        return list(self.G.successors(node_id)) + list(self.G.predecessors(node_id))

    def nodes_by_type(self, node_type: str) -> list[dict[str, Any]]:
        return [
            data
            for _, data in self.G.nodes(data=True)
            if data.get("type") == node_type
        ]

    # ── Serialization ────────────────────────────────────────────────────

    def to_dict(self) -> dict[str, list[Any]]:
        nodes = [dict(data) for _, data in self.G.nodes(data=True)]
        edges = [dict(data) for _, _, _, data in self.G.edges(keys=True, data=True)]
        return {"nodes": nodes, "edges": edges}

    def to_knowledge_graph(self) -> KnowledgeGraph:
        data = self.to_dict()
        return KnowledgeGraph(
            nodes=[Node(**n) for n in data["nodes"]],
            edges=[Edge(**e) for e in data["edges"]],
        )

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, default=str)

    # ── Validation ───────────────────────────────────────────────────────

    def validate(self) -> list[str]:
        errors: list[str] = []
        node_ids = set(self.G.nodes)
        for u, v, data in self.G.edges(data=True):
            if u not in node_ids:
                errors.append(f"Edge references missing source node: {u}")
            if v not in node_ids:
                errors.append(f"Edge references missing target node: {v}")
            if data.get("carcinogen") and data["carcinogen"] not in node_ids:
                errors.append(
                    f"Edge '{u}→{v}' references carcinogen '{data['carcinogen']}' "
                    f"which is not in the graph"
                )
        return errors
