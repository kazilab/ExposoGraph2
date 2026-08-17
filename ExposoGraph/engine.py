"""NetworkX-backed graph engine for building and querying the knowledge graph."""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import networkx as nx

from .config import GraphMode
from .grounding import prepare_knowledge_graph
from .models import Edge, KnowledgeGraph, Node

_PACKAGE_DIR = Path(__file__).resolve().parent
_DEFAULT_GRAPH_DATA_PATH = _PACKAGE_DIR / "map" / "graph-data.json"
_DEFAULT_TISSUE_EXPRESSION_PATH = _PACKAGE_DIR / "data" / "tissue_expression_data_raw.json"
_DEFAULT_INTERACTION_PARAMETERS_PATH = _PACKAGE_DIR / "data" / "interaction_parameters.json"

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

    def load_reference_graph(
        self,
        *,
        graph_data_path: str | Path | None = None,
        tissue_expression_path: str | Path | None = None,
        interaction_parameters_path: str | Path | None = None,
    ) -> list[str]:
        """Load the bundled reference graph from ``map/graph-data.json``.

        This is the canonical way to instantiate the reference knowledge
        graph in Python (``graph-data.js`` remains the separate artifact
        consumed by the Streamlit/D3 viewer -- see ``exporter.to_graph_data_js``).

        After the base graph is loaded, two sources are (re)applied on top of
        it, in order:

        1. Tissue expression data from ``data/tissue_expression_data_raw.json`` is
           applied to the relevant enzyme nodes -- see
           :meth:`_apply_tissue_expression` for details. This *overwrites*
           whatever ``tissue_weights`` the bundled graph-data.json baked
           directly into those node attributes, so the freshly-sourced
           values become the sole source of truth.
        2. Interaction kinetics (both ``competitive_inhibition`` and
           ``phase2_conjugation`` blocks) from
           ``data/interaction_parameters.json`` are applied to the matching
           enzyme/substrate edges -- see :meth:`_apply_interaction_parameters`
           for details. Same overwrite semantics: the JSON file is the
           trusted source for ``Edge.kinetics``, not graph-data.json.

        Returns the combined warning messages from all three steps.
        """
        from .exporter import parse_graph_artifact  # local import avoids an import cycle

        resolved_graph_path = Path(graph_data_path) if graph_data_path else _DEFAULT_GRAPH_DATA_PATH
        kg = parse_graph_artifact(resolved_graph_path)
        warnings = self.load(kg)
        warnings.extend(self._apply_tissue_expression(tissue_expression_path))
        warnings.extend(self._apply_interaction_parameters(interaction_parameters_path))
        return warnings

    def _apply_tissue_expression(self, path: str | Path | None = None) -> list[str]:
        """(Re)apply ``tissue_expression_data_raw.json`` to the relevant enzyme nodes.

        The raw file covers 10 tissues (the original 8 plus
        ``Skin_NotSunExposed``/``Skin_SunExposed``) and 76 genes -- a superset
        of the older, pre-normalized ``tissue_expression_data.json`` (8
        tissues, 59 genes), which remains bundled only for the separate
        GTEx lookup helpers in ``tissue_subgraphs.py``, not for this method.

        For every ``Enzyme`` node with an entry in the source file's
        ``expression`` table, the node's attributes are set to:

        - ``tissue_weights_raw``: the raw per-tissue expression values,
          taken directly from the source file.
        - ``tissue_weights``: the same values normalized by dividing by
          the highest raw value across that enzyme's tissues (so the
          most-expressing tissue is always ``1.0``). This overwrites
          any ``tissue_weights`` the node already had (e.g. baked in by
          the bundled graph-data.json), which is no longer trusted as a
          data source once this method has run.

        Enzyme nodes with no entry in the source file are left with
        neither attribute (any pre-existing ``tissue_weights`` on them is
        also cleared, since it can no longer be attributed to this
        source of truth) and are reported as a warning.

        Returns a list of warning messages for enzyme nodes present in
        the graph but absent from the tissue expression source file.
        """
        resolved_path = Path(path) if path else _DEFAULT_TISSUE_EXPRESSION_PATH
        # This file's "expression" table is raw nTPM values with no
        # normalization applied -- the divide-by-max step below is this
        # method's own responsibility, unchanged regardless of source file.
        expression: dict[str, dict[str, float]] = json.loads(
            resolved_path.read_text(encoding="utf-8")
        )["expression"]

        warnings: list[str] = []
        enzyme_ids = [
            node_id for node_id, data in self.G.nodes(data=True) if data.get("type") == "Enzyme"
        ]
        for enzyme_id in enzyme_ids:
            node_data = self.G.nodes[enzyme_id]
            raw = expression.get(enzyme_id)
            if raw is None:
                node_data.pop("tissue_weights", None)
                node_data.pop("tissue_weights_raw", None)
                warnings.append(f"No tissue expression data for enzyme: {enzyme_id}")
                continue

            max_raw = max(raw.values()) if raw else 0.0
            normalized = (
                {tissue: value / max_raw for tissue, value in raw.items()}
                if max_raw
                else dict.fromkeys(raw, 0.0)
            )
            node_data["tissue_weights_raw"] = raw
            node_data["tissue_weights"] = normalized

        return warnings

    #: Sibling top-level blocks of ``interaction_parameters.json`` that share
    #: the identical ``<enzyme>.substrates.<substrate>`` shape and are both
    #: processed by :meth:`_apply_interaction_parameters`.
    _INTERACTION_PARAMETER_BLOCKS: tuple[str, ...] = (
        "competitive_inhibition",
        "phase2_conjugation",
    )

    def _apply_interaction_parameters(self, path: str | Path | None = None) -> list[str]:
        """(Re)apply enzyme/substrate kinetics from ``interaction_parameters.json``
        onto the matching edges.

        Two sibling top-level blocks share an identical shape and are both
        processed here, in this order: ``competitive_inhibition`` (Phase I
        bioactivation/detoxification competing for the same CYP active site)
        and ``phase2_conjugation`` (Phase II conjugation enzymes competing
        for the same transferase active site / co-substrate pool). See
        ``_INTERACTION_PARAMETER_BLOCKS``.

        For every ``<block>.<enzyme>.substrates.<substrate>`` entry in the
        source file, the entry's own ``graph_node_id`` field -- not the
        ``<substrate>`` JSON key itself -- names the graph node id it
        corresponds to. Every entry carries this field explicitly, including
        the majority where it's simply equal to the JSON key (an exact-match
        ``Carcinogen``, ``Metabolite``, or ``Substrate`` node id); naming
        mismatches (case differences, aliases such as
        ``trichloroethylene`` -> ``TCE``, or a substrate resolving to its
        parent carcinogen such as ``BPDE`` -> ``BaP``) are therefore resolved
        in the JSON itself, where a human reviewing the data can audit them,
        rather than through a runtime alias/exclusion table in this module.

        The entry's remaining fields (Km_uM, Vmax_relative, Ki_uM, product,
        product_carcinogenic, ...) are set, unchanged, as ``kinetics`` on
        every existing edge whose ``source`` is that enzyme and whose
        ``carcinogen`` attribute equals ``graph_node_id``. This *overwrites*
        any ``kinetics`` already present on those edges (e.g. baked in by the
        bundled graph-data.json), which is no longer trusted as a data source
        once this method has run -- the same overwrite semantics as
        :meth:`_apply_tissue_expression`. The two blocks' enzyme keys are
        disjoint (Phase I CYPs vs. Phase II transferases), so there is no
        cross-block collision risk in this overwrite step.

        Entries with a missing ``graph_node_id``, a ``graph_node_id`` that is
        not a node in the graph, or no matching edge are reported as warnings
        rather than raising, mirroring the tissue-expression warning pattern.

        Returns a list of warning messages.
        """
        resolved_path = Path(path) if path else _DEFAULT_INTERACTION_PARAMETERS_PATH
        source_data = json.loads(resolved_path.read_text(encoding="utf-8"))

        warnings: list[str] = []
        pending: dict[tuple[str, str], dict[str, Any]] = {}
        pending_block: dict[tuple[str, str], str] = {}
        for block_name in self._INTERACTION_PARAMETER_BLOCKS:
            block = source_data.get(block_name, {})
            for enzyme_id, enzyme_block in block.items():
                if enzyme_id == "_description" or not isinstance(enzyme_block, dict):
                    continue
                for substrate_key, params in enzyme_block.get("substrates", {}).items():
                    resolved = params.get("graph_node_id")
                    if not resolved:
                        warnings.append(
                            f"No graph_node_id declared for {block_name} "
                            f"substrate: {enzyme_id}/{substrate_key}"
                        )
                        continue
                    if resolved not in self.G:
                        warnings.append(
                            f"graph_node_id {resolved!r} for {block_name} substrate "
                            f"{enzyme_id}/{substrate_key} is not a node in the graph"
                        )
                        continue
                    kinetics = {k: v for k, v in params.items() if k != "graph_node_id"}
                    pending[(enzyme_id, resolved)] = kinetics
                    pending_block[(enzyme_id, resolved)] = block_name

        applied: set[tuple[str, str]] = set()
        for source_id, _target_id, edge_data in self.G.edges(data=True):
            key = (source_id, edge_data.get("carcinogen"))
            if key in pending:
                edge_data["kinetics"] = dict(pending[key])
                applied.add(key)

        for enzyme_id, resolved in pending:
            if (enzyme_id, resolved) not in applied:
                block_name = pending_block[(enzyme_id, resolved)]
                warnings.append(
                    f"No edge found for {block_name} pair: {enzyme_id} -> {resolved}"
                )

        return warnings

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

    # ── Domain-specific filters ──────────────────────────────────────────
    #
    # These read whatever tissue/group vocabulary is actually present in
    # the loaded graph (never a hardcoded list), and the subgraph-returning
    # methods all share one output shape -- {"nodes": [...], "edges": [...]}
    # of plain attribute dicts, matching :meth:`to_dict` -- so callers can
    # compose/pass them around uniformly. This is also the intended
    # foundation for eventually replacing the ad hoc per-filter JS in
    # ExposoGraph/map/index.html (applyCarcinogenFilter/applyTissueFilter)
    # with a single server/engine-side filtering path.

    def get_tissues(self) -> list[str]:
        """Return every tissue name present in any node's ``tissue_weights``.

        Currently only Enzyme nodes carry ``tissue_weights`` (populated by
        :meth:`_apply_tissue_expression` from ``tissue_expression_data_raw.json``
        at load time), but this scans all node types so it stays correct if
        that ever changes. Sorted, deduplicated; empty if no node has
        ``tissue_weights`` set.
        """
        tissues: set[str] = set()
        for _, data in self.G.nodes(data=True):
            weights = data.get("tissue_weights")
            if weights:
                tissues.update(weights.keys())
        return sorted(tissues)

    def get_carcinogen_groups(self) -> list[str]:
        """Return every distinct ``group`` label among Carcinogen nodes.

        e.g. ``"Aldehyde"``, ``"UV_Radiation"``, ``"PFAS"``. Scoped to
        ``type == "Carcinogen"`` specifically -- other node types (Enzyme,
        Gene, ...) use the same ``group`` attribute for unrelated groupings
        (e.g. DNA-repair pathway families), which this deliberately excludes.
        Sorted, deduplicated; empty if no Carcinogen node has a ``group``.
        """
        groups: set[str] = set()
        for _, data in self.G.nodes(data=True):
            if data.get("type") == "Carcinogen" and data.get("group"):
                groups.add(data["group"])
        return sorted(groups)

    def carcinogens_by_group(self, group: str) -> list[dict[str, Any]]:
        """Return every Carcinogen node whose ``group`` equals *group*.

        See :meth:`get_carcinogen_groups` for the available values. Empty
        list if *group* matches no Carcinogen node.
        """
        return [
            data
            for _, data in self.G.nodes(data=True)
            if data.get("type") == "Carcinogen" and data.get("group") == group
        ]

    def node_neighborhood(self, node_id: str) -> dict[str, list[Any]]:
        """Return the 1-hop neighborhood subgraph of *node_id*.

        Nodes: *node_id* itself plus every direct successor/predecessor
        (see :meth:`neighbors`). Edges: every edge directly incident to
        *node_id* (either direction, all parallel edges). Returns
        ``{"nodes": [], "edges": []}`` if *node_id* isn't in the graph.
        """
        if node_id not in self.G:
            return {"nodes": [], "edges": []}
        neighbor_ids = {node_id, *self.neighbors(node_id)}
        nodes = [dict(self.G.nodes[n]) for n in neighbor_ids]
        edges = [
            dict(data)
            for u, v, data in self.G.edges(data=True)
            if u == node_id or v == node_id
        ]
        return {"nodes": nodes, "edges": edges}

    def subgraph_by_node_type(self, node_type: str) -> dict[str, list[Any]]:
        """Return the node-only subgraph of every node with ``type == node_type``.

        Equivalent to ``{"nodes": nodes_by_type(node_type), "edges": []}``
        -- deliberately no edges, even if two matching nodes happen to be
        directly connected. See :meth:`subgraph_by_edge_type` for the
        edge-driven counterpart.
        """
        return {"nodes": self.nodes_by_type(node_type), "edges": []}

    def subgraph_by_edge_type(self, edge_type: str) -> dict[str, list[Any]]:
        """Return every edge with ``type == edge_type``, plus their adjacent nodes.

        Nodes are deduplicated across all matching edges' endpoints.
        """
        edges: list[dict[str, Any]] = []
        node_ids: set[str] = set()
        for u, v, data in self.G.edges(data=True):
            if data.get("type") == edge_type:
                edges.append(dict(data))
                node_ids.add(u)
                node_ids.add(v)
        nodes = [dict(self.G.nodes[n]) for n in node_ids]
        return {"nodes": nodes, "edges": edges}

    def enzymes_by_tissue_threshold(self, tissue: str, threshold: float) -> dict[str, list[Any]]:
        """Return Enzyme nodes at/above *threshold* for *tissue*, plus their neighborhood.

        Qualifying enzymes: ``type == "Enzyme"`` and
        ``tissue_weights[tissue] >= threshold`` (enzymes missing that
        tissue key never qualify, regardless of *threshold*'s sign).
        The result also includes every node/edge directly adjacent to a
        qualifying enzyme -- mirroring ExposoGraph/map/index.html's
        ``applyTissueFilter`` expansion -- not just the enzymes themselves.
        """
        qualifying_ids = {
            node_id
            for node_id, data in self.G.nodes(data=True)
            if data.get("type") == "Enzyme"
            and data.get("tissue_weights")
            and data["tissue_weights"].get(tissue, float("-inf")) >= threshold
        }
        node_ids = set(qualifying_ids)
        edges: list[dict[str, Any]] = []
        for u, v, data in self.G.edges(data=True):
            if u in qualifying_ids or v in qualifying_ids:
                edges.append(dict(data))
                node_ids.add(u)
                node_ids.add(v)
        nodes = [dict(self.G.nodes[n]) for n in node_ids]
        return {"nodes": nodes, "edges": edges}

    def filtered_subgraph(
        self,
        *,
        group: str | None = None,
        tissue: str | None = None,
        tissue_threshold: float | None = None,
        edge_type: str | None = None,
        node_type: str | None = None,
    ) -> dict[str, list[Any]]:
        """Return the subgraph at the intersection of up to four filter axes.

        Each axis is optional and left-``None`` axes impose no restriction
        (the identity element for the intersection); passing nothing
        returns the full graph. The four axes:

        - *group*: an edge must touch a Carcinogen node in this group,
          either directly (source/target) or via its ``carcinogen``
          field. See :meth:`get_carcinogen_groups`.
        - *tissue* (+ optional *tissue_threshold*, default ``-inf`` i.e.
          any value): an edge must touch an Enzyme node with
          ``tissue_weights[tissue] >= tissue_threshold``. See
          :meth:`get_tissues`. *tissue_threshold* without *tissue* raises
          ``ValueError``.
        - *edge_type*: the edge's own ``type`` must match.
        - *node_type*: an edge must touch a node with this ``type``.

        IMPORTANT -- this is a strict logical AND across whichever axes
        are given: an edge survives only if it satisfies *every* given
        axis simultaneously, evaluated per-edge. This is stricter than
        (and not a drop-in replacement for) ExposoGraph/map/index.html's
        current filter buttons, which apply each filter independently as
        an OR-based highlight rather than intersecting them -- e.g.
        requesting ``group="PAH"`` and ``node_type="Enzyme"`` here returns
        only PAH-carcinogen edges that *also* touch an Enzyme node, not
        "all PAH-related nodes" unioned with "all Enzyme nodes". The
        returned node set is the union of the surviving edges' endpoints
        only -- a carcinogen matching *group* with no edge satisfying the
        other axes will not appear. If that per-edge-AND semantics isn't
        what's wanted for a given caller, filter node lists from the other
        methods above directly instead.
        """
        if tissue_threshold is not None and tissue is None:
            raise ValueError("tissue_threshold requires tissue to also be given")

        group_node_ids: set[str] | None = None
        if group is not None:
            group_node_ids = {
                node_id
                for node_id, data in self.G.nodes(data=True)
                if data.get("type") == "Carcinogen" and data.get("group") == group
            }

        tissue_node_ids: set[str] | None = None
        if tissue is not None:
            threshold = tissue_threshold if tissue_threshold is not None else float("-inf")
            tissue_node_ids = {
                node_id
                for node_id, data in self.G.nodes(data=True)
                if data.get("type") == "Enzyme"
                and data.get("tissue_weights")
                and data["tissue_weights"].get(tissue, float("-inf")) >= threshold
            }

        type_node_ids: set[str] | None = None
        if node_type is not None:
            type_node_ids = {
                node_id for node_id, data in self.G.nodes(data=True) if data.get("type") == node_type
            }

        def _touches(u: str, v: str, carcinogen: str | None, candidates: set[str]) -> bool:
            return u in candidates or v in candidates or (carcinogen is not None and carcinogen in candidates)

        surviving_edges: list[dict[str, Any]] = []
        node_ids: set[str] = set()
        for u, v, data in self.G.edges(data=True):
            if edge_type is not None and data.get("type") != edge_type:
                continue
            carcinogen = data.get("carcinogen")
            if group_node_ids is not None and not _touches(u, v, carcinogen, group_node_ids):
                continue
            if tissue_node_ids is not None and not _touches(u, v, None, tissue_node_ids):
                continue
            if type_node_ids is not None and not (u in type_node_ids or v in type_node_ids):
                continue
            surviving_edges.append(dict(data))
            node_ids.add(u)
            node_ids.add(v)

        nodes = [dict(self.G.nodes[n]) for n in node_ids]
        return {"nodes": nodes, "edges": surviving_edges}

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
