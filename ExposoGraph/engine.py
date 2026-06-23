"""NetworkX-backed graph engine for building and querying the knowledge graph."""

from __future__ import annotations

import json
import logging
from typing import Any
from pathlib import Path
import networkx as nx
import re
from .config import GraphMode
from .grounding import prepare_knowledge_graph
from .models import (
    Edge,
    KnowledgeGraph,
    Node,
    GenotypeModifiersContainer,
    CompetitiveInhibitionSubstrate,
    Phase2ConjugationSubstrate,
    EnzymeInduction,
)

logger = logging.getLogger(__name__)


def load_from_json(f_name):
    with open(f_name, "r", encoding="utf-8") as file:
        return json.load(file)


class GraphEngine:
    """Thin wrapper around a NetworkX MultiDiGraph that speaks our domain model."""

    def __init__(self) -> None:
        self.G: nx.MultiDiGraph = nx.MultiDiGraph()

    # ── Mutations ────────────────────────────────────────────────────────
    def _impute_node(
        self, node_id: str, node_type: str = "unknown", visible: bool = True
    ) -> None:
        """Helper to initialize a missing node with default attributes."""
        if node_id not in self.G:
            logger.warning(f"Node_id {node_id} is missing. Imputing node into graph.")
            self.G.add_node(node_id, type=node_type, visible=visible)

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

    def modify_node(
        self,
        node_id: str,
        update_obj_name: str,
        update_obj_value,
        default_type: str = "unknown",
        default_visible: bool = True,
    ):
        # Pass type and visible configurations down to the imputation helper
        self._impute_node(node_id, node_type=default_type, visible=default_visible)
        self.G.nodes[node_id][update_obj_name] = update_obj_value

    def modify_edge(
        self,
        source_id: str,
        target_id: str,
        key: str,
        value,
        default_type: list | str = "unknown",
        default_visible: bool = True,
    ):
        logger.debug(
            f"modifying edge for source_id: {source_id}, target_id: {target_id}, key: {key} "
        )
        edge_needs_creation = 0
        s = source_id
        t = target_id

        # Parse type array/list defaults safely
        if isinstance(default_type, (list, tuple)) and len(default_type) >= 2:
            source_type = default_type[0]
            target_type = default_type[1]
        else:
            source_type = default_type
            target_type = default_type

        # 1. Check and impute endpoint nodes with structural index overrides
        if source_id not in self.G:
            logger.warning(f"source_id is not in G: {source_id}. Imputing node.")
            self._impute_node(source_id, node_type=source_type, visible=default_visible)

        if target_id not in self.G:
            logger.warning(f"target_id is not in G: {target_id}. Imputing node.")
            self._impute_node(target_id, node_type=target_type, visible=default_visible)
            t = target_id

        # 2. Check and flag missing edge
        if (source_id, target_id) not in self.G.edges and (
            target_id,
            source_id,
        ) not in self.G.edges:
            logger.warning("edge hasn't been created yet. Imputing edge.")
            edge_needs_creation = 1

        if (
            source_id,
            target_id,
        ) not in self.G.edges and edge_needs_creation == 0:  # edge exists, but wrong order
            s = target_id
            t = source_id

        # 3. Apply assignment
        if edge_needs_creation:
            self.G.add_edge(s, t)

        self.G[s][t][key] = value

    # ── Getters ──────────────────────────────────────────────────
    def get_node(self, node_id: str):
        return dict(self.G[node_id])

    def get_edge(self, source_id: str, target_id: str):
        if (source_id, target_id) in self.G.edges:
            return dict(self.G[source_id][target_id])
        if (target_id, source_id) in self.G.edges:
            return dict(self.G[target_id][source_id])

    # ── Filtering ──────────────────────────────────────────────────

    def filter_by_criteria(
        self,
        *,
        node_type: str | None = None,
        edge_type: str | None = None,
        carcinogen_class: str | None = None,
        tissue: str | None = None,
        min_tissue_weight: float | None = None,
    ) -> GraphEngine:
        """Filter the network graph across multiple criteria using subgraph intersections."""
        # Track all active subgraph overlays generated by the selected criteria
        active_subgraphs: list[nx.MultiDiGraph] = []

        # ── 1. Node Type Filtering ───────────────────────────────────────
        if node_type:
            matching_nodes = [
                n
                for n, d in self.G.nodes(data=True)
                if d.get("type") in list(node_type)  # == node_type
            ]
            # Subgraph containing specified nodes and any edges between them
            active_subgraphs.append(self.G.subgraph(matching_nodes).copy())

        # ── 2. Edge Type Filtering ───────────────────────────────────────
        if edge_type:
            matching_edges = [
                (u, v, k)
                for u, v, k, d in self.G.edges(keys=True, data=True)
                if d.get("type") == edge_type
            ]
            # Subgraph containing exclusively these edges and their attached endpoints
            active_subgraphs.append(self.G.edge_subgraph(matching_edges).copy())

        # ── 3. Carcinogen Class (Neighbors of Neighbors) ─────────────────
        if carcinogen_class:
            # Level 0: Core carcinogens matching the specified group field
            layer_0 = {
                n
                for n, d in self.G.nodes(data=True)
                if d.get("type") == "Carcinogen" and d.get("group") == carcinogen_class
            }

            # Level 1: Neighbors of those core carcinogens
            layer_1 = set(layer_0)
            for node_id in layer_0:
                layer_1.update(self.G.successors(node_id))
                layer_1.update(self.G.predecessors(node_id))

            # Level 2: Neighbors of neighbors
            layer_2 = set(layer_1)
            for node_id in layer_1:
                layer_2.update(self.G.successors(node_id))
                layer_2.update(self.G.predecessors(node_id))

            active_subgraphs.append(self.G.subgraph(layer_2).copy())

        # ── 4. Tissue & Weight Filtering ─────────────────────────────────
        if tissue:
            threshold = min_tissue_weight if min_tissue_weight is not None else 0.0

            # Find nodes where specified tissue field exists and is above threshold
            core_tissue_nodes = set()
            for n, d in self.G.nodes(data=True):
                # Assumes your tissue weight is stored in a dictionary attribute or directly
                tissue_data = d.get("tissue_weights")
                if (
                    isinstance(tissue_data, dict)
                    and tissue_data.get(tissue, 0.0) >= threshold
                ):
                    core_tissue_nodes.add(n)
                # elif d.get("tissue") == tissue and d.get("weight", 0.0) >= threshold:
                #     core_tissue_nodes.add(n)

            # Gather those core nodes and all of their immediate neighbors
            tissue_neighborhood = set(core_tissue_nodes)
            for node_id in core_tissue_nodes:
                tissue_neighborhood.update(self.G.successors(node_id))
                tissue_neighborhood.update(self.G.predecessors(node_id))

            active_subgraphs.append(self.G.subgraph(tissue_neighborhood).copy())

        # ── 5. Intersection Logic ────────────────────────────────────────
        # Catch-all fallback: If NO filters are actively chosen,
        # return a full copy of the current state immediately!
        if not active_subgraphs:
            result_engine = GraphEngine()
            result_engine.G = self.G.copy()
            return result_engine

        # If we have active filters, perform the strict intersection
        # Seed the intersection arrays with the very first active filter layer
        final_nodes = set(active_subgraphs[0].nodes)
        final_edges = set(active_subgraphs[0].edges(keys=True))

        # Core logic adjustment: Intersect only across the filters that were actually executed
        for sub in active_subgraphs[1:]:
            final_nodes &= set(sub.nodes)
            final_edges &= set(sub.edges(keys=True))

        # Build and populate the final isolated GraphEngine container
        filtered_engine = GraphEngine()

        # Add intersecting nodes with their original dictionary properties
        for node_id in final_nodes:
            filtered_engine.G.add_node(node_id, **self.G.nodes[node_id])

        # Add intersecting multi-edges with their original properties
        for u, v, k in final_edges:
            if u in final_nodes and v in final_nodes:
                filtered_engine.G.add_edge(u, v, key=k, **self.G.edges[u, v, k])

        return filtered_engine

    # ── Bulk operations ──────────────────────────────────────────────────
    def load_base_data(self, file_path: Path) -> None:
        raw_data = load_from_json(file_path)
        for node in raw_data.get("nodes"):
            self.add_node(Node(**node))
        for edge in raw_data.get("edges"):
            self.add_edge(Edge(**edge))

    def load_interaction_params(self, interaction_file_path: Path) -> None:
        raw_data = load_from_json(interaction_file_path)


        genotype_modifiers_metadata = raw_data["genotype_modifiers"]
        genotype_modifiers_metadata.pop("_description")

        # fix later, dumb to run through twice
        genotype_class = GenotypeModifiersContainer(
            **genotype_modifiers_metadata
        ).model_dump()

        for enzyme, gene_mod in genotype_class.items():
            self.modify_node(enzyme, "genotype_modifiers", gene_mod)

        competitive_inhibition_metadata = raw_data["competitive_inhibition"]
        competitive_inhibition_metadata.pop("_description")
        for enzyme_name, inter_metadata in competitive_inhibition_metadata.items():
            inter_metadata.pop("_description")
            for substrate_name, metadata_dict in inter_metadata["substrates"].items():
                self.modify_edge(
                    enzyme_name,
                    substrate_name,
                    "competitive_inhibition",
                    CompetitiveInhibitionSubstrate(**metadata_dict),
                    default_type=["Enzyme", "Substrate"],
                    default_visible=False,
                )

        phase2_conjugation_metadata = raw_data["phase2_conjugation_metadata"]
        phase2_conjugation_metadata.pop("_description")
        for enzyme_name, inter_metadata in phase2_conjugation_metadata.items():
            inter_metadata.pop("_description")
            for substrate_name, metadata_dict in inter_metadata["substrates"].items():
                self.modify_edge(
                    enzyme_name,
                    substrate_name,
                    "phase2_conjugation",
                    Phase2ConjugationSubstrate(**metadata_dict),
                    default_type=["Enzyme", "Substrate"],
                    default_visible=False,
                )

        enzyme_induction_metadata = raw_data["enzyme_induction_metadata"]
        enzyme_induction_metadata.pop("_description")
        for lifestyle_name, inter_metadata in enzyme_induction_metadata.items():
            inter_metadata.pop("_primary_mechanism")
            for enzyme, metadata_dict in inter_metadata.items():
                self.modify_edge(
                    lifestyle_name,
                    substrate_name,
                    "enzyme_induction",
                    EnzymeInduction(**metadata_dict),
                    default_type=["Lifestyle", "Enzyme"],
                    default_visible=False,
                )

    def load_exposure_params(self, exposure_file_path: Path ) -> None:

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

    def get_node(self, node_id: str) -> dict[str, Any] | None:
        if node_id in self.G:
            return dict(self.G.nodes[node_id])
        return None

    def neighbors(self, node_id: str) -> list[str]:
        if node_id not in self.G:
            return []
        return list(self.G.successors(node_id)) + list(self.G.predecessors(node_id))

    def nodes_by_type(self, node_type: str) -> list[dict[str, Any]]:
        return [
            data for _, data in self.G.nodes(data=True) if data.get("type") == node_type
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
