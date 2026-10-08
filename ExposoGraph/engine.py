"""NetworkX-backed graph engine for building and querying the knowledge graph."""

from __future__ import annotations

import json
import logging
import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import networkx as nx

from .config import GraphMode
from .grounding import prepare_knowledge_graph
from .interaction_schema import GSHConsumer, InductionRule
from .models import Edge, KnowledgeGraph, Node
from .parameter_provider import (
    INTERACTION_BLOCK_MARKER,
    INTERACTION_ENZYME_MARKER,
    INTERACTION_SUBSTRATE_MARKER,
    InteractionParameterProvider,
    JSONInteractionParameterProvider,
    KGInteractionParameterProvider,
)

_PACKAGE_DIR = Path(__file__).resolve().parent
_DEFAULT_GRAPH_DATA_PATH = _PACKAGE_DIR / "map" / "graph-data.json"
_DEFAULT_TISSUE_EXPRESSION_PATH = _PACKAGE_DIR / "data" / "tissue_expression_data_raw.json"
_DEFAULT_INTERACTION_PARAMETERS_PATH = _PACKAGE_DIR / "data" / "interaction_parameters.json"
_DEFAULT_FLUX_KINETIC_PARAMETERS_PATH = _PACKAGE_DIR / "data" / "kinetic_parameters.json"
_DEFAULT_PROXY_FLUX_PARAMETERS_PATH = _PACKAGE_DIR / "data" / "proxy_flux_parameters.json"
_DEFAULT_PROXY_FLUX_PROVENANCE_PATH = _PACKAGE_DIR / "data" / "proxy_flux_provenance.json"

# Residual activity retained by GSTM1/GSTT1 homozygous deletion carriers
# (see kinetic_parameters.json genotype_modifiers.special_cases).
_GST_NULL_RESIDUAL_ACTIVITY = 0.05

_DEFAULT_EXPOSURE_DB_PATH = _PACKAGE_DIR / "data" / "exposure_database.json"

_FALLBACK_QIVIVE_TISSUES: dict[str, dict[str, float]] = {
    "liver": {"mppgl_mg_per_g": 40.0, "organ_weight_g": 1500.0},
    "lung": {"mppgl_mg_per_g": 20.0, "organ_weight_g": 1000.0},
    "kidney": {"mppgl_mg_per_g": 12.0, "organ_weight_g": 300.0},
    "intestine": {"mppgl_mg_per_g": 35.0, "organ_weight_g": 900.0},
}

_STEADY_STATE_DEFAULTS: dict[str, float] = {
    "body_weight_kg": 70.0,
    "volume_l_per_kg": 0.7,
    "absorption_fraction": 1.0,
    "exposure_frequency_per_day": 1.0,
    "cardiac_output_l_per_day": 7200.0,
    "background_clearance_rate_per_day": 0.05,
    "reactive_intermediate_loss_rate_per_day": 1.0,
    "detoxified_metabolite_loss_rate_per_day": 1.0,
    "flux_rate_scale_per_day": 1.0,
    "rate_reference_concentration_uM": 1.0,
}

_FALLBACK_STEADY_STATE_TISSUES: dict[str, dict[str, float]] = {
    "liver": {
        "organ_weight_g": 1500.0,
        "tissue_partition_coefficient": 1.0,
        "tissue_blood_flow_fraction": 0.25,
    },
    "lung": {
        "organ_weight_g": 1000.0,
        "tissue_partition_coefficient": 0.8,
        "tissue_blood_flow_fraction": 1.0,
    },
    "kidney": {
        "organ_weight_g": 300.0,
        "tissue_partition_coefficient": 1.1,
        "tissue_blood_flow_fraction": 0.2,
    },
    "intestine": {
        "organ_weight_g": 900.0,
        "tissue_partition_coefficient": 0.9,
        "tissue_blood_flow_fraction": 0.12,
    },
    "bladder": {
        "organ_weight_g": 150.0,
        "tissue_partition_coefficient": 0.7,
        "tissue_blood_flow_fraction": 0.02,
    },
    "breast": {
        "organ_weight_g": 500.0,
        "tissue_partition_coefficient": 1.4,
        "tissue_blood_flow_fraction": 0.03,
    },
    "colon": {
        "organ_weight_g": 600.0,
        "tissue_partition_coefficient": 0.9,
        "tissue_blood_flow_fraction": 0.08,
    },
    "prostate": {
        "organ_weight_g": 30.0,
        "tissue_partition_coefficient": 0.8,
        "tissue_blood_flow_fraction": 0.01,
    },
    "esophagus": {
        "organ_weight_g": 40.0,
        "tissue_partition_coefficient": 0.8,
        "tissue_blood_flow_fraction": 0.01,
    },
    "skin": {
        "organ_weight_g": 3300.0,
        "tissue_partition_coefficient": 1.2,
        "tissue_blood_flow_fraction": 0.05,
    },
}

def _positive_context_float(context: Mapping[str, Any], key: str, fallback: float) -> float:
    """Read a positive numeric context value with a conservative fallback."""
    try:
        value = float(context.get(key, fallback))
    except (TypeError, ValueError):
        return fallback
    if value <= 0 or not math.isfinite(value):
        return fallback
    return value


def _bounded_fraction_context(context: Mapping[str, Any], key: str, fallback: float) -> float:
    """Read a fraction constrained to the open interval used by PBPK rates."""
    try:
        value = float(context.get(key, fallback))
    except (TypeError, ValueError):
        return fallback
    if not math.isfinite(value):
        return fallback
    return min(max(value, 1e-6), 1.0)

# Fields a flux-parameter term carries that are represented on the
# FluxReaction record itself (or as provenance pointers) rather than inside
# its verbatim ``params`` payload. ``gene`` and ``note`` stay in ``params``
# because the flux proxy-term helpers consume them directly.
_FLUX_RESERVED_TERM_FIELDS = frozenset(
    {"graph_node_id", "rate_law", "equation", "confidence", "notes", "provenance_ref", "sources"}
)

# Audit metadata ``_apply_flux_edge_kinetics`` adds to each baked edge
# payload on top of the term's parameter dict. ``get_edge_flux_reactions``
# strips exactly these keys (and no parameter keys) when reconstructing a
# FluxReaction from an edge payload, so the two records stay identical.
_FLUX_PAYLOAD_META_FIELDS = frozenset(
    {
        "rate_law",
        "role",
        "confidence",
        "provenance_ref",
        "source",
        "dispatch_status",
        "graph_group",
        "sources",
        "notes",
        "enzyme_id",
    }
)

# Pathway-block names that map unambiguously onto the coarse flux role
# vocabulary. Idiosyncratic blocks (e.g. ``ethanol_oxidation``) are recorded
# as "other" with their verbatim pathway name until the flux JSONs carry
# explicit per-term role annotations.
_FLUX_ROLE_KINDS = {
    "activation": "activation",
    "activation_terms": "activation",
    "detoxification": "detoxification",
    "detox_terms": "detoxification",
    "repair_terms": "repair",
}

# Bridge mapping from flux-parameter class keys to the Carcinogen node
# ``group`` labels used in graph-data.json. Transitional: once the flux JSONs
# carry explicit ``graph_group`` / ``carcinogen_node_ids`` annotations, the
# JSON-declared values take precedence and this mapping is only a fallback.
_FLUX_CLASS_GRAPH_GROUPS = {
    "PAH": "PAHs",
    "Aflatoxin": "Mycotoxins",
    "Aldehyde": "Aldehydes",
    "Nitrosamine": "Nitrosamines",
    "NDMA": "Nitrosamines",
    "NDEA": "Nitrosamines",
    "HCA": "HCAs",
    "AromaticAmines": "Aromatic Amines",
    "EstrogenMetabolites": "Estrogen",
    "Benzene": "Benzene",
    "VinylChloride": "Vinyl Chloride",
    "ChlorinatedSolvent": "Chlorinated Solvents",
    "UV_Radiation": "UV Radiation",
    "Dioxin": "Dioxins",
    "HeavyMetal": "Heavy Metals",
}

# Edge types that can carry flux-reaction scope between a Carcinogen (or its
# metabolites) and an Enzyme; used by get_flux_reaction_coverage.
_FLUX_SCOPE_EDGE_TYPES = frozenset({"SUBSTRATE_OF", "DETOXIFIED_BY", "REPAIRED_BY"})


@dataclass(frozen=True)
class FluxReaction:
    """One reaction term of the flux contract, per reaction (not per class).

    Built by :meth:`GraphEngine._apply_flux_parameters` from
    ``kinetic_parameters.json`` and ``proxy_flux_parameters.json``. The
    engine-side index is the interim bridge: group-level flux is an
    aggregation over these records, and once the graph carries role-typed
    scope edges and per-term kinetics, scope can be derived from edges
    instead of the parameter rosters.
    """

    carcinogen_class: str
    term_key: str
    pathway: str
    role: str
    rate_law: str
    params: dict[str, Any]
    confidence: str
    provenance_ref: str
    enzyme_id: str | None
    graph_group: str | None
    source: str
    # Audit metadata retained so the edge-baking step
    # (``_apply_flux_edge_kinetics``) can attach the full term payload --
    # including PMID lists and notes -- to ``Edge.kinetics["flux_terms"]``.
    sources: list[str] | None = None
    notes: str | None = None

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
        self._interaction_parameters: dict[str, Any] | None = None
        self._parameter_provider: JSONInteractionParameterProvider | None = None
        self._kg_parameter_provider: InteractionParameterProvider | None = None
        self._flux_reactions_by_class: dict[str, list[FluxReaction]] | None = None
        self._flux_class_sources: dict[str, set[str]] | None = None
        self._flux_shadowed_reactions: list[FluxReaction] = []
        self._flux_aggregation_by_class: dict[str, dict[str, Any]] | None = None
        self._flux_metadata: dict[str, Any] | None = None
        self._flux_genotype_modifiers: dict[str, Any] | None = None
        self._flux_proxy_exposure_defaults: dict[str, dict[str, Any]] = {}
        self._flux_proxy_class_cfg: dict[str, dict[str, Any]] = {}
        self._flux_group_class_md_cache: dict[str, dict[str, Any]] | None = None
        self._flux_group_class_md_graph: Any = None
        self._flux_class_signal_cache: dict[str, dict[str, Any]] | None = None
        self._flux_class_signal_graph: Any = None
        self._flux_proxy_provenance: dict[str, Any] | None = None
        self._exposure_database: dict[str, Any] | None = None
        self._tissue_expression_raw: dict[str, dict[str, float]] | None = None
        self._tissue_expression_normalized: dict[str, dict[str, float]] | None = None

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
        self._interaction_parameters = None
        self._parameter_provider = None
        self._kg_parameter_provider = None
        self._flux_reactions_by_class = None
        self._flux_class_sources = None
        self._flux_shadowed_reactions = []
        self._flux_aggregation_by_class = None
        self._flux_metadata = None
        self._flux_genotype_modifiers = None
        self._flux_proxy_exposure_defaults = {}
        self._flux_proxy_class_cfg = {}
        self._flux_group_class_md_cache = None
        self._flux_group_class_md_graph = None
        self._flux_class_signal_cache = None
        self._flux_class_signal_graph = None
        self._flux_proxy_provenance = None
        self._exposure_database = None
        self._tissue_expression_raw = None
        self._tissue_expression_normalized = None

    def load_reference_graph(
        self,
        *,
        graph_data_path: str | Path | None = None,
        tissue_expression_path: str | Path | None = None,
        interaction_parameters_path: str | Path | None = None,
        kinetic_parameters_path: str | Path | None = None,
        proxy_flux_parameters_path: str | Path | None = None,
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
        3. Flux parameters from ``data/kinetic_parameters.json`` and
           ``data/proxy_flux_parameters.json`` build the engine-side
           flux-reaction index -- see :meth:`_apply_flux_parameters`. This
           does not mutate the graph; it makes the quantitative flux
           contract queryable through the engine
           (:meth:`get_flux_reactions` and friends).
        4. Substrate-bound flux terms are then baked onto their
           substrate→enzyme edges under ``Edge.kinetics["flux_terms"]`` --
           see :meth:`_apply_flux_edge_kinetics`. Like the interaction
           overlay, the parameter JSONs (not graph-data.json) remain the
           trusted source; the edge payload is derived at load time.

        Returns the combined warning messages from all steps.
        """
        from .exporter import parse_graph_artifact  # local import avoids an import cycle

        resolved_graph_path = Path(graph_data_path) if graph_data_path else _DEFAULT_GRAPH_DATA_PATH
        kg = parse_graph_artifact(resolved_graph_path)
        warnings = self.load(kg)
        warnings.extend(self._apply_tissue_expression(tissue_expression_path))
        warnings.extend(self._apply_interaction_parameters(interaction_parameters_path))
        warnings.extend(self._apply_flux_parameters(kinetic_parameters_path, proxy_flux_parameters_path))
        warnings.extend(self._apply_flux_edge_kinetics())
        return warnings

    def _ensure_tissue_expression(self) -> dict[str, dict[str, float]]:
        """Lazily load and normalize ``tissue_expression_data_raw.json``.

        Retains both the raw per-tissue nTPM values and the per-enzyme
        divide-by-max normalization (most-expressing tissue = 1.0) that
        :meth:`_apply_tissue_expression` bakes onto enzyme nodes. The
        normalized table is the single tissue-weight source served to flux
        models via :meth:`get_tissue_expression`.
        """
        if self._tissue_expression_raw is None or self._tissue_expression_normalized is None:
            try:
                doc = json.loads(_DEFAULT_TISSUE_EXPRESSION_PATH.read_text(encoding="utf-8"))
                expression: dict[str, dict[str, float]] = doc["expression"]
            except (OSError, json.JSONDecodeError, KeyError) as exc:
                # This file is the authoritative flux tissue-weight table;
                # silently degrading every lookup to a fallback weight would
                # compute a full profile on invented numbers.
                raise RuntimeError(
                    f"Tissue expression source {_DEFAULT_TISSUE_EXPRESSION_PATH} is the "
                    "authoritative tissue-weight table and must be readable with an "
                    f"'expression' block: {exc}"
                ) from exc
            normalized: dict[str, dict[str, float]] = {}
            for gene, raw in expression.items():
                max_raw = max(raw.values()) if raw else 0.0
                normalized[gene] = (
                    {tissue: value / max_raw for tissue, value in raw.items()}
                    if max_raw
                    else dict.fromkeys(raw, 0.0)
                )
            self._tissue_expression_raw = expression
            self._tissue_expression_normalized = normalized
        return self._tissue_expression_normalized

    def get_tissue_expression(self, gene: str) -> dict[str, float] | None:
        """Return normalized per-tissue expression weights for a gene.

        Weights come from ``tissue_expression_data_raw.json`` (GTEx v8
        nTPM via Human Protein Atlas), normalized per enzyme by its
        most-expressing tissue. Returns ``None`` when the gene has no
        entry in the expression source.
        """
        return self._ensure_tissue_expression().get(gene)

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
        if resolved_path == _DEFAULT_TISSUE_EXPRESSION_PATH:
            # Share the retained (and flux-served) normalized table so node
            # baking and get_tissue_expression can never diverge.
            self._ensure_tissue_expression()
            expression = self._tissue_expression_raw or {}
            normalized_by_gene = self._tissue_expression_normalized or {}
        else:
            # This file's "expression" table is raw nTPM values with no
            # normalization applied -- the divide-by-max step below is this
            # method's own responsibility, unchanged regardless of source file.
            expression = json.loads(resolved_path.read_text(encoding="utf-8"))["expression"]
            normalized_by_gene = {}
            for gene, raw in expression.items():
                max_raw = max(raw.values()) if raw else 0.0
                normalized_by_gene[gene] = (
                    {tissue: value / max_raw for tissue, value in raw.items()}
                    if max_raw
                    else dict.fromkeys(raw, 0.0)
                )

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

            node_data["tissue_weights_raw"] = raw
            node_data["tissue_weights"] = normalized_by_gene[enzyme_id]

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

        Each enzyme block declares ``target_node_id``. Each substrate entry
        declares ``source_node_id``. Those are the ends of the reaction edge:
        ``source_node_id -[SUBSTRATE_OF|DETOXIFIED_BY]-> target_node_id``.
        ``PRODUCES`` edges are not carriers for these kinetics.

        The entry's remaining fields are set, unchanged, as ``kinetics`` on
        that edge, plus three self-description markers
        (``interaction_enzyme``, ``interaction_substrate``,
        ``interaction_block``). This overwrites any ``kinetics`` already on
        the edge. Flux bindings already stored under ``flux_terms`` are kept.

        A second entry for the same source and target is reported and
        skipped; the earlier entry is kept. A missing id, an id that is not
        a node, or no reaction edge is a warning, not an error.

        Returns a list of warning messages.
        """
        resolved_path = Path(path) if path else _DEFAULT_INTERACTION_PARAMETERS_PATH
        source_data = json.loads(resolved_path.read_text(encoding="utf-8"))

        # Retain the parsed document and a typed provider over it so the
        # engine can serve the non-edge parameter blocks (enzyme_induction,
        # gsh_depletion, genotype_modifiers, interaction_rules) through
        # getters instead of every consumer reading data/*.json directly.
        # This is the interim bridge until those parameters are carried as
        # graph node/edge attributes (see the interaction-parameter getters
        # below the mutation section).
        self._interaction_parameters = source_data
        if resolved_path.name == "interaction_parameters.json":
            self._parameter_provider = JSONInteractionParameterProvider(resolved_path.parent)
        else:
            # Custom-named overlay documents have no provider-side
            # counterpart file; fall back to the bundled default provider.
            self._parameter_provider = None

        warnings: list[str] = []
        pending: dict[tuple[str, str], dict[str, Any]] = {}
        pending_block: dict[tuple[str, str], str] = {}
        reaction_types = {"SUBSTRATE_OF", "DETOXIFIED_BY"}
        pending_substrate: dict[tuple[str, str], str] = {}
        for block_name in self._INTERACTION_PARAMETER_BLOCKS:
            block = source_data.get(block_name, {})
            for enzyme_key, enzyme_block in block.items():
                if enzyme_key == "_description" or not isinstance(enzyme_block, dict):
                    continue
                if "substrates" not in enzyme_block:
                    continue
                target_id = enzyme_block.get("target_node_id")
                if not target_id:
                    warnings.append(
                        f"No target_node_id declared for {block_name} enzyme: {enzyme_key}"
                    )
                    continue
                if target_id not in self.G:
                    warnings.append(
                        f"target_node_id {target_id!r} for {block_name} enzyme "
                        f"{enzyme_key} is not a node in the graph"
                    )
                    continue
                for substrate_key, params in enzyme_block.get("substrates", {}).items():
                    if not isinstance(params, dict):
                        continue
                    source_id = params.get("source_node_id")
                    if not source_id:
                        warnings.append(
                            f"No source_node_id declared for {block_name} "
                            f"substrate: {enzyme_key}/{substrate_key}"
                        )
                        continue
                    if source_id not in self.G:
                        warnings.append(
                            f"source_node_id {source_id!r} for {block_name} substrate "
                            f"{enzyme_key}/{substrate_key} is not a node in the graph"
                        )
                        continue
                    key = (target_id, source_id)
                    if key in pending:
                        warnings.append(
                            f"Duplicate source_node_id {source_id!r} for "
                            f"{block_name} target {target_id}: {substrate_key} conflicts with "
                            f"{pending_substrate[key]}; kept the earlier entry"
                        )
                        continue
                    kinetics = {k: v for k, v in params.items() if k != "source_node_id"}
                    # Self-description markers so the graph-walk provider
                    # (KGInteractionParameterProvider) can reconstruct
                    # typed records from the edge kinetics without
                    # re-reading the JSON document: the substrate key and
                    # source block are not recoverable from node identity
                    # alone.
                    kinetics[INTERACTION_ENZYME_MARKER] = enzyme_key
                    kinetics[INTERACTION_SUBSTRATE_MARKER] = substrate_key
                    kinetics[INTERACTION_BLOCK_MARKER] = block_name
                    pending[key] = kinetics
                    pending_block[key] = block_name
                    pending_substrate[key] = substrate_key

        applied: set[tuple[str, str]] = set()
        for source_id, target_id, edge_data in self.G.edges(data=True):
            if edge_data.get("type") not in reaction_types:
                continue
            key = (target_id, source_id)
            if key not in pending:
                continue
            previous = edge_data.get("kinetics")
            edge_data["kinetics"] = dict(pending[key])
            # The engine-owned flux-term namespace (see
            # ``_apply_flux_edge_kinetics``) survives this overwrite:
            # re-applying interaction parameters must not clobber
            # flux bindings baked in an earlier pass.
            if isinstance(previous, dict) and "flux_terms" in previous:
                edge_data["kinetics"]["flux_terms"] = previous["flux_terms"]
            applied.add(key)

        for target_id, source_id in pending:
            if (target_id, source_id) not in applied:
                block_name = pending_block[(target_id, source_id)]
                warnings.append(
                    f"No source-to-target reaction edge for {block_name} pair: "
                    f"{source_id} -> {target_id}"
                )

        return warnings


    def _flux_term_enzyme_id(self, term: Mapping[str, Any], term_key: str) -> str | None:
        """Resolve a flux term's enzyme node.

        Exact term keys win. Alias keys (``CYP2E1_liver``) use
        ``target_node_id`` when that target is an Enzyme, then ``gene``.
        """
        if term_key in self.G:
            return term_key
        target = term.get("target_node_id")
        if (
            isinstance(target, str)
            and target in self.G
            and self.G.nodes[target].get("type") == "Enzyme"
        ):
            return target
        gene = term.get("gene")
        if isinstance(gene, str) and gene in self.G:
            return gene
        return None

    def _apply_flux_parameters(
        self,
        kinetic_path: str | Path | None = None,
        proxy_path: str | Path | None = None,
    ) -> list[str]:
        """(Re)build the flux-reaction side index from the two flux parameter JSONs.

        Reads ``data/kinetic_parameters.json`` (mechanistic Km/Vmax classes)
        and ``data/proxy_flux_parameters.json`` (semi-quantitative proxy
        classes) and builds ``self._flux_reactions_by_class`` -- an
        engine-side index, NOT node attributes: non-enzyme proxy terms
        (``general_ROS``, ``AHRR_feedback``) have no node to live on, and node
        storage would leak into ``to_dict()``/exports. GraphEngine is still
        the contract because the public getters live on it.

        Precedence: classes present in both files (ChlorinatedSolvent,
        Dioxin, HeavyMetal) keep only their proxy entries, matching current
        flux-engine dispatch behavior; the kinetic entries are shadowed
        (visible via ``get_flux_reaction_coverage`` source reporting) and
        retained in ``_flux_shadowed_reactions`` so
        :meth:`_apply_flux_edge_kinetics` can still bake their graph-edge
        bindings.

        Enzyme linkage: a term whose key exactly matches a node id resolves
        to it. Alias keys (``CYP2E1_liver``, ``ADH1B_star1``) resolve to
        ``target_node_id`` when that target is an Enzyme node, otherwise to
        a ``gene`` that names a node, otherwise ``None``.

        Returns warning messages, mirroring the other overlay methods.
        """
        kinetic_resolved = Path(kinetic_path) if kinetic_path else _DEFAULT_FLUX_KINETIC_PARAMETERS_PATH
        proxy_resolved = Path(proxy_path) if proxy_path else _DEFAULT_PROXY_FLUX_PARAMETERS_PATH

        warnings: list[str] = []
        index: dict[str, list[FluxReaction]] = {}
        sources: dict[str, set[str]] = {}
        aggregations: dict[str, dict[str, Any]] = {}

        def _record(
            cls: str,
            term_key: str,
            block: str,
            term: Mapping[str, Any],
            class_data: Mapping[str, Any],
            source: str,
            provenance_prefix: str,
        ) -> FluxReaction:
            enzyme_id = self._flux_term_enzyme_id(term, term_key)
            declared_group = class_data.get("graph_group")
            return FluxReaction(
                carcinogen_class=cls,
                term_key=term_key,
                pathway=block,
                role=_FLUX_ROLE_KINDS.get(block, "other"),
                rate_law=str(
                    term.get("rate_law")
                    or term.get("equation")
                    or class_data.get("kinetics_model")
                    or ""
                ),
                params={k: v for k, v in term.items() if k not in _FLUX_RESERVED_TERM_FIELDS},
                confidence=str(term.get("confidence") or class_data.get("confidence_overall") or ""),
                provenance_ref=str(
                    term.get("provenance_ref") or f"{provenance_prefix}.{block}.{term_key}"
                ),
                enzyme_id=enzyme_id,
                graph_group=str(declared_group) if declared_group else _FLUX_CLASS_GRAPH_GROUPS.get(cls),
                source=source,
                sources=list(term["sources"]) if isinstance(term.get("sources"), list) else None,
                notes=str(term["notes"]) if term.get("notes") is not None else None,
            )

        # kinetic_parameters.json -- mechanistic classes
        try:
            kinetic_doc = json.loads(kinetic_resolved.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            warnings.append(f"Could not read flux kinetic parameters at {kinetic_resolved}: {exc}")
            kinetic_doc = {}
        for cls, cls_data in kinetic_doc.get("carcinogen_classes", {}).items():
            if not isinstance(cls_data, dict):
                continue
            reactions = []
            for block, terms in cls_data.get("pathways", {}).items():
                if not isinstance(terms, dict):
                    continue
                for term_key, term in terms.items():
                    if term_key.startswith("_") or not isinstance(term, dict):
                        continue
                    reactions.append(
                        _record(cls, term_key, block, term, cls_data, "kinetic_parameters", f"carcinogen_classes.{cls}.pathways")
                    )
            index[cls] = reactions
            sources[cls] = {"kinetic_parameters"}
            aggregations[cls] = dict(cls_data.get("aggregation", {})) if isinstance(cls_data.get("aggregation"), dict) else {}

        # proxy_flux_parameters.json -- semi-quantitative classes; proxy wins
        try:
            proxy_doc = json.loads(proxy_resolved.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            warnings.append(f"Could not read proxy flux parameters at {proxy_resolved}: {exc}")
            proxy_doc = {}
        for cls, cls_data in proxy_doc.get("classes", {}).items():
            if not isinstance(cls_data, dict):
                continue
            reactions = []
            for block, terms in cls_data.items():
                if not block.endswith("_terms") or not isinstance(terms, dict):
                    continue
                for term_key, term in terms.items():
                    if term_key.startswith("_") or not isinstance(term, dict):
                        continue
                    reactions.append(
                        _record(cls, term_key, block, term, cls_data, "proxy_flux_parameters", f"classes.{cls}")
                    )
            if cls in index:
                # Proxy wins for dispatch, but the kinetic entries are
                # retained as ``shadowed`` reactions so the graph-edge
                # baking step can still attach their bindings (dispatch
                # choice and graph knowledge are separate concerns).
                self._flux_shadowed_reactions.extend(index[cls])
                sources[cls].add("proxy_flux_parameters")
            else:
                sources[cls] = {"proxy_flux_parameters"}
            index[cls] = reactions
            if isinstance(cls_data.get("exposure_default"), dict):
                self._flux_proxy_exposure_defaults[cls] = dict(cls_data["exposure_default"])
            self._flux_proxy_class_cfg[cls] = dict(cls_data)

        self._flux_reactions_by_class = index
        self._flux_class_sources = sources
        self._flux_aggregation_by_class = aggregations
        self._flux_metadata = dict(kinetic_doc.get("metadata", {})) if isinstance(kinetic_doc.get("metadata"), dict) else {}
        self._flux_genotype_modifiers = (
            dict(kinetic_doc.get("genotype_modifiers", {}))
            if isinstance(kinetic_doc.get("genotype_modifiers"), dict)
            else {}
        )
        return warnings

    def _apply_flux_edge_kinetics(self) -> list[str]:
        """Bake edge-bound flux terms onto their source→target edges.

        Iterates the flux-reaction index built by
        :meth:`_apply_flux_parameters` and, for every term carrying both
        a ``source_node_id`` and a ``target_node_id`` binding, attaches the
        term's verbatim parameter payload to the graph edge between those
        two nodes, under ``Edge.kinetics["flux_terms"]``, nested as
        ``flux_terms[carcinogen_class][pathway][term_key]``. The binding is
        endpoint-generic by design: ``source_node_id`` / ``target_node_id``
        cover substrate-to-enzyme reaction edges and non-enzymatic terms -- formation /
        driver terms anchored on ``carcinogen --FORMS_ADDUCT--> lesion``
        edges, receptor-feedback terms on ``INHIBITS`` edges, and so on.
        The nesting is
        required because several edges carry multiple terms -- diplotype
        variants (``ALDH2_star1`` / ``ALDH2_star1_star2`` with different
        kinetics) and competing reaction channels catalyzed by one enzyme
        (``CYP3A4`` activation vs. ``CYP3A4_AFQ1`` detoxification), and one
        edge (``Formaldehyde → ADH5``) carries terms from two different
        carcinogen classes. ``kinetic_parameters.json`` stays the single
        source of truth: the payload is baked at load time, not copied into
        ``graph-data.json``.

        Terms with bindings whose edge does not exist (e.g. the known
        ``PhIP → CYP1B1`` / ``PhIP → NAT2`` scope gaps) are skipped and
        reported as warnings rather than raising. Binding-less terms
        (substrate-only, driver/proxy, non-enzymatic) are ignored here by
        design. Dual-source classes (ChlorinatedSolvent, Dioxin, HeavyMetal)
        dispatch on their proxy entries, but their shadowed kinetic entries
        are baked here too -- dispatch choice and graph knowledge are
        separate concerns -- with ``source: kinetic_parameters`` recording
        where each payload came from and ``dispatch_status``
        (``active`` | ``shadowed``) declaring whether the term is a member
        of the active dispatch roster or carried graph knowledge. The
        graph walk (:meth:`get_edge_flux_reactions`) skips shadowed
        payloads, making the active-vs-carried status explicit on each
        payload for the graph-native roster cutover (the GSTT1/TCE
        conjugation term is the canonical case: its kinetic payload is
        carried, while the published model treats GSTT1-mediated TCE
        conjugation as bioactivation).

        Where an edge also carries interaction kinetics from
        :meth:`_apply_interaction_parameters`, both coexist: the interaction
        block keeps its flat keys (``Km_uM``/``Ki_uM``/...) and the flux
        terms live under the ``flux_terms`` namespace.

        Each baked payload is self-contained: the term's parameter dict
        (Km/Vmax/CLint, genotype labels, ...) plus the audit metadata the
        reaction index carries (``rate_law``, ``role``, ``confidence``,
        ``provenance_ref``, ``source``, ``graph_group``, ``sources`` PMID
        list, ``notes``), so a consumer walking the graph can compute and
        audit the reaction from the edge alone.

        Returns a list of warning messages, mirroring the other overlay
        methods.
        """
        warnings: list[str] = []
        if self._flux_reactions_by_class is None:
            self._apply_flux_parameters()
        baked: list[tuple[FluxReaction, bool]] = [
            (reaction, False)
            for reactions in self._flux_reactions_by_class.values()
            for reaction in reactions
        ]
        # Shadowed reactions (dual-source classes): carried for graph
        # knowledge, excluded from the active dispatch roster.
        baked.extend((reaction, True) for reaction in self._flux_shadowed_reactions)
        for reaction, shadowed in baked:
            source = reaction.params.get("source_node_id")
            target = reaction.params.get("target_node_id")
            if not source or not target:
                continue
            if not self.G.has_edge(source, target):
                warnings.append(
                    f"No edge for flux binding: {reaction.carcinogen_class}/{reaction.term_key} "
                    f"({source} -> {target})"
                )
                continue
            edge_view = self.G[source][target]
            edge_data = next(iter(edge_view.values()))
            kinetics = edge_data.setdefault("kinetics", {})
            flux_terms = kinetics.setdefault("flux_terms", {})
            by_class = flux_terms.setdefault(reaction.carcinogen_class, {})
            by_pathway = by_class.setdefault(reaction.pathway, {})
            payload = dict(reaction.params)
            payload["rate_law"] = reaction.rate_law
            payload["role"] = reaction.role
            payload["confidence"] = reaction.confidence
            payload["provenance_ref"] = reaction.provenance_ref
            payload["source"] = reaction.source
            payload["dispatch_status"] = "shadowed" if shadowed else "active"
            if reaction.graph_group is not None:
                payload["graph_group"] = reaction.graph_group
            # The roster's *resolved* enzyme id (exact key / Enzyme-typed
            # ``target_node_id`` / gene), which can differ from the binding
            # annotation (alias pseudo-terms like ``CYP2E1_liver``).
            # Baked so the edge payload reconstructs the index record exactly.
            if reaction.enzyme_id is not None:
                payload["enzyme_id"] = reaction.enzyme_id
            if reaction.sources is not None:
                payload["sources"] = reaction.sources
            if reaction.notes is not None:
                payload["notes"] = reaction.notes
            by_pathway[reaction.term_key] = payload

        # Class-level signal blocks with declared bindings (Dioxin's AhR
        # occupancy Hill term, bound to ``TCDD --AGONIZES--> AHR``): bake
        # the block onto its edge under a dedicated namespace so the
        # class's signaling parameters are served graph-first, mirroring
        # how term payloads ride their binding edges. Signal blocks are
        # class config, not roster terms -- a separate namespace keeps
        # the roster walk and coverage accounting untouched.
        for cls_name, class_cfg in self._flux_proxy_class_cfg.items():
            signal = class_cfg.get("signal")
            if not isinstance(signal, dict):
                continue
            signal_source = signal.get("source_node_id")
            signal_target = signal.get("target_node_id")
            if not signal_source or not signal_target:
                continue
            if not self.G.has_edge(signal_source, signal_target):
                warnings.append(
                    f"No edge for flux signal binding: {cls_name} "
                    f"({signal_source} -> {signal_target})"
                )
                continue
            edge_view = self.G[signal_source][signal_target]
            edge_data = next(iter(edge_view.values()))
            kinetics = edge_data.setdefault("kinetics", {})
            kinetics.setdefault("flux_class_signal", {})[cls_name] = dict(signal)
        # The bake mutates the graph in place, so any previously built
        # signal cache (from an earlier get_flux_class_config call) is
        # stale; drop it unconditionally.
        self._flux_class_signal_cache = None
        self._flux_class_signal_graph = None
        return warnings

    # ── Interaction-parameter access ─────────────────────────────────────

    def get_interaction_parameters(self) -> dict[str, Any]:
        """Return the parsed ``interaction_parameters.json`` document.

        The engine retains the document it applied in
        :meth:`_apply_interaction_parameters` so downstream modules can source
        interaction data (enzyme induction, GSH depletion, genotype modifiers,
        interaction rules) through the engine instead of reading
        ``data/interaction_parameters.json`` directly. When no overlay has been
        applied, the bundled default document is read on first use.

        This is the documented interim bridge for the four parameter blocks
        that are not yet carried as graph node/edge attributes. The returned
        mapping is the engine's live copy -- treat it as read-only.
        """
        if self._interaction_parameters is None:
            self._interaction_parameters = json.loads(
                _DEFAULT_INTERACTION_PARAMETERS_PATH.read_text(encoding="utf-8")
            )
        return self._interaction_parameters

    def _has_baked_interaction_kinetics(self) -> bool:
        """Return True if edge kinetics carry an interaction-parameter bake.

        ``G.number_of_edges()`` alone is not a safe gate: a graph merged or
        loaded by other means can have edges without the interaction bake,
        and a KG provider over it would silently serve an empty (or
        partial) kinetic roster. The bake's block marker is the actual
        invariant ``get_parameter_provider`` needs.
        """
        for _source_id, _target_id, edge_data in self.G.edges(data=True):
            kinetics = edge_data.get("kinetics")
            if (
                isinstance(kinetics, dict)
                and kinetics.get(INTERACTION_BLOCK_MARKER) == "competitive_inhibition"
            ):
                return True
        return False

    def get_parameter_provider(self) -> InteractionParameterProvider:
        """Return the typed interaction-parameter provider owned by the engine.

        With the reference graph loaded, this is the hybrid
        ``KGInteractionParameterProvider``: kinetic records (competitive
        interactions, per-enzyme and per-carcinogen reactions) are
        reconstructed from the edge kinetics this engine baked, while
        induction rules, GSH consumers, and parameter evidence fall
        through to the JSON document the engine applied. A bare engine (no
        graph loaded) returns the JSON provider unchanged. Consumers
        needing typed records should take the provider from here rather
        than constructing their own, so the whole application reads one
        copy of the parameter data.
        """
        if self._parameter_provider is None:
            self._parameter_provider = JSONInteractionParameterProvider()
        if self._has_baked_interaction_kinetics():
            if self._kg_parameter_provider is None:
                self._kg_parameter_provider = KGInteractionParameterProvider(
                    self, data_dir=self._parameter_provider.data_dir
                )
            return self._kg_parameter_provider
        return self._parameter_provider

    def get_induction_rules(
        self,
        exposure_context: str | None = None,
        tissue: str | None = None,
    ) -> list[InductionRule]:
        """Return enzyme-induction rules as typed ``InductionRule`` records.

        ``exposure_context`` optionally filters to one lifestyle section
        (e.g. ``"smoking"``, ``"chronic_alcohol"``, ``"TCDD_dioxin"``);
        ``tissue`` mirrors the provider's tissue filter.
        """
        rules = self.get_parameter_provider().get_induction_rules(tissue=tissue)
        if exposure_context is not None:
            rules = [rule for rule in rules if rule.exposure_context == exposure_context]
        return rules

    def get_gsh_consumers(self, tissue: str | None = None) -> list[GSHConsumer]:
        """Return GSH-depletion consumer records, optionally tissue-filtered."""
        return self.get_parameter_provider().get_gsh_consumers(tissue=tissue)

    def get_gsh_parameters(self) -> dict[str, Any]:
        """Return the ``gsh_depletion`` block (scalars, consumers, biology model)."""
        return self.get_interaction_parameters().get("gsh_depletion", {})

    def get_genotype_modifiers(self, gene: str | None = None) -> dict[str, Any]:
        """Return genotype-modifier tables, or one gene's table when ``gene`` is given."""
        block = self.get_interaction_parameters().get("genotype_modifiers", {})
        if gene is None:
            return block
        return block.get(gene, {})

    def get_interaction_rules(self) -> dict[str, Any]:
        """Return the ``interaction_rules`` block (thresholds, synergies, antagonisms)."""
        return self.get_interaction_parameters().get("interaction_rules", {})

    # ── Flux-reaction access ─────────────────────────────────────────────

    def _ensure_flux_index(self) -> None:
        """Lazily build the flux-reaction index from the bundled default files."""
        if self._flux_reactions_by_class is None or self._flux_class_sources is None:
            self._apply_flux_parameters()

    def get_flux_reactions(
        self,
        carcinogen_class: Any,
        role: str | None = None,
    ) -> list[FluxReaction]:
        """Return the flux-reaction terms for a carcinogen class.

        ``carcinogen_class`` accepts either the plain class key (``"PAH"``)
        or a ``CarcinogenClass`` enum member. ``role`` optionally filters on
        the coarse role vocabulary (``"activation"`` / ``"detoxification"`` /
        ``"repair"``); ``"other"`` matches terms whose JSON pathway block does
        not map onto that vocabulary yet (e.g. Aldehyde's
        ``ethanol_oxidation``). For dual-source classes
        (ChlorinatedSolvent, Dioxin, HeavyMetal) the proxy entries are
        returned, matching current flux-engine dispatch.
        """
        self._ensure_flux_index()
        cls = getattr(carcinogen_class, "value", carcinogen_class)
        reactions = list(self._flux_reactions_by_class.get(cls, []))
        if role is not None:
            reactions = [reaction for reaction in reactions if reaction.role == role]
        return reactions

    def get_edge_flux_reactions(self, carcinogen_class: Any) -> list[FluxReaction]:
        """Return flux-reaction terms with their parameters read from the graph.

        Graph-walking counterpart to :meth:`get_flux_reactions`: walks every
        edge of the loaded graph, collects the
        ``Edge.kinetics["flux_terms"]`` payloads baked by
        :meth:`_apply_flux_edge_kinetics` for *carcinogen_class*, and returns
        the class's reaction roster with each edge-anchored term's record
        reconstructed from its edge payload. Terms carried on the class's
        CarcinogenGroup node
        (``flux_class_metadata[cls].class_level_terms``) are likewise
        reconstructed from the group-node payload.

        When the walk's coverage is complete -- every side-index term has a
        graph representation -- the returned roster is *graph-native*:
        membership comes from the walk alone and the list is ordered by
        ``(pathway, term_key)`` for determinism, independent of edge or
        parameter-file order. The graph-native order may change the
        presentation order of per-term/per-enzyme entries relative to the
        side-index order; callers that require a specific display order
        should sort explicitly. The records themselves are identical
        (field-for-field) to ``get_flux_reactions`` by construction --
        only the *source of the values* and the order differ.

        Fallbacks, all visible in the returned records' ``source``/params:

        - If the walk does not cover every side-index term (a partially
          baked or foreign graph), the roster keeps index membership and
          order, with walked records overlaying their index counterparts
          -- the graph never silently shrinks a roster below what the
          index serves.
        - Terms whose bindings have no edge yet and no group carrier
          keep their side-index record in that fallback mode.
        - With a bare engine (no graph loaded) there are no edges or
          carriers to walk, so this degrades to ``get_flux_reactions``
          unchanged.
        """
        cls = getattr(carcinogen_class, "value", carcinogen_class)
        roster = self.get_flux_reactions(cls)
        if not self.G.number_of_edges():
            return roster
        walked: dict[tuple[str, str], FluxReaction] = {}
        for _source_id, _target_id, edge_data in self.G.edges(data=True):
            kinetics = edge_data.get("kinetics")
            if not isinstance(kinetics, dict):
                continue
            flux_terms = kinetics.get("flux_terms")
            if not isinstance(flux_terms, dict):
                continue
            class_terms = flux_terms.get(cls)
            if not isinstance(class_terms, dict):
                continue
            for pathway, terms in class_terms.items():
                if not isinstance(terms, dict):
                    continue
                for term_key, payload in terms.items():
                    if not isinstance(payload, dict):
                        continue
                    if payload.get("dispatch_status") == "shadowed":
                        # Dispatch-shadowed payload (dual-source class: the
                        # kinetic parameterization carried as graph
                        # knowledge while the proxy model won dispatch).
                        # Not a member of the active roster; skipping keeps
                        # the walk semantically correct for a graph-native
                        # roster enumeration.
                        continue
                    params = {
                        key: value
                        for key, value in payload.items()
                        if key not in _FLUX_PAYLOAD_META_FIELDS
                    }
                    walked[(pathway, term_key)] = FluxReaction(
                        carcinogen_class=cls,
                        term_key=term_key,
                        pathway=pathway,
                        role=str(payload.get("role", "other")),
                        rate_law=str(payload.get("rate_law", "")),
                        params=params,
                        confidence=str(payload.get("confidence", "")),
                        provenance_ref=str(payload.get("provenance_ref", "")),
                        enzyme_id=payload.get("enzyme_id"),
                        graph_group=payload.get("graph_group"),
                        source=str(payload.get("source", "")),
                        sources=list(payload["sources"]) if isinstance(payload.get("sources"), list) else None,
                        notes=str(payload["notes"]) if payload.get("notes") is not None else None,
                    )
        # Graph-carried class-level terms (CarcinogenGroup
        # ``flux_class_metadata[cls].class_level_terms``): rebuild their
        # records from the group-node payload, mirroring the index build's
        # field resolution so the records are identical by construction --
        # only the source of the values differs (``source`` marks the
        # carrier so the reconstruction is auditable).
        carried_entry = self._flux_group_class_md().get(cls) or {}
        carried = carried_entry.get("class_level_terms")
        if isinstance(carried, dict) and carried:
            class_data = self._flux_proxy_class_cfg.get(cls, {})
            for section, terms in carried.items():
                if not isinstance(terms, dict):
                    continue
                for term_key, term in terms.items():
                    if not isinstance(term, dict):
                        continue
                    graph_node_id = term.get("graph_node_id")
                    enzyme_id: str | None = None
                    if graph_node_id:
                        enzyme_id = str(graph_node_id)
                    elif term_key in self.G:
                        enzyme_id = term_key
                    gene = term.get("gene")
                    if enzyme_id is None and isinstance(gene, str) and gene in self.G:
                        enzyme_id = gene
                    declared_group = class_data.get("graph_group")
                    walked[(section, term_key)] = FluxReaction(
                        carcinogen_class=cls,
                        term_key=term_key,
                        pathway=section,
                        role=_FLUX_ROLE_KINDS.get(section, "other"),
                        rate_law=str(
                            term.get("rate_law")
                            or term.get("equation")
                            or class_data.get("kinetics_model")
                            or ""
                        ),
                        params={
                            key: value
                            for key, value in term.items()
                            if key not in _FLUX_RESERVED_TERM_FIELDS
                        },
                        confidence=str(
                            term.get("confidence") or class_data.get("confidence_overall") or ""
                        ),
                        provenance_ref=str(
                            term.get("provenance_ref") or f"classes.{cls}.{section}.{term_key}"
                        ),
                        enzyme_id=enzyme_id,
                        graph_group=str(declared_group) if declared_group else _FLUX_CLASS_GRAPH_GROUPS.get(cls),
                        source="graph_class_level_carrier",
                        sources=list(term["sources"]) if isinstance(term.get("sources"), list) else None,
                        notes=str(term["notes"]) if term.get("notes") is not None else None,
                    )
        if not walked:
            return roster
        # Graph-native completeness gate: the walk serves as the roster on
        # its own only when it covers every side-index term. The reference
        # graph is at full carriage (verified by the roster diff), so the
        # loaded path is graph-native; a partially baked or foreign graph
        # falls back to the seeded mapping below rather than silently
        # shrinking the roster.
        index_keys = {(reaction.pathway, reaction.term_key) for reaction in roster}
        if index_keys <= walked.keys():
            return [walked[key] for key in sorted(walked)]
        return [walked.get((reaction.pathway, reaction.term_key), reaction) for reaction in roster]

    def get_flux_classes(self) -> list[str]:
        """Return the union of kinetic + proxy flux classes.

        Ordering follows file insertion order: kinetic classes in
        ``kinetic_parameters.json`` order first, then proxy-only classes in
        ``proxy_flux_parameters.json`` order (dual-source classes keep their
        kinetic position). Consumers that need the ``CarcinogenClass`` enum
        order (e.g. ``compute_full_profile`` output ordering) should apply it
        themselves rather than relying on this list.
        """
        self._ensure_flux_index()
        return list(self._flux_reactions_by_class)

    def _flux_group_class_md(self) -> dict[str, dict[str, Any]]:
        """Map flux class name -> its CarcinogenGroup node's metadata entry.

        Built lazily from the loaded reference graph and rebuilt whenever
        the underlying graph object changes (e.g. a reload), so class-level
        reads can prefer the graph over the bundled parameter files.
        Empty on a bare engine (no graph loaded) -- callers fall back to
        the JSON-sourced side index.
        """
        if self._flux_group_class_md_cache is None or self._flux_group_class_md_graph is not self.G:
            md: dict[str, dict[str, Any]] = {}
            for _, data in self.G.nodes(data=True):
                if data.get("type") != "CarcinogenGroup":
                    continue
                for cls, entry in (data.get("flux_class_metadata") or {}).items():
                    if isinstance(entry, dict):
                        md[str(cls)] = entry
            self._flux_group_class_md_cache = md
            self._flux_group_class_md_graph = self.G
        return self._flux_group_class_md_cache

    def _flux_class_signals(self) -> dict[str, dict[str, Any]]:
        """Map flux class name -> its edge-carried signal block.

        Built lazily from the loaded reference graph's baked
        ``kinetics["flux_class_signal"]`` payloads (see
        :meth:`_apply_flux_edge_kinetics` -- currently Dioxin's AhR
        occupancy Hill term on ``TCDD --AGONIZES--> AHR``) and rebuilt
        whenever the underlying graph object changes. Empty on a bare
        engine (no graph loaded) -- callers fall back to the proxy JSON's
        class config.
        """
        if self._flux_class_signal_cache is None or self._flux_class_signal_graph is not self.G:
            signals: dict[str, dict[str, Any]] = {}
            for _source, _target, data in self.G.edges(data=True):
                kinetics = data.get("kinetics")
                if not isinstance(kinetics, dict):
                    continue
                carried = kinetics.get("flux_class_signal")
                if not isinstance(carried, dict):
                    continue
                for cls_name, block in carried.items():
                    if isinstance(block, dict) and block:
                        signals[str(cls_name)] = block
            self._flux_class_signal_cache = signals
            self._flux_class_signal_graph = self.G
        return self._flux_class_signal_cache

    def get_flux_aggregation(self, carcinogen_class: Any) -> dict[str, Any]:
        """Return the per-class flux aggregation spec, graph-first.

        Prefers the aggregation block carried on the class's CarcinogenGroup
        node (``flux_class_metadata[cls].aggregation`` in the loaded
        reference graph) and falls back to the kinetic JSON's side index on
        a bare engine. The aggregation block names how a class's reaction
        terms combine: per-role vmax fields, derived terms, detox
        fractions, total scaling, rounding, and unit notes. Classes
        without a block (proxy classes, or mechanistic classes pending
        annotation) return an empty dict.
        """
        self._ensure_flux_index()
        cls = getattr(carcinogen_class, "value", carcinogen_class)
        entry = self._flux_group_class_md().get(cls)
        if entry and isinstance(entry.get("aggregation"), dict) and entry["aggregation"]:
            return entry["aggregation"]
        return self._flux_aggregation_by_class.get(cls, {})

    def get_flux_metadata(self) -> dict[str, Any]:
        """Return the metadata block of kinetic_parameters.json.

        Carries the doc-level parameter context shared across classes --
        exposure defaults (uM), qivive and steady-state tissue defaults --
        so consumers (e.g. the Aldehyde ethanol-oxidation substrate
        default) read it through the engine instead of the JSON file.
        """
        self._ensure_flux_index()
        return self._flux_metadata or {}

    def get_genotype_modifier(self, diplotype: str, gene: str) -> float:
        """Return Vmax scaling factor (0.0-2.0) based on metaboliser phenotype.

        Standard scale: PM=0.0, IM=0.5, NM=1.0, RM=1.5, UM=2.0.
        Special cases for ALDH2 heterozygotes, GSTM1/GSTT1 copy-number states,
        and cohort-specific CYP aliases used by the extension modules.

        Args:
            diplotype: Phenotype label (for example ``"NM"``, ``"null"``,
                or ``"*1/*2"``).
            gene: Gene name (for example ``"GSTM1"`` or ``"ALDH2"``).

        Returns:
            Scaling factor between 0.0 and 2.0.
        """
        self._ensure_flux_index()
        modifiers = self._flux_genotype_modifiers or {}
        special = modifiers["special_cases"]
        std = modifiers["standard_scale"]

        diplotype_lower = diplotype.lower().strip()
        gene_upper = gene.upper().strip()

        # Gene-specific manuscript/reference aliases that are more precise than
        # the generic PM/IM/NM/RM/UM scale.
        if gene_upper == "CYP1A2":
            if diplotype_lower in ("*1f/*1f", "1f/1f", "cyp1a2*1f/*1f", "um_1f_1f"):
                return 1.5
            if diplotype_lower in ("*1a/*1f", "1a/1f", "*1f/*1a", "1f/1a"):
                # Heterozygous *1F: intermediate inducibility (Sachse et al. 1999;
                # Ghotbi et al. 2007). Splits the difference between *1A/*1A (1.0)
                # and *1F/*1F (1.5).
                return 1.25
            if diplotype_lower in ("*1a/*1a", "1a/1a"):
                return 1.0
            if diplotype_lower in ("*1k", "*1k/*1k", "1k/1k"):
                return 0.5
            if diplotype_lower in ("pm", "poor", "poor metabolizer", "poor_metabolizer"):
                return 0.3

        if gene_upper == "NAT2":
            if diplotype_lower in ("slow", "slow acetylator", "slow_acetylator", "sa"):
                return 0.2
            if diplotype_lower in ("intermediate", "intermediate acetylator", "intermediate_acetylator"):
                return 0.5
            if diplotype_lower in ("rapid", "rapid acetylator", "rapid_acetylator", "ra"):
                return 1.0

        if gene_upper == "NAT1":
            # NAT1 phenotype assignments: *4 is the reference (rapid), *10 has
            # been associated with modestly increased acetylation activity in
            # several reports (Bell et al. 1995; Lin et al. 1998), and *14 is a
            # well-established slow allele (Hughes et al. 1998).
            if diplotype_lower in ("*4/*4", "4/4", "rapid", "ra"):
                return 1.0
            if diplotype_lower in ("*4/*10", "*10/*4", "4/10", "10/4"):
                return 1.05
            if diplotype_lower in ("*10/*10", "10/10"):
                return 1.1
            if diplotype_lower in ("*4/*14", "*14/*4", "4/14", "14/4"):
                return 0.75
            if diplotype_lower in ("*10/*14", "*14/*10", "10/14", "14/10"):
                return 0.8
            if diplotype_lower in ("*14/*14", "14/14", "slow", "sa"):
                return 0.5

        if gene_upper == "CYP2D6":
            if diplotype_lower in ("*1/*1", "*1/*2", "*2/*2"):
                return 1.0
            if diplotype_lower in ("*1/*4", "*1/*5", "*2/*4", "*10/*10", "im"):
                return 0.5
            if diplotype_lower in ("*4/*4", "*5/*5", "*4/*5", "pm", "poor"):
                return 0.0
            if "x2" in diplotype_lower or diplotype_lower in ("um", "ultrarapid"):
                return 2.0

        # Special cases
        if gene_upper == "ALDH2" and diplotype in ("*1/*2", "heterozygote"):
            return float(special["ALDH2_star1_star2"]["activity_fraction"])

        if gene_upper == "ALDH2" and diplotype in ("*2/*2", "PM_ALDH2"):
            return float(
                special.get("ALDH2_star2_homozygous", {}).get("activity_fraction", 0.001)
            )

        if gene_upper == "GSTM1" and diplotype_lower in (
            "null",
            "null/null",
            "deletion",
            "deleted",
            "0",
            "0/0",
        ):
            return _GST_NULL_RESIDUAL_ACTIVITY

        if gene_upper == "GSTT1" and diplotype_lower in (
            "null",
            "null/null",
            "deletion",
            "deleted",
            "0",
            "0/0",
        ):
            return _GST_NULL_RESIDUAL_ACTIVITY

        if gene_upper in {"GSTM1", "GSTT1"} and diplotype_lower in (
            "present",
            "active",
            "wt",
            "wildtype",
            "*1/*1",
            "1/1",
        ):
            return 1.0

        # Extension cohort aliases.
        if gene_upper == "CYP2E1" and diplotype_lower in ("um_c1c1", "*1c/*1c", "c1/c1"):
            # Source interaction model defines CYP2E1*1C/*1C as a 140% activity state.
            return 1.4

        if gene_upper == "CYP1A1":
            if diplotype_lower in ("*1/*2a", "*1/2a", "wt/*2a", "*2a carrier"):
                # CYP1A1*2A is an inducibility/risk allele rather than a well-calibrated
                # kinetic phenotype, so use a conservative step-up above NM.
                return 1.25
            if diplotype_lower in ("*2a/*2a", "2a/2a"):
                return 1.5

        # Gene-specific phenotype scales (override the generic PM=0/IM=0.5/NM=1.0
        # standard for genes whose null/slow phenotype retains substantial residual
        # activity in vivo).
        if gene_upper == "EPHX1":
            # Hassett 1994; Smith 1997: Y113H/Y113H ("slow") retains ~30-50% epoxide
            # hydrolase activity, not zero.
            if diplotype_lower in ("pm", "slow"):
                return 0.4
            if diplotype_lower in ("im", "intermediate"):
                return 0.7
            if diplotype_lower in ("rm", "rapid", "fast"):
                return 1.3

        if gene_upper == "NQO1":
            # Siegel 1999; Ross 2004: NQO1*2 (Pro187Ser) homozygotes retain ~3-5%
            # activity due to ubiquitin-mediated degradation; heterozygotes ~50%.
            if diplotype_lower in ("pm", "*2/*2"):
                return 0.05
            if diplotype_lower in ("im", "*1/*2"):
                return 0.5

        if gene_upper == "GSTP1":
            # Watson 1998; Hu 1997: Ile105Val (Val/Val) retains ~30-50% activity
            # for many PAH-diol epoxide substrates (reduced thermal stability and
            # affinity, not loss of function).
            if diplotype_lower in ("pm", "val/val"):
                return 0.4
            if diplotype_lower in ("im", "ile/val"):
                return 0.7

        if gene_upper == "CYP1B1":
            # Bailey 1998; Shimada 1999: CYP1B1*3 (L432V) carriers have modestly
            # elevated catalysis (~25-50%), not the standard 2x ultrarapid scale.
            if diplotype_lower in ("rm", "*1/*3", "leu/val"):
                return 1.25
            if diplotype_lower in ("um", "*3/*3", "val/val"):
                return 1.5

        # Standard phenotype scale
        phenotype_map = {
            "pm": std["PM"],
            "poor": std["PM"],
            "im": std["IM"],
            "intermediate": std["IM"],
            "nm": std["NM"],
            "normal": std["NM"],
            "wt": std["NM"],
            "wildtype": std["NM"],
            "*1/*1": std["NM"],
            "rm": std["RM"],
            "rapid": std["RM"],
            "um": std["UM"],
            "ultrarapid": std["UM"],
            "null": std["PM"],
            "0/0": 0.0,
        }

        key = diplotype_lower
        if key in phenotype_map:
            return float(phenotype_map[key])

        # Try numeric
        try:
            val = float(diplotype)
            if 0.0 <= val <= 2.0:
                return val
        except (ValueError, TypeError):
            pass

        warnings.warn(
            f"Unrecognized diplotype '{diplotype}' for gene '{gene}'; defaulting to NM (1.0)",
            stacklevel=2,
        )
        return 1.0

    def get_proxy_genotype_modifier(self, diplotype: str, gene: str | None) -> float:
        """Return a silent genotype scaling factor for proxy-model terms."""
        if not gene:
            return 1.0

        gene_upper = gene.upper().strip()
        diplotype_lower = str(diplotype).lower().strip()
        self._ensure_flux_index()
        modifiers = self._flux_genotype_modifiers or {}
        special = modifiers["special_cases"]
        std = modifiers["standard_scale"]

        if gene_upper == "CYP1A2":
            if diplotype_lower in ("*1f/*1f", "1f/1f", "cyp1a2*1f/*1f", "um_1f_1f"):
                return 1.5
            if diplotype_lower in ("*1a/*1f", "1a/1f", "*1f/*1a", "1f/1a"):
                return 1.25
            if diplotype_lower in ("*1a/*1a", "1a/1a"):
                return 1.0
            if diplotype_lower in ("*1k", "*1k/*1k", "1k/1k"):
                return 0.5
            if diplotype_lower in ("pm", "poor", "poor metabolizer", "poor_metabolizer"):
                return 0.3

        if gene_upper == "NAT2":
            if diplotype_lower in ("slow", "slow acetylator", "slow_acetylator", "sa"):
                return 0.2
            if diplotype_lower in ("intermediate", "intermediate acetylator", "intermediate_acetylator"):
                return 0.5
            if diplotype_lower in ("rapid", "rapid acetylator", "rapid_acetylator", "ra"):
                return 1.0

        if gene_upper == "NAT1":
            if diplotype_lower in ("*4/*4", "4/4", "rapid", "ra"):
                return 1.0
            if diplotype_lower in ("*4/*10", "*10/*4", "4/10", "10/4"):
                return 1.05
            if diplotype_lower in ("*10/*10", "10/10"):
                return 1.1
            if diplotype_lower in ("*4/*14", "*14/*4", "4/14", "14/4"):
                return 0.75
            if diplotype_lower in ("*10/*14", "*14/*10", "10/14", "14/10"):
                return 0.8
            if diplotype_lower in ("*14/*14", "14/14", "slow", "sa"):
                return 0.5

        if gene_upper == "CYP2D6":
            if diplotype_lower in ("*1/*1", "*1/*2", "*2/*2"):
                return 1.0
            if diplotype_lower in ("*1/*4", "*1/*5", "*2/*4", "*10/*10", "im"):
                return 0.5
            if diplotype_lower in ("*4/*4", "*5/*5", "*4/*5", "pm", "poor"):
                return 0.0
            if "x2" in diplotype_lower or diplotype_lower in ("um", "ultrarapid"):
                return 2.0

        if gene_upper == "ALDH2" and diplotype in ("*1/*2", "heterozygote"):
            return float(special["ALDH2_star1_star2"]["activity_fraction"])
        if gene_upper == "ALDH2" and diplotype in ("*2/*2", "PM_ALDH2"):
            return float(
                special.get("ALDH2_star2_homozygous", {}).get("activity_fraction", 0.001)
            )
        if gene_upper in {"GSTM1", "GSTT1"} and diplotype_lower in (
            "null",
            "null/null",
            "deletion",
            "deleted",
            "0",
            "0/0",
        ):
            return _GST_NULL_RESIDUAL_ACTIVITY
        if gene_upper in {"GSTM1", "GSTT1"} and diplotype_lower in (
            "present",
            "active",
            "wt",
            "wildtype",
            "*1/*1",
            "1/1",
        ):
            return 1.0
        if gene_upper == "CYP2E1" and diplotype_lower in ("um_c1c1", "*1c/*1c", "c1/c1"):
            return 1.4
        if gene_upper == "CYP1A1":
            if diplotype_lower in ("*1/*2a", "*1/2a", "wt/*2a", "*2a carrier"):
                return 1.25
            if diplotype_lower in ("*2a/*2a", "2a/2a"):
                return 1.5

        if gene_upper == "EPHX1":
            if diplotype_lower in ("pm", "slow"):
                return 0.4
            if diplotype_lower in ("im", "intermediate"):
                return 0.7
            if diplotype_lower in ("rm", "rapid", "fast"):
                return 1.3
        if gene_upper == "NQO1":
            if diplotype_lower in ("pm", "*2/*2"):
                return 0.05
            if diplotype_lower in ("im", "*1/*2"):
                return 0.5
        if gene_upper == "GSTP1":
            if diplotype_lower in ("pm", "val/val"):
                return 0.4
            if diplotype_lower in ("im", "ile/val"):
                return 0.7
        if gene_upper == "CYP1B1":
            if diplotype_lower in ("rm", "*1/*3", "leu/val"):
                return 1.25
            if diplotype_lower in ("um", "*3/*3", "val/val"):
                return 1.5

        phenotype_map = {
            "pm": std["PM"],
            "poor": std["PM"],
            "im": std["IM"],
            "intermediate": std["IM"],
            "nm": std["NM"],
            "normal": std["NM"],
            "wt": std["NM"],
            "wildtype": std["NM"],
            "*1/*1": std["NM"],
            "rm": std["RM"],
            "rapid": std["RM"],
            "um": std["UM"],
            "ultrarapid": std["UM"],
            "null": std["PM"],
            "0/0": 0.0,
        }
        if diplotype_lower in phenotype_map:
            return float(phenotype_map[diplotype_lower])
        try:
            val = float(diplotype)
            if 0.0 <= val <= 40.0:
                return val
        except (ValueError, TypeError):
            pass
        return 1.0


        # ── Tissue weight (GTEx integration) ──────────────────────────────────────

        _FLUX_TISSUE_TO_GTEX: dict[str, str] = {
        "liver": "Liver",
        "lung": "Lung",
        "prostate": "Prostate",
        "bladder": "Bladder",
        "colon": "Colon",
        "breast": "Breast",
        "kidney": "Kidney",
        "esophagus": "Esophagus",
        }

        _TISSUE_ALIASES: dict[str, str] = {
        "liver": "liver",
        "hepatic": "liver",
        "lung": "lung",
        "pulmonary": "lung",
        "prostate": "prostate",
        "breast": "breast",
        "colon": "colon",
        "colorectal": "colon",
        "kidney": "kidney",
        "renal": "kidney",
        "bladder": "bladder",
        "lymphocyte": "lymphocyte",
        "blood": "lymphocyte",
        "esophagus": "esophagus",
        "esophageal": "esophagus",
        "stomach": "stomach",
        "gastric": "stomach",
        "intestine": "intestine",
        "nasal_mucosa": "nasal_mucosa",
        "brain": "brain",
        "heart": "heart",
        "muscle": "muscle",
        "adipose": "adipose",
        "skin": "skin",
        "placenta": "placenta",
        }

    def get_qivive_context(
        self,
        tissue_key: str,
        overrides: Mapping[str, float] | None = None,
    ) -> dict[str, float]:
        """Return MPPGL/organ-weight context for optional QIVIVE flux scaling.

        ``tissue_key`` is the caller-normalized tissue key; falls back to
        reference physiology when the kinetic metadata omits a tissue.
        """
        self._ensure_flux_index()
        metadata = self.get_flux_metadata()
        qivive_defaults = metadata.get("qivive_defaults", {})
        tissue_defaults = qivive_defaults.get("tissues", {})
        source = tissue_defaults.get(tissue_key, _FALLBACK_QIVIVE_TISSUES.get(tissue_key, _FALLBACK_QIVIVE_TISSUES["liver"]))
        context = {
            "mppgl_mg_per_g": float(source["mppgl_mg_per_g"]),
            "organ_weight_g": float(source["organ_weight_g"]),
        }
        if overrides:
            if "mppgl_mg_per_g" in overrides:
                context["mppgl_mg_per_g"] = float(overrides["mppgl_mg_per_g"])
            if "organ_weight_g" in overrides:
                context["organ_weight_g"] = float(overrides["organ_weight_g"])
        context["scale"] = round(context["mppgl_mg_per_g"] * context["organ_weight_g"], 6)
        return context

    def get_steady_state_context(
        self,
        tissue_key: str,
        overrides: Mapping[str, float] | None = None,
    ) -> dict[str, float]:
        """Return validated defaults for the flux-coupled steady-state solver.

        ``tissue_key`` is the caller-normalized tissue key.
        """
        self._ensure_flux_index()
        metadata = self.get_flux_metadata()
        configured = metadata.get("steady_state_defaults", {})

        qivive_context = self.get_qivive_context(tissue_key)
        fallback_tissue = _FALLBACK_STEADY_STATE_TISSUES.get(
            tissue_key,
            _FALLBACK_STEADY_STATE_TISSUES["liver"],
        )
        configured_tissues = configured.get("tissues", {})
        configured_tissue = configured_tissues.get(tissue_key, {})
        tissue_source: dict[str, Any] = {
            **fallback_tissue,
            **configured_tissue,
            "organ_weight_g": configured_tissue.get(
                "organ_weight_g",
                qivive_context.get("organ_weight_g", fallback_tissue["organ_weight_g"]),
            ),
        }

        context: dict[str, float] = {}
        for key, fallback in _STEADY_STATE_DEFAULTS.items():
            context[key] = _positive_context_float(configured, key, fallback)
        context["absorption_fraction"] = _bounded_fraction_context(
            configured,
            "absorption_fraction",
            _STEADY_STATE_DEFAULTS["absorption_fraction"],
        )
        context["organ_weight_g"] = _positive_context_float(
            tissue_source,
            "organ_weight_g",
            fallback_tissue["organ_weight_g"],
        )
        context["tissue_partition_coefficient"] = _positive_context_float(
            tissue_source,
            "tissue_partition_coefficient",
            fallback_tissue["tissue_partition_coefficient"],
        )
        context["tissue_blood_flow_fraction"] = _bounded_fraction_context(
            tissue_source,
            "tissue_blood_flow_fraction",
            fallback_tissue["tissue_blood_flow_fraction"],
        )

        if overrides:
            for key, value in overrides.items():
                if key in {"absorption_fraction", "tissue_blood_flow_fraction"}:
                    context[key] = _bounded_fraction_context(overrides, key, context[key])
                else:
                    context[key] = _positive_context_float(overrides, key, context.get(key, 1.0))

        context["central_volume_l"] = round(
            context["body_weight_kg"] * context["volume_l_per_kg"],
            6,
        )
        context["tissue_volume_l"] = round(context["organ_weight_g"] / 1000.0, 6)
        context["tissue_blood_flow_l_per_day"] = round(
            context["cardiac_output_l_per_day"] * context["tissue_blood_flow_fraction"],
            6,
        )
        return context

    def get_default_concentration(self, carcinogen_class: Any) -> float:
        """Return default environmental exposure concentration in uM.

        Kinetic classes read ``metadata.exposure_defaults_uM`` through the
        ``exposure_default_field`` named in their aggregation block; proxy
        classes through their retained ``exposure_default`` spec (a kinetic
        parameters field or an exposure-database scenario). Other classes
        fall back to 0.1 uM.
        """
        self._ensure_flux_index()
        cls = getattr(carcinogen_class, "value", carcinogen_class)
        metadata = self.get_flux_metadata()
        defaults = metadata.get("exposure_defaults_uM", {})
        field_name = self.get_flux_aggregation(cls).get("exposure_default_field")
        if field_name is not None:
            value = defaults.get(field_name)
            if value is not None:
                return float(value)

        spec = self._flux_proxy_exposure_defaults.get(cls)
        entry = self._flux_group_class_md().get(cls)
        if entry and isinstance(entry.get("exposure_defaults"), dict):
            spec = entry["exposure_defaults"]
        if isinstance(spec, dict):
            source = spec.get("source")
            if source == "kinetic_parameters":
                value = defaults.get(spec.get("field"))
                if value is not None:
                    return float(value)
            elif source == "exposure_database":
                exposure_db = self._ensure_exposure_database()
                scenario = (
                    exposure_db.get("carcinogen_classes", {})
                    .get(spec.get("class"), {})
                    .get("exposure_scenarios", {})
                    .get(spec.get("scenario"), {})
                )
                value = scenario.get(spec.get("field"))
                if value is not None:
                    return float(value)
        return 0.1

    def _ensure_exposure_database(self) -> dict[str, Any]:
        """Lazily load the packaged exposure database for proxy defaults."""
        if self._exposure_database is None:
            try:
                self._exposure_database = json.loads(
                    _DEFAULT_EXPOSURE_DB_PATH.read_text(encoding="utf-8")
                )
            except (OSError, json.JSONDecodeError) as exc:
                raise RuntimeError(
                    f"Could not read exposure database at {_DEFAULT_EXPOSURE_DB_PATH}: {exc}"
                ) from exc
        return self._exposure_database


    def get_flux_class_config(self, carcinogen_class: Any) -> dict[str, Any]:
        """Return a class's proxy config block, class-level fields graph-first.

        Empty dict for classes that only have kinetic-parameters entries
        (this emptiness is load-bearing -- callers use it to distinguish
        measured-kinetics from proxy classes). For proxy classes the
        class-level fields (``model_kind``, ``exposure_default``) are
        served from the class's CarcinogenGroup node when the reference
        graph is loaded, falling back to the proxy JSON on a bare engine;
        non-enzyme class-level terms (driver/proxy terms with no edge
        binding, e.g. ``general_ROS``) are likewise served graph-first
        from the group node's ``flux_class_metadata[cls].class_level_terms``
        carrier, falling back to the proxy JSON; enzyme terms keep
        coming from the proxy JSON until the anchoring backlog gives
        them graph edges. A class's ``signal`` block (receptor/
        signaling-model classes, e.g. Dioxin) is served graph-first
        from its baked edge payload when the block carries a binding
        (``TCDD --AGONIZES--> AHR``), falling back to the proxy JSON on a
        bare engine. Also serves unit notes and signal configuration for
        the proxy flux classes.
        """
        self._ensure_flux_index()
        cls = getattr(carcinogen_class, "value", carcinogen_class)
        cfg = dict(self._flux_proxy_class_cfg.get(cls, {}))
        if cfg:
            entry = self._flux_group_class_md().get(cls)
            if entry:
                if entry.get("model_kind"):
                    cfg["model_kind"] = entry["model_kind"]
                if isinstance(entry.get("exposure_defaults"), dict):
                    cfg["exposure_default"] = dict(entry["exposure_defaults"])
                carried = entry.get("class_level_terms")
                if isinstance(carried, dict):
                    for section, terms in carried.items():
                        if not isinstance(terms, dict) or not terms:
                            continue
                        merged = dict(cfg.get(section, {}))
                        for term_key, term in terms.items():
                            if isinstance(term, dict):
                                merged[term_key] = dict(term)
                        cfg[section] = merged
        signal_block = self._flux_class_signals().get(cls)
        if isinstance(signal_block, dict) and signal_block:
            cfg["signal"] = dict(signal_block)
        return cfg

    def get_flux_provenance_entry(self, ref: str) -> dict[str, Any]:
        """Resolve a dotted ``provenance_ref`` into proxy_flux_provenance.json."""
        node: Any = self._ensure_flux_proxy_provenance()
        for part in ref.split("."):
            node = node[part]
        return cast(dict[str, Any], node)

    def _ensure_flux_proxy_provenance(self) -> dict[str, Any]:
        """Lazily load the packaged proxy flux provenance doc."""
        if self._flux_proxy_provenance is None:
            try:
                self._flux_proxy_provenance = json.loads(
                    _DEFAULT_PROXY_FLUX_PROVENANCE_PATH.read_text(encoding="utf-8")
                )
            except (OSError, json.JSONDecodeError) as exc:
                raise RuntimeError(
                    f"Could not read proxy flux provenance at {_DEFAULT_PROXY_FLUX_PROVENANCE_PATH}: {exc}"
                ) from exc
        return self._flux_proxy_provenance

    def get_flux_reaction_coverage(self, carcinogen_class: Any) -> dict[str, Any]:
        """Compare a class's flux-reaction roster against the graph.

        Two anchoring criteria are reported per roster term:

        - **Binding edges (primary).** A term carrying both
          ``source_node_id`` and ``target_node_id`` is
          edge-anchored when an edge between those two nodes carries its
          baked flux
          payload -- ``kinetics.flux_terms[cls][pathway][term_key]``, the
          exact payload :meth:`_apply_flux_edge_kinetics` bakes and
          :meth:`get_edge_flux_reactions` reads. Edge presence alone
          (``binding_edge_present``) is necessary but not sufficient:
          terms whose binding edge exists but carries no payload are
          ``missing_flux_payloads``; terms with no bindings at all are
          ``unbound_terms`` -- they need substrate-oriented annotations,
          and most proxy/semi-quantitative rosters land there by design.
          Terms whose binding edge does not exist are the genuine
          curation backlog (``missing_binding_edges``).
        - **Class-carcinogen scope edges (secondary, legacy).** Whether
          any scope edge (``SUBSTRATE_OF`` / ``DETOXIFIED_BY`` /
          ``REPAIRED_BY``, either direction) connects the class's
          Carcinogen nodes to the enzyme. Metabolite-anchored bindings
          (e.g. ``AFB1_epoxide → GSTA1``) legitimately have no such edge,
          so this is a diagnostic, not the backlog.

        For dual-source classes (ChlorinatedSolvent, Dioxin, HeavyMetal)
        the roster carries the proxy entries; the shadowed kinetic
        entries are reported separately in ``shadowed_reactions`` with
        the same anchoring fields.

        ``missing_scope_edges`` is kept as an alias of
        ``missing_binding_edges`` for the established vocabulary.
        ``edge_anchored_count`` is an independent axis from the
        unresolved/unbound/missing taxonomy: an alias pseudo-term (e.g.
        ``CYP2E1_with_b5``) can be edge-anchored via its binding while
        still reporting an unresolved legacy ``enzyme_id``.
        Requires the reference graph to be loaded to say anything
        meaningful; on a bare engine every term reports unanchored.
        This getter is the curation-backlog driver, not a pass/fail
        check.
        """
        self._ensure_flux_index()
        cls = getattr(carcinogen_class, "value", carcinogen_class)
        reactions = self._flux_reactions_by_class.get(cls, [])
        graph_group = next((r.graph_group for r in reactions if r.graph_group), None) or _FLUX_CLASS_GRAPH_GROUPS.get(cls)
        carcinogen_ids = []
        if graph_group:
            carcinogen_ids = [
                node_id
                for node_id, data in self.G.nodes(data=True)
                if data.get("type") == "Carcinogen" and data.get("group") == graph_group
            ]

        def _binding_probe(reaction: FluxReaction) -> tuple[bool, bool, list[str]] | None:
            """``(edge_present, payload_present, edge_types)`` for the term's
            binding nodes; ``None`` when the term carries no bindings.

            ``payload_present`` requires an edge between the binding nodes
            to carry ``kinetics.flux_terms[cls][pathway][term_key]`` -- the
            exact payload the edge-walk reader consumes. Edge presence
            alone is necessary but not sufficient.
            """
            substrate = reaction.params.get("source_node_id")
            enzyme = reaction.params.get("target_node_id")
            if not substrate or not enzyme:
                return None
            if not self.G.has_edge(substrate, enzyme):
                return (False, False, [])
            edge_types: set[str] = set()
            payload_present = False
            for data in self.G.get_edge_data(substrate, enzyme).values():
                if data.get("type"):
                    edge_types.add(str(data.get("type")))
                flux_terms = (data.get("kinetics") or {}).get("flux_terms", {})
                if reaction.term_key in flux_terms.get(cls, {}).get(reaction.pathway, {}):
                    payload_present = True
            return (True, payload_present, sorted(edge_types))

        def _class_carcinogen_edge_types(enzyme_id: str | None) -> list[str]:
            """Legacy scope-edge diagnostic between class carcinogens and the enzyme."""
            if not enzyme_id or not carcinogen_ids:
                return []
            edge_types: set[str] = set()
            for carcinogen_id in carcinogen_ids:
                for source_id, target_id in ((carcinogen_id, enzyme_id), (enzyme_id, carcinogen_id)):
                    if not self.G.has_edge(source_id, target_id):
                        continue
                    for data in self.G.get_edge_data(source_id, target_id).values():
                        if data.get("type") in _FLUX_SCOPE_EDGE_TYPES:
                            edge_types.add(str(data.get("type")))
            return sorted(edge_types)

        # Terms whose config blocks are carried on the class's
        # CarcinogenGroup node (``flux_class_metadata[cls].class_level_terms``
        # -- graph-first serving with the proxy JSON as fallback). These
        # stay class-level by design: they have no edge binding, so they are
        # excluded from the unresolved/unbound/missing backlog below and
        # reported separately as a resolved carrier state.
        carried_terms: set[str] = set()
        carried_entry = self._flux_group_class_md().get(cls) or {}
        for _section, _terms in (carried_entry.get("class_level_terms") or {}).items():
            if isinstance(_terms, dict):
                carried_terms.update(_terms)

        def _row(reaction: FluxReaction) -> dict[str, Any]:
            resolved = reaction.enzyme_id is not None and reaction.enzyme_id in self.G
            probe = _binding_probe(reaction)
            edge_present, payload_present, binding_types = probe if probe else (False, False, [])
            return {
                "term_key": reaction.term_key,
                "pathway": reaction.pathway,
                "role": reaction.role,
                "rate_law": reaction.rate_law,
                "enzyme_id": reaction.enzyme_id,
                "source_node_id": reaction.params.get("source_node_id"),
                "target_node_id": reaction.params.get("target_node_id"),
                "resolved": resolved,
                "binding_edge_present": edge_present,
                "binding_edge_types": binding_types,
                "edge_anchored": payload_present,
                "group_carried": reaction.term_key in carried_terms,
                "class_carcinogen_edge_types": _class_carcinogen_edge_types(reaction.enzyme_id),
            }

        rows: list[dict[str, Any]] = []
        unresolved: list[str] = []
        unbound: list[str] = []
        missing_binding: list[str] = []
        missing_payload: list[str] = []
        for reaction in reactions:
            row = _row(reaction)
            if not row["resolved"] and reaction.term_key not in carried_terms:
                unresolved.append(reaction.term_key)
            rows.append(row)

        # Classify the resolved terms: unbound = no binding annotations;
        # missing_binding = bindings whose edge doesn't exist;
        # missing_payload = binding edge exists but carries no baked
        # flux_terms payload for this term. Group-carried terms are a
        # resolved carrier state and are excluded from all three.
        for reaction, row in zip(reactions, rows):
            if not row["resolved"] or reaction.term_key in carried_terms:
                continue
            if not reaction.params.get("source_node_id") or not reaction.params.get("target_node_id"):
                unbound.append(reaction.term_key)
            elif not row["binding_edge_present"]:
                missing_binding.append(reaction.term_key)
            elif not row["edge_anchored"]:
                missing_payload.append(reaction.term_key)

        shadowed = [r for r in self._flux_shadowed_reactions if r.carcinogen_class == cls]
        shadowed_rows = [_row(reaction) for reaction in shadowed]

        group_carried = [r.term_key for r in reactions if r.term_key in carried_terms]

        return {
            "carcinogen_class": cls,
            "graph_group": graph_group,
            "carcinogen_node_ids": carcinogen_ids,
            "sources": sorted(self._flux_class_sources.get(cls, set())),
            "reaction_count": len(reactions),
            "resolved_enzyme_count": len(reactions) - len(unresolved),
            "edge_anchored_count": sum(1 for row in rows if row["edge_anchored"]),
            "unresolved_terms": unresolved,
            "unbound_terms": unbound,
            "group_carried_terms": group_carried,
            "missing_binding_edges": missing_binding,
            "missing_scope_edges": missing_binding,
            "missing_flux_payloads": missing_payload,
            "shadowed_reactions": shadowed_rows,
            "reactions": rows,
        }

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

        e.g. ``"Aldehydes"``, ``"UV Radiation"``, ``"PFAS"``. Scoped to
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

    #: Node types hidden from viewer-facing subgraphs by default (they
    #: remain in the underlying graph for callers that bypass this method,
    #: e.g. the flux engine) -- see :meth:`subgraph_by_node_types`.
    _DEFAULT_EXCLUDED_VIEWER_NODE_TYPES: tuple[str, ...] = ("Substrate", "CarcinogenGroup")

    def subgraph_by_node_types(
        self,
        node_types: Sequence[str] | None = None,
        *,
        exclude_types: Sequence[str] = _DEFAULT_EXCLUDED_VIEWER_NODE_TYPES,
    ) -> dict[str, list[Any]]:
        """Return the subgraph restricted to *node_types*, minus *exclude_types*.

        Generalizes :meth:`subgraph_by_node_type` to a multi-select: *node_types*
        of ``None`` or empty imposes no restriction (every type is included) --
        the identity element for this axis, mirroring the ``None``-means-
        unrestricted convention :meth:`filtered_subgraph` already uses.
        *exclude_types* is always subtracted afterward regardless of whether it
        appears in *node_types*, defaulting to hiding ``Substrate`` nodes from
        viewer-facing subgraphs -- pass ``exclude_types=()`` to see them.

        Edges are included only when *both* endpoints survive the type
        filter; an edge's separate ``carcinogen`` context field is not
        itself a node-type constraint.
        """
        excluded = set(exclude_types)
        allowed = set(node_types) if node_types else None
        nodes = [
            dict(data)
            for _, data in self.G.nodes(data=True)
            if data.get("type") not in excluded
            and (allowed is None or data.get("type") in allowed)
        ]
        kept_ids = {node["id"] for node in nodes}
        edges = [
            dict(data)
            for u, v, data in self.G.edges(data=True)
            if u in kept_ids and v in kept_ids
        ]
        return {"nodes": nodes, "edges": edges}

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
        requesting ``group="PAHs"`` and ``node_type="Enzyme"`` here returns
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

    # ── Path traversal ───────────────────────────────────────────────────

    def _maximal_simple_paths(
        self,
        start: str,
        *,
        forward: bool,
        edge_types: Sequence[str] | None,
    ) -> list[dict[str, list[Any]]]:
        """Enumerate every *maximal* simple directed path touching *start*.

        Shared DFS core for :meth:`paths_from_carcinogen` (``forward=True``,
        walking successors) and :meth:`paths_to_carcinogen` (``forward=False``,
        walking predecessors). A path is "simple" in the graph-theory sense
        (no repeated node) and "maximal" in that it is extended greedily
        until every next step would either revisit an already-visited node
        in *this* path or there is no further edge to take -- at which point
        the path is recorded and that branch stops. Simple-path semantics is
        what makes this well-defined at all: this graph is not acyclic (e.g.
        ``MeHg_GSH``/``MethylmercuryCompounds`` form a 2-cycle via paired
        ``ACTIVATES``/``DETOXIFIES`` edges), so an unbounded walk could
        otherwise recurse forever. A would-revisit step is treated as a dead
        end for that branch, not an error.

        *edge_types* optionally restricts which edges are walkable (e.g. to
        exclude ``PATHWAY`` edges into generic KEGG pathway-annotation nodes
        rather than genuine biochemical transformation steps); ``None``
        (the default) walks every edge type.

        Each returned path is ``{"nodes": [...], "edges": [...]}`` with
        nodes/edges in traversal order from *start* onward -- ``edges[i]``
        connects ``nodes[i]`` to ``nodes[i + 1]`` -- regardless of *forward*,
        so callers never need to know which direction produced a given
        path. Nodes/edges are full attribute dicts (matching :meth:`to_dict`),
        not bare ids. A graph can easily have many branching maximal paths
        from one node, and the same edge can legitimately appear in several
        returned paths -- callers that want a deduplicated edge set across
        many paths should union on each edge dict's ``(source, target)``.
        Excludes the degenerate zero-edge "path" of *start* alone.
        """
        allowed_types = set(edge_types) if edge_types is not None else None
        paths: list[dict[str, list[Any]]] = []

        def _step(node: str) -> list[tuple[str, dict[str, Any]]]:
            view = (
                self.G.out_edges(node, data=True)
                if forward
                else self.G.in_edges(node, data=True)
            )
            steps = []
            for u, v, data in view:
                if allowed_types is not None and data.get("type") not in allowed_types:
                    continue
                neighbor = v if forward else u
                steps.append((neighbor, data))
            return steps

        def _record(node_ids: list[str], edge_data: list[dict[str, Any]]) -> None:
            if not edge_data:
                return
            ordered_nodes = node_ids if forward else list(reversed(node_ids))
            ordered_edges = edge_data if forward else list(reversed(edge_data))
            paths.append(
                {
                    "nodes": [dict(self.G.nodes[n]) for n in ordered_nodes],
                    "edges": [dict(e) for e in ordered_edges],
                }
            )

        def _dfs(node: str, visited: set[str], node_ids: list[str], edge_data: list[dict[str, Any]]) -> None:
            extended = False
            for neighbor, data in _step(node):
                if neighbor in visited:
                    continue
                extended = True
                visited.add(neighbor)
                node_ids.append(neighbor)
                edge_data.append(data)
                _dfs(neighbor, visited, node_ids, edge_data)
                edge_data.pop()
                node_ids.pop()
                visited.discard(neighbor)
            if not extended:
                _record(node_ids, edge_data)

        _dfs(start, {start}, [start], [])
        return paths

    def paths_from_carcinogen(
        self,
        carcinogen_id: str,
        *,
        edge_types: Sequence[str] | None = None,
    ) -> list[dict[str, list[Any]]]:
        """Return every maximal simple directed path starting at *carcinogen_id*.

        *carcinogen_id* must name an existing node with ``type == "Carcinogen"``
        (raises ``ValueError`` otherwise, mirroring :meth:`add_edge`'s
        validation style). Each path begins at *carcinogen_id* and walks
        forward (successors) until it dead-ends -- see
        :meth:`_maximal_simple_paths` for the exact semantics (simple-path
        cycle handling, *edge_types* filtering, and the returned shape).

        Complements :meth:`paths_to_carcinogen`. A path returned here can
        itself end at a different Carcinogen node (e.g. a precursor chain
        that resolves into another named carcinogen), which is expected,
        not an error.
        """
        if carcinogen_id not in self.G:
            raise ValueError(f"Unknown node: {carcinogen_id}")
        if self.G.nodes[carcinogen_id].get("type") != "Carcinogen":
            raise ValueError(f"Node {carcinogen_id!r} is not a Carcinogen node")
        return self._maximal_simple_paths(carcinogen_id, forward=True, edge_types=edge_types)

    def paths_to_carcinogen(
        self,
        carcinogen_id: str,
        *,
        edge_types: Sequence[str] | None = None,
    ) -> list[dict[str, list[Any]]]:
        """Return every maximal simple directed path ending at *carcinogen_id*.

        *carcinogen_id* must name an existing node with ``type == "Carcinogen"``
        (raises ``ValueError`` otherwise). Each path walks backward
        (predecessors) from *carcinogen_id* and is then reversed so it reads
        in natural source-to-target order, ending at *carcinogen_id* -- see
        :meth:`_maximal_simple_paths` for the exact semantics.

        This deliberately includes length-1 paths: an enzyme with a direct
        ``DETOXIFIES`` (or any other) edge straight into *carcinogen_id* --
        e.g. ``GSTM1 -> ArsenicInorganic`` -- is itself a complete, maximal
        path here, not just a building block of a longer one. Complements
        :meth:`paths_from_carcinogen`.
        """
        if carcinogen_id not in self.G:
            raise ValueError(f"Unknown node: {carcinogen_id}")
        if self.G.nodes[carcinogen_id].get("type") != "Carcinogen":
            raise ValueError(f"Node {carcinogen_id!r} is not a Carcinogen node")
        return self._maximal_simple_paths(carcinogen_id, forward=False, edge_types=edge_types)

    def carcinogen_group_paths_subgraph(
        self,
        groups: Sequence[str],
        *,
        edge_types: Sequence[str] | None = None,
    ) -> dict[str, list[Any]]:
        """Return the union of every forward path out of each carcinogen in *groups*.

        For every group in *groups* (see :meth:`get_carcinogen_groups`), every
        matching Carcinogen node (:meth:`carcinogens_by_group`) contributes
        its full set of maximal simple directed paths
        (:meth:`paths_from_carcinogen`); the result is the union of all of
        those, at the class/group level rather than per-individual-carcinogen
        -- selecting a group pulls in every carcinogen in it. A carcinogen
        with no outgoing edges still contributes its own node.

        Nodes are deduplicated by id. Edges are deduplicated by
        ``(source, target)`` -- the same structural edge recurring across
        multiple paths (or multiple carcinogens) is kept once, using the
        attributes from wherever it was first encountered. Groups matching no
        Carcinogen node contribute nothing (not an error); empty *groups*
        returns an empty subgraph -- unlike :meth:`subgraph_by_node_types`,
        this axis has no "no restriction" identity element, since selecting
        nothing here means "show no carcinogen paths", not "show every path".
        See :meth:`map_viewer_subgraph` for how the two axes combine.
        """
        node_ids: set[str] = set()
        nodes: list[dict[str, Any]] = []
        edges_by_key: dict[tuple[str, str], dict[str, Any]] = {}

        def _add_node(node_id: str) -> None:
            if node_id not in node_ids and node_id in self.G:
                node_ids.add(node_id)
                nodes.append(dict(self.G.nodes[node_id]))

        for group in groups:
            for carcinogen in self.carcinogens_by_group(group):
                carcinogen_id = carcinogen["id"]
                _add_node(carcinogen_id)
                for path in self.paths_from_carcinogen(carcinogen_id, edge_types=edge_types):
                    for node in path["nodes"]:
                        _add_node(node["id"])
                    for edge in path["edges"]:
                        key = (edge["source"], edge["target"])
                        edges_by_key.setdefault(key, edge)

        return {"nodes": nodes, "edges": list(edges_by_key.values())}

    # ── Map viewer filtering ────────────────────────────────────────────────

    def map_viewer_subgraph(
        self,
        *,
        node_types: Sequence[str] | None = None,
        carcinogen_groups: Sequence[str] | None = None,
        exclude_types: Sequence[str] = _DEFAULT_EXCLUDED_VIEWER_NODE_TYPES,
    ) -> dict[str, list[Any]]:
        """Return the Reference Map viewer's subgraph: node-type ∩ carcinogen-path axes.

        Each axis is independently optional: a ``None``/empty *node_types*
        imposes no restriction on the type axis, and a ``None``/empty
        *carcinogen_groups* imposes no restriction on the carcinogen-path
        axis (unlike :meth:`carcinogen_group_paths_subgraph` called directly,
        where empty *groups* means "match nothing" -- here it means "this
        axis doesn't apply"). Passing neither returns the full graph minus
        *exclude_types* (``Substrate`` by default, dropped unconditionally
        regardless of either axis).

        *node_types* restricts to nodes of exactly those types -- excluded
        types are absent from the result entirely, never merely
        de-emphasized. *carcinogen_groups* restricts to nodes reachable via
        :meth:`carcinogen_group_paths_subgraph`.

        The two axes are intersected on node id (this is the multi-select
        "intersection of filters" behavior, not "most recently changed filter
        wins"). Edges are the induced subgraph over the surviving node ids --
        every real graph edge between two surviving nodes, not merely the
        specific path edges that produced carcinogen-axis membership.

        This never dims -- it only removes. Tissue-threshold de-emphasis is a
        separate, non-removing overlay applied afterward; see
        :meth:`dim_by_tissue_threshold`.
        """
        type_subgraph = self.subgraph_by_node_types(node_types, exclude_types=exclude_types)
        kept_ids = {node["id"] for node in type_subgraph["nodes"]}

        if carcinogen_groups:
            carcinogen_subgraph = self.carcinogen_group_paths_subgraph(carcinogen_groups)
            kept_ids &= {node["id"] for node in carcinogen_subgraph["nodes"]}

        nodes = [dict(self.G.nodes[node_id]) for node_id in kept_ids]
        edges = [
            dict(data)
            for u, v, data in self.G.edges(data=True)
            if u in kept_ids and v in kept_ids
        ]
        return {"nodes": nodes, "edges": edges}

    def dim_by_tissue_threshold(
        self,
        subgraph: dict[str, list[Any]],
        tissue: str,
        threshold: float,
    ) -> dict[str, list[Any]]:
        """Return *subgraph* with under-expressed Enzyme nodes/edges marked ``_dimmed``.

        Every ``Enzyme`` node in *subgraph* whose ``tissue_weights[tissue]``
        is below *threshold* (or has no entry for *tissue* at all) is
        annotated ``_dimmed: True``; every edge in *subgraph* directly
        incident to a dimmed enzyme (as ``source`` or ``target`` -- not via
        the separate ``carcinogen`` context field) is dimmed too. Every other
        node and edge is annotated ``_dimmed: False``. Non-Enzyme nodes are
        never dimmed by this, regardless of *tissue*/*threshold*.

        This never removes anything: the returned dict has exactly the same
        node and edge sets as *subgraph*, just annotated -- composes cleanly
        after any node-removing filter such as :meth:`map_viewer_subgraph`.
        Non-mutating: *subgraph* and its contents are not modified; every
        returned node/edge is a fresh copy.
        """
        dimmed_ids: set[str] = set()
        nodes: list[dict[str, Any]] = []
        for node in subgraph["nodes"]:
            node = dict(node)
            is_dimmed = False
            if node.get("type") == "Enzyme":
                weights = node.get("tissue_weights") or {}
                weight = weights.get(tissue)
                is_dimmed = weight is None or weight < threshold
            node["_dimmed"] = is_dimmed
            if is_dimmed:
                dimmed_ids.add(node["id"])
            nodes.append(node)

        edges: list[dict[str, Any]] = []
        for edge in subgraph["edges"]:
            edge = dict(edge)
            edge["_dimmed"] = edge.get("source") in dimmed_ids or edge.get("target") in dimmed_ids
            edges.append(edge)

        return {"nodes": nodes, "edges": edges}

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
