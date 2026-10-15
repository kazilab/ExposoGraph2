"""Streamlined GraphEngine, for comparison with ``engine.py``.

Not imported by the package. ``graph-data.json`` is the topology contract
and is read directly here. ``data/*.json`` is read only while baking.
After that, getters read the graph: endpoint payloads ride edges, class
rosters and class config ride the class node, and document-level blocks
ride ``G.graph``. This module does not parse artifacts through ``exporter``
or ground records through ``grounding``.

Intentionally not a drop-in. The old per-filter getter names
(``get_node``, ``nodes_by_type``, ``subgraph_by_*``) collapse onto ``get``,
``nodes``, ``values``, and ``subgraph``.
"""

from __future__ import annotations

import json
import logging
import math
import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

import networkx as nx

from .interaction_schema import GSHConsumer, InductionRule, TissueContext
from .models import Edge, KnowledgeGraph, Node
from .parameter_provider import (
    INTERACTION_BLOCK_MARKER,
    INTERACTION_ENZYME_MARKER,
    INTERACTION_SUBSTRATE_MARKER,
)

logger = logging.getLogger(__name__)

_PKG = Path(__file__).resolve().parent
DATA_DIR = _PKG / "data"
GRAPH_PATH = _PKG / "map" / "graph-data.json"

# (document, root key, provenance prefix). Both files use pathways.<block>.
_FLUX_DOCS: tuple[tuple[str, str, str], ...] = (
    ("kinetic_parameters", "carcinogen_classes", "carcinogen_classes.{cls}.pathways"),
    ("proxy_flux_parameters", "classes", "classes.{cls}.pathways"),
)
_ROLE = {"activation": "activation", "detoxification": "detoxification", "repair": "repair"}
_TERM_META = {"graph_node_id", "rate_law", "confidence", "notes", "provenance_ref", "sources"}
_NESTING = {"carcinogen_class", "pathway", "term_key"}
_REACTION_EDGES = {"SUBSTRATE_OF", "DETOXIFIED_BY"}
_FRACTIONS = {"absorption_fraction", "tissue_blood_flow_fraction"}
_VIEWER_EXCLUDED = ("Substrate", "CarcinogenGroup")

Subgraph = dict[str, list[dict[str, Any]]]


def dig(data: Any, path: str | Sequence[str], default: Any = None) -> Any:
    """Walk a dotted (or segment-tuple) path through nested mappings."""
    for seg in path.split(".") if isinstance(path, str) else path:
        if not isinstance(data, Mapping) or seg not in data:
            return default
        data = data[seg]
    return data


def _entries(block: Any) -> dict[str, dict]:
    """Dict-valued, non-underscore children of a JSON block."""
    return {k: v for k, v in (block or {}).items()
            if not k.startswith("_") and isinstance(v, dict)}


def _nest(root: dict, *keys: str) -> dict:
    for key in keys:
        root = root.setdefault(key, {})
    return root


def _key(cls: Any) -> str:
    return getattr(cls, "value", cls)


def _clean(key: str, value: Any, fallback: Any) -> Any:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return fallback
    if not math.isfinite(value):
        return fallback
    if key in _FRACTIONS:
        return min(max(value, 1e-6), 1.0)
    return value if value > 0 else fallback


def _normalize(raw: Mapping[str, float]) -> dict[str, float]:
    top = max(raw.values(), default=0.0)
    return {t: (v / top if top else 0.0) for t, v in raw.items()}


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _activity_fraction(catalog: Mapping[str, Any], gene: str, diplotype: str) -> float:
    """Same lookup as ``genotypes.activity_fraction``, against a baked catalog."""
    key = diplotype.lower().strip()
    gene_block = catalog.get("genes", {}).get(gene.upper().strip(), {})
    for phenotype in gene_block.get("phenotypes", {}).values():
        aliases = {str(alias).lower() for alias in phenotype.get("aliases", [])}
        if key in aliases:
            return float(phenotype["activity_fraction"])
    for rule in gene_block.get("contains", []):
        if str(rule["token"]).lower() in key:
            return float(rule["activity_fraction"])
    scale_key = catalog.get("label_aliases", {}).get(key)
    if scale_key is not None:
        return float(catalog["standard_scale"][scale_key])
    try:
        value = float(diplotype)
    except (TypeError, ValueError):
        value = None
    if value is not None and 0.0 <= value <= 2.0:
        return value
    warnings.warn(
        f"Unrecognized diplotype '{diplotype}' for gene '{gene}'; defaulting to NM (1.0)",
        stacklevel=2,
    )
    return 1.0


@dataclass(frozen=True)
class FluxReaction:
    """One flux term. Its baked graph payload is ``asdict`` minus the nesting keys."""

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
    sources: list[str] | None = None
    notes: str | None = None


class GraphEngine:
    """MultiDiGraph wrapper; JSON sources are baked on load, getters read paths."""

    def __init__(self, data_dir: str | Path = DATA_DIR) -> None:
        self.G: nx.MultiDiGraph = nx.MultiDiGraph()
        self.data_dir = Path(data_dir)
        self.clear()

    # ── Sources ──────────────────────────────────────────────────────────

    def _doc(self, name: str) -> dict[str, Any]:
        """Load-time read of ``data/<name>.json``. Getters do not call this."""
        if name not in self._docs:
            path = self.data_dir / f"{name}.json"
            self._docs[name] = json.loads(path.read_text(encoding="utf-8"))
        return self._docs[name]

    # ── Mutations ────────────────────────────────────────────────────────

    def add_node(self, node: Node) -> None:
        self.G.add_node(node.id, **node.model_dump(exclude_none=True, mode="json"))

    def add_edge(self, edge: Edge) -> None:
        for role, ref in (("source", edge.source), ("target", edge.target),
                          ("carcinogen context", edge.carcinogen)):
            if ref and ref not in self.G:
                raise ValueError(f"Missing {role} node: {ref}")
        self.G.add_edge(edge.source, edge.target,
                        **edge.model_dump(exclude_none=True, mode="json"))

    def remove_node(self, node_id: str) -> None:
        if node_id in self.G:
            self.G.remove_node(node_id)

    def remove_edge(self, source: str, target: str, key: str | None = None) -> None:
        if self.G.has_edge(source, target, key):
            self.G.remove_edge(source, target, key)

    def clear(self) -> None:
        self.G.clear()
        self._docs: dict[str, dict[str, Any]] = {}
        self._flux: tuple[dict, list, dict] | None = None

    def _add_records(self, nodes: Iterable[Node | dict[str, Any]],
                     edges: Iterable[Edge | dict[str, Any]]) -> list[str]:
        """Add records as given. No grounding and no strict-mode filtering."""
        warnings: list[str] = []
        for raw in nodes:
            self.add_node(raw if isinstance(raw, Node) else Node(**raw))
        for raw in edges:
            edge = raw if isinstance(raw, Edge) else Edge(**raw)
            try:
                self.add_edge(edge)
            except ValueError as exc:
                warnings.append(str(exc))
                logger.warning("Skipped edge during merge: %s", exc)
        return warnings

    def load(self, kg: KnowledgeGraph) -> list[str]:
        self.clear()
        return self.merge(kg)

    def merge(self, kg: KnowledgeGraph) -> list[str]:
        return self._add_records(kg.nodes, kg.edges)

    def load_reference_graph(self, graph_path: str | Path = GRAPH_PATH) -> list[str]:
        """Load ``graph-data.json`` and bake every ``data/*.json`` overlay onto it."""
        path = Path(graph_path)
        if path.suffix.lower() != ".json":
            raise ValueError(f"engine_revised loads JSON only, not {path.name}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        self.clear()
        warnings = self._add_records(payload.get("nodes", []), payload.get("edges", []))
        self._bake_catalogs()
        for bake in (self._bake_tissue_expression, self._bake_interaction_kinetics,
                     self._bake_flux):
            warnings += bake()
        return warnings

    def _bake_catalogs(self) -> None:
        """Copy document-level quantitative blocks onto the graph. Load only."""
        self.G.graph["flux_metadata"] = self._doc("kinetic_parameters").get("metadata", {})
        self.G.graph["flux_provenance"] = self._doc("proxy_flux_provenance")
        self.G.graph["exposure_scenarios"] = self._doc("exposure_database_revised").get(
            "exposure_scenarios", {})
        self.G.graph["genotypes"] = self._doc("genotypes")

    # ── Bakes ────────────────────────────────────────────────────────────

    def _edges(self, source: Any, target: Any, types: set[str] | None = None) -> list[dict]:
        if not (source and target and self.G.has_edge(source, target)):
            return []
        return [d for d in self.G[source][target].values()
                if types is None or d.get("type") in types]

    def _expression(self) -> tuple[dict, dict]:
        raw = self._doc("tissue_expression_data_raw")["expression"]
        return raw, {gene: _normalize(values) for gene, values in raw.items()}

    def _bake_tissue_expression(self) -> list[str]:
        """Join the expression table on node id."""
        raw, norm = self._expression()
        warnings = []
        for node_id, data in self.G.nodes(data=True):
            if node_id in raw:
                data["tissue_weights_raw"], data["tissue_weights"] = raw[node_id], norm[node_id]
                continue
            if data.get("type") != "Enzyme":
                continue
            data.pop("tissue_weights", None)
            data.pop("tissue_weights_raw", None)
            warnings.append(f"No tissue expression data for enzyme: {node_id}")
        return warnings

    def _bake_interaction_kinetics(self) -> list[str]:
        doc, warnings, seen = self._doc("interaction_parameters"), [], set()
        self.G.graph["enzyme_induction"] = doc.get("enzyme_induction", {})
        self.G.graph["gsh_depletion"] = doc.get("gsh_depletion", {})
        for block in ("competitive_inhibition", "phase2_conjugation"):
            for enzyme, enzyme_block in _entries(doc.get(block)).items():
                target = enzyme_block.get("target_node_id")
                for substrate, params in _entries(enzyme_block.get("substrates")).items():
                    source = params.get("source_node_id")
                    edges = self._edges(source, target, _REACTION_EDGES)
                    if not edges or (source, target) in seen:
                        warnings.append(f"Unbaked {block} entry {enzyme}/{substrate}: "
                                        f"{source} -> {target}")
                        continue
                    seen.add((source, target))
                    kinetics = {k: v for k, v in params.items() if k != "source_node_id"}
                    kinetics |= {INTERACTION_ENZYME_MARKER: enzyme,
                                 INTERACTION_SUBSTRATE_MARKER: substrate,
                                 INTERACTION_BLOCK_MARKER: block}
                    for edge in edges:
                        if "flux_terms" in (edge.get("kinetics") or {}):
                            kinetics["flux_terms"] = edge["kinetics"]["flux_terms"]
                        edge["kinetics"] = dict(kinetics)
        return warnings

    def _bake_flux(self) -> list[str]:
        """Terms and signals with endpoints ride their edge; the rest ride the class node."""
        active, shadowed, configs = self._flux_index()
        warnings = []
        tagged = [(r, "active") for rs in active.values() for r in rs]
        tagged += [(r, "shadowed") for r in shadowed]
        for reaction, status in tagged:
            payload = {k: v for k, v in asdict(reaction).items() if k not in _NESTING}
            payload["dispatch_status"] = status
            if reaction.carcinogen_class in self.G:
                bucket = "flux_shadowed" if status == "shadowed" else "flux_active"
                self.G.nodes[reaction.carcinogen_class].setdefault(bucket, []).append(
                    {"pathway": reaction.pathway, "term_key": reaction.term_key, **payload})
            ends = (reaction.params.get("source_node_id"), reaction.params.get("target_node_id"))
            if all(ends):
                edges = self._edges(*ends)
                if not edges:
                    warnings.append(f"No edge for flux binding: {reaction.carcinogen_class}/"
                                    f"{reaction.term_key} ({ends[0]} -> {ends[1]})")
                    continue
                carrier = _nest(edges[0], "kinetics", "flux_terms")
            elif reaction.carcinogen_class in self.G:
                carrier = _nest(self.G.nodes[reaction.carcinogen_class], "flux_terms")
            else:
                continue
            _nest(carrier, reaction.carcinogen_class, reaction.pathway)[reaction.term_key] = payload
        self.G.graph["flux_class_order"] = list(active)
        for cls, by_doc in configs.items():
            if cls not in self.G:
                continue
            node = self.G.nodes[cls]
            proxy = by_doc.get("proxy_flux_parameters")
            if proxy:
                node["flux_class_config"] = dict(proxy)
            aggregation = (by_doc.get("kinetic_parameters") or {}).get("aggregation")
            if isinstance(aggregation, dict):
                node["flux_aggregation"] = dict(aggregation)
            node["flux_sources"] = sorted(by_doc)
            signal = dig(by_doc, "proxy_flux_parameters.signal")
            if signal:
                for edge in self._edges(signal.get("source_node_id"), signal.get("target_node_id"))[:1]:
                    _nest(edge, "kinetics", "flux_class_signal")[cls] = dict(signal)
        return warnings

    # ── Flux index (pure function of the flux JSONs) ─────────────────────

    def _enzyme_id(self, term_key: str, term: Mapping[str, Any]) -> str | None:
        target, gene = term.get("target_node_id"), term.get("gene")
        if term_key in self.G:
            return term_key
        if target in self.G and self.G.nodes[target].get("type") == "Enzyme":
            return target
        return gene if isinstance(gene, str) and gene in self.G else None

    def _reaction(self, cls: str, block: str, key: str, term: dict, cfg: dict,
                  source: str, prefix: str) -> FluxReaction:
        return FluxReaction(
            carcinogen_class=cls, term_key=key, pathway=block,
            role=_ROLE.get(block, "other"),
            rate_law=str(term.get("rate_law") or cfg.get("kinetics_model") or ""),
            params={k: v for k, v in term.items() if k not in _TERM_META},
            confidence=str(term.get("confidence") or cfg.get("confidence_overall") or ""),
            provenance_ref=str(term.get("provenance_ref") or f"{prefix}.{block}.{key}"),
            enzyme_id=self._enzyme_id(key, term),
            graph_group=str(cfg["graph_group"]) if cfg.get("graph_group") else None,
            source=source,
            sources=list(term["sources"]) if isinstance(term.get("sources"), list) else None,
            notes=str(term["notes"]) if term.get("notes") is not None else None,
        )

    def _reactions(self, cls: str, cfg: dict, name: str, prefix: str) -> list[FluxReaction]:
        return [
            self._reaction(cls, block, key, term, cfg, name, prefix.format(cls=cls))
            for block, terms in _entries(cfg.get("pathways")).items()
            for key, term in _entries(terms).items()
        ]

    def _flux_index(self) -> tuple[dict[str, list[FluxReaction]], list[FluxReaction], dict]:
        """Active roster follows ``dispatch``, not file order."""
        if self._flux is None:
            active: dict[str, list[FluxReaction]] = {}
            shadowed: list[FluxReaction] = []
            configs: dict[str, dict[str, dict]] = {}
            kinetic_shadowed: set[str] = set()
            for name, root, prefix in _FLUX_DOCS:
                for cls, cfg in _entries(self._doc(name).get(root)).items():
                    reactions = self._reactions(cls, cfg, name, prefix)
                    configs.setdefault(cls, {})[name] = cfg
                    if name == "kinetic_parameters":
                        active[cls] = reactions
                        if cfg.get("dispatch") == "shadowed":
                            kinetic_shadowed.add(cls)
                    elif cfg.get("dispatch") == "shadowed":
                        shadowed.extend(reactions)
                    elif cls in kinetic_shadowed:
                        shadowed.extend(active.get(cls, []))
                        active[cls] = reactions
                    elif cls not in active:
                        active[cls] = reactions
            self._flux = (active, shadowed, configs)
        return self._flux

    def _baked_reaction(self, cls: str, record: Mapping[str, Any]) -> FluxReaction:
        payload = {k: v for k, v in record.items() if k not in {"pathway", "term_key", "dispatch_status"}}
        return FluxReaction(carcinogen_class=cls, pathway=record["pathway"],
                            term_key=record["term_key"], **payload)

    def get_flux_classes(self) -> list[str]:
        return [cls for cls in self.G.graph.get("flux_class_order", [])
                if self.G.nodes[cls].get("flux_active")]

    def get_flux_reactions(self, cls: Any, role: str | None = None) -> list[FluxReaction]:
        cls = _key(cls)
        records = self.G.nodes[cls].get("flux_active", []) if cls in self.G else []
        reactions = [self._baked_reaction(cls, record) for record in records]
        return [r for r in reactions if role is None or r.role == role]

    def get_edge_flux_reactions(self, cls: Any) -> list[FluxReaction]:
        """Active terms on edges the constrained walk from ``cls`` crossed.

        The class node is the start. A term counts only when its binding
        edge was crossed. Terms on edges the walk does not reach, including
        a different start such as ethanol, are not returned.
        """
        cls = _key(cls)
        if cls not in self.G:
            return []
        found: list[FluxReaction] = []
        seen: set[tuple[str, str]] = set()
        for edge in self.cone(cls)["edges"]:
            blocks = ((edge.get("kinetics") or {}).get("flux_terms") or {}).get(cls) or {}
            if not isinstance(blocks, dict):
                continue
            for pathway, terms in blocks.items():
                if not isinstance(terms, dict):
                    continue
                for term_key, payload in terms.items():
                    if not isinstance(payload, dict) or payload.get("dispatch_status") == "shadowed":
                        continue
                    if (pathway, term_key) in seen:
                        continue
                    seen.add((pathway, term_key))
                    found.append(self._baked_reaction(
                        cls, {"pathway": pathway, "term_key": term_key, **payload}))
        return found

    def get_flux_reaction_coverage(self, cls: Any) -> dict[str, Any]:
        """Binding-edge coverage for one class. Curation report, not a pass/fail check."""
        cls = _key(cls)
        node = self.G.nodes[cls] if cls in self.G else {}
        active = [self._baked_reaction(cls, record) for record in node.get("flux_active", [])]
        shadowed = [self._baked_reaction(cls, record) for record in node.get("flux_shadowed", [])]

        def row(reaction: FluxReaction) -> dict[str, Any]:
            source, target = reaction.params.get("source_node_id"), reaction.params.get("target_node_id")
            edges = self._edges(source, target)
            payload = any(
                reaction.term_key in ((d.get("kinetics") or {}).get("flux_terms", {})
                                      .get(cls, {}).get(reaction.pathway, {}))
                for d in edges
            )
            return {
                "term_key": reaction.term_key,
                "pathway": reaction.pathway,
                "role": reaction.role,
                "rate_law": reaction.rate_law,
                "enzyme_id": reaction.enzyme_id,
                "source_node_id": source,
                "target_node_id": target,
                "resolved": reaction.enzyme_id in self.G if reaction.enzyme_id else False,
                "binding_edge_present": bool(edges),
                "binding_edge_types": sorted({d.get("type") for d in edges if d.get("type")}),
                "edge_anchored": payload,
            }

        rows = [row(r) for r in active]
        unresolved, unbound, missing_binding, missing_payload = [], [], [], []
        for reaction, item in zip(active, rows):
            if not item["resolved"]:
                unresolved.append(reaction.term_key)
                continue
            if not (item["source_node_id"] and item["target_node_id"]):
                unbound.append(reaction.term_key)
            elif not item["binding_edge_present"]:
                missing_binding.append(reaction.term_key)
            elif not item["edge_anchored"]:
                missing_payload.append(reaction.term_key)
        return {
            "carcinogen_class": cls,
            "sources": list(node.get("flux_sources", [])),
            "reaction_count": len(rows),
            "resolved_enzyme_count": len(rows) - len(unresolved),
            "edge_anchored_count": sum(item["edge_anchored"] for item in rows),
            "unresolved_terms": unresolved,
            "unbound_terms": unbound,
            "missing_binding_edges": missing_binding,
            "missing_flux_payloads": missing_payload,
            "shadowed_reactions": [row(r) for r in shadowed],
            "reactions": rows,
        }

    def get_flux_class_config(self, cls: Any) -> dict[str, Any]:
        """Proxy-class block baked on the class node."""
        cls = _key(cls)
        return dict(self.G.nodes[cls].get("flux_class_config", {})) if cls in self.G else {}

    def get_flux_aggregation(self, cls: Any) -> dict[str, Any]:
        cls = _key(cls)
        return dict(self.G.nodes[cls].get("flux_aggregation", {})) if cls in self.G else {}

    def get_flux_metadata(self) -> dict[str, Any]:
        return dict(self.G.graph.get("flux_metadata") or {})

    def get_flux_provenance_entry(self, ref: str) -> dict[str, Any]:
        entry = dig(self.G.graph.get("flux_provenance") or {}, ref)
        if entry is None:
            raise KeyError(ref)
        return entry

    def get_default_concentration(self, cls: Any) -> float:
        defaults = self.get_flux_metadata().get("exposure_defaults_uM", {})
        spec = self.get_flux_class_config(cls).get("exposure_default") or {}
        fields = [self.get_flux_aggregation(cls).get("exposure_default_field")]
        if spec.get("source") == "kinetic_parameters":
            fields.append(spec.get("field"))
        for value in [defaults.get(f) for f in fields if f]:
            if value is not None:
                return float(value)
        if spec.get("source") == "exposure_database":
            value = dig(self.G.graph.get("exposure_scenarios") or {},
                        (spec.get("scenario_id"), spec.get("field")))
            if value is not None:
                return float(value)
        return 0.1

    # ── Physiology contexts ──────────────────────────────────────────────

    def get_qivive_context(self, tissue: str, overrides: Mapping[str, float] | None = None) -> dict[str, float]:
        tissues = dig(self.get_flux_metadata(), "qivive_defaults.tissues")
        merged = {**tissues.get(tissue, tissues["liver"]), **(overrides or {})}
        ctx = {k: float(merged[k]) for k in ("mppgl_mg_per_g", "organ_weight_g")}
        ctx["scale"] = round(ctx["mppgl_mg_per_g"] * ctx["organ_weight_g"], 6)
        return ctx

    def get_steady_state_context(self, tissue: str, overrides: Mapping[str, float] | None = None) -> dict[str, float]:
        cfg = self.get_flux_metadata()["steady_state_defaults"]
        base = {k: v for k, v in cfg.items() if k not in ("description", "tissues")}
        base |= cfg["tissues"].get(tissue, cfg["tissues"]["liver"])
        ctx = {k: _clean(k, v, None) for k, v in base.items()}
        if None in ctx.values():
            raise ValueError(f"Invalid steady_state_defaults for {tissue}: {ctx}")
        for k, v in (overrides or {}).items():
            ctx[k] = _clean(k, v, ctx.get(k, 1.0))
        ctx["central_volume_l"] = round(ctx["body_weight_kg"] * ctx["volume_l_per_kg"], 6)
        ctx["tissue_volume_l"] = round(ctx["organ_weight_g"] / 1000.0, 6)
        ctx["tissue_blood_flow_l_per_day"] = round(
            ctx["cardiac_output_l_per_day"] * ctx["tissue_blood_flow_fraction"], 6)
        return ctx

    # ── Interaction parameters ───────────────────────────────────────────

    def get_induction_rules(self, exposure_context: str | None = None,
                            tissue: str | None = None) -> list[InductionRule]:
        rules = []
        for context, exposure_data in (self.G.graph.get("enzyme_induction") or {}).items():
            if not isinstance(exposure_data, dict):
                continue
            for enzyme, rule_data in exposure_data.items():
                if enzyme.startswith("_") or not isinstance(rule_data, dict):
                    continue
                if exposure_context not in (None, context):
                    continue
                specificity = rule_data.get("tissue_specificity")
                if tissue and specificity and tissue.lower() not in str(specificity).lower():
                    continue
                rules.append(InductionRule(
                    exposure_context=context,
                    enzyme=enzyme,
                    fold_induction=_optional_float(rule_data.get("fold_induction")),
                    range_min=_optional_float(rule_data.get("range_min")),
                    range_max=_optional_float(rule_data.get("range_max")),
                    mechanism=rule_data.get("mechanism"),
                    tissue_context=TissueContext(tissue=specificity),
                    metadata=dict(rule_data),
                ))
        return rules

    def get_gsh_consumers(self, tissue: str | None = None) -> list[GSHConsumer]:
        consumers = (self.G.graph.get("gsh_depletion") or {}).get("consumers", {})
        parsed = []
        for name, data in consumers.items():
            data = data if isinstance(data, dict) else {}
            context = TissueContext(tissue=data.get("tissue"))
            if tissue and context.tissue and context.tissue.lower() != tissue.lower():
                continue
            parsed.append(GSHConsumer(
                name=name,
                enzyme=data.get("enzyme"),
                substrate_class=data.get("substrate_class"),
                gsh_per_umol_substrate=_optional_float(data.get("gsh_per_umol_substrate")),
                tissue_context=context,
                metadata=dict(data),
            ))
        return parsed

    def get_genotype_modifiers(self, gene: str | None = None) -> dict[str, Any]:
        """Baked genotype catalog, or one gene's block."""
        genes = (self.G.graph.get("genotypes") or {}).get("genes", {})
        return genes if gene is None else genes.get(gene, {})

    def get_genotype_modifier(self, diplotype: str, gene: str) -> float:
        """Vmax activity fraction from the catalog baked on the graph."""
        return _activity_fraction(self.G.graph.get("genotypes") or {}, gene, diplotype)

    def get_proxy_genotype_modifier(self, diplotype: str, gene: str | None) -> float:
        return 1.0 if not gene else self.get_genotype_modifier(diplotype, gene)

    def get_tissue_expression(self, gene: str) -> dict[str, float] | None:
        return self.get(gene, key="tissue_weights")

    # ── Queries ──────────────────────────────────────────────────────────

    @property
    def node_count(self) -> int:
        return self.G.number_of_nodes()

    @property
    def edge_count(self) -> int:
        return self.G.number_of_edges()

    def get(self, source: str, target: str | None = None, key: str | Sequence[str] | None = None,
            *, edge_key: Any = None, default: Any = None) -> Any:
        """Node (``target=None``) or edge attributes, optionally drilled to ``key``."""
        if target is None:
            data = self.G.nodes[source] if source in self.G else None
        else:
            view = self.G[source][target] if self.G.has_edge(source, target) else {}
            data = view.get(edge_key) if edge_key is not None else next(iter(view.values()), None)
        if data is None:
            return None if key is None else default
        return dict(data) if key is None else dig(data, key, default)

    def nodes(self, where: Callable[[dict], bool] | None = None, **attrs: Any) -> list[dict[str, Any]]:
        """Node dicts matching every ``attr=value`` and the optional predicate."""
        return [dict(d) for _, d in self.G.nodes(data=True)
                if all(d.get(k) == v for k, v in attrs.items()) and (where is None or where(d))]

    def values(self, attr: str, **attrs: Any) -> list[Any]:
        """Sorted distinct values (or mapping keys) of ``attr`` over matching nodes."""
        found: set[Any] = set()
        for d in self.nodes(**attrs):
            value = d.get(attr)
            found.update(value if isinstance(value, Mapping) else [value] if value else [])
        return sorted(found)

    def neighbors(self, node_id: str) -> list[str]:
        return list(nx.all_neighbors(self.G, node_id)) if node_id in self.G else []

    def subgraph(self, node_ids: Iterable[str] | None = None,
                 edge_where: Callable[[str, str, dict], bool] | None = None) -> Subgraph:
        """Induced subgraph on ``node_ids``, or matching edges plus their endpoints."""
        if edge_where is None:
            keep = set(self.G) if node_ids is None else set(node_ids) & set(self.G)
            edges = [dict(d) for u, v, d in self.G.edges(data=True) if u in keep and v in keep]
        else:
            matched = [(u, v, d) for u, v, d in self.G.edges(data=True) if edge_where(u, v, d)]
            keep = {n for u, v, _ in matched for n in (u, v)}
            edges = [dict(d) for _, _, d in matched]
        return {"nodes": [dict(self.G.nodes[n]) for n in keep], "edges": edges}

    def _tissue_ids(self, tissue: str, threshold: float = float("-inf")) -> set[str]:
        return {d["id"] for d in self.nodes(type="Enzyme")
                if (d.get("tissue_weights") or {}).get(tissue, float("-inf")) >= threshold}

    def _steps(self, node: str, arrival: str | None, *,
               reverse: bool = False, edge_types: Sequence[str] | None = None):
        """Outgoing steps that keep the substrate or agonist we arrived with.

        An enzyme ``PRODUCES`` edge is taken only when its ``substrate_node_id``
        is the node we entered from. A receptor ``INDUCES`` or ``INHIBITS``
        edge is taken only when its ``agonist_node_id`` is that same node.
        Reverse walks apply the same pairing from the downstream side.
        """
        ntype = self.G.nodes[node].get("type") if node in self.G else None
        step = self.G.in_edges if reverse else self.G.out_edges
        for u, v, d in step(node, data=True):
            other = u if reverse else v
            if edge_types is not None and d.get("type") not in edge_types:
                continue
            et = d.get("type")
            if not reverse:
                if et == "PRODUCES" and ntype == "Enzyme" and d.get("substrate_node_id") != arrival:
                    continue
                if et in ("INDUCES", "INHIBITS") and ntype == "Receptor" and d.get("agonist_node_id") != arrival:
                    continue
            elif ntype == "Enzyme" and et == "SUBSTRATE_OF":
                if not any(p.get("type") == "PRODUCES" and p.get("target") == arrival
                           and p.get("substrate_node_id") == other
                           for _, _, p in self.G.out_edges(node, data=True)):
                    continue
            elif ntype == "Receptor" and et == "AGONIZES":
                if not any(p.get("type") in ("INDUCES", "INHIBITS") and p.get("target") == arrival
                           and p.get("agonist_node_id") == other
                           for _, _, p in self.G.out_edges(node, data=True)):
                    continue
            yield other, d

    def cone(self, *starts: str) -> Subgraph:
        """Nodes and edges reached by a constrained walk from ``starts``.

        The same enzyme may be entered from more than one substrate. Each
        entry is expanded once. Edges are the ones the walk crossed, not
        every edge among the reached nodes.
        """
        queued = {(s, None) for s in starts if s in self.G}
        seen: set[tuple[str, str | None]] = set()
        edges: list[dict] = []
        taken: set[tuple[str, str, int]] = set()
        while queued:
            node, arrival = queued.pop()
            if (node, arrival) in seen:
                continue
            seen.add((node, arrival))
            for nxt, data in self._steps(node, arrival):
                key = (data.get("source", node), data.get("target", nxt), id(data))
                if key not in taken:
                    taken.add(key)
                    edges.append(dict(data))
                queued.add((nxt, node))
        node_ids = {n for n, _ in seen}
        return {"nodes": [dict(self.G.nodes[n]) for n in node_ids], "edges": edges}

    def filtered_subgraph(self, *, group: str | None = None, tissue: str | None = None,
                          tissue_threshold: float | None = None, edge_type: str | None = None,
                          node_type: str | None = None) -> Subgraph:
        """Constrained cone for a carcinogen group, then the other axes.

        A group is the union of walks from its carcinogens, not every edge
        tagged with those ids. Tissue, edge type, and node type still filter
        that cone. With no group, those axes apply to the whole graph.
        """
        if tissue_threshold is not None and tissue is None:
            raise ValueError("tissue_threshold requires tissue to also be given")
        if group is not None:
            starts = [d["id"] for d in self.nodes(type="Carcinogen", group=group)]
            walked = self.cone(*starts)
            edges = walked["edges"]
        else:
            edges = [dict(d) for _, _, d in self.G.edges(data=True)]
        if edge_type is not None:
            edges = [e for e in edges if e.get("type") == edge_type]
        if tissue is not None:
            ids_t = self._tissue_ids(tissue, tissue_threshold if tissue_threshold is not None else float("-inf"))
            edges = [e for e in edges if {e.get("source"), e.get("target")} & ids_t]
        if node_type is not None:
            ids_n = {d["id"] for d in self.nodes(type=node_type)}
            edges = [e for e in edges if {e.get("source"), e.get("target")} & ids_n]
        keep = {n for e in edges for n in (e.get("source"), e.get("target")) if n in self.G}
        if group is not None:
            keep |= {d["id"] for d in self.nodes(type="Carcinogen", group=group)}
        return {"nodes": [dict(self.G.nodes[n]) for n in keep], "edges": edges}

    def map_viewer_subgraph(self, *, node_types: Sequence[str] | None = None,
                            carcinogen_groups: Sequence[str] | None = None,
                            exclude_types: Sequence[str] = _VIEWER_EXCLUDED) -> Subgraph:
        """Node-type filter over the constrained cone, not the induced descendant set."""
        if carcinogen_groups:
            starts = [d["id"] for g in carcinogen_groups for d in self.nodes(type="Carcinogen", group=g)]
            walked = self.cone(*starts)
            def allowed(node_id: str) -> bool:
                kind = self.G.nodes[node_id].get("type")
                return kind not in exclude_types and (not node_types or kind in node_types)

            edges = [e for e in walked["edges"] if allowed(e["source"]) and allowed(e["target"])]
            keep = {n for e in edges for n in (e["source"], e["target"])}
            keep |= {s for s in starts if allowed(s)}
            return {"nodes": [dict(self.G.nodes[n]) for n in keep], "edges": edges}
        keep = {d["id"] for d in self.nodes(
            where=lambda d: d.get("type") not in exclude_types
            and (not node_types or d.get("type") in node_types))}
        return self.subgraph(keep)

    def dim_by_tissue_threshold(self, subgraph: Subgraph, tissue: str, threshold: float) -> Subgraph:
        def dimmed(n: dict) -> bool:
            weight = (n.get("tissue_weights") or {}).get(tissue)
            return n.get("type") == "Enzyme" and (weight is None or weight < threshold)

        nodes = [{**n, "_dimmed": dimmed(n)} for n in subgraph["nodes"]]
        off = {n["id"] for n in nodes if n["_dimmed"]}
        edges = [{**e, "_dimmed": bool({e.get("source"), e.get("target")} & off)}
                 for e in subgraph["edges"]]
        return {"nodes": nodes, "edges": edges}

    def paths(self, source_node_id: str, *, reverse: bool = False,
              edge_types: Sequence[str] | None = None) -> list[Subgraph]:
        """Every maximal simple path from a node, without fanning through a shared enzyme or receptor."""
        found: list[Subgraph] = []

        def walk(node: str, arrival: str | None, seen: list[str], edges: list[dict]) -> None:
            nxt = [(nbr, data) for nbr, data in self._steps(
                node, arrival, reverse=reverse, edge_types=edge_types) if nbr not in seen]
            for neighbor, data in nxt:
                walk(neighbor, node, seen + [neighbor], edges + [data])
            if not nxt and edges:
                order = slice(None, None, -1 if reverse else 1)
                found.append({"nodes": [dict(self.G.nodes[n]) for n in seen[order]],
                              "edges": [dict(e) for e in edges[order]]})

        walk(source_node_id, None, [source_node_id], [])
        return found

    # ── Serialization ────────────────────────────────────────────────────

    def to_dict(self) -> Subgraph:
        return self.subgraph()

    def to_knowledge_graph(self) -> KnowledgeGraph:
        data = self.to_dict()
        return KnowledgeGraph(nodes=[Node(**n) for n in data["nodes"]],
                              edges=[Edge(**e) for e in data["edges"]])

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, default=str)

    def validate(self) -> list[str]:
        return [f"Edge {u}→{v} references missing node {ref!r}"
                for u, v, d in self.G.edges(data=True)
                for ref in (u, v, d.get("carcinogen")) if ref and ref not in self.G]
