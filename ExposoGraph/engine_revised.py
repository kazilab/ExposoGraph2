"""Streamlined GraphEngine, for comparison with ``engine.py``.

Not imported by the package. ``graph-data.json`` is the topology contract
and is read directly here. ``data/*.json`` supplies quantitative values and
provenance. This module does not parse artifacts through ``exporter`` or
ground records through ``grounding``. Load bakes the join onto the graph:
a payload with ``source_node_id`` / ``target_node_id`` rides that edge;
anything else rides the node whose id is the class id. Getters are then a
path lookup or a filtered walk.

Intentionally not a drop-in. The old per-filter getter names
(``get_node``, ``nodes_by_type``, ``subgraph_by_*``) collapse onto ``get``,
``nodes``, ``values``, and ``subgraph``.
"""

from __future__ import annotations

import json
import logging
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import networkx as nx

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

    def doc(self, name: str) -> dict[str, Any]:
        """Parsed ``data/<name>.json``, read once per engine."""
        if name not in self._docs:
            path = self.data_dir / f"{name}.json"
            self._docs[name] = json.loads(path.read_text(encoding="utf-8"))
        return self._docs[name]

    def param(self, name: str, path: str | Sequence[str], default: Any = None) -> Any:
        """Path lookup into a source document, e.g. ``param("interaction_parameters", "gsh_depletion")``."""
        return dig(self.doc(name), path, default)

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
        self._kg_provider: InteractionParameterProvider | None = None

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
        for bake in (self._bake_tissue_expression, self._bake_interaction_kinetics,
                     self._bake_flux):
            warnings += bake()
        return warnings

    # ── Bakes ────────────────────────────────────────────────────────────

    def _edges(self, source: Any, target: Any, types: set[str] | None = None) -> list[dict]:
        if not (source and target and self.G.has_edge(source, target)):
            return []
        return [d for d in self.G[source][target].values()
                if types is None or d.get("type") in types]

    def _expression(self) -> tuple[dict, dict]:
        raw = self.doc("tissue_expression_data_raw")["expression"]
        return raw, {gene: _normalize(values) for gene, values in raw.items()}

    def _bake_tissue_expression(self) -> list[str]:
        """Join the expression table on node id, not node type."""
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
        doc, warnings, seen = self.doc("interaction_parameters"), [], set()
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
        for cls, by_doc in configs.items():
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
                for cls, cfg in _entries(self.doc(name).get(root)).items():
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

    def get_flux_classes(self) -> list[str]:
        return list(self._flux_index()[0])

    def get_flux_reactions(self, cls: Any, role: str | None = None) -> list[FluxReaction]:
        reactions = self._flux_index()[0].get(_key(cls), [])
        return [r for r in reactions if role is None or r.role == role]

    def get_edge_flux_reactions(self, cls: Any) -> list[FluxReaction]:
        """Active roster read back from the graph carriers; index order if incomplete."""
        cls = _key(cls)
        carriers = [d for _, _, d in self.G.edges(data=True)] + list(self.G.nodes.values())
        walked = {
            (pathway, key): FluxReaction(carcinogen_class=cls, pathway=pathway, term_key=key,
                                         **{k: v for k, v in p.items() if k != "dispatch_status"})
            for data in carriers
            for pathway, terms in (dig(data, ("kinetics", "flux_terms", cls))
                                   or dig(data, ("flux_terms", cls)) or {}).items()
            for key, p in terms.items() if p.get("dispatch_status") == "active"
        }
        roster = self.get_flux_reactions(cls)
        if {(r.pathway, r.term_key) for r in roster} <= walked.keys():
            return [walked[k] for k in sorted(walked)]
        return [walked.get((r.pathway, r.term_key), r) for r in roster]

    def get_flux_reaction_coverage(self, cls: Any) -> dict[str, Any]:
        """Binding-edge coverage for one class. Curation report, not a pass/fail check."""
        cls = _key(cls)
        active, shadowed, configs = self._flux_index()

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

        rows = [row(r) for r in active.get(cls, [])]
        unresolved, unbound, missing_binding, missing_payload = [], [], [], []
        for reaction, item in zip(active.get(cls, []), rows):
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
            "sources": sorted(configs.get(cls, {})),
            "reaction_count": len(rows),
            "resolved_enzyme_count": len(rows) - len(unresolved),
            "edge_anchored_count": sum(item["edge_anchored"] for item in rows),
            "unresolved_terms": unresolved,
            "unbound_terms": unbound,
            "missing_binding_edges": missing_binding,
            "missing_flux_payloads": missing_payload,
            "shadowed_reactions": [row(r) for r in shadowed if r.carcinogen_class == cls],
            "reactions": rows,
        }

    def get_flux_class_config(self, cls: Any) -> dict[str, Any]:
        """Proxy-class block; empty for measured-kinetics-only classes."""
        return dict(dig(self._flux_index()[2], (_key(cls), "proxy_flux_parameters"), {}))

    def get_flux_aggregation(self, cls: Any) -> dict[str, Any]:
        return dig(self._flux_index()[2], (_key(cls), "kinetic_parameters", "aggregation"), {})

    def get_flux_metadata(self) -> dict[str, Any]:
        return self.param("kinetic_parameters", "metadata", {})

    def get_flux_provenance_entry(self, ref: str) -> dict[str, Any]:
        entry = self.param("proxy_flux_provenance", ref)
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
            value = self.param("exposure_database_revised",
                               ("exposure_scenarios", spec.get("scenario_id"), spec.get("field")))
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

    def get_parameter_provider(self) -> InteractionParameterProvider:
        baked = any(dig(d, ("kinetics", INTERACTION_BLOCK_MARKER)) == "competitive_inhibition"
                    for _, _, d in self.G.edges(data=True))
        if not baked:
            return JSONInteractionParameterProvider(self.data_dir)
        if self._kg_provider is None:
            self._kg_provider = KGInteractionParameterProvider(self, data_dir=self.data_dir)
        return self._kg_provider

    def get_induction_rules(self, exposure_context: str | None = None,
                            tissue: str | None = None) -> list[InductionRule]:
        rules = self.get_parameter_provider().get_induction_rules(tissue=tissue)
        return [r for r in rules if exposure_context in (None, r.exposure_context)]

    def get_gsh_consumers(self, tissue: str | None = None) -> list[GSHConsumer]:
        return self.get_parameter_provider().get_gsh_consumers(tissue=tissue)

    def get_genotype_modifiers(self, gene: str | None = None) -> dict[str, Any]:
        """Genotype catalog, or one gene's block."""
        from .genotypes import load_genotypes

        genes = load_genotypes()["genes"]
        return genes if gene is None else genes.get(gene, {})

    def get_genotype_modifier(self, diplotype: str, gene: str) -> float:
        """Vmax activity fraction from ``genotypes.json``."""
        from .genotypes import activity_fraction

        return activity_fraction(gene, diplotype)

    def get_proxy_genotype_modifier(self, diplotype: str, gene: str | None) -> float:
        return 1.0 if not gene else self.get_genotype_modifier(diplotype, gene)

    def get_tissue_expression(self, gene: str) -> dict[str, float] | None:
        return self._expression()[1].get(gene)

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

    def filtered_subgraph(self, *, group: str | None = None, tissue: str | None = None,
                          tissue_threshold: float | None = None, edge_type: str | None = None,
                          node_type: str | None = None) -> Subgraph:
        """Edges satisfying every given axis (per-edge AND), plus their endpoints."""
        if tissue_threshold is not None and tissue is None:
            raise ValueError("tissue_threshold requires tissue to also be given")
        axes: list[Callable[[str, str, dict], bool]] = []
        if edge_type is not None:
            axes.append(lambda u, v, d: d.get("type") == edge_type)
        if group is not None:
            ids = {d["id"] for d in self.nodes(type="Carcinogen", group=group)}
            axes.append(lambda u, v, d: bool({u, v, d.get("carcinogen")} & ids))
        if tissue is not None:
            ids_t = self._tissue_ids(tissue, tissue_threshold if tissue_threshold is not None else float("-inf"))
            axes.append(lambda u, v, d: bool({u, v} & ids_t))
        if node_type is not None:
            ids_n = {d["id"] for d in self.nodes(type=node_type)}
            axes.append(lambda u, v, d: bool({u, v} & ids_n))
        return self.subgraph(edge_where=lambda u, v, d: all(a(u, v, d) for a in axes))

    def map_viewer_subgraph(self, *, node_types: Sequence[str] | None = None,
                            carcinogen_groups: Sequence[str] | None = None,
                            exclude_types: Sequence[str] = _VIEWER_EXCLUDED) -> Subgraph:
        """Induced subgraph: node-type axis ∩ nodes reachable from the groups' carcinogens."""
        keep = {d["id"] for d in self.nodes(
            where=lambda d: d.get("type") not in exclude_types
            and (not node_types or d.get("type") in node_types))}
        if carcinogen_groups:
            starts = [d["id"] for g in carcinogen_groups for d in self.nodes(type="Carcinogen", group=g)]
            keep &= set(starts).union(*(nx.descendants(self.G, s) for s in starts))
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

    def paths(self, carcinogen_id: str, *, reverse: bool = False,
              edge_types: Sequence[str] | None = None) -> list[Subgraph]:
        """Every maximal simple directed path from (or, ``reverse``, into) a Carcinogen."""
        if self.get(carcinogen_id, key="type") != "Carcinogen":
            raise ValueError(f"Node {carcinogen_id!r} is not a Carcinogen node")
        step = self.G.in_edges if reverse else self.G.out_edges
        found: list[Subgraph] = []

        def walk(node: str, seen: list[str], edges: list[dict]) -> None:
            nxt = [(u if reverse else v, d) for u, v, d in step(node, data=True)
                   if (edge_types is None or d.get("type") in edge_types)
                   and (u if reverse else v) not in seen]
            for neighbor, data in nxt:
                walk(neighbor, seen + [neighbor], edges + [data])
            if not nxt and edges:
                order = slice(None, None, -1 if reverse else 1)
                found.append({"nodes": [dict(self.G.nodes[n]) for n in seen[order]],
                              "edges": [dict(e) for e in edges[order]]})

        walk(carcinogen_id, [carcinogen_id], [])
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
