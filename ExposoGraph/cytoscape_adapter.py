"""Cytoscape-oriented graph export helpers for the advanced viewer."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Mapping

import networkx as nx

from .config import GraphVisibility, normalize_graph_visibility
from .engine import GraphEngine
from .graph_filters import filter_knowledge_graph, filtered_engine
from .models import KnowledgeGraph


class ViewerLayoutMode(str, Enum):
    COSE = "cose"
    BREADTHFIRST = "breadthfirst"
    CIRCLE = "circle"
    PRESET = "preset"


_NODE_COLORS: dict[str, str] = {
    "Carcinogen": "#e05565",
    "Enzyme": "#4f98a3",
    "Gene": "#3d8b8b",
    "Metabolite": "#e8945a",
    "DNA_Adduct": "#a86fdf",
    "Pathway": "#5591c7",
    "Tissue": "#c2855a",
}

_NODE_SHAPES: dict[str, str] = {
    "Carcinogen": "diamond",
    "Enzyme": "ellipse",
    "Gene": "ellipse",
    "Metabolite": "ellipse",
    "DNA_Adduct": "hexagon",
    "Pathway": "round-rectangle",
    "Tissue": "tag",
}

_NODE_SIZES: dict[str, int] = {
    "Carcinogen": 52,
    "Enzyme": 40,
    "Gene": 36,
    "Metabolite": 34,
    "DNA_Adduct": 38,
    "Pathway": 58,
    "Tissue": 34,
}

_EDGE_COLORS: dict[str, str] = {
    "ACTIVATES": "#e05565",
    "DETOXIFIES": "#6daa45",
    "TRANSPORTS": "#5591c7",
    "FORMS_ADDUCT": "#a86fdf",
    "REPAIRS": "#e8af34",
    "PATHWAY": "#707a8a",
    "EXPRESSED_IN": "#c2855a",
    "INDUCES": "#d4a843",
    "INHIBITS": "#8b4a6b",
    "ENCODES": "#3d8b8b",
    "CUSTOM": "#9ea9bd",
}

_EDGE_WIDTHS: dict[str, int] = {
    "ACTIVATES": 3,
    "DETOXIFIES": 3,
    "TRANSPORTS": 3,
    "FORMS_ADDUCT": 4,
    "REPAIRS": 3,
    "PATHWAY": 2,
    "EXPRESSED_IN": 2,
    "INDUCES": 3,
    "INHIBITS": 3,
    "ENCODES": 2,
    "CUSTOM": 2,
}


@dataclass(frozen=True)
class CytoscapeBundle:
    elements: list[dict[str, Any]]
    stylesheet: list[dict[str, Any]]
    layout: dict[str, Any]
    metadata: dict[str, Any]
    positions: dict[str, dict[str, float]]

    def to_dict(self) -> dict[str, Any]:
        return {
            "elements": self.elements,
            "stylesheet": self.stylesheet,
            "layout": self.layout,
            "metadata": self.metadata,
            "positions": self.positions,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CytoscapeBundle":
        return cls(
            elements=list(data.get("elements", [])),
            stylesheet=list(data.get("stylesheet", [])),
            layout=dict(data.get("layout", {})),
            metadata=dict(data.get("metadata", {})),
            positions={
                str(node_id): {
                    "x": float(coords["x"]),
                    "y": float(coords["y"]),
                }
                for node_id, coords in dict(data.get("positions", {})).items()
                if isinstance(coords, Mapping) and "x" in coords and "y" in coords
            },
        )


def normalize_viewer_layout_mode(value: str | ViewerLayoutMode | None) -> ViewerLayoutMode:
    if isinstance(value, ViewerLayoutMode):
        return value
    if value is None:
        return ViewerLayoutMode.COSE

    cleaned = value.strip().lower()
    aliases = {
        "cose": ViewerLayoutMode.COSE,
        "force": ViewerLayoutMode.COSE,
        "force-directed": ViewerLayoutMode.COSE,
        "breadthfirst": ViewerLayoutMode.BREADTHFIRST,
        "hierarchical": ViewerLayoutMode.BREADTHFIRST,
        "circle": ViewerLayoutMode.CIRCLE,
        "preset": ViewerLayoutMode.PRESET,
        "saved": ViewerLayoutMode.PRESET,
        "fixed": ViewerLayoutMode.PRESET,
    }
    return aliases.get(cleaned, ViewerLayoutMode.COSE)


def _slugify_class(value: str | None) -> str:
    if not value:
        return "unknown"
    chars = [
        ch.lower() if ch.isalnum() else "-"
        for ch in value.strip()
    ]
    collapsed = "".join(chars).strip("-")
    while "--" in collapsed:
        collapsed = collapsed.replace("--", "-")
    return collapsed or "unknown"


def _coerce_engine(
    graph: GraphEngine | KnowledgeGraph,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> GraphEngine:
    normalized = (
        visibility
        if isinstance(visibility, GraphVisibility)
        else normalize_graph_visibility(visibility)
    )
    if isinstance(graph, GraphEngine):
        return filtered_engine(graph, normalized)

    filtered_graph = filter_knowledge_graph(graph, normalized)
    engine = GraphEngine()
    for node in filtered_graph.nodes:
        engine.add_node(node)
    for edge in filtered_graph.edges:
        engine.add_edge(edge)
    return engine


def compute_viewer_positions(
    graph: GraphEngine | KnowledgeGraph,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    seed: int = 42,
    width: int = 1400,
    height: int = 980,
) -> dict[str, dict[str, float]]:
    engine = _coerce_engine(graph, visibility=visibility)
    if engine.node_count == 0:
        return {}
    if engine.node_count == 1:
        only = next(iter(engine.G.nodes()))
        return {str(only): {"x": width / 2, "y": height / 2}}

    layout_graph = nx.Graph()
    layout_graph.add_nodes_from(engine.G.nodes())
    layout_graph.add_edges_from((u, v) for u, v in engine.G.edges())
    spring_k = 1.8 / math.sqrt(max(layout_graph.number_of_nodes(), 1))
    raw_positions = nx.spring_layout(layout_graph, seed=seed, k=spring_k, iterations=250)

    scale_x = width * 0.36
    scale_y = height * 0.34
    center_x = width / 2
    center_y = height / 2
    return {
        str(node_id): {
            "x": round(center_x + float(coords[0]) * scale_x, 2),
            "y": round(center_y + float(coords[1]) * scale_y, 2),
        }
        for node_id, coords in raw_positions.items()
    }


def load_viewer_positions(path: str | Path) -> dict[str, dict[str, float]]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        str(node_id): {
            "x": float(coords["x"]),
            "y": float(coords["y"]),
        }
        for node_id, coords in raw.items()
        if isinstance(coords, Mapping) and "x" in coords and "y" in coords
    }


def load_cytoscape_bundle(path: str | Path) -> CytoscapeBundle:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return CytoscapeBundle.from_dict(raw)


def write_viewer_positions(
    graph: GraphEngine | KnowledgeGraph,
    path: str | Path,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    seed: int = 42,
    positions: Mapping[str, Mapping[str, float]] | None = None,
) -> Path:
    path = Path(path)
    payload = (
        {
            str(node_id): {"x": float(coords["x"]), "y": float(coords["y"])}
            for node_id, coords in positions.items()
        }
        if positions is not None
        else compute_viewer_positions(graph, visibility=visibility, seed=seed)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def viewer_stylesheet() -> list[dict[str, Any]]:
    return [
        {
            "selector": "node",
            "style": {
                "label": "data(label)",
                "background-color": "data(color)",
                "shape": "data(shape)",
                "width": "data(size)",
                "height": "data(size)",
                "border-width": 2,
                "border-color": "#08111f",
                "color": "#e6edf7",
                "font-size": 11,
                "font-family": "Inter, system-ui, sans-serif",
                "text-wrap": "wrap",
                "text-max-width": 120,
                "text-valign": "bottom",
                "text-margin-y": 12,
                "text-outline-color": "#08111f",
                "text-outline-width": 2,
                "background-opacity": 0.94,
                "overlay-opacity": 0,
                "z-index": 10,
            },
        },
        {
            "selector": "edge",
            "style": {
                "curve-style": "bezier",
                "line-color": "data(color)",
                "target-arrow-color": "data(color)",
                "target-arrow-shape": "triangle",
                "arrow-scale": 1.1,
                "width": "data(width)",
                "opacity": 0.72,
                "line-style": "solid",
                "z-index": 2,
            },
        },
        {
            "selector": ".type-carcinogen",
            "style": {
                "font-size": 12,
                "font-weight": 700,
                "text-max-width": 132,
            },
        },
        {
            "selector": ".type-pathway",
            "style": {
                "font-size": 10,
                "font-weight": 600,
                "background-opacity": 0.72,
                "text-max-width": 140,
            },
        },
        {
            "selector": ".edge-pathway",
            "style": {
                "opacity": 0.28,
                "width": 2,
                "line-style": "dashed",
                "target-arrow-shape": "none",
            },
        },
        {
            "selector": ".edge-forms-adduct",
            "style": {
                "width": 4,
                "opacity": 0.88,
            },
        },
        {
            "selector": ".edge-repairs",
            "style": {
                "line-style": "dotted",
            },
        },
        {
            "selector": ".edge-transports",
            "style": {
                "curve-style": "unbundled-bezier",
            },
        },
        {
            "selector": ".match-canonical, .match-alias",
            "style": {
                "border-style": "solid",
                "border-width": 3,
            },
        },
        {
            "selector": ".match-unmatched",
            "style": {
                "border-style": "dashed",
                "border-color": "#f5c06b",
            },
        },
        {
            "selector": ".match-custom",
            "style": {
                "border-style": "double",
                "border-color": "#d9dee7",
            },
        },
        {
            "selector": ".origin-llm",
            "style": {
                "border-color": "#6ec8ff",
            },
        },
        {
            "selector": ".origin-user",
            "style": {
                "border-color": "#f3e9a7",
            },
        },
        {
            "selector": ".origin-seeded",
            "style": {
                "border-color": "#7de3a0",
            },
        },
        {
            "selector": ".dimmed",
            "style": {
                "opacity": 0.08,
                "text-opacity": 0.12,
            },
        },
        {
            "selector": ".connected",
            "style": {
                "opacity": 1,
                "text-opacity": 1,
                "z-index": 25,
            },
        },
        {
            "selector": ".selected",
            "style": {
                "border-width": 4,
                "border-color": "#ffffff",
                "shadow-blur": 18,
                "shadow-color": "#ffffff",
                "shadow-opacity": 0.35,
                "z-index": 30,
            },
        },
        {
            "selector": "edge.selected",
            "style": {
                "width": 5,
                "opacity": 0.96,
                "target-arrow-color": "#ffffff",
                "line-color": "#ffffff",
                "z-index": 32,
            },
        },
    ]


def viewer_layout(
    layout_mode: ViewerLayoutMode | str = ViewerLayoutMode.COSE,
) -> dict[str, Any]:
    normalized = normalize_viewer_layout_mode(layout_mode)
    if normalized == ViewerLayoutMode.PRESET:
        return {"name": "preset", "fit": True, "padding": 40, "animate": False}
    if normalized == ViewerLayoutMode.BREADTHFIRST:
        return {
            "name": "breadthfirst",
            "fit": True,
            "padding": 40,
            "animate": True,
            "spacingFactor": 1.15,
        }
    if normalized == ViewerLayoutMode.CIRCLE:
        return {"name": "circle", "fit": True, "padding": 50, "animate": True}
    return {
        "name": "cose",
        "fit": True,
        "padding": 45,
        "animate": True,
        "nodeRepulsion": 120000,
        "idealEdgeLength": 120,
        "edgeElasticity": 140,
    }


def _node_classes(data: Mapping[str, Any]) -> str:
    classes = [
        "node",
        f"type-{_slugify_class(str(data.get('type', 'Node')))}",
        f"origin-{_slugify_class(str(data.get('origin', 'imported')))}",
        f"match-{_slugify_class(str(data.get('match_status', 'unknown')))}",
    ]
    group = data.get("group")
    if group:
        classes.append(f"group-{_slugify_class(str(group))}")
    return " ".join(classes)


def _edge_classes(data: Mapping[str, Any]) -> str:
    classes = [
        "edge",
        f"edge-{_slugify_class(str(data.get('type', 'EDGE')))}",
        f"origin-{_slugify_class(str(data.get('origin', 'imported')))}",
        f"match-{_slugify_class(str(data.get('match_status', 'unknown')))}",
    ]
    return " ".join(classes)


def _node_data(data: Mapping[str, Any]) -> dict[str, Any]:
    node_type = str(data.get("type", "Node"))
    return {
        **dict(data),
        "kind": "node",
        "color": _NODE_COLORS.get(node_type, "#76c3ff"),
        "shape": _NODE_SHAPES.get(node_type, "ellipse"),
        "size": _NODE_SIZES.get(node_type, 36),
    }


def _edge_data(data: Mapping[str, Any], *, edge_id: str) -> dict[str, Any]:
    edge_type = str(data.get("type", "EDGE"))
    label = data.get("label") or edge_type.replace("_", " ").title()
    return {
        **dict(data),
        "id": edge_id,
        "kind": "edge",
        "label": label,
        "color": _EDGE_COLORS.get(edge_type, "#8ea4bb"),
        "width": _EDGE_WIDTHS.get(edge_type, 2),
    }


def build_cytoscape_elements(
    graph: GraphEngine | KnowledgeGraph,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    positions: Mapping[str, Mapping[str, float]] | None = None,
) -> list[dict[str, Any]]:
    engine = _coerce_engine(graph, visibility=visibility)
    elements: list[dict[str, Any]] = []
    edge_counts: dict[str, int] = {}

    for node_id, data in engine.G.nodes(data=True):
        element: dict[str, Any] = {
            "data": _node_data(data),
            "classes": _node_classes(data),
            "selectable": True,
        }
        coords = positions.get(str(node_id)) if positions is not None else None
        if coords is not None and "x" in coords and "y" in coords:
            element["position"] = {
                "x": float(coords["x"]),
                "y": float(coords["y"]),
            }
        elements.append(element)

    for source, target, key, data in engine.G.edges(keys=True, data=True):
        edge_key = str(key)
        count = edge_counts.get(edge_key, 0)
        edge_counts[edge_key] = count + 1
        edge_id = edge_key if count == 0 else f"{edge_key}:{count}"
        element = {
            "data": _edge_data(data, edge_id=edge_id),
            "classes": _edge_classes(data),
            "selectable": True,
        }
        elements.append(element)

    return elements


def build_cytoscape_metadata(
    graph: GraphEngine | KnowledgeGraph,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    layout_mode: ViewerLayoutMode | str = ViewerLayoutMode.COSE,
) -> dict[str, Any]:
    engine = _coerce_engine(graph, visibility=visibility)
    node_type_counts: dict[str, int] = {}
    edge_type_counts: dict[str, int] = {}
    carcinogen_groups: dict[str, int] = {}

    for _, data in engine.G.nodes(data=True):
        node_type = str(data.get("type", "Node"))
        node_type_counts[node_type] = node_type_counts.get(node_type, 0) + 1
        if node_type == "Carcinogen" and data.get("group"):
            group = str(data["group"])
            carcinogen_groups[group] = carcinogen_groups.get(group, 0) + 1

    for _, _, data in engine.G.edges(data=True):
        edge_type = str(data.get("type", "EDGE"))
        edge_type_counts[edge_type] = edge_type_counts.get(edge_type, 0) + 1

    normalized_visibility = (
        visibility.value
        if isinstance(visibility, GraphVisibility)
        else normalize_graph_visibility(visibility).value
    )
    normalized_layout = normalize_viewer_layout_mode(layout_mode).value
    return {
        "visibility": normalized_visibility,
        "layout_mode": normalized_layout,
        "node_count": engine.node_count,
        "edge_count": engine.edge_count,
        "node_type_counts": dict(sorted(node_type_counts.items())),
        "edge_type_counts": dict(sorted(edge_type_counts.items())),
        "carcinogen_groups": dict(sorted(carcinogen_groups.items())),
    }


def build_cytoscape_bundle(
    graph: GraphEngine | KnowledgeGraph,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    positions: Mapping[str, Mapping[str, float]] | None = None,
    layout_mode: ViewerLayoutMode | str = ViewerLayoutMode.COSE,
) -> CytoscapeBundle:
    normalized_layout = normalize_viewer_layout_mode(layout_mode)
    resolved_positions = (
        {
            str(node_id): {
                "x": float(coords["x"]),
                "y": float(coords["y"]),
            }
            for node_id, coords in positions.items()
        }
        if positions is not None
        else compute_viewer_positions(graph, visibility=visibility)
    )
    return CytoscapeBundle(
        elements=build_cytoscape_elements(
            graph,
            visibility=visibility,
            positions=resolved_positions or None,
        ),
        stylesheet=viewer_stylesheet(),
        layout=viewer_layout(normalized_layout),
        metadata=build_cytoscape_metadata(
            graph,
            visibility=visibility,
            layout_mode=normalized_layout,
        ),
        positions=resolved_positions,
    )


def write_cytoscape_bundle(
    graph: GraphEngine | KnowledgeGraph,
    path: str | Path,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    positions: Mapping[str, Mapping[str, float]] | None = None,
    layout_mode: ViewerLayoutMode | str = ViewerLayoutMode.COSE,
) -> Path:
    path = Path(path)
    bundle = build_cytoscape_bundle(
        graph,
        visibility=visibility,
        positions=positions,
        layout_mode=layout_mode,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bundle.to_dict(), indent=2), encoding="utf-8")
    return path
