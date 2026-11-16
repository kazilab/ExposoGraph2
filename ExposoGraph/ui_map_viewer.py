"""Tab: interactive Reference Map viewer.

Filtering (node type, carcinogen group, tissue threshold) is decided here in
Python via `GraphEngine.map_viewer_subgraph`/`GraphEngine.dim_by_tissue_threshold`
-- never in JavaScript. The Streamlit widgets below recompute the filtered
subgraph on every interaction and pass it (via `_map_payload`), alongside
the full viewer graph, to the map component, which lays out the full graph once
and only shows, hides, or dims its nodes and edges as the filters change. It is
drawn with the force-directed
`ExposoGraph/map/map.{html,css,js}` D3 renderer. That renderer no longer
contains any filtering logic of its own; see
`docs/design/refactor_ui_map_viewer_plan.md` for the full rationale.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import streamlit as st

from ._app_shared import (
    EDGE_COLORS,
    EDGE_LABELS,
    NODE_COLORS,
    NODE_LABELS,
    NODE_SIZES,
)
from .engine import GraphEngine
from .validation.models import NodeType

_MAP_DIR = Path(__file__).resolve().parent / "map"
_MAP_HTML = _MAP_DIR / "map.html"
_MAP_CSS = _MAP_DIR / "map.css"
_MAP_JS = _MAP_DIR / "map.js"
# D3 v7.9.2 (ISC, see map/vendor/d3.LICENSE), vendored from the npm package
# so the map doesn't depend on a CDN.
_D3_JS = _MAP_DIR / "vendor" / "d3.v7.9.2.min.js"

# Registered once at import; Streamlit mounts it into the app page (in a shadow
# root) instead of an iframe. D3 is prepended so map.js can use the global.
_MAP_COMPONENT = st.components.v2.component(
    "exposograph_reference_map",
    html=_MAP_HTML.read_text(encoding="utf-8"),
    css=_MAP_CSS.read_text(encoding="utf-8"),
    js=_D3_JS.read_text(encoding="utf-8") + "\n" + _MAP_JS.read_text(encoding="utf-8"),
)

# Substrate and CarcinogenGroup nodes are part of the underlying knowledge
# graph but must never appear in this viewer (they're intermediate bookkeeping
# and structural metadata, not something a reader needs to see on the map) --
# so they're never offered as a selectable node-type option, and
# GraphEngine excludes them by default regardless.
_SELECTABLE_NODE_TYPES = [
    t.value
    for t in NodeType
    if t not in (NodeType.SUBSTRATE, NodeType.CARCINOGEN_GROUP)
]

_NO_TISSUE_OPTION = "(none — no tissue dimming)"


def _clean_for_js(obj: Any) -> Any:
    """Strip None values so the JS side sees clean objects."""
    if isinstance(obj, dict):
        return {k: _clean_for_js(v) for k, v in obj.items() if v is not None}
    if isinstance(obj, list):
        return [_clean_for_js(i) for i in obj]
    return obj


def _tissue_note(meta: dict[str, Any]) -> str:
    """How the nodes' tissue weights were derived, from the expression table's metadata."""
    labels = meta["tissue_labels"]
    mapped = [
        f"{labels[tissue]} {'averages' if len(sources) > 1 else 'uses'} "
        + " and ".join(sources)
        for tissue, sources in meta["source_tissue_map"].items()
        if sources != [labels[tissue]]
    ]
    return " ".join(
        [
            f"{meta['source']}, {meta['unit']}; weights are each gene's value "
            "divided by its maximum across the bundled tissues.",
            "; ".join(mapped) + ".",
            f"No value in the source table for {', '.join(meta['hpa_unavailable_graph_enzymes'])}.",
        ]
    )


def _edge_key(edge: dict[str, Any]) -> str:
    """How map.js identifies an edge (it has no id of its own)."""
    return f"{edge['source']}|{edge['target']}|{edge['type']}"


def _map_payload(
    graph: dict[str, list[Any]], shown: dict[str, list[Any]], tissue_note: str
) -> dict[str, Any]:
    """Data for the map component.

    ``graph`` is the full viewer graph, which the map lays out once; ``shown``
    is the filtered subgraph of it, whose nodes/edges the map shows (dimming
    those flagged ``_dimmed``) while hiding the rest. ``tissue_note`` is shown
    with a node's tissue weights in its details.
    """
    style = {
        "node_colors": NODE_COLORS,
        "node_labels": NODE_LABELS,
        "node_sizes": NODE_SIZES,
        "edge_colors": EDGE_COLORS,
        "edge_labels": EDGE_LABELS,
    }
    return {
        "graph": _clean_for_js(graph),
        "shown": {
            "nodes": [n["id"] for n in shown["nodes"]],
            "edges": [_edge_key(e) for e in shown["edges"]],
        },
        "dimmed": {
            "nodes": [n["id"] for n in shown["nodes"] if n.get("_dimmed")],
            "edges": [_edge_key(e) for e in shown["edges"] if e.get("_dimmed")],
        },
        "style": style,
        "notes": {"tissue": tissue_note},
    }


def render(engine: GraphEngine) -> None:
    """Render the Reference Map viewer from the session engine."""
    st.markdown("#### Reference Map")
    st.caption(
        "Force-directed graph rendered from the loaded reference graph. "
        "Filters below run against the knowledge graph in Python"
    )

    if engine.node_count == 0:
        st.info("No data yet — use LLM Extract or Manual Entry to add nodes.")
        return

    carcinogen_groups = engine.get_carcinogen_groups()
    tissues = engine.get_tissues()

    col_types, col_groups = st.columns(2)
    with col_types:
        selected_node_types = st.multiselect(
            "Node types",
            options=_SELECTABLE_NODE_TYPES,
            default=[],
            help="Leave empty to show every node type (except Substrate).",
            key="map_viewer_node_types",
        )
    with col_groups:
        selected_groups = st.multiselect(
            "Carcinogen categories",
            options=carcinogen_groups,
            default=[],
            help=(
                "Restricts the map to nodes reachable via a directed path "
                "from any carcinogen in the selected categories. Leave empty "
                "for no such restriction."
            ),
            key="map_viewer_carcinogen_groups",
        )

    col_tissue, col_threshold = st.columns(2)
    with col_tissue:
        tissue_options = [_NO_TISSUE_OPTION, *tissues]
        selected_tissue = st.selectbox(
            "Tissue (dim low-expression enzymes)",
            options=tissue_options,
            index=0,
            key="map_viewer_tissue",
        )
    with col_threshold:
        threshold = st.slider(
            "Minimum nTPM weight",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.05,
            disabled=selected_tissue == _NO_TISSUE_OPTION,
            help="Enzyme nodes (and their corresponding edges) below this weight are grayed out.",
            key="map_viewer_tissue_threshold",
        )

    with st.spinner("Rendering map..."):
        subgraph = engine.map_viewer_subgraph(
            node_types=selected_node_types or None,
            carcinogen_groups=selected_groups or None,
        )
        if selected_tissue != _NO_TISSUE_OPTION:
            subgraph = engine.dim_by_tissue_threshold(
                subgraph, selected_tissue, threshold
            )

        payload = _map_payload(
            engine.map_viewer_subgraph(),
            subgraph,
            _tissue_note(engine.get_tissue_expression_metadata()),
        )

    st.caption(
        f"{len(subgraph['nodes'])} nodes / {len(subgraph['edges'])} edges shown."
    )
    _MAP_COMPONENT(key="reference_map", data=payload)
