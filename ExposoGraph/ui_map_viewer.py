"""Tab: interactive Reference Map viewer.

Filtering (node type, carcinogen group, tissue threshold) is decided here in
Python via `GraphEngine.map_viewer_subgraph`/`GraphEngine.dim_by_tissue_threshold`
-- never in JavaScript. The Streamlit widgets below recompute the filtered
subgraph on every interaction and hand it to `_subgraph_to_html`, which renders
it with the force-directed `ExposoGraph/map/map.{html,css,js}`
D3 renderer. That renderer no longer contains any filtering logic of its own;
see `docs/design/refactor_ui_map_viewer_plan.md` for the full rationale.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import streamlit as st
import streamlit.components.v1 as components

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
_FONTS_URL = (
    "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700"
    "&family=JetBrains+Mono:wght@400;500&display=swap"
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


def _subgraph_to_html(graph: dict[str, list[Any]]) -> str:
    """Assemble the map document from ``map/map.{html,css,js}`` with *graph*
    embedded as ``GRAPH_DATA`` (plus ``MAP_STYLE``)."""
    style = {
        "node_colors": NODE_COLORS,
        "node_labels": NODE_LABELS,
        "node_sizes": NODE_SIZES,
        "edge_colors": EDGE_COLORS,
        "edge_labels": EDGE_LABELS,
    }
    data_script = (
        f"const GRAPH_DATA = {json.dumps(_clean_for_js(graph), indent=2)};\n"
        f"const MAP_STYLE = {json.dumps(style, indent=2)};"
    )
    return (
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n'
        f"<script>\n{_D3_JS.read_text(encoding='utf-8')}\n</script>\n"
        f'<link href="{_FONTS_URL}" rel="stylesheet">\n'
        f"<style>\n{_MAP_CSS.read_text(encoding='utf-8')}</style>\n"
        "</head>\n<body>\n"
        f"{_MAP_HTML.read_text(encoding='utf-8')}"
        f"<script>\n{data_script}\n</script>\n"
        f"<script>\n{_MAP_JS.read_text(encoding='utf-8')}</script>\n"
        "</body>\n</html>\n"
    )


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

        html = _subgraph_to_html(subgraph)

    st.caption(
        f"{len(subgraph['nodes'])} nodes / {len(subgraph['edges'])} edges shown."
    )
    components.html(html, height=920, scrolling=True)
