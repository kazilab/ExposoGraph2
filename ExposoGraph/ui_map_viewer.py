"""Tab: interactive Reference Map viewer.

Filtering (node type, carcinogen group, tissue threshold) is decided here in
Python via `GraphEngine.map_viewer_subgraph`/`GraphEngine.dim_by_tissue_threshold`
-- never in JavaScript. The Streamlit widgets below recompute the filtered
subgraph on every interaction and hand it to `exporter.subgraph_to_html_string`,
which renders it against the existing `ExposoGraph/map/index.html` force-directed
D3 template. That template no longer contains any filtering logic of its own;
see `docs/design/refactor_ui_map_viewer_plan.md` for the full rationale.
"""

from __future__ import annotations

from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from .engine import GraphEngine
from .exporter import subgraph_to_html_string
from .models import NodeType

_MAP_DIR = Path(__file__).resolve().parent / "map"
_MAP_TEMPLATE = _MAP_DIR / "index.html"

# Substrate nodes are part of the underlying knowledge graph but must never
# appear in this viewer (they're intermediate bookkeeping, not something a
# reader needs to see on the map) -- so they're never offered as a selectable
# node-type option, and GraphEngine excludes them by default regardless.
_SELECTABLE_NODE_TYPES = [t.value for t in NodeType if t is not NodeType.SUBSTRATE]

_NO_TISSUE_OPTION = "(none — no tissue dimming)"


def render(engine: GraphEngine) -> None:
    """Render the Reference Map viewer, filtered/dimmed via `engine`."""
    st.markdown("#### Reference Map")
    st.caption(
        "Force-directed graph rendered from the loaded reference graph. "
        "Filters below run against the knowledge graph in Python; `Substrate` "
        "nodes are never shown."
    )

    if not _MAP_TEMPLATE.exists():
        st.error(f"Bundled map template is missing: {_MAP_TEMPLATE.name}")
        return

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
            "Carcinogen group(s)",
            options=carcinogen_groups,
            default=[],
            help=(
                "Restricts the map to nodes reachable via a directed path "
                "from any carcinogen in the selected group(s). Leave empty "
                "for no carcinogen-path restriction."
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
            help="Enzyme nodes (and their edges) below this weight are grayed out, not removed.",
            key="map_viewer_tissue_threshold",
        )

    with st.spinner("Rendering map..."):
        subgraph = engine.map_viewer_subgraph(
            node_types=selected_node_types or None,
            carcinogen_groups=selected_groups or None,
        )
        if selected_tissue != _NO_TISSUE_OPTION:
            subgraph = engine.dim_by_tissue_threshold(subgraph, selected_tissue, threshold)

        html = subgraph_to_html_string(subgraph, template_path=_MAP_TEMPLATE)

    st.caption(f"{len(subgraph['nodes'])} nodes / {len(subgraph['edges'])} edges shown.")
    components.html(html, height=920, scrolling=True)
