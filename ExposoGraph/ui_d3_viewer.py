"""Tab: embedded interactive HTML viewer.

This wires the existing HTML export path (`exporter.to_interactive_html_string`)
into the Streamlit app, so users can view the standalone D3-style HTML
representation without downloading it first.
"""

from __future__ import annotations

import streamlit as st
import streamlit.components.v1 as components

from ._app_shared import slugify_project_name
from .config import GraphVisibility
from .engine import GraphEngine
from .exporter import to_interactive_html_string
from .graph_filters import graph_visibility_label


def render(engine: GraphEngine) -> None:
    st.markdown("#### D3 HTML Viewer")

    if engine.node_count == 0:
        st.info("No data yet — use LLM Extract or Manual Entry to add nodes.")
        return

    visibility = st.session_state.get("graph_visibility", GraphVisibility.ALL.value)
    st.caption(f"Visibility: {graph_visibility_label(visibility)}")

    with st.spinner("Rendering interactive HTML viewer..."):
        html = to_interactive_html_string(engine, visibility=visibility)

    components.html(html, height=840, scrolling=True)

    st.download_button(
        "Download HTML",
        data=html,
        file_name=f"{slugify_project_name(st.session_state.project_name)}_d3_viewer.html",
        mime="text/html",
        use_container_width=True,
    )

