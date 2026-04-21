"""Tab: embedded D3 viewer for the checked-in reference map bundle."""

from __future__ import annotations

from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from .exporter import bundle_to_html_string

_MAP_DIR = Path(__file__).resolve().parent / "map"


def render() -> None:
    """Render the bundled D3 map shipped in ``ExposoGraph/map``."""
    st.markdown("#### Reference Map")
    st.caption("Bundled D3 viewer from `ExposoGraph/map/index.html`.")

    index_path = _MAP_DIR / "index.html"
    graph_data_path = _MAP_DIR / "graph-data.js"
    missing = [str(path.name) for path in (index_path, graph_data_path) if not path.exists()]
    if missing:
        st.error(f"Bundled map assets are missing: {', '.join(missing)}")
        return

    html = bundle_to_html_string(index_path, graph_data_path)
    components.html(html, height=920, scrolling=True)
