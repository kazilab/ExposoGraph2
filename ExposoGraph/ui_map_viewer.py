"""Tab: embedded D3 viewer for the checked-in reference map bundle."""

import streamlit as st
import streamlit.components.v1 as components
from pathlib import Path
from .exporter import bundle_to_html_string
from .engine import GraphEngine

_MAP_DIR = Path(__file__).resolve().parent / "map"

def render() -> None:
    """Render the D3 map frame using the backend graph engine intersection filter."""
    st.markdown("#### Reference Map")

    if "engine" not in st.session_state:
        st.error("Graph Engine is not initialized.")
        return
    master_engine: GraphEngine = st.session_state.engine

    selected_node_type = st.session_state.get("node_type_select", None)
    selected_edge_type = st.session_state.get("edge_type_select", None)
    selected_carcinogen = st.session_state.get("carcinogen_group_select", None)
    selected_tissue = st.session_state.get("tissue_select", None)
    tissue_threshold = st.session_state.get("tissue_threshold", None)

    # 2. Compute the pristine intersecting sub-network
    pruned_engine = master_engine.filter_by_criteria(
        node_type=selected_node_type,
        edge_type=selected_edge_type,
        carcinogen_class=selected_carcinogen,
        tissue=selected_tissue,
        min_tissue_weight=tissue_threshold,
    )

    # 3. Export the data layout dictionary matching what your frontend consumes
    d3_graph_data = pruned_engine.to_dict()

    # 4. Generate visual wrapper code
    index_path = _MAP_DIR / "index.html"
    if not index_path.exists():
        st.error("Bundled template index.html is missing.")
        return

    html = bundle_to_html_string(index_path, d3_graph_data)
    components.html(html, height=920, scrolling=True)
