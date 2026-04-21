"""Streamlit UI for ExposoGraph.

Run with:  streamlit run ExposoGraph/app.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

# Ensure the package is importable when run via `streamlit run ExposoGraph/app.py`
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from ExposoGraph import (
    ui_d3_viewer,
    ui_data,
    ui_extract,
    ui_manual,
    ui_map_viewer,
    ui_preview,
    ui_sidebar,
)
from ExposoGraph._app_shared import get_engine, get_repository
from ExposoGraph.branding import (
    APP_NAME,
    APP_TAGLINE,
    APP_VERSION,
    CONTACT_EMAIL,
    COPYRIGHT_HOLDER,
    DEVELOPED_BY,
)

# ── Page config ──────────────────────────────────────────────────────────

st.set_page_config(
    page_title=f"{APP_NAME} v{APP_VERSION}",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    [data-testid="stSidebar"] { min-width: 340px; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    </style>
    """,
    unsafe_allow_html=True,
)

# ── Session state bootstrap ──────────────────────────────────────────────

engine = get_engine()
repository = get_repository()
if "extract_text" not in st.session_state:
    st.session_state.extract_text = ""
if "project_name" not in st.session_state:
    st.session_state.project_name = "knowledge_graph"
if "revision_note" not in st.session_state:
    st.session_state.revision_note = ""

# ── Sidebar ──────────────────────────────────────────────────────────────

ui_sidebar.render(engine, repository)

# ── Main area tabs ───────────────────────────────────────────────────────

st.markdown(f"## {APP_NAME}")
st.caption(f"{APP_TAGLINE} · Version {APP_VERSION}")
st.caption(f"{DEVELOPED_BY} · {CONTACT_EMAIL} · Copyright {COPYRIGHT_HOLDER}")

tab_map, tab_extract, tab_manual, tab_preview, tab_d3_viewer, tab_data = st.tabs(
    ["Reference Map", "LLM Extract", "Manual Entry", "Graph Preview", "D3 HTML Viewer", "Raw Data"]
)

with tab_map:
    ui_map_viewer.render()

with tab_extract:
    ui_extract.render(engine)

with tab_manual:
    ui_manual.render(engine)

with tab_preview:
    ui_preview.render(engine)

with tab_d3_viewer:
    ui_d3_viewer.render(engine)

with tab_data:
    ui_data.render(engine)
