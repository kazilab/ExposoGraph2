"""Sidebar: import, revision history, gene panels, activity scores, export."""

from __future__ import annotations

import streamlit as st

from ._app_shared import (
    APP_MODE,
    PERSISTENCE_ENABLED,
    PROJECTS_DIR,
    REPOSITORY_PATH,
    load_into_engine,
    parse_uploaded_graph,
    relative_path,
    revision_label,
    saved_project_paths,
    slugify_project_name,
)
from .branding import (
    APP_NAME,
    APP_TAGLINE,
    APP_VERSION,
    CONTACT_EMAIL,
    COPYRIGHT_HOLDER,
    DEVELOPED_BY,
)
from .config import AppMode, GraphVisibility
from .engine import GraphEngine
from .exporter import (
    parse_graph_artifact,
    to_interactive_html,
    to_interactive_html_string,
    to_plotly_html_string,
)
from .graph_filters import filtered_engine, graph_visibility_label
from .reference_data import (
    ACTIVITY_SCORES,
    build_full_panel,
    build_tier1_panel,
    build_tier2_panel,
    get_activity_score_metadata,
    get_activity_score_references,
    get_activity_scores,
)
from .storage import GraphRepository


def render(engine: GraphEngine, repository: GraphRepository | None) -> None:
    """Render the full sidebar."""
    with st.sidebar:
        st.markdown(f"### {APP_NAME}")
        st.caption(APP_TAGLINE)
        st.caption(f"Version {APP_VERSION} · {DEVELOPED_BY}")
        st.caption(f"Contact: {CONTACT_EMAIL} · Copyright {COPYRIGHT_HOLDER}")
        st.caption(
            f"**{engine.node_count}** nodes · **{engine.edge_count}** edges"
        )
        st.divider()
        _render_graph_filters(engine)
        active_visibility = st.selectbox(
            "Graph visibility",
            options=[visibility.value for visibility in GraphVisibility],
            key="graph_visibility",
            format_func=graph_visibility_label,
            help="Controls Preview, Raw Data, and Export views without changing the stored graph.",
        )
        visible_engine = filtered_engine(engine, active_visibility)
        if active_visibility != GraphVisibility.ALL.value:
            st.caption(
                f"Current view: **{visible_engine.node_count}** nodes · "
                f"**{visible_engine.edge_count}** edges"
            )
        if APP_MODE == AppMode.STATELESS:
            st.info(
                "Mode: stateless. User graphs are not saved on the server. "
                "Download the interactive HTML file to keep your work."
            )
        else:
            st.caption(f"Mode: local persistence (`{relative_path(REPOSITORY_PATH)}`)")

        st.divider()
        _render_import(engine)
        st.divider()
        _render_project_name()
        if PERSISTENCE_ENABLED and repository is not None:
            _render_revision_history(engine, repository)
            _render_html_snapshots(engine)
            st.divider()
        _render_gene_panels(engine)
        st.divider()
        _render_activity_scores()
        st.divider()
        _render_export(engine)
        st.divider()
        _render_validation(engine)
        st.divider()
        if st.button("Clear graph", type="secondary", use_container_width=True):
            engine.clear()
            st.rerun()


def _render_import(engine: GraphEngine) -> None:
    st.markdown("##### Import Existing Data")
    replace_import = st.checkbox(
        "Replace current graph on import",
        value=engine.node_count == 0,
        help="When enabled, imported data clears the current graph first.",
    )
    st.caption("Upload a saved JSON, HTML, or graph-data.js export. There is no canonical viewer file in the repo.")
    uploaded = st.file_uploader(
        "Upload graph",
        type=["json", "html", "js"],
        label_visibility="collapsed",
    )
    if uploaded is not None:
        try:
            kg = parse_uploaded_graph(uploaded.name, uploaded.read().decode("utf-8"))
            warnings = load_into_engine(engine, kg, replace=replace_import)
            message = f"Loaded {len(kg.nodes)} nodes, {len(kg.edges)} edges"
            if warnings:
                message += f" with {len(warnings)} warning(s)"
            st.success(message)
            st.rerun()
        except Exception as exc:
            st.error(f"Upload failed: {exc}")


def _render_project_name() -> None:
    st.text_input(
        "Project name",
        key="project_name",
        help="Used as the human-readable project name and the base name for exports.",
    )


def _render_revision_history(engine: GraphEngine, repository: GraphRepository) -> None:
    st.text_input(
        "Revision note",
        key="revision_note",
        help="Optional note describing what changed in this revision.",
    )

    st.divider()

    st.markdown("##### Revision History")
    replace_import = st.session_state.get("Replace current graph on import", engine.node_count == 0)
    visibility = st.session_state.get("graph_visibility", GraphVisibility.ALL.value)
    visible_engine = filtered_engine(engine, visibility)
    st.caption(f"Save scope: {graph_visibility_label(visibility)}")

    col_r1, col_r2 = st.columns(2)
    with col_r1:
        if st.button("Save revision", use_container_width=True):
            if visible_engine.node_count == 0:
                st.error("Current revision view is empty")
            else:
                try:
                    saved = repository.save_engine(
                        graph_key=slugify_project_name(st.session_state.project_name),
                        graph_name=st.session_state.project_name.strip() or "knowledge_graph",
                        engine=engine,
                        visibility=visibility,
                        note=st.session_state.revision_note.strip() or None,
                    )
                    st.success(
                        f"Saved revision r{saved.revision_number} "
                        f"({graph_visibility_label(saved.visibility)}) to `{REPOSITORY_PATH.name}`"
                    )
                    st.session_state.revision_note = ""
                    st.rerun()
                except Exception as exc:
                    st.error(f"Revision save failed: {exc}")
    with col_r2:
        db_graphs = repository.list_graphs()
        selected_db_graph = st.selectbox(
            "Stored project",
            options=[""] + [graph.graph_key for graph in db_graphs],
            format_func=lambda key: next(
                (
                    f"{graph.graph_name} (r{graph.revision_number})"
                    for graph in db_graphs
                    if graph.graph_key == key
                ),
                "Select project",
            ),
            label_visibility="collapsed",
        )
        if st.button("Load latest", use_container_width=True):
            if not selected_db_graph:
                st.error("Select a stored project first")
            else:
                revision = repository.get_latest_revision(selected_db_graph)
                if revision is None:
                    st.error("No revisions found for that project")
                else:
                    warnings = load_into_engine(
                        engine,
                        revision.to_knowledge_graph(),
                        replace=replace_import,
                    )
                    message = f"Loaded {revision.graph_name} r{revision.revision_number}"
                    if warnings:
                        message += f" with {len(warnings)} warning(s)"
                    st.success(message)
                    st.rerun()

    if db_graphs:
        selected_history_graph = st.selectbox(
            "Revision history",
            options=[""] + [graph.graph_key for graph in db_graphs],
            format_func=lambda key: next(
                (graph.graph_name for graph in db_graphs if graph.graph_key == key),
                "Select project",
            ),
        )
        if selected_history_graph:
            revisions = repository.list_revisions(selected_history_graph)
            revisions_by_id = {revision.revision_id: revision for revision in revisions}
            selected_revision_id = st.selectbox(
                "Revision",
                options=[revision.revision_id for revision in revisions],
                format_func=lambda revision_id: revision_label(revision_id, revisions_by_id),
            )
            selected_revision = repository.get_revision(selected_revision_id)
            if selected_revision is not None:
                st.caption(f"Database: `{relative_path(REPOSITORY_PATH)}`")
                col_db1, col_db2 = st.columns(2)
                with col_db1:
                    if st.button("Load selected revision", use_container_width=True):
                        warnings = load_into_engine(
                            engine,
                            selected_revision.to_knowledge_graph(),
                            replace=replace_import,
                        )
                        message = (
                            f"Loaded {selected_revision.graph_name} "
                            f"r{selected_revision.revision_number}"
                        )
                        if warnings:
                            message += f" with {len(warnings)} warning(s)"
                        st.success(message)
                        st.rerun()
                with col_db2:
                    st.download_button(
                        "Download revision HTML",
                        data=selected_revision.html,
                        file_name=(
                            f"{selected_revision.graph_key}_r"
                            f"{selected_revision.revision_number}.html"
                        ),
                        mime="text/html",
                        use_container_width=True,
                    )

    st.divider()


def _render_html_snapshots(engine: GraphEngine) -> None:
    st.markdown("##### HTML Files")
    visibility = st.session_state.get("graph_visibility", GraphVisibility.ALL.value)
    visible_engine = filtered_engine(engine, visibility)
    st.caption(f"Snapshot scope: {graph_visibility_label(visibility)}")
    col_s1, col_s2 = st.columns(2)
    with col_s1:
        if st.button("Save HTML snapshot", use_container_width=True):
            if visible_engine.node_count == 0:
                st.error("Current snapshot view is empty")
            else:
                project_path = PROJECTS_DIR / f"{slugify_project_name(st.session_state.project_name)}.html"
                to_interactive_html(
                    engine,
                    project_path,
                    visibility=visibility,
                )
                st.success(f"Saved `{project_path.name}`")
    with col_s2:
        saved_projects = saved_project_paths()
        selected_project = st.selectbox(
            "Saved",
            options=[""] + [path.name for path in saved_projects],
            label_visibility="collapsed",
        )
        replace_import = st.session_state.get("Replace current graph on import", engine.node_count == 0)
        if st.button("Load HTML snapshot", use_container_width=True):
            if not selected_project:
                st.error("Select a saved project first")
            else:
                project_path = PROJECTS_DIR / selected_project
                kg = parse_graph_artifact(project_path)
                warnings = load_into_engine(engine, kg, replace=replace_import)
                message = f"Loaded `{project_path.name}`"
                if warnings:
                    message += f" with {len(warnings)} warning(s)"
                st.success(message)
                st.rerun()


def _render_gene_panels(engine: GraphEngine) -> None:
    st.markdown("##### Reference Gene Panels")
    col_p1, col_p2, col_p3 = st.columns(3)
    with col_p1:
        if st.button("Tier 1 (13)", use_container_width=True, help="Core CYP/GST/NAT/UGT enzymes"):
            kg = build_tier1_panel()
            engine.merge(kg)
            st.success(f"Loaded {len(kg.nodes)} Tier 1 genes")
            st.rerun()
    with col_p2:
        if st.button("Tier 2 (23)", use_container_width=True, help="Extended hormone, transport, and DNA-repair panel"):
            kg = build_tier2_panel()
            engine.merge(kg)
            st.success(f"Loaded {len(kg.nodes)} Tier 2 genes")
            st.rerun()
    with col_p3:
        if st.button("All (36)", use_container_width=True, help="Full Tier 1 + Tier 2 panel"):
            kg = build_full_panel()
            engine.merge(kg)
            st.success(f"Loaded {len(kg.nodes)} genes")
            st.rerun()


def _render_activity_scores() -> None:
    st.markdown("##### Activity Score Lookup")
    activity_gene = st.selectbox("Gene", [""] + sorted(ACTIVITY_SCORES.keys()))
    if activity_gene:
        st.dataframe(get_activity_scores(activity_gene), use_container_width=True, hide_index=True)
        activity_meta = get_activity_score_metadata(activity_gene) or {}
        if activity_meta:
            st.caption(
                f"Evidence basis: {activity_meta.get('evidence_basis', 'Unspecified')}. "
                f"{activity_meta.get('note', '')}"
            )
        activity_refs = get_activity_score_references(activity_gene) or []
        if activity_refs:
            with st.expander("Activity score references", expanded=False):
                ref_rows = [
                    {
                        "source_db": ref.get("source_db", ""),
                        "record_id": ref.get("record_id", ""),
                        "pmid": ref.get("pmid", ""),
                        "citation": ref.get("citation", ""),
                        "url": ref.get("url", ""),
                    }
                    for ref in activity_refs
                ]
                st.dataframe(ref_rows, use_container_width=True, hide_index=True)


def _render_export(engine: GraphEngine) -> None:
    st.markdown("##### Export")
    visibility = st.session_state.get("graph_visibility", GraphVisibility.ALL.value)
    export_engine = filtered_engine(engine, visibility)
    html_export = ""
    html_export_error = ""
    plotly_html_export = ""
    plotly_html_error = ""
    if export_engine.node_count > 0:
        try:
            html_export = to_interactive_html_string(export_engine)
        except Exception as exc:
            html_export_error = str(exc)
        try:
            plotly_html_export = to_plotly_html_string(export_engine)
        except Exception as exc:
            plotly_html_error = str(exc)
    st.caption(
        "Export scope: "
        f"{graph_visibility_label(visibility)} "
        f"({export_engine.node_count} nodes, {export_engine.edge_count} edges)"
    )
    if PERSISTENCE_ENABLED:
        st.caption("Use HTML Files below to save a named local snapshot.")
    else:
        st.caption("Server-side file export disabled in stateless mode.")
    if export_engine.node_count > 0:
        col_e1, col_e2 = st.columns(2)
        if not html_export_error:
            with col_e1:
                st.download_button(
                    "Download App HTML",
                    data=html_export,
                    file_name=f"{slugify_project_name(st.session_state.project_name)}.html",
                    mime="text/html",
                    use_container_width=True,
                    help="Parseable ExposoGraph HTML export for reload/revisions.",
                )
        if not plotly_html_error:
            with col_e2:
                st.download_button(
                    "Download Plotly HTML",
                    data=plotly_html_export,
                    file_name=f"{slugify_project_name(st.session_state.project_name)}_plotly.html",
                    mime="text/html",
                    use_container_width=True,
                    help="Standalone Plotly viewer optimized for interactive exploration.",
                )
    if html_export_error:
        st.caption(f"HTML export unavailable: {html_export_error}")
    if plotly_html_error:
        st.caption(f"Plotly HTML export unavailable: {plotly_html_error}")


def _render_validation(engine: GraphEngine) -> None:
    errors = engine.validate()
    if errors:
        st.warning(f"{len(errors)} validation issue(s)")
        for e in errors:
            st.caption(f"⚠ {e}")
    elif engine.node_count > 0:
        st.success("Graph valid")

def _render_graph_filters(engine: GraphEngine) -> None:
    """Render universal graph filters synchronized with Python intersection logic."""
    st.markdown("##### Active Graph Filters")
    
    # ── 1. Carcinogen Class Selection ──────────────────────────────────────
    # Dynamically look up all unique node group designations from memory
    all_nodes = [data for _, data in engine.G.nodes(data=True)]
    carcinogen_groups = sorted(list({
        data.get("group") for data in all_nodes 
        if data.get("type") == "Carcinogen" and data.get("group")
    }))
    
    selected_carcinogen = st.selectbox(
        "Carcinogen Class",
        options=["All Groups"] + carcinogen_groups,
        format_func=lambda x: x.replace("_", " ") if x != "All Groups" else x,
        help="Returns neighbors-of-neighbors of carcinogens in this group."
    )
    # Map the "All Groups" string back to None for the intersection logic
    st.session_state["carcinogen_group_select"] = (
        None if selected_carcinogen == "All Groups" else selected_carcinogen
    )

    # ── 2. Node & Edge Type Toggles ────────────────────────────────────────
    # Extract unique categories present across the graph
    node_types = sorted(list({data.get("type") for _, data in engine.G.nodes(data=True) if data.get("type")}))
    edge_types = sorted(list({data.get("type") for _, _, _, data in engine.G.edges(keys=True, data=True) if data.get("type")}))

    selected_node_type = st.selectbox("Node Types", options=["All Nodes"] + node_types)
    st.session_state["node_type_select"] = (
        None if selected_node_type == "All Nodes" else selected_node_type
    )

    selected_edge_type = st.selectbox("Edge Types", options=["All Edges"] + edge_types)
    st.session_state["edge_type_select"] = (
        None if selected_edge_type == "All Edges" else selected_edge_type
    )

    # ── 3. Tissue & Weight Metrics ─────────────────────────────────────────
    # Scrape unique tissues from data blocks
    discovered_tissues = set()
    for _, d in engine.G.nodes(data=True):
        t_field = d.get("tissue") or d.get("tissues")
        if isinstance(t_field, dict):
            discovered_tissues.update(t_field.keys())
        elif isinstance(t_field, str):
            discovered_tissues.add(t_field)

    selected_tissue = st.selectbox("Filter by Tissue Environment", options=["All Tissues"] + sorted(list(discovered_tissues)))
    
    if selected_tissue != "All Tissues":
        st.session_state["tissue_select"] = selected_tissue
        st.session_state["tissue_threshold"] = st.slider(
            "Minimum Tissue weight", 
            min_value=0.0, 
            max_value=1.0, 
            value=0.05, 
            step=0.05
        )
    else:
        st.session_state["tissue_select"] = None
        st.session_state["tissue_threshold"] = None
