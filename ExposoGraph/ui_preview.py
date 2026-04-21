"""Tab 3: Interactive Graph Preview."""

from __future__ import annotations

import streamlit as st

try:
    from streamlit_agraph import Config, agraph
    from streamlit_agraph import Edge as AEdge
    from streamlit_agraph import Node as ANode
except ModuleNotFoundError as exc:  # pragma: no cover - depends on optional extra
    Config = None
    AEdge = None
    ANode = None
    agraph = None
    _AGRAPH_IMPORT_ERROR: ModuleNotFoundError | None = exc
else:
    _AGRAPH_IMPORT_ERROR = None

from ._app_shared import (
    EDGE_COLORS,
    EDGE_SEARCH_FIELDS,
    NODE_COLORS,
    NODE_SEARCH_FIELDS,
    edge_tooltip,
    matches_query,
    node_tooltip,
)
from .config import GraphVisibility
from .engine import GraphEngine
from .graph_filters import filter_knowledge_graph, graph_visibility_label


def render(engine: GraphEngine) -> None:
    """Render the Graph Preview tab."""
    st.markdown("#### Graph Preview")

    if _AGRAPH_IMPORT_ERROR is not None:
        st.error(
            "Graph Preview requires the optional dependency `streamlit-agraph`. "
            "Install it with `pip install streamlit-agraph` or "
            "`pip install ExposoGraph[streamlit]`."
        )
        return

    if engine.node_count == 0:
        st.info("No data yet — use LLM Extract or Manual Entry to add nodes.")
        return

    visibility = st.session_state.get("graph_visibility", GraphVisibility.ALL.value)
    visible_graph = filter_knowledge_graph(engine.to_knowledge_graph(), visibility)
    all_nodes = [
        (node.id, node.model_dump(exclude_none=True, mode="json"))
        for node in visible_graph.nodes
    ]
    all_edges = [
        (
            edge.source,
            edge.target,
            f"{edge.source}-{edge.type.value}-{edge.target}",
            edge.model_dump(exclude_none=True, mode="json"),
        )
        for edge in visible_graph.edges
    ]
    node_type_options = sorted({data.get("type", "?") for _, data in all_nodes})
    edge_type_options = sorted({data.get("type", "?") for _, _, _, data in all_edges})

    st.caption(f"Visibility: {graph_visibility_label(visibility)}")

    col_f1, col_f2, col_f3 = st.columns(3)
    with col_f1:
        preview_query = st.text_input(
            "Search graph",
            key="preview_query",
            placeholder="label, id, detail, PMID, tissue, evidence...",
        ).strip()
    with col_f2:
        preview_node_types = st.multiselect(
            "Node types",
            node_type_options,
            default=node_type_options,
            key="preview_node_types",
        )
    with col_f3:
        preview_edge_types = st.multiselect(
            "Edge types",
            edge_type_options,
            default=edge_type_options,
            key="preview_edge_types",
        )

    allowed_node_ids = {
        node_id
        for node_id, data in all_nodes
        if data.get("type", "?") in preview_node_types
    }
    edge_match_node_ids: set[str] = set()
    if preview_query:
        for u, v, _, data in all_edges:
            if (
                data.get("type", "?") in preview_edge_types
                and u in allowed_node_ids
                and v in allowed_node_ids
                and matches_query(data, preview_query, EDGE_SEARCH_FIELDS)
            ):
                edge_match_node_ids.update((u, v))

    if preview_query:
        visible_node_ids = {
            node_id
            for node_id, data in all_nodes
            if node_id in allowed_node_ids
            and (
                matches_query(data, preview_query, NODE_SEARCH_FIELDS)
                or node_id in edge_match_node_ids
            )
        }
    else:
        visible_node_ids = allowed_node_ids

    filtered_nodes = [
        (node_id, data)
        for node_id, data in all_nodes
        if node_id in visible_node_ids
    ]
    filtered_edges = [
        (u, v, edge_key, data)
        for u, v, edge_key, data in all_edges
        if data.get("type", "?") in preview_edge_types
        and u in visible_node_ids
        and v in visible_node_ids
    ]

    st.caption(
        f"Showing {len(filtered_nodes)} of {engine.node_count} nodes and "
        f"{len(filtered_edges)} of {engine.edge_count} edges"
    )

    if not filtered_nodes:
        st.info("No nodes match the current preview filters.")
        return

    agraph_nodes = []
    agraph_edges = []

    for _, data in filtered_nodes:
        ntype = data.get("type", "Enzyme")
        color = NODE_COLORS.get(ntype, "#888")
        agraph_nodes.append(
            ANode(
                id=data["id"],
                label=data.get("label", data["id"]),
                color=color,
                size=25 if ntype == "Carcinogen" else 18,
                title=node_tooltip(data),
            )
        )

    for u, v, _, data in filtered_edges:
        etype = data.get("type", "ACTIVATES")
        color = EDGE_COLORS.get(etype, "#555")
        agraph_edges.append(
            AEdge(
                source=u,
                target=v,
                label=data.get("label", etype),
                color=color,
                title=edge_tooltip(data),
                width=2,
            )
        )

    config = Config(
        width=1200,
        height=700,
        directed=True,
        physics=True,
        hierarchical=False,
        nodeHighlightBehavior=True,
        highlightColor="#5cc6d0",
        collapsible=False,
        node={"labelProperty": "label"},
        link={"labelProperty": "label", "renderLabel": False},
    )

    agraph(nodes=agraph_nodes, edges=agraph_edges, config=config)

    st.caption(
        f"{len(filtered_nodes)} nodes · {len(filtered_edges)} edges · "
        "Drag to rearrange, scroll to zoom"
    )
