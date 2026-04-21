"""Tab 4: Raw Data View."""

from __future__ import annotations

import streamlit as st

from ._app_shared import (
    EDGE_COLORS,
    EDGE_SEARCH_FIELDS,
    NODE_COLORS,
    NODE_SEARCH_FIELDS,
    matches_query,
)
from .config import GraphVisibility
from .engine import GraphEngine
from .graph_filters import filter_knowledge_graph, graph_visibility_label


def render(engine: GraphEngine) -> None:
    """Render the Raw Data tab."""
    st.markdown("#### Current Graph Data")

    if engine.node_count == 0:
        st.info("Graph is empty.")
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
        data_query = st.text_input(
            "Search raw data",
            key="data_query",
            placeholder="Search nodes, edges, evidence, PMIDs...",
        ).strip()
    with col_f2:
        data_node_types = st.multiselect(
            "Node types",
            node_type_options,
            default=node_type_options,
            key="data_node_types",
        )
    with col_f3:
        data_edge_types = st.multiselect(
            "Edge types",
            edge_type_options,
            default=edge_type_options,
            key="data_edge_types",
        )

    filtered_nodes = [
        (node_id, data)
        for node_id, data in all_nodes
        if data.get("type", "?") in data_node_types
        and matches_query(data, data_query, NODE_SEARCH_FIELDS)
    ]
    filtered_edges = [
        (u, v, edge_key, data)
        for u, v, edge_key, data in all_edges
        if data.get("type", "?") in data_edge_types
        and matches_query(data, data_query, EDGE_SEARCH_FIELDS)
    ]

    node_rows = [
        {
            "id": data.get("id"),
            "label": data.get("label"),
            "type": data.get("type"),
            "detail": data.get("detail", ""),
            "match_status": data.get("match_status", ""),
            "origin": data.get("origin", ""),
            "status": (data.get("curation") or {}).get("status", ""),
            "confidence": (data.get("curation") or {}).get("confidence", ""),
            "provenance_count": len(data.get("provenance") or []),
            "source_db": data.get("source_db", ""),
            "tissue": data.get("tissue", ""),
        }
        for _, data in filtered_nodes
    ]
    edge_rows = [
        {
            "source": u,
            "target": v,
            "type": data.get("type"),
            "label": data.get("label", ""),
            "match_status": data.get("match_status", ""),
            "origin": data.get("origin", ""),
            "status": (data.get("curation") or {}).get("status", ""),
            "confidence": (data.get("curation") or {}).get("confidence", ""),
            "provenance_count": len(data.get("provenance") or []),
            "carcinogen": data.get("carcinogen", ""),
            "pmid": data.get("pmid", ""),
        }
        for u, v, _, data in filtered_edges
    ]

    st.caption(
        f"Filtered to {len(filtered_nodes)} nodes and {len(filtered_edges)} edges"
    )
    col_n, col_e = st.columns(2)
    with col_n:
        st.markdown(f"**Nodes ({len(filtered_nodes)})**")
        if node_rows:
            st.dataframe(node_rows, use_container_width=True, hide_index=True)
        for _, data in filtered_nodes:
            ntype = data.get("type", "?")
            color = NODE_COLORS.get(ntype, "#888")
            with st.expander(f":{color[1:]}[●] {data.get('label', data.get('id'))} — {ntype}"):
                st.json(data)
                if st.button("Delete", key=f"del_node_{data['id']}"):
                    engine.remove_node(data["id"])
                    st.rerun()

    with col_e:
        st.markdown(f"**Edges ({len(filtered_edges)})**")
        if edge_rows:
            st.dataframe(edge_rows, use_container_width=True, hide_index=True)
        for i, (u, v, edge_key, data) in enumerate(filtered_edges):
            etype = data.get("type", "?")
            color = EDGE_COLORS.get(etype, "#555")
            lbl = data.get("label", etype)
            with st.expander(f":{color[1:]}[→] {u} → {v}  ({lbl})"):
                st.json(data)
                if st.button("Delete", key=f"del_edge_{i}_{edge_key}"):
                    engine.remove_edge(u, v, edge_key)
                    st.rerun()
