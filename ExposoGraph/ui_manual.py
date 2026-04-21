"""Tab 2: Manual Node & Edge Entry."""

from __future__ import annotations

import streamlit as st

from ._app_shared import build_curation_record, build_provenance_record
from .engine import GraphEngine
from .models import (
    CurationConfidence,
    CurationStatus,
    Edge,
    EdgeType,
    Node,
    NodeType,
)


def render(engine: GraphEngine) -> None:
    """Render the Manual Entry tab."""
    st.markdown("#### Add Nodes & Edges Manually")

    node_tab, edge_tab = st.tabs(["Add Node", "Add Edge"])

    with node_tab:
        _render_add_node(engine)

    with edge_tab:
        _render_add_edge(engine)


def _render_add_node(engine: GraphEngine) -> None:
    with st.form("add_node", clear_on_submit=True):
        st.markdown("**New Node**")
        c1, c2 = st.columns(2)
        with c1:
            n_label = st.text_input("Label *", placeholder="e.g. CYP1A1")
            n_type = st.selectbox("Type *", [t.value for t in NodeType])
        with c2:
            n_id = st.text_input("ID (auto from label)", placeholder="CYP1A1")
            n_detail = st.text_input("Detail", placeholder="Short description")

        c3, c4, c5 = st.columns(3)
        with c3:
            n_group = st.text_input("Group / Class", placeholder="PAH, DNA Repair (BER), …")
            n_iarc = st.selectbox("IARC", [None, "Group 1", "Group 2A", "Group 2B", "Group 3"])
        with c4:
            n_phase = st.selectbox("Phase", [None, "I", "II", "III"])
            n_role = st.selectbox(
                "Role",
                [None, "Activation", "Detoxification", "Mixed", "Transport", "Repair"],
            )
        with c5:
            n_reactivity = st.selectbox("Reactivity", [None, "High", "Intermediate", "Low"])

        with st.expander("Provenance & Curation", expanded=False):
            a1, a2 = st.columns(2)
            with a1:
                n_source_db = st.text_input("Source DB", placeholder="NCBI Gene; GTEx; ClinPGx")
                n_record_id = st.text_input("Record ID / Accession", placeholder="e.g. 1543 or CYP1A1")
                n_pmid = st.text_input("PMID", placeholder="e.g. 29134625")
                n_citation = st.text_input("Citation", placeholder="Short citation or paper title")
                n_variant = st.text_input("Variant", placeholder="e.g. CYP1A1*2C")
            with a2:
                n_tissue = st.text_input("Tissue", placeholder="e.g. lung, liver, bladder")
                n_url = st.text_input("URL", placeholder="https://...")
                n_phenotype = st.text_input("Phenotype", placeholder="e.g. poor metabolizer")
                n_activity_score = st.text_input("Activity score", placeholder="e.g. 1.5")
            n_evidence = st.text_area("Evidence", placeholder="Short provenance or evidence note", height=80)

            curation_col1, curation_col2 = st.columns(2)
            with curation_col1:
                n_status = st.selectbox(
                    "Status",
                    [status.value for status in CurationStatus],
                    index=0,
                )
                n_confidence = st.selectbox(
                    "Confidence",
                    [""] + [confidence.value for confidence in CurationConfidence],
                )
                n_curator = st.text_input("Curator", placeholder="Name or initials")
            with curation_col2:
                n_reviewed_by = st.text_input("Reviewed by", placeholder="Reviewer name")
                n_reviewed_at = st.text_input("Reviewed at", placeholder="YYYY-MM-DD")
                n_curation_notes = st.text_area(
                    "Curation notes",
                    placeholder="Why this record is draft/reviewed/approved",
                    height=80,
                )

        if st.form_submit_button("Add Node", type="primary"):
            if not n_label:
                st.error("Label is required")
            elif n_activity_score:
                try:
                    activity_score = float(n_activity_score)
                except ValueError:
                    st.error("Activity score must be numeric")
                    activity_score = None
            else:
                activity_score = None
            if n_label and (not n_activity_score or activity_score is not None):
                provenance = build_provenance_record(
                    source_db=n_source_db,
                    record_id=n_record_id,
                    evidence=n_evidence,
                    pmid=n_pmid,
                    tissue=n_tissue,
                    citation=n_citation,
                    url=n_url,
                )
                curation = build_curation_record(
                    status=n_status,
                    confidence=n_confidence,
                    curator=n_curator,
                    reviewed_by=n_reviewed_by,
                    reviewed_at=n_reviewed_at,
                    notes=n_curation_notes,
                )
                node = Node(
                    id=n_id or n_label.replace(" ", "_"),
                    label=n_label,
                    type=NodeType(n_type),
                    detail=n_detail,
                    group=n_group or None,
                    iarc=n_iarc,
                    phase=n_phase,
                    role=n_role,
                    reactivity=n_reactivity,
                    source_db=n_source_db or None,
                    evidence=n_evidence or None,
                    pmid=n_pmid or None,
                    tissue=n_tissue or None,
                    variant=n_variant or None,
                    phenotype=n_phenotype or None,
                    activity_score=activity_score,
                    provenance=provenance,
                    curation=curation,
                )
                engine.add_node(node)
                st.success(f"Added node **{node.label}** ({node.type.value})")
                st.rerun()


def _render_add_edge(engine: GraphEngine) -> None:
    node_ids = sorted(engine.G.nodes) if engine.node_count > 0 else []

    with st.form("add_edge", clear_on_submit=True):
        st.markdown("**New Edge**")
        if not node_ids:
            st.info("Add at least two nodes first.")

        c1, c2, c3 = st.columns(3)
        with c1:
            e_source = st.selectbox("Source *", node_ids or [""])
        with c2:
            e_target = st.selectbox("Target *", node_ids or [""])
        with c3:
            e_type = st.selectbox("Type *", [t.value for t in EdgeType])

        c4, c5 = st.columns(2)
        with c4:
            e_label = st.text_input("Label", placeholder="e.g. epoxidation")
        with c5:
            carcinogen_ids = [
                n for n in node_ids
                if (nd := engine.get_node(n)) is not None and nd.get("type") == "Carcinogen"
            ]
            e_carcin = st.selectbox("Carcinogen (context)", [""] + carcinogen_ids)

        with st.expander("Provenance & Curation", expanded=False):
            e_source_db = st.text_input("Source DB", placeholder="NCBI Gene; CTD; ClinPGx")
            e_record_id = st.text_input("Record ID / Accession", placeholder="e.g. rs2228001 or pathway accession")
            e_pmid = st.text_input("PMID", placeholder="e.g. 25786862")
            e_tissue = st.text_input("Tissue", placeholder="e.g. liver, lung")
            e_citation = st.text_input("Citation", placeholder="Short citation or paper title")
            e_url = st.text_input("URL", placeholder="https://...")
            e_evidence = st.text_area("Evidence", placeholder="Short provenance or evidence note", height=80)
            curation_col1, curation_col2 = st.columns(2)
            with curation_col1:
                e_status = st.selectbox(
                    "Status",
                    [status.value for status in CurationStatus],
                    index=0,
                )
                e_confidence = st.selectbox(
                    "Confidence",
                    [""] + [confidence.value for confidence in CurationConfidence],
                )
                e_curator = st.text_input("Curator", placeholder="Name or initials")
            with curation_col2:
                e_reviewed_by = st.text_input("Reviewed by", placeholder="Reviewer name")
                e_reviewed_at = st.text_input("Reviewed at", placeholder="YYYY-MM-DD")
                e_curation_notes = st.text_area(
                    "Curation notes",
                    placeholder="Why this edge is draft/reviewed/approved",
                    height=80,
                )

        if st.form_submit_button("Add Edge", type="primary"):
            if not e_source or not e_target:
                st.error("Source and target are required")
            else:
                provenance = build_provenance_record(
                    source_db=e_source_db,
                    record_id=e_record_id,
                    evidence=e_evidence,
                    pmid=e_pmid,
                    tissue=e_tissue,
                    citation=e_citation,
                    url=e_url,
                )
                curation = build_curation_record(
                    status=e_status,
                    confidence=e_confidence,
                    curator=e_curator,
                    reviewed_by=e_reviewed_by,
                    reviewed_at=e_reviewed_at,
                    notes=e_curation_notes,
                )
                edge = Edge(
                    source=e_source,
                    target=e_target,
                    type=EdgeType(e_type),
                    label=e_label or None,
                    carcinogen=e_carcin or None,
                    source_db=e_source_db or None,
                    evidence=e_evidence or None,
                    pmid=e_pmid or None,
                    tissue=e_tissue or None,
                    provenance=provenance,
                    curation=curation,
                )
                try:
                    engine.add_edge(edge)
                    st.success(f"Added edge **{e_source}** → **{e_target}** ({e_type})")
                    st.rerun()
                except ValueError as exc:
                    st.error(str(exc))
