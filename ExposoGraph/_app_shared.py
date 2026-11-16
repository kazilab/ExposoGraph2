"""Shared constants, helpers, and session-state bootstrap for the Streamlit UI."""

from __future__ import annotations

from typing import Any

import streamlit as st

from .engine import GraphEngine
from .validation.models import (
    CurationConfidence,
    CurationRecord,
    CurationStatus,
    ProvenanceRecord,
)

# ── Node/edge styles ──────────────────────────────────────────────────────
# Single source for the colors, legend labels, and sizes used by the
# Reference Map (passed into the map renderer) and the Raw Data tab.

NODE_COLORS = {
    "Carcinogen": "#e05565",
    "Enzyme": "#4f98a3",
    "Receptor": "#c74fb0",
    "Gene": "#3d8b8b",
    "Metabolite": "#e8945a",
    "DNA_Adduct": "#a86fdf",
    "Pathway": "#5591c7",
    "Tissue": "#c2855a",
}
NODE_LABELS = {
    "Carcinogen": "Carcinogens (IARC-classified)",
    "Enzyme": "Enzymes (Phase I/II/III + Repair)",
    "Receptor": "Receptors (Nuclear/Ligand-Activated Sensors)",
    "Gene": "Genes (Receptors & Variant Loci)",
    "Metabolite": "Metabolites (Intermediates & Conjugates)",
    "DNA_Adduct": "DNA Adducts (Lesion Types)",
    "Pathway": "Pathways (KEGG References)",
    "Tissue": "Tissues (Expression Context)",
}
NODE_SIZES = {
    "Carcinogen": 12,
    "Enzyme": 9,
    "Receptor": 8,
    "Gene": 8,
    "Metabolite": 7,
    "DNA_Adduct": 8,
    "Pathway": 10,
    "Tissue": 9,
}
EDGE_COLORS = {
    # Legacy edge types (retained for backward compatibility; rarely/never
    # emitted into graph-data.json anymore -- see EdgeType in models.py)
    "ACTIVATES": "#e05565",
    "DETOXIFIES": "#6daa45",
    "TRANSPORTS": "#5591c7",
    "REPAIRS": "#e8af34",
    "EXPRESSED_IN": "#c2855a",
    "ENCODES": "#3d8b8b",
    "CUSTOM": "#9ea9bd",
    # Directional-schema edge types actually emitted into graph-data.json
    "FORMS_ADDUCT": "#a86fdf",
    "PRODUCES": "#d9704a",
    "DETOXIFIED_BY": "#6daa45",
    "TRANSPORTED_BY": "#5591c7",
    "REPAIRED_BY": "#e8af34",
    "SUBSTRATE_OF": "#9b8b3f",
    "AGONIZES": "#d1638f",
    "INDUCES": "#4fa382",
    "INHIBITS": "#a84a3c",
    "TRANSFORMS_SPONTANEOUSLY": "#7d8ba1",
    "SCHEMATIC_LINKAGE": "#a3a3a3",
    "PATHWAY": "#555b6a",
}
EDGE_LABELS = {
    "ACTIVATES": "Activates (procarcinogen → reactive metabolite)",
    "DETOXIFIES": "Detoxifies (conjugation / inactivation)",
    "TRANSPORTS": "Transports (efflux of conjugates)",
    "REPAIRS": "Repairs (DNA damage removal)",
    "EXPRESSED_IN": "Expressed in (gene → tissue)",
    "ENCODES": "Encodes (variant locus → gene or enzyme)",
    "CUSTOM": "Custom curated relation",
    "FORMS_ADDUCT": "Forms Adduct (metabolite → DNA lesion)",
    "PRODUCES": "Produces (enzyme → downstream metabolite/carcinogen)",
    "DETOXIFIED_BY": "Detoxified by (compound → detoxifying enzyme)",
    "TRANSPORTED_BY": "Transported by (compound → transporter enzyme)",
    "REPAIRED_BY": "Repaired by (DNA lesion → repair enzyme)",
    "SUBSTRATE_OF": "Substrate of (compound → metabolizing enzyme)",
    "AGONIZES": "Agonizes (ligand → receptor)",
    "INDUCES": "Induces (upregulates enzyme/receptor)",
    "INHIBITS": "Inhibits (downregulates enzyme/receptor)",
    "TRANSFORMS_SPONTANEOUSLY": "Transforms spontaneously (non-enzymatic)",
    "SCHEMATIC_LINKAGE": "Mechanism consolidated/summarized (dashed)",
    "PATHWAY": "Pathway membership",
}

# ── Search field tuples ───────────────────────────────────────────────────

NODE_SEARCH_FIELDS = (
    "id",
    "label",
    "type",
    "detail",
    "group",
    "iarc",
    "phase",
    "role",
    "source_db",
    "evidence",
    "pmid",
    "tissue",
    "variant",
    "phenotype",
    "provenance",
    "curation",
)
EDGE_SEARCH_FIELDS = (
    "source",
    "target",
    "type",
    "label",
    "carcinogen",
    "source_db",
    "evidence",
    "pmid",
    "tissue",
    "provenance",
    "curation",
)

# ── Session state bootstrap ──────────────────────────────────────────────


def get_engine() -> GraphEngine:
    if "engine" not in st.session_state:
        st.session_state.engine = GraphEngine()
    engine: GraphEngine = st.session_state.engine
    return engine


def start_engine(engine: GraphEngine) -> None:
    """Populate *engine* from the bundled reference graph, once per session.

    Uses :meth:`GraphEngine.load_reference_graph`, which loads
    ``map/graph-data.json`` (the canonical knowledge graph) and layers the
    tissue-expression and interaction-parameter overlays on top of it -- this
    is the single source of truth for every tab, not just the Reference Map.
    """
    if "graph_initialized" not in st.session_state:
        try:
            engine.load_reference_graph()
            st.session_state.graph_initialized = True
        except Exception as e:
            st.warning(f"Could not pre-load reference map: {e}")


# ── Pure helpers ──────────────────────────────────────────────────────────


def matches_query(data: dict[str, Any], query: str, fields: tuple[str, ...]) -> bool:
    if not query:
        return True

    def _flatten(value: object) -> list[str]:
        if value is None:
            return []
        if isinstance(value, dict):
            result: list[str] = []
            for nested in value.values():
                result.extend(_flatten(nested))
            return result
        if isinstance(value, list):
            result_list: list[str] = []
            for nested in value:
                result_list.extend(_flatten(nested))
            return result_list
        return [str(value)]

    values: list[str] = []
    for field in fields:
        values.extend(_flatten(data.get(field)))
    haystack = " ".join(values)
    return query.lower() in haystack.lower()


def build_provenance_record(
    *,
    source_db: str,
    record_id: str,
    evidence: str,
    pmid: str,
    tissue: str,
    citation: str,
    url: str,
) -> list[ProvenanceRecord]:
    record = ProvenanceRecord(
        source_db=source_db or None,
        record_id=record_id or None,
        evidence=evidence or None,
        pmid=pmid or None,
        tissue=tissue or None,
        citation=citation or None,
        url=url or None,
    )
    if not record.model_dump(exclude_none=True):
        return []
    return [record]


def build_curation_record(
    *,
    status: str,
    confidence: str,
    curator: str,
    reviewed_by: str,
    reviewed_at: str,
    notes: str,
) -> CurationRecord | None:
    if not any([status, confidence, curator, reviewed_by, reviewed_at, notes]):
        return None
    return CurationRecord(
        status=CurationStatus(status) if status else CurationStatus.DRAFT,
        confidence=CurationConfidence(confidence) if confidence else None,
        curator=curator or None,
        reviewed_by=reviewed_by or None,
        reviewed_at=reviewed_at or None,
        notes=notes or None,
    )
