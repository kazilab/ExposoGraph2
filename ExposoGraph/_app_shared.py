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

# ── Color maps ────────────────────────────────────────────────────────────

NODE_COLORS = {
    "Carcinogen": "#e05565",
    "Enzyme": "#4f98a3",
    "Gene": "#3d8b8b",
    "Metabolite": "#e8945a",
    "DNA_Adduct": "#a86fdf",
    "Pathway": "#5591c7",
    "Tissue": "#c2855a",
}
EDGE_COLORS = {
    "ACTIVATES": "#e05565",
    "DETOXIFIES": "#6daa45",
    "TRANSPORTS": "#5591c7",
    "FORMS_ADDUCT": "#a86fdf",
    "REPAIRS": "#e8af34",
    "PATHWAY": "#888888",
    "EXPRESSED_IN": "#c2855a",
    "INDUCES": "#d4a843",
    "INHIBITS": "#8b4a6b",
    "ENCODES": "#3d8b8b",
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
