"""Shared constants, helpers, and session-state bootstrap for the Streamlit UI."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any

import streamlit as st

from .config import GraphVisibility, get_app_mode, persistence_enabled
from .engine import GraphEngine
from .exporter import parse_graph_data_text
from .graph_filters import graph_visibility_label
from .llm_extractor import EXAMPLE_INPUT
from .models import (
    CurationConfidence,
    CurationRecord,
    CurationStatus,
    KnowledgeGraph,
    ProvenanceRecord,
)
from .storage import GraphRepository, GraphRevisionSummary

# ── Path constants ────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
REPOSITORY_PATH = DATA_DIR / "ExposoGraph.sqlite3"
PROJECTS_DIR = PROJECT_ROOT / "saved_graphs"
APP_MODE = get_app_mode()
PERSISTENCE_ENABLED = persistence_enabled(APP_MODE)
if PERSISTENCE_ENABLED:
    DATA_DIR.mkdir(exist_ok=True)
    PROJECTS_DIR.mkdir(exist_ok=True)

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
    if "graph_initialized" not in st.session_state:
        try:
            default_js_path = (
                Path(__file__).resolve().parent / "map" / "graph_data_final.json"
            )
            if default_js_path.exists():
                engine.load_from_json(default_js_path)
                st.session_state.graph_initialized = True
        except Exception as e:
            st.warning(f"Could not pre-load reference map: {e}")


@st.cache_resource
def get_repository() -> GraphRepository | None:
    if not PERSISTENCE_ENABLED:
        return None
    return GraphRepository(REPOSITORY_PATH)


def get_pending_extraction() -> KnowledgeGraph | None:
    raw = st.session_state.get("pending_extraction")
    if raw is None:
        return None
    return KnowledgeGraph(**raw)


def load_example_text() -> None:
    st.session_state["extract_text"] = EXAMPLE_INPUT


# ── Pure helpers ──────────────────────────────────────────────────────────


def get_secret(key: str, default: str = "") -> str:
    """Read a secret from Streamlit secrets (Cloud) or env vars (local)."""
    try:
        return str(st.secrets[key])
    except (KeyError, FileNotFoundError):
        return os.environ.get(key, default)


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


def relative_path(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def load_into_engine(
    engine: GraphEngine, kg: KnowledgeGraph, *, replace: bool
) -> list[str]:
    if replace:
        return engine.load(kg)
    return engine.merge(kg)


def slugify_project_name(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    slug = slug.strip("._")
    return slug or "knowledge_graph"


def saved_project_paths() -> list[Path]:
    return sorted(list(PROJECTS_DIR.glob("*.html")) + list(PROJECTS_DIR.glob("*.json")))


def revision_label(
    revision_id: int, revisions_by_id: dict[int, GraphRevisionSummary]
) -> str:
    revision = revisions_by_id[revision_id]
    note = f" - {revision.note}" if revision.note else ""
    visibility = (
        f" | {graph_visibility_label(revision.visibility)}"
        if revision.visibility != GraphVisibility.ALL
        else ""
    )
    return (
        f"r{revision.revision_number} | {revision.node_count}n/{revision.edge_count}e | "
        f"{revision.created_at[:19]}{visibility}{note}"
    )


def parse_uploaded_graph(name: str, raw_text: str) -> KnowledgeGraph:
    suffix = Path(name).suffix.lower()
    if suffix == ".json":
        return KnowledgeGraph(**json.loads(raw_text))
    if suffix in {".html", ".js"}:
        return parse_graph_data_text(raw_text)
    raise ValueError(f"Unsupported upload type: {suffix}")


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


def annotation_lines(data: dict[str, Any]) -> list[str]:
    lines: list[str] = []
    for label, key in [
        ("Source", "source_db"),
        ("Evidence", "evidence"),
        ("PMID", "pmid"),
        ("Tissue", "tissue"),
        ("Variant", "variant"),
        ("Phenotype", "phenotype"),
    ]:
        value = data.get(key)
        if value:
            lines.append(f"{label}: {value}")
    provenance = data.get("provenance") or []
    if provenance:
        record_ids = []
        for record in provenance:
            record_id = record.get("record_id")
            if record_id and record_id not in record_ids:
                record_ids.append(record_id)
        if record_ids:
            lines.append(f"Record IDs: {'; '.join(record_ids)}")
        lines.append(f"Provenance records: {len(provenance)}")
    curation = data.get("curation") or {}
    for label, key in [
        ("Status", "status"),
        ("Confidence", "confidence"),
        ("Curator", "curator"),
        ("Reviewed by", "reviewed_by"),
    ]:
        value = curation.get(key)
        if value:
            lines.append(f"{label}: {value}")
    if data.get("activity_score") is not None:
        lines.append(f"Activity score: {data['activity_score']}")
    return lines


def node_tooltip(data: dict[str, Any]) -> str:
    parts = [data.get("detail", "")]
    parts.extend(annotation_lines(data))
    return "\n".join(part for part in parts if part)


def edge_tooltip(data: dict[str, Any]) -> str:
    parts = [data.get("label", data.get("type", ""))]
    parts.extend(annotation_lines(data))
    return "\n".join(part for part in parts if part)
