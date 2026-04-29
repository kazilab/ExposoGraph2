"""Evidence classification and reporting utilities."""

from __future__ import annotations

from collections import Counter
from typing import Any, Iterable

MODEL_PROXY_KEYWORDS = ("equivalent", "index", "composite", "score")
SCENARIO_PROXY_IDS = {"blood_ethanol", "saliva_acetaldehyde"}


def classify_source_status(record: dict[str, Any]) -> str:
    """Assign a conservative coverage classification to a resolved record."""
    explicit = record.get("coverage_status") or record.get("source_status")
    if explicit in {"nhanes_native", "literature_biomarker", "scenario_proxy", "model_proxy", "derived_proxy", "unsupported"}:
        return explicit
    biomarker = str(record.get("biomarker") or record.get("biomarker_id") or "").lower()
    if biomarker in SCENARIO_PROXY_IDS:
        return "scenario_proxy"
    if any(keyword in biomarker for keyword in MODEL_PROXY_KEYWORDS):
        return "model_proxy"
    if record.get("nhanes_variables") or record.get("nhanes_variable"):
        return "nhanes_native"
    if record.get("references") or record.get("evidence"):
        return "literature_biomarker"
    return "unsupported"


def coverage_report(records: Iterable[dict[str, Any]]) -> dict[str, Any]:
    records = list(records)
    statuses = Counter(classify_source_status(record) for record in records)
    return {
        "total_mapping_entries": len(records),
        "entries_with_nhanes_evidence": statuses.get("nhanes_native", 0),
        "entries_with_literature_proxy_only": statuses.get("literature_biomarker", 0) + statuses.get("scenario_proxy", 0) + statuses.get("model_proxy", 0) + statuses.get("derived_proxy", 0),
        "entries_with_brenda_kinetic_support": sum(1 for record in records if record.get("Km_uM") is not None or record.get("target_enzyme")),
        "entries_with_pubchem_id": sum(1 for record in records if record.get("pubchem_cid")),
        "entries_with_comptox_id": sum(1 for record in records if record.get("comptox_dtxsid")),
        "entries_missing_Km": sum(1 for record in records if record.get("Km_uM") is None),
        "entries_missing_source_references": sum(1 for record in records if not record.get("references") and not record.get("evidence")),
        "entries_with_model_proxy_status": statuses.get("model_proxy", 0),
        "status_counts": dict(statuses),
    }
