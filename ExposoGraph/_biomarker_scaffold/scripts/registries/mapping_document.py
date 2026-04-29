"""Helpers for traceable, update-compatible biomarker mapping documents."""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, UTC
from typing import Any

REQUIRED_ENTRY_FIELDS = ("biomarker", "matrix", "reference_range", "reference_units", "source_status")
REQUIRED_METADATA_FIELDS = ("schema_version", "generated_at", "forward_update_compatible", "customizable")
REQUIRED_TRACE_FIELDS = ("created_at", "created_index")


def _entry_identity(entry: dict[str, Any]) -> str:
    biomarker = str(entry.get("biomarker") or entry.get("biomarker_id") or "unknown")
    lifestyle = str(entry.get("lifestyle_factor") or entry.get("lifestyle") or "default")
    return str(entry.get("entry_id") or f"{biomarker}::{lifestyle}")


def build_mapping_document(
    entries: list[dict[str, Any]],
    metadata: dict[str, Any] | None = None,
    update_list: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Build a validation-ready mapping document with traceability fields."""
    normalized_entries: list[dict[str, Any]] = []
    for idx, raw in enumerate(entries):
        entry = deepcopy(raw)
        entry.setdefault("entry_id", _entry_identity(entry))
        entry.setdefault("trace", {})
        entry["trace"].setdefault("created_index", idx)
        entry["trace"].setdefault("created_at", datetime.now(UTC).isoformat())
        normalized_entries.append(entry)

    out_metadata = {
        "schema_version": "1.1.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "forward_update_compatible": True,
        "customizable": True,
    }
    if metadata:
        out_metadata.update(metadata)

    return {
        "_metadata": out_metadata,
        "_update_list": list(update_list or []),
        "entries": normalized_entries,
    }


def normalize_mapping_document(
    document: dict[str, Any],
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a normalized mapping document preserving known top-level fields."""
    base_metadata = dict(document.get("_metadata", {})) if isinstance(document.get("_metadata"), dict) else {}
    if metadata:
        base_metadata.update(metadata)
    updates = document.get("_update_list", [])
    if not isinstance(updates, list):
        updates = []
    entries = document.get("entries", [])
    if not isinstance(entries, list):
        entries = []
    normalized = build_mapping_document(entries=entries, metadata=base_metadata, update_list=updates)
    for key, value in document.items():
        if key in {"_metadata", "_update_list", "entries"}:
            continue
        normalized[key] = deepcopy(value)
    return normalized


def apply_update_list(document: dict[str, Any], updates: list[dict[str, Any]]) -> dict[str, Any]:
    """Apply forward-compatible updates (upsert/remove) by entry_id."""
    out = deepcopy(document)
    entries = [dict(item) for item in out.get("entries", []) if isinstance(item, dict)]
    index = {str(item.get("entry_id") or _entry_identity(item)): i for i, item in enumerate(entries)}

    for update in updates:
        op = update.get("op", "upsert")
        target = str(update.get("entry_id") or "")
        payload = dict(update.get("entry", {}))
        if not target and payload:
            target = _entry_identity(payload)
        if not target:
            continue

        if op == "remove":
            if target in index:
                entries.pop(index[target])
                index = {str(item.get("entry_id") or _entry_identity(item)): i for i, item in enumerate(entries)}
            continue

        payload.setdefault("entry_id", target)
        payload.setdefault("trace", {})
        payload["trace"].setdefault("updated_at", datetime.now(UTC).isoformat())
        if target in index:
            merged = dict(entries[index[target]])
            merged.update(payload)
            entries[index[target]] = merged
        else:
            entries.append(payload)
            index[target] = len(entries) - 1

    out["entries"] = entries
    out.setdefault("_update_list", [])
    out["_update_list"].extend(updates)
    return out


def validate_mapping_document(document: dict[str, Any]) -> list[str]:
    """Return validation issues for mapping document structure/content."""
    issues: list[str] = []
    metadata = document.get("_metadata")
    if not isinstance(metadata, dict):
        issues.append("Document must include an '_metadata' object.")
    else:
        for field in REQUIRED_METADATA_FIELDS:
            if metadata.get(field) in (None, ""):
                issues.append(f"Document metadata missing required field: {field}")
    update_list = document.get("_update_list")
    if not isinstance(update_list, list):
        issues.append("Document must include an '_update_list' list.")

    entries = document.get("entries")
    if not isinstance(entries, list):
        return ["Document must include an 'entries' list."]

    seen: set[tuple[str, str]] = set()
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            issues.append(f"entries[{i}] is not an object.")
            continue
        missing = [field for field in REQUIRED_ENTRY_FIELDS if entry.get(field) in (None, "")]
        if missing:
            issues.append(f"entries[{i}] missing required fields: {', '.join(missing)}")
        trace = entry.get("trace")
        if not isinstance(trace, dict):
            issues.append(f"entries[{i}] missing trace object")
        else:
            missing_trace = [field for field in REQUIRED_TRACE_FIELDS if trace.get(field) in (None, "")]
            if missing_trace:
                issues.append(f"entries[{i}] trace missing fields: {', '.join(missing_trace)}")
        if entry.get("entry_id") in (None, ""):
            issues.append(f"entries[{i}] missing required fields: entry_id")
        pair = (
            str(entry.get("biomarker", "")),
            str(entry.get("lifestyle_factor") or entry.get("lifestyle") or "default"),
        )
        if not entry.get("allow_duplicate_combo"):
            if pair in seen:
                issues.append(f"Duplicate biomarker/lifestyle combination not allowed: {pair[0]}::{pair[1]}")
            seen.add(pair)
    return issues

