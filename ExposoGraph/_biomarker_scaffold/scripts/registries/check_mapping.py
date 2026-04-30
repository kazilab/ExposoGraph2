"""Validate and optionally normalize ``biomarker_mapping.json``."""

from __future__ import annotations

import argparse
import copy
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .build_mapping import _entry_id
from .loader import load_json_mapping, write_json_mapping


def validate_biomarker_mapping_document(
    document: dict[str, Any], *, fix: bool = False
) -> tuple[dict[str, Any], list[str]]:
    """Validate a biomarker mapping document.

    When ``fix`` is true, missing ``entry_id`` and minimal trace fields are
    populated in a copied document. The returned error list is empty only when
    the original or fixed document satisfies the scaffold contract.
    """
    normalized = copy.deepcopy(document)
    errors: list[str] = []

    metadata = normalized.setdefault("_metadata", {})
    if not isinstance(metadata, dict):
        errors.append("_metadata must be an object")
        return normalized, errors
    metadata.setdefault("forward_update_compatible", True)
    metadata.setdefault("customizable", True)

    update_list = normalized.setdefault("_update_list", [])
    if not isinstance(update_list, list):
        errors.append("_update_list must be a list")

    entries = normalized.get("entries")
    if not isinstance(entries, list):
        errors.append("entries must be a list")
        return normalized, errors

    generated_at = str(
        metadata.get("generated_at") or datetime.now(timezone.utc).isoformat()
    )
    seen: set[str] = set()
    for index, raw_entry in enumerate(entries):
        if not isinstance(raw_entry, dict):
            errors.append(f"entries[{index}] must be an object")
            continue
        try:
            entry_id = _entry_id(raw_entry)
        except ValueError as exc:
            errors.append(str(exc))
            continue

        if entry_id in seen:
            errors.append(f"duplicate biomarker entry_id: {entry_id}")
        seen.add(entry_id)

        if fix:
            raw_entry.setdefault("entry_id", entry_id)
            trace = raw_entry.setdefault("trace", {})
            if isinstance(trace, dict):
                trace.setdefault("created_index", index)
                trace.setdefault("created_at", generated_at)
                trace.setdefault("source_registry", "")
        elif raw_entry.get("entry_id") != entry_id:
            errors.append(f"{entry_id} has a missing or inconsistent entry_id")

        trace_value = raw_entry.get("trace")
        if not isinstance(trace_value, dict):
            errors.append(f"{entry_id} trace must be an object")
        elif not trace_value.get("source_registry"):
            errors.append(f"{entry_id} trace.source_registry is required")

    return normalized, errors


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate ExposoGraph biomarker_mapping.json trace/update fields."
    )
    parser.add_argument(
        "--mapping",
        default="ExposoGraph/data/biomarker_mapping.json",
        help="Path to biomarker_mapping.json.",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Normalize missing forward-compatible fields in memory.",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write normalized content back to --mapping. Requires --fix.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.write and not args.fix:
        raise SystemExit("--write requires --fix")

    mapping_path = Path(args.mapping)
    document = load_json_mapping(mapping_path)
    normalized, errors = validate_biomarker_mapping_document(document, fix=args.fix)
    if args.write and not errors:
        write_json_mapping(mapping_path, normalized)

    summary = {
        "mapping": str(mapping_path),
        "entry_count": len(normalized.get("entries") or []),
        "valid": not errors,
        "errors": errors,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    raise SystemExit(0 if not errors else 1)


if __name__ == "__main__":
    main()
