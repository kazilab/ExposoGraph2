"""Build ``biomarker_mapping.json`` from split YAML registry files."""

from __future__ import annotations

import argparse
import copy
import json
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .loader import load_json_mapping, load_registry_document, write_json_mapping

_PRIVATE_SOURCE_KEY = "_source_document"


def _repo_root_from_source(source_path: Path) -> Path:
    """Infer the repository root from the documented scaffold path."""
    marker = Path("ExposoGraph") / "_biomarker_scaffold" / "data" / "registries"
    parts = source_path.resolve().parts
    marker_parts = marker.parts
    for index in range(0, len(parts) - len(marker_parts) + 1):
        if parts[index : index + len(marker_parts)] == marker_parts:
            return Path(*parts[:index])
    return Path.cwd()


def _resolve_registry_path(repo_root: Path, path_value: str) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    return repo_root / path


def _normalise_source_documents(
    manifest_path: Path, manifest: dict[str, Any]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    metadata = dict(manifest.get("_metadata") or manifest.get("metadata") or {})
    raw_sources = (
        manifest.get("source_documents")
        or manifest.get("sources")
        or metadata.get("source_documents")
        or []
    )
    if not isinstance(raw_sources, list):
        raise ValueError("Biomarker manifest source_documents must be a list")

    default_registry = "ExposoGraph/_biomarker_scaffold/data/registries/biomarkers_master.yaml"
    metadata.setdefault("source_registry_kind", "manifest")
    metadata.setdefault("source_registry", default_registry)
    metadata.setdefault("manifest_type", "biomarker_source_manifest")

    sources: list[dict[str, Any]] = []
    for raw_source in raw_sources:
        if not isinstance(raw_source, dict):
            raise ValueError(f"Invalid source document entry in {manifest_path}: {raw_source!r}")
        source = dict(raw_source)
        source_path = source.get("path")
        if not isinstance(source_path, str) or not source_path:
            raise ValueError(f"Source document is missing a path: {raw_source!r}")
        sources.append(source)
    return metadata, sources


def _load_manifest_entries(
    source_path: str | Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    manifest_path = Path(source_path)
    manifest = load_registry_document(manifest_path)
    repo_root = _repo_root_from_source(manifest_path)
    metadata, sources = _normalise_source_documents(manifest_path, manifest)

    entries: list[dict[str, Any]] = []
    normalized_sources: list[dict[str, Any]] = []
    for source in sources:
        registry_path = _resolve_registry_path(repo_root, str(source["path"]))
        source_document = load_registry_document(registry_path)
        source_metadata = dict(source_document.get("_metadata") or {})
        source_entries = source_document.get("entries") or []
        if not isinstance(source_entries, list):
            raise ValueError(f"Registry entries must be a list in {registry_path}")

        normalized_source = {**source_metadata, **source}
        normalized_source["entry_count"] = len(source_entries)
        normalized_sources.append(normalized_source)

        for raw_entry in source_entries:
            if not isinstance(raw_entry, dict):
                raise ValueError(f"Invalid biomarker row in {registry_path}: {raw_entry!r}")
            entry = copy.deepcopy(raw_entry)
            entry[_PRIVATE_SOURCE_KEY] = normalized_source
            entries.append(entry)

    metadata["source_documents"] = normalized_sources
    metadata["source_registry_entry_count"] = len(entries)
    return metadata, entries


def _entry_id(entry: dict[str, Any]) -> str:
    existing = entry.get("entry_id")
    if isinstance(existing, str) and existing:
        return existing
    biomarker = str(entry.get("biomarker", "")).strip()
    lifestyle_factor = str(entry.get("lifestyle_factor", "")).strip()
    if not biomarker or not lifestyle_factor:
        raise ValueError(f"Biomarker entry is missing key fields: {entry!r}")
    return f"{biomarker}::{lifestyle_factor}"


def _created_at(entry: dict[str, Any], fallback: str) -> str:
    trace = entry.get("trace")
    if isinstance(trace, dict):
        created_at = trace.get("created_at")
        if isinstance(created_at, str) and created_at:
            return created_at
    return fallback


def _public_entry(entry: dict[str, Any], index: int, generated_at: str) -> dict[str, Any]:
    source_document = entry.get(_PRIVATE_SOURCE_KEY)
    if not isinstance(source_document, dict):
        raise ValueError(f"Biomarker entry is missing source document context: {entry!r}")

    row = {
        key: copy.deepcopy(value)
        for key, value in entry.items()
        if key not in {_PRIVATE_SOURCE_KEY, "trace"}
    }
    row["entry_id"] = _entry_id(row)
    row["trace"] = {
        "created_index": index,
        "created_at": _created_at(entry, generated_at),
        "source_registry": source_document["path"],
        "source_family": source_document.get("source_family", ""),
        "registry_phase": source_document.get("registry_phase", ""),
        "registry_tier": source_document.get("registry_tier", ""),
        "source_note": source_document.get("source_note", ""),
    }
    return row


def load_biomarker_mapping_manifest(
    source_path: str | Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load the source manifest and return normalized metadata plus entries."""
    metadata, sourced_entries = _load_manifest_entries(source_path)
    entries = [
        {key: copy.deepcopy(value) for key, value in entry.items() if key != _PRIVATE_SOURCE_KEY}
        for entry in sourced_entries
    ]
    return metadata, entries


def build_biomarker_mapping_document(source_path: str | Path) -> dict[str, Any]:
    """Build the flattened JSON biomarker mapping document from a manifest."""
    metadata, sourced_entries = _load_manifest_entries(source_path)
    generated_at = datetime.now(timezone.utc).isoformat()
    metadata = copy.deepcopy(metadata)
    metadata.setdefault("schema_version", "1.2.0")
    metadata["generated_at"] = generated_at
    metadata.setdefault("forward_update_compatible", True)
    metadata.setdefault("customizable", True)
    metadata.setdefault(
        "description",
        "Manifest for rebuilding biomarker_mapping.json from split source YAML files.",
    )
    metadata.setdefault(
        "context",
        "v6 Methods / Exposure integration module: 'The complete mapping of Tier 2 "
        "lifestyle multipliers to their primary biomarker data sources, published "
        "concentration ranges, and derived [S]/Km ratios.'",
    )
    metadata.setdefault("version", "1.0.0")
    metadata.setdefault("created", "2026-04-20")
    metadata.setdefault("registry_phase", "split_manifest")
    metadata.setdefault("registry_tier", "mixed")
    metadata.setdefault(
        "source_note",
        "Master manifest for the split biomarker registry. Edit the owning source "
        "file, then rebuild and compare.",
    )
    entries = [
        _public_entry(entry, index=index, generated_at=generated_at)
        for index, entry in enumerate(sourced_entries)
    ]
    return {
        "_metadata": metadata,
        "_update_list": [],
        "entries": entries,
    }


def _compare_key(entry: dict[str, Any]) -> str:
    return _entry_id(entry)


def _compare_entry(entry: dict[str, Any]) -> dict[str, Any]:
    return {
        key: copy.deepcopy(value)
        for key, value in entry.items()
        if key != "trace"
    }


def _entries_by_id(entries: Iterable[Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for raw_entry in entries:
        if not isinstance(raw_entry, dict):
            raise ValueError(f"Invalid biomarker entry: {raw_entry!r}")
        result[_compare_key(raw_entry)] = _compare_entry(raw_entry)
    return result


def compare_biomarker_mapping_documents(
    old_document: dict[str, Any], new_document: dict[str, Any]
) -> dict[str, Any]:
    """Compare biomarker rows while ignoring scaffold trace metadata."""
    old_entries = _entries_by_id(old_document.get("entries") or [])
    new_entries = _entries_by_id(new_document.get("entries") or [])

    old_ids = set(old_entries)
    new_ids = set(new_entries)
    added = sorted(new_ids - old_ids)
    removed = sorted(old_ids - new_ids)
    changed = sorted(
        entry_id
        for entry_id in old_ids & new_ids
        if old_entries[entry_id] != new_entries[entry_id]
    )

    return {
        "mapped_biomarkers_unchanged": not added and not removed and not changed,
        "old_count": len(old_entries),
        "new_count": len(new_entries),
        "added_count": len(added),
        "removed_count": len(removed),
        "changed_count": len(changed),
        "added": added,
        "removed": removed,
        "changed": changed,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build ExposoGraph biomarker_mapping.json from split YAML registries."
    )
    parser.add_argument(
        "--source",
        default="ExposoGraph/_biomarker_scaffold/data/registries/biomarkers_master.yaml",
        help="Path to the biomarker source manifest.",
    )
    parser.add_argument(
        "--out",
        default="ExposoGraph/data/biomarker_mapping.json",
        help="Output JSON mapping path.",
    )
    parser.add_argument(
        "--old",
        default=None,
        help="Optional previous JSON mapping path to compare against.",
    )
    parser.add_argument(
        "--compare-only",
        action="store_true",
        help="Compare --old and --out without rebuilding --out.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.compare_only:
        if not args.old:
            raise SystemExit("--compare-only requires --old")
        old_document = load_json_mapping(args.old)
        new_document = load_json_mapping(args.out)
        report = compare_biomarker_mapping_documents(old_document, new_document)
        print(json.dumps(report, indent=2, sort_keys=True))
        raise SystemExit(0 if report["mapped_biomarkers_unchanged"] else 1)

    document = build_biomarker_mapping_document(args.source)
    write_json_mapping(args.out, document)
    if args.old:
        old_document = load_json_mapping(args.old)
        report = compare_biomarker_mapping_documents(old_document, document)
        print(json.dumps(report, indent=2, sort_keys=True))
        raise SystemExit(0 if report["mapped_biomarkers_unchanged"] else 1)
    print(f"Wrote {len(document['entries'])} biomarker rows to {args.out}")


if __name__ == "__main__":
    main()
