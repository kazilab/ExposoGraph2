"""Build and compare biomarker_mapping.json documents from split YAML sources."""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from .loader import load_json_mapping
from .mapping_document import build_mapping_document, validate_mapping_document

DEFAULT_MANIFEST_PATH = (
    Path(__file__).resolve().parents[3]
    / "_biomarker_scaffold"
    / "data"
    / "registries"
    / "biomarkers_master.yaml"
)
DEFAULT_SOURCE_PATH = DEFAULT_MANIFEST_PATH
DEFAULT_OUTPUT_PATH = Path(__file__).resolve().parents[3] / "data" / "biomarker_mapping.json"
DEFAULT_OLD_PATH = Path(__file__).resolve().parents[3] / "data" / "biomarker_mapping_old.json"
REPO_ROOT = Path(__file__).resolve().parents[4]

_IGNORED_ENTRY_KEYS = {"entry_id", "trace", "biomarker_id"}
_SOURCE_DOCUMENT_KEYS = {"biomarkers", "entries"}
_MANIFEST_KEYS = {"sources"}


def _coerce_source_entry(raw_entry: dict[str, Any]) -> dict[str, Any]:
    entry = deepcopy(raw_entry)
    biomarker = entry.get("biomarker") or entry.get("biomarker_id")
    if not biomarker:
        raise ValueError("Each source entry must define biomarker or biomarker_id")
    entry.setdefault("biomarker", biomarker)
    entry.pop("biomarker_id", None)
    if "lifestyle" in entry and "lifestyle_factor" not in entry:
        entry["lifestyle_factor"] = entry.pop("lifestyle")
    entry.setdefault("lifestyle_factor", "default")
    if "references" not in entry or entry["references"] is None:
        entry["references"] = []
    return entry


def _load_yaml_document(path: str | Path) -> Any:
    path = Path(path)
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _registry_path_text(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(resolved)


def _split_source_document(
    path: Path,
    raw: Any,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    metadata: dict[str, Any] = {}
    entries: Any = raw
    if isinstance(raw, dict):
        if "sources" in raw:
            raise ValueError(f"Manifest files must be loaded through the manifest loader: {path}")
        metadata = {k: v for k, v in raw.items() if k not in _SOURCE_DOCUMENT_KEYS}
        if "biomarkers" in raw:
            entries = raw["biomarkers"]
        elif "entries" in raw:
            entries = raw["entries"]
        else:
            entries = []
    if entries is None:
        entries = []
    if not isinstance(entries, list):
        raise ValueError(f"Source YAML must contain a list of biomarker records: {path}")

    out: list[dict[str, Any]] = []
    for idx, item in enumerate(entries):
        if not isinstance(item, dict):
            raise ValueError(f"Source entry {idx} in {path} is not a mapping")
        out.append(item)
    return metadata, out


def _resolve_source_spec(
    manifest_path: Path,
    spec: str | dict[str, Any],
) -> tuple[Path, dict[str, Any]]:
    if isinstance(spec, str):
        source_path = Path(spec)
        source_meta: dict[str, Any] = {}
    elif isinstance(spec, dict):
        source_path = Path(spec.get("path") or "")
        if not source_path:
            raise ValueError("Each manifest source entry must include a path")
        source_meta = {k: v for k, v in spec.items() if k != "path"}
    else:
        raise ValueError("Manifest source entries must be strings or objects")

    if not source_path.is_absolute():
        source_path = manifest_path.parent / source_path
    return source_path, source_meta


def _trace_metadata(source_path: Path, source_meta: dict[str, Any]) -> dict[str, Any]:
    trace: dict[str, Any] = {"source_registry": _registry_path_text(source_path)}
    for key in ("source_family", "registry_phase", "registry_tier", "source_note"):
        value = source_meta.get(key)
        if value not in (None, ""):
            trace[key] = value
    return trace


def _prepare_entries(
    entries: list[dict[str, Any]],
    source_path: Path,
    source_meta: dict[str, Any],
) -> list[dict[str, Any]]:
    trace_template = _trace_metadata(source_path, source_meta)
    prepared_entries: list[dict[str, Any]] = []
    for entry in entries:
        prepared = _coerce_source_entry(entry)
        trace = dict(prepared.get("trace", {}))
        for key, value in trace_template.items():
            trace.setdefault(key, value)
        prepared["trace"] = trace
        prepared_entries.append(prepared)
    return prepared_entries


def _summarize_source(
    source_path: Path,
    source_meta: dict[str, Any],
    entries: list[dict[str, Any]],
) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "path": _registry_path_text(source_path),
        "entry_count": len(entries),
    }
    for key in ("source_family", "registry_phase", "registry_tier", "source_note"):
        value = source_meta.get(key)
        if value not in (None, ""):
            summary[key] = value
    return summary


def _load_source_file(
    path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    raw = _load_yaml_document(path)
    metadata, entries = _split_source_document(path, raw)
    return metadata, entries


def _load_single_source_bundle(
    path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    source_metadata, entries = _load_source_file(path)
    prepared_entries = _prepare_entries(entries, path, source_metadata)
    metadata = dict(source_metadata)
    metadata["source_registry_kind"] = "single"
    metadata["source_registry"] = _registry_path_text(path)
    metadata["source_registry_entry_count"] = len(prepared_entries)
    metadata["source_documents"] = [
        _summarize_source(path, source_metadata, prepared_entries),
    ]
    return metadata, prepared_entries


def _load_manifest_bundle(
    path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    raw = _load_yaml_document(path)
    if not isinstance(raw, dict) or "sources" not in raw:
        return _load_single_source_bundle(path)

    manifest_metadata = {k: v for k, v in raw.items() if k not in _MANIFEST_KEYS}
    source_documents: list[dict[str, Any]] = []
    entries: list[dict[str, Any]] = []
    for spec in raw["sources"]:
        source_path, source_meta = _resolve_source_spec(path, spec)
        source_metadata, source_entries = _load_source_file(source_path)
        combined_meta = dict(source_metadata)
        combined_meta.update(source_meta)
        prepared_entries = _prepare_entries(source_entries, source_path, combined_meta)
        entries.extend(prepared_entries)
        source_documents.append(_summarize_source(source_path, combined_meta, prepared_entries))

    manifest_metadata["source_registry_kind"] = "manifest"
    manifest_metadata["source_registry"] = _registry_path_text(path)
    manifest_metadata["source_registry_entry_count"] = len(entries)
    manifest_metadata["source_documents"] = source_documents
    return manifest_metadata, entries


def load_biomarker_mapping_source(
    path: str | Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load a biomarker mapping YAML source or manifest document."""
    path = Path(path)
    raw = _load_yaml_document(path)
    if isinstance(raw, dict) and "sources" in raw:
        return _load_manifest_bundle(path)
    return _load_single_source_bundle(path)


def load_biomarker_mapping_manifest(
    path: str | Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load a biomarker mapping manifest and all of its referenced sources."""
    path = Path(path)
    return _load_manifest_bundle(path)


def build_biomarker_mapping_document(
    source_path: str | Path = DEFAULT_SOURCE_PATH,
    *,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a normalized biomarker mapping document from a YAML source."""
    source_metadata, entries = load_biomarker_mapping_source(source_path)
    build_metadata = dict(source_metadata)
    if metadata:
        build_metadata.update(metadata)
    return build_mapping_document(entries=entries, metadata=build_metadata)


def _entry_key(entry: dict[str, Any]) -> str:
    biomarker = str(entry.get("biomarker") or entry.get("biomarker_id") or "").strip()
    lifestyle = str(entry.get("lifestyle_factor") or entry.get("lifestyle") or "default").strip()
    if biomarker:
        return f"{biomarker}::{lifestyle}"
    return str(entry.get("entry_id") or f"entry::{lifestyle}")


def _compareable_entry(entry: dict[str, Any]) -> dict[str, Any]:
    payload = deepcopy(entry)
    for key in _IGNORED_ENTRY_KEYS:
        payload.pop(key, None)
    return payload


def _index_entries(entries: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for idx, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"Entry {idx} is not a mapping")
        key = _entry_key(entry)
        if key in index:
            raise ValueError(f"Duplicate biomarker/lifestyle combination in mapping: {key}")
        index[key] = _compareable_entry(entry)
    return index


def compare_biomarker_mapping_documents(
    old_document: dict[str, Any],
    new_document: dict[str, Any],
) -> dict[str, Any]:
    """Compare two mapping documents using stable biomarker/lifestyle identities."""
    old_entries = old_document.get("entries", [])
    new_entries = new_document.get("entries", [])
    if not isinstance(old_entries, list):
        raise ValueError("Old document must contain an entries list")
    if not isinstance(new_entries, list):
        raise ValueError("New document must contain an entries list")

    old_index = _index_entries([entry for entry in old_entries if isinstance(entry, dict)])
    new_index = _index_entries([entry for entry in new_entries if isinstance(entry, dict)])

    old_keys = set(old_index)
    new_keys = set(new_index)
    shared_keys = sorted(old_keys & new_keys)
    added_keys = sorted(new_keys - old_keys)
    removed_keys = sorted(old_keys - new_keys)

    changed_entries: list[dict[str, Any]] = []
    for key in shared_keys:
        old_payload = old_index[key]
        new_payload = new_index[key]
        if old_payload == new_payload:
            continue
        changed_fields = sorted(
            field
            for field in set(old_payload) | set(new_payload)
            if old_payload.get(field) != new_payload.get(field)
        )
        changed_entries.append({"entry_id": key, "fields_changed": changed_fields})

    return {
        "old_entry_count": len(old_index),
        "new_entry_count": len(new_index),
        "shared_entry_count": len(shared_keys),
        "added_count": len(added_keys),
        "removed_count": len(removed_keys),
        "changed_count": len(changed_entries),
        "added_entries": added_keys,
        "removed_entries": removed_keys,
        "changed_entries": changed_entries,
        "mapped_biomarkers_unchanged": not added_keys and not removed_keys and not changed_entries,
    }


def _print_summary(report: dict[str, Any]) -> None:
    print(
        "Comparison summary: "
        f"{report['shared_entry_count']} shared, "
        f"{report['added_count']} added, "
        f"{report['removed_count']} removed, "
        f"{report['changed_count']} changed"
    )
    if report["mapped_biomarkers_unchanged"]:
        print("Mapped biomarkers unchanged between the two documents.")
    else:
        if report["added_entries"]:
            print("Added:")
            for key in report["added_entries"]:
                print(f"- {key}")
        if report["removed_entries"]:
            print("Removed:")
            for key in report["removed_entries"]:
                print(f"- {key}")
        if report["changed_entries"]:
            print("Changed:")
            for item in report["changed_entries"]:
                fields = ", ".join(item["fields_changed"])
                print(f"- {item['entry_id']}: {fields}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild biomarker_mapping.json from a split YAML manifest and compare it "
            "with the old snapshot."
        )
    )
    parser.add_argument(
        "--source",
        default=str(DEFAULT_SOURCE_PATH),
        help="YAML source manifest or source file path.",
    )
    parser.add_argument(
        "--out",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Output biomarker_mapping.json path.",
    )
    parser.add_argument(
        "--old",
        default=str(DEFAULT_OLD_PATH),
        help="Preserved old biomarker_mapping.json path used for comparison.",
    )
    parser.add_argument(
        "--compare-only",
        action="store_true",
        help="Skip rebuilding and compare --old and --out as existing JSON documents.",
    )
    parser.add_argument(
        "--report",
        help="Optional JSON path for the comparison report.",
    )
    args = parser.parse_args(argv)

    if args.compare_only:
        old_path = Path(args.old)
        new_path = Path(args.out)
        old_document = load_json_mapping(old_path)
        new_document = load_json_mapping(new_path)
        report = compare_biomarker_mapping_documents(old_document, new_document)
        _print_summary(report)
        if args.report:
            report_path = Path(args.report)
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_text = json.dumps(report, indent=2, sort_keys=True) + "\n"
            report_path.write_text(report_text, encoding="utf-8")
        return 0

    source_path = Path(args.source)
    output_path = Path(args.out)
    old_path = Path(args.old)

    document = build_biomarker_mapping_document(source_path)
    issues = validate_mapping_document(document)
    if issues:
        print(f"FAILED: {len(issues)} issue(s) found")
        for issue in issues:
            print(f"- {issue}")
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(document, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(f"Wrote rebuilt biomarker mapping to {output_path}")

    if old_path.exists():
        old_document = load_json_mapping(old_path)
        report = compare_biomarker_mapping_documents(old_document, document)
        _print_summary(report)
        if args.report:
            report_path = Path(args.report)
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_text = json.dumps(report, indent=2, sort_keys=True) + "\n"
            report_path.write_text(report_text, encoding="utf-8")
    else:
        print(f"Old snapshot not found at {old_path}; comparison skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
