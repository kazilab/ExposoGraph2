"""Resolve biomarker mapping records across measurement and external evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .evidence import classify_source_status, coverage_report
from .loader import load_json_mapping, load_yaml_registry

MATRIX_PREFIXES = (
    "urinary_",
    "blood_",
    "serum_",
    "plasma_",
    "whole_blood_",
    "saliva_",
    "exhaled_breath_",
)


def normalize_alias(name: str) -> str:
    value = name.strip().lower().replace(" ", "_").replace("-", "_")
    value = "_".join(part for part in value.split("_") if part)
    for prefix in MATRIX_PREFIXES:
        if value.startswith(prefix):
            return value[len(prefix):]
    return value


def _iter_mapping_entries(model_mapping: dict[str, Any]) -> list[dict[str, Any]]:
    if isinstance(model_mapping.get("entries"), list):
        entries: list[dict[str, Any]] = []
        for idx, raw_entry in enumerate(model_mapping["entries"]):
            if not isinstance(raw_entry, dict):
                continue
            entry = dict(raw_entry)
            biomarker = entry.get("biomarker") or entry.get("biomarker_id")
            lifestyle = entry.get("lifestyle_factor") or entry.get("lifestyle") or "default"
            entry.setdefault("biomarker", biomarker)
            if biomarker:
                entry.setdefault("biomarker_id", f"{biomarker}::{lifestyle}")
            else:
                entry.setdefault("biomarker_id", f"entry_{idx}")
            entry.setdefault(
                "trace",
                {
                    "source_format": "entries_list",
                    "source_index": idx,
                    "lifestyle_factor": lifestyle,
                },
            )
            entries.append(entry)
        return entries

    entries: list[dict[str, Any]] = []
    for key, value in model_mapping.items():
        if key.startswith("_"):
            continue
        if isinstance(value, dict):
            entry = dict(value)
            entry.setdefault("biomarker", entry.get("biomarker_id", key))
            entry.setdefault("biomarker_id", entry.get("biomarker", key))
            entry.setdefault("trace", {"source_format": "keyed_object", "mapping_key": key})
            entries.append(entry)
    return entries


def _index_measurement_registry(measurement_registry: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for item in measurement_registry:
        biomarker = item.get("biomarker") or item.get("biomarker_id")
        if biomarker:
            index[str(biomarker)] = item
            index[normalize_alias(str(biomarker))] = item
        variable = item.get("nhanes_variable") or item.get("variable")
        if variable:
            index[str(variable).upper()] = item
    return index


def _index_mapping_entries(model_mapping: dict[str, Any]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for entry in _iter_mapping_entries(model_mapping):
        biomarker_id = str(entry.get("biomarker_id") or entry.get("biomarker") or "")
        if not biomarker_id:
            continue
        index[biomarker_id] = entry
        index[normalize_alias(biomarker_id)] = entry
        biomarker = entry.get("biomarker")
        if biomarker:
            index[str(biomarker)] = entry
            index[normalize_alias(str(biomarker))] = entry
    return index


def _index_external_records(records: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for item in records:
        for key in ("biomarker_id", "biomarker", "canonical_name", "identifier"):
            value = item.get(key)
            if not value:
                continue
            index[str(value)] = item
            index[normalize_alias(str(value))] = item
    return index


def resolve_biomarker(
    biomarker_id: str,
    measurement_registry: list[dict[str, Any]],
    model_mapping: dict[str, Any],
    external_sources: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    measurement_index = _index_measurement_registry(measurement_registry)
    measurement = (
        measurement_index.get(biomarker_id)
        or measurement_index.get(str(biomarker_id).upper())
        or measurement_index.get(normalize_alias(biomarker_id))
        or {}
    )

    mapping_index = _index_mapping_entries(model_mapping)
    canonical_id = (
        measurement.get("biomarker")
        or measurement.get("biomarker_id")
        or biomarker_id
    )
    mapping_entry = dict(
        mapping_index.get(str(canonical_id))
        or mapping_index.get(normalize_alias(str(canonical_id)))
        or model_mapping.get(str(canonical_id), {})
        or {}
    )
    mapping_entry.setdefault("biomarker", str(canonical_id))
    mapping_entry.setdefault("biomarker_id", str(canonical_id))

    resolved = {**measurement, **mapping_entry}
    resolved["biomarker_id"] = resolved.get("biomarker_id") or resolved.get("biomarker") or biomarker_id
    resolved["biomarker"] = resolved.get("biomarker") or resolved["biomarker_id"]

    external_sources = external_sources or {}
    for source_name, source_index in external_sources.items():
        candidate = (
            source_index.get(str(canonical_id))
            or source_index.get(normalize_alias(str(canonical_id)))
            or source_index.get(biomarker_id)
            or source_index.get(normalize_alias(biomarker_id))
        )
        if not candidate:
            continue
        resolved.setdefault("external_sources", {})[source_name] = candidate
        for key, value in candidate.items():
            if key not in resolved or resolved[key] in (None, "", []):
                resolved[key] = value

    resolved["coverage_status"] = classify_source_status(resolved)
    resolved.setdefault("confidence_metadata", {})
    if measurement:
        resolved["confidence_metadata"]["measurement_match"] = "direct_or_alias"
    return resolved


def resolve_all(model_mapping: dict[str, Any], measurement_registry: list[dict[str, Any]] | None = None, external_sources: dict[str, dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    measurement_registry = measurement_registry or []
    return [resolve_biomarker(entry["biomarker_id"], measurement_registry, model_mapping, external_sources) for entry in _iter_mapping_entries(model_mapping)]


def _load_measurement_file(path: str | None) -> list[dict[str, Any]]:
    if not path:
        return []
    p = Path(path)
    if p.suffix.lower() in {".yaml", ".yml", ".json"}:
        data = load_yaml_registry(p) if p.suffix.lower() in {".yaml", ".yml"} else json.loads(p.read_text(encoding="utf-8"))
        return data if isinstance(data, list) else data.get("records", [])
    if p.suffix.lower() == ".csv":
        import pandas as pd

        return pd.read_csv(p).to_dict(orient="records")
    raise ValueError(f"Unsupported measurement registry file: {path}")


def _load_external_sources(paths: list[str] | None) -> dict[str, dict[str, Any]]:
    if not paths:
        return {}
    out: dict[str, dict[str, Any]] = {}
    for raw_path in paths:
        p = Path(raw_path)
        source_name = p.stem
        data = _load_measurement_file(raw_path)
        out[source_name] = _index_external_records(data)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Resolve biomarker mapping records and generate a coverage report.")
    parser.add_argument("--mapping", required=True, help="Path to biomarker_mapping.json")
    parser.add_argument("--nhanes", help="Optional NHANES long CSV or registry file")
    parser.add_argument(
        "--external",
        action="append",
        default=[],
        help="Optional external source file (JSON/YAML/CSV). Repeat for multiple sources.",
    )
    parser.add_argument("--out", required=True, help="Output JSON path for resolved records")
    parser.add_argument("--report", help="Optional JSON path for coverage report")
    args = parser.parse_args(argv)

    mapping = load_json_mapping(args.mapping)
    measurements = _load_measurement_file(args.nhanes)
    external_sources = _load_external_sources(args.external)
    resolved = resolve_all(mapping, measurements, external_sources=external_sources)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(resolved, indent=2, sort_keys=True), encoding="utf-8")
    report = coverage_report(resolved)
    report_path = Path(args.report) if args.report else out.with_suffix(".coverage.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
