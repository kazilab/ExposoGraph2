"""Registry loading helpers."""

from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path
from typing import Any

import yaml

from .schema import BiomarkerRecord, EvidenceRecord


class RegistryLoadError(ValueError):
    """Raised when a registry or mapping file cannot be parsed."""


def _read_structured_file(path: str | Path) -> Any:
    path = Path(path)
    if not path.exists():
        raise RegistryLoadError(f"Registry file does not exist: {path}")
    try:
        text = path.read_text(encoding="utf-8")
        if path.suffix.lower() in {".yaml", ".yml"}:
            return yaml.safe_load(text) or []
        if path.suffix.lower() == ".json":
            return json.loads(text)
        # Try YAML first because JSON is a YAML subset.
        return yaml.safe_load(text) or []
    except (OSError, json.JSONDecodeError, yaml.YAMLError) as exc:
        raise RegistryLoadError(f"Failed to parse {path}: {exc}") from exc


def load_yaml_registry(path: str | Path) -> list[dict[str, Any]]:
    data = _read_structured_file(path)
    if isinstance(data, dict) and "biomarkers" in data:
        data = data["biomarkers"]
    if not isinstance(data, list):
        raise RegistryLoadError(f"Expected list of registry records in {path}")
    if not all(isinstance(item, dict) for item in data):
        raise RegistryLoadError(f"Every registry item in {path} must be an object")
    return data


def load_json_mapping(path: str | Path) -> dict[str, Any]:
    data = _read_structured_file(path)
    if not isinstance(data, dict):
        raise RegistryLoadError(f"Expected mapping object in {path}")
    return data


def _build_evidence(raw_items: Any) -> list[EvidenceRecord]:
    if raw_items is None:
        return []
    if not isinstance(raw_items, list):
        raise RegistryLoadError("evidence must be a list")
    known = {field.name for field in fields(EvidenceRecord)}
    out: list[EvidenceRecord] = []
    for item in raw_items:
        if isinstance(item, str):
            out.append(EvidenceRecord(source=item))
            continue
        if not isinstance(item, dict):
            raise RegistryLoadError("each evidence item must be a string or object")
        kwargs = {k: v for k, v in item.items() if k in known and k != "raw"}
        extras = {k: v for k, v in item.items() if k not in known}
        kwargs.setdefault("raw", extras)
        out.append(EvidenceRecord(**kwargs))
    return out


def _record_from_dict(item: dict[str, Any]) -> BiomarkerRecord:
    known = {field.name for field in fields(BiomarkerRecord)}
    kwargs = {k: v for k, v in item.items() if k in known and k not in {"raw", "evidence"}}
    extras = {k: v for k, v in item.items() if k not in known}
    kwargs["evidence"] = _build_evidence(item.get("evidence", []))
    kwargs["raw"] = extras
    missing = [name for name in ["biomarker_id", "canonical_name", "matrix", "chemical_class", "source_status"] if name not in kwargs or kwargs[name] in (None, "")]
    if missing:
        ident = item.get("biomarker_id", "<unknown>")
        raise RegistryLoadError(f"Record {ident} missing required fields: {', '.join(missing)}")
    return BiomarkerRecord(**kwargs)


def load_biomarker_records(path: str | Path) -> list[BiomarkerRecord]:
    return [_record_from_dict(item) for item in load_yaml_registry(path)]
