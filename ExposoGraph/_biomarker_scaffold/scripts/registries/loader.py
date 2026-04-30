"""Load and write biomarker mapping registry documents."""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any, cast


def load_json_mapping(path: str | Path) -> dict[str, Any]:
    """Load a JSON biomarker mapping document."""
    source = Path(path)
    with source.open("r", encoding="utf-8") as handle:
        loaded = json.load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"Expected JSON object in {source}")
    return loaded


def load_yaml_mapping(path: str | Path) -> dict[str, Any]:
    """Load a YAML registry document."""
    source = Path(path)
    yaml_module = importlib.import_module("yaml")
    with source.open("r", encoding="utf-8") as handle:
        loaded = cast(Any, yaml_module).safe_load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"Expected YAML mapping in {source}")
    return loaded


def load_registry_document(path: str | Path) -> dict[str, Any]:
    """Load a registry document from ``.json``, ``.yaml``, or ``.yml``."""
    source = Path(path)
    if source.suffix.lower() == ".json":
        return load_json_mapping(source)
    if source.suffix.lower() in {".yaml", ".yml"}:
        return load_yaml_mapping(source)
    raise ValueError(f"Unsupported registry file type: {source.suffix}")


def write_json_mapping(path: str | Path, document: dict[str, Any]) -> None:
    """Write a JSON biomarker mapping document with stable formatting."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

