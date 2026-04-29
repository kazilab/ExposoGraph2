"""Registry package for biomarker schemas, validation, and resolution."""

from __future__ import annotations

from .loader import load_biomarker_records, load_json_mapping, load_yaml_registry
from .mapping_document import (
    apply_update_list,
    build_mapping_document,
    normalize_mapping_document,
    validate_mapping_document,
)
from .resolver import resolve_all, resolve_biomarker
from .validator import validate_biomarker_record, validate_registry

_BUILD_EXPORTS = {
    "build_biomarker_mapping_document",
    "compare_biomarker_mapping_documents",
    "load_biomarker_mapping_manifest",
    "load_biomarker_mapping_source",
}


def __getattr__(name: str):
    if name in _BUILD_EXPORTS:
        from .build_mapping import (
            build_biomarker_mapping_document,
            compare_biomarker_mapping_documents,
            load_biomarker_mapping_manifest,
            load_biomarker_mapping_source,
        )

        exports = {
            "build_biomarker_mapping_document": build_biomarker_mapping_document,
            "compare_biomarker_mapping_documents": compare_biomarker_mapping_documents,
            "load_biomarker_mapping_manifest": load_biomarker_mapping_manifest,
            "load_biomarker_mapping_source": load_biomarker_mapping_source,
        }
        return exports[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "load_yaml_registry",
    "load_json_mapping",
    "load_biomarker_records",
    "build_mapping_document",
    "normalize_mapping_document",
    "apply_update_list",
    "validate_mapping_document",
    "validate_biomarker_record",
    "validate_registry",
    "resolve_biomarker",
    "resolve_all",
    "build_biomarker_mapping_document",
    "compare_biomarker_mapping_documents",
    "load_biomarker_mapping_manifest",
    "load_biomarker_mapping_source",
]
