"""Validation helpers for biomarker registries."""

from __future__ import annotations

from collections import Counter
from typing import Iterable, Mapping, Any

from .schema import BiomarkerRecord, RegistryValidationReport, ValidationIssue

_REQUIRED = ["biomarker_id", "matrix", "chemical_class"]
_PROXY_STATUSES = {"scenario_proxy", "model_proxy", "derived_proxy", "literature_biomarker"}


def _get(record: BiomarkerRecord | Mapping[str, Any], field: str) -> Any:
    return getattr(record, field) if isinstance(record, BiomarkerRecord) else record.get(field)


def validate_biomarker_record(record: BiomarkerRecord | Mapping[str, Any]) -> list[ValidationIssue]:
    biomarker_id = _get(record, "biomarker_id")
    issues: list[ValidationIssue] = []
    for field in _REQUIRED:
        if _get(record, field) in (None, ""):
            issues.append(ValidationIssue("error", biomarker_id, field, f"Missing required field: {field}"))
    if _get(record, "source_status") in (None, ""):
        issues.append(ValidationIssue("error", biomarker_id, "source_status", "Missing required field: source_status"))
    if not _get(record, "evidence"):
        issues.append(ValidationIssue("warning", biomarker_id, "evidence", "Missing evidence records"))
    if not (_get(record, "pubchem_cid") or _get(record, "comptox_dtxsid") or _get(record, "cas")):
        issues.append(ValidationIssue("warning", biomarker_id, "external_identifiers", "Missing external chemical identifiers"))
    if _get(record, "Km_uM") is None:
        issues.append(ValidationIssue("warning", biomarker_id, "Km_uM", "Missing enzyme kinetic Km value"))
    if _get(record, "reference_units") is None:
        issues.append(ValidationIssue("warning", biomarker_id, "reference_units", "Missing unit conversion or reference units metadata"))
    status = _get(record, "coverage_status") or _get(record, "source_status")
    if status in _PROXY_STATUSES:
        issues.append(ValidationIssue("info", biomarker_id, "source_status", f"Entry is classified as {status}"))
    return issues


def validate_registry(records: Iterable[BiomarkerRecord | Mapping[str, Any]]) -> RegistryValidationReport:
    records = list(records)
    report = RegistryValidationReport()
    ids = [_get(record, "biomarker_id") for record in records if _get(record, "biomarker_id")]
    for duplicated in sorted([item for item, count in Counter(ids).items() if count > 1]):
        report.issues.append(ValidationIssue("error", duplicated, "biomarker_id", "Duplicate biomarker_id"))
    for record in records:
        report.issues.extend(validate_biomarker_record(record))
    return report
