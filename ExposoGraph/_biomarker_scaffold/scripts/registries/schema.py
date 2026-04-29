"""Typed registry schemas for biomarker evidence records.

The registry layer intentionally separates measurement evidence, chemical identity,
mechanistic interpretation, and model parameters. Unknown/extra fields are kept in
``raw`` so source-specific metadata is not lost during validation or resolution.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

CoverageStatus = Literal[
    "nhanes_native",
    "literature_biomarker",
    "scenario_proxy",
    "model_proxy",
    "derived_proxy",
    "unsupported",
]
IssueLevel = Literal["error", "warning", "info"]


@dataclass(slots=True)
class EvidenceRecord:
    """A single source-backed evidence item."""

    source: str
    source_type: str | None = None
    identifier: str | None = None
    url: str | None = None
    citation: str | None = None
    confidence: float | None = None
    notes: str | None = None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class BiomarkerIdentity:
    biomarker_id: str
    canonical_name: str
    display_name: str | None = None
    matrix: str | None = None
    chemical_class: str | None = None
    source_status: str | None = None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class MeasurementSource:
    source: str
    cycle: str | None = None
    class_name: str | None = None
    variables: list[dict[str, Any]] = field(default_factory=list)
    units: str | None = None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ChemicalIdentity:
    parent_compound: str | None = None
    pubchem_cid: str | None = None
    comptox_dtxsid: str | None = None
    cas: str | None = None
    mw_g_mol: float | None = None
    synonyms: list[str] = field(default_factory=list)
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class MechanisticMapping:
    target_tissue: str | None = None
    target_enzyme: str | None = None
    Km_uM: float | None = None
    partition_coefficient: float | None = None
    S_over_Km_central: float | None = None
    S_over_Km_range: list[float] | None = None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ModelParameters:
    tier2_multiplier: float | None = None
    reference_range: list[float] | None = None
    reference_units: str | None = None
    lifestyle: str | None = None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class BiomarkerRecord:
    biomarker_id: str
    canonical_name: str
    matrix: str
    chemical_class: str
    source_status: str
    evidence: list[EvidenceRecord] = field(default_factory=list)
    display_name: str | None = None
    parent_compound: str | None = None
    pubchem_cid: str | None = None
    comptox_dtxsid: str | None = None
    cas: str | None = None
    mw_g_mol: float | None = None
    nhanes_variables: list[dict[str, Any]] = field(default_factory=list)
    target_tissue: str | None = None
    target_enzyme: str | None = None
    Km_uM: float | None = None
    partition_coefficient: float | None = None
    S_over_Km_central: float | None = None
    S_over_Km_range: list[float] | None = None
    tier2_multiplier: float | None = None
    reference_range: list[float] | None = None
    reference_units: str | None = None
    references: list[str] = field(default_factory=list)
    provenance_note: str | None = None
    coverage_status: CoverageStatus | None = None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ValidationIssue:
    level: IssueLevel
    biomarker_id: str | None
    field: str | None
    message: str


@dataclass(slots=True)
class RegistryValidationReport:
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.level == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.level == "warning"]

    @property
    def infos(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.level == "info"]

    @property
    def ok(self) -> bool:
        return not self.errors

    def as_dict(self) -> dict[str, Any]:
        return {
            "ok": self.ok,
            "errors": [issue.__dict__ for issue in self.errors],
            "warnings": [issue.__dict__ for issue in self.warnings],
            "infos": [issue.__dict__ for issue in self.infos],
        }
