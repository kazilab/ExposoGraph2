"""Typed schema objects for mechanism-resolved interaction parameters.

Phase 3 defines data containers only. The schema intentionally carries unknown
biology and SME review metadata without turning either into active behavior.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from typing import Any


JsonDict = dict[str, Any]


class ValueEnum(str, Enum):
    """String enum with JSON-friendly values."""

    def __str__(self) -> str:
        return self.value


class ReactionRole(ValueEnum):
    BIOACTIVATION = "bioactivation"
    DETOXIFICATION = "detoxification"
    CLEARANCE = "clearance"
    DUAL_ROLE = "dual_role"
    PROBE_ONLY = "probe_only"
    UNKNOWN = "unknown"


class RiskDirectionIfFluxDecreases(ValueEnum):
    DECREASE = "decrease"
    INCREASE = "increase"
    MIXED = "mixed"
    NEUTRAL = "neutral"
    UNKNOWN = "unknown"


class EvidenceGrade(ValueEnum):
    CURATED = "curated"
    INFERRED_FROM_LOCAL_METADATA = "inferred_from_local_metadata"
    PLACEHOLDER = "placeholder"
    UNKNOWN = "unknown"


class RiskEndpoint(ValueEnum):
    DNA_ADDUCT = "DNA_adduct"
    REACTIVE_EPOXIDE = "reactive_epoxide"
    QUINONE_REDOX = "quinone_redox"
    ALDEHYDE = "aldehyde"
    ROS = "ROS"
    GSH_CONJUGATION = "GSH_conjugation"
    NONE = "none"
    UNKNOWN = "unknown"


class SMEReviewStatus(ValueEnum):
    CURATED = "curated"
    CANDIDATE = "candidate"
    PENDING_TEAM_AGREEMENT = "pending_team_agreement"
    DEFERRED_3_0 = "deferred_3_0"
    UNKNOWN = "unknown"


class ReleaseTarget(ValueEnum):
    V2_0 = "v2_0"
    V3_0 = "v3_0"
    FUTURE = "future"
    UNKNOWN = "unknown"


def enum_from_value(enum_cls: type[ValueEnum], value: Any, default: ValueEnum) -> ValueEnum:
    """Parse a local JSON value into an enum without inventing semantics."""
    if isinstance(value, enum_cls):
        return value
    if value is None:
        return default
    try:
        return enum_cls(str(value))
    except ValueError:
        return default


def to_jsonable(value: Any) -> Any:
    """Convert dataclasses/enums recursively into JSON-compatible values."""
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return {field.name: to_jsonable(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, list):
        return [to_jsonable(item) for item in value]
    if isinstance(value, tuple):
        return [to_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    return value


@dataclass
class SerializableRecord:
    """Mixin for simple dict serialization."""

    def to_dict(self) -> JsonDict:
        return to_jsonable(self)


@dataclass
class EvidenceRecord(SerializableRecord):
    source: str | None = None
    grade: EvidenceGrade = EvidenceGrade.UNKNOWN
    confidence: str | None = None
    provenance_ref: str | None = None
    notes: str | None = None
    metadata: JsonDict | None = None


@dataclass
class ParameterUncertainty(SerializableRecord):
    confidence: str | None = None
    lower: float | None = None
    upper: float | None = None
    unit: str | None = None
    method: str | None = None
    notes: str | None = None


@dataclass
class AssumptionWarning(SerializableRecord):
    code: str
    message: str
    severity: str = "warning"
    field: str | None = None
    review_status: SMEReviewStatus = SMEReviewStatus.UNKNOWN


@dataclass
class SMEReviewNote(SerializableRecord):
    status: SMEReviewStatus = SMEReviewStatus.UNKNOWN
    source: str | None = None
    release_target: ReleaseTarget = ReleaseTarget.UNKNOWN
    tissue_context: str | None = None
    endpoint_context: str | None = None
    team_agreement_status: str | None = None
    notes: str | None = None
    assumption_warnings: list[AssumptionWarning] | None = None


@dataclass
class TissueContext(SerializableRecord):
    tissue: str | None = None
    species: str | None = None
    endpoint: RiskEndpoint = RiskEndpoint.UNKNOWN
    exposure_context: str | None = None
    metadata: JsonDict | None = None


@dataclass
class KineticParameterSet(SerializableRecord):
    km_uM: float | None = None
    ki_uM: float | None = None
    vmax_relative: float | None = None
    hill_coefficient: float | None = None
    relative_priority: int | None = None
    assumed_ki: bool | None = None
    product: str | None = None
    product_hazard: JsonDict | None = None
    evidence: EvidenceRecord | None = None
    uncertainty: ParameterUncertainty | None = None
    metadata: JsonDict | None = None


@dataclass
class MetabolicReaction(SerializableRecord):
    enzyme: str
    substrate: str
    product: str | None = None
    kinetic_parameters: KineticParameterSet | None = None
    reaction_role: ReactionRole = ReactionRole.UNKNOWN
    risk_direction_if_flux_decreases: RiskDirectionIfFluxDecreases = RiskDirectionIfFluxDecreases.UNKNOWN
    risk_endpoints: list[RiskEndpoint] | None = None
    tissue_context: TissueContext | None = None
    evidence: EvidenceRecord | list[EvidenceRecord] | None = None
    warnings: list[AssumptionWarning] | None = None
    sme_notes: list[SMEReviewNote] | None = None
    metadata: JsonDict | None = None


@dataclass
class CompetitiveInteraction(SerializableRecord):
    enzyme: str
    substrate: str
    kinetic_parameters: KineticParameterSet | None = None
    reaction_role: ReactionRole = ReactionRole.UNKNOWN
    risk_direction_if_flux_decreases: RiskDirectionIfFluxDecreases = RiskDirectionIfFluxDecreases.UNKNOWN
    evidence: EvidenceRecord | None = None
    warnings: list[AssumptionWarning] | None = None
    sme_notes: list[SMEReviewNote] | None = None
    metadata: JsonDict | None = None


@dataclass
class GSHConsumer(SerializableRecord):
    name: str
    enzyme: str | None = None
    substrate_class: str | None = None
    gsh_per_umol_substrate: float | None = None
    tissue_context: TissueContext | None = None
    evidence: EvidenceRecord | None = None
    warnings: list[AssumptionWarning] | None = None
    metadata: JsonDict | None = None


@dataclass
class InductionRule(SerializableRecord):
    exposure_context: str
    enzyme: str
    fold_induction: float | None = None
    range_min: float | None = None
    range_max: float | None = None
    mechanism: str | None = None
    tissue_context: TissueContext | None = None
    evidence: EvidenceRecord | None = None
    warnings: list[AssumptionWarning] | None = None
    metadata: JsonDict | None = None


@dataclass
class ProviderResult(SerializableRecord):
    source: str
    items: list[Any]
    warnings: list[AssumptionWarning] | None = None
    metadata: JsonDict | None = None

