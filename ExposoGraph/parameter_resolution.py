"""Schema-only containers for future kinetic parameter resolution.

Phase 3 deliberately does not implement Ki resolution, IC50 conversion,
affinity conversion, RDKit fingerprints, ECFP4, Tanimoto fallbacks, or Km-as-Ki
policy decisions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .interaction_schema import EvidenceRecord, ParameterUncertainty, SerializableRecord, TissueContext, ValueEnum


class ParameterSourceKind(ValueEnum):
    MEASURED = "measured"
    CURATED = "curated"
    ASSUMED = "assumed"
    LOCAL_METADATA = "local_metadata"
    LITERATURE = "literature"
    MODEL_OUTPUT = "model_output"
    UNKNOWN = "unknown"


class ParameterResolutionMethod(ValueEnum):
    DIRECT_LOOKUP = "direct_lookup"
    MEASURED_VALUE = "measured_value"
    ASSUMED_EQUAL_KM = "assumed_equal_km"
    IC50_CONVERSION_NOT_IMPLEMENTED = "ic50_conversion_not_implemented"
    KM_AS_KI_POLICY_NOT_IMPLEMENTED = "km_as_ki_policy_not_implemented"
    AFFINITY_CONVERSION_NOT_IMPLEMENTED = "affinity_conversion_not_implemented"
    UNRESOLVED = "unresolved"
    UNKNOWN = "unknown"


class AffinityFallbackStatus(ValueEnum):
    INACTIVE = "inactive"
    NOT_APPLICABLE = "not_applicable"
    PENDING = "pending"
    UNAVAILABLE = "unavailable"
    UNKNOWN = "unknown"


@dataclass
class ParameterResolutionWarning(SerializableRecord):
    code: str
    message: str
    severity: str = "warning"
    source_field: str | None = None


@dataclass
class ResolvedParameter(SerializableRecord):
    name: str
    value: float | None = None
    unit: str | None = None
    source_kind: ParameterSourceKind = ParameterSourceKind.UNKNOWN
    resolution_method: ParameterResolutionMethod = ParameterResolutionMethod.UNRESOLVED
    evidence: EvidenceRecord | None = None
    uncertainty: ParameterUncertainty | None = None
    fallback_status: AffinityFallbackStatus = AffinityFallbackStatus.NOT_APPLICABLE
    warnings: list[ParameterResolutionWarning] | None = None
    metadata: dict[str, Any] | None = None


@dataclass
class KiResolverRequest(SerializableRecord):
    enzyme: str
    inhibitor: str
    target_substrate: str | None = None
    tissue_context: TissueContext | None = None
    endpoint_context: str | None = None
    requested_unit: str = "uM"
    allow_ic50_conversion: bool = False
    allow_affinity_fallback: bool = False
    allow_km_as_ki_policy: bool = False
    metadata: dict[str, Any] | None = None
