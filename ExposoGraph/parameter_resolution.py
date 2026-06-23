"""Typed containers for local kinetic parameter resolution.

The legacy Ki lookup contract remains available for competitive interactions.
Reversible-inhibition resolution records direct and derived parameters without
wiring those calculations into live platform execution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .interaction_schema import (
    ApplicabilityDomain,
    ConcentrationBasis,
    EvidenceRecord,
    InhibitionMode,
    ParameterUncertainty,
    SerializableRecord,
    TissueContext,
    ValueEnum,
)


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
    IC50_CHENG_PRUSOFF_COMPETITIVE = "ic50_cheng_prusoff_competitive"
    IC50_CHENG_PRUSOFF_UNCOMPETITIVE = "ic50_cheng_prusoff_uncompetitive"
    IC50_PURE_NONCOMPETITIVE = "ic50_pure_noncompetitive"
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


class InhibitionResolutionStatus(ValueEnum):
    RESOLVED_DIRECT = "resolved_direct"
    RESOLVED_DERIVED = "resolved_derived"
    SENSITIVITY_ONLY = "sensitivity_only"
    QUALITATIVE_ONLY = "qualitative_only"
    REVIEW_REQUIRED = "review_required"
    INVALID = "invalid"


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


@dataclass
class ReversibleInhibitionResolutionRequest(SerializableRecord):
    mode: InhibitionMode = InhibitionMode.UNKNOWN
    enzyme: str | None = None
    inhibitor: str | None = None
    target_substrate: str | None = None
    km_uM: float | None = None
    ki_free_enzyme_uM: float | None = None
    ki_enzyme_substrate_uM: float | None = None
    ic50_uM: float | None = None
    inhibitor_concentration_uM: float | None = None
    substrate_concentration_uM: float | None = None
    assay_substrate_concentration_uM: float | None = None
    vmax: float = 1.0
    concentration_basis: ConcentrationBasis = ConcentrationBasis.UNKNOWN
    parameter_concentration_basis: ConcentrationBasis = ConcentrationBasis.UNKNOWN
    applicability_domain: ApplicabilityDomain = ApplicabilityDomain.NOT_ASSESSABLE
    evidence: EvidenceRecord | None = None
    metadata: dict[str, Any] | None = None


@dataclass
class ResolvedReversibleInhibition(SerializableRecord):
    status: InhibitionResolutionStatus
    mode: InhibitionMode
    applicability_domain: ApplicabilityDomain
    ki_free_enzyme_uM: ResolvedParameter | None = None
    ki_enzyme_substrate_uM: ResolvedParameter | None = None
    km_uM: ResolvedParameter | None = None
    inhibitor_concentration_uM: ResolvedParameter | None = None
    substrate_concentration_uM: ResolvedParameter | None = None
    substrate_to_km_ratio: float | None = None
    concentration_basis: ConcentrationBasis = ConcentrationBasis.UNKNOWN
    parameter_concentration_basis: ConcentrationBasis = ConcentrationBasis.UNKNOWN
    kernel_result: Any | None = None
    warnings: list[ParameterResolutionWarning] | None = None
    assumptions: list[str] | None = None
    evidence: EvidenceRecord | None = None
    metadata: dict[str, Any] | None = None
