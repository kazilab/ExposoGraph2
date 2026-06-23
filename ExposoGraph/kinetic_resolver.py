"""Centralized kinetic-parameter resolver for Phase 5.

The resolver is local-data-only. It preserves the legacy competitive Ki lookup
contract and exposes a reversible-inhibition context resolver for local
parameter resolution and the live interaction integration seam.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isclose, isfinite
from typing import Any

from .interaction_schema import (
    ApplicabilityDomain,
    CompetitiveInteraction,
    ConcentrationBasis,
    EvidenceGrade,
    InhibitionMode,
    ParameterUncertainty,
    ReleaseTarget,
    SerializableRecord,
    enum_from_value,
)
from .parameter_provider import JSONInteractionParameterProvider
from .parameter_resolution import (
    AffinityFallbackStatus,
    InhibitionResolutionStatus,
    ParameterResolutionMethod,
    ParameterResolutionWarning,
    ParameterSourceKind,
    ResolvedParameter,
    ResolvedReversibleInhibition,
    ReversibleInhibitionResolutionRequest,
)
from .reversible_inhibition import compute_reversible_inhibition


MODULE3_KM_STATIC_FOR_2_0 = True
_QUANTITATIVE_CONCENTRATION_BASES = {
    ConcentrationBasis.UNBOUND,
    ConcentrationBasis.INTRACELLULAR,
    ConcentrationBasis.TISSUE_EFFECTIVE,
    ConcentrationBasis.MODEL_DERIVED,
}
_BLOCKING_WARNING_CODES = {
    "UNKNOWN_INHIBITION_MODE",
    "IC50_CONVERSION_REQUIRES_MODE",
    "IC50_CONVERSION_REQUIRES_ASSAY_SUBSTRATE",
    "IC50_CONVERSION_REQUIRES_KM",
    "MIXED_INHIBITION_REQUIRES_TWO_CONSTANTS",
    "CONCENTRATION_BASIS_MISMATCH",
    "INHIBITOR_CONCENTRATION_MISSING",
    "SUBSTRATE_CONTEXT_MISSING",
    "KI_MISSING",
    "KM_MISSING",
    "OUTSIDE_REVERSIBLE_INHIBITION_DOMAIN",
    "MODE_PARAMETER_MISMATCH",
    "PURE_NONCOMPETITIVE_REQUIRES_EQUAL_CONSTANTS",
}


@dataclass
class KiResolutionContext(SerializableRecord):
    """Context carried with a Phase 5 Ki request."""

    enzyme: str | None = None
    inhibitor: str | None = None
    target_substrate: str | None = None
    tissue: str | None = None
    endpoint: str | None = None
    assay_context: dict[str, Any] | None = None
    release_target: ReleaseTarget = ReleaseTarget.V2_0
    allow_ic50_conversion: bool = False
    allow_km_proxy: bool = True
    allow_affinity_fallback: bool = False
    metadata: dict[str, Any] | None = None


class KineticParameterResolver:
    """Resolve kinetic parameters from local provider data with explicit fallbacks."""

    def __init__(self, provider: JSONInteractionParameterProvider | None = None) -> None:
        self.provider = provider or JSONInteractionParameterProvider()

    def get_ki(
        self,
        enzyme: str,
        inhibitor: str,
        target_substrate: str | None = None,
        context: KiResolutionContext | dict[str, Any] | None = None,
    ) -> ResolvedParameter:
        request = _coerce_context(
            context,
            enzyme=enzyme,
            inhibitor=inhibitor,
            target_substrate=target_substrate,
        )
        interaction = self._find_interaction(enzyme, inhibitor)
        warnings: list[ParameterResolutionWarning] = []

        if interaction is None or interaction.kinetic_parameters is None:
            warnings.append(_warning("ki_missing", "No local interaction record with Ki was found.", "Ki_uM"))
            return self._unresolved(
                enzyme=enzyme,
                inhibitor=inhibitor,
                request=request,
                warnings=warnings,
                affinity_requested=request.allow_affinity_fallback,
            )

        kinetics = interaction.kinetic_parameters
        evidence = interaction.evidence or kinetics.evidence
        local_fields = (kinetics.metadata or {}).get("local_fields", {})
        provenance = evidence.metadata if evidence else {}
        metadata = _metadata(enzyme, inhibitor, target_substrate, local_fields, provenance, request)

        if kinetics.ki_uM is not None:
            if _is_curated_ki(evidence, provenance):
                metadata["resolution_order"] = "exact_curated_ki"
                metadata["is_curated_ki"] = True
                return ResolvedParameter(
                    name="Ki",
                    value=kinetics.ki_uM,
                    unit="uM",
                    source_kind=ParameterSourceKind.CURATED,
                    resolution_method=ParameterResolutionMethod.MEASURED_VALUE,
                    evidence=evidence,
                    uncertainty=ParameterUncertainty(
                        confidence=evidence.confidence if evidence else "curated",
                        unit="uM",
                        method="curated_ki_lookup",
                    ),
                    fallback_status=AffinityFallbackStatus.NOT_APPLICABLE,
                    warnings=None,
                    metadata=metadata,
                )
            warnings.append(
                _warning(
                    "ki_provenance_incomplete",
                    "Local Ki is present but curated Ki provenance is incomplete or non-curated.",
                    "Ki_uM",
                )
            )
            metadata["resolution_order"] = "non_curated_local_ki"
            metadata["is_curated_ki"] = False
            return ResolvedParameter(
                name="Ki",
                value=kinetics.ki_uM,
                unit="uM",
                source_kind=ParameterSourceKind.LOCAL_METADATA,
                resolution_method=ParameterResolutionMethod.DIRECT_LOOKUP,
                evidence=evidence,
                uncertainty=ParameterUncertainty(confidence="unknown", unit="uM", method="local_ki_lookup"),
                fallback_status=AffinityFallbackStatus.NOT_APPLICABLE,
                warnings=warnings,
                metadata=metadata,
            )

        warnings.append(_warning("ki_missing", "No curated Ki_uM is present in local data.", "Ki_uM"))
        if request.allow_ic50_conversion:
            warnings.append(
                _warning(
                    "ic50_conversion_unavailable",
                    "IC50-to-Ki conversion was requested, but no sufficient assay context and local IC50 were available.",
                    "IC50_uM",
                )
            )

        if request.allow_km_proxy and kinetics.km_uM is not None:
            warnings.append(
                _warning(
                    "km_used_as_ki_proxy",
                    "Km_uM was used as a low-confidence Ki proxy because curated Ki is missing.",
                    "Km_uM",
                )
            )
            metadata["resolution_order"] = "km_proxy"
            metadata["is_curated_ki"] = False
            metadata["proxy_source_field"] = "Km_uM"
            return ResolvedParameter(
                name="Ki",
                value=kinetics.km_uM,
                unit="uM",
                source_kind=ParameterSourceKind.ASSUMED,
                resolution_method=ParameterResolutionMethod.ASSUMED_EQUAL_KM,
                evidence=evidence,
                uncertainty=ParameterUncertainty(
                    confidence="low",
                    unit="uM",
                    method="km_as_ki_proxy",
                    notes="Phase 5 conservative fallback; not curated Ki.",
                ),
                fallback_status=AffinityFallbackStatus.NOT_APPLICABLE,
                warnings=warnings,
                metadata=metadata,
            )

        return self._unresolved(
            enzyme=enzyme,
            inhibitor=inhibitor,
            request=request,
            warnings=warnings,
            evidence=evidence,
            metadata=metadata,
            affinity_requested=request.allow_affinity_fallback,
        )

    def resolve_reversible_inhibition(
        self,
        request: ReversibleInhibitionResolutionRequest | dict[str, Any],
    ) -> ResolvedReversibleInhibition:
        """Resolve a reversible-inhibition context without live engine integration."""

        return _resolve_reversible_context(_coerce_reversible_request(request))

    def _find_interaction(self, enzyme: str, inhibitor: str):
        for interaction in self.provider.get_competitive_interactions(enzyme):
            if _normalize(interaction.substrate) == _normalize(inhibitor):
                return interaction
        return None

    def _unresolved(
        self,
        *,
        enzyme: str,
        inhibitor: str,
        request: KiResolutionContext,
        warnings: list[ParameterResolutionWarning],
        evidence=None,
        metadata: dict[str, Any] | None = None,
        affinity_requested: bool = False,
    ) -> ResolvedParameter:
        if affinity_requested:
            warnings.append(
                _warning(
                    "affinity_fallback_unavailable",
                    "Affinity fallback was requested but is unavailable in Phase 5.",
                    "allow_affinity_fallback",
                )
            )
        warnings.append(_warning("no_parameter_resolved", "No Ki parameter was resolved.", "Ki_uM"))
        result_metadata = metadata or _metadata(enzyme, inhibitor, request.target_substrate, {}, {}, request)
        result_metadata["resolution_order"] = "unresolved"
        result_metadata["is_curated_ki"] = False
        return ResolvedParameter(
            name="Ki",
            value=None,
            unit="uM",
            source_kind=ParameterSourceKind.LOCAL_METADATA,
            resolution_method=ParameterResolutionMethod.UNRESOLVED,
            evidence=evidence,
            uncertainty=ParameterUncertainty(confidence="unknown", unit="uM", method="unresolved"),
            fallback_status=(
                AffinityFallbackStatus.UNAVAILABLE
                if affinity_requested
                else AffinityFallbackStatus.NOT_APPLICABLE
            ),
            warnings=warnings,
            metadata=result_metadata,
        )


def get_ki(
    enzyme: str,
    inhibitor: str,
    target_substrate: str | None = None,
    context: KiResolutionContext | dict[str, Any] | None = None,
    provider: JSONInteractionParameterProvider | None = None,
) -> ResolvedParameter:
    """Resolve Ki through the centralized Phase 5 local-data path."""

    return KineticParameterResolver(provider=provider).get_ki(
        enzyme,
        inhibitor,
        target_substrate=target_substrate,
        context=context,
    )


def resolve_reversible_inhibition(
    request: ReversibleInhibitionResolutionRequest | dict[str, Any],
    provider: JSONInteractionParameterProvider | None = None,
) -> ResolvedReversibleInhibition:
    """Resolve a reversible-inhibition context and call the reversible-inhibition kernel if valid."""

    return KineticParameterResolver(provider=provider).resolve_reversible_inhibition(request)


def request_from_competitive_interaction(
    interaction: CompetitiveInteraction,
    *,
    inhibitor_concentration_uM: float,
    substrate_concentration_uM: float,
    concentration_basis: ConcentrationBasis | str = ConcentrationBasis.UNKNOWN,
    parameter_concentration_basis: ConcentrationBasis | str = ConcentrationBasis.UNKNOWN,
    target_substrate: str | None = None,
    mode: InhibitionMode | str = InhibitionMode.COMPETITIVE,
    assay_substrate_concentration_uM: float | None = None,
    vmax: float = 1.0,
) -> ReversibleInhibitionResolutionRequest:
    """Adapt an existing competitive record into the reversible resolver contract."""

    kinetics = interaction.kinetic_parameters
    evidence = interaction.evidence or (kinetics.evidence if kinetics else None)
    local_fields = (kinetics.metadata or {}).get("local_fields", {}) if kinetics else {}
    metadata = {
        "adapter": "competitive_interaction_v1",
        "source_kind": "local_json_provider",
        "live_integration": False,
        "local_fields_used": dict(local_fields),
    }
    active_mode = enum_from_value(InhibitionMode, mode, InhibitionMode.UNKNOWN)
    ki_value = kinetics.ki_uM if kinetics else None
    return ReversibleInhibitionResolutionRequest(
        mode=active_mode,
        enzyme=interaction.enzyme,
        inhibitor=interaction.substrate,
        target_substrate=target_substrate,
        km_uM=kinetics.km_uM if kinetics else None,
        ki_free_enzyme_uM=ki_value if active_mode is InhibitionMode.COMPETITIVE else None,
        ki_enzyme_substrate_uM=None,
        inhibitor_concentration_uM=inhibitor_concentration_uM,
        substrate_concentration_uM=substrate_concentration_uM,
        assay_substrate_concentration_uM=assay_substrate_concentration_uM,
        vmax=vmax,
        concentration_basis=enum_from_value(ConcentrationBasis, concentration_basis, ConcentrationBasis.UNKNOWN),
        parameter_concentration_basis=enum_from_value(
            ConcentrationBasis,
            parameter_concentration_basis,
            ConcentrationBasis.UNKNOWN,
        ),
        applicability_domain=ApplicabilityDomain.NOT_ASSESSABLE,
        evidence=evidence,
        metadata=metadata,
    )


def _resolve_reversible_context(request: ReversibleInhibitionResolutionRequest) -> ResolvedReversibleInhibition:
    warnings: list[ParameterResolutionWarning] = []
    assumptions: list[str] = []
    metadata = dict(request.metadata or {})
    metadata.setdefault("resolver", "kinetic_resolver.resolve_reversible_inhibition")
    metadata.setdefault("live_integration", False)

    mode = enum_from_value(InhibitionMode, request.mode, InhibitionMode.UNKNOWN)
    concentration_basis = enum_from_value(
        ConcentrationBasis,
        request.concentration_basis,
        ConcentrationBasis.UNKNOWN,
    )
    parameter_basis = enum_from_value(
        ConcentrationBasis,
        request.parameter_concentration_basis,
        ConcentrationBasis.UNKNOWN,
    )
    applicability = enum_from_value(
        ApplicabilityDomain,
        request.applicability_domain,
        ApplicabilityDomain.NOT_ASSESSABLE,
    )

    km_value, invalid = _numeric_value(request.km_uM, "km_uM", warnings, positive=True)
    ki_e_value, ki_e_invalid = _numeric_value(
        request.ki_free_enzyme_uM,
        "ki_free_enzyme_uM",
        warnings,
        positive=True,
    )
    ki_es_value, ki_es_invalid = _numeric_value(
        request.ki_enzyme_substrate_uM,
        "ki_enzyme_substrate_uM",
        warnings,
        positive=True,
    )
    ic50_value, ic50_invalid = _numeric_value(request.ic50_uM, "ic50_uM", warnings, positive=True)
    inhibitor_value, inhibitor_invalid = _numeric_value(
        request.inhibitor_concentration_uM,
        "inhibitor_concentration_uM",
        warnings,
        nonnegative=True,
    )
    substrate_value, substrate_invalid = _numeric_value(
        request.substrate_concentration_uM,
        "substrate_concentration_uM",
        warnings,
        nonnegative=True,
    )
    assay_substrate_value, assay_invalid = _numeric_value(
        request.assay_substrate_concentration_uM,
        "assay_substrate_concentration_uM",
        warnings,
        positive=True,
    )
    vmax_value, vmax_invalid = _numeric_value(request.vmax, "vmax", warnings, nonnegative=True)
    any_invalid = any(
        [
            invalid,
            ki_e_invalid,
            ki_es_invalid,
            ic50_invalid,
            inhibitor_invalid,
            substrate_invalid,
            assay_invalid,
            vmax_invalid,
        ]
    )

    km_param = _parameter(
        "Km",
        km_value,
        "km_uM",
        ParameterResolutionMethod.MEASURED_VALUE,
        request.evidence,
        ParameterSourceKind.MEASURED,
    )
    ki_e_param = _parameter(
        "Ki_E",
        ki_e_value,
        "ki_free_enzyme_uM",
        ParameterResolutionMethod.MEASURED_VALUE,
        request.evidence,
        ParameterSourceKind.MEASURED,
    )
    ki_es_param = _parameter(
        "Ki_ES",
        ki_es_value,
        "ki_enzyme_substrate_uM",
        ParameterResolutionMethod.MEASURED_VALUE,
        request.evidence,
        ParameterSourceKind.MEASURED,
    )
    inhibitor_param = _parameter(
        "I",
        inhibitor_value,
        "inhibitor_concentration_uM",
        ParameterResolutionMethod.DIRECT_LOOKUP,
        request.evidence,
        ParameterSourceKind.MODEL_OUTPUT,
    )
    substrate_param = _parameter(
        "S",
        substrate_value,
        "substrate_concentration_uM",
        ParameterResolutionMethod.DIRECT_LOOKUP,
        request.evidence,
        ParameterSourceKind.MODEL_OUTPUT,
    )

    derived = False
    if mode is InhibitionMode.UNKNOWN:
        warnings.append(
            _resolution_warning(
                "UNKNOWN_INHIBITION_MODE",
                "Unknown inhibition mode cannot be resolved as competitive by default.",
                "mode",
            )
        )
        if ic50_value is not None:
            warnings.append(
                _resolution_warning(
                    "IC50_CONVERSION_REQUIRES_MODE",
                    "IC50 conversion requires a supported inhibition mode.",
                    "ic50_uM",
                )
            )

    if ic50_value is not None and mode is InhibitionMode.COMPETITIVE and ki_e_value is None:
        if km_value is None:
            warnings.append(_resolution_warning("IC50_CONVERSION_REQUIRES_KM", "Competitive IC50 conversion requires Km.", "km_uM"))
        if assay_substrate_value is None:
            warnings.append(
                _resolution_warning(
                    "IC50_CONVERSION_REQUIRES_ASSAY_SUBSTRATE",
                    "Competitive IC50 conversion requires assay substrate concentration.",
                    "assay_substrate_concentration_uM",
                )
            )
        if km_value is not None and assay_substrate_value is not None:
            ki_e_value = ic50_value / (1.0 + assay_substrate_value / km_value)
            ki_e_param = _parameter(
                "Ki_E",
                ki_e_value,
                "ic50_uM",
                ParameterResolutionMethod.IC50_CHENG_PRUSOFF_COMPETITIVE,
                request.evidence,
                ParameterSourceKind.ASSUMED,
            )
            assumptions.append("competitive_ic50_cheng_prusoff")
            derived = True

    if ic50_value is not None and mode is InhibitionMode.UNCOMPETITIVE and ki_es_value is None:
        if km_value is None:
            warnings.append(_resolution_warning("IC50_CONVERSION_REQUIRES_KM", "Uncompetitive IC50 conversion requires Km.", "km_uM"))
        if assay_substrate_value is None:
            warnings.append(
                _resolution_warning(
                    "IC50_CONVERSION_REQUIRES_ASSAY_SUBSTRATE",
                    "Uncompetitive IC50 conversion requires assay substrate concentration.",
                    "assay_substrate_concentration_uM",
                )
            )
        if km_value is not None and assay_substrate_value is not None:
            ki_es_value = ic50_value / (1.0 + km_value / assay_substrate_value)
            ki_es_param = _parameter(
                "Ki_ES",
                ki_es_value,
                "ic50_uM",
                ParameterResolutionMethod.IC50_CHENG_PRUSOFF_UNCOMPETITIVE,
                request.evidence,
                ParameterSourceKind.ASSUMED,
            )
            assumptions.append("uncompetitive_ic50_cheng_prusoff")
            derived = True

    if ic50_value is not None and mode is InhibitionMode.PURE_NONCOMPETITIVE:
        if ki_e_value is None and ki_es_value is None:
            ki_e_value = ic50_value
            ki_es_value = ic50_value
            ki_e_param = _parameter(
                "Ki_E",
                ki_e_value,
                "ic50_uM",
                ParameterResolutionMethod.IC50_PURE_NONCOMPETITIVE,
                request.evidence,
                ParameterSourceKind.ASSUMED,
            )
            ki_es_param = _parameter(
                "Ki_ES",
                ki_es_value,
                "ic50_uM",
                ParameterResolutionMethod.IC50_PURE_NONCOMPETITIVE,
                request.evidence,
                ParameterSourceKind.ASSUMED,
            )
            assumptions.append("pure_noncompetitive_ic50_equals_ki")
            derived = True

    if ic50_value is not None and mode is InhibitionMode.MIXED and (ki_e_value is None or ki_es_value is None):
        warnings.append(
            _resolution_warning(
                "MIXED_INHIBITION_REQUIRES_TWO_CONSTANTS",
                "A single IC50 cannot infer both mixed-inhibition binding constants.",
                "ic50_uM",
            )
        )

    if mode is InhibitionMode.PURE_NONCOMPETITIVE:
        if ki_e_value is None and ki_es_value is not None:
            ki_e_value = ki_es_value
            ki_e_param = _parameter(
                "Ki_E",
                ki_e_value,
                "ki_enzyme_substrate_uM",
                ParameterResolutionMethod.DIRECT_LOOKUP,
                request.evidence,
                ParameterSourceKind.ASSUMED,
            )
            assumptions.append("pure_noncompetitive_equal_ki_applied_to_free_enzyme_arm")
        if ki_es_value is None and ki_e_value is not None:
            ki_es_value = ki_e_value
            ki_es_param = _parameter(
                "Ki_ES",
                ki_es_value,
                "ki_free_enzyme_uM",
                ParameterResolutionMethod.DIRECT_LOOKUP,
                request.evidence,
                ParameterSourceKind.ASSUMED,
            )
            assumptions.append("pure_noncompetitive_equal_ki_applied_to_enzyme_substrate_arm")
        if ki_e_value is not None and ki_es_value is not None and not isclose(ki_e_value, ki_es_value):
            warnings.append(
                _resolution_warning(
                    "PURE_NONCOMPETITIVE_REQUIRES_EQUAL_CONSTANTS",
                    "Pure non-competitive inhibition requires equal Ki values for both arms.",
                    "ki_free_enzyme_uM",
                )
            )

    if mode is InhibitionMode.COMPETITIVE and ki_es_value is not None:
        warnings.append(
            _resolution_warning(
                "MODE_PARAMETER_MISMATCH",
                "Competitive inhibition must not include a Ki_ES arm for quantitative resolution.",
                "ki_enzyme_substrate_uM",
            )
        )
    if mode is InhibitionMode.UNCOMPETITIVE and ki_e_value is not None:
        warnings.append(
            _resolution_warning(
                "MODE_PARAMETER_MISMATCH",
                "Uncompetitive inhibition must not include a Ki_E arm for quantitative resolution.",
                "ki_free_enzyme_uM",
            )
        )

    _add_missing_mode_warnings(mode, ki_e_value, ki_es_value, warnings)
    if km_value is None:
        warnings.append(_resolution_warning("KM_MISSING", "Km is required for reversible-inhibition resolution.", "km_uM"))
    if inhibitor_value is None:
        warnings.append(
            _resolution_warning(
                "INHIBITOR_CONCENTRATION_MISSING",
                "Local inhibitor concentration is required for quantitative resolution.",
                "inhibitor_concentration_uM",
            )
        )
    if substrate_value is None:
        warnings.append(
            _resolution_warning(
                "SUBSTRATE_CONTEXT_MISSING",
                "Substrate concentration is required for quantitative resolution.",
                "substrate_concentration_uM",
            )
        )

    if not _basis_permits_quantitative(concentration_basis, parameter_basis):
        warnings.append(
            _resolution_warning(
                "CONCENTRATION_BASIS_MISMATCH",
                "Nominal, total, external, or unknown concentration basis cannot be treated as local quantitative concentration.",
                "concentration_basis",
            )
        )

    substrate_to_km_ratio = None
    if km_value is not None and substrate_value is not None:
        substrate_to_km_ratio = substrate_value / km_value

    warning_codes = {warning.code for warning in warnings}
    status = _resolution_status(any_invalid, warning_codes, derived)
    if ApplicabilityDomain.OUTSIDE_DOMAIN is applicability and status in {
        InhibitionResolutionStatus.RESOLVED_DIRECT,
        InhibitionResolutionStatus.RESOLVED_DERIVED,
    }:
        warnings.append(
            _resolution_warning(
                "OUTSIDE_REVERSIBLE_INHIBITION_DOMAIN",
                "Caller supplied an outside-domain applicability state.",
                "applicability_domain",
            )
        )
        status = InhibitionResolutionStatus.REVIEW_REQUIRED
        warning_codes = {warning.code for warning in warnings}

    kernel_result = None
    applicability = _applicability_for_status(status, applicability, warning_codes)
    if status in {InhibitionResolutionStatus.RESOLVED_DIRECT, InhibitionResolutionStatus.RESOLVED_DERIVED}:
        try:
            ki_e_kernel = ki_e_value if mode in {InhibitionMode.COMPETITIVE, InhibitionMode.PURE_NONCOMPETITIVE, InhibitionMode.MIXED} else None
            ki_es_kernel = ki_es_value if mode in {InhibitionMode.UNCOMPETITIVE, InhibitionMode.PURE_NONCOMPETITIVE, InhibitionMode.MIXED} else None
            kernel_result = compute_reversible_inhibition(
                mode=mode,
                substrate_concentration=substrate_value,
                km=km_value,
                vmax=vmax_value,
                inhibitor_concentration=inhibitor_value,
                ki_free_enzyme=ki_e_kernel,
                ki_enzyme_substrate=ki_es_kernel,
            )
        except ValueError as exc:
            warnings.append(
                _resolution_warning(
                    "OUTSIDE_REVERSIBLE_INHIBITION_DOMAIN",
                    str(exc),
                    "reversible_inhibition_kernel",
                )
            )
            status = InhibitionResolutionStatus.INVALID
            applicability = ApplicabilityDomain.OUTSIDE_DOMAIN

    return ResolvedReversibleInhibition(
        status=status,
        mode=mode,
        applicability_domain=applicability,
        ki_free_enzyme_uM=ki_e_param if ki_e_value is not None else None,
        ki_enzyme_substrate_uM=ki_es_param if ki_es_value is not None else None,
        km_uM=km_param if km_value is not None else None,
        inhibitor_concentration_uM=inhibitor_param if inhibitor_value is not None else None,
        substrate_concentration_uM=substrate_param if substrate_value is not None else None,
        substrate_to_km_ratio=substrate_to_km_ratio,
        concentration_basis=concentration_basis,
        parameter_concentration_basis=parameter_basis,
        kernel_result=kernel_result,
        warnings=warnings,
        assumptions=assumptions,
        evidence=request.evidence,
        metadata=metadata,
    )


def _coerce_context(
    context: KiResolutionContext | dict[str, Any] | None,
    *,
    enzyme: str,
    inhibitor: str,
    target_substrate: str | None,
) -> KiResolutionContext:
    if isinstance(context, KiResolutionContext):
        return KiResolutionContext(
            enzyme=context.enzyme or enzyme,
            inhibitor=context.inhibitor or inhibitor,
            target_substrate=context.target_substrate or target_substrate,
            tissue=context.tissue,
            endpoint=context.endpoint,
            assay_context=context.assay_context,
            release_target=context.release_target,
            allow_ic50_conversion=context.allow_ic50_conversion,
            allow_km_proxy=context.allow_km_proxy,
            allow_affinity_fallback=context.allow_affinity_fallback,
            metadata=context.metadata,
        )
    if isinstance(context, dict):
        return KiResolutionContext(
            enzyme=str(context.get("enzyme") or enzyme),
            inhibitor=str(context.get("inhibitor") or inhibitor),
            target_substrate=context.get("target_substrate") or target_substrate,
            tissue=context.get("tissue"),
            endpoint=context.get("endpoint"),
            assay_context=context.get("assay_context"),
            release_target=context.get("release_target", ReleaseTarget.V2_0),
            allow_ic50_conversion=bool(context.get("allow_ic50_conversion", False)),
            allow_km_proxy=bool(context.get("allow_km_proxy", True)),
            allow_affinity_fallback=bool(context.get("allow_affinity_fallback", False)),
            metadata=context.get("metadata"),
        )
    return KiResolutionContext(enzyme=enzyme, inhibitor=inhibitor, target_substrate=target_substrate)


def _coerce_reversible_request(
    request: ReversibleInhibitionResolutionRequest | dict[str, Any],
) -> ReversibleInhibitionResolutionRequest:
    if isinstance(request, ReversibleInhibitionResolutionRequest):
        source = request.to_dict()
        evidence = request.evidence
    else:
        source = dict(request)
        evidence = source.get("evidence")
    return ReversibleInhibitionResolutionRequest(
        mode=enum_from_value(InhibitionMode, source.get("mode"), InhibitionMode.UNKNOWN),
        enzyme=source.get("enzyme"),
        inhibitor=source.get("inhibitor"),
        target_substrate=source.get("target_substrate"),
        km_uM=source.get("km_uM"),
        ki_free_enzyme_uM=source.get("ki_free_enzyme_uM"),
        ki_enzyme_substrate_uM=source.get("ki_enzyme_substrate_uM"),
        ic50_uM=source.get("ic50_uM"),
        inhibitor_concentration_uM=source.get("inhibitor_concentration_uM"),
        substrate_concentration_uM=source.get("substrate_concentration_uM"),
        assay_substrate_concentration_uM=source.get("assay_substrate_concentration_uM"),
        vmax=source.get("vmax", 1.0),
        concentration_basis=enum_from_value(
            ConcentrationBasis,
            source.get("concentration_basis"),
            ConcentrationBasis.UNKNOWN,
        ),
        parameter_concentration_basis=enum_from_value(
            ConcentrationBasis,
            source.get("parameter_concentration_basis"),
            ConcentrationBasis.UNKNOWN,
        ),
        applicability_domain=enum_from_value(
            ApplicabilityDomain,
            source.get("applicability_domain"),
            ApplicabilityDomain.NOT_ASSESSABLE,
        ),
        evidence=evidence,
        metadata=source.get("metadata"),
    )


def _numeric_value(
    value: Any,
    field_name: str,
    warnings: list[ParameterResolutionWarning],
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> tuple[float | None, bool]:
    if value is None:
        return None, False
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        warnings.append(_invalid_numeric_warning(field_name))
        return None, True
    if not isfinite(numeric):
        warnings.append(_invalid_numeric_warning(field_name))
        return None, True
    if positive and numeric <= 0.0:
        warnings.append(_invalid_numeric_warning(field_name))
        return None, True
    if nonnegative and numeric < 0.0:
        warnings.append(_invalid_numeric_warning(field_name))
        return None, True
    if numeric == 0.0:
        return 0.0, False
    return float(numeric), False


def _parameter(
    name: str,
    value: float | None,
    source_field: str,
    method: ParameterResolutionMethod,
    evidence,
    source_kind: ParameterSourceKind,
) -> ResolvedParameter | None:
    if value is None:
        return None
    return ResolvedParameter(
        name=name,
        value=value,
        unit="uM",
        source_kind=source_kind,
        resolution_method=method,
        evidence=evidence,
        uncertainty=ParameterUncertainty(confidence="resolved", unit="uM", method=str(method)),
        fallback_status=AffinityFallbackStatus.NOT_APPLICABLE,
        warnings=None,
        metadata={"source_field": source_field},
    )


def _add_missing_mode_warnings(
    mode: InhibitionMode,
    ki_e_value: float | None,
    ki_es_value: float | None,
    warnings: list[ParameterResolutionWarning],
) -> None:
    if mode is InhibitionMode.COMPETITIVE and ki_e_value is None:
        warnings.append(_resolution_warning("KI_MISSING", "Competitive inhibition requires Ki_E.", "ki_free_enzyme_uM"))
    if mode is InhibitionMode.UNCOMPETITIVE and ki_es_value is None:
        warnings.append(_resolution_warning("KI_MISSING", "Uncompetitive inhibition requires Ki_ES.", "ki_enzyme_substrate_uM"))
    if mode is InhibitionMode.PURE_NONCOMPETITIVE and (ki_e_value is None or ki_es_value is None):
        warnings.append(
            _resolution_warning(
                "KI_MISSING",
                "Pure non-competitive inhibition requires an explicit Ki for both equivalent arms.",
                "ki_free_enzyme_uM",
            )
        )
    if mode is InhibitionMode.MIXED and (ki_e_value is None or ki_es_value is None):
        warnings.append(
            _resolution_warning(
                "KI_MISSING",
                "Mixed inhibition requires direct Ki_E and Ki_ES constants.",
                "ki_free_enzyme_uM",
            )
        )


def _basis_permits_quantitative(
    concentration_basis: ConcentrationBasis,
    parameter_basis: ConcentrationBasis,
) -> bool:
    if concentration_basis not in _QUANTITATIVE_CONCENTRATION_BASES:
        return False
    if parameter_basis in {ConcentrationBasis.NOMINAL, ConcentrationBasis.TOTAL}:
        return False
    return True


def _resolution_status(
    any_invalid: bool,
    warning_codes: set[str],
    derived: bool,
) -> InhibitionResolutionStatus:
    if any_invalid:
        return InhibitionResolutionStatus.INVALID
    if warning_codes & _BLOCKING_WARNING_CODES:
        return InhibitionResolutionStatus.REVIEW_REQUIRED
    if derived:
        return InhibitionResolutionStatus.RESOLVED_DERIVED
    return InhibitionResolutionStatus.RESOLVED_DIRECT


def _applicability_for_status(
    status: InhibitionResolutionStatus,
    requested: ApplicabilityDomain,
    warning_codes: set[str],
) -> ApplicabilityDomain:
    if status is InhibitionResolutionStatus.INVALID:
        return ApplicabilityDomain.OUTSIDE_DOMAIN
    if "CONCENTRATION_BASIS_MISMATCH" in warning_codes:
        return ApplicabilityDomain.OUTSIDE_DOMAIN
    if status is InhibitionResolutionStatus.REVIEW_REQUIRED:
        return ApplicabilityDomain.NOT_ASSESSABLE
    if requested in {ApplicabilityDomain.IN_DOMAIN, ApplicabilityDomain.CONDITIONALLY_IN_DOMAIN}:
        return requested
    return ApplicabilityDomain.IN_DOMAIN


def _invalid_numeric_warning(field_name: str) -> ParameterResolutionWarning:
    return _resolution_warning(
        "OUTSIDE_REVERSIBLE_INHIBITION_DOMAIN",
        f"{field_name} must be finite and within the reversible-inhibition domain.",
        field_name,
        severity="error",
    )


def _is_curated_ki(evidence, provenance: dict[str, Any] | None) -> bool:
    if evidence and evidence.grade is EvidenceGrade.CURATED:
        return True
    return bool(provenance and provenance.get("ki_status") == "curated")


def _metadata(
    enzyme: str,
    inhibitor: str,
    target_substrate: str | None,
    local_fields: dict[str, Any],
    provenance: dict[str, Any] | None,
    request: KiResolutionContext,
) -> dict[str, Any]:
    return {
        "enzyme": enzyme,
        "inhibitor": inhibitor,
        "target_substrate": target_substrate,
        "source_kind": "local_json_provider",
        "local_fields_used": dict(local_fields),
        "provenance_fields_used": dict(provenance or {}),
        "context": {
            "tissue": request.tissue,
            "endpoint": request.endpoint,
            "assay_context": request.assay_context,
            "release_target": str(request.release_target),
            "allow_ic50_conversion": request.allow_ic50_conversion,
            "allow_km_proxy": request.allow_km_proxy,
            "allow_affinity_fallback": request.allow_affinity_fallback,
        },
        "module3_km_static_for_2_0": MODULE3_KM_STATIC_FOR_2_0,
    }


def _warning(code: str, message: str, source_field: str) -> ParameterResolutionWarning:
    return ParameterResolutionWarning(code=code, message=message, source_field=source_field)


def _resolution_warning(
    code: str,
    message: str,
    source_field: str,
    *,
    severity: str = "warning",
) -> ParameterResolutionWarning:
    return ParameterResolutionWarning(code=code, message=message, severity=severity, source_field=source_field)


def _normalize(value: object) -> str:
    return "".join(character.lower() for character in str(value) if character.isalnum())
