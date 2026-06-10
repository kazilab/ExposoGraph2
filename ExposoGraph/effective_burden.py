"""Standalone effective carcinogenic burden integration for Phase 8.

This module combines already-computed mechanism ratios into a transparent
semi-mechanistic relative burden. It does not import interaction_engine, modify
public risk outputs, compute kinetic factors, or implement attribution logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Any

from .endpoint_toxic_flux import EndpointToxicFluxResult
from .gsh_redox_capacity import (
    GSHRedoxCapacityInput,
    GSHRedoxCapacityResult,
    compute_gsh_redox_capacity,
)
from .interaction_schema import (
    AssumptionWarning,
    EvidenceGrade,
    EvidenceRecord,
    JsonDict,
    SMEReviewStatus,
    SerializableRecord,
)


EffectiveBurdenWarning = AssumptionWarning


@dataclass
class SusceptibilityModifier(SerializableRecord):
    """Caller-supplied relative susceptibility modifier."""

    modifier_ratio: float | None = 1.0
    label: str | None = None
    evidence: EvidenceRecord | None = None
    warnings: list[AssumptionWarning] | None = None
    metadata: JsonDict | None = None


@dataclass
class GSHBurdenCouplingInput(SerializableRecord):
    """Inputs for standalone Phase 8 GSH consumption coupling."""

    gsh_relevant: bool = False
    base_gsh_consumption_load: float | None = 0.0
    upstream_activation_burden_ratio: float | None = None
    d_factor: float | None = None
    k_factor: float | None = None
    tissue: str | None = None
    synthesis_capacity: float | None = None
    turnover_capacity: float | None = None
    baseline_capacity: float | None = None
    evidence: EvidenceRecord | None = None
    metadata: JsonDict | None = None


@dataclass
class GSHBurdenCouplingResult(SerializableRecord):
    """Result of scaling GSH consumption by upstream activation burden."""

    gsh_relevant: bool
    base_gsh_consumption_load: float
    upstream_activation_burden_ratio: float
    gsh_consumption_load_scaled: float
    scaling_source: str
    d_factor: float | None
    k_factor: float | None
    gsh_redox_capacity_result: GSHRedoxCapacityResult | None
    warnings: list[AssumptionWarning]
    evidence: EvidenceRecord | None = None
    metadata: JsonDict | None = None

    @property
    def gsh_fraction(self) -> float | None:
        if self.gsh_redox_capacity_result is None:
            return None
        return self.gsh_redox_capacity_result.gsh_fraction

    @property
    def redox_capacity_ratio(self) -> float | None:
        if self.gsh_redox_capacity_result is None:
            return None
        return self.gsh_redox_capacity_result.redox_capacity_ratio

    @property
    def detox_penalty_multiplier(self) -> float:
        if self.gsh_redox_capacity_result is None:
            return 1.0
        return self.gsh_redox_capacity_result.detox_penalty_multiplier


@dataclass
class EffectiveBurdenInput(SerializableRecord):
    """Inputs for the standalone Phase 8 effective burden ratio."""

    activation_burden_ratio: float | None = None
    detox_failure_ratio: float | None = None
    susceptibility_modifier: float | SusceptibilityModifier | None = None
    gsh_relevant: bool = False
    gsh_detox_penalty_ratio: float | None = None
    endpoint_toxic_flux_result: EndpointToxicFluxResult | None = None
    gsh_redox_capacity_result: GSHRedoxCapacityResult | None = None
    gsh_coupling: GSHBurdenCouplingInput | None = None
    evidence: EvidenceRecord | None = None
    metadata: JsonDict | None = None


@dataclass
class EffectiveBurdenResult(SerializableRecord):
    """Standalone relative effective carcinogenic burden result."""

    activation_burden_ratio: float
    detox_failure_ratio: float
    gsh_detox_penalty_ratio: float
    susceptibility_modifier: float
    effective_carcinogenic_burden_ratio: float
    gsh_relevant: bool
    gsh_consumption_load: float | None
    gsh_consumption_load_scaled: float | None
    gsh_fraction: float | None
    redox_capacity_ratio: float | None
    model_boundary: str
    warnings: list[AssumptionWarning]
    evidence: EvidenceRecord | None = None
    gsh_coupling_result: GSHBurdenCouplingResult | None = None
    metadata: JsonDict | None = None


def compute_effective_carcinogenic_burden(
    burden_input: EffectiveBurdenInput | None = None,
    *,
    activation_burden_ratio: float | None = None,
    detox_failure_ratio: float | None = None,
    susceptibility_modifier: float | SusceptibilityModifier | None = None,
    gsh_relevant: bool = False,
    gsh_detox_penalty_ratio: float | None = None,
    endpoint_toxic_flux_result: EndpointToxicFluxResult | None = None,
    gsh_redox_capacity_result: GSHRedoxCapacityResult | None = None,
    gsh_coupling: GSHBurdenCouplingInput | None = None,
    evidence: EvidenceRecord | None = None,
    metadata: JsonDict | None = None,
) -> EffectiveBurdenResult:
    """Combine precomputed mechanism ratios into one relative burden ratio."""

    request = burden_input or EffectiveBurdenInput(
        activation_burden_ratio=activation_burden_ratio,
        detox_failure_ratio=detox_failure_ratio,
        susceptibility_modifier=susceptibility_modifier,
        gsh_relevant=gsh_relevant,
        gsh_detox_penalty_ratio=gsh_detox_penalty_ratio,
        endpoint_toxic_flux_result=endpoint_toxic_flux_result,
        gsh_redox_capacity_result=gsh_redox_capacity_result,
        gsh_coupling=gsh_coupling,
        evidence=evidence,
        metadata=metadata,
    )
    warnings: list[AssumptionWarning] = []
    ratio_sources: JsonDict = {}

    endpoint_result = request.endpoint_toxic_flux_result
    if endpoint_result is not None:
        activation = _resolve_nonnegative_ratio(
            endpoint_result.activation_burden_ratio,
            1.0,
            "activation_burden_ratio",
            warnings,
            missing_is_warning=False,
        )
        detox = _resolve_nonnegative_ratio(
            endpoint_result.detox_failure_ratio,
            1.0,
            "detox_failure_ratio",
            warnings,
            missing_is_warning=False,
        )
        warnings.extend(endpoint_result.warnings)
        ratio_sources["activation_burden_ratio"] = "endpoint_toxic_flux_result"
        ratio_sources["detox_failure_ratio"] = "endpoint_toxic_flux_result"
    else:
        activation = _resolve_nonnegative_ratio(
            request.activation_burden_ratio,
            1.0,
            "activation_burden_ratio",
            warnings,
        )
        detox = _resolve_nonnegative_ratio(
            request.detox_failure_ratio,
            1.0,
            "detox_failure_ratio",
            warnings,
        )
        ratio_sources["activation_burden_ratio"] = "explicit_or_neutral_default"
        ratio_sources["detox_failure_ratio"] = "explicit_or_neutral_default"

    susceptibility = _resolve_susceptibility_modifier(request.susceptibility_modifier, warnings)
    gsh_penalty, gsh_context = _resolve_gsh_penalty(request, warnings)

    effective = activation * detox * gsh_penalty * susceptibility
    if not isfinite(effective):
        warnings.append(
            _warning(
                "effective_burden_invalid_bounded",
                "Effective burden ratio was invalid; using neutral bounded output.",
                field="effective_carcinogenic_burden_ratio",
            )
        )
        effective = 1.0
    if effective < 0.0:
        warnings.append(
            _warning(
                "effective_burden_negative_clamped",
                "Effective burden ratio was negative and was clamped to zero.",
                field="effective_carcinogenic_burden_ratio",
            )
        )
        effective = 0.0

    result_metadata: JsonDict = {
        "phase": "phase_8_effective_carcinogenic_burden",
        "model_family": "semi_mechanistic_relative_burden",
        "relative_burden_formula": (
            "ActivationBurdenRatio * DetoxFailureRatio * "
            "GSHDetoxPenaltyRatio * SusceptibilityModifier"
        ),
        "uses_endpoint_toxic_flux_result": endpoint_result is not None,
        "uses_gsh_redox_capacity_result": gsh_context["gsh_redox_capacity_result"] is not None,
        "ratio_sources": ratio_sources,
        "public_risk_output": "not_produced_or_modified",
        "engine_integration": False,
        "kinetic_factor_computation": False,
        "validated_pbpk_ode_model": False,
        "clinical_risk_model": False,
        "shapley_attribution": False,
        "phase9_behavior": False,
    }
    if request.metadata:
        result_metadata["input_metadata"] = dict(request.metadata)
    if metadata:
        result_metadata["caller_metadata"] = dict(metadata)

    return EffectiveBurdenResult(
        activation_burden_ratio=round(activation, 6),
        detox_failure_ratio=round(detox, 6),
        gsh_detox_penalty_ratio=round(gsh_penalty, 6),
        susceptibility_modifier=round(susceptibility, 6),
        effective_carcinogenic_burden_ratio=round(effective, 6),
        gsh_relevant=request.gsh_relevant,
        gsh_consumption_load=gsh_context["gsh_consumption_load"],
        gsh_consumption_load_scaled=gsh_context["gsh_consumption_load_scaled"],
        gsh_fraction=gsh_context["gsh_fraction"],
        redox_capacity_ratio=gsh_context["redox_capacity_ratio"],
        model_boundary="standalone_internal_semi_mechanistic_relative_burden_ratio",
        warnings=warnings,
        evidence=request.evidence or _default_evidence(),
        gsh_coupling_result=gsh_context["gsh_coupling_result"],
        metadata=result_metadata,
    )


def couple_gsh_consumption_to_activation_burden(
    coupling_input: GSHBurdenCouplingInput | None = None,
    *,
    gsh_relevant: bool = False,
    base_gsh_consumption_load: float | None = 0.0,
    upstream_activation_burden_ratio: float | None = None,
    d_factor: float | None = None,
    k_factor: float | None = None,
    tissue: str | None = None,
    synthesis_capacity: float | None = None,
    turnover_capacity: float | None = None,
    baseline_capacity: float | None = None,
    evidence: EvidenceRecord | None = None,
    metadata: JsonDict | None = None,
) -> GSHBurdenCouplingResult:
    """Scale GSH consumption by caller-supplied upstream activation burden."""

    request = coupling_input or GSHBurdenCouplingInput(
        gsh_relevant=gsh_relevant,
        base_gsh_consumption_load=base_gsh_consumption_load,
        upstream_activation_burden_ratio=upstream_activation_burden_ratio,
        d_factor=d_factor,
        k_factor=k_factor,
        tissue=tissue,
        synthesis_capacity=synthesis_capacity,
        turnover_capacity=turnover_capacity,
        baseline_capacity=baseline_capacity,
        evidence=evidence,
        metadata=metadata,
    )
    warnings: list[AssumptionWarning] = []
    base_load = _resolve_nonnegative_ratio(
        request.base_gsh_consumption_load,
        0.0,
        "base_gsh_consumption_load",
        warnings,
        missing_is_warning=True,
    )

    if not request.gsh_relevant:
        return GSHBurdenCouplingResult(
            gsh_relevant=False,
            base_gsh_consumption_load=round(base_load, 6),
            upstream_activation_burden_ratio=1.0,
            gsh_consumption_load_scaled=round(base_load, 6),
            scaling_source="not_gsh_relevant_neutral",
            d_factor=request.d_factor,
            k_factor=request.k_factor,
            gsh_redox_capacity_result=None,
            warnings=warnings,
            evidence=request.evidence or _default_evidence(),
            metadata=_gsh_coupling_metadata("not_gsh_relevant_neutral", request.metadata, metadata),
        )

    scaling_source = "neutral_fallback"
    upstream = 1.0
    if request.upstream_activation_burden_ratio is not None:
        upstream = _resolve_nonnegative_ratio(
            request.upstream_activation_burden_ratio,
            1.0,
            "upstream_activation_burden_ratio",
            warnings,
            missing_is_warning=False,
        )
        scaling_source = "explicit_upstream_activation_burden_ratio"
    elif request.d_factor is not None and request.k_factor is not None:
        d_value = _resolve_nonnegative_ratio(
            request.d_factor,
            1.0,
            "d_factor",
            warnings,
            missing_is_warning=False,
        )
        k_value = _resolve_nonnegative_ratio(
            request.k_factor,
            1.0,
            "k_factor",
            warnings,
            missing_is_warning=False,
        )
        upstream = d_value * k_value
        scaling_source = "d_times_k_approximation"
    else:
        warnings.append(
            _warning(
                "gsh_upstream_activation_missing_neutral",
                "GSH-relevant substrate lacked explicit upstream activation or supplied D/K factors; using neutral scaling.",
                field="upstream_activation_burden_ratio",
            )
        )

    scaled_load = base_load * upstream
    gsh_result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(
            tissue=request.tissue,
            consumption_load=scaled_load,
            synthesis_capacity=request.synthesis_capacity,
            turnover_capacity=request.turnover_capacity,
            baseline_capacity=request.baseline_capacity,
            evidence=request.evidence,
            metadata={
                "phase8_gsh_coupling": True,
                "scaling_source": scaling_source,
                "base_gsh_consumption_load": base_load,
                "upstream_activation_burden_ratio": upstream,
            },
        )
    )
    warnings.extend(gsh_result.warnings)

    return GSHBurdenCouplingResult(
        gsh_relevant=True,
        base_gsh_consumption_load=round(base_load, 6),
        upstream_activation_burden_ratio=round(upstream, 6),
        gsh_consumption_load_scaled=round(scaled_load, 6),
        scaling_source=scaling_source,
        d_factor=request.d_factor,
        k_factor=request.k_factor,
        gsh_redox_capacity_result=gsh_result,
        warnings=warnings,
        evidence=request.evidence or _default_evidence(),
        metadata=_gsh_coupling_metadata(scaling_source, request.metadata, metadata),
    )


def effective_burden_from_endpoint_and_gsh(
    endpoint_toxic_flux_result: EndpointToxicFluxResult,
    gsh_redox_capacity_result: GSHRedoxCapacityResult | None = None,
    *,
    susceptibility_modifier: float | SusceptibilityModifier | None = None,
    gsh_relevant: bool = False,
    gsh_coupling: GSHBurdenCouplingInput | None = None,
    metadata: JsonDict | None = None,
) -> EffectiveBurdenResult:
    """Build effective burden from accepted Phase 6 and optional Phase 7 outputs."""

    return compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            endpoint_toxic_flux_result=endpoint_toxic_flux_result,
            gsh_redox_capacity_result=gsh_redox_capacity_result,
            susceptibility_modifier=susceptibility_modifier,
            gsh_relevant=gsh_relevant,
            gsh_coupling=gsh_coupling,
            metadata=metadata,
        )
    )


def _resolve_gsh_penalty(
    request: EffectiveBurdenInput,
    warnings: list[AssumptionWarning],
) -> tuple[float, JsonDict]:
    context: JsonDict = {
        "gsh_consumption_load": None,
        "gsh_consumption_load_scaled": None,
        "gsh_fraction": None,
        "redox_capacity_ratio": None,
        "gsh_redox_capacity_result": None,
        "gsh_coupling_result": None,
    }

    if not request.gsh_relevant:
        return 1.0, context

    if request.gsh_coupling is not None:
        coupling_result = couple_gsh_consumption_to_activation_burden(request.gsh_coupling)
        context["gsh_coupling_result"] = coupling_result
        context["gsh_redox_capacity_result"] = coupling_result.gsh_redox_capacity_result
        context["gsh_consumption_load"] = coupling_result.base_gsh_consumption_load
        context["gsh_consumption_load_scaled"] = coupling_result.gsh_consumption_load_scaled
        context["gsh_fraction"] = coupling_result.gsh_fraction
        context["redox_capacity_ratio"] = coupling_result.redox_capacity_ratio
        warnings.extend(coupling_result.warnings)
        return coupling_result.detox_penalty_multiplier, context

    if request.gsh_redox_capacity_result is not None:
        result = request.gsh_redox_capacity_result
        context["gsh_redox_capacity_result"] = result
        context["gsh_consumption_load"] = result.consumption_load
        context["gsh_consumption_load_scaled"] = result.consumption_load
        context["gsh_fraction"] = result.gsh_fraction
        context["redox_capacity_ratio"] = result.redox_capacity_ratio
        warnings.extend(result.warnings)
        return _resolve_nonnegative_ratio(
            result.detox_penalty_multiplier,
            1.0,
            "gsh_detox_penalty_ratio",
            warnings,
            missing_is_warning=False,
        ), context

    if request.gsh_detox_penalty_ratio is not None:
        penalty = _resolve_nonnegative_ratio(
            request.gsh_detox_penalty_ratio,
            1.0,
            "gsh_detox_penalty_ratio",
            warnings,
            missing_is_warning=False,
        )
        return penalty, context

    warnings.append(
        _warning(
            "gsh_detox_penalty_missing_neutral",
            "GSH-relevant burden lacked GSH result, coupling input, or explicit penalty; using neutral GSH penalty.",
            field="gsh_detox_penalty_ratio",
        )
    )
    return 1.0, context


def _resolve_susceptibility_modifier(
    modifier: float | SusceptibilityModifier | None,
    warnings: list[AssumptionWarning],
) -> float:
    if isinstance(modifier, SusceptibilityModifier):
        warnings.extend(modifier.warnings or [])
        return _resolve_nonnegative_ratio(
            modifier.modifier_ratio,
            1.0,
            "susceptibility_modifier",
            warnings,
        )
    return _resolve_nonnegative_ratio(modifier, 1.0, "susceptibility_modifier", warnings)


def _resolve_nonnegative_ratio(
    value: float | None,
    default: float,
    field: str,
    warnings: list[AssumptionWarning],
    *,
    missing_is_warning: bool = True,
) -> float:
    if value is None:
        if missing_is_warning:
            warnings.append(
                _warning(
                    f"{field}_missing_neutral_default",
                    f"Missing {field}; using neutral default {default}.",
                    field=field,
                )
            )
        return float(default)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        warnings.append(
            _warning(
                f"{field}_invalid_neutral_default",
                f"Invalid {field}; using neutral default {default}.",
                field=field,
            )
        )
        return float(default)
    if not isfinite(numeric):
        warnings.append(
            _warning(
                f"{field}_invalid_neutral_default",
                f"Invalid {field}; using neutral default {default}.",
                field=field,
            )
        )
        return float(default)
    if numeric < 0.0:
        warnings.append(
            _warning(
                f"{field}_negative_clamped",
                f"Negative {field} was clamped to zero for bounded relative-burden output.",
                field=field,
            )
        )
        return 0.0
    return numeric


def _gsh_coupling_metadata(
    scaling_source: str,
    input_metadata: JsonDict | None,
    caller_metadata: JsonDict | None,
) -> JsonDict:
    metadata: JsonDict = {
        "phase": "phase_8_gsh_activation_coupling",
        "scaling_source": scaling_source,
        "gsh_relevant_neutral_when_false": True,
        "explicit_upstream_ratio_preferred": True,
        "d_times_k_supplied_factors_only": True,
        "internal_d_or_k_computation": False,
        "double_counts_endpoint_toxic_flux": False,
        "engine_integration": False,
    }
    if input_metadata:
        metadata["input_metadata"] = dict(input_metadata)
    if caller_metadata:
        metadata["caller_metadata"] = dict(caller_metadata)
    return metadata


def _default_evidence() -> EvidenceRecord:
    return EvidenceRecord(
        source="ExposoGraph Phase 8 local integration layer",
        grade=EvidenceGrade.PLACEHOLDER,
        confidence="local_2_0_semi_mechanistic_relative_burden",
        notes="Standalone internal relative burden integration; not a public risk-output change.",
        metadata={"phase": "phase_8"},
    )


def _warning(code: str, message: str, *, field: str | None = None) -> AssumptionWarning:
    return AssumptionWarning(
        code=code,
        message=message,
        field=field,
        review_status=SMEReviewStatus.UNKNOWN,
    )
