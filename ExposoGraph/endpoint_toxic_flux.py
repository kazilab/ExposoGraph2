"""Endpoint-specific toxic-flux interpretation for Phase 6.

This module accepts precomputed flux ratios and interprets them through the
accepted Phase 4 reaction-role semantics. It does not calculate kinetic
parameters, call competitive-flux helpers, or integrate with public risk-engine
outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .interaction_schema import (
    AssumptionWarning,
    JsonDict,
    ReactionRole,
    RiskDirectionIfFluxDecreases,
    RiskEndpoint,
    SMEReviewNote,
    SMEReviewStatus,
    SerializableRecord,
)
from .reaction_role_semantics import (
    ReactionRoleAnnotation,
    ReactionRoleRegistry,
    get_default_reaction_role_registry,
    interpret_competitive_flux_ratio,
)


ToxicFluxWarning = AssumptionWarning


@dataclass
class EndpointFluxInput(SerializableRecord):
    """Precomputed competitive-flux fact plus endpoint context."""

    enzyme: str
    substrate: str
    flux_ratio: float
    tissue: str | None = None
    endpoint: str | RiskEndpoint | None = None
    annotation: ReactionRoleAnnotation | None = None
    metadata: JsonDict | None = None


@dataclass
class EndpointToxicFluxResult(SerializableRecord):
    """Standalone endpoint burden interpretation for a flux ratio."""

    enzyme: str
    substrate: str
    tissue: str | None
    endpoint: str | None
    flux_ratio: float
    endpoint_toxic_flux_ratio: float
    burden_multiplier: float
    activation_burden_ratio: float
    detox_failure_ratio: float
    reaction_role: ReactionRole
    risk_direction_if_flux_decreases: RiskDirectionIfFluxDecreases
    annotation_record_id: str | None
    warnings: list[AssumptionWarning]
    sme_notes: list[SMEReviewNote]
    metadata: JsonDict | None = None


def interpret_endpoint_toxic_flux(
    endpoint_input: EndpointFluxInput,
    *,
    registry: ReactionRoleRegistry | None = None,
) -> EndpointToxicFluxResult:
    """Interpret a precomputed flux ratio in endpoint-specific burden terms."""

    annotation = endpoint_input.annotation
    if annotation is None:
        active_registry = registry or get_default_reaction_role_registry()
        annotation = active_registry.lookup(
            endpoint_input.enzyme,
            endpoint_input.substrate,
            tissue=endpoint_input.tissue,
            endpoint=endpoint_input.endpoint,
        )

    return interpret_competitive_endpoint_flux(
        endpoint_input.flux_ratio,
        annotation,
        enzyme=endpoint_input.enzyme,
        substrate=endpoint_input.substrate,
        tissue=endpoint_input.tissue,
        endpoint=endpoint_input.endpoint,
        metadata=endpoint_input.metadata,
    )


def interpret_competitive_endpoint_flux(
    flux_ratio: float,
    annotation: ReactionRoleAnnotation,
    *,
    enzyme: str | None = None,
    substrate: str | None = None,
    tissue: str | None = None,
    endpoint: str | RiskEndpoint | None = None,
    metadata: JsonDict | None = None,
) -> EndpointToxicFluxResult:
    """Map a competitive flux ratio to endpoint toxic-burden channels."""

    if flux_ratio <= 0:
        raise ValueError("flux_ratio must be greater than zero")

    phase4_result = interpret_competitive_flux_ratio(flux_ratio, annotation)
    warnings = list(phase4_result.warnings)
    sme_notes = list(phase4_result.sme_notes)
    channel = _interpretation_channel(annotation)

    if channel == "neutral":
        burden_multiplier = 1.0
        activation_burden_ratio = 1.0
        detox_failure_ratio = 1.0
        _ensure_neutral_warning(annotation, warnings)
    elif channel == "detox_failure":
        burden_multiplier = 1.0 / flux_ratio
        activation_burden_ratio = 1.0
        detox_failure_ratio = burden_multiplier
    else:
        burden_multiplier = flux_ratio
        activation_burden_ratio = burden_multiplier
        detox_failure_ratio = 1.0

    result_metadata = {
        "phase": "phase_6_endpoint_specific_toxic_flux",
        "raw_flux_ratio": flux_ratio,
        "phase4_burden_multiplier": phase4_result.burden_multiplier,
        "interpretation_channel": channel,
        "precomputed_flux_ratio_required": True,
        "engine_integration": False,
        "public_risk_output": "not_produced",
        "reaction_role_metadata": dict(annotation.metadata or {}),
    }
    if metadata:
        result_metadata["caller_metadata"] = dict(metadata)

    return EndpointToxicFluxResult(
        enzyme=enzyme or annotation.enzyme,
        substrate=substrate or annotation.substrate,
        tissue=tissue,
        endpoint=_endpoint_value(endpoint),
        flux_ratio=flux_ratio,
        endpoint_toxic_flux_ratio=burden_multiplier,
        burden_multiplier=burden_multiplier,
        activation_burden_ratio=activation_burden_ratio,
        detox_failure_ratio=detox_failure_ratio,
        reaction_role=annotation.reaction_role,
        risk_direction_if_flux_decreases=annotation.risk_direction_if_flux_decreases,
        annotation_record_id=annotation.record_id,
        warnings=warnings,
        sme_notes=sme_notes,
        metadata=result_metadata,
    )


def endpoint_toxic_flux_from_registry_lookup(
    enzyme: str,
    substrate: str,
    flux_ratio: float,
    *,
    tissue: str | None = None,
    endpoint: str | RiskEndpoint | None = None,
    registry: ReactionRoleRegistry | None = None,
    metadata: JsonDict | None = None,
) -> EndpointToxicFluxResult:
    """Lookup a Phase 4 annotation and interpret a precomputed flux ratio."""

    active_registry = registry or get_default_reaction_role_registry()
    annotation = active_registry.lookup(enzyme, substrate, tissue=tissue, endpoint=endpoint)
    return interpret_endpoint_toxic_flux(
        EndpointFluxInput(
            enzyme=enzyme,
            substrate=substrate,
            flux_ratio=flux_ratio,
            tissue=tissue,
            endpoint=endpoint,
            annotation=annotation,
            metadata=metadata,
        )
    )


def _interpretation_channel(annotation: ReactionRoleAnnotation) -> str:
    if not annotation.active:
        return "neutral"
    if annotation.review_status in {
        SMEReviewStatus.CANDIDATE,
        SMEReviewStatus.PENDING_TEAM_AGREEMENT,
        SMEReviewStatus.DEFERRED_3_0,
    }:
        return "neutral"
    if annotation.reaction_role in {ReactionRole.UNKNOWN, ReactionRole.PROBE_ONLY, ReactionRole.DUAL_ROLE}:
        return "neutral"
    if annotation.risk_direction_if_flux_decreases in {
        RiskDirectionIfFluxDecreases.UNKNOWN,
        RiskDirectionIfFluxDecreases.NEUTRAL,
        RiskDirectionIfFluxDecreases.MIXED,
    }:
        return "neutral"
    if annotation.reaction_role in {ReactionRole.DETOXIFICATION, ReactionRole.CLEARANCE}:
        return "detox_failure"
    if annotation.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.INCREASE:
        return "detox_failure"
    if annotation.reaction_role is ReactionRole.BIOACTIVATION:
        return "activation"
    if annotation.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.DECREASE:
        return "activation"
    return "neutral"


def _ensure_neutral_warning(
    annotation: ReactionRoleAnnotation,
    warnings: list[AssumptionWarning],
) -> None:
    codes = {warning.code for warning in warnings}
    if not annotation.active and "endpoint_role_inactive" not in codes:
        warnings.append(
            _warning(
                "endpoint_role_inactive",
                "Endpoint toxic-flux interpretation is neutral because the annotation is inactive.",
            )
        )
    if annotation.review_status in {SMEReviewStatus.CANDIDATE, SMEReviewStatus.PENDING_TEAM_AGREEMENT}:
        if "endpoint_role_pending" not in codes:
            warnings.append(
                _warning(
                    "endpoint_role_pending",
                    "Endpoint toxic-flux interpretation is neutral because the annotation is not accepted curated behavior.",
                )
            )
    if annotation.review_status is SMEReviewStatus.DEFERRED_3_0 and "endpoint_role_deferred" not in codes:
        warnings.append(
            _warning(
                "endpoint_role_deferred",
                "Endpoint toxic-flux interpretation is neutral because the annotation is deferred.",
            )
        )
    if annotation.reaction_role is ReactionRole.DUAL_ROLE and "endpoint_dual_role_neutral" not in codes:
        warnings.append(
            _warning(
                "endpoint_dual_role_neutral",
                "Dual-role reactions require explicit context-specific resolution; using neutral burden.",
            )
        )
    if annotation.reaction_role is ReactionRole.PROBE_ONLY and "endpoint_probe_only_neutral" not in codes:
        warnings.append(
            _warning(
                "endpoint_probe_only_neutral",
                "Probe-only reactions are not interpreted as endpoint toxic burden.",
            )
        )
    if annotation.reaction_role is ReactionRole.UNKNOWN and "endpoint_role_unknown" not in codes:
        warnings.append(
            _warning(
                "endpoint_role_unknown",
                "Endpoint toxic-flux interpretation is neutral because the reaction role is unknown.",
            )
        )
    if annotation.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.UNKNOWN:
        if "endpoint_direction_unknown" not in codes:
            warnings.append(
                _warning(
                    "endpoint_direction_unknown",
                    "Endpoint toxic-flux interpretation is neutral because the risk direction is unknown.",
                )
            )
    if annotation.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.MIXED:
        if "endpoint_direction_mixed" not in codes:
            warnings.append(
                _warning(
                    "endpoint_direction_mixed",
                    "Mixed endpoint direction requires context-specific resolution; using neutral burden.",
                )
            )
    if annotation.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.NEUTRAL:
        if "endpoint_direction_neutral" not in codes:
            warnings.append(
                _warning(
                    "endpoint_direction_neutral",
                    "Neutral endpoint direction leaves toxic burden unchanged.",
                )
            )


def _warning(code: str, message: str) -> AssumptionWarning:
    return AssumptionWarning(
        code=code,
        message=message,
        field="endpoint_toxic_flux_ratio",
        review_status=SMEReviewStatus.UNKNOWN,
    )


def _endpoint_value(endpoint: str | RiskEndpoint | None) -> str | None:
    if isinstance(endpoint, RiskEndpoint):
        return endpoint.value
    return endpoint
