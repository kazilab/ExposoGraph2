"""Centralized kinetic-parameter resolver for Phase 5.

The resolver is a local-data-only seam for Ki lookup. It prefers curated Ki,
keeps Km proxy use visible and low confidence, and guards unavailable IC50 or
affinity paths without wiring anything into the interaction engine.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .interaction_schema import EvidenceGrade, ParameterUncertainty, ReleaseTarget, SerializableRecord
from .parameter_provider import JSONInteractionParameterProvider
from .parameter_resolution import (
    AffinityFallbackStatus,
    ParameterResolutionMethod,
    ParameterResolutionWarning,
    ParameterSourceKind,
    ResolvedParameter,
)


MODULE3_KM_STATIC_FOR_2_0 = True


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
    """Resolve Ki values from local provider data with explicit fallbacks."""

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
    """Resolve Ki through the centralized Phase 5 local-data seam."""

    return KineticParameterResolver(provider=provider).get_ki(
        enzyme,
        inhibitor,
        target_substrate=target_substrate,
        context=context,
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


def _normalize(value: object) -> str:
    return "".join(character.lower() for character in str(value) if character.isalnum())
