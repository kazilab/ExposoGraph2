"""Corrected/versioned GSH redox-capacity model for Phase 7.

This standalone module implements a quasi-steady, semi-mechanistic relative
capacity layer. It does not integrate with interaction_engine public outputs
or couple GSH state to Phase 8/9 mechanisms.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .interaction_schema import (
    AssumptionWarning,
    EvidenceGrade,
    EvidenceRecord,
    JsonDict,
    SerializableRecord,
    SMEReviewStatus,
    TissueContext,
    ValueEnum,
)


class GSHModelVersion(ValueEnum):
    """Named GSH model versions exposed by the Phase 7 standalone layer."""

    PHASE7_QUASI_STEADY_RELATIVE_CAPACITY = "phase7_quasi_steady_relative_capacity"
    LEGACY_DETACHED_GSH_PENALTY = "legacy_detached_gsh_penalty"


@dataclass
class GSHTissuePreset(SerializableRecord):
    """Local deterministic tissue preset for relative GSH capacity."""

    tissue: str
    relative_synthesis_capacity: float
    relative_turnover_capacity: float
    baseline_capacity: float = 1.0
    evidence: EvidenceRecord | None = None
    warnings: list[AssumptionWarning] | None = None
    metadata: JsonDict | None = None


@dataclass
class GSHRedoxCapacityInput(SerializableRecord):
    """Inputs for the standalone Phase 7 redox-capacity calculation."""

    tissue: str | None = None
    consumption_load: float | None = 0.0
    synthesis_capacity: float | None = None
    turnover_capacity: float | None = None
    baseline_capacity: float | None = None
    model_version: GSHModelVersion = GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY
    evidence: EvidenceRecord | None = None
    metadata: JsonDict | None = None


@dataclass
class GSHRedoxCapacityResult(SerializableRecord):
    """Standalone relative GSH redox-capacity result."""

    tissue: str
    model_version: GSHModelVersion
    gsh_fraction: float
    redox_capacity_ratio: float
    detox_penalty_multiplier: float
    synthesis_capacity: float
    consumption_load: float
    turnover_capacity: float
    baseline_capacity: float
    clamped: bool
    warnings: list[AssumptionWarning]
    evidence: EvidenceRecord | None = None
    tissue_context: TissueContext | None = None
    metadata: JsonDict | None = None


GSHRedoxWarning = AssumptionWarning


def get_default_gsh_tissue_presets() -> dict[str, GSHTissuePreset]:
    """Return local 2.0 presets for relative GSH capacity calculations."""

    return {
        "default": _preset(
            "default",
            synthesis=1.0,
            turnover=1.0,
            baseline=1.0,
            notes="Default 2.0 semi-mechanistic relative-capacity preset.",
        ),
        "liver": _preset(
            "liver",
            synthesis=1.25,
            turnover=1.1,
            baseline=1.0,
            notes="Liver 2.0 local preset with modestly higher relative replenishment.",
        ),
        "lung": _preset(
            "lung",
            synthesis=0.85,
            turnover=0.9,
            baseline=1.0,
            notes="Lung 2.0 local preset for relative-capacity comparisons.",
        ),
        "bone_marrow": _preset(
            "bone_marrow",
            synthesis=0.75,
            turnover=0.8,
            baseline=1.0,
            notes="Bone marrow 2.0 local preset for hematopoietic context comparisons.",
        ),
        "hematopoietic": _preset(
            "hematopoietic",
            synthesis=0.75,
            turnover=0.8,
            baseline=1.0,
            notes="Hematopoietic alias for the bone marrow 2.0 local preset.",
        ),
    }


def compute_quasi_steady_gsh_fraction(
    synthesis_capacity: float,
    consumption_load: float,
) -> float:
    """Compute bounded quasi-steady GSH fraction from synthesis and load."""

    synthesis = max(0.0, float(synthesis_capacity))
    consumption = max(0.0, float(consumption_load))
    if consumption == 0.0:
        return 1.0 if synthesis > 0.0 else 0.0
    if synthesis == 0.0:
        return 0.0
    return _clamp01(synthesis / (synthesis + consumption))


def compute_gsh_redox_capacity(
    gsh_input: GSHRedoxCapacityInput | None = None,
    *,
    tissue: str | None = None,
    consumption_load: float | None = None,
    synthesis_capacity: float | None = None,
    turnover_capacity: float | None = None,
    baseline_capacity: float | None = None,
    presets: dict[str, GSHTissuePreset] | None = None,
    metadata: JsonDict | None = None,
) -> GSHRedoxCapacityResult:
    """Compute Phase 7 quasi-steady relative GSH redox capacity."""

    request = gsh_input or GSHRedoxCapacityInput(
        tissue=tissue,
        consumption_load=consumption_load,
        synthesis_capacity=synthesis_capacity,
        turnover_capacity=turnover_capacity,
        baseline_capacity=baseline_capacity,
        metadata=metadata,
    )
    active_presets = presets or get_default_gsh_tissue_presets()
    preset, preset_warnings = _lookup_preset(request.tissue or tissue, active_presets)
    warnings = list(preset.warnings or []) + preset_warnings

    synthesis, synthesis_warned = _resolve_nonnegative(
        request.synthesis_capacity if request.synthesis_capacity is not None else synthesis_capacity,
        preset.relative_synthesis_capacity,
        "synthesis_capacity",
        warnings,
    )
    consumption, consumption_warned = _resolve_nonnegative(
        request.consumption_load if request.consumption_load is not None else consumption_load,
        0.0,
        "consumption_load",
        warnings,
    )
    turnover, turnover_warned = _resolve_nonnegative(
        request.turnover_capacity if request.turnover_capacity is not None else turnover_capacity,
        preset.relative_turnover_capacity,
        "turnover_capacity",
        warnings,
    )
    baseline, baseline_warned = _resolve_nonnegative(
        request.baseline_capacity if request.baseline_capacity is not None else baseline_capacity,
        preset.baseline_capacity,
        "baseline_capacity",
        warnings,
    )

    if baseline == 0.0:
        warnings.append(
            _warning(
                "baseline_capacity_zero",
                "Baseline capacity is zero; redox capacity ratio is bounded at zero.",
                field="baseline_capacity",
            )
        )

    raw_fraction = _raw_quasi_steady_fraction(synthesis, consumption)
    gsh_fraction = _clamp01(raw_fraction)
    redox_capacity_ratio = _clamp01(gsh_fraction * baseline)
    detox_penalty_multiplier = _detox_penalty(redox_capacity_ratio)
    clamped = (
        raw_fraction != gsh_fraction
        or synthesis_warned
        or consumption_warned
        or turnover_warned
        or baseline_warned
        or redox_capacity_ratio != gsh_fraction * baseline
    )

    result_metadata: JsonDict = {
        "phase": "phase_7_gsh_redox_capacity",
        "model_family": "semi_mechanistic_relative_capacity",
        "quasi_steady_expression": "synthesis_capacity / (synthesis_capacity + consumption_load)",
        "turnover_at_current_fraction": turnover * gsh_fraction,
        "tissue_preset": preset.tissue,
        "local_deterministic_2_0_preset": True,
        "validated_pbpk_ode_gsh_gssg_nrf2": False,
        "phase8_coupling": False,
        "shapley_attribution": False,
        "public_risk_output": "not_produced",
        "engine_integration": False,
    }
    if request.metadata:
        result_metadata["input_metadata"] = dict(request.metadata)
    if metadata:
        result_metadata["caller_metadata"] = dict(metadata)

    return GSHRedoxCapacityResult(
        tissue=preset.tissue,
        model_version=request.model_version,
        gsh_fraction=round(gsh_fraction, 6),
        redox_capacity_ratio=round(redox_capacity_ratio, 6),
        detox_penalty_multiplier=round(detox_penalty_multiplier, 6),
        synthesis_capacity=synthesis,
        consumption_load=consumption,
        turnover_capacity=turnover,
        baseline_capacity=baseline,
        clamped=clamped,
        warnings=warnings,
        evidence=request.evidence or preset.evidence,
        tissue_context=TissueContext(tissue=preset.tissue, metadata={"preset": preset.tissue}),
        metadata=result_metadata,
    )


def legacy_gsh_depletion_fraction(
    synthesis_capacity: float,
    consumption_load: float,
) -> float:
    """Explicit legacy detached GSH penalty comparator, not Phase 7 default."""

    synthesis = max(0.0, float(synthesis_capacity))
    consumption = max(0.0, float(consumption_load))
    if consumption == 0.0:
        return 1.0
    if synthesis == 0.0:
        return 0.0
    return _clamp01(1.0 - (consumption / synthesis))


def _preset(
    tissue: str,
    *,
    synthesis: float,
    turnover: float,
    baseline: float,
    notes: str,
) -> GSHTissuePreset:
    warning = _warning(
        "local_2_0_gsh_preset",
        "GSH tissue preset is a local 2.0 semi-mechanistic relative-capacity preset, not a curated tissue constant.",
        field="tissue",
    )
    return GSHTissuePreset(
        tissue=tissue,
        relative_synthesis_capacity=synthesis,
        relative_turnover_capacity=turnover,
        baseline_capacity=baseline,
        evidence=EvidenceRecord(
            source="ExposoGraph Phase 7 local preset",
            grade=EvidenceGrade.PLACEHOLDER,
            confidence="local_2_0_preset",
            notes=notes,
            metadata={"local_deterministic_preset": True},
        ),
        warnings=[warning],
        metadata={"model_boundary": "semi_mechanistic_relative_capacity"},
    )


def _lookup_preset(
    tissue: str | None,
    presets: dict[str, GSHTissuePreset],
) -> tuple[GSHTissuePreset, list[AssumptionWarning]]:
    if tissue is None:
        return presets["default"], []
    key = _normalize_key(tissue)
    alias = "bone_marrow" if key in {"bonemarrow", "marrow"} else key
    for preset_key, preset in presets.items():
        if _normalize_key(preset_key) == alias or _normalize_key(preset.tissue) == alias:
            return preset, []
    return presets["default"], [
        _warning(
            "unknown_tissue_default_preset",
            "Unknown tissue used default GSH redox-capacity preset instead of invented tissue constants.",
            field="tissue",
        )
    ]


def _resolve_nonnegative(
    value: float | None,
    default: float,
    field: str,
    warnings: list[AssumptionWarning],
) -> tuple[float, bool]:
    if value is None:
        return float(default), False
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        warnings.append(
            _warning(
                f"{field}_invalid",
                f"Invalid {field} value; using bounded default instead.",
                field=field,
            )
        )
        return float(default), True
    if numeric < 0.0:
        warnings.append(
            _warning(
                f"{field}_negative_clamped",
                f"Negative {field} was clamped to zero for bounded relative-capacity output.",
                field=field,
            )
        )
        return 0.0, True
    return numeric, False


def _raw_quasi_steady_fraction(synthesis_capacity: float, consumption_load: float) -> float:
    if consumption_load <= 0.0:
        return 1.0 if synthesis_capacity > 0.0 else 0.0
    if synthesis_capacity <= 0.0:
        return 0.0
    return synthesis_capacity / (synthesis_capacity + consumption_load)


def _detox_penalty(redox_capacity_ratio: float) -> float:
    if redox_capacity_ratio <= 0.0:
        return 20.0
    return min(20.0, 1.0 / redox_capacity_ratio)


def _clamp01(value: float) -> float:
    return min(1.0, max(0.0, float(value)))


def _warning(code: str, message: str, *, field: str | None = None) -> AssumptionWarning:
    return AssumptionWarning(
        code=code,
        message=message,
        field=field,
        review_status=SMEReviewStatus.UNKNOWN,
    )


def _normalize_key(value: object) -> str:
    return "".join(character.lower() for character in str(value) if character.isalnum())
