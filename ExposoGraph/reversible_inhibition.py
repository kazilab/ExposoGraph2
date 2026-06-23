"""Isolated reversible-inhibition equation kernel."""

from __future__ import annotations

from dataclasses import dataclass
from math import isclose, isfinite
from typing import Any

from .interaction_schema import AssumptionWarning, InhibitionMode, SerializableRecord, enum_from_value


EQUATION_ID = "reversible_inhibition.mixed.v1"
EQUATION_VERSION = 1


@dataclass
class ReversibleInhibitionResult(SerializableRecord):
    equation_id: str
    equation_version: int
    mode: InhibitionMode
    alpha_free_enzyme: float
    alpha_enzyme_substrate: float
    substrate_to_km_ratio: float
    flux_modifier: float
    baseline_rate: float
    inhibited_rate: float
    apparent_km: float
    apparent_vmax_fraction: float
    warnings: list[AssumptionWarning]


def compute_reversible_inhibition(
    *,
    mode: InhibitionMode | str,
    substrate_concentration: float,
    km: float,
    vmax: float,
    inhibitor_concentration: float,
    ki_free_enzyme: float | None = None,
    ki_enzyme_substrate: float | None = None,
) -> ReversibleInhibitionResult:
    """Compute the parent mixed reversible-inhibition equation.

    The binding constants define alpha values for the free-enzyme arm
    (``alpha_E``) and enzyme-substrate arm (``alpha_ES``). An absent arm is
    represented as ``None`` and contributes exactly ``1.0``.
    """

    active_mode = enum_from_value(InhibitionMode, mode, InhibitionMode.UNKNOWN)
    substrate = _finite_float(substrate_concentration, "substrate_concentration")
    km_value = _finite_float(km, "km")
    vmax_value = _finite_float(vmax, "vmax")
    inhibitor = _finite_float(inhibitor_concentration, "inhibitor_concentration")
    ki_e = _optional_finite_float(ki_free_enzyme, "ki_free_enzyme")
    ki_es = _optional_finite_float(ki_enzyme_substrate, "ki_enzyme_substrate")

    if active_mode is InhibitionMode.UNKNOWN:
        raise ValueError("unknown inhibition mode cannot be used for quantitative calculation")
    if substrate < 0.0:
        raise ValueError("substrate_concentration must be non-negative")
    if inhibitor < 0.0:
        raise ValueError("inhibitor_concentration must be non-negative")
    if km_value <= 0.0:
        raise ValueError("km must be greater than zero")
    if vmax_value < 0.0:
        raise ValueError("vmax must be non-negative")

    _validate_mode_contract(active_mode, ki_e, ki_es)

    alpha_e = _alpha(inhibitor, ki_e, "ki_free_enzyme")
    alpha_es = _alpha(inhibitor, ki_es, "ki_enzyme_substrate")
    baseline_rate = _baseline_rate(vmax_value, substrate, km_value)
    inhibited_rate = _inhibited_rate(vmax_value, substrate, km_value, alpha_e, alpha_es)
    denominator = alpha_e * km_value + alpha_es * substrate
    flux_modifier = (km_value + substrate) / denominator
    apparent_vmax_fraction = 1.0 / alpha_es
    apparent_km = alpha_e * km_value / alpha_es

    warnings: list[AssumptionWarning] = []
    if substrate == 0.0:
        warnings.append(
            AssumptionWarning(
                code="zero_substrate_absolute_rates_zero",
                message="Zero substrate gives zero absolute rates; normalized modifier uses the parent equation.",
                field="substrate_concentration",
            )
        )

    return ReversibleInhibitionResult(
        equation_id=EQUATION_ID,
        equation_version=EQUATION_VERSION,
        mode=active_mode,
        alpha_free_enzyme=_clean_float(alpha_e),
        alpha_enzyme_substrate=_clean_float(alpha_es),
        substrate_to_km_ratio=_clean_float(substrate / km_value),
        flux_modifier=_clean_float(flux_modifier),
        baseline_rate=_clean_float(baseline_rate),
        inhibited_rate=_clean_float(inhibited_rate),
        apparent_km=_clean_float(apparent_km),
        apparent_vmax_fraction=_clean_float(apparent_vmax_fraction),
        warnings=warnings,
    )


def _baseline_rate(vmax: float, substrate: float, km: float) -> float:
    return vmax * substrate / (km + substrate)


def _inhibited_rate(
    vmax: float,
    substrate: float,
    km: float,
    alpha_free_enzyme: float,
    alpha_enzyme_substrate: float,
) -> float:
    return vmax * substrate / (alpha_free_enzyme * km + alpha_enzyme_substrate * substrate)


def _validate_mode_contract(
    mode: InhibitionMode,
    ki_free_enzyme: float | None,
    ki_enzyme_substrate: float | None,
) -> None:
    if mode is InhibitionMode.COMPETITIVE:
        _require_active_ki(ki_free_enzyme, "ki_free_enzyme")
        if ki_enzyme_substrate is not None:
            raise ValueError("competitive inhibition requires ki_enzyme_substrate to be absent")
        return

    if mode is InhibitionMode.UNCOMPETITIVE:
        if ki_free_enzyme is not None:
            raise ValueError("uncompetitive inhibition requires ki_free_enzyme to be absent")
        _require_active_ki(ki_enzyme_substrate, "ki_enzyme_substrate")
        return

    if mode is InhibitionMode.PURE_NONCOMPETITIVE:
        _require_active_ki(ki_free_enzyme, "ki_free_enzyme")
        _require_active_ki(ki_enzyme_substrate, "ki_enzyme_substrate")
        if not isclose(float(ki_free_enzyme), float(ki_enzyme_substrate), rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("pure_noncompetitive inhibition requires equal Ki values")
        return

    if mode is InhibitionMode.MIXED:
        _require_active_ki(ki_free_enzyme, "ki_free_enzyme")
        _require_active_ki(ki_enzyme_substrate, "ki_enzyme_substrate")
        return

    raise ValueError("unsupported inhibition mode")


def _alpha(inhibitor_concentration: float, ki: float | None, field_name: str) -> float:
    if ki is None:
        return 1.0
    _require_active_ki(ki, field_name)
    return 1.0 + inhibitor_concentration / ki


def _require_active_ki(value: float | None, field_name: str) -> None:
    if value is None:
        raise ValueError(f"{field_name} is required for this inhibition mode")
    if value <= 0.0:
        raise ValueError(f"{field_name} must be greater than zero when the binding arm is active")


def _finite_float(value: float, field_name: str) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be a finite numeric value") from exc
    if not isfinite(numeric):
        raise ValueError(f"{field_name} must be a finite numeric value")
    return numeric


def _optional_finite_float(value: float | None, field_name: str) -> float | None:
    if value is None:
        return None
    return _finite_float(value, field_name)


def _clean_float(value: float) -> float:
    if not isfinite(value):
        raise ValueError("reversible inhibition result contains a non-finite value")
    if value == 0.0:
        return 0.0
    return float(value)
