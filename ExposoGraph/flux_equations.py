"""Pure equation helpers for Module 3 flux calculations.

These helpers keep the paper-facing kinetic formulas separate from the larger
workflow engine. They intentionally contain no data loading, file IO, or public
workflow policy.
"""

from __future__ import annotations

import math


def michaelis_menten(S: float, Vmax: float, Km: float) -> float:
    """Calculate reaction velocity using Michaelis-Menten kinetics.

    v = Vmax * S / (Km + S)
    """
    if Km <= 0:
        raise ValueError(f"Km must be positive, got {Km}")
    if S < 0:
        raise ValueError(f"Substrate concentration cannot be negative, got {S}")
    return Vmax * S / (Km + S)


def hill_equation(S: float, Vmax: float, K50: float, n: float) -> float:
    """Calculate reaction velocity using Hill cooperative kinetics.

    v = Vmax * S^n / (K50^n + S^n)
    """
    if K50 <= 0:
        raise ValueError(f"K50 must be positive, got {K50}")
    if n <= 0:
        raise ValueError(f"Hill coefficient n must be positive, got {n}")
    if S < 0:
        raise ValueError(f"Substrate concentration cannot be negative, got {S}")
    Sn = S ** n
    K50n = K50 ** n
    return float(Vmax * Sn / (K50n + Sn))


def finite_nonnegative(value: float, field_name: str) -> float:
    """Normalize a numeric equation input while preserving legacy errors."""
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be finite and non-negative") from exc
    if not math.isfinite(numeric) or numeric < 0.0:
        raise ValueError(f"{field_name} must be finite and non-negative")
    if numeric == 0.0:
        return 0.0
    return numeric


def scaled_vmax(
    base_vmax: float,
    genotype_modifier: float,
    tissue_weight: float,
    *,
    vmax_relative: float = 1.0,
    relative_capacity_scale: float = 1.0,
) -> float:
    """Apply relative capacity, genotype, and tissue scaling to Vmax."""
    return (
        base_vmax
        * vmax_relative
        * relative_capacity_scale
        * genotype_modifier
        * tissue_weight
    )


def saturating_flux(S: float, scale: float, Km: float) -> float:
    """Evaluate the simple saturating proxy form scale * S / (Km + S)."""
    return scale * S / (Km + S)


def activation_detox_ratio(total_activation: float, total_detox: float) -> float:
    """Return the activation/detoxification ratio with legacy edge behavior."""
    activation = float(total_activation)
    detox = float(total_detox)
    if detox > 0:
        return activation / detox
    if activation > 0:
        return 999.0
    return 1.0


def susceptibility_score_log2(net_ratio: float) -> float:
    """Return the log2 susceptibility score used by Module 3 outputs."""
    if net_ratio <= 0 or not math.isfinite(net_ratio):
        return 0.0
    return round(math.log2(net_ratio), 4)
