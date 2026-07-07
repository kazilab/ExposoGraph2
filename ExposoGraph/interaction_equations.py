"""Pure equation helpers for Module 5 interaction-risk calculations.

The interaction engine still owns public workflows and data assembly. This
module names the small equations that are useful for review, testing, and paper
methods text.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

RoundFn = Callable[[float, int], float]


def explicit_dk_activation_scale(
    base_key: str,
    rate_key: str,
    exposure_profile: Mapping[str, Any],
) -> tuple[float | None, dict[str, float] | None]:
    """Return the explicit D x K activation scale when both factors exist."""
    for key in (base_key, rate_key):
        value = exposure_profile.get(key)
        if not isinstance(value, Mapping):
            continue
        if "d_factor" not in value or "k_factor" not in value:
            continue
        try:
            d_factor = max(0.0, float(value["d_factor"]))
            k_factor = max(0.0, float(value["k_factor"]))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(d_factor) or not math.isfinite(k_factor):
            continue
        return d_factor * k_factor, {"d_factor": d_factor, "k_factor": k_factor}
    return None, None


def gsh_detox_components(
    carcinogen: str,
    gsh_fraction: float,
    genotypes: Mapping[str, str],
    gsh_detox_map: Mapping[str, str | None],
    *,
    redox_detox_penalty: float | None = None,
    round_fn: RoundFn = round,
) -> tuple[float, float, float]:
    """Return susceptibility, GSH-pool penalty, and combined GSH multiplier."""
    if gsh_detox_map.get(carcinogen) is None:
        return 1.0, 1.0, 1.0

    genotype_factor = 1.0
    gstm1 = str(genotypes.get("GSTM1", "active")).lower()
    gstp1 = genotypes.get("GSTP1", "Ile105Ile")

    if carcinogen == "PAH" and gstm1 in {"null", "null/null", "deletion", "0"}:
        genotype_factor = 2.5
    elif carcinogen == "PAH" and gstp1 == "Val105Val":
        genotype_factor = 1.5

    if redox_detox_penalty is not None:
        gsh_penalty = max(1.0, float(redox_detox_penalty))
    else:
        if gsh_fraction >= 0.30:
            gsh_penalty = 1.0
        elif gsh_fraction >= 0.20:
            gsh_penalty = 1.0 + (0.30 - gsh_fraction) / 0.10 * 0.5
        else:
            gsh_penalty = 1.5 + (0.20 - gsh_fraction) / 0.20 * 3.0

    rounded_gsh_penalty = round_fn(gsh_penalty, 3)
    combined_penalty = round_fn(genotype_factor * gsh_penalty, 3)
    return genotype_factor, rounded_gsh_penalty, combined_penalty


def gsh_upstream_activation_scale(
    *,
    direct_activation_ratio: float | None,
    direct_status: str | None,
    direct_review_required: bool | None,
    explicit_dk_scale: float | None,
    explicit_dk_details: Mapping[str, float] | None,
    round_fn: RoundFn = round,
) -> tuple[float, str, dict[str, Any], str | None]:
    """Resolve the GSH upstream activation scale and fallback warning code."""
    direct_details: dict[str, Any] = {}
    direct_ratio = float("nan") if direct_activation_ratio is None else float(direct_activation_ratio)
    if direct_status is not None:
        direct_details = {
            "direct_activation_available": False,
            "inhibition_status": direct_status,
            "review_required": bool(direct_review_required),
            "activation_burden_ratio": (
                round_fn(direct_ratio, 6) if math.isfinite(direct_ratio) else None
            ),
        }

    has_direct_activation = (
        direct_status == "mechanism_resolved"
        and not direct_review_required
        and math.isfinite(direct_ratio)
        and direct_ratio >= 0.0
        and not math.isclose(direct_ratio, 1.0)
    )
    if has_direct_activation:
        return (
            max(0.0, direct_ratio),
            "direct_activation_burden_ratio",
            {
                "activation_burden_ratio": round_fn(direct_ratio, 6),
                "inhibition_status": direct_status,
            },
            None,
        )

    if explicit_dk_scale is not None:
        return (
            explicit_dk_scale,
            "d_times_k_approximation",
            {
                **direct_details,
                **dict(explicit_dk_details or {}),
                "internal_d_or_k_computation": False,
            },
            "gsh_d_times_k_fallback_used",
        )

    return (
        1.0,
        "neutral_missing_upstream_activation",
        {**direct_details, "upstream_activation_burden_ratio": 1.0},
        "gsh_upstream_activation_missing_neutral",
    )


def final_mechanism_multiplier(
    induction_multiplier: float,
    inhibition_burden_multiplier: float,
    gsh_penalty: float,
    *,
    round_fn: RoundFn = round,
) -> float:
    """Combine induction, inhibition burden, and GSH terms once."""
    return round_fn(induction_multiplier * inhibition_burden_multiplier * gsh_penalty, 6)


def adjusted_relative_risk(
    baseline_relative_risk: float,
    final_multiplier: float,
    *,
    round_fn: RoundFn = round,
) -> float:
    """Apply the final mechanism multiplier to baseline relative risk."""
    return round_fn(baseline_relative_risk * final_multiplier, 3)


def pairwise_synergy_factor(
    left_individual_risk: float,
    right_individual_risk: float,
    left_adjusted_risk: float,
    right_adjusted_risk: float,
    *,
    round_fn: RoundFn = round,
) -> float:
    """Return adjusted-pair risk divided by independent-pair risk."""
    independent_total = left_individual_risk + right_individual_risk
    adjusted_total = left_adjusted_risk + right_adjusted_risk
    return round_fn(adjusted_total / independent_total if independent_total > 0 else 1.0, 3)
