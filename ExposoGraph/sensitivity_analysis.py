"""Global sensitivity analysis for ExposoGraph 2.0.

Provides a self-contained Sobol-index estimator (no scipy dependency) for use
in place of the one-at-a-time ±20% perturbations described in Manuscript v6
Supplementary Figure S2:

    "future versions will implement global sensitivity analysis (Sobol
    indices) to quantify parameter interactions and identify the subset of
    parameters to which risk scores are most sensitive".

The estimator follows the Saltelli 2010 scheme, drawing two independent
quasi-pseudo-random sample matrices ``A`` and ``B`` of size ``N x k`` and
``k`` hybrid matrices ``A_B[i]`` in which column ``i`` of ``A`` is replaced
with column ``i`` of ``B``. For each parameter it returns:

* ``S1`` - the first-order Sobol index, the fraction of output variance
  explained by parameter ``i`` alone.
* ``ST`` - the total-order Sobol index, the fraction of variance that would
  remain if every parameter other than ``i`` were fixed.

Dependencies: only ``math`` and ``random`` from the Python stdlib.

A thin convenience wrapper :func:`sobol_for_synergy_score` evaluates the
indices against the interaction-engine synergy output for a supplied
carcinogen mixture.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Any, Callable, Iterable

ParamBounds = dict[str, tuple[float, float]]


@dataclass
class SobolIndex:
    """First- and total-order Sobol indices for a single parameter."""

    name: str
    S1: float
    ST: float
    S1_lower: float | None = None
    S1_upper: float | None = None
    ST_lower: float | None = None
    ST_upper: float | None = None


@dataclass
class SobolResult:
    """Full result of a Sobol sensitivity analysis."""

    indices: list[SobolIndex]
    n_samples: int
    output_mean: float
    output_variance: float
    total_model_calls: int


# ── Sampling helpers ──────────────────────────────────────────────────────


def _sample_matrix(
    rng: random.Random,
    bounds: ParamBounds,
    n: int,
) -> list[dict[str, float]]:
    """Return ``n`` rows of uniform samples drawn within ``bounds``."""
    matrix: list[dict[str, float]] = []
    for _ in range(n):
        row: dict[str, float] = {}
        for name, (lo, hi) in bounds.items():
            row[name] = rng.uniform(float(lo), float(hi))
        matrix.append(row)
    return matrix


def _evaluate(
    model_fn: Callable[[dict[str, float]], float],
    samples: Iterable[dict[str, float]],
) -> list[float]:
    return [float(model_fn(row)) for row in samples]


def _variance(values: list[float]) -> float:
    if not values:
        return 0.0
    mean = sum(values) / len(values)
    return sum((v - mean) ** 2 for v in values) / len(values)


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _bootstrap_indices(
    A_values: list[float],
    B_values: list[float],
    Ab_values_by_param: dict[str, list[float]],
    *,
    n_resamples: int,
    rng: random.Random,
) -> dict[str, tuple[float, float, float, float]]:
    """Return bootstrap 95% CI for S1 and ST per parameter."""
    n = len(A_values)
    cis: dict[str, tuple[float, float, float, float]] = {}

    for name, ab_values in Ab_values_by_param.items():
        s1_samples: list[float] = []
        st_samples: list[float] = []
        for _ in range(n_resamples):
            idx = [rng.randrange(n) for _ in range(n)]
            A_r = [A_values[i] for i in idx]
            B_r = [B_values[i] for i in idx]
            Ab_r = [ab_values[i] for i in idx]
            var_y = _variance(A_r + B_r)
            if var_y <= 0:
                s1_samples.append(0.0)
                st_samples.append(0.0)
                continue
            s1 = sum(B_r[i] * (Ab_r[i] - A_r[i]) for i in range(n)) / n / var_y
            st = sum((A_r[i] - Ab_r[i]) ** 2 for i in range(n)) / (2 * n) / var_y
            s1_samples.append(s1)
            st_samples.append(st)
        s1_samples.sort()
        st_samples.sort()

        def _percentile(samples: list[float], pct: float) -> float:
            if not samples:
                return 0.0
            idx = max(0, min(len(samples) - 1, int(round(pct * (len(samples) - 1)))))
            return samples[idx]

        cis[name] = (
            _percentile(s1_samples, 0.025),
            _percentile(s1_samples, 0.975),
            _percentile(st_samples, 0.025),
            _percentile(st_samples, 0.975),
        )
    return cis


# ── Core estimator ────────────────────────────────────────────────────────


def sobol_indices(
    model_fn: Callable[[dict[str, float]], float],
    bounds: ParamBounds,
    *,
    n_samples: int = 256,
    seed: int = 0,
    bootstrap_resamples: int = 0,
) -> SobolResult:
    """Estimate first- and total-order Sobol indices for ``model_fn``.

    Parameters
    ----------
    model_fn:
        Callable mapping ``{param_name: value}`` -> scalar output.
    bounds:
        Uniform sampling bounds per parameter as ``{name: (low, high)}``.
    n_samples:
        Base sample size ``N`` for the ``A`` and ``B`` matrices. Total model
        calls are ``N * (k + 2)`` where ``k`` is the number of parameters.
    seed:
        Seed for the stdlib RNG used to draw samples.
    bootstrap_resamples:
        If > 0, compute bootstrap 95% CIs for each Sobol index.
    """
    if n_samples < 4:
        raise ValueError("n_samples must be >= 4 for stable variance estimation")
    if not bounds:
        raise ValueError("bounds must contain at least one parameter")

    rng = random.Random(seed)
    param_names = list(bounds.keys())

    A = _sample_matrix(rng, bounds, n_samples)
    B = _sample_matrix(rng, bounds, n_samples)

    A_values = _evaluate(model_fn, A)
    B_values = _evaluate(model_fn, B)

    Ab_values_by_param: dict[str, list[float]] = {}
    total_calls = 2 * n_samples

    for name in param_names:
        Ab = [dict(row) for row in A]
        for i, row in enumerate(Ab):
            row[name] = B[i][name]
        Ab_values = _evaluate(model_fn, Ab)
        Ab_values_by_param[name] = Ab_values
        total_calls += n_samples

    full_output = A_values + B_values
    var_y = _variance(full_output)
    mean_y = _mean(full_output)

    ci_map: dict[str, tuple[float, float, float, float]] = {}
    if bootstrap_resamples > 0:
        ci_rng = random.Random(seed + 1)
        ci_map = _bootstrap_indices(
            A_values,
            B_values,
            Ab_values_by_param,
            n_resamples=bootstrap_resamples,
            rng=ci_rng,
        )

    indices: list[SobolIndex] = []
    for name in param_names:
        ab_values = Ab_values_by_param[name]
        if var_y <= 0:
            s1 = 0.0
            st = 0.0
        else:
            s1 = sum(
                B_values[i] * (ab_values[i] - A_values[i])
                for i in range(n_samples)
            ) / n_samples / var_y
            st = sum(
                (A_values[i] - ab_values[i]) ** 2
                for i in range(n_samples)
            ) / (2 * n_samples) / var_y

        ci = ci_map.get(name)
        indices.append(
            SobolIndex(
                name=name,
                S1=round(s1, 4),
                ST=round(st, 4),
                S1_lower=round(ci[0], 4) if ci else None,
                S1_upper=round(ci[1], 4) if ci else None,
                ST_lower=round(ci[2], 4) if ci else None,
                ST_upper=round(ci[3], 4) if ci else None,
            )
        )

    indices.sort(key=lambda entry: -entry.ST)

    return SobolResult(
        indices=indices,
        n_samples=n_samples,
        output_mean=round(mean_y, 4),
        output_variance=round(var_y, 6),
        total_model_calls=total_calls,
    )


# ── Convenience wrapper for interaction-engine synergy scores ─────────────


def sobol_for_synergy_score(
    mixture: dict[str, float],
    *,
    km_range: tuple[float, float] = (0.5, 1.5),
    expression_range: tuple[float, float] = (0.7, 1.3),
    n_samples: int = 128,
    seed: int = 0,
    carcinogen_pair: tuple[str, str] | None = None,
    tissue: str = "Liver",
    lifestyle: dict[str, Any] | None = None,
) -> SobolResult:
    """Run Sobol sensitivity on the interaction-engine synergy score.

    Parameters perturbed: per-substrate Km scale factors and per-enzyme
    expression scale factors passed into
    :func:`ExposoGraph.interaction_engine.compute_interaction_matrix` via its
    ``param_perturbations`` / ``expression_perturbations`` hooks.

    Returns Sobol indices for the scalar synergy score of a chosen carcinogen
    pair; defaults to the first pair yielded by the matrix.
    """
    from .interaction_engine import compute_interaction_matrix

    substrates = list(mixture.keys())
    enzymes = ("CYP1A1", "CYP1A2", "CYP2E1", "CYP3A4")

    bounds: ParamBounds = {}
    for name in substrates:
        bounds[f"km::{name}"] = km_range
    for name in enzymes:
        bounds[f"expr::{name}"] = expression_range

    def _model(params: dict[str, float]) -> float:
        param_perturbations: dict[str, dict[str, float]] = {}
        expression_perturbations: dict[str, float] = {}
        for key, value in params.items():
            prefix, sep, payload = key.partition("::")
            if not sep:
                continue
            if prefix == "km":
                param_perturbations[payload] = {"km_scale": float(value)}
            elif prefix == "expr":
                expression_perturbations[payload] = float(value)
        try:
            result = compute_interaction_matrix(
                mixture,
                tissue=tissue,
                lifestyle=lifestyle,
                param_perturbations=param_perturbations or None,
                expression_perturbations=expression_perturbations or None,
            )
        except Exception:
            return 1.0
        if not result.synergy_matrix:
            return 1.0
        if carcinogen_pair is not None:
            left, right = carcinogen_pair
            pair_key = f"{left}_x_{right}"
            if pair_key in result.synergy_matrix:
                return float(result.synergy_matrix[pair_key])
            reversed_key = f"{right}_x_{left}"
            if reversed_key in result.synergy_matrix:
                return float(result.synergy_matrix[reversed_key])
        first_value = next(iter(result.synergy_matrix.values()))
        return float(first_value)

    return sobol_indices(_model, bounds, n_samples=n_samples, seed=seed)


__all__ = [
    "ParamBounds",
    "SobolIndex",
    "SobolResult",
    "sobol_for_synergy_score",
    "sobol_indices",
]
