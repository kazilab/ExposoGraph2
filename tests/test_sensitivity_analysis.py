"""Smoke tests for the global sensitivity (Sobol indices) estimator."""

import pytest

from ExposoGraph import (
    SobolIndex,
    SobolResult,
    sobol_for_synergy_score,
    sobol_indices,
)


def test_sobol_indices_linear_single_parameter_is_fully_dominant():
    # y = a  ->  S1[a] ~= 1, ST[a] ~= 1
    result = sobol_indices(
        lambda params: params["a"],
        {"a": (0.0, 1.0)},
        n_samples=128,
        seed=42,
    )
    assert isinstance(result, SobolResult)
    assert len(result.indices) == 1
    idx = result.indices[0]
    assert idx.name == "a"
    assert idx.S1 == pytest.approx(1.0, abs=0.1)
    assert idx.ST == pytest.approx(1.0, abs=0.1)


def test_sobol_indices_identify_dominant_additive_parameter():
    # y = 10*a + b with a, b ~ U(0,1) -> a carries >>90% of variance
    result = sobol_indices(
        lambda params: 10.0 * params["a"] + params["b"],
        {"a": (0.0, 1.0), "b": (0.0, 1.0)},
        n_samples=256,
        seed=0,
    )
    ranked = {idx.name: idx for idx in result.indices}
    assert ranked["a"].S1 > ranked["b"].S1
    assert ranked["a"].ST > ranked["b"].ST
    # First-order for 'a' dominates total variance
    assert ranked["a"].S1 > 0.8


def test_sobol_indices_rejects_degenerate_inputs():
    with pytest.raises(ValueError):
        sobol_indices(lambda p: 1.0, {}, n_samples=128)
    with pytest.raises(ValueError):
        sobol_indices(lambda p: p["x"], {"x": (0, 1)}, n_samples=2)


def test_sobol_indices_returns_expected_call_count():
    # N*(k+2) model evaluations for k=2, N=32 -> 128 calls
    result = sobol_indices(
        lambda params: params["a"] * params["b"],
        {"a": (0.0, 1.0), "b": (0.0, 1.0)},
        n_samples=32,
        seed=1,
    )
    assert result.total_model_calls == 32 * (2 + 2)
    assert result.n_samples == 32


def test_sobol_indices_bootstrap_attaches_ci():
    result = sobol_indices(
        lambda params: params["a"] + 0.5 * params["b"],
        {"a": (0.0, 1.0), "b": (0.0, 1.0)},
        n_samples=64,
        seed=0,
        bootstrap_resamples=50,
    )
    for idx in result.indices:
        assert isinstance(idx, SobolIndex)
        assert idx.S1_lower is not None
        assert idx.S1_upper is not None
        assert idx.ST_lower is not None
        assert idx.ST_upper is not None
        assert idx.S1_lower <= idx.S1_upper
        assert idx.ST_lower <= idx.ST_upper


def test_sobol_indices_reproducible_with_seed():
    bounds = {"a": (0.0, 1.0), "b": (0.0, 2.0)}
    r1 = sobol_indices(lambda p: p["a"] * p["b"], bounds, n_samples=64, seed=7)
    r2 = sobol_indices(lambda p: p["a"] * p["b"], bounds, n_samples=64, seed=7)
    names_1 = [i.name for i in r1.indices]
    names_2 = [i.name for i in r2.indices]
    assert names_1 == names_2
    for a, b in zip(r1.indices, r2.indices):
        assert a.S1 == pytest.approx(b.S1)
        assert a.ST == pytest.approx(b.ST)


def test_sobol_for_synergy_score_runs_and_returns_bounded_indices():
    # Reduced sample count to keep the smoke test under a second.
    result = sobol_for_synergy_score(
        {"PAH": 3.0, "HCA": 1.5, "benzene": 4.0},
        n_samples=8,
        seed=0,
        tissue="Liver",
    )
    assert isinstance(result, SobolResult)
    assert result.indices, "expected at least one perturbed parameter"
    # Parameter names use "km::<substrate>" / "expr::<enzyme>" prefixes.
    names = {idx.name for idx in result.indices}
    assert any(n.startswith("km::") for n in names)
    assert any(n.startswith("expr::") for n in names)
