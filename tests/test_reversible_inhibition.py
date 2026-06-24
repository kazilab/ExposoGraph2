import json
from math import isfinite

import pytest

from ExposoGraph.interaction_schema import InhibitionMode
from ExposoGraph.reversible_inhibition import (
    EQUATION_ID,
    EQUATION_VERSION,
    compute_reversible_inhibition,
)


def _result(**overrides):
    params = {
        "mode": InhibitionMode.MIXED,
        "substrate_concentration": 20.0,
        "km": 10.0,
        "vmax": 100.0,
        "inhibitor_concentration": 5.0,
        "ki_free_enzyme": 10.0,
        "ki_enzyme_substrate": 20.0,
    }
    params.update(overrides)
    return compute_reversible_inhibition(**params)


def _assert_no_nonfinite(value):
    if isinstance(value, float):
        assert isfinite(value)
    elif isinstance(value, dict):
        for item in value.values():
            _assert_no_nonfinite(item)
    elif isinstance(value, list):
        for item in value:
            _assert_no_nonfinite(item)


def test_zero_inhibitor_returns_unit_modifier():
    result = _result(inhibitor_concentration=0.0)

    assert result.equation_id == EQUATION_ID
    assert result.equation_version == EQUATION_VERSION
    assert result.alpha_free_enzyme == 1.0
    assert result.alpha_enzyme_substrate == 1.0
    assert result.flux_modifier == pytest.approx(1.0)
    assert result.baseline_rate == pytest.approx(result.inhibited_rate)


def test_valid_positive_inputs_produce_bounded_inhibition_modifier():
    result = _result()

    assert 0.0 < result.flux_modifier <= 1.0
    assert result.inhibited_rate <= result.baseline_rate


def test_competitive_inhibition_weakens_at_high_substrate():
    low_substrate = _result(
        mode=InhibitionMode.COMPETITIVE,
        substrate_concentration=1.0,
        ki_free_enzyme=5.0,
        ki_enzyme_substrate=None,
    )
    high_substrate = _result(
        mode=InhibitionMode.COMPETITIVE,
        substrate_concentration=1_000_000.0,
        ki_free_enzyme=5.0,
        ki_enzyme_substrate=None,
    )

    assert low_substrate.flux_modifier < high_substrate.flux_modifier
    assert high_substrate.flux_modifier == pytest.approx(1.0, rel=1e-4)


def test_pure_noncompetitive_inhibition_is_substrate_independent():
    low_substrate = _result(
        mode=InhibitionMode.PURE_NONCOMPETITIVE,
        substrate_concentration=1.0,
        ki_free_enzyme=5.0,
        ki_enzyme_substrate=5.0,
    )
    high_substrate = _result(
        mode=InhibitionMode.PURE_NONCOMPETITIVE,
        substrate_concentration=1_000_000.0,
        ki_free_enzyme=5.0,
        ki_enzyme_substrate=5.0,
    )

    assert low_substrate.flux_modifier == pytest.approx(high_substrate.flux_modifier)
    assert low_substrate.flux_modifier == pytest.approx(1.0 / low_substrate.alpha_free_enzyme)


def test_uncompetitive_inhibition_limits_at_low_and_high_substrate():
    low_substrate = _result(
        mode=InhibitionMode.UNCOMPETITIVE,
        substrate_concentration=1e-9,
        ki_free_enzyme=None,
        ki_enzyme_substrate=5.0,
    )
    high_substrate = _result(
        mode=InhibitionMode.UNCOMPETITIVE,
        substrate_concentration=1_000_000.0,
        ki_free_enzyme=None,
        ki_enzyme_substrate=5.0,
    )

    assert low_substrate.flux_modifier == pytest.approx(1.0)
    assert high_substrate.flux_modifier == pytest.approx(
        1.0 / high_substrate.alpha_enzyme_substrate,
        rel=1e-4,
    )


def test_parent_equation_reduces_to_competitive_when_es_arm_absent():
    result = _result(
        mode=InhibitionMode.COMPETITIVE,
        ki_free_enzyme=5.0,
        ki_enzyme_substrate=None,
    )
    alpha_e = 1.0 + 5.0 / 5.0
    expected = (10.0 + 20.0) / (alpha_e * 10.0 + 20.0)

    assert result.alpha_free_enzyme == pytest.approx(alpha_e)
    assert result.alpha_enzyme_substrate == 1.0
    assert result.flux_modifier == pytest.approx(expected)


def test_mixed_reduces_to_pure_noncompetitive_when_affinities_match():
    mixed = _result(mode=InhibitionMode.MIXED, ki_free_enzyme=5.0, ki_enzyme_substrate=5.0)
    pure = _result(
        mode=InhibitionMode.PURE_NONCOMPETITIVE,
        ki_free_enzyme=5.0,
        ki_enzyme_substrate=5.0,
    )

    assert mixed.alpha_free_enzyme == pytest.approx(mixed.alpha_enzyme_substrate)
    assert mixed.flux_modifier == pytest.approx(pure.flux_modifier)
    assert mixed.apparent_km == pytest.approx(pure.apparent_km)
    assert mixed.apparent_vmax_fraction == pytest.approx(pure.apparent_vmax_fraction)


def test_parent_equation_reduces_to_uncompetitive_when_free_enzyme_arm_absent():
    result = _result(
        mode=InhibitionMode.UNCOMPETITIVE,
        ki_free_enzyme=None,
        ki_enzyme_substrate=5.0,
    )
    alpha_es = 1.0 + 5.0 / 5.0
    expected = (10.0 + 20.0) / (10.0 + alpha_es * 20.0)

    assert result.alpha_free_enzyme == 1.0
    assert result.alpha_enzyme_substrate == pytest.approx(alpha_es)
    assert result.flux_modifier == pytest.approx(expected)


def test_apparent_kinetic_parameters_follow_alpha_definitions():
    result = _result(ki_free_enzyme=4.0, ki_enzyme_substrate=8.0)

    assert result.apparent_vmax_fraction == pytest.approx(1.0 / result.alpha_enzyme_substrate)
    assert result.apparent_km == pytest.approx(
        result.alpha_free_enzyme * 10.0 / result.alpha_enzyme_substrate
    )


def test_unit_scaled_equivalent_concentration_ratios_match():
    micromolar = _result(
        substrate_concentration=5.0,
        km=10.0,
        inhibitor_concentration=2.0,
        ki_free_enzyme=4.0,
        ki_enzyme_substrate=8.0,
    )
    nanomolar_scaled = _result(
        substrate_concentration=5000.0,
        km=10000.0,
        inhibitor_concentration=2000.0,
        ki_free_enzyme=4000.0,
        ki_enzyme_substrate=8000.0,
    )

    assert micromolar.flux_modifier == pytest.approx(nanomolar_scaled.flux_modifier)
    assert micromolar.apparent_vmax_fraction == pytest.approx(
        nanomolar_scaled.apparent_vmax_fraction
    )


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"inhibitor_concentration": -1.0}, "inhibitor_concentration"),
        ({"substrate_concentration": -1.0}, "substrate_concentration"),
        ({"km": 0.0}, "km"),
        ({"km": -1.0}, "km"),
        ({"ki_free_enzyme": 0.0}, "ki_free_enzyme"),
        ({"substrate_concentration": float("nan")}, "finite"),
        ({"mode": InhibitionMode.UNKNOWN}, "unknown inhibition mode"),
        ({"mode": "not-a-mode"}, "unknown inhibition mode"),
        ({"mode": InhibitionMode.MIXED, "ki_enzyme_substrate": None}, "ki_enzyme_substrate"),
        (
            {
                "mode": InhibitionMode.PURE_NONCOMPETITIVE,
                "ki_free_enzyme": 5.0,
                "ki_enzyme_substrate": 6.0,
            },
            "equal Ki",
        ),
    ],
)
def test_invalid_inputs_fail_deterministically(overrides, message):
    with pytest.raises(ValueError, match=message):
        _result(**overrides)


def test_zero_substrate_is_deterministic_and_serializable():
    result = _result(substrate_concentration=0.0)
    payload = result.to_dict()

    assert result.baseline_rate == 0.0
    assert result.inhibited_rate == 0.0
    assert result.flux_modifier == pytest.approx(1.0 / result.alpha_free_enzyme)
    assert payload["warnings"][0]["code"] == "zero_substrate_absolute_rates_zero"
    _assert_no_nonfinite(payload)
    json.dumps(payload)


def test_serialized_result_has_no_nan_or_infinity():
    payload = _result().to_dict()

    _assert_no_nonfinite(payload)
    assert "Infinity" not in json.dumps(payload)
    assert "NaN" not in json.dumps(payload)


def test_reversible_inhibition_kernel_is_not_routed_into_live_engines():
    import ExposoGraph.flux_engine as flux_engine
    import ExposoGraph.interaction_engine as interaction_engine
    import ExposoGraph.unified_api as unified_api

    assert "reversible_inhibition" not in interaction_engine.__dict__
    assert "compute_reversible_inhibition" not in interaction_engine.__dict__
    assert "reversible_inhibition" not in flux_engine.__dict__
    assert "compute_reversible_inhibition" not in flux_engine.__dict__
    assert "reversible_inhibition" not in unified_api.__dict__
    assert "compute_reversible_inhibition" not in unified_api.__dict__
