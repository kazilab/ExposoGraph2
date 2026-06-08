import pytest

import ExposoGraph.gsh_redox_capacity as gsh_redox_capacity
from ExposoGraph.gsh_redox_capacity import (
    GSHModelVersion,
    GSHRedoxCapacityInput,
    compute_gsh_redox_capacity,
    compute_quasi_steady_gsh_fraction,
    get_default_gsh_tissue_presets,
    legacy_gsh_depletion_fraction,
)


def _warning_codes(result):
    return {warning.code for warning in result.warnings}


def test_zero_consumption_gives_full_gsh_fraction():
    result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(tissue="liver", consumption_load=0.0)
    )

    assert result.gsh_fraction == pytest.approx(1.0)
    assert result.redox_capacity_ratio == pytest.approx(1.0)
    assert result.detox_penalty_multiplier == pytest.approx(1.0)


def test_high_consumption_reduces_gsh_fraction():
    low_load = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(tissue="liver", consumption_load=0.25)
    )
    high_load = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(tissue="liver", consumption_load=10.0)
    )

    assert high_load.gsh_fraction < low_load.gsh_fraction
    assert high_load.redox_capacity_ratio < low_load.redox_capacity_ratio


def test_zero_synthesis_with_positive_consumption_gives_zero_fraction():
    result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(
            tissue="default",
            synthesis_capacity=0.0,
            consumption_load=1.0,
        )
    )

    assert result.gsh_fraction == pytest.approx(0.0)
    assert result.redox_capacity_ratio == pytest.approx(0.0)
    assert result.detox_penalty_multiplier > 1.0


def test_quasi_steady_helper_bounds_expected_edges():
    assert compute_quasi_steady_gsh_fraction(1.0, 0.0) == pytest.approx(1.0)
    assert compute_quasi_steady_gsh_fraction(0.0, 2.0) == pytest.approx(0.0)
    assert 0.0 <= compute_quasi_steady_gsh_fraction(1.0, 3.0) <= 1.0


def test_final_fraction_and_capacity_are_clamped_to_unit_interval():
    result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(
            tissue="default",
            synthesis_capacity=1.0,
            consumption_load=0.0,
            baseline_capacity=2.0,
        )
    )

    assert 0.0 <= result.gsh_fraction <= 1.0
    assert 0.0 <= result.redox_capacity_ratio <= 1.0
    assert result.clamped is True


def test_negative_and_invalid_inputs_warn_and_return_bounded_output():
    negative = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(
            tissue="default",
            synthesis_capacity=-1.0,
            consumption_load=-2.0,
            turnover_capacity=-3.0,
        )
    )
    invalid = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(tissue="default", consumption_load="not-a-number")
    )

    assert 0.0 <= negative.gsh_fraction <= 1.0
    assert 0.0 <= invalid.gsh_fraction <= 1.0
    assert {
        "synthesis_capacity_negative_clamped",
        "consumption_load_negative_clamped",
        "turnover_capacity_negative_clamped",
    }.issubset(_warning_codes(negative))
    assert "consumption_load_invalid" in _warning_codes(invalid)


def test_tissue_presets_are_deterministic_and_include_required_contexts():
    first = get_default_gsh_tissue_presets()
    second = get_default_gsh_tissue_presets()

    assert {"default", "liver", "lung", "bone_marrow", "hematopoietic"}.issubset(first)
    assert first["liver"].to_dict() == second["liver"].to_dict()
    assert first["bone_marrow"].relative_synthesis_capacity == pytest.approx(
        first["hematopoietic"].relative_synthesis_capacity
    )


def test_unknown_tissue_falls_back_to_default_with_warning():
    result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(tissue="unmapped tissue", consumption_load=0.5)
    )

    assert result.tissue == "default"
    assert "unknown_tissue_default_preset" in _warning_codes(result)


def test_detox_penalty_increases_as_gsh_fraction_decreases():
    moderate = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(synthesis_capacity=1.0, consumption_load=1.0)
    )
    depleted = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(synthesis_capacity=1.0, consumption_load=9.0)
    )

    assert depleted.gsh_fraction < moderate.gsh_fraction
    assert depleted.detox_penalty_multiplier > moderate.detox_penalty_multiplier


def test_legacy_behavior_is_explicitly_named_legacy():
    assert "legacy" in legacy_gsh_depletion_fraction.__name__
    assert "legacy" in GSHModelVersion.LEGACY_DETACHED_GSH_PENALTY.value
    assert legacy_gsh_depletion_fraction(1.0, 0.25) == pytest.approx(0.75)


def test_result_serialization_is_json_friendly():
    result_dict = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(tissue="lung", consumption_load=0.5)
    ).to_dict()

    assert result_dict["model_version"] == "phase7_quasi_steady_relative_capacity"
    assert result_dict["tissue"] == "lung"
    assert isinstance(result_dict["warnings"], list)
    assert result_dict["metadata"]["quasi_steady_expression"] == (
        "synthesis_capacity / (synthesis_capacity + consumption_load)"
    )


def test_no_adjusted_risk_field_is_produced():
    result_dict = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(consumption_load=0.5)
    ).to_dict()

    assert "adjusted_risk" not in result_dict
    assert "adjusted_risk" not in result_dict["metadata"]


def test_no_phase8_coupling_to_activation_burden_is_implemented():
    result_dict = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(consumption_load=0.5)
    ).to_dict()

    assert result_dict["metadata"]["phase8_coupling"] is False
    assert "activation_burden_ratio" not in result_dict
    assert "endpoint_toxic_flux_ratio" not in result_dict
    assert "susceptibility_modifier" not in result_dict
    assert "endpoint_toxic_flux" not in gsh_redox_capacity.__dict__


def test_no_interaction_engine_integration_occurs():
    assert "interaction_engine" not in gsh_redox_capacity.__dict__
    assert "compute_interaction_matrix" not in gsh_redox_capacity.__dict__
    assert "gsh_depletion_model" not in gsh_redox_capacity.__dict__


def test_no_pbpk_ode_overclaim_in_model_metadata():
    result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(consumption_load=0.5)
    )
    metadata = result.metadata

    assert metadata["model_family"] == "semi_mechanistic_relative_capacity"
    assert metadata["validated_pbpk_ode_gsh_gssg_nrf2"] is False
    assert "fully_validated" not in " ".join(str(value).lower() for value in metadata.values())
