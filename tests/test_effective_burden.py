import math

import pytest

import ExposoGraph.effective_burden as effective_burden
from ExposoGraph.effective_burden import (
    EffectiveBurdenInput,
    GSHBurdenCouplingInput,
    SusceptibilityModifier,
    compute_effective_carcinogenic_burden,
    couple_gsh_consumption_to_activation_burden,
    effective_burden_from_endpoint_and_gsh,
)
from ExposoGraph.endpoint_toxic_flux import interpret_competitive_endpoint_flux
from ExposoGraph.gsh_redox_capacity import GSHRedoxCapacityInput, compute_gsh_redox_capacity
from ExposoGraph.interaction_schema import ReactionRole, RiskDirectionIfFluxDecreases, SMEReviewStatus
from ExposoGraph.reaction_role_semantics import ReactionRoleAnnotation


def _annotation(role, direction):
    return ReactionRoleAnnotation(
        enzyme="CYPX",
        substrate="substrate",
        reaction_role=role,
        risk_direction_if_flux_decreases=direction,
        review_status=SMEReviewStatus.CURATED,
        record_id="phase8_test_annotation",
    )


def _warning_codes(result):
    return {warning.code for warning in result.warnings}


def test_effective_burden_multiplies_mechanism_factors():
    result = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            activation_burden_ratio=2.0,
            detox_failure_ratio=3.0,
            gsh_relevant=True,
            gsh_detox_penalty_ratio=4.0,
            susceptibility_modifier=0.5,
        )
    )

    assert result.effective_carcinogenic_burden_ratio == pytest.approx(12.0)
    assert result.activation_burden_ratio == pytest.approx(2.0)
    assert result.detox_failure_ratio == pytest.approx(3.0)
    assert result.gsh_detox_penalty_ratio == pytest.approx(4.0)
    assert result.susceptibility_modifier == pytest.approx(0.5)


def test_neutral_defaults_return_effective_burden_one_with_warnings():
    result = compute_effective_carcinogenic_burden()

    assert result.effective_carcinogenic_burden_ratio == pytest.approx(1.0)
    assert result.activation_burden_ratio == pytest.approx(1.0)
    assert result.detox_failure_ratio == pytest.approx(1.0)
    assert result.gsh_detox_penalty_ratio == pytest.approx(1.0)
    assert result.susceptibility_modifier == pytest.approx(1.0)
    assert {
        "activation_burden_ratio_missing_neutral_default",
        "detox_failure_ratio_missing_neutral_default",
        "susceptibility_modifier_missing_neutral_default",
    }.issubset(_warning_codes(result))


def test_phase6_endpoint_toxic_flux_result_feeds_activation_and_detox_ratios():
    endpoint_result = interpret_competitive_endpoint_flux(
        0.5,
        _annotation(ReactionRole.DETOXIFICATION, RiskDirectionIfFluxDecreases.INCREASE),
    )

    burden = effective_burden_from_endpoint_and_gsh(endpoint_result)

    assert burden.activation_burden_ratio == pytest.approx(1.0)
    assert burden.detox_failure_ratio == pytest.approx(2.0)
    assert burden.effective_carcinogenic_burden_ratio == pytest.approx(2.0)
    assert burden.metadata["uses_endpoint_toxic_flux_result"] is True


def test_non_gsh_relevant_substrates_keep_gsh_penalty_neutral():
    result = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            activation_burden_ratio=2.0,
            detox_failure_ratio=2.0,
            gsh_relevant=False,
            gsh_detox_penalty_ratio=10.0,
            susceptibility_modifier=1.0,
        )
    )

    assert result.gsh_detox_penalty_ratio == pytest.approx(1.0)
    assert result.effective_carcinogenic_burden_ratio == pytest.approx(4.0)


def test_gsh_relevant_substrates_scale_consumption_by_upstream_activation():
    result = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(
            gsh_relevant=True,
            base_gsh_consumption_load=1.5,
            upstream_activation_burden_ratio=3.0,
            tissue="default",
        )
    )

    assert result.scaling_source == "explicit_upstream_activation_burden_ratio"
    assert result.base_gsh_consumption_load == pytest.approx(1.5)
    assert result.upstream_activation_burden_ratio == pytest.approx(3.0)
    assert result.gsh_consumption_load_scaled == pytest.approx(4.5)
    assert result.gsh_redox_capacity_result.consumption_load == pytest.approx(4.5)


def test_explicit_upstream_activation_ratio_is_preferred_over_d_times_k():
    result = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(
            gsh_relevant=True,
            base_gsh_consumption_load=1.0,
            upstream_activation_burden_ratio=2.0,
            d_factor=10.0,
            k_factor=10.0,
        )
    )

    assert result.scaling_source == "explicit_upstream_activation_burden_ratio"
    assert result.upstream_activation_burden_ratio == pytest.approx(2.0)
    assert result.gsh_consumption_load_scaled == pytest.approx(2.0)


def test_d_times_k_approximation_is_used_only_when_explicit_ratio_absent():
    result = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(
            gsh_relevant=True,
            base_gsh_consumption_load=1.0,
            d_factor=2.0,
            k_factor=3.0,
        )
    )

    assert result.scaling_source == "d_times_k_approximation"
    assert result.upstream_activation_burden_ratio == pytest.approx(6.0)
    assert result.gsh_consumption_load_scaled == pytest.approx(6.0)
    assert result.metadata["internal_d_or_k_computation"] is False


def test_missing_upstream_activation_for_gsh_relevant_substrate_warns_and_uses_neutral_fallback():
    result = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(gsh_relevant=True, base_gsh_consumption_load=2.0)
    )

    assert result.scaling_source == "neutral_fallback"
    assert result.upstream_activation_burden_ratio == pytest.approx(1.0)
    assert result.gsh_consumption_load_scaled == pytest.approx(2.0)
    assert "gsh_upstream_activation_missing_neutral" in {warning.code for warning in result.warnings}


def test_gsh_coupling_increases_detox_penalty_when_activation_burden_increases():
    low_activation = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(
            gsh_relevant=True,
            base_gsh_consumption_load=1.0,
            upstream_activation_burden_ratio=1.0,
            tissue="default",
        )
    )
    high_activation = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(
            gsh_relevant=True,
            base_gsh_consumption_load=1.0,
            upstream_activation_burden_ratio=5.0,
            tissue="default",
        )
    )

    assert high_activation.gsh_fraction < low_activation.gsh_fraction
    assert high_activation.detox_penalty_multiplier > low_activation.detox_penalty_multiplier


def test_effective_burden_can_use_phase8_gsh_coupling_result():
    result = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            activation_burden_ratio=2.0,
            detox_failure_ratio=1.0,
            susceptibility_modifier=1.0,
            gsh_relevant=True,
            gsh_coupling=GSHBurdenCouplingInput(
                gsh_relevant=True,
                base_gsh_consumption_load=1.0,
                upstream_activation_burden_ratio=3.0,
            ),
        )
    )

    assert result.gsh_consumption_load == pytest.approx(1.0)
    assert result.gsh_consumption_load_scaled == pytest.approx(3.0)
    assert result.gsh_detox_penalty_ratio > 1.0
    assert result.effective_carcinogenic_burden_ratio > 2.0


def test_phase7_gsh_result_can_feed_gsh_detox_penalty_without_rerunning_gsh():
    gsh_result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(tissue="default", consumption_load=3.0)
    )

    result = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            activation_burden_ratio=1.0,
            detox_failure_ratio=1.0,
            susceptibility_modifier=1.0,
            gsh_relevant=True,
            gsh_redox_capacity_result=gsh_result,
        )
    )

    assert result.gsh_detox_penalty_ratio == pytest.approx(gsh_result.detox_penalty_multiplier)
    assert result.gsh_fraction == pytest.approx(gsh_result.gsh_fraction)
    assert result.metadata["uses_gsh_redox_capacity_result"] is True


def test_negative_and_invalid_ratios_warn_and_return_bounded_output():
    result = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            activation_burden_ratio=-2.0,
            detox_failure_ratio=math.nan,
            gsh_relevant=True,
            gsh_detox_penalty_ratio=-3.0,
            susceptibility_modifier=-1.0,
        )
    )

    assert result.effective_carcinogenic_burden_ratio == pytest.approx(0.0)
    assert result.activation_burden_ratio == pytest.approx(0.0)
    assert result.detox_failure_ratio == pytest.approx(1.0)
    assert result.gsh_detox_penalty_ratio == pytest.approx(0.0)
    assert result.susceptibility_modifier == pytest.approx(0.0)
    assert {
        "activation_burden_ratio_negative_clamped",
        "detox_failure_ratio_invalid_neutral_default",
        "gsh_detox_penalty_ratio_negative_clamped",
        "susceptibility_modifier_negative_clamped",
    }.issubset(_warning_codes(result))


def test_susceptibility_modifier_is_applied_and_default_is_neutral():
    with_modifier = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            activation_burden_ratio=2.0,
            detox_failure_ratio=2.0,
            susceptibility_modifier=SusceptibilityModifier(modifier_ratio=1.5, label="test_modifier"),
        )
    )
    default_modifier = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(activation_burden_ratio=2.0, detox_failure_ratio=2.0)
    )

    assert with_modifier.effective_carcinogenic_burden_ratio == pytest.approx(6.0)
    assert default_modifier.effective_carcinogenic_burden_ratio == pytest.approx(4.0)
    assert with_modifier.susceptibility_modifier == pytest.approx(1.5)


def test_no_adjusted_risk_public_output_is_produced():
    result_dict = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(activation_burden_ratio=2.0, detox_failure_ratio=1.0)
    ).to_dict()

    assert "adjusted_risk" not in result_dict
    assert "adjusted_risk" not in result_dict["metadata"]
    assert result_dict["metadata"]["public_risk_output"] == "not_produced_or_modified"


def test_no_phase9_shapley_state_or_output_is_produced():
    result_dict = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(activation_burden_ratio=1.0, detox_failure_ratio=1.0)
    ).to_dict()

    assert result_dict["metadata"]["shapley_attribution"] is False
    assert result_dict["metadata"]["phase9_behavior"] is False
    assert "shapley_values" not in result_dict
    assert "attribution_state" not in result_dict


def test_no_interaction_engine_integration_occurs():
    assert "interaction_engine" not in effective_burden.__dict__
    assert "compute_interaction_matrix" not in effective_burden.__dict__
    assert "competitive_inhibition_flux" not in effective_burden.__dict__
    assert "gsh_depletion_model" not in effective_burden.__dict__
    assert "get_ki" not in effective_burden.__dict__


def test_no_pbpk_ode_overclaim_appears_in_metadata():
    result = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(activation_burden_ratio=1.0, detox_failure_ratio=1.0)
    )
    metadata = result.metadata
    values = " ".join(str(value).lower() for value in metadata.values())

    assert metadata["model_family"] == "semi_mechanistic_relative_burden"
    assert metadata["validated_pbpk_ode_model"] is False
    assert metadata["clinical_risk_model"] is False
    assert "fully_validated" not in values
    assert "clinical prediction" not in values


def test_result_serialization_is_json_friendly():
    result_dict = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            activation_burden_ratio=2.0,
            detox_failure_ratio=1.5,
            susceptibility_modifier=1.25,
            gsh_relevant=True,
            gsh_coupling=GSHBurdenCouplingInput(
                gsh_relevant=True,
                base_gsh_consumption_load=1.0,
                d_factor=2.0,
                k_factor=2.0,
            ),
        )
    ).to_dict()

    assert result_dict["effective_carcinogenic_burden_ratio"] > 0.0
    assert result_dict["gsh_coupling_result"]["scaling_source"] == "d_times_k_approximation"
    assert isinstance(result_dict["warnings"], list)
    assert result_dict["metadata"]["relative_burden_formula"] == (
        "ActivationBurdenRatio * DetoxFailureRatio * GSHDetoxPenaltyRatio * SusceptibilityModifier"
    )
