import pytest

import ExposoGraph.mechanism_attribution as mechanism_attribution
from ExposoGraph.mechanism_attribution import (
    MechanismName,
    MechanismState,
    MechanismStateValue,
    compute_explicit_mechanism_interactions,
    compute_mechanism_attribution,
    compute_shapley_main_effects,
    evaluate_mechanism_states,
    generate_mechanism_states,
    validate_eight_state_values,
)


def _state_values(model):
    return {state: model(state) for state in generate_mechanism_states()}


def _effect_by_mechanism(effects):
    return {effect.mechanism: effect.effect for effect in effects}


def _interaction_by_key(terms):
    return {tuple(mechanism.value for mechanism in term.mechanisms): term.effect for term in terms}


def test_all_eight_mechanism_states_are_generated_deterministically():
    states = generate_mechanism_states()

    assert [state.key for state in states] == [
        "none",
        "induction",
        "competition",
        "gsh",
        "induction+competition",
        "induction+gsh",
        "competition+gsh",
        "induction+competition+gsh",
    ]
    assert len(set(states)) == 8


def test_missing_one_state_raises_clear_value_error():
    values = _state_values(lambda state: 1.0)
    values.pop(MechanismState.from_active((MechanismName.GSH,)))

    with pytest.raises(ValueError, match="Missing Phase 9 mechanism state values: gsh"):
        validate_eight_state_values(values)


def test_non_finite_state_value_raises_clear_value_error():
    values = _state_values(lambda state: 1.0)
    values[MechanismState.from_active((MechanismName.INDUCTION,))] = float("nan")

    with pytest.raises(ValueError, match="State value for induction must be finite"):
        compute_mechanism_attribution(values)


def test_additive_state_values_give_shapley_equal_to_singletons_and_zero_interactions():
    values = _state_values(
        lambda state: 10.0
        + (2.0 if state.induction else 0.0)
        + (3.0 if state.competition else 0.0)
        + (4.0 if state.gsh else 0.0)
    )

    result = compute_mechanism_attribution(values)
    shapley = _effect_by_mechanism(result.shapley_main_effects)
    singletons = _interaction_by_key(result.singleton_effects)
    pairs = _interaction_by_key(result.pairwise_interactions)

    assert shapley[MechanismName.INDUCTION] == pytest.approx(2.0)
    assert shapley[MechanismName.COMPETITION] == pytest.approx(3.0)
    assert shapley[MechanismName.GSH] == pytest.approx(4.0)
    assert singletons[("induction",)] == pytest.approx(2.0)
    assert singletons[("competition",)] == pytest.approx(3.0)
    assert singletons[("gsh",)] == pytest.approx(4.0)
    assert all(value == pytest.approx(0.0) for value in pairs.values())
    assert result.three_way_interaction.effect == pytest.approx(0.0)


def test_pairwise_interaction_case_reports_correct_pairwise_term():
    values = _state_values(
        lambda state: 1.0
        + (2.0 if state.induction else 0.0)
        + (3.0 if state.competition else 0.0)
        + (4.0 if state.gsh else 0.0)
        + (5.0 if state.induction and state.competition else 0.0)
    )

    result = compute_mechanism_attribution(values)
    pairs = _interaction_by_key(result.pairwise_interactions)

    assert pairs[("induction", "competition")] == pytest.approx(5.0)
    assert pairs[("induction", "gsh")] == pytest.approx(0.0)
    assert pairs[("competition", "gsh")] == pytest.approx(0.0)
    assert result.three_way_interaction.effect == pytest.approx(0.0)


def test_three_way_interaction_case_reports_correct_three_way_term():
    values = _state_values(
        lambda state: 1.0 + (11.0 if state.induction and state.competition and state.gsh else 0.0)
    )

    result = compute_mechanism_attribution(values)

    assert all(term.effect == pytest.approx(0.0) for term in result.singleton_effects)
    assert all(term.effect == pytest.approx(0.0) for term in result.pairwise_interactions)
    assert result.three_way_interaction.effect == pytest.approx(11.0)


def test_shapley_main_effects_sum_to_total_effect_within_tolerance():
    values = _state_values(
        lambda state: 2.0
        + (1.0 if state.induction else 0.0)
        + (2.0 if state.competition else 0.0)
        + (3.0 if state.gsh else 0.0)
        + (4.0 if state.induction and state.gsh else 0.0)
        + (5.0 if state.induction and state.competition and state.gsh else 0.0)
    )

    result = compute_mechanism_attribution(values)

    assert sum(effect.effect for effect in result.shapley_main_effects) == pytest.approx(result.total_effect)
    assert result.shapley_residual == pytest.approx(0.0)


def test_explicit_interaction_reconstruction_has_zero_residual_within_tolerance():
    values = _state_values(
        lambda state: 3.0
        + (1.0 if state.induction else 0.0)
        + (2.0 if state.competition else 0.0)
        + (3.0 if state.gsh else 0.0)
        + (4.0 if state.induction and state.competition else 0.0)
        + (5.0 if state.induction and state.gsh else 0.0)
        + (6.0 if state.competition and state.gsh else 0.0)
        + (7.0 if state.induction and state.competition and state.gsh else 0.0)
    )

    result = compute_mechanism_attribution(values)
    reconstructed = (
        sum(term.effect for term in result.singleton_effects)
        + sum(term.effect for term in result.pairwise_interactions)
        + result.three_way_interaction.effect
    )

    assert reconstructed == pytest.approx(result.total_effect)
    assert result.interaction_reconstruction_residual == pytest.approx(0.0)


def test_pairwise_interactions_are_reported_for_all_three_pairs():
    result = compute_mechanism_attribution(_state_values(lambda state: float(len(state.active_mechanisms))))

    assert [tuple(mechanism.value for mechanism in term.mechanisms) for term in result.pairwise_interactions] == [
        ("induction", "competition"),
        ("induction", "gsh"),
        ("competition", "gsh"),
    ]


def test_three_way_interaction_is_reported():
    result = compute_mechanism_attribution(_state_values(lambda state: float(len(state.active_mechanisms))))

    assert tuple(mechanism.value for mechanism in result.three_way_interaction.mechanisms) == (
        "induction",
        "competition",
        "gsh",
    )
    assert result.three_way_interaction.term_type == "three_way_interaction"


def test_residuals_are_zero_within_tolerance_for_valid_complete_input():
    result = compute_mechanism_attribution(_state_values(lambda state: 1.0 + len(state.active_mechanisms)))

    assert result.residuals_are_zero_within_tolerance is True
    assert result.shapley_residual == pytest.approx(0.0)
    assert result.interaction_reconstruction_residual == pytest.approx(0.0)
    assert result.warnings == []


def test_evaluator_helper_calls_caller_provided_evaluator_for_exactly_eight_states():
    calls = []

    def evaluator(state):
        calls.append(state.key)
        return float(len(calls))

    evaluated = evaluate_mechanism_states(evaluator)

    assert calls == [state.key for state in generate_mechanism_states()]
    assert len(evaluated) == 8
    assert all(isinstance(item, MechanismStateValue) for item in evaluated)


def test_shapley_and_explicit_helpers_accept_serializable_state_value_records():
    records = [MechanismStateValue(state=state, value=float(index), state_key=state.key) for index, state in enumerate(generate_mechanism_states())]

    shapley = compute_shapley_main_effects(records)
    singleton_terms, pairwise_terms, three_way = compute_explicit_mechanism_interactions(records)

    assert len(shapley) == 3
    assert len(singleton_terms) == 3
    assert len(pairwise_terms) == 3
    assert three_way.order == 3


def test_no_adjusted_risk_public_output_is_produced():
    result_dict = compute_mechanism_attribution(_state_values(lambda state: 1.0)).to_dict()

    assert "adjusted_risk" not in result_dict
    assert "adjusted_risk" not in result_dict["metadata"]


def test_no_interaction_engine_integration_occurs():
    assert "interaction_engine" not in mechanism_attribution.__dict__
    assert "compute_interaction_matrix" not in mechanism_attribution.__dict__
    assert "competitive_inhibition_flux" not in mechanism_attribution.__dict__
    assert "compute_effective_carcinogenic_burden" not in mechanism_attribution.__dict__
    assert "compute_gsh_redox_capacity" not in mechanism_attribution.__dict__
    assert "interpret_competitive_endpoint_flux" not in mechanism_attribution.__dict__
    assert "get_ki" not in mechanism_attribution.__dict__


def test_no_phase10_model_card_or_uncertainty_output_behavior_is_implemented():
    result_dict = compute_mechanism_attribution(_state_values(lambda state: 1.0)).to_dict()

    assert "model_card" not in result_dict
    assert "uncertainty_summary" not in result_dict
    assert result_dict["metadata"]["phase10_behavior"] is False
    assert result_dict["metadata"]["uncertainty_propagation"] is False


def test_output_serialization_is_json_friendly_if_to_dict_is_available():
    result_dict = compute_mechanism_attribution(_state_values(lambda state: 1.0 + len(state.active_mechanisms))).to_dict()

    assert result_dict["mechanisms"] == ["induction", "competition", "gsh"]
    assert result_dict["state_values"][0]["state_key"] == "none"
    assert result_dict["shapley_main_effects"][0]["mechanism"] == "induction"
    assert result_dict["pairwise_interactions"][0]["mechanisms"] == ["induction", "competition"]


def test_metadata_states_model_output_attribution_not_biological_causality():
    result = compute_mechanism_attribution(_state_values(lambda state: 1.0 + len(state.active_mechanisms)))

    assert result.metadata["attribution_type"] == "model_output_attribution"
    assert result.metadata["not_biological_causality"] is True
    assert result.metadata["residual_policy"] == "no_unexplained_residual_term"
