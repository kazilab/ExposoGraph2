import json
import math

import pytest

from ExposoGraph.endpoint_toxic_flux import interpret_competitive_endpoint_flux
from ExposoGraph.interaction_engine import _interaction_matrix_to_compat_dict, cli_main, competitive_inhibition_flux, compute_interaction_matrix
from ExposoGraph.interaction_schema import ConcentrationBasis, InhibitionMode, ReactionRole, RiskDirectionIfFluxDecreases
from ExposoGraph.reaction_role_semantics import ReactionRoleAnnotation
from ExposoGraph.unified_api import patient_risk_query


def _finite_walk(value):
    if isinstance(value, float):
        assert math.isfinite(value)
    elif isinstance(value, dict):
        for item in value.values():
            _finite_walk(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _finite_walk(item)


def _bio(enzyme, substrate, substrates):
    return competitive_inhibition_flux(enzyme, substrates).substrates[substrate].biological_output


def test_live_modified_flux_reaches_reaction_role_and_endpoint_outputs():
    ndma = _bio("CYP2E1", "NDMA", {"NDMA": 1.0, "ethanol": 2000.0})
    benzene = _bio("CYP2E1", "benzene", {"benzene": 10.0, "ethanol": 2000.0})

    assert ndma["reaction_role_interpretation"]["role"] == "bioactivation"
    assert ndma["reaction_role_interpretation"]["directional_interpretation"] == "toxic_product_formation_may_decrease"
    assert ndma["endpoint_toxic_flux"]["activation_burden_ratio"] < 1.0

    assert benzene["reaction_role_interpretation"]["role"] == "detoxification"
    assert benzene["reaction_role_interpretation"]["directional_interpretation"] == "burden_may_increase"
    assert benzene["endpoint_toxic_flux"]["detox_failure_ratio"] > 1.0
    assert benzene["effective_burden"]["effective_carcinogenic_burden_ratio"] > 1.0


def test_elimination_and_protective_repair_roles_increase_burden_when_inhibited():
    clearance = ReactionRoleAnnotation(
        enzyme="CLEAR",
        substrate="parent",
        reaction_role=ReactionRole.CLEARANCE,
        risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.INCREASE,
    )
    repair = ReactionRoleAnnotation(
        enzyme="REPAIR",
        substrate="damage",
        reaction_role=ReactionRole.PROTECTIVE_REPAIR,
        risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.INCREASE,
    )

    assert interpret_competitive_endpoint_flux(0.5, clearance).detox_failure_ratio == pytest.approx(2.0)
    assert interpret_competitive_endpoint_flux(0.5, repair).detox_failure_ratio == pytest.approx(2.0)


def test_unknown_role_and_unresolved_inhibition_withhold_directional_interpretation():
    unknown = _bio("CYP2E1", "acetaminophen", {"acetaminophen": 10.0, "ethanol": 2000.0})
    unresolved = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={
            "benzene": {
                "mode": InhibitionMode.UNKNOWN,
                "inhibitor_concentration_uM": 2.0,
                "concentration_basis": ConcentrationBasis.UNBOUND,
                "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
            }
        },
    ).substrates["benzene"].biological_output

    assert unknown["reaction_role_interpretation"]["role"] == "unknown"
    assert unknown["reaction_role_interpretation"]["directional_interpretation"] == "withheld_unknown_role"
    assert unknown["reaction_role_interpretation"]["review_required"] is True

    assert unresolved["kinetic_effect"]["mechanism_state"] == "mechanism_unresolved"
    assert unresolved["endpoint_toxic_flux"]["endpoint_toxic_flux_ratio"] is None
    assert unresolved["effective_burden"]["effective_carcinogenic_burden_ratio"] is None


def test_gsh_capacity_changes_only_for_biologically_annotated_pathways():
    acetaminophen = _bio("CYP2E1", "acetaminophen", {"acetaminophen": 10.0, "ethanol": 2000.0})
    ndma = _bio("CYP2E1", "NDMA", {"NDMA": 1.0, "ethanol": 2000.0})

    assert acetaminophen["gsh_capacity_effect"]["gsh_relevant"] is True
    assert acetaminophen["gsh_capacity_effect"]["gsh_fraction"] is not None
    assert acetaminophen["gsh_capacity_effect"]["detox_penalty_multiplier"] >= 1.0

    assert ndma["gsh_capacity_effect"]["gsh_relevant"] is False
    assert ndma["gsh_capacity_effect"]["gsh_fraction"] is None
    assert ndma["gsh_capacity_effect"]["detox_penalty_multiplier"] == pytest.approx(1.0)


def test_mechanism_attribution_uses_live_engine_generated_eight_states():
    result = compute_interaction_matrix(
        {"benzene": 1.0, "ethanol": 1.0},
        enable_induction=False,
        enable_gsh_depletion=False,
    )

    attribution = result.mechanism_attribution
    assert attribution["state_calculation_source"] == "interaction_engine.compute_interaction_matrix"
    assert len(attribution["state_values"]) == 8
    assert attribution["metadata"]["caller_metadata"]["engine_generated_states"] is True
    assert attribution["metadata"]["caller_metadata"]["disabled_inhibition_modifier"] == pytest.approx(1.0)
    assert attribution["mechanism_state_distinctions"]["mechanism_unresolved"]


def test_multiple_explicit_inhibition_contexts_remain_review_required_in_biological_output():
    result = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={
            "benzene": [
                {
                    "mode": InhibitionMode.COMPETITIVE,
                    "ki_free_enzyme_uM": 5.0,
                    "inhibitor_concentration_uM": 2.0,
                    "concentration_basis": ConcentrationBasis.UNBOUND,
                    "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
                },
                {
                    "mode": InhibitionMode.UNCOMPETITIVE,
                    "ki_enzyme_substrate_uM": 6.0,
                    "inhibitor_concentration_uM": 1.5,
                    "concentration_basis": ConcentrationBasis.UNBOUND,
                    "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
                },
            ]
        },
    )
    benzene = result.substrates["benzene"]

    assert benzene.kinetic_modifier is None
    assert benzene.biological_output["kinetic_effect"]["review_required"] is True
    assert "MULTIPLE_INHIBITORS_NOT_JOINTLY_RESOLVED" in benzene.biological_output["kinetic_effect"]["warnings"]


def test_model_transparency_and_serialization_expose_biological_output_assumptions_without_nan():
    output = _bio("CYP2E1", "benzene", {"benzene": 10.0, "ethanol": 2000.0})

    assert output["kinetic_effect"]["equation_id"] == "reversible_inhibition.mixed.v1"
    assert output["model_transparency"]["live_engine_integration"] is True
    assert output["model_transparency"]["model_card_summary"]["validation_summary"]["live_biological_output_integration"] is True
    json.dumps(output, allow_nan=False)
    _finite_walk(output)


def test_unified_api_and_cli_outputs_are_additive_and_backward_compatible(tmp_path):
    profile = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Liver",
        include_tissue_report=False,
    )
    payload = profile.biological_output_integration
    assert "substrate_outputs" in payload
    assert "mechanism_attribution" in payload
    json.dumps(payload, allow_nan=False)

    result = compute_interaction_matrix({"benzene": 1.0, "ethanol": 1.0})
    compat = _interaction_matrix_to_compat_dict(result)
    benzene = compat["competitive_effects"]["CYP2E1"]["benzene"]
    assert {
        "single_flux",
        "competitive_flux",
        "flux_change_fraction",
        "inhibition_term",
        "activated_product_flux",
        "Km_uM",
        "concentration_uM",
        "product",
        "product_carcinogenic",
    }.issubset(benzene)
    assert "biological_output" in benzene

    output_path = tmp_path / "interaction.json"
    assert cli_main(["--profile", "smoker", "--output-json", str(output_path)]) == 0
    cli_payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert "individual_risks" in cli_payload
    assert "mechanism_attribution" in cli_payload
    json.dumps(cli_payload, allow_nan=False)
