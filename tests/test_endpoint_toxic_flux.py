import pytest

import ExposoGraph.endpoint_toxic_flux as endpoint_toxic_flux
from ExposoGraph.endpoint_toxic_flux import (
    EndpointFluxInput,
    endpoint_toxic_flux_from_registry_lookup,
    interpret_competitive_endpoint_flux,
    interpret_endpoint_toxic_flux,
)
from ExposoGraph.interaction_schema import (
    ReactionRole,
    RiskDirectionIfFluxDecreases,
    RiskEndpoint,
    SMEReviewStatus,
)
from ExposoGraph.reaction_role_semantics import ReactionRoleAnnotation


def _annotation(role, direction, **kwargs):
    return ReactionRoleAnnotation(
        enzyme=kwargs.pop("enzyme", "CYPX"),
        substrate=kwargs.pop("substrate", "substrate"),
        reaction_role=role,
        risk_direction_if_flux_decreases=direction,
        **kwargs,
    )


def _warning_codes(result):
    return {warning.code for warning in result.warnings}


def test_bioactivation_decrease_direction_lowers_activation_burden():
    annotation = _annotation(
        ReactionRole.BIOACTIVATION,
        RiskDirectionIfFluxDecreases.DECREASE,
        record_id="bioactivation_case",
    )

    result = interpret_competitive_endpoint_flux(
        0.25,
        annotation,
        enzyme="CYP1A1",
        substrate="benzene",
        endpoint=RiskEndpoint.DNA_ADDUCT,
    )

    assert result.flux_ratio == pytest.approx(0.25)
    assert result.endpoint_toxic_flux_ratio == pytest.approx(0.25)
    assert result.activation_burden_ratio == pytest.approx(0.25)
    assert result.detox_failure_ratio == pytest.approx(1.0)
    assert result.reaction_role is ReactionRole.BIOACTIVATION
    assert "adjusted_risk" not in result.to_dict()


def test_detoxification_increase_direction_raises_detox_failure_burden():
    annotation = _annotation(
        ReactionRole.DETOXIFICATION,
        RiskDirectionIfFluxDecreases.INCREASE,
        record_id="detox_case",
    )

    result = interpret_competitive_endpoint_flux(0.25, annotation)

    assert result.endpoint_toxic_flux_ratio == pytest.approx(4.0)
    assert result.activation_burden_ratio == pytest.approx(1.0)
    assert result.detox_failure_ratio == pytest.approx(4.0)
    assert result.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.INCREASE


def test_clearance_role_raises_detox_failure_burden_when_flux_decreases():
    annotation = _annotation(
        ReactionRole.CLEARANCE,
        RiskDirectionIfFluxDecreases.INCREASE,
        record_id="clearance_case",
    )

    result = interpret_competitive_endpoint_flux(0.5, annotation)

    assert result.endpoint_toxic_flux_ratio == pytest.approx(2.0)
    assert result.detox_failure_ratio == pytest.approx(2.0)
    assert result.activation_burden_ratio == pytest.approx(1.0)


def test_unknown_role_returns_neutral_with_warning():
    annotation = _annotation(
        ReactionRole.UNKNOWN,
        RiskDirectionIfFluxDecreases.UNKNOWN,
        record_id="unknown_case",
    )

    result = interpret_competitive_endpoint_flux(0.25, annotation)

    assert result.endpoint_toxic_flux_ratio == pytest.approx(1.0)
    assert result.activation_burden_ratio == pytest.approx(1.0)
    assert result.detox_failure_ratio == pytest.approx(1.0)
    assert _warning_codes(result) >= {"reaction_role_unknown", "endpoint_role_unknown"}


def test_dual_role_or_mixed_direction_is_neutral_without_context_resolution():
    annotation = _annotation(
        ReactionRole.DUAL_ROLE,
        RiskDirectionIfFluxDecreases.MIXED,
        record_id="dual_mixed_case",
    )

    result = interpret_competitive_endpoint_flux(0.25, annotation)

    assert result.endpoint_toxic_flux_ratio == pytest.approx(1.0)
    assert result.activation_burden_ratio == pytest.approx(1.0)
    assert result.detox_failure_ratio == pytest.approx(1.0)
    assert {"endpoint_dual_role_neutral", "endpoint_direction_mixed"}.issubset(_warning_codes(result))


def test_probe_only_is_neutral_with_warning():
    annotation = _annotation(
        ReactionRole.PROBE_ONLY,
        RiskDirectionIfFluxDecreases.NEUTRAL,
        record_id="probe_case",
    )

    result = interpret_competitive_endpoint_flux(0.2, annotation)

    assert result.endpoint_toxic_flux_ratio == pytest.approx(1.0)
    assert "endpoint_probe_only_neutral" in _warning_codes(result)


def test_inactive_pending_tce_annotation_is_neutral_and_warned():
    result = endpoint_toxic_flux_from_registry_lookup("CYP2E1", "TCE", 0.5)

    assert result.annotation_record_id == "reaction_role_tce_cyp2e1_candidate_pending"
    assert result.endpoint_toxic_flux_ratio == pytest.approx(1.0)
    assert result.activation_burden_ratio == pytest.approx(1.0)
    assert result.detox_failure_ratio == pytest.approx(1.0)
    assert {"reaction_role_inactive", "endpoint_role_inactive", "endpoint_role_pending"}.issubset(
        _warning_codes(result)
    )


def test_benzene_cyp2e1_record_maps_flux_decrease_to_higher_burden():
    result = endpoint_toxic_flux_from_registry_lookup(
        "CYP2E1",
        "benzene",
        0.5,
        endpoint=RiskEndpoint.DNA_ADDUCT,
    )

    assert result.annotation_record_id == "reaction_role_benzene_cyp2e1_v2_direction_adjustment"
    assert result.reaction_role is ReactionRole.DETOXIFICATION
    assert result.endpoint_toxic_flux_ratio == pytest.approx(2.0)
    assert result.detox_failure_ratio == pytest.approx(2.0)
    assert result.activation_burden_ratio == pytest.approx(1.0)


def test_benzene_cyp2f1_and_cyp2a13_context_specific_behavior_is_preserved():
    for enzyme in ("CYP2F1", "CYP2A13"):
        result = endpoint_toxic_flux_from_registry_lookup(
            enzyme,
            "benzene",
            0.5,
            tissue="bone marrow",
            endpoint=RiskEndpoint.DNA_ADDUCT,
        )

        assert result.annotation_record_id == f"reaction_role_benzene_{enzyme.lower()}_bone_marrow_direction_adjustment"
        assert result.endpoint_toxic_flux_ratio == pytest.approx(2.0)
        assert result.detox_failure_ratio == pytest.approx(2.0)


def test_benzene_cyp2f1_and_cyp2a13_outside_context_do_not_globally_flip():
    for enzyme in ("CYP2F1", "CYP2A13"):
        result = endpoint_toxic_flux_from_registry_lookup(enzyme, "benzene", 0.5)

        assert result.annotation_record_id == f"reaction_role_benzene_{enzyme.lower()}_outside_context_unknown"
        assert result.endpoint_toxic_flux_ratio == pytest.approx(1.0)
        assert result.detox_failure_ratio == pytest.approx(1.0)
        assert "endpoint_role_unknown" in _warning_codes(result)


def test_product_carcinogenic_metadata_alone_does_not_drive_interpretation():
    hazardous = _annotation(
        ReactionRole.UNKNOWN,
        RiskDirectionIfFluxDecreases.UNKNOWN,
        metadata={"product_carcinogenic": True},
    )
    nonhazardous = _annotation(
        ReactionRole.UNKNOWN,
        RiskDirectionIfFluxDecreases.UNKNOWN,
        metadata={"product_carcinogenic": False},
    )

    hazardous_result = interpret_competitive_endpoint_flux(0.2, hazardous)
    nonhazardous_result = interpret_competitive_endpoint_flux(0.2, nonhazardous)

    assert hazardous_result.endpoint_toxic_flux_ratio == pytest.approx(1.0)
    assert nonhazardous_result.endpoint_toxic_flux_ratio == pytest.approx(1.0)
    assert hazardous_result.reaction_role is nonhazardous_result.reaction_role is ReactionRole.UNKNOWN
    assert (
        hazardous_result.risk_direction_if_flux_decreases
        is nonhazardous_result.risk_direction_if_flux_decreases
        is RiskDirectionIfFluxDecreases.UNKNOWN
    )


def test_no_adjusted_risk_field_is_produced():
    annotation = _annotation(
        ReactionRole.BIOACTIVATION,
        RiskDirectionIfFluxDecreases.DECREASE,
    )

    result_dict = interpret_competitive_endpoint_flux(0.75, annotation).to_dict()

    assert "adjusted_risk" not in result_dict
    assert "adjusted_risk" not in result_dict["metadata"]


def test_endpoint_toxic_flux_has_no_interaction_engine_integration():
    assert "interaction_engine" not in endpoint_toxic_flux.__dict__
    assert "compute_interaction_matrix" not in endpoint_toxic_flux.__dict__
    assert "competitive_inhibition_flux" not in endpoint_toxic_flux.__dict__
    assert "get_ki" not in endpoint_toxic_flux.__dict__


def test_endpoint_flux_input_uses_registry_lookup_when_annotation_is_absent():
    result = interpret_endpoint_toxic_flux(
        EndpointFluxInput(
            enzyme="CYP2E1",
            substrate="benzene",
            flux_ratio=0.5,
            tissue="liver",
            endpoint=RiskEndpoint.DNA_ADDUCT,
        )
    )

    assert result.annotation_record_id == "reaction_role_benzene_cyp2e1_v2_direction_adjustment"
    assert result.flux_ratio == pytest.approx(0.5)
    assert result.endpoint_toxic_flux_ratio == pytest.approx(2.0)


def test_result_serialization_is_json_friendly():
    annotation = _annotation(
        ReactionRole.BIOACTIVATION,
        RiskDirectionIfFluxDecreases.DECREASE,
        record_id="serialize_case",
    )

    result_dict = interpret_competitive_endpoint_flux(
        0.5,
        annotation,
        enzyme="CYP1A1",
        substrate="benzene",
        endpoint=RiskEndpoint.DNA_ADDUCT,
    ).to_dict()

    assert result_dict["reaction_role"] == "bioactivation"
    assert result_dict["risk_direction_if_flux_decreases"] == "decrease"
    assert result_dict["endpoint"] == "DNA_adduct"
    assert result_dict["metadata"]["precomputed_flux_ratio_required"] is True
