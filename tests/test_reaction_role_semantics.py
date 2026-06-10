import pytest

from ExposoGraph.interaction_schema import (
    ReactionRole,
    ReleaseTarget,
    RiskDirectionIfFluxDecreases,
    SMEReviewStatus,
)
from ExposoGraph.reaction_role_semantics import (
    ReactionRoleAnnotation,
    get_default_reaction_role_registry,
    get_spyros_sme_reaction_role_records,
    interpret_competitive_flux_ratio,
)


def test_required_spyros_records_are_present_with_statuses():
    records = {record.record_id: record for record in get_spyros_sme_reaction_role_records()}

    required = {
        "spyros_benzene_cyp2e1_v2_flip",
        "spyros_benzene_cyp2f1_bone_marrow_flip",
        "spyros_benzene_cyp2a13_bone_marrow_flip",
        "spyros_benzene_cyp1a1_no_default_flip",
        "spyros_ndma_cyp2e1_no_flip",
        "spyros_vinyl_chloride_cyp2e1_no_flip",
        "spyros_hca_cyp1a1_no_flip",
        "spyros_hca_cyp1a2_no_flip",
        "spyros_tce_cyp2e1_candidate_pending",
    }

    assert required.issubset(records)
    assert records["spyros_benzene_cyp2e1_v2_flip"].release_target is ReleaseTarget.V2_0
    assert records["spyros_tce_cyp2e1_candidate_pending"].active is False
    assert (
        records["spyros_tce_cyp2e1_candidate_pending"].review_status
        is SMEReviewStatus.PENDING_TEAM_AGREEMENT
    )


def test_benzene_cyp2e1_record_expresses_sign_flip():
    registry = get_default_reaction_role_registry()
    record = registry.lookup("cyp2e1", "BENZENE")

    assert record.record_id == "spyros_benzene_cyp2e1_v2_flip"
    assert record.reaction_role is ReactionRole.DETOXIFICATION
    assert record.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.INCREASE

    interpretation = interpret_competitive_flux_ratio(0.5, record)
    assert interpretation.burden_multiplier == pytest.approx(2.0)


def test_benzene_cyp2f1_and_cyp2a13_flip_only_in_bone_marrow_context():
    registry = get_default_reaction_role_registry()

    for enzyme in ("CYP2F1", "CYP2A13"):
        global_record = registry.lookup(enzyme, "benzene")
        marrow_record = registry.lookup(enzyme, "benzene", tissue="bone marrow")
        endpoint_record = registry.lookup(enzyme, "benzene", endpoint="hematopoietic benzene toxicity")

        assert global_record.risk_direction_if_flux_decreases is not RiskDirectionIfFluxDecreases.INCREASE
        assert interpret_competitive_flux_ratio(0.5, global_record).burden_multiplier == pytest.approx(1.0)
        assert marrow_record.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.INCREASE
        assert endpoint_record.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.INCREASE
        assert interpret_competitive_flux_ratio(0.5, marrow_record).burden_multiplier == pytest.approx(2.0)


def test_no_flip_records_map_flux_decrease_to_lower_burden():
    registry = get_default_reaction_role_registry()
    cases = [
        ("CYP1A1", "benzene", "spyros_benzene_cyp1a1_no_default_flip"),
        ("CYP2E1", "NDMA", "spyros_ndma_cyp2e1_no_flip"),
        ("CYP2E1", "vinyl chloride", "spyros_vinyl_chloride_cyp2e1_no_flip"),
        ("CYP1A1", "HCA", "spyros_hca_cyp1a1_no_flip"),
        ("CYP1A2", "heterocyclic amine", "spyros_hca_cyp1a2_no_flip"),
    ]

    for enzyme, substrate, record_id in cases:
        record = registry.lookup(enzyme, substrate)
        interpretation = interpret_competitive_flux_ratio(0.5, record)

        assert record.record_id == record_id
        assert record.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.DECREASE
        assert interpretation.burden_multiplier == pytest.approx(0.5)


def test_tce_cyp2e1_is_pending_candidate_not_active_curated_behavior():
    registry = get_default_reaction_role_registry()
    record = registry.lookup("CYP2E1", "TCE")
    interpretation = interpret_competitive_flux_ratio(0.5, record)

    assert record.record_id == "spyros_tce_cyp2e1_candidate_pending"
    assert record.review_status is SMEReviewStatus.PENDING_TEAM_AGREEMENT
    assert record.active is False
    assert interpretation.burden_multiplier == pytest.approx(1.0)
    assert any(warning.code == "reaction_role_inactive" for warning in interpretation.warnings)


def test_broader_detox_arm_records_are_deferred_or_absent_from_active_2_0_behavior():
    registry = get_default_reaction_role_registry()

    for enzyme in ("NAT2", "ALDH", "GST", "EPHX1"):
        deferred = registry.lookup(enzyme, "broader detox-arm cases")
        assert deferred.review_status is SMEReviewStatus.DEFERRED_3_0
        assert deferred.release_target is ReleaseTarget.V3_0
        assert deferred.active is False

    nat2_specific = registry.lookup("NAT2", "N_OH_PhIP")
    assert nat2_specific.record_id == "unknown_reaction_role"
    assert nat2_specific.active is False


def test_unknown_reaction_role_returns_neutral_interpretation_with_warning():
    registry = get_default_reaction_role_registry()
    record = registry.lookup("CYP9Z9", "unknown substrate")
    interpretation = interpret_competitive_flux_ratio(0.25, record)

    assert record.reaction_role is ReactionRole.UNKNOWN
    assert record.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.UNKNOWN
    assert interpretation.burden_multiplier == pytest.approx(1.0)
    assert any(warning.code == "reaction_role_unmatched" for warning in interpretation.warnings)


def test_product_carcinogenic_metadata_alone_does_not_drive_sign_rule():
    hazardous = ReactionRoleAnnotation(
        enzyme="CYPX",
        substrate="same",
        metadata={"product_carcinogenic": True},
    )
    nonhazardous = ReactionRoleAnnotation(
        enzyme="CYPX",
        substrate="same",
        metadata={"product_carcinogenic": False},
    )

    assert hazardous.reaction_role is nonhazardous.reaction_role is ReactionRole.UNKNOWN
    assert (
        hazardous.risk_direction_if_flux_decreases
        is nonhazardous.risk_direction_if_flux_decreases
        is RiskDirectionIfFluxDecreases.UNKNOWN
    )
    assert interpret_competitive_flux_ratio(0.5, hazardous).burden_multiplier == pytest.approx(1.0)
    assert interpret_competitive_flux_ratio(0.5, nonhazardous).burden_multiplier == pytest.approx(1.0)


def test_interpretation_helper_maps_decrease_and_increase_directions_differently():
    lower_burden = ReactionRoleAnnotation(
        enzyme="E1",
        substrate="S",
        reaction_role=ReactionRole.BIOACTIVATION,
        risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.DECREASE,
    )
    higher_burden = ReactionRoleAnnotation(
        enzyme="E2",
        substrate="S",
        reaction_role=ReactionRole.DETOXIFICATION,
        risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.INCREASE,
    )

    assert interpret_competitive_flux_ratio(0.25, lower_burden).burden_multiplier == pytest.approx(0.25)
    assert interpret_competitive_flux_ratio(0.25, higher_burden).burden_multiplier == pytest.approx(4.0)


def test_registry_lookup_is_deterministic_and_case_insensitive():
    registry = get_default_reaction_role_registry()
    ids = [registry.lookup("cyp2e1", "benzene").record_id for _ in range(3)]

    assert ids == ["spyros_benzene_cyp2e1_v2_flip"] * 3
    assert registry.lookup("CYP2E1", "benzene").record_id == registry.lookup("cyp2e1", "BENZENE").record_id
    assert registry.lookup("CYP2E1", "vinyl_chloride").record_id == registry.lookup(
        "cyp2e1",
        "vinyl chloride",
    ).record_id


def test_reaction_role_semantics_has_no_interaction_engine_integration():
    import ExposoGraph.reaction_role_semantics as semantics

    assert "interaction_engine" not in semantics.__dict__
    assert "compute_interaction_matrix" not in semantics.__dict__
