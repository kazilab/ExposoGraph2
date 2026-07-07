import json

import pytest

from ExposoGraph import interaction_engine as engine
from ExposoGraph.effective_burden import (
    GSHBurdenCouplingInput,
    couple_gsh_consumption_to_activation_burden,
)
from ExposoGraph.gsh_redox_capacity import GSHModelVersion
from ExposoGraph.interaction_engine import (
    _interaction_matrix_to_compat_dict,
    compute_interaction_matrix,
)
from ExposoGraph.unified_api import patient_risk_query


def _warning_codes(records):
    codes = set()
    for record in records or []:
        if hasattr(record, "code"):
            codes.add(record.code)
        elif isinstance(record, dict):
            codes.add(record.get("code"))
    return {code for code in codes if code}


def _inhibition_burden(
    *,
    activation_burden_ratio,
    status="mechanism_resolved",
    review_required=False,
):
    return engine._InhibitionBurdenResolution(
        burden_multiplier=1.0,
        activation_burden_ratio=activation_burden_ratio,
        detox_failure_ratio=1.0,
        endpoint_toxic_flux_ratio=1.0,
        status=status,
        review_required=review_required,
    )


def _assert_adjusted_risk_formula(result, carcinogen):
    resolved = result.mechanism_resolved_risks[carcinogen]
    assert resolved.final_mechanism_multiplier == pytest.approx(
        round(
            resolved.induction_multiplier
            * resolved.inhibition_burden_multiplier
            * resolved.matrix_gsh_penalty,
            6,
        )
    )
    assert result.interaction_adjusted_risks[carcinogen] == pytest.approx(
        round(result.individual_risks[carcinogen] * resolved.final_mechanism_multiplier, 3)
    )
    return resolved


def test_gsh_neutral_when_no_gsh_relevant_exposure():
    result = compute_interaction_matrix(
        {"HCA": 2.0},
        enable_induction=False,
        enable_competition=False,
        include_biological_outputs=True,
    )

    assert result.gsh_status.consumption_rate_umol_h_g == pytest.approx(0.0)
    assert result.gsh_status.detox_penalty_multiplier == pytest.approx(1.0)
    assert result.gsh_status.individual_contributions == {}
    hca = _assert_adjusted_risk_formula(result, "HCA")
    assert hca.gsh_pool_penalty == pytest.approx(1.0)
    assert hca.matrix_gsh_penalty == pytest.approx(1.0)


def test_gsh_scaled_by_activation_burden_when_available():
    low = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 1.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={"PAH": _inhibition_burden(activation_burden_ratio=2.0)},
    )
    high = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 1.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={"PAH": _inhibition_burden(activation_burden_ratio=5.0)},
    )

    contribution = high.individual_contributions["PAH_GSTM1"]
    assert contribution["upstream_scaling_source"] == "direct_activation_burden_ratio"
    assert contribution["upstream_activation_scale"] == pytest.approx(5.0)
    assert contribution["gsh_consumption_umol_h_g"] == pytest.approx(5.0)
    assert high.fraction_normal < low.fraction_normal


def test_gsh_d_times_k_fallback_warns():
    coupling = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(
            gsh_relevant=True,
            base_gsh_consumption_load=1.0,
            d_factor=2.0,
            k_factor=3.0,
            tissue="Liver",
        )
    )
    assert coupling.scaling_source == "d_times_k_approximation"
    assert coupling.upstream_activation_burden_ratio == pytest.approx(6.0)
    assert "gsh_d_times_k_fallback_used" in _warning_codes(coupling.warnings)

    matrix_status = engine._compute_matrix_gsh_redox_status(
        {
            "PAH_umol_h_g": {
                "flux_umol_h_g": 1.0,
                "d_factor": 2.0,
                "k_factor": 3.0,
            }
        },
        genotypes={},
        tissue="Liver",
        inhibition_burdens={
            "PAH": _inhibition_burden(
                activation_burden_ratio=1.0,
                status="mechanism_unresolved",
                review_required=True,
            )
        },
    )
    contribution = matrix_status.individual_contributions["PAH_GSTM1"]
    assert contribution["upstream_scaling_source"] == "d_times_k_approximation"
    assert contribution["upstream_scaling_provenance"]["d_factor"] == pytest.approx(2.0)
    assert contribution["upstream_scaling_provenance"]["k_factor"] == pytest.approx(3.0)
    assert "gsh_d_times_k_fallback_used" in _warning_codes(matrix_status.warnings)


def test_gsh_neutral_fallback_warns_when_no_upstream_burden():
    status = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 2.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={
            "PAH": _inhibition_burden(
                activation_burden_ratio=1.0,
                status="mechanism_unresolved",
                review_required=True,
            )
        },
    )

    contribution = status.individual_contributions["PAH_GSTM1"]
    assert contribution["upstream_scaling_source"] == "neutral_missing_upstream_activation"
    assert contribution["upstream_activation_scale"] == pytest.approx(1.0)
    assert "gsh_upstream_activation_missing_neutral" in _warning_codes(status.warnings)


def test_diagnostic_gsh_output_not_authoritative():
    result = compute_interaction_matrix(
        {"PAH": 8.0, "HCA": 2.0, "NDMA": 2.0, "ethanol": 12.0, "acrolein": 8.0},
        genotypes={"GSTM1": "null"},
        lifestyle={"smoking": True, "chronic_alcohol": True},
        include_biological_outputs=True,
    )

    selected = result.competitive_effects["CYP1A1"].substrates["PhIP"].biological_output
    assert selected["selected_authoritative_effect"] is True
    assert selected["gsh_capacity_effect"]["diagnostic_only"] is True
    assert selected["gsh_capacity_effect"]["included_in_authoritative_adjusted_risk"] is False
    assert selected["effective_burden"]["includes_diagnostic_gsh_capacity"] is False
    _assert_adjusted_risk_formula(result, "HCA")


def test_gsh_penalty_applied_once():
    result = compute_interaction_matrix(
        {"PAH": 8.0, "acrolein": 8.0},
        genotypes={"GSTM1": "null"},
        include_biological_outputs=True,
    )

    pah = _assert_adjusted_risk_formula(result, "PAH")
    assert pah.matrix_gsh_penalty == pytest.approx(
        pah.gsh_pool_penalty * pah.susceptibility_modifier
    )
    assert pah.provenance["gsh"]["matrix_gsh_penalty_applied_once"] is True
    assert pah.provenance["gsh"]["diagnostic_gsh_capacity_included"] is False


def test_gstm1_null_susceptibility_not_double_counted():
    result = compute_interaction_matrix(
        {"PAH": 8.0, "acrolein": 8.0},
        genotypes={"GSTM1": "null"},
        lifestyle={"smoking": True},
        include_biological_outputs=True,
    )

    pah = _assert_adjusted_risk_formula(result, "PAH")
    assert pah.susceptibility_modifier > 1.0
    assert pah.susceptibility_applied_in == "matrix_gsh_penalty"
    assert pah.final_mechanism_multiplier == pytest.approx(
        round(
            pah.induction_multiplier
            * pah.inhibition_burden_multiplier
            * pah.gsh_pool_penalty
            * pah.susceptibility_modifier,
            6,
        )
    )


def test_high_mixed_exposure_reports_redox_capacity_and_warnings():
    result = compute_interaction_matrix(
        {
            "PAH": 8.0,
            "acrolein": 8.0,
            "chromium_VI": 4.0,
            "arsenic": 4.0,
            "ethanol": 12.0,
        },
        genotypes={"GSTM1": "null"},
        include_biological_outputs=True,
    )

    assert result.gsh_status.model_version == GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY.value
    assert result.gsh_status.redox_capacity_ratio is not None
    assert result.gsh_status.detox_penalty_multiplier > 1.0
    assert result.gsh_status.individual_contributions
    assert _warning_codes(result.gsh_status.warnings)
    json.dumps(_interaction_matrix_to_compat_dict(result)["gsh_status"], allow_nan=False)


def test_adjusted_risk_output_contains_gsh_provenance():
    result = compute_interaction_matrix(
        {"PAH": 8.0, "acrolein": 8.0},
        genotypes={"GSTM1": "null"},
        include_biological_outputs=True,
    )

    pah = _assert_adjusted_risk_formula(result, "PAH")
    gsh = pah.provenance["gsh"]
    assert pah.provenance["gsh_source"] == "gsh_redox_capacity.compute_gsh_redox_capacity"
    assert gsh["model_version"] == GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY.value
    assert gsh["redox_capacity_ratio"] == result.gsh_status.redox_capacity_ratio
    assert gsh["detox_penalty_multiplier"] == result.gsh_status.detox_penalty_multiplier
    assert "warnings" in gsh
    json.dumps(pah.to_dict(), allow_nan=False)


def test_platform_integration_module5_gsh_case():
    result = compute_interaction_matrix(
        {"PAH": 4.0, "HCA": 1.0, "acrolein": 4.0},
        genotypes={"GSTM1": "NM"},
        include_biological_outputs=True,
    )
    payload = _interaction_matrix_to_compat_dict(result)

    assert payload["module5_model_card"]["gsh_model_version"] == (
        GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY.value
    )
    assert payload["gsh_status"]["model_version"] == GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY.value
    assert payload["mechanism_resolved_risks"]["PAH"]["provenance"]["gsh"]
    json.dumps(payload, allow_nan=False)

    profile = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Liver",
        include_tissue_report=False,
    )
    integration = profile.biological_output_integration
    assert integration["module5_model_card"]["gsh_model_version"] == (
        GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY.value
    )
    assert "mechanism_resolved_risks" in integration
    json.dumps(integration, allow_nan=False)
