import json
import math

import pytest

from ExposoGraph import interaction_engine as engine
from ExposoGraph.interaction_engine import (
    CompetitiveInhibitionResult,
    SubstrateFluxChange,
    _interaction_matrix_to_compat_dict,
    _resolve_endpoint_inhibition_burden,
    compute_interaction_matrix,
)
from ExposoGraph.gsh_redox_capacity import GSHModelVersion
from ExposoGraph.unified_api import _build_biological_output_integration, patient_risk_query


def _finite_walk(value):
    if isinstance(value, float):
        assert math.isfinite(value)
    elif isinstance(value, dict):
        for item in value.values():
            _finite_walk(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _finite_walk(item)


def _principal_numbers(result):
    return {
        "individual_risks": dict(result.individual_risks),
        "interaction_adjusted_risks": dict(result.interaction_adjusted_risks),
        "synergy_matrix": dict(result.synergy_matrix),
        "total_independent_risk": result.total_independent_risk,
        "total_interaction_risk": result.total_interaction_risk,
        "interaction_factor": result.interaction_factor,
        "gsh_status": {
            "baseline_gsh_mM": result.gsh_status.baseline_gsh_mM,
            "steady_state_gsh_mM": result.gsh_status.steady_state_gsh_mM,
            "fraction_normal": result.gsh_status.fraction_normal,
            "consumption_rate_umol_h_g": result.gsh_status.consumption_rate_umol_h_g,
            "synthesis_rate_umol_h_g": result.gsh_status.synthesis_rate_umol_h_g,
            "net_rate_umol_h_g": result.gsh_status.net_rate_umol_h_g,
            "tipping_point_multiplier": result.gsh_status.tipping_point_multiplier,
            "time_to_depletion_h": result.gsh_status.time_to_depletion_h,
            "individual_contributions": deepcopy_dict(result.gsh_status.individual_contributions),
            "model_version": result.gsh_status.model_version,
            "redox_capacity_ratio": result.gsh_status.redox_capacity_ratio,
            "detox_penalty_multiplier": result.gsh_status.detox_penalty_multiplier,
            "warnings": deepcopy_dict(result.gsh_status.warnings),
            "metadata": deepcopy_dict(result.gsh_status.metadata),
        },
        "induction_effects": dict(result.induction_effects.enzyme_folds),
        "competitive_effects": {
            enzyme: {
                substrate: {
                    "single_flux": flux.single_flux,
                    "competitive_flux": flux.competitive_flux,
                    "flux_change_fraction": flux.flux_change_fraction,
                    "inhibition_term": flux.inhibition_term,
                    "activated_product_flux": flux.activated_product_flux,
                    "Km_uM": flux.Km_uM,
                    "concentration_uM": flux.concentration_uM,
                    "baseline_flux": flux.baseline_flux,
                    "kinetic_modifier": flux.kinetic_modifier,
                    "modified_flux": flux.modified_flux,
                }
                for substrate, flux in enzyme_result.substrates.items()
            }
            for enzyme, enzyme_result in result.competitive_effects.items()
        },
    }


def deepcopy_dict(value):
    return json.loads(json.dumps(value))


def _inhibition_burden(
    *,
    activation_burden_ratio: float,
    status: str = "mechanism_resolved",
    review_required: bool = False,
):
    return engine._InhibitionBurdenResolution(
        burden_multiplier=1.0,
        activation_burden_ratio=activation_burden_ratio,
        detox_failure_ratio=1.0,
        endpoint_toxic_flux_ratio=1.0,
        status=status,
        review_required=review_required,
    )


def _flux(
    *,
    competitive_flux: float,
    status: str = "resolved_direct",
    modifier_applied_once: bool = True,
) -> SubstrateFluxChange:
    return SubstrateFluxChange(
        single_flux=1.0,
        competitive_flux=competitive_flux,
        flux_change_fraction=round(competitive_flux - 1.0, 4),
        inhibition_term=1.0,
        activated_product_flux=0.0,
        Km_uM=1.0,
        concentration_uM=1.0,
        product="phenol",
        product_carcinogenic=False,
        baseline_flux=1.0,
        kinetic_modifier=competitive_flux if modifier_applied_once else None,
        modified_flux=competitive_flux if modifier_applied_once else None,
        kinetic_resolution_status=status,
        inhibition_mode="competitive",
        centralized_resolver_used=True,
        modifier_applied_once=modifier_applied_once,
    )


def _selected_biological_output(result, enzyme, substrate):
    output = result.competitive_effects[enzyme].substrates[substrate].biological_output
    assert output is not None
    assert output["selected_authoritative_effect"] is True
    return output


def _inhibition_provenance(result, carcinogen):
    return result.mechanism_resolved_risks[carcinogen].provenance["inhibition"]


def test_pulmonary_benzene_nnk_selects_resolved_cyp2a13_candidate():
    exposure = {"benzene": 1.0, "NNK": 1.0}
    result = compute_interaction_matrix(
        exposure,
        tissue="Lung",
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )
    without_outputs = compute_interaction_matrix(
        exposure,
        tissue="Lung",
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=False,
    )

    cyp2a13 = result.competitive_effects["CYP2A13"].substrates["benzene"]
    cyp2f1 = result.competitive_effects["CYP2F1"].substrates["benzene"]
    selected = result.mechanism_resolved_risks["benzene"]
    selected_output = cyp2a13.biological_output
    diagnostic_output = cyp2f1.biological_output

    assert engine._kinetic_mechanism_state(cyp2a13.kinetic_resolution_status) == "mechanism_resolved"
    assert engine._kinetic_mechanism_state(cyp2f1.kinetic_resolution_status) == "mechanism_absent"
    assert cyp2a13.flux_change_fraction < 0.0
    assert cyp2a13.kinetic_modifier != pytest.approx(1.0)
    assert cyp2f1.flux_change_fraction == pytest.approx(0.0)

    inhibition = selected.provenance["inhibition"]
    assert inhibition["enzyme"] == "CYP2A13"
    assert selected.inhibition_status == "mechanism_resolved"
    assert selected_output is not None
    assert selected_output["selected_authoritative_effect"] is True
    assert selected_output["kinetic_effect"]["mechanism_state"] == "mechanism_resolved"
    assert selected_output["kinetic_effect"]["provenance"]["enzyme"] == "CYP2A13"
    assert selected_output["interpretation_substrate"] == "benzene"

    assert diagnostic_output is not None
    assert diagnostic_output["selected_authoritative_effect"] is False
    assert diagnostic_output["kinetic_effect"]["mechanism_state"] == "mechanism_absent"
    assert selected.inhibition_burden_multiplier == pytest.approx(1.0)
    assert selected.review_required is True
    assert selected_output["endpoint_toxic_flux"]["burden_multiplier"] == pytest.approx(1.0)
    assert selected_output["reaction_role_interpretation"]["review_required"] is True
    assert selected_output["reaction_role_interpretation"]["warnings"]
    assert _principal_numbers(result) == _principal_numbers(without_outputs)


def _stub_recursive_attribution(monkeypatch):
    monkeypatch.setattr(
        engine,
        "_compute_live_mechanism_attribution",
        lambda *args, **kwargs: {"stubbed_for_selected_serialization_test": True},
    )


def _spy_endpoint_and_effective_calls(monkeypatch):
    endpoint_calls = []
    effective_calls = []
    real_endpoint = engine.interpret_endpoint_toxic_flux
    real_effective = engine.compute_effective_carcinogenic_burden

    def endpoint_spy(endpoint_input, *args, **kwargs):
        endpoint_calls.append(endpoint_input)
        return real_endpoint(endpoint_input, *args, **kwargs)

    def effective_spy(*args, **kwargs):
        request = args[0] if args else kwargs.get("burden_input")
        endpoint_result = getattr(request, "endpoint_toxic_flux_result", None)
        if endpoint_result is None:
            endpoint_result = kwargs.get("endpoint_toxic_flux_result")
        effective_calls.append(endpoint_result)
        return real_effective(*args, **kwargs)

    monkeypatch.setattr(engine, "interpret_endpoint_toxic_flux", endpoint_spy)
    monkeypatch.setattr(engine, "compute_effective_carcinogenic_burden", effective_spy)
    return endpoint_calls, effective_calls


def _spy_role_registry_lookup(monkeypatch):
    lookup_calls = []
    real_factory = engine.get_default_reaction_role_registry

    class LookupSpy:
        def __init__(self, registry):
            self._registry = registry

        def lookup(self, enzyme, substrate, *args, **kwargs):
            tissue = kwargs.get("tissue")
            if tissue is None and args:
                tissue = args[0]
            endpoint = kwargs.get("endpoint")
            lookup_calls.append(
                {
                    "enzyme": enzyme,
                    "substrate": substrate,
                    "tissue": tissue,
                    "endpoint": endpoint,
                }
            )
            return self._registry.lookup(enzyme, substrate, *args, **kwargs)

        def __getattr__(self, name):
            return getattr(self._registry, name)

    monkeypatch.setattr(engine, "get_default_reaction_role_registry", lambda: LookupSpy(real_factory()))
    return lookup_calls


def _selected_role_lookup_calls(lookup_calls, enzyme, interpretation_substrate, tissue="Liver"):
    return [
        call
        for call in lookup_calls
        if call["enzyme"] == enzyme
        and call["substrate"] == interpretation_substrate
        and call["tissue"] == tissue
    ]


def _selected_endpoint_calls(endpoint_calls, enzyme, interpretation_substrate):
    return [
        call
        for call in endpoint_calls
        if call.enzyme == enzyme and call.substrate == interpretation_substrate
    ]


def _selected_effective_calls(effective_calls, enzyme, interpretation_substrate):
    return [
        endpoint_result
        for endpoint_result in effective_calls
        if endpoint_result is not None
        and endpoint_result.enzyme == enzyme
        and endpoint_result.substrate == interpretation_substrate
    ]


def _force_benzene_unresolved_competition(monkeypatch):
    real_competitive_inhibition_flux = engine.competitive_inhibition_flux

    def fake_competitive_inhibition_flux(enzyme, substrates, *args, **kwargs):
        if enzyme == "CYP2E1":
            return CompetitiveInhibitionResult(
                "CYP2E1",
                {
                    "benzene": _flux(
                        competitive_flux=1.0,
                        status="review_required",
                        modifier_applied_once=False,
                    )
                },
            )
        return real_competitive_inhibition_flux(enzyme, substrates, *args, **kwargs)

    monkeypatch.setattr(engine, "competitive_inhibition_flux", fake_competitive_inhibition_flux)


def test_detoxification_inhibition_increases_adjusted_risk():
    result = compute_interaction_matrix(
        {"benzene": 1.0, "ethanol": 10.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    flux = result.competitive_effects["CYP2E1"].substrates["benzene"]
    resolved = result.mechanism_resolved_risks["benzene"]

    assert flux.flux_change_fraction < 0.0
    assert resolved.activation_burden_ratio == pytest.approx(1.0)
    assert resolved.detox_failure_ratio > 1.0
    assert resolved.inhibition_burden_multiplier == pytest.approx(
        resolved.detox_failure_ratio
    )
    assert result.interaction_adjusted_risks["benzene"] > result.individual_risks["benzene"]


@pytest.mark.parametrize(
    ("exposure", "carcinogen", "enzyme", "substrate"),
    [
        ({"benzene": 1.0, "ethanol": 10.0}, "benzene", "CYP2E1", "benzene"),
        ({"NDMA": 1.0, "ethanol": 10.0}, "NDMA", "CYP2E1", "NDMA"),
    ],
)
def test_selected_effect_serialization_reuses_authoritative_resolution(
    exposure,
    carcinogen,
    enzyme,
    substrate,
):
    result = compute_interaction_matrix(
        exposure,
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    resolved = result.mechanism_resolved_risks[carcinogen]
    inhibition = _inhibition_provenance(result, carcinogen)
    selected = _selected_biological_output(result, enzyme, substrate)

    assert selected["interpretation_substrate"] == inhibition["interpretation_substrate"]
    assert selected["kinetic_effect"]["mechanism_state"] == resolved.inhibition_status
    assert selected["endpoint_toxic_flux"]["activation_burden_ratio"] == pytest.approx(
        resolved.activation_burden_ratio
    )
    assert selected["endpoint_toxic_flux"]["detox_failure_ratio"] == pytest.approx(
        resolved.detox_failure_ratio
    )
    assert selected["effective_burden"]["effective_carcinogenic_burden_ratio"] == pytest.approx(
        resolved.inhibition_burden_multiplier
    )
    assert selected["endpoint_toxic_flux"]["warnings"] == inhibition["endpoint_toxic_flux"][
        "warnings"
    ]
    assert selected["effective_burden"]["includes_diagnostic_gsh_capacity"] is False


@pytest.mark.parametrize(
    (
        "exposure",
        "enzyme",
        "flux_substrate",
        "interpretation_substrate",
        "expected_status",
        "force_unresolved",
    ),
    [
        (
            {"benzene": 1.0, "ethanol": 10.0},
            "CYP2E1",
            "benzene",
            "benzene",
            "mechanism_resolved",
            False,
        ),
        (
            {"HCA": 1.0, "PAH": 0.0},
            "CYP1A1",
            "PhIP",
            "HCA",
            "mechanism_absent",
            False,
        ),
        (
            {"benzene": 1.0, "ethanol": 10.0},
            "CYP2E1",
            "benzene",
            "benzene",
            "mechanism_unresolved",
            True,
        ),
    ],
)
def test_selected_role_lookup_occurs_once_for_each_selected_status(
    monkeypatch,
    exposure,
    enzyme,
    flux_substrate,
    interpretation_substrate,
    expected_status,
    force_unresolved,
):
    _stub_recursive_attribution(monkeypatch)
    lookup_calls = _spy_role_registry_lookup(monkeypatch)
    if force_unresolved:
        _force_benzene_unresolved_competition(monkeypatch)

    result = compute_interaction_matrix(
        exposure,
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    selected = _selected_biological_output(result, enzyme, flux_substrate)
    assert selected["kinetic_effect"]["mechanism_state"] == expected_status
    assert selected["interpretation_substrate"] == interpretation_substrate
    assert len(_selected_role_lookup_calls(lookup_calls, enzyme, interpretation_substrate)) == 1


def test_selected_role_lookup_not_called_when_no_effect_is_selected(monkeypatch):
    _stub_recursive_attribution(monkeypatch)
    lookup_calls = _spy_role_registry_lookup(monkeypatch)

    no_effect = compute_interaction_matrix(
        {"benzene": 1.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )
    disabled = compute_interaction_matrix(
        {"benzene": 1.0, "ethanol": 10.0},
        enable_induction=False,
        enable_competition=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    assert no_effect.mechanism_resolved_risks["benzene"].inhibition_status == "mechanism_absent"
    assert disabled.mechanism_resolved_risks["benzene"].inhibition_status == "mechanism_disabled"
    assert _selected_role_lookup_calls(lookup_calls, "CYP2E1", "benzene") == []


@pytest.mark.parametrize(
    ("exposure", "carcinogen", "enzyme", "flux_substrate", "interpretation_substrate"),
    [
        ({"benzene": 1.0, "ethanol": 10.0}, "benzene", "CYP2E1", "benzene", "benzene"),
        ({"PAH": 3.0, "HCA": 2.0}, "HCA", "CYP1A1", "PhIP", "HCA"),
    ],
)
def test_selected_endpoint_and_effective_calls_are_exclusive(
    monkeypatch,
    exposure,
    carcinogen,
    enzyme,
    flux_substrate,
    interpretation_substrate,
):
    _stub_recursive_attribution(monkeypatch)
    endpoint_calls, effective_calls = _spy_endpoint_and_effective_calls(monkeypatch)

    result = compute_interaction_matrix(
        exposure,
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    selected = _selected_biological_output(result, enzyme, flux_substrate)
    resolved = result.mechanism_resolved_risks[carcinogen]

    assert selected["interpretation_substrate"] == interpretation_substrate
    assert selected["effective_burden"]["effective_carcinogenic_burden_ratio"] == pytest.approx(
        resolved.inhibition_burden_multiplier
    )
    assert len(_selected_endpoint_calls(endpoint_calls, enzyme, interpretation_substrate)) == 1
    assert len(_selected_effective_calls(effective_calls, enzyme, interpretation_substrate)) == 1


def test_selected_mechanism_absent_output_is_neutral_selected_and_exclusive(monkeypatch):
    _stub_recursive_attribution(monkeypatch)
    endpoint_calls, effective_calls = _spy_endpoint_and_effective_calls(monkeypatch)

    result = compute_interaction_matrix(
        {"HCA": 1.0, "PAH": 0.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    resolved = result.mechanism_resolved_risks["HCA"]
    selected = _selected_biological_output(result, "CYP1A1", "PhIP")

    assert resolved.inhibition_status == "mechanism_absent"
    assert resolved.inhibition_burden_multiplier == pytest.approx(1.0)
    assert resolved.provenance["inhibition"]["flux_substrate"] == "PhIP"
    assert resolved.provenance["inhibition"]["interpretation_substrate"] == "HCA"
    assert resolved.provenance["inhibition"]["reaction_role"] == "bioactivation"
    assert resolved.provenance["inhibition"]["annotation_record_id"] == "spyros_hca_cyp1a1_no_flip"
    assert selected["diagnostic_role"] == "selected_authoritative_inhibition_effect"
    assert selected["interpretation_substrate"] == "HCA"
    assert selected["reaction_role_interpretation"]["role"] == "bioactivation"
    assert (
        selected["reaction_role_interpretation"]["risk_direction_if_flux_decreases"]
        == "decrease"
    )
    assert (
        selected["reaction_role_interpretation"]["annotation_record_id"]
        == "spyros_hca_cyp1a1_no_flip"
    )
    assert selected["kinetic_effect"]["mechanism_state"] == "mechanism_absent"
    assert selected["endpoint_toxic_flux"]["selected_authoritative_effect"] is True
    assert selected["endpoint_toxic_flux"]["diagnostic_only"] is False
    assert selected["endpoint_toxic_flux"]["status"] == "mechanism_absent"
    assert selected["endpoint_toxic_flux"]["reaction_role"] == "bioactivation"
    assert selected["endpoint_toxic_flux"]["risk_direction_if_flux_decreases"] == "decrease"
    assert (
        selected["endpoint_toxic_flux"]["annotation_record_id"]
        == "spyros_hca_cyp1a1_no_flip"
    )
    assert selected["endpoint_toxic_flux"]["endpoint_toxic_flux_ratio"] == pytest.approx(1.0)
    assert selected["effective_burden"]["selected_authoritative_effect"] is True
    assert selected["effective_burden"]["diagnostic_only"] is False
    assert selected["effective_burden"]["reaction_role"] == "bioactivation"
    assert (
        selected["effective_burden"]["annotation_record_id"]
        == "spyros_hca_cyp1a1_no_flip"
    )
    assert selected["effective_burden"]["effective_carcinogenic_burden_ratio"] == pytest.approx(1.0)
    assert selected["effective_burden"]["includes_diagnostic_gsh_capacity"] is False
    assert not _selected_endpoint_calls(endpoint_calls, "CYP1A1", "HCA")
    assert not _selected_effective_calls(effective_calls, "CYP1A1", "HCA")


def test_selected_unresolved_output_is_neutral_reviewed_and_exclusive(monkeypatch):
    _stub_recursive_attribution(monkeypatch)
    _force_benzene_unresolved_competition(monkeypatch)
    endpoint_calls, effective_calls = _spy_endpoint_and_effective_calls(monkeypatch)

    result = compute_interaction_matrix(
        {"benzene": 1.0, "ethanol": 10.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    resolved = result.mechanism_resolved_risks["benzene"]
    selected = _selected_biological_output(result, "CYP2E1", "benzene")

    assert resolved.inhibition_status == "mechanism_unresolved"
    assert resolved.inhibition_burden_multiplier == pytest.approx(1.0)
    assert resolved.review_required is True
    assert selected["reaction_role_interpretation"]["role"] == "detoxification"
    assert selected["reaction_role_interpretation"]["review_required"] is True
    assert selected["endpoint_toxic_flux"]["endpoint_toxic_flux_ratio"] == pytest.approx(1.0)
    assert selected["endpoint_toxic_flux"]["review_required"] is True
    assert selected["effective_burden"]["effective_carcinogenic_burden_ratio"] == pytest.approx(1.0)
    assert selected["effective_burden"]["review_required"] is True
    assert not _selected_endpoint_calls(endpoint_calls, "CYP2E1", "benzene")
    assert not _selected_effective_calls(effective_calls, "CYP2E1", "benzene")


def test_review_required_role_is_preserved_in_neutral_selected_state(monkeypatch):
    _stub_recursive_attribution(monkeypatch)
    endpoint_calls, effective_calls = _spy_endpoint_and_effective_calls(monkeypatch)

    result = compute_interaction_matrix(
        {"benzene": 1.0},
        tissue="Lung",
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    resolved = result.mechanism_resolved_risks["benzene"]
    selected = _selected_biological_output(result, "CYP2A13", "benzene")
    role = selected["reaction_role_interpretation"]

    assert resolved.inhibition_status == "mechanism_absent"
    assert resolved.inhibition_burden_multiplier == pytest.approx(1.0)
    assert resolved.review_required is True
    assert resolved.provenance["inhibition"]["annotation_record_id"] == (
        "spyros_benzene_cyp2a13_outside_context_unknown"
    )
    assert role["role"] == "unknown"
    assert role["review_required"] is True
    assert role["annotation_record_id"] == "spyros_benzene_cyp2a13_outside_context_unknown"
    assert role["warnings"] or role["sme_notes"]
    assert selected["endpoint_toxic_flux"]["reaction_role"] == "unknown"
    assert selected["endpoint_toxic_flux"]["annotation_record_id"] == role["annotation_record_id"]
    assert selected["endpoint_toxic_flux"]["review_required"] is True
    assert selected["endpoint_toxic_flux"]["endpoint_toxic_flux_ratio"] == pytest.approx(1.0)
    assert selected["effective_burden"]["annotation_record_id"] == role["annotation_record_id"]
    assert selected["effective_burden"]["review_required"] is True
    assert selected["effective_burden"]["effective_carcinogenic_burden_ratio"] == pytest.approx(1.0)
    assert not _selected_endpoint_calls(endpoint_calls, "CYP2A13", "benzene")
    assert not _selected_effective_calls(effective_calls, "CYP2A13", "benzene")


def test_bioactivation_inhibition_decreases_adjusted_risk():
    result = compute_interaction_matrix(
        {"NDMA": 1.0, "ethanol": 10.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    flux = result.competitive_effects["CYP2E1"].substrates["NDMA"]
    resolved = result.mechanism_resolved_risks["NDMA"]

    assert flux.flux_change_fraction < 0.0
    assert resolved.activation_burden_ratio < 1.0
    assert resolved.detox_failure_ratio == pytest.approx(1.0)
    assert resolved.inhibition_burden_multiplier == pytest.approx(
        resolved.activation_burden_ratio
    )
    assert result.interaction_adjusted_risks["NDMA"] < result.individual_risks["NDMA"]


def test_hca_selected_phip_flux_uses_hca_role_mapping_in_output():
    result = compute_interaction_matrix(
        {"PAH": 3.0, "HCA": 2.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    flux = result.competitive_effects["CYP1A1"].substrates["PhIP"]
    resolved = result.mechanism_resolved_risks["HCA"]
    selected = _selected_biological_output(result, "CYP1A1", "PhIP")

    assert flux.flux_change_fraction < 0.0
    assert selected["interpretation_substrate"] == "HCA"
    assert selected["endpoint_toxic_flux"]["substrate"] == "HCA"
    assert selected["endpoint_toxic_flux"]["reaction_role"] == "bioactivation"
    assert (
        selected["endpoint_toxic_flux"]["annotation_record_id"]
        == "spyros_hca_cyp1a1_no_flip"
    )
    assert selected["reaction_role_interpretation"]["role"] == "bioactivation"
    assert (
        selected["reaction_role_interpretation"]["annotation_record_id"]
        == "spyros_hca_cyp1a1_no_flip"
    )
    assert selected["effective_burden"]["effective_carcinogenic_burden_ratio"] == pytest.approx(
        resolved.inhibition_burden_multiplier
    )
    assert resolved.inhibition_burden_multiplier < 1.0
    assert resolved.inhibition_burden_multiplier != pytest.approx(1.0)


def test_unknown_and_unresolved_inhibition_are_neutral_and_warned():
    unknown = _resolve_endpoint_inhibition_burden(
        "benzene",
        {
            "CYP2A13": CompetitiveInhibitionResult(
                "CYP2A13",
                {"benzene": _flux(competitive_flux=0.5)},
            )
        },
        tissue="Lung",
        enable_competition=True,
    )
    unresolved = _resolve_endpoint_inhibition_burden(
        "benzene",
        {
            "CYP2E1": CompetitiveInhibitionResult(
                "CYP2E1",
                {
                    "benzene": _flux(
                        competitive_flux=1.0,
                        status="review_required",
                        modifier_applied_once=False,
                    )
                },
            )
        },
        tissue="Liver",
        enable_competition=True,
    )

    assert unknown.burden_multiplier == pytest.approx(1.0)
    assert unknown.review_required is True
    assert {"context_missing", "reaction_role_unknown"}.issubset(unknown.warnings)

    assert unresolved.burden_multiplier == pytest.approx(1.0)
    assert unresolved.review_required is True
    assert unresolved.status == "mechanism_unresolved"


def test_no_competition_and_disabled_competition_are_neutral_without_changing_other_factors():
    no_competition = compute_interaction_matrix(
        {"benzene": 1.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )
    no_comp_resolved = no_competition.mechanism_resolved_risks["benzene"]
    assert no_comp_resolved.inhibition_status == "mechanism_absent"
    assert no_comp_resolved.inhibition_burden_multiplier == pytest.approx(1.0)
    assert no_competition.interaction_adjusted_risks["benzene"] == pytest.approx(
        no_competition.individual_risks["benzene"]
    )

    exposure = {"PAH": 8.0, "HCA": 2.0, "NDMA": 2.0, "ethanol": 12.0, "acrolein": 8.0}
    lifestyle = {"smoking": True, "chronic_alcohol": True}
    genotypes = {"GSTM1": "null"}
    enabled = compute_interaction_matrix(
        exposure,
        genotypes=genotypes,
        lifestyle=lifestyle,
        include_biological_outputs=True,
    )
    disabled = compute_interaction_matrix(
        exposure,
        genotypes=genotypes,
        lifestyle=lifestyle,
        enable_competition=False,
        include_biological_outputs=True,
    )

    assert disabled.induction_effects.enzyme_folds == enabled.induction_effects.enzyme_folds
    assert disabled.gsh_status.fraction_normal == enabled.gsh_status.fraction_normal
    assert all(
        resolved.inhibition_burden_multiplier == pytest.approx(1.0)
        for resolved in disabled.mechanism_resolved_risks.values()
    )
    assert disabled.mechanism_resolved_risks["PAH"].matrix_gsh_penalty > 1.0


def test_biological_output_flag_does_not_change_principal_numbers():
    exposure = {"benzene": 1.0, "NDMA": 1.0, "PAH": 4.0, "ethanol": 10.0}
    with_outputs = compute_interaction_matrix(
        exposure,
        lifestyle={"chronic_alcohol": True},
        include_biological_outputs=True,
    )
    without_outputs = compute_interaction_matrix(
        exposure,
        lifestyle={"chronic_alcohol": True},
        include_biological_outputs=False,
    )

    assert with_outputs.mechanism_attribution is not None
    assert without_outputs.mechanism_attribution is None
    assert _principal_numbers(with_outputs) == _principal_numbers(without_outputs)
    assert with_outputs.mechanism_resolved_risks
    assert without_outputs.mechanism_resolved_risks == {}
    assert any(
        flux.biological_output is not None
        for enzyme_result in with_outputs.competitive_effects.values()
        for flux in enzyme_result.substrates.values()
    )
    assert all(
        flux.biological_output is None
        for enzyme_result in without_outputs.competitive_effects.values()
        for flux in enzyme_result.substrates.values()
    )
    assert "mechanism_resolved_risks" in _interaction_matrix_to_compat_dict(with_outputs)
    assert "mechanism_resolved_risks" not in _interaction_matrix_to_compat_dict(without_outputs)
    assert "mechanism_resolved_risks" in _build_biological_output_integration(with_outputs)
    assert "mechanism_resolved_risks" not in _build_biological_output_integration(without_outputs)


def test_exact_once_non_neutral_induction_gsh_and_susceptibility_factors():
    result = compute_interaction_matrix(
        {"PAH": 8.0, "HCA": 2.0, "NDMA": 2.0, "ethanol": 12.0, "acrolein": 8.0},
        genotypes={"GSTM1": "null"},
        lifestyle={"smoking": True, "chronic_alcohol": True},
        include_biological_outputs=True,
    )

    pah = result.mechanism_resolved_risks["PAH"]
    assert pah.induction_multiplier > 1.0
    assert pah.susceptibility_modifier > 1.0
    assert pah.matrix_gsh_penalty == pytest.approx(
        pah.gsh_pool_penalty * pah.susceptibility_modifier
    )
    assert pah.final_mechanism_multiplier == pytest.approx(
        pah.induction_multiplier * pah.inhibition_burden_multiplier * pah.matrix_gsh_penalty
    )
    assert result.interaction_adjusted_risks["PAH"] == pytest.approx(
        round(result.individual_risks["PAH"] * pah.final_mechanism_multiplier, 3)
    )

    hca_selected = _selected_biological_output(result, "CYP1A1", "PhIP")
    diagnostic_gsh = hca_selected["gsh_capacity_effect"]["detox_penalty_multiplier"]
    hca = result.mechanism_resolved_risks["HCA"]
    assert hca_selected["gsh_capacity_effect"]["included_in_authoritative_adjusted_risk"] is False
    assert hca_selected["effective_burden"]["includes_diagnostic_gsh_capacity"] is False
    assert hca.final_mechanism_multiplier == pytest.approx(
        round(
            hca.induction_multiplier
            * hca.inhibition_burden_multiplier
            * hca.matrix_gsh_penalty,
            6,
        )
    )
    if diagnostic_gsh != pytest.approx(1.0):
        assert hca.final_mechanism_multiplier != pytest.approx(
            hca.induction_multiplier
            * hca.inhibition_burden_multiplier
            * hca.matrix_gsh_penalty
            * diagnostic_gsh
        )
    else:
        assert hca_selected["effective_burden"]["gsh_detox_penalty_ratio"] == pytest.approx(1.0)


def test_matrix_gsh_status_uses_phase7_redox_model_as_authoritative_pool():
    result = compute_interaction_matrix(
        {"PAH": 8.0, "ethanol": 12.0},
        genotypes={"GSTM1": "null"},
        lifestyle={"smoking": True, "chronic_alcohol": True},
        include_biological_outputs=True,
    )

    assert (
        result.gsh_status.model_version
        == GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY.value
    )
    assert result.gsh_status.redox_capacity_ratio == pytest.approx(
        result.gsh_status.fraction_normal
    )
    assert result.gsh_status.detox_penalty_multiplier >= 1.0

    pah = result.mechanism_resolved_risks["PAH"]
    assert pah.gsh_pool_penalty == pytest.approx(
        round(result.gsh_status.detox_penalty_multiplier, 3)
    )
    assert pah.matrix_gsh_penalty == pytest.approx(
        pah.gsh_pool_penalty * pah.susceptibility_modifier
    )
    assert pah.provenance["gsh"]["matrix_gsh_penalty_applied_once"] is True
    assert pah.provenance["gsh"]["diagnostic_gsh_capacity_included"] is False
    assert pah.provenance["gsh"]["legacy_matrix_gsh_behavior"] is False


def test_matrix_gsh_direct_activation_burden_scales_shared_pool():
    low_activation = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 1.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={"PAH": _inhibition_burden(activation_burden_ratio=2.0)},
    )
    high_activation = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 1.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={"PAH": _inhibition_burden(activation_burden_ratio=5.0)},
    )

    assert high_activation.fraction_normal < low_activation.fraction_normal
    assert (
        high_activation.detox_penalty_multiplier
        > low_activation.detox_penalty_multiplier
    )
    assert high_activation.consumption_rate_umol_h_g == pytest.approx(
        low_activation.consumption_rate_umol_h_g * 2.5
    )
    assert (
        high_activation.individual_contributions["PAH_GSTM1"]["upstream_scaling_source"]
        == "direct_activation_burden_ratio"
    )


def test_matrix_gsh_direct_activation_is_preferred_over_explicit_d_times_k():
    status = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": {"flux_umol_h_g": 1.0, "d_factor": 4.0, "k_factor": 5.0}},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={"PAH": _inhibition_burden(activation_burden_ratio=2.0)},
    )

    contribution = status.individual_contributions["PAH_GSTM1"]
    assert contribution["upstream_activation_scale"] == pytest.approx(2.0)
    assert contribution["gsh_consumption_umol_h_g"] == pytest.approx(2.0)
    assert contribution["upstream_scaling_source"] == "direct_activation_burden_ratio"


def test_matrix_gsh_uses_d_times_k_when_unresolved_activation_is_not_direct():
    status = engine._compute_matrix_gsh_redox_status(
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
    neutral = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 1.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={},
    )

    contribution = status.individual_contributions["PAH_GSTM1"]
    assert contribution["upstream_activation_scale"] == pytest.approx(6.0)
    assert contribution["gsh_consumption_umol_h_g"] == pytest.approx(6.0)
    assert contribution["upstream_scaling_source"] == "d_times_k_approximation"
    assert contribution["upstream_scaling_provenance"]["direct_activation_available"] is False
    assert contribution["upstream_scaling_provenance"]["inhibition_status"] == "mechanism_unresolved"
    assert status.fraction_normal < neutral.fraction_normal
    assert status.detox_penalty_multiplier > neutral.detox_penalty_multiplier
    assert "gsh_upstream_activation_missing_neutral" not in {
        warning["code"] for warning in status.warnings
    }


def test_matrix_gsh_missing_upstream_burden_with_unresolved_record_warns_and_uses_neutral():
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
    assert contribution["upstream_activation_scale"] == pytest.approx(1.0)
    assert contribution["gsh_consumption_umol_h_g"] == pytest.approx(2.0)
    assert contribution["upstream_scaling_source"] == "neutral_missing_upstream_activation"
    assert contribution["upstream_scaling_provenance"]["direct_activation_available"] is False
    assert contribution["upstream_scaling_provenance"]["inhibition_status"] == "mechanism_unresolved"
    assert "gsh_upstream_activation_missing_neutral" in {
        warning["code"] for warning in status.warnings
    }


def test_non_gsh_relevant_pathway_does_not_receive_shared_pool_penalty():
    result = compute_interaction_matrix(
        {"HCA": 4.0, "acetaminophen_umol_h_g": 200.0},
        include_biological_outputs=True,
    )

    assert result.gsh_status.detox_penalty_multiplier > 1.0
    hca = result.mechanism_resolved_risks["HCA"]
    assert hca.gsh_pool_penalty == pytest.approx(1.0)
    assert hca.matrix_gsh_penalty == pytest.approx(1.0)
    assert hca.final_mechanism_multiplier == pytest.approx(
        hca.induction_multiplier * hca.inhibition_burden_multiplier
    )


def test_factor_application_and_derived_fields_use_resolved_risks():
    result = compute_interaction_matrix(
        {"benzene": 1.0, "NDMA": 1.0, "ethanol": 10.0},
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )

    for carcinogen, resolved in result.mechanism_resolved_risks.items():
        assert resolved.final_mechanism_multiplier == pytest.approx(
            resolved.induction_multiplier
            * resolved.inhibition_burden_multiplier
            * resolved.matrix_gsh_penalty
        )
        assert result.interaction_adjusted_risks[carcinogen] == pytest.approx(
            round(result.individual_risks[carcinogen] * resolved.final_mechanism_multiplier, 3)
        )

    benzene_inhibition = result.mechanism_resolved_risks["benzene"].provenance["inhibition"]
    assert benzene_inhibition["modifier_applied_once"] is True

    assert result.total_interaction_risk == pytest.approx(
        round(sum(result.interaction_adjusted_risks.values()), 3)
    )
    assert result.interaction_factor == pytest.approx(
        round(result.total_interaction_risk / result.total_independent_risk, 4)
    )
    for pair, synergy in result.synergy_matrix.items():
        left, right = pair.split("_x_")
        independent_total = result.individual_risks[left] + result.individual_risks[right]
        adjusted_total = (
            result.interaction_adjusted_risks[left] + result.interaction_adjusted_risks[right]
        )
        assert synergy == pytest.approx(round(adjusted_total / independent_total, 3))


def test_public_serializers_keep_existing_fields_and_finite_values():
    result = compute_interaction_matrix({"benzene": 1.0, "ethanol": 10.0})
    payload = _interaction_matrix_to_compat_dict(result)

    assert {
        "individual_risks",
        "interaction_adjusted_risks",
        "synergy_matrix",
        "competitive_effects",
        "total_interaction_risk",
        "interaction_factor",
        "mechanism_attribution",
    }.issubset(payload)
    assert "mechanism_resolved_risks" in payload
    json.dumps(payload, allow_nan=False)
    _finite_walk(payload)

    profile = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        include_tissue_report=False,
    )
    assert "mechanism_resolved_risks" in profile.biological_output_integration
    json.dumps(profile.biological_output_integration, allow_nan=False)
    _finite_walk(profile.biological_output_integration)
