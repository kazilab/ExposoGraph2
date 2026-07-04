"""Smoke tests for the multi-carcinogen interaction engine."""

import json
from dataclasses import asdict
from pathlib import Path

import pytest

from ExposoGraph import (
    InteractionMatrixResult,
    SynergyConfidenceInterval,
    SynergyDecomposition,
    assumed_ki_pairs,
    compute_interaction_matrix,
    decompose_synergy,
    get_interaction_expansion_backlog,
    get_interaction_source_catalog,
    get_parameter_provenance,
    monte_carlo_synergy_ci,
)


MECHANISM_STATE_KEYS = {
    "none",
    "induction",
    "competition",
    "gsh",
    "induction+competition",
    "induction+gsh",
    "competition+gsh",
    "induction+competition+gsh",
}


def test_interaction_matrix_captures_synergy_between_two_classes():
    result = compute_interaction_matrix(
        {"PAH": 1.0, "HCA": 0.5},
        tissue="Liver",
    )

    assert isinstance(result, InteractionMatrixResult)
    assert "PAH" in result.individual_risks
    assert "HCA" in result.individual_risks
    assert result.total_independent_risk > 0
    assert result.total_interaction_risk > 0
    assert result.synergy_matrix is not None


def test_all_mechanisms_disabled_yields_unit_synergy():
    result = compute_interaction_matrix(
        {"PAH": 3.0, "NNK": 4.0, "HCA": 2.0},
        lifestyle={"smoking": True},
        enable_induction=False,
        enable_competition=False,
        enable_gsh_depletion=False,
    )
    assert all(v == pytest.approx(1.0) for v in result.synergy_matrix.values())
    assert result.interaction_factor == pytest.approx(1.0)


def test_decompose_synergy_returns_eight_state_shapley_terms():
    decomp = decompose_synergy(
        {"PAH": 3.0, "HCA": 2.0, "benzene": 6.0},
        lifestyle={"smoking": True},
    )
    assert decomp  # at least one pair
    for pair, dec in decomp.items():
        assert isinstance(dec, SynergyDecomposition)
        assert dec.pair == pair
        assert dec.decomposition_basis == "eight_state_shapley"
        assert set(dec.state_values) == MECHANISM_STATE_KEYS
        assert set(dec.main_effects) == {"induction", "competition", "gsh"}
        assert set(dec.pairwise_interactions) == {
            "induction+competition",
            "induction+gsh",
            "competition+gsh",
        }
        assert isinstance(dec.three_way_interaction, float)
        assert dec.shapley_decomposition["state_count"] == 8
        assert dec.shapley_decomposition["source"] == "compute_interaction_matrix"

        reconstructed = (
            1.0
            + dec.delta_ind
            + dec.delta_comp
            + dec.delta_gsh
            + sum(dec.pairwise_interactions.values())
            + dec.three_way_interaction
            + dec.reconstruction_residual
        )
        assert reconstructed == pytest.approx(dec.composite, abs=1e-6)
        assert dec.residual == pytest.approx(0.0, abs=1e-9)
        assert dec.reconstruction_residual == pytest.approx(0.0, abs=1e-9)
        assert dec.shapley_residual == pytest.approx(0.0, abs=1e-9)
        assert dec.residuals_are_zero_within_tolerance is True
        assert dec.residual_policy == "numerical_reconstruction_check_only"
        assert dec.compatibility_fields["policy"] == "compatibility_only"
        json.dumps(asdict(dec), allow_nan=False)


def test_decompose_synergy_induction_dominates_for_smoker():
    decomp = decompose_synergy(
        {"PAH": 3.0, "HCA": 2.0},
        lifestyle={"smoking": True, "pack_years": 25},
    )
    # Smoking induces CYP1A2/CYP1A1 → induction delta should be the dominant
    # positive contributor for at least one pair involving HCA or PAH.
    any_induction_driven = any(
        dec.main_effects["induction"]
        >= max(abs(dec.main_effects["competition"]), abs(dec.main_effects["gsh"]))
        for dec in decomp.values()
    )
    assert any_induction_driven


def test_decompose_synergy_uses_authoritative_interaction_matrix_values():
    exposure = {"PAH": 3.0, "HCA": 2.0, "benzene": 6.0, "ethanol": 8.0}
    lifestyle = {"smoking": True, "chronic_alcohol": True}
    matrix = compute_interaction_matrix(
        exposure,
        lifestyle=lifestyle,
        include_biological_outputs=False,
    )
    decomp = decompose_synergy(exposure, lifestyle=lifestyle)

    for pair, score in matrix.synergy_matrix.items():
        assert pair in decomp
        assert decomp[pair].composite == pytest.approx(score)
        assert decomp[pair].state_values["induction+competition+gsh"] == pytest.approx(score)


def test_decompose_synergy_biological_output_flag_does_not_change_principal_terms():
    exposure = {"PAH": 3.0, "HCA": 2.0, "benzene": 6.0, "ethanol": 8.0}
    lifestyle = {"smoking": True, "chronic_alcohol": True}
    with_outputs = decompose_synergy(
        exposure,
        lifestyle=lifestyle,
        include_biological_outputs=True,
    )
    without_outputs = decompose_synergy(
        exposure,
        lifestyle=lifestyle,
        include_biological_outputs=False,
    )

    def principal(entry):
        return {
            "composite": entry.composite,
            "main_effects": entry.main_effects,
            "pairwise_interactions": entry.pairwise_interactions,
            "three_way_interaction": entry.three_way_interaction,
            "reconstruction_residual": entry.reconstruction_residual,
            "state_values": entry.state_values,
        }

    assert {pair: principal(dec) for pair, dec in with_outputs.items()} == {
        pair: principal(dec) for pair, dec in without_outputs.items()
    }


def test_figure3_notebook_generator_uses_final_decomposition_fields():
    source = Path("tools/create_figure3_notebook.py").read_text(encoding="utf-8")

    assert "dec.dominant_mechanism" in source
    assert "main_effect_induction" in source
    assert "interaction_induction_competition" in source
    assert "reconstruction_residual" in source
    assert '"residual": "" if dec is None else dec.residual' not in source


def test_s6_decomposition_rows_use_final_decomposition_fields():
    source = Path("tools/generate_interaction_results_s6.py").read_text(encoding="utf-8")

    assert "def _decomposition_row" in source
    assert "residual_policy" in source
    assert '"residual_policy": dec.residual_policy' in source
    assert "compatibility_policy" in source
    for field in {
        "main_effect_induction",
        "main_effect_competition",
        "main_effect_gsh",
        "interaction_induction_competition",
        "interaction_induction_gsh",
        "interaction_competition_gsh",
        "interaction_three_way",
        "reconstruction_residual",
    }:
        assert field in source
    assert '"residual": dec.residual' not in source


def test_monte_carlo_synergy_ci_returns_ordered_bounds():
    intervals = monte_carlo_synergy_ci(
        {"PAH": 3.0, "NNK": 4.0, "HCA": 2.0},
        lifestyle={"smoking": True},
        n_iterations=30,
        seed=13,
    )
    assert intervals
    for ci in intervals.values():
        assert isinstance(ci, SynergyConfidenceInterval)
        lower, upper = ci.composite_ci95
        assert lower <= ci.composite_mean <= upper


def test_parameter_provenance_catalog_lists_expected_pairs():
    prov = get_parameter_provenance()
    assert "_metadata" in prov
    assert "pairs" in prov
    # Every (enzyme, substrate) should have a ki_status entry
    for enzyme, substrates in prov["pairs"].items():
        for sub_name, entry in substrates.items():
            assert entry.get("ki_status") in {"curated", "assumed_equal_km"}


def test_assumed_ki_pairs_enumerate_correctly():
    pairs = assumed_ki_pairs()
    assert pairs
    prov = get_parameter_provenance()
    for enzyme, sub_name in pairs:
        entry = prov["pairs"][enzyme][sub_name]
        assert entry["ki_status"] == "assumed_equal_km"


def test_interaction_source_catalog_lists_brenda_and_rendic_guengerich():
    catalog = get_interaction_source_catalog()
    sources = {entry["source"] for entry in catalog}
    assert "BRENDA enzyme database" in sources
    assert any("Rendic S, Guengerich FP." in source for source in sources)


def test_interaction_expansion_backlog_separates_present_vs_pending_pairs():
    backlog = get_interaction_expansion_backlog()
    assert backlog["reported_gap_pairs"] == 70
    assert backlog["already_parameterized_pairs_from_requested_list"] == 27
    assert backlog["remaining_pairs_to_parameterize"] == 43
    assert "CYP1B1" in backlog["already_parameterized_enzymes_from_requested_list"]
    assert "CYP3A5" in backlog["groups"]["phase_I_cyps"]["pending_enzymes"]
    assert backlog["groups"]["ahr_regulatory"]["model_mode"] == "induction_not_competition"
    assert backlog["groups"]["ahr_regulatory"]["pending_regulators"] == [
        "AHR",
        "ARNT",
    ]


def test_interaction_backlog_scientific_validity_triage_is_structured():
    triage = get_interaction_expansion_backlog()["scientific_validity_triage"]
    assert triage["requested_priority_pairs_count"] == 36
    assert triage["tier_a_requested_pairs"] == 21
    assert triage["tier_b_requested_pairs"] == 15
    assert triage["current_coverage_of_requested_pairs"] == {
        "present_curated": 5,
        "present_assumed_equal_km": 14,
        "missing": 17,
    }
    assert len(triage["green_direct_competition_candidates"]) == 26
    assert len(triage["yellow_provisional_candidates"]) == 7
    assert len(triage["red_exclude_or_rename_pairs"]) == 3

    cyp1b1_e2 = next(
        item
        for item in triage["green_direct_competition_candidates"]
        if item["enzyme"] == "CYP1B1" and item["requested_substrate"] == "E2"
    )
    assert cyp1b1_e2["canonical_substrate"] == "estradiol"
    assert cyp1b1_e2["recommended_action"] == "rename_then_add_as_direct_substrate"

    cyp2e1_acetaldehyde = next(
        item
        for item in triage["red_exclude_or_rename_pairs"]
        if item["enzyme"] == "CYP2E1"
        and item["requested_substrate"] == "acetaldehyde"
    )
    assert cyp2e1_acetaldehyde["recommended_action"] == (
        "exclude_from_competition_matrix"
    )

    gstp1_chromium = next(
        item
        for item in triage["yellow_provisional_candidates"]
        if item["enzyme"] == "GSTP1"
        and item["requested_substrate"] == "chromium(VI)"
    )
    assert gstp1_chromium["recommended_action"] == (
        "avoid_until_direct_substrate_evidence_is_curated"
    )


def test_pulmonary_benzene_activates_cyp2a13_and_cyp2f1():
    """Low-dose inhalation benzene uses pulmonary CYP2A13/CYP2F1 at lung."""
    prov = get_parameter_provenance()
    assert "CYP2A13" in prov["pairs"]
    assert "CYP2F1" in prov["pairs"]
    assert "benzene" in prov["pairs"]["CYP2A13"]
    assert "benzene" in prov["pairs"]["CYP2F1"]


def test_pulmonary_tissue_benzene_interaction_differs_from_liver():
    """The pulmonary pathway should produce a different competition profile.

    NDMA is included so that hepatic CYP2E1 has >= 2 substrates and produces a
    competition multiplier on the liver side, giving the test something to
    compare against the pulmonary CYP2A13/CYP2F1 pathway at lung.
    """
    mixture = {"benzene": 4.0, "NNK": 3.0, "NDMA": 1.5}
    liver_result = compute_interaction_matrix(mixture, tissue="Liver")
    lung_result = compute_interaction_matrix(mixture, tissue="Lung")
    assert liver_result.synergy_matrix
    assert lung_result.synergy_matrix
    diffs = []
    for pair in liver_result.synergy_matrix:
        if pair in lung_result.synergy_matrix:
            diffs.append(
                abs(liver_result.synergy_matrix[pair] - lung_result.synergy_matrix[pair])
            )
    assert any(d > 1e-6 for d in diffs), (
        "Expected pulmonary CYP2A13/CYP2F1 pathway to alter at least one "
        "synergy entry vs the hepatic CYP2E1 pathway."
    )
