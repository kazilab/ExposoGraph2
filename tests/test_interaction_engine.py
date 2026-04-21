"""Smoke tests for the multi-carcinogen interaction engine."""

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


def test_decompose_synergy_returns_mechanism_deltas():
    decomp = decompose_synergy(
        {"PAH": 3.0, "HCA": 2.0, "benzene": 6.0},
        lifestyle={"smoking": True},
    )
    assert decomp  # at least one pair
    for pair, dec in decomp.items():
        assert isinstance(dec, SynergyDecomposition)
        assert dec.pair == pair
        # composite == 1 + sum(deltas) + residual by construction
        reconstructed = 1.0 + dec.delta_comp + dec.delta_gsh + dec.delta_ind + dec.residual
        assert reconstructed == pytest.approx(dec.composite, abs=1e-3)


def test_decompose_synergy_induction_dominates_for_smoker():
    decomp = decompose_synergy(
        {"PAH": 3.0, "HCA": 2.0},
        lifestyle={"smoking": True, "pack_years": 25},
    )
    # Smoking induces CYP1A2/CYP1A1 → induction delta should be the dominant
    # positive contributor for at least one pair involving HCA or PAH.
    any_induction_driven = any(
        dec.delta_ind >= max(abs(dec.delta_comp), abs(dec.delta_gsh))
        for dec in decomp.values()
    )
    assert any_induction_driven


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
