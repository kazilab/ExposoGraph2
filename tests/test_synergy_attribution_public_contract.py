import json
from pathlib import Path

import pytest

from ExposoGraph.interaction_engine import (
    _interaction_matrix_to_compat_dict,
    compute_interaction_matrix,
    decompose_synergy,
)
from ExposoGraph.unified_api import patient_risk_query


EXPOSURE = {"PAH": 4.0, "HCA": 1.0, "acrolein": 4.0, "ethanol": 8.0}
GENOTYPES = {"GSTM1": "NM"}
LIFESTYLE = {"smoking": True}


@pytest.fixture(scope="module")
def module5_result():
    return compute_interaction_matrix(
        EXPOSURE,
        genotypes=GENOTYPES,
        lifestyle=LIFESTYLE,
        include_biological_outputs=True,
    )


@pytest.fixture(scope="module")
def decomposition():
    return decompose_synergy(EXPOSURE, genotypes=GENOTYPES, lifestyle=LIFESTYLE)


def _effect_map(effects):
    return {str(effect["mechanism"]): float(effect["effect"]) for effect in effects}


def test_shapley_terms_sum_to_total_interaction_effect(module5_result):
    attribution = module5_result.mechanism_attribution
    effects = _effect_map(attribution["shapley_main_effects"])

    assert attribution["decomposition_basis"] == "eight_state_shapley"
    assert set(effects) == {"induction", "competition", "gsh"}
    assert sum(effects.values()) + attribution["shapley_residual"] == pytest.approx(
        attribution["total_effect"]
    )
    assert attribution["residuals_are_zero_within_tolerance"] is True


def test_induction_competition_gsh_terms_present(module5_result, decomposition):
    attribution = module5_result.mechanism_attribution

    assert set(attribution["mechanisms"]) == {"induction", "competition", "gsh"}
    assert len(attribution["state_values"]) == 8
    assert {_term["mechanisms"][0] for _term in attribution["singleton_effects"]} == {
        "induction",
        "competition",
        "gsh",
    }

    first_pair = next(iter(decomposition.values()))
    assert set(first_pair.main_effects) == {"induction", "competition", "gsh"}
    assert set(first_pair.state_values) == {
        "none",
        "induction",
        "competition",
        "gsh",
        "induction+competition",
        "induction+gsh",
        "competition+gsh",
        "induction+competition+gsh",
    }


def test_pairwise_terms_present_or_explicitly_absent(decomposition):
    assert decomposition

    for record in decomposition.values():
        assert record.decomposition_basis == "eight_state_shapley"
        assert set(record.pairwise_interactions) == {
            "induction+competition",
            "induction+gsh",
            "competition+gsh",
        }
        assert "pairwise_interactions" in record.compatibility_fields["primary_fields"]


def test_residual_zero_or_documented(module5_result, decomposition):
    attribution = module5_result.mechanism_attribution

    assert attribution["metadata"]["residual_policy"] == "no_unexplained_residual_term"
    assert attribution["shapley_residual"] == pytest.approx(0.0)
    assert attribution["interaction_reconstruction_residual"] == pytest.approx(0.0)

    for record in decomposition.values():
        assert record.residual_policy == "numerical_reconstruction_check_only"
        assert record.reconstruction_residual == pytest.approx(0.0)
        assert record.shapley_residual == pytest.approx(0.0)
        assert record.residuals_are_zero_within_tolerance is True


def test_synergy_heatmap_values_match_interaction_matrix(module5_result, decomposition):
    assert set(decomposition) == set(module5_result.synergy_matrix)

    for pair, record in decomposition.items():
        assert record.composite == pytest.approx(module5_result.synergy_matrix[pair])
        left, right = pair.split("_x_")
        independent_total = (
            module5_result.individual_risks[left] + module5_result.individual_risks[right]
        )
        adjusted_total = (
            module5_result.interaction_adjusted_risks[left]
            + module5_result.interaction_adjusted_risks[right]
        )
        assert module5_result.synergy_matrix[pair] == pytest.approx(
            round(adjusted_total / independent_total, 3)
        )


def test_mechanism_attribution_serializes(module5_result, decomposition):
    payload = _interaction_matrix_to_compat_dict(module5_result)

    json.dumps(payload["mechanism_attribution"], allow_nan=False)
    for record in decomposition.values():
        json.dumps(
            {
                "pair": record.pair,
                "main_effects": record.main_effects,
                "pairwise_interactions": record.pairwise_interactions,
                "three_way_interaction": record.three_way_interaction,
                "residual": record.residual,
                "shapley_decomposition": record.shapley_decomposition,
            },
            allow_nan=False,
        )


def test_figure3_can_consume_mechanism_attribution_output(module5_result, decomposition):
    figure3_script = Path("tools/create_figure3_notebook.py").read_text(encoding="utf-8")
    assert "decompose_synergy" in figure3_script
    assert "result.synergy_matrix" in figure3_script

    rows = []
    for pair, record in decomposition.items():
        rows.append(
            {
                "pair": pair,
                "heatmap_value": module5_result.synergy_matrix[pair],
                "dominant_mechanism": record.dominant_mechanism,
                "induction": record.main_effects["induction"],
                "competition": record.main_effects["competition"],
                "gsh": record.main_effects["gsh"],
                "pairwise_terms": record.pairwise_interactions,
                "three_way_term": record.three_way_interaction,
                "residual": record.reconstruction_residual,
            }
        )

    assert rows
    json.dumps(rows, allow_nan=False)


def test_platform_output_contains_synergy_fields(module5_result):
    payload = _interaction_matrix_to_compat_dict(module5_result)

    assert payload["synergy_matrix"] == module5_result.synergy_matrix
    assert payload["interaction_factor"] == module5_result.interaction_factor
    assert payload["mechanism_attribution"]["decomposition_basis"] == "eight_state_shapley"
    assert payload["module5_model_card"]["synergy_decomposition_basis"] == "eight_state_shapley"

    profile = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Liver",
        include_tissue_report=False,
    )
    integration = profile.biological_output_integration
    assert "synergy_matrix" in integration
    assert "interaction_factor" in integration
    assert integration["synergy_reporting"]["authoritative_attribution"] == "mechanism_attribution"
    assert integration["synergy_reporting"]["pairwise_heatmap"] == "synergy_matrix_descriptive"
    assert integration["module5_model_card"]["synergy_decomposition_basis"] == "eight_state_shapley"
