from __future__ import annotations

import json
from dataclasses import asdict

import pytest

from ExposoGraph import get_default_v2_provider
from ExposoGraph.flux_equations import (
    activation_detox_ratio,
    hill_equation,
    michaelis_menten,
    scaled_vmax,
    susceptibility_score_log2,
)
from ExposoGraph.flux_engine import compute_pathway_flux
from ExposoGraph.interaction_engine import (
    _interaction_matrix_to_compat_dict,
    compute_interaction_matrix,
    decompose_synergy,
)
from ExposoGraph.interaction_equations import (
    adjusted_relative_risk,
    final_mechanism_multiplier,
    gsh_detox_components,
    gsh_upstream_activation_scale,
    pairwise_synergy_factor,
)


def _module3_payload(result) -> dict:
    payload = asdict(result)
    payload["risk_classification"] = result.risk_classification.value
    payload["tissue_weight_source"] = result.tissue_weight_source.value
    return payload


def _module5_principal_payload(result) -> dict:
    payload = _interaction_matrix_to_compat_dict(result)
    return {
        "individual_risks": payload["individual_risks"],
        "interaction_adjusted_risks": payload["interaction_adjusted_risks"],
        "synergy_matrix": payload["synergy_matrix"],
        "total_independent_risk": payload["total_independent_risk"],
        "total_interaction_risk": payload["total_interaction_risk"],
        "interaction_factor": payload["interaction_factor"],
        "module5_model_card": payload["module5_model_card"],
    }


def test_engine_equation_helpers_import() -> None:
    assert michaelis_menten(1.0, 2.0, 3.0) == pytest.approx(0.5)
    assert hill_equation(1.0, 2.0, 1.0, 2.0) == pytest.approx(1.0)
    assert scaled_vmax(2.0, 0.5, 0.25) == pytest.approx(0.25)
    assert activation_detox_ratio(4.0, 2.0) == pytest.approx(2.0)
    assert susceptibility_score_log2(2.0) == pytest.approx(1.0)
    assert final_mechanism_multiplier(2.0, 0.5, 1.5) == pytest.approx(1.5)
    assert adjusted_relative_risk(10.0, 1.5) == pytest.approx(15.0)
    assert pairwise_synergy_factor(10.0, 5.0, 12.0, 6.0) == pytest.approx(1.2)
    assert gsh_detox_components("PAH", 0.25, {"GSTM1": "active"}, {"PAH": "GSTM1"})[1] > 1.0
    assert gsh_upstream_activation_scale(
        direct_activation_ratio=2.0,
        direct_status="mechanism_resolved",
        direct_review_required=False,
        explicit_dk_scale=None,
        explicit_dk_details=None,
    )[1] == "direct_activation_burden_ratio"


def test_refactor_preserves_module3_outputs() -> None:
    direct = compute_pathway_flux(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )
    via_provider = get_default_v2_provider().compute_module3_pathway_flux(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )

    assert _module3_payload(via_provider) == _module3_payload(direct)
    json.dumps(_module3_payload(direct), sort_keys=True, allow_nan=False)


def test_refactor_preserves_module5_outputs() -> None:
    exposure = {"benzene": 1.0, "ethanol": 10.0, "PAH": 4.0, "acrolein": 4.0}
    kwargs = {
        "genotypes": {"GSTM1": "NM"},
        "tissue": "Liver",
        "include_biological_outputs": True,
    }
    direct = compute_interaction_matrix(exposure, **kwargs)
    via_provider = get_default_v2_provider().compute_module5_interaction_matrix(exposure, **kwargs)

    assert _module5_principal_payload(via_provider) == _module5_principal_payload(direct)
    json.dumps(_module5_principal_payload(direct), sort_keys=True, allow_nan=False)


def test_refactor_preserves_figure3_input_shape() -> None:
    exposure = {"PAH": 4.0, "HCA": 1.0, "acrolein": 4.0, "ethanol": 8.0}
    result = compute_interaction_matrix(
        exposure,
        genotypes={"GSTM1": "NM"},
        lifestyle={"smoking": True},
        include_biological_outputs=True,
    )
    decomposition = decompose_synergy(
        exposure,
        genotypes={"GSTM1": "NM"},
        lifestyle={"smoking": True},
    )

    rows = [
        {
            "pair": pair,
            "heatmap_value": result.synergy_matrix[pair],
            "dominant_mechanism": record.dominant_mechanism,
            "induction": record.main_effects["induction"],
            "competition": record.main_effects["competition"],
            "gsh": record.main_effects["gsh"],
            "pairwise_terms": record.pairwise_interactions,
            "three_way_term": record.three_way_interaction,
            "residual": record.reconstruction_residual,
        }
        for pair, record in decomposition.items()
    ]

    assert rows
    assert {row["pair"] for row in rows} == set(result.synergy_matrix)
    json.dumps(rows, sort_keys=True, allow_nan=False)

