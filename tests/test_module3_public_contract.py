"""Public-contract checks for the Module 3 single-carcinogen workflow."""

from __future__ import annotations

import json
import math
from dataclasses import asdict

import pytest


def _reference_module3_result():
    from ExposoGraph import CarcinogenClass, compute_pathway_flux

    return compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )


def test_module3_imports():
    import ExposoGraph
    from ExposoGraph import flux_engine

    assert hasattr(ExposoGraph, "compute_pathway_flux")
    assert hasattr(ExposoGraph, "PathwayFluxResult")
    assert hasattr(flux_engine, "compute_pathway_flux")


def test_module3_public_api_available():
    from ExposoGraph import (
        CarcinogenClass,
        FluxTissueWeightSource,
        PathwayFluxResult,
        compute_pathway_flux,
    )

    assert callable(compute_pathway_flux)
    assert CarcinogenClass.PAH.value == "PAH"
    assert FluxTissueWeightSource.CURATED.value == "curated"
    assert PathwayFluxResult.__name__ == "PathwayFluxResult"


def test_module3_single_carcinogen_default_case():
    from ExposoGraph import FluxTissueWeightSource, PathwayFluxResult, RiskClassification

    result = _reference_module3_result()

    assert isinstance(result, PathwayFluxResult)
    assert result.carcinogen_class == "PAH"
    assert result.tissue == "Lung"
    assert result.substrate_concentration_uM > 0.0
    assert result.total_activation > 0.0
    assert result.total_detox > 0.0
    assert math.isfinite(result.net_ratio)
    assert result.risk_classification in set(RiskClassification)
    assert result.tissue_weight_source == FluxTissueWeightSource.CURATED
    assert result.induction_factors_used == {}


def test_module3_static_km_behavior_preserved_if_applicable():
    from ExposoGraph import kinetic_resolver
    from ExposoGraph.kinetic_resolver import MODULE3_KM_STATIC_FOR_2_0, get_ki

    result = get_ki("CYP2E1", "benzene")

    assert MODULE3_KM_STATIC_FOR_2_0 is True
    assert result.metadata["module3_km_static_for_2_0"] is True
    assert "update_module3_km_from_affinity" not in dir(kinetic_resolver)
    assert "dynamic_module3_km" not in dir(kinetic_resolver)


def test_module3_tissue_expression_vmax_scaling_preserved_if_applicable():
    from ExposoGraph import CarcinogenClass, compute_pathway_flux, get_flux_tissue_weight

    lung_weight = get_flux_tissue_weight("CYP1A1", "Lung")
    liver_weight = get_flux_tissue_weight("CYP1A1", "Liver")
    lung = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )
    liver = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Liver",
    )

    assert lung_weight != liver_weight
    assert lung.total_activation != pytest.approx(liver.total_activation)
    assert lung.net_ratio != pytest.approx(liver.net_ratio)


def test_module3_output_shape_stable():
    result = _reference_module3_result()
    payload = asdict(result)

    expected_keys = {
        "carcinogen_class",
        "tissue",
        "substrate_concentration_uM",
        "genotypes_used",
        "activation_enzymes",
        "detox_enzymes",
        "total_activation",
        "total_detox",
        "net_ratio",
        "susceptibility_score_log2",
        "risk_classification",
        "tissue_weight_source",
        "model_kind",
        "parameter_source",
        "warnings",
        "induction_factors_used",
        "qivive_applied",
        "steady_state_concentrations_uM",
    }

    assert expected_keys <= set(payload)
    assert payload["activation_enzymes"]
    assert payload["detox_enzymes"]
    json.dumps(payload, default=str, sort_keys=True, allow_nan=False)


def test_module3_does_not_require_module5_interaction_context():
    import inspect

    from ExposoGraph import compute_pathway_flux

    signature = inspect.signature(compute_pathway_flux)
    required = [
        name
        for name, parameter in signature.parameters.items()
        if parameter.default is inspect._empty
        and parameter.kind
        in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
    ]

    result = _reference_module3_result()

    assert required == ["carcinogen_class", "genotypes", "tissue"]
    assert result.induction_factors_used == {}
    assert "INDUCTION_FACTORS_APPLIED" not in result.warnings


def test_module3_integration_harness_case():
    from ExposoGraph import PathwayFluxResult

    result = _reference_module3_result()

    assert isinstance(result, PathwayFluxResult)
    assert result.carcinogen_class == "PAH"
    assert result.total_activation > 0.0
    assert result.total_detox > 0.0
    assert math.isfinite(result.net_ratio)
