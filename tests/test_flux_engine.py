"""Smoke tests for the quantitative flux engine."""

import math
from dataclasses import asdict

import pytest

from ExposoGraph import (
    CarcinogenClass,
    FluxTissueWeightSource,
    PathwayFluxResult,
    RiskClassification,
    compute_full_profile,
    compute_pathway_flux,
    flux_engine,
)


def test_compute_pathway_flux_returns_finite_ratio_for_reference_pah_genotype():
    result = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )

    assert isinstance(result, PathwayFluxResult)
    assert result.total_activation > 0
    assert result.total_detox > 0
    assert math.isfinite(result.net_ratio)
    assert result.risk_classification in set(RiskClassification)
    assert result.tissue_weight_source == FluxTissueWeightSource.CURATED.value


def test_measured_flux_result_exposes_measured_parameter_source():
    result = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTP1": "NM"},
        tissue="Lung",
    )

    payload = asdict(result)
    assert payload["model_kind"] == "measured_kinetics"
    assert payload["parameter_source"] == "kinetic_parameters.json"
    assert all(
        enzyme["parameter_source"] == "kinetic_parameters.json"
        for enzyme in payload["activation_enzymes"] + payload["detox_enzymes"]
    )
    assert all(
        enzyme["provenance_ref"] == ""
        for enzyme in payload["activation_enzymes"] + payload["detox_enzymes"]
    )


def test_proxy_flux_result_exposes_provenance_metadata():
    result = compute_pathway_flux(
        CarcinogenClass.NDEA,
        {"CYP2A13": "NM", "CYP2E1": "NM", "MGMT": "NM"},
        tissue="Liver",
    )

    payload = asdict(result)
    assert payload["model_kind"] == "semi_quantitative_proxy"
    assert payload["parameter_source"] == "proxy_flux_parameters.json"

    cyp2a13 = next(
        enzyme for enzyme in payload["activation_enzymes"] if enzyme["enzyme"] == "CYP2A13"
    )
    assert cyp2a13["parameter_source"] == "proxy_flux_parameters.json"
    assert cyp2a13["provenance_ref"] == "classes.NDEA.activation_terms.CYP2A13"
    assert cyp2a13["provenance_sources"]
    assert cyp2a13["parameter_basis"]


def test_compute_full_profile_covers_all_supported_classes():
    profile = compute_full_profile(
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM", "NAT2": "NM"},
        tissue="Liver",
    )

    assert profile.total_classes_modeled == len(flux_engine._DISPATCH)
    assert "PAH" in profile.per_class_results
    assert "AromaticAmines" in profile.per_class_results
    assert "EstrogenMetabolites" in profile.per_class_results
    assert "NDEA" in profile.per_class_results
    assert "VinylChloride" in profile.per_class_results
    assert "ChlorinatedSolvent" in profile.per_class_results
    assert "UV_Radiation" in profile.per_class_results
    assert "Dioxin" in profile.per_class_results
    assert "HeavyMetal" in profile.per_class_results


@pytest.mark.parametrize(
    ("carcinogen_class", "tissue", "genotypes"),
    [
        (
            CarcinogenClass.CHLORINATED_SOLVENT,
            "Kidney",
            {"CYP2E1": "NM", "GSTT1": "NM"},
        ),
        (
            CarcinogenClass.DIOXIN,
            "Lung",
            {"CYP1A1": "NM", "CYP1B1": "NM"},
        ),
        (
            CarcinogenClass.HEAVY_METAL,
            "Liver",
            {"AS3MT": "NM"},
        ),
        (
            CarcinogenClass.AROMATIC_AMINES,
            "Bladder",
            {"CYP1A2": "NM", "NAT1": "NM", "NAT2": "NM", "GSTM1": "NM"},
        ),
        (
            CarcinogenClass.ESTROGEN_METABOLITES,
            "Breast",
            {"CYP1B1": "NM", "COMT": "NM", "SULT1E1": "NM", "UGT2B7": "NM"},
        ),
        (
            CarcinogenClass.NDEA,
            "Liver",
            {"CYP2A13": "NM", "CYP2E1": "NM", "MGMT": "NM"},
        ),
        (
            CarcinogenClass.VINYL_CHLORIDE,
            "Liver",
            {"CYP2E1": "NM", "GSTT1": "NM", "EPHX1": "NM", "ALDH2": "NM"},
        ),
        (
            CarcinogenClass.UV_RADIATION,
            "Skin",
            {"XPC": "NM", "ERCC2": "NM", "OGG1": "NM", "POLH": "NM"},
        ),
    ],
)
def test_compute_pathway_flux_supports_packaged_wave2_classes(
    carcinogen_class,
    tissue,
    genotypes,
):
    result = compute_pathway_flux(carcinogen_class, genotypes, tissue=tissue)

    assert isinstance(result, PathwayFluxResult)
    assert result.risk_classification != RiskClassification.INSUFFICIENT_DATA
    assert result.total_activation > 0
    assert result.total_detox > 0
    assert math.isfinite(result.net_ratio)
    assert "ESTIMATED_PARAMS" in result.warnings
