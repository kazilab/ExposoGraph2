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
    genotype_modifier,
    qivive_intrinsic_clearance,
    solve_flux_steady_state,
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
    assert result.susceptibility_score_log2 == round(math.log2(result.net_ratio), 4)
    assert result.steady_state_concentrations_uM["reactive_intermediate_uM"] >= 0
    assert result.steady_state_model["model"] == "one_tissue_perfusion_limited_pbpk_steady_state"
    assert result.steady_state_model["time_to_steady_state_days"] > 0
    assert result.steady_state_concentration_proxy_uM["reactive_intermediate_proxy_uM"] >= 0
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


def test_flux_engine_recognizes_manuscript_genotype_aliases():
    assert genotype_modifier("*1F/*1F", "CYP1A2") == 1.5
    assert genotype_modifier("PM", "CYP1A2") == 0.3
    assert genotype_modifier("slow", "NAT2") == 0.2
    assert genotype_modifier("null", "GSTM1") == 0.05
    assert genotype_modifier("null", "GSTT1") == 0.05


def test_induction_factors_scale_activation_vmax_without_changing_inputs():
    baseline = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "CYP1B1": "NM", "GSTM1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        tissue="Lung",
    )
    induced = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "CYP1B1": "NM", "GSTM1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        tissue="Lung",
        induction_factors={"CYP1A1": 3.0},
    )

    assert induced.total_activation > baseline.total_activation
    assert induced.net_ratio > baseline.net_ratio
    assert induced.induction_factors_used == {"CYP1A1": 3.0}
    cyp1a1 = next(enzyme for enzyme in induced.activation_enzymes if enzyme.enzyme == "CYP1A1")
    assert cyp1a1.induction_modifier == 3.0
    assert "INDUCTION_FACTORS_APPLIED" in induced.warnings


def test_qivive_scaling_exposes_upscaled_flux_context_and_preserves_ratio():
    baseline = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "CYP1B1": "NM", "GSTM1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        tissue="Lung",
    )
    scaled = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "CYP1B1": "NM", "GSTM1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        tissue="Lung",
        qivive=True,
    )

    assert scaled.qivive_applied is True
    assert scaled.qivive_context["mppgl_mg_per_g"] > 0
    assert scaled.qivive_context["organ_weight_g"] > 0
    assert scaled.total_activation > baseline.total_activation
    assert scaled.net_ratio == baseline.net_ratio
    assert "QIVIVE_SCALE_APPLIED" in scaled.warnings


def test_qivive_intrinsic_clearance_helper_uses_mppgl_and_organ_weight():
    assert qivive_intrinsic_clearance(
        10.0,
        2.0,
        microsomal_protein_mg_per_g_tissue=40.0,
        organ_weight_g=1500.0,
    ) == 300000.0


def test_flux_steady_state_solver_uses_explicit_pbpk_context():
    protected = solve_flux_steady_state(
        substrate_conc_uM=1.0,
        activation_flux=2.0,
        detox_flux=6.0,
        tissue="Liver",
    )
    impaired_detox = solve_flux_steady_state(
        substrate_conc_uM=1.0,
        activation_flux=2.0,
        detox_flux=0.2,
        tissue="Liver",
    )

    assert protected.model["central_volume_l"] == 49.0
    assert protected.model["tissue_blood_flow_l_per_day"] > 0
    assert protected.model["extraction_ratio"] > 0
    assert protected.model["time_to_steady_state_days"] > 0
    assert protected.concentrations_uM["central_substrate_uM"] >= 0
    assert impaired_detox.concentrations_uM["reactive_intermediate_uM"] > (
        protected.concentrations_uM["reactive_intermediate_uM"]
    )
