"""Smoke tests for the exposure integration engine."""

import pytest

from ExposoGraph import (
    ExposureWeightedRisk,
    LifetimeCancerRisk,
    compute_exposure_weighted_risk,
    compute_lifetime_cancer_risk,
    get_database,
    get_exposure_scenarios,
)


def test_lifetime_cancer_risk_produces_finite_score_above_threshold():
    risk = compute_lifetime_cancer_risk(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "NM"},
        daily_dose_mg_kg=0.001,
    )

    assert isinstance(risk, LifetimeCancerRisk)
    assert risk.carcinogen_class == "PAH"
    assert risk.lecr > 0
    assert risk.slope_factor > 0


def test_lifetime_cancer_risk_normalizes_microgram_slope_factor_units():
    risk = compute_lifetime_cancer_risk(
        "AflatoxinB1",
        {},
        daily_dose_mg_kg=0.000001,
        flux_ratio_override=1.0,
    )

    assert risk.slope_factor == pytest.approx(1000.0)
    assert risk.lecr == pytest.approx(0.001)


def test_exposure_scenarios_are_tier_classified():
    scenarios = {item.scenario_id: item for item in get_exposure_scenarios("PAH")}

    assert scenarios["general_population"].tier == 1
    assert scenarios["smoker"].tier == 2
    assert scenarios["occupational_high"].tier == 3


def test_exposure_database_curation_fixes_are_present():
    db = get_database()
    classes = db["carcinogen_classes"]

    assert classes["VinylChloride"]["class_id"] == 8

    heavy_metal_scenarios = {
        item.scenario_id: item for item in get_exposure_scenarios("HeavyMetals")
    }
    assert heavy_metal_scenarios["smoker"].tissue_conc_uM == pytest.approx(0.004)
    assert heavy_metal_scenarios["occupational_cadmium"].tissue_conc_uM == pytest.approx(0.04)

    nitroso = classes["DietaryNitroso"]
    assert nitroso["bmdl10"]["value"] == 10
    high_processed = nitroso["exposure_scenarios"]["high_processed_meat"]
    assert "total nitrosamines" in high_processed["note"]
    assert "Scientific Reports 2025" in high_processed["source"]

    tce_limits = classes["ChlorinatedSolvents"]["regulatory_limits"]
    assert "100 ppm TWA" in tce_limits["OSHA_PEL_TCE"]
    assert "0.2 ppm" in tce_limits["EPA_TSCA_interim_ECEL_TCE"]


def test_biomarker_measurement_overrides_scenario_exposure_multiplier():
    risk = compute_exposure_weighted_risk(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "NM"},
        "Lung",
        exposure_scenario="general_population",
        biomarker_measurements={"urinary_1_hydroxypyrene": 0.175},
    )

    assert isinstance(risk, ExposureWeightedRisk)
    assert risk.exposure_tier == 3
    assert risk.exposure_multiplier == 3.0
    assert risk.tissue_conc_uM == pytest.approx(0.162)
    assert risk.biomarker_dose_estimate["biomarker"] == "urinary_1_hydroxypyrene"
