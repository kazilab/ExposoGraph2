"""Smoke tests for the exposure integration engine."""

from ExposoGraph import (
    LifetimeCancerRisk,
    compute_lifetime_cancer_risk,
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
