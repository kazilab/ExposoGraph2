"""Smoke tests for the COSMIC SBS mutational-signature validation endpoint."""

import pytest

from ExposoGraph import (
    SignatureAttributionScore,
    SignatureAttributionSummary,
    SignaturePrediction,
    carcinogen_class_to_signatures,
    get_mutational_signature_catalog,
    predict_expected_signatures,
    score_signature_attribution,
    signature_summary_to_dict,
)


def test_catalog_has_expected_signatures_and_class_map():
    catalog = get_mutational_signature_catalog()
    sigs = catalog.get("signatures", {})
    class_map = catalog.get("carcinogen_class_map", {})
    # Core signatures must be present
    for required in ("SBS4", "SBS22", "SBS24", "SBS29", "SBS16", "SBS11"):
        assert required in sigs, f"Missing signature {required}"
    for required in ("PAH", "Aflatoxin", "Ethanol", "AlkylatingAgent"):
        assert required in class_map


def test_carcinogen_class_to_signatures_is_safe_for_unknown():
    out = carcinogen_class_to_signatures("NotARealClass")
    assert out == {"primary": [], "secondary": []}


def test_carcinogen_class_to_signatures_returns_primaries():
    out = carcinogen_class_to_signatures("PAH")
    assert isinstance(out["primary"], list)
    assert isinstance(out["secondary"], list)
    # SBS4 is canonical tobacco/PAH primary
    assert "SBS4" in out["primary"]


def test_predict_expected_signatures_empty_risk_returns_empty():
    out = predict_expected_signatures({})
    assert out == []


def test_predict_expected_signatures_respects_primary_dominance():
    # Pure PAH risk -> SBS4 (primary) should outweigh SBS24 (aflatoxin primary)
    predictions = predict_expected_signatures({"PAH": 5.0})
    assert predictions, "Expected at least one prediction"
    assert all(isinstance(p, SignaturePrediction) for p in predictions)
    top = predictions[0]
    assert top.signature == "SBS4"
    assert top.weight > 0
    # Weights are normalized to sum to 1 by default
    assert abs(sum(p.weight for p in predictions) - 1.0) < 1e-3


def test_predict_secondary_weight_scales_contribution():
    # Weighted lower than primary attribution.
    default = predict_expected_signatures(
        {"PAH": 5.0},
        normalize=False,
        secondary_weight=0.35,
    )
    full_weighted = predict_expected_signatures(
        {"PAH": 5.0},
        normalize=False,
        secondary_weight=1.0,
    )
    # Primary signatures should show equal contribution; secondaries should
    # grow when secondary_weight rises.
    default_map = {p.signature: p.weight for p in default}
    full_map = {p.signature: p.weight for p in full_weighted}
    shared = set(default_map) & set(full_map)
    assert shared
    grew = [s for s in shared if full_map[s] > default_map[s] + 1e-6]
    assert grew, "At least one signature should grow with heavier secondary weighting"


def test_predict_is_robust_to_unknown_classes_and_zero_scores():
    out = predict_expected_signatures(
        {"NotAClass": 99.0, "PAH": 0.0, "Aflatoxin": 4.0}
    )
    # All expected signatures are from Aflatoxin primaries (e.g. SBS24)
    sigs = {p.signature for p in out}
    assert "SBS24" in sigs


def test_score_signature_attribution_perfect_match_has_zero_residual():
    predicted = {"SBS4": 0.8, "SBS24": 0.2}
    summary = score_signature_attribution(predicted, predicted)
    assert isinstance(summary, SignatureAttributionSummary)
    assert summary.total_abs_residual == pytest.approx(0.0, abs=1e-4)
    assert summary.mean_abs_residual == pytest.approx(0.0, abs=1e-4)
    assert summary.spearman_like_agreement == pytest.approx(1.0)


def test_score_signature_attribution_accepts_prediction_list():
    predictions = predict_expected_signatures({"PAH": 3.0})
    summary = score_signature_attribution(
        predictions,
        {p.signature: p.weight for p in predictions},
    )
    assert summary.total_abs_residual == pytest.approx(0.0, abs=1e-4)


def test_score_signature_attribution_reports_per_signature_residuals():
    predicted = {"SBS4": 0.7, "SBS24": 0.3}
    observed = {"SBS4": 0.4, "SBS24": 0.6}
    summary = score_signature_attribution(predicted, observed)
    residuals = {e.signature: e.residual for e in summary.per_signature}
    assert residuals["SBS4"] == pytest.approx(0.3, abs=1e-3)
    assert residuals["SBS24"] == pytest.approx(-0.3, abs=1e-3)
    assert all(isinstance(e, SignatureAttributionScore) for e in summary.per_signature)


def test_summary_serializes_to_dict():
    summary = score_signature_attribution({"SBS4": 1.0}, {"SBS4": 1.0})
    payload = signature_summary_to_dict(summary)
    assert set(payload) == {
        "per_signature",
        "total_abs_residual",
        "mean_abs_residual",
        "spearman_like_agreement",
    }
    assert payload["per_signature"][0]["signature"] == "SBS4"
