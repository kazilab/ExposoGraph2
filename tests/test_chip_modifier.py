"""Smoke tests for the CHIP somatic modifier layer."""

import pytest

from ExposoGraph import (
    AGGREGATE_CHIP_MODIFIER,
    CHIP_DRIVER_GENES,
    VAF_THRESHOLD,
    CHIPEffect,
    CHIPStatus,
    apply_chip_to_risk_map,
    apply_chip_to_risk_score,
    compute_chip_modifier,
    derive_chip_status,
    parse_chip_calls,
    patient_risk_query,
)


def test_parse_chip_calls_with_empty_list_is_negative():
    status = parse_chip_calls([])
    assert isinstance(status, CHIPStatus)
    assert status.present is False
    assert status.mutations == []
    assert status.max_vaf == 0.0


def test_parse_chip_calls_filters_sub_threshold_vaf():
    status = parse_chip_calls([
        {"gene": "DNMT3A", "vaf": 0.005},
        {"gene": "TP53", "vaf": 0.01},
    ])
    assert status.present is False
    assert all(m.variant_allele_fraction < VAF_THRESHOLD for m in status.mutations)


def test_parse_chip_calls_accepts_driver_over_threshold():
    status = parse_chip_calls([{"gene": "DNMT3A", "vaf": 0.08}])
    assert status.present is True
    assert status.max_vaf == pytest.approx(0.08)
    assert status.mutations[0].gene in CHIP_DRIVER_GENES


def test_compute_chip_modifier_returns_aggregate_for_single_driver():
    status = parse_chip_calls([{"gene": "TP53", "vaf": 0.06}])
    effect = compute_chip_modifier(status)
    assert isinstance(effect, CHIPEffect)
    assert effect.chip_positive is True
    assert effect.aggregate_modifier == pytest.approx(AGGREGATE_CHIP_MODIFIER)
    assert "TP53" in effect.driver_genes


def test_compute_chip_modifier_combines_multiple_drivers_multiplicatively():
    status = parse_chip_calls([
        {"gene": "DNMT3A", "vaf": 0.04},
        {"gene": "TP53", "vaf": 0.07},
    ])
    effect = compute_chip_modifier(status)
    assert effect.chip_positive is True
    # Multiple drivers: per-class PAH multiplier should exceed either single effect
    pah_mult = effect.per_class_multipliers.get("PAH", 1.0)
    assert pah_mult > 1.0


def test_apply_chip_to_risk_score_uses_class_multiplier():
    status = parse_chip_calls([{"gene": "TP53", "vaf": 0.08}])
    effect = compute_chip_modifier(status)
    base_score = 25.0
    adjusted = apply_chip_to_risk_score(base_score, effect, carcinogen_class="PAH")
    assert adjusted > base_score


def test_apply_chip_to_risk_map_scales_all_entries():
    status = parse_chip_calls([{"gene": "DNMT3A", "vaf": 0.08}])
    effect = compute_chip_modifier(status)
    base = {"PAH": 100.0, "HCA": 50.0, "benzene": 20.0}
    adjusted = apply_chip_to_risk_map(base, effect)
    assert set(adjusted) == set(base)
    assert all(adjusted[k] >= base[k] for k in base)


def test_derive_chip_status_binary_flag_only():
    status = derive_chip_status(None, is_present=True)
    assert status.present is True
    assert status.mutations == []


def test_chip_absent_returns_identity_effect():
    effect = compute_chip_modifier(None)
    assert effect.chip_positive is False
    assert effect.aggregate_modifier == 1.0
    assert effect.per_class_multipliers == {}


def test_patient_risk_query_accepts_chip_list():
    profile = patient_risk_query(
        genotypes={"GSTM1": "null"},
        tissue="Liver",
        lifestyle={"smoking": True},
        chip_status=[{"gene": "TP53", "vaf": 0.08}],
    )
    assert profile.chip_status is not None
    assert profile.chip_status.present is True
    assert profile.chip_effect is not None
    assert profile.chip_effect.aggregate_modifier == pytest.approx(AGGREGATE_CHIP_MODIFIER)
    if profile.interactions is not None:
        # CHIP-adjusted risk should uplift at least one carcinogen class
        base = profile.interactions.interaction_adjusted_risks
        lifted = [k for k in base if profile.chip_adjusted_risks.get(k, 0.0) > base[k]]
        assert lifted, "CHIP should elevate at least one risk entry"


def test_patient_risk_query_without_chip_is_unchanged():
    profile = patient_risk_query(
        genotypes={"GSTM1": "null"},
        tissue="Liver",
        lifestyle={"smoking": True},
    )
    assert profile.chip_status is None
    assert profile.chip_effect is None
    assert profile.chip_adjusted_risks == {}
