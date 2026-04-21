"""Smoke tests for the epigenomic (promoter-methylation) modifier layer."""

import pytest

from ExposoGraph import (
    AGGREGATE_METHYLATION_MODIFIER,
    EPIGENETIC_GENE_EFFECTS,
    HYPER_BETA_THRESHOLD,
    HYPO_BETA_THRESHOLD,
    EpigeneticEffect,
    MethylationStatus,
    apply_methylation_to_risk_map,
    apply_methylation_to_risk_score,
    compute_epigenetic_modifier,
    derive_methylation_status,
    epigenetic_effect_to_dict,
    get_epigenetic_gene_effects,
    methylation_status_to_dict,
    parse_methylation_calls,
    patient_risk_query,
)


def test_parse_methylation_calls_empty_is_absent():
    status = parse_methylation_calls([])
    assert isinstance(status, MethylationStatus)
    assert status.present is False
    assert status.marks == []


def test_parse_methylation_calls_promotes_curated_hyper():
    status = parse_methylation_calls(
        [{"gene": "MGMT", "beta_value": 0.72}]
    )
    assert status.present is True
    assert status.marks[0].gene == "MGMT"
    assert status.marks[0].resolved_status() == "hyper"


def test_parse_methylation_calls_ignores_normal_beta_values():
    status = parse_methylation_calls(
        [{"gene": "MGMT", "beta_value": 0.35}]
    )
    assert status.present is False
    assert status.marks[0].resolved_status() == "normal"


def test_parse_methylation_calls_respects_explicit_status():
    status = parse_methylation_calls(
        [{"gene": "GSTP1", "status": "hyper"}]
    )
    assert status.present is True
    assert status.marks[0].resolved_status() == "hyper"


def test_parse_methylation_calls_requires_curated_direction():
    # CYP1A1 maps to "hyper" -> protective. Providing "hypo" should not flag.
    status = parse_methylation_calls(
        [{"gene": "CYP1A1", "beta_value": 0.10}]
    )
    # hypo on CYP1A1 isn't a curated trigger (curated direction is hyper)
    assert status.present is False


def test_ahrr_hypomethylation_is_biomarker_only():
    status = parse_methylation_calls(
        [{"gene": "AHRR", "beta_value": 0.05}]
    )
    # AHRR hypo is the canonical smoking biomarker (curated direction=hypo)
    # so detrimental_genes includes it, but the effect has no multipliers.
    assert status.present is True
    effect = compute_epigenetic_modifier(status)
    assert effect.methylation_positive is False
    assert effect.per_class_multipliers == {}
    assert any("biomarker-only" in note for note in effect.notes)


def test_mgmt_hypermethylation_amplifies_nitrosamine_risk():
    status = parse_methylation_calls(
        [{"gene": "MGMT", "beta_value": 0.80}]
    )
    effect = compute_epigenetic_modifier(status)
    assert isinstance(effect, EpigeneticEffect)
    assert effect.methylation_positive is True
    assert effect.per_class_multipliers["Nitrosamine"] > 1.0
    assert effect.per_class_multipliers["AlkylatingAgent"] > 1.0
    assert "MGMT" in effect.affected_genes


def test_cyp1a1_hypermethylation_is_protective_for_pah():
    status = parse_methylation_calls(
        [{"gene": "CYP1A1", "beta_value": 0.75}]
    )
    effect = compute_epigenetic_modifier(status)
    assert effect.methylation_positive is True
    # Bioactivator silencing lowers PAH multiplier below 1.0
    assert effect.per_class_multipliers["PAH"] < 1.0


def test_multiple_genes_combine_multiplicatively_for_pah():
    # GSTP1 hyper: PAH 1.5; OGG1 hyper: PAH 1.3 -> product 1.95
    status = parse_methylation_calls(
        [
            {"gene": "GSTP1", "beta_value": 0.80},
            {"gene": "OGG1", "beta_value": 0.70},
        ]
    )
    effect = compute_epigenetic_modifier(status)
    assert effect.per_class_multipliers["PAH"] == pytest.approx(1.95, abs=1e-3)


def test_aggregate_modifier_applies_when_harmful_present():
    status = parse_methylation_calls(
        [{"gene": "MGMT", "beta_value": 0.80}]
    )
    effect = compute_epigenetic_modifier(status)
    assert effect.aggregate_modifier == pytest.approx(AGGREGATE_METHYLATION_MODIFIER)


def test_aggregate_modifier_is_identity_when_only_protective():
    status = parse_methylation_calls(
        [{"gene": "CYP1A1", "beta_value": 0.75}]
    )
    effect = compute_epigenetic_modifier(status)
    # All multipliers < 1 => no harmful marks => aggregate stays at 1.0
    assert effect.aggregate_modifier == pytest.approx(1.0)


def test_apply_methylation_to_risk_score_scales_class():
    status = parse_methylation_calls(
        [{"gene": "MGMT", "beta_value": 0.80}]
    )
    effect = compute_epigenetic_modifier(status)
    base = 10.0
    adjusted = apply_methylation_to_risk_score(base, effect, "Nitrosamine")
    assert adjusted > base


def test_apply_methylation_to_risk_map_is_identity_when_absent():
    effect = compute_epigenetic_modifier(None)
    risks = {"PAH": 100.0, "HCA": 50.0}
    out = apply_methylation_to_risk_map(risks, effect)
    assert out == {"PAH": 100.0, "HCA": 50.0}


def test_derive_methylation_status_binary_flag():
    status = derive_methylation_status(None, is_present=True)
    assert status.present is True
    assert status.marks == []


def test_status_and_effect_serialize_to_dict():
    status = parse_methylation_calls([{"gene": "MGMT", "beta_value": 0.80}])
    effect = compute_epigenetic_modifier(status)
    s_dict = methylation_status_to_dict(status)
    e_dict = epigenetic_effect_to_dict(effect)
    assert s_dict["present"] is True
    assert "MGMT" in s_dict["detrimental_genes"]
    assert e_dict["methylation_positive"] is True
    assert "Nitrosamine" in e_dict["per_class_multipliers"]


def test_threshold_constants_are_within_expected_bounds():
    assert 0.0 < HYPO_BETA_THRESHOLD < HYPER_BETA_THRESHOLD < 1.0


def test_curated_table_accessor_returns_deep_copy():
    table = get_epigenetic_gene_effects()
    assert "MGMT" in table
    table["MGMT"]["carcinogen_class_multipliers"]["Nitrosamine"] = 99.0
    # Original table is unaffected
    assert EPIGENETIC_GENE_EFFECTS["MGMT"]["carcinogen_class_multipliers"]["Nitrosamine"] != 99.0


def test_patient_risk_query_accepts_methylation_list():
    profile = patient_risk_query(
        genotypes={"GSTM1": "null"},
        tissue="Liver",
        lifestyle={"smoking": True},
        methylation_status=[{"gene": "MGMT", "beta_value": 0.80}],
    )
    assert profile.methylation_status is not None
    assert profile.methylation_status.present is True
    assert profile.methylation_effect is not None
    assert profile.methylation_effect.methylation_positive is True
    if profile.interactions is not None:
        base = profile.interactions.interaction_adjusted_risks
        lifted = [
            k for k in base if profile.methylation_adjusted_risks.get(k, 0.0) > base[k]
        ]
        assert lifted, "Methylation should elevate at least one risk entry"


def test_patient_risk_query_without_methylation_is_unchanged():
    profile = patient_risk_query(
        genotypes={"GSTM1": "null"},
        tissue="Liver",
        lifestyle={"smoking": True},
    )
    assert profile.methylation_status is None
    assert profile.methylation_effect is None
    assert profile.methylation_adjusted_risks == {}
