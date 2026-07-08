from __future__ import annotations

import json
from pathlib import Path

import pytest

from ExposoGraph import interaction_engine as engine
from ExposoGraph.effective_burden import (
    GSHBurdenCouplingInput,
    couple_gsh_consumption_to_activation_burden,
)
from ExposoGraph.gsh_redox_capacity import GSHModelVersion
from ExposoGraph.interaction_engine import (
    _interaction_matrix_to_compat_dict,
    compute_interaction_matrix,
)
from ExposoGraph.interaction_equations import (
    adjusted_relative_risk,
    final_mechanism_multiplier,
    gsh_detox_components,
    gsh_upstream_activation_scale,
    pairwise_synergy_factor,
)
from ExposoGraph.interaction_schema import ConcentrationBasis, EvidenceGrade, EvidenceRecord, InhibitionMode
from ExposoGraph.kinetic_resolver import get_ki, resolve_reversible_inhibition
from ExposoGraph.model_transparency import summarize_model_boundaries
from ExposoGraph.parameter_resolution import (
    InhibitionResolutionStatus,
    ParameterResolutionMethod,
    ParameterSourceKind,
    ReversibleInhibitionResolutionRequest,
)


def _warning_codes(records) -> set[str]:
    codes = set()
    for record in records or []:
        if hasattr(record, "code"):
            codes.add(record.code)
        elif isinstance(record, dict):
            codes.add(record.get("code"))
    return {code for code in codes if code}


def _inhibition_burden(
    *,
    activation_burden_ratio: float,
    status: str = "mechanism_resolved",
    review_required: bool = False,
):
    return engine._InhibitionBurdenResolution(
        burden_multiplier=1.0,
        activation_burden_ratio=activation_burden_ratio,
        detox_failure_ratio=1.0,
        endpoint_toxic_flux_ratio=1.0,
        status=status,
        review_required=review_required,
    )


def _evidence() -> EvidenceRecord:
    return EvidenceRecord(
        source="local-test",
        grade=EvidenceGrade.CURATED,
        confidence="high",
        provenance_ref="unit://open-issue-resolution",
    )


def _request(mode: InhibitionMode, **overrides) -> ReversibleInhibitionResolutionRequest:
    values = {
        "enzyme": "CYP2E1",
        "target_substrate": "benzene",
        "inhibitor": "ethanol",
        "mode": mode,
        "km_uM": 10.0,
        "inhibitor_concentration_uM": 2.0,
        "substrate_concentration_uM": 20.0,
        "vmax": 1.0,
        "concentration_basis": ConcentrationBasis.UNBOUND,
        "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
        "evidence": _evidence(),
    }
    values.update(overrides)
    return ReversibleInhibitionResolutionRequest(**values)


def test_issue5_gsh_activation_scaled() -> None:
    low = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 1.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={"PAH": _inhibition_burden(activation_burden_ratio=2.0)},
    )
    high = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 1.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={"PAH": _inhibition_burden(activation_burden_ratio=5.0)},
    )

    assert high.individual_contributions["PAH_GSTM1"]["upstream_scaling_source"] == (
        "direct_activation_burden_ratio"
    )
    assert high.consumption_rate_umol_h_g > low.consumption_rate_umol_h_g
    assert high.fraction_normal < low.fraction_normal


def test_issue5_d_times_k_fallback_warned() -> None:
    coupling = couple_gsh_consumption_to_activation_burden(
        GSHBurdenCouplingInput(
            gsh_relevant=True,
            base_gsh_consumption_load=1.0,
            d_factor=2.0,
            k_factor=3.0,
            tissue="Liver",
        )
    )
    assert coupling.scaling_source == "d_times_k_approximation"
    assert "gsh_d_times_k_fallback_used" in _warning_codes(coupling.warnings)


def test_issue5_neutral_fallback_warned() -> None:
    status = engine._compute_matrix_gsh_redox_status(
        {"PAH_umol_h_g": 2.0},
        genotypes={},
        tissue="Liver",
        inhibition_burdens={
            "PAH": _inhibition_burden(
                activation_burden_ratio=1.0,
                status="mechanism_unresolved",
                review_required=True,
            )
        },
    )

    assert status.individual_contributions["PAH_GSTM1"]["upstream_scaling_source"] == (
        "neutral_missing_upstream_activation"
    )
    assert "gsh_upstream_activation_missing_neutral" in _warning_codes(status.warnings)


def test_issue5_diagnostic_output_not_authoritative() -> None:
    result = compute_interaction_matrix(
        {"PAH": 8.0, "HCA": 2.0, "NDMA": 2.0, "ethanol": 12.0, "acrolein": 8.0},
        genotypes={"GSTM1": "null"},
        lifestyle={"smoking": True, "chronic_alcohol": True},
        include_biological_outputs=True,
    )

    selected = result.competitive_effects["CYP1A1"].substrates["PhIP"].biological_output
    assert selected["selected_authoritative_effect"] is True
    assert selected["gsh_capacity_effect"]["diagnostic_only"] is True
    assert selected["gsh_capacity_effect"]["included_in_authoritative_adjusted_risk"] is False
    assert selected["effective_burden"]["includes_diagnostic_gsh_capacity"] is False


def test_issue7_module5_mechanism_resolved_contract() -> None:
    result = compute_interaction_matrix(
        {"benzene": 1.0, "ethanol": 10.0, "PAH": 4.0, "acrolein": 4.0},
        genotypes={"GSTM1": "NM"},
        include_biological_outputs=True,
    )

    assert result.mechanism_resolved_risks
    assert result.mechanism_attribution["decomposition_basis"] == "eight_state_shapley"
    assert result.gsh_status.model_version == GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY.value
    assert result.synergy_matrix


def test_issue7_transparency_model_card_visible() -> None:
    result = compute_interaction_matrix(
        {"benzene": 1.0, "ethanol": 10.0},
        include_biological_outputs=True,
    )
    payload = _interaction_matrix_to_compat_dict(result)

    assert payload["module5_model_card"]["mechanism_model_version"] == (
        "module5_mechanism_resolved_v2"
    )
    assert "ki_resolver_statuses" in payload["module5_model_card"]
    json.dumps(payload, allow_nan=False)


def test_issue7_equation_helpers_cover_mechanism_components() -> None:
    scale, source, details, warning = gsh_upstream_activation_scale(
        direct_activation_ratio=None,
        direct_status=None,
        direct_review_required=None,
        explicit_dk_scale=6.0,
        explicit_dk_details={"d_factor": 2.0, "k_factor": 3.0},
    )

    assert callable(gsh_detox_components)
    assert callable(final_mechanism_multiplier)
    assert callable(adjusted_relative_risk)
    assert callable(pairwise_synergy_factor)
    assert scale == pytest.approx(6.0)
    assert source == "d_times_k_approximation"
    assert details["internal_d_or_k_computation"] is False
    assert warning == "gsh_d_times_k_fallback_used"


def test_issue9_centralized_ki_resolution_path() -> None:
    resolved = get_ki("CYP2E1", "ethanol", target_substrate="benzene")

    assert resolved.source_kind is ParameterSourceKind.CURATED
    assert resolved.resolution_method is ParameterResolutionMethod.MEASURED_VALUE
    assert resolved.metadata["is_curated_ki"] is True


def test_issue9_km_proxy_warned_low_confidence() -> None:
    resolved = get_ki("CYP2E1", "benzene", target_substrate="ethanol")

    assert resolved.source_kind is ParameterSourceKind.ASSUMED
    assert resolved.resolution_method is ParameterResolutionMethod.ASSUMED_EQUAL_KM
    assert resolved.uncertainty.confidence == "low"
    assert "km_used_as_ki_proxy" in _warning_codes(resolved.warnings)


def test_issue9_ic50_context_review_required() -> None:
    result = resolve_reversible_inhibition(
        _request(
            InhibitionMode.COMPETITIVE,
            ki_free_enzyme_uM=None,
            ic50_uM=30.0,
            assay_substrate_concentration_uM=None,
        )
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert "IC50_CONVERSION_REQUIRES_ASSAY_SUBSTRATE" in _warning_codes(result.warnings)


def test_issue9_unknown_mode_review_required() -> None:
    result = resolve_reversible_inhibition(
        _request(InhibitionMode.UNKNOWN, ki_free_enzyme_uM=None, ic50_uM=30.0)
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert "UNKNOWN_INHIBITION_MODE" in _warning_codes(result.warnings)


def test_issue13_v2_v3_scope_documented() -> None:
    docs = (
        Path("docs/module3_module5_user_paths.rst").read_text(encoding="utf-8")
        + "\n"
        + Path("docs/module5_mechanism_resolved_model.rst").read_text(encoding="utf-8")
    )

    assert "Module 3 remains in v2.0 and is not deprecated" in docs
    assert "deferred to v3 or later" in docs
    assert "does not replace Module 3" in docs


def test_issue13_pending_cases_neutral_visible() -> None:
    result = compute_interaction_matrix(
        {"benzene": 1.0},
        tissue="Lung",
        enable_induction=False,
        enable_gsh_depletion=False,
        include_biological_outputs=True,
    )
    resolved = result.mechanism_resolved_risks["benzene"]

    assert resolved.review_required is True
    assert resolved.inhibition_burden_multiplier == pytest.approx(1.0)
    assert resolved.provenance["inhibition"]["reaction_role_interpretation"]["role"] == "unknown"


def test_issue13_no_global_endpoint_sign_rule() -> None:
    text = " ".join(statement.statement for statement in summarize_model_boundaries())

    assert "product_carcinogenic is evidence and metadata only" in text
    assert "not the final sign rule" in text

