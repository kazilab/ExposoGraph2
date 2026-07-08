import json
import math

from ExposoGraph.interaction_engine import (
    _interaction_matrix_to_compat_dict,
    competitive_inhibition_flux,
    compute_interaction_matrix,
)
from ExposoGraph.interaction_schema import (
    ConcentrationBasis,
    EvidenceGrade,
    EvidenceRecord,
    InhibitionMode,
)
from ExposoGraph.kinetic_resolver import (
    get_ki,
    resolve_reversible_inhibition,
)
from ExposoGraph.parameter_resolution import (
    InhibitionResolutionStatus,
    ParameterResolutionMethod,
    ParameterSourceKind,
    ReversibleInhibitionResolutionRequest,
)
from ExposoGraph.unified_api import patient_risk_query


def _evidence() -> EvidenceRecord:
    return EvidenceRecord(
        source="unit-test",
        grade=EvidenceGrade.CURATED,
        confidence="high",
        provenance_ref="unit://kinetic-public-contract",
        notes="Synthetic kinetic public-contract evidence.",
    )


def _request(
    mode: InhibitionMode,
    **overrides,
) -> ReversibleInhibitionResolutionRequest:
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


def _warning_codes(record) -> set[str]:
    return {warning.code for warning in (record.warnings or [])}


def test_direct_ki_resolution_preserves_curated_provenance():
    resolved = get_ki("CYP2E1", "ethanol", target_substrate="benzene")

    assert resolved.value == 13000.0
    assert resolved.unit == "uM"
    assert resolved.source_kind is ParameterSourceKind.CURATED
    assert resolved.resolution_method is ParameterResolutionMethod.MEASURED_VALUE
    assert resolved.warnings is None
    assert resolved.evidence is not None
    assert resolved.evidence.provenance_ref
    assert resolved.metadata["is_curated_ki"] is True


def test_missing_ki_is_review_required_and_not_quantified():
    result = resolve_reversible_inhibition(
        _request(InhibitionMode.COMPETITIVE, ki_free_enzyme_uM=None, ic50_uM=None)
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.kernel_result is None
    assert "KI_MISSING" in _warning_codes(result)


def test_km_as_ki_proxy_is_low_confidence_warned_and_provenanced():
    resolved = get_ki("CYP2E1", "benzene", target_substrate="ethanol")

    assert resolved.source_kind is ParameterSourceKind.ASSUMED
    assert resolved.resolution_method is ParameterResolutionMethod.ASSUMED_EQUAL_KM
    assert resolved.uncertainty is not None
    assert resolved.uncertainty.confidence == "low"
    assert "ki_missing" in _warning_codes(resolved)
    assert "km_used_as_ki_proxy" in _warning_codes(resolved)
    assert resolved.metadata["proxy_source_field"] == "Km_uM"
    assert resolved.evidence is not None
    assert resolved.evidence.provenance_ref


def test_ic50_without_required_context_is_review_required():
    result = resolve_reversible_inhibition(
        _request(
            InhibitionMode.COMPETITIVE,
            ki_free_enzyme_uM=None,
            ic50_uM=30.0,
            assay_substrate_concentration_uM=None,
        )
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.kernel_result is None
    assert "IC50_CONVERSION_REQUIRES_ASSAY_SUBSTRATE" in _warning_codes(result)


def test_unknown_inhibition_mode_is_review_required():
    result = resolve_reversible_inhibition(
        _request(InhibitionMode.UNKNOWN, ki_free_enzyme_uM=None, ic50_uM=30.0)
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.kernel_result is None
    assert "UNKNOWN_INHIBITION_MODE" in _warning_codes(result)
    assert "IC50_CONVERSION_REQUIRES_MODE" in _warning_codes(result)


def test_multiple_inhibition_contexts_are_not_silently_combined():
    flux = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={
            "benzene": [
                {
                    "mode": InhibitionMode.COMPETITIVE,
                    "ki_free_enzyme_uM": 5.0,
                    "inhibitor_concentration_uM": 2.0,
                    "concentration_basis": ConcentrationBasis.UNBOUND,
                    "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
                },
                {
                    "mode": InhibitionMode.UNCOMPETITIVE,
                    "ki_enzyme_substrate_uM": 6.0,
                    "inhibitor_concentration_uM": 3.0,
                    "concentration_basis": ConcentrationBasis.UNBOUND,
                    "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
                },
            ]
        },
    )

    change = flux.substrates["benzene"]
    assert change.kinetic_modifier is None
    assert change.kinetic_resolution_status == "review_required"
    assert change.modifier_applied_once is False
    assert "MULTIPLE_INHIBITORS_NOT_JOINTLY_RESOLVED" in change.kinetic_warning_codes
    kinetic_effect = change.biological_output["kinetic_effect"]
    assert kinetic_effect["review_required"] is True
    assert "MULTIPLE_INHIBITORS_NOT_JOINTLY_RESOLVED" in kinetic_effect["warnings"]


def test_competitive_mode_resolves_modifier():
    result = resolve_reversible_inhibition(
        _request(InhibitionMode.COMPETITIVE, ki_free_enzyme_uM=5.0)
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.kernel_result is not None
    assert result.kernel_result.mode is InhibitionMode.COMPETITIVE
    assert math.isfinite(result.kernel_result.flux_modifier)
    assert result.kernel_result.flux_modifier < 1.0


def test_pure_noncompetitive_mode_resolves_modifier():
    result = resolve_reversible_inhibition(
        _request(InhibitionMode.PURE_NONCOMPETITIVE, ki_free_enzyme_uM=5.0)
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.kernel_result is not None
    assert result.kernel_result.mode is InhibitionMode.PURE_NONCOMPETITIVE
    assert math.isfinite(result.kernel_result.flux_modifier)
    assert "pure_noncompetitive_equal_ki_applied_to_enzyme_substrate_arm" in result.assumptions


def test_uncompetitive_mode_resolves_modifier():
    result = resolve_reversible_inhibition(
        _request(InhibitionMode.UNCOMPETITIVE, ki_enzyme_substrate_uM=5.0)
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.kernel_result is not None
    assert result.kernel_result.mode is InhibitionMode.UNCOMPETITIVE
    assert math.isfinite(result.kernel_result.flux_modifier)


def test_mixed_mode_resolves_when_both_constants_are_present():
    result = resolve_reversible_inhibition(
        _request(
            InhibitionMode.MIXED,
            ki_free_enzyme_uM=5.0,
            ki_enzyme_substrate_uM=8.0,
        )
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.kernel_result is not None
    assert result.kernel_result.mode is InhibitionMode.MIXED
    assert math.isfinite(result.kernel_result.flux_modifier)


def test_kinetic_resolution_serializes_status_warnings_and_provenance():
    flux = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "ethanol": 2000.0},
    )

    change = flux.substrates["benzene"]
    output = change.biological_output
    kinetic_effect = output["kinetic_effect"]
    assert kinetic_effect["status"] in {"resolved_direct", "resolved_derived"}
    assert "warnings" in kinetic_effect
    assert kinetic_effect["centralized_resolver_used"] is True
    assert kinetic_effect["modifier_applied_once"] is True
    assert kinetic_effect["provenance"]["aggregate_resolution"]["competitors"]
    competitor = kinetic_effect["provenance"]["aggregate_resolution"]["competitors"][0]
    assert competitor["source_kind"]
    assert competitor["resolution_method"]
    assert "metadata" in competitor
    json.dumps(output, allow_nan=False)


def test_kinetic_resolution_is_visible_in_platform_outputs():
    matrix = compute_interaction_matrix({"benzene": 1.0, "ethanol": 1.0})
    compat = _interaction_matrix_to_compat_dict(matrix)

    effects = compat["competitive_effects"]["CYP2E1"]["benzene"]
    kinetic_effect = effects["biological_output"]["kinetic_effect"]
    assert "status" in kinetic_effect
    assert "warnings" in kinetic_effect
    assert "provenance" in kinetic_effect
    assert "ki_resolver_statuses" in compat["module5_model_card"]
    json.dumps(compat, allow_nan=False)

    risk = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Liver",
        include_tissue_report=False,
    )
    integration = risk.biological_output_integration
    assert "ki_resolver_statuses" in integration["module5_model_card"]
    if integration["substrate_outputs"]:
        first_output = next(iter(integration["substrate_outputs"].values()))["kinetic_effect"]
        assert "status" in first_output
        assert "warnings" in first_output
        assert "provenance" in first_output
    json.dumps(integration, allow_nan=False)
