import math

from ExposoGraph.interaction_schema import (
    ApplicabilityDomain,
    ConcentrationBasis,
    EvidenceGrade,
    EvidenceRecord,
    InhibitionMode,
)
from ExposoGraph.kinetic_resolver import (
    request_from_competitive_interaction,
    resolve_reversible_inhibition,
)
from ExposoGraph.parameter_provider import JSONInteractionParameterProvider
from ExposoGraph.parameter_resolution import (
    InhibitionResolutionStatus,
    ParameterResolutionMethod,
    ReversibleInhibitionResolutionRequest,
)


def _base_request(mode, **overrides):
    data = {
        "mode": mode,
        "enzyme": "CYPX",
        "inhibitor": "inhibitor-x",
        "target_substrate": "substrate-y",
        "km_uM": 10.0,
        "inhibitor_concentration_uM": 2.0,
        "substrate_concentration_uM": 20.0,
        "vmax": 1.0,
        "concentration_basis": ConcentrationBasis.UNBOUND,
        "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
        "evidence": EvidenceRecord(
            source="fixture",
            grade=EvidenceGrade.CURATED,
            confidence="high",
            provenance_ref="fixture#kinetics",
        ),
    }
    data.update(overrides)
    return ReversibleInhibitionResolutionRequest(**data)


def _warning_codes(result):
    return {warning.code for warning in result.warnings}


def _assert_no_nonfinite_numbers(value):
    if isinstance(value, float):
        assert math.isfinite(value)
        return
    if isinstance(value, dict):
        for item in value.values():
            _assert_no_nonfinite_numbers(item)
        return
    if isinstance(value, list):
        for item in value:
            _assert_no_nonfinite_numbers(item)


def test_direct_competitive_ki_resolution_calls_reversible_inhibition_kernel():
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.COMPETITIVE, ki_free_enzyme_uM=5.0)
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.applicability_domain is ApplicabilityDomain.IN_DOMAIN
    assert result.ki_free_enzyme_uM.value == 5.0
    assert result.ki_enzyme_substrate_uM is None
    assert result.substrate_to_km_ratio == 2.0
    assert result.kernel_result.equation_id == "reversible_inhibition.mixed.v1"
    assert result.kernel_result.mode is InhibitionMode.COMPETITIVE


def test_direct_pure_noncompetitive_ki_resolution_uses_explicit_equal_arm_assumption():
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.PURE_NONCOMPETITIVE, ki_free_enzyme_uM=5.0)
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.ki_free_enzyme_uM.value == 5.0
    assert result.ki_enzyme_substrate_uM.value == 5.0
    assert "pure_noncompetitive_equal_ki_applied_to_enzyme_substrate_arm" in result.assumptions
    assert result.kernel_result.alpha_free_enzyme == result.kernel_result.alpha_enzyme_substrate


def test_direct_uncompetitive_ki_es_resolution():
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.UNCOMPETITIVE, ki_enzyme_substrate_uM=5.0)
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.ki_free_enzyme_uM is None
    assert result.ki_enzyme_substrate_uM.value == 5.0
    assert result.kernel_result.mode is InhibitionMode.UNCOMPETITIVE


def test_direct_mixed_resolution_requires_two_constants_and_calls_kernel():
    result = resolve_reversible_inhibition(
        _base_request(
            InhibitionMode.MIXED,
            ki_free_enzyme_uM=5.0,
            ki_enzyme_substrate_uM=8.0,
        )
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.ki_free_enzyme_uM.value == 5.0
    assert result.ki_enzyme_substrate_uM.value == 8.0
    assert result.kernel_result.mode is InhibitionMode.MIXED


def test_competitive_ic50_conversion_uses_assay_substrate_and_km():
    result = resolve_reversible_inhibition(
        _base_request(
            InhibitionMode.COMPETITIVE,
            ic50_uM=30.0,
            assay_substrate_concentration_uM=10.0,
        )
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DERIVED
    assert result.ki_free_enzyme_uM.value == 15.0
    assert result.ki_free_enzyme_uM.resolution_method is ParameterResolutionMethod.IC50_CHENG_PRUSOFF_COMPETITIVE
    assert "competitive_ic50_cheng_prusoff" in result.assumptions
    assert result.kernel_result is not None


def test_uncompetitive_ic50_conversion_uses_assay_substrate_and_km():
    result = resolve_reversible_inhibition(
        _base_request(
            InhibitionMode.UNCOMPETITIVE,
            ic50_uM=30.0,
            assay_substrate_concentration_uM=10.0,
        )
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DERIVED
    assert result.ki_enzyme_substrate_uM.value == 15.0
    assert result.ki_enzyme_substrate_uM.resolution_method is ParameterResolutionMethod.IC50_CHENG_PRUSOFF_UNCOMPETITIVE
    assert "uncompetitive_ic50_cheng_prusoff" in result.assumptions
    assert result.kernel_result is not None


def test_pure_noncompetitive_ic50_equals_ki_only_when_mode_is_explicit():
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.PURE_NONCOMPETITIVE, ic50_uM=30.0)
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DERIVED
    assert result.ki_free_enzyme_uM.value == 30.0
    assert result.ki_enzyme_substrate_uM.value == 30.0
    assert result.ki_free_enzyme_uM.resolution_method is ParameterResolutionMethod.IC50_PURE_NONCOMPETITIVE
    assert "pure_noncompetitive_ic50_equals_ki" in result.assumptions


def test_mixed_single_ic50_is_review_required_not_inferred():
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.MIXED, ic50_uM=30.0)
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.kernel_result is None
    assert "MIXED_INHIBITION_REQUIRES_TWO_CONSTANTS" in _warning_codes(result)


def test_unknown_mode_does_not_default_to_competitive():
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.UNKNOWN, ic50_uM=30.0)
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.mode is InhibitionMode.UNKNOWN
    assert result.kernel_result is None
    assert "UNKNOWN_INHIBITION_MODE" in _warning_codes(result)
    assert "IC50_CONVERSION_REQUIRES_MODE" in _warning_codes(result)


def test_missing_ki_produces_controlled_review_status():
    result = resolve_reversible_inhibition(_base_request(InhibitionMode.COMPETITIVE))

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.kernel_result is None
    assert "KI_MISSING" in _warning_codes(result)


def test_missing_km_blocks_ic50_conversion_when_required():
    result = resolve_reversible_inhibition(
        _base_request(
            InhibitionMode.COMPETITIVE,
            km_uM=None,
            ic50_uM=30.0,
            assay_substrate_concentration_uM=10.0,
        )
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.kernel_result is None
    assert "IC50_CONVERSION_REQUIRES_KM" in _warning_codes(result)
    assert "KM_MISSING" in _warning_codes(result)


def test_missing_assay_substrate_blocks_mode_aware_ic50_conversion():
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.COMPETITIVE, ic50_uM=30.0)
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.kernel_result is None
    assert "IC50_CONVERSION_REQUIRES_ASSAY_SUBSTRATE" in _warning_codes(result)


def test_incompatible_concentration_basis_blocks_quantitative_resolution():
    result = resolve_reversible_inhibition(
        _base_request(
            InhibitionMode.COMPETITIVE,
            ki_free_enzyme_uM=5.0,
            concentration_basis=ConcentrationBasis.NOMINAL,
        )
    )

    assert result.status is InhibitionResolutionStatus.REVIEW_REQUIRED
    assert result.applicability_domain is ApplicabilityDomain.OUTSIDE_DOMAIN
    assert result.kernel_result is None
    assert "CONCENTRATION_BASIS_MISMATCH" in _warning_codes(result)


def test_valid_concentration_basis_permits_quantitative_resolution():
    result = resolve_reversible_inhibition(
        _base_request(
            InhibitionMode.COMPETITIVE,
            ki_free_enzyme_uM=5.0,
            concentration_basis=ConcentrationBasis.TISSUE_EFFECTIVE,
            parameter_concentration_basis=ConcentrationBasis.UNKNOWN,
        )
    )

    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.applicability_domain is ApplicabilityDomain.IN_DOMAIN
    assert result.kernel_result is not None


def test_resolution_preserves_parameter_provenance():
    evidence = EvidenceRecord(
        source="parameter_provenance.json",
        grade=EvidenceGrade.CURATED,
        confidence="curated",
        provenance_ref="parameter_provenance.json#pairs/CYPX/inhibitor-x",
    )
    result = resolve_reversible_inhibition(
        _base_request(InhibitionMode.COMPETITIVE, ki_free_enzyme_uM=5.0, evidence=evidence)
    )

    assert result.evidence.provenance_ref == "parameter_provenance.json#pairs/CYPX/inhibitor-x"
    assert result.ki_free_enzyme_uM.evidence.provenance_ref == evidence.provenance_ref


def test_invalid_inputs_do_not_produce_nan_or_infinity():
    result = resolve_reversible_inhibition(
        _base_request(
            InhibitionMode.COMPETITIVE,
            ic50_uM=-30.0,
            assay_substrate_concentration_uM=10.0,
        )
    )

    assert result.status is InhibitionResolutionStatus.INVALID
    assert result.kernel_result is None
    assert "OUTSIDE_REVERSIBLE_INHIBITION_DOMAIN" in _warning_codes(result)
    _assert_no_nonfinite_numbers(result.to_dict())


def test_legacy_competitive_parameter_record_can_be_adapted_without_live_integration():
    provider = JSONInteractionParameterProvider()
    interaction = next(
        item
        for item in provider.get_competitive_interactions("CYP2E1")
        if item.substrate == "ethanol"
    )

    request = request_from_competitive_interaction(
        interaction,
        inhibitor_concentration_uM=100.0,
        substrate_concentration_uM=75.0,
        target_substrate="benzene",
        concentration_basis=ConcentrationBasis.UNBOUND,
        parameter_concentration_basis=ConcentrationBasis.UNBOUND,
    )
    result = resolve_reversible_inhibition(request)

    assert request.ki_free_enzyme_uM == 13000.0
    assert result.status is InhibitionResolutionStatus.RESOLVED_DIRECT
    assert result.evidence.provenance_ref == "parameter_provenance.json#pairs/CYP2E1/ethanol"
    assert result.metadata["adapter"] == "competitive_interaction_v1"
    assert result.metadata["live_integration"] is False
    assert result.kernel_result.mode is InhibitionMode.COMPETITIVE
