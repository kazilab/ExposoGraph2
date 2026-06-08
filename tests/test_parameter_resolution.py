import ExposoGraph.parameter_resolution as parameter_resolution
from ExposoGraph.interaction_schema import EvidenceRecord, ParameterUncertainty
from ExposoGraph.parameter_resolution import (
    AffinityFallbackStatus,
    KiResolverRequest,
    ParameterResolutionMethod,
    ParameterResolutionWarning,
    ParameterSourceKind,
    ResolvedParameter,
)


def test_resolved_parameter_can_represent_unresolved_ki_explicitly():
    resolved = ResolvedParameter(
        name="Ki",
        value=None,
        unit="uM",
        source_kind=ParameterSourceKind.LOCAL_METADATA,
        resolution_method=ParameterResolutionMethod.UNRESOLVED,
        evidence=EvidenceRecord(source="parameter_provenance.json"),
        uncertainty=ParameterUncertainty(confidence="low", unit="uM"),
        fallback_status=AffinityFallbackStatus.INACTIVE,
        warnings=[
            ParameterResolutionWarning(
                code="ki_missing",
                message="No Ki value was available in local metadata.",
                source_field="Ki_uM",
            )
        ],
    )

    payload = resolved.to_dict()

    assert payload["value"] is None
    assert payload["source_kind"] == "local_metadata"
    assert payload["resolution_method"] == "unresolved"
    assert payload["fallback_status"] == "inactive"
    assert payload["warnings"][0]["code"] == "ki_missing"


def test_ki_resolver_request_is_a_contract_not_a_resolver_policy():
    request = KiResolverRequest(
        enzyme="CYP2E1",
        inhibitor="benzene",
        target_substrate="toluene",
        allow_ic50_conversion=False,
        allow_affinity_fallback=False,
        allow_km_as_ki_policy=False,
        metadata={"phase": "3"},
    )

    assert request.enzyme == "CYP2E1"
    assert request.inhibitor == "benzene"
    assert request.target_substrate == "toluene"
    assert request.allow_ic50_conversion is False
    assert request.allow_affinity_fallback is False
    assert request.allow_km_as_ki_policy is False
    assert request.metadata["phase"] == "3"


def test_resolution_enums_name_deferred_conversion_paths_without_implementing_them():
    assert ParameterResolutionMethod.IC50_CONVERSION_NOT_IMPLEMENTED.value == "ic50_conversion_not_implemented"
    assert ParameterResolutionMethod.KM_AS_KI_POLICY_NOT_IMPLEMENTED.value == "km_as_ki_policy_not_implemented"
    assert ParameterResolutionMethod.AFFINITY_CONVERSION_NOT_IMPLEMENTED.value == "affinity_conversion_not_implemented"
    assert AffinityFallbackStatus.PENDING.value == "pending"
    assert AffinityFallbackStatus.UNAVAILABLE.value == "unavailable"


def test_phase3_module_exposes_no_ki_conversion_or_fallback_helpers():
    forbidden_helpers = {
        "convert_ic50_to_ki",
        "convert_affinity_to_ki",
        "resolve_ki",
        "get_ki",
        "infer_ki_from_km",
        "tanimoto_affinity_fallback",
    }

    assert forbidden_helpers.isdisjoint(set(dir(parameter_resolution)))
