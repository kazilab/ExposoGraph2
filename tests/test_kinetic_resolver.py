import json

from ExposoGraph.kinetic_resolver import (
    MODULE3_KM_STATIC_FOR_2_0,
    KiResolutionContext,
    KineticParameterResolver,
    get_ki,
)
import ExposoGraph.kinetic_resolver as kinetic_resolver
from ExposoGraph.parameter_provider import JSONInteractionParameterProvider
from ExposoGraph.parameter_resolution import (
    AffinityFallbackStatus,
    ParameterResolutionMethod,
    ParameterSourceKind,
)


def _local_provider_km(enzyme, substrate):
    provider = JSONInteractionParameterProvider()
    interaction = next(
        item
        for item in provider.get_competitive_interactions(enzyme)
        if item.substrate == substrate
    )

    assert interaction.kinetic_parameters.ki_uM is None
    return interaction.kinetic_parameters.km_uM


def test_curated_ki_is_preferred_over_km_proxy():
    result = get_ki("CYP2E1", "ethanol", target_substrate="benzene")

    assert result.value == 13000.0
    assert result.unit == "uM"
    assert result.source_kind is ParameterSourceKind.CURATED
    assert result.resolution_method is ParameterResolutionMethod.MEASURED_VALUE
    assert result.evidence is not None
    assert result.evidence.provenance_ref == "parameter_provenance.json#pairs/CYP2E1/ethanol"
    assert result.metadata["resolution_order"] == "exact_curated_ki"
    assert result.metadata["is_curated_ki"] is True
    assert result.warnings is None


def test_km_proxy_fallback_is_low_confidence_and_warned():
    expected_km_proxy = _local_provider_km("CYP2E1", "benzene")
    result = get_ki("CYP2E1", "benzene", target_substrate="ethanol")

    assert result.value == expected_km_proxy
    assert result.source_kind is ParameterSourceKind.ASSUMED
    assert result.resolution_method is ParameterResolutionMethod.ASSUMED_EQUAL_KM
    assert result.uncertainty.confidence == "low"
    assert result.metadata["proxy_source_field"] == "Km_uM"
    assert result.metadata["is_curated_ki"] is False
    assert {warning.code for warning in result.warnings} >= {"ki_missing", "km_used_as_ki_proxy"}


def test_km_proxy_can_be_disabled():
    result = get_ki(
        "CYP2E1",
        "benzene",
        context=KiResolutionContext(allow_km_proxy=False),
    )

    assert result.value is None
    assert result.resolution_method is ParameterResolutionMethod.UNRESOLVED
    assert any(warning.code == "ki_missing" for warning in result.warnings)
    assert any(warning.code == "no_parameter_resolved" for warning in result.warnings)


def test_unknown_pair_is_unresolved_without_inventing_default():
    result = get_ki("CYP9Z9", "missing")

    assert result.value is None
    assert result.resolution_method is ParameterResolutionMethod.UNRESOLVED
    assert result.metadata["resolution_order"] == "unresolved"
    assert any(warning.code == "ki_missing" for warning in result.warnings)
    assert any(warning.code == "no_parameter_resolved" for warning in result.warnings)


def test_ic50_conversion_is_guarded_without_assay_context():
    expected_km_proxy = _local_provider_km("CYP2E1", "benzene")
    result = get_ki(
        "CYP2E1",
        "benzene",
        context=KiResolutionContext(allow_ic50_conversion=True),
    )

    assert result.value == expected_km_proxy
    assert result.resolution_method is ParameterResolutionMethod.ASSUMED_EQUAL_KM
    assert any(warning.code == "ic50_conversion_unavailable" for warning in result.warnings)
    assert "IC50_uM" not in result.metadata["local_fields_used"]


def test_affinity_fallback_request_remains_unavailable_and_no_helpers_exist():
    result = get_ki(
        "CYP9Z9",
        "missing",
        context=KiResolutionContext(allow_affinity_fallback=True),
    )

    assert result.value is None
    assert result.fallback_status is AffinityFallbackStatus.UNAVAILABLE
    assert any(warning.code == "affinity_fallback_unavailable" for warning in result.warnings)
    assert "convert_affinity_to_ki" not in dir(kinetic_resolver)
    assert "tanimoto_affinity_fallback" not in dir(kinetic_resolver)
    assert "convert_ic50_to_ki" not in dir(kinetic_resolver)
    assert not any("rdkit" in name.lower() for name in dir(kinetic_resolver))
    assert not any("tanimoto" in name.lower() for name in dir(kinetic_resolver))


def test_resolver_uses_json_provider_with_synthetic_local_fixture(tmp_path):
    interaction_data = {
        "competitive_inhibition": {
            "CYPX": {
                "substrates": {
                    "curated": {
                        "Km_uM": 40,
                        "Ki_uM": 2.5,
                        "assumed_ki": False,
                        "provenance_ref": "parameter_provenance.json#pairs/CYPX/curated",
                    },
                    "noncurated": {
                        "Km_uM": 33,
                        "Ki_uM": 9.5,
                        "assumed_ki": False,
                        "provenance_ref": "parameter_provenance.json#pairs/CYPX/noncurated",
                    },
                }
            }
        }
    }
    provenance_data = {
        "pairs": {
            "CYPX": {
                "curated": {
                    "ki_status": "curated",
                    "ki_value_uM": 2.5,
                    "ki_reference": "local curated Ki",
                    "km_confidence": "high",
                },
                "noncurated": {
                    "ki_status": "assumed_equal_km",
                    "ki_value_uM": None,
                    "ki_reference": None,
                    "km_confidence": "low",
                },
            }
        }
    }
    (tmp_path / "interaction_parameters.json").write_text(json.dumps(interaction_data), encoding="utf-8")
    (tmp_path / "parameter_provenance.json").write_text(json.dumps(provenance_data), encoding="utf-8")
    provider = JSONInteractionParameterProvider(data_dir=tmp_path)
    resolver = KineticParameterResolver(provider=provider)

    curated = resolver.get_ki("cypx", "CURATED")
    noncurated = resolver.get_ki("CYPX", "noncurated")

    assert curated.value == 2.5
    assert curated.source_kind is ParameterSourceKind.CURATED
    assert noncurated.value == 9.5
    assert noncurated.source_kind is ParameterSourceKind.LOCAL_METADATA
    assert noncurated.resolution_method is ParameterResolutionMethod.DIRECT_LOOKUP
    assert any(warning.code == "ki_provenance_incomplete" for warning in noncurated.warnings)


def test_module3_km_remains_static_for_2_0():
    result = get_ki("CYP2E1", "benzene")

    assert MODULE3_KM_STATIC_FOR_2_0 is True
    assert result.metadata["module3_km_static_for_2_0"] is True
    assert "update_module3_km_from_affinity" not in dir(kinetic_resolver)
    assert "dynamic_module3_km" not in dir(kinetic_resolver)
