import ast
import math
from pathlib import Path

import pytest

import ExposoGraph.flux_engine as flux_engine
import ExposoGraph.interaction_engine as interaction_engine
from ExposoGraph.interaction_schema import ConcentrationBasis, InhibitionMode
from ExposoGraph.parameter_resolution import InhibitionResolutionStatus


def _finite_walk(value):
    if isinstance(value, float):
        assert math.isfinite(value)
    elif isinstance(value, dict):
        for item in value.values():
            _finite_walk(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _finite_walk(item)


def test_live_competitive_inhibition_routes_through_centralized_resolver(monkeypatch):
    calls = []
    real_resolver = interaction_engine.resolve_reversible_inhibition

    def spy(request):
        calls.append(request)
        return real_resolver(request)

    monkeypatch.setattr(interaction_engine, "resolve_reversible_inhibition", spy)

    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "ethanol": 2000.0},
    )
    benzene = result.substrates["benzene"]

    assert calls
    assert benzene.centralized_resolver_used is True
    assert benzene.kinetic_resolution_status == InhibitionResolutionStatus.RESOLVED_DIRECT.value
    assert benzene.kinetic_modifier is not None
    assert benzene.modified_flux == benzene.competitive_flux


def test_live_competitive_result_matches_legacy_formula_within_tolerance():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "ethanol": 2000.0},
    )
    benzene = result.substrates["benzene"]

    legacy_single = 10.0 / (75.0 + 10.0)
    legacy_term = 2000.0 / 13000.0
    legacy_competitive = 10.0 / (75.0 * (1.0 + legacy_term) + 10.0)

    assert benzene.single_flux == pytest.approx(legacy_single, abs=1e-6)
    assert benzene.competitive_flux == pytest.approx(legacy_competitive, abs=1e-6)
    assert benzene.inhibition_term == pytest.approx(legacy_term, abs=1e-4)
    assert benzene.discrepancy_classification == "numerical_precision"


def test_flux_engine_applies_kinetic_modifier_exactly_once():
    application = flux_engine.apply_kinetic_modifier_once(10.0, 0.25)

    assert application.baseline_flux == 10.0
    assert application.kinetic_modifier == 0.25
    assert application.modified_flux == 2.5
    assert application.applied_once is True


def test_live_modifier_is_not_applied_twice():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "ethanol": 2000.0},
    )
    benzene = result.substrates["benzene"]

    assert benzene.modifier_applied_once is True
    assert benzene.competitive_flux == pytest.approx(benzene.single_flux * benzene.kinetic_modifier, abs=1e-6)
    assert benzene.competitive_flux != pytest.approx(
        benzene.single_flux * benzene.kinetic_modifier * benzene.kinetic_modifier,
        abs=1e-6,
    )


@pytest.mark.parametrize(
    ("mode", "context"),
    [
        (
            InhibitionMode.PURE_NONCOMPETITIVE,
            {"ki_free_enzyme_uM": 5.0},
        ),
        (
            InhibitionMode.UNCOMPETITIVE,
            {"ki_enzyme_substrate_uM": 5.0},
        ),
        (
            InhibitionMode.MIXED,
            {"ki_free_enzyme_uM": 5.0, "ki_enzyme_substrate_uM": 8.0},
        ),
    ],
)
def test_live_reversible_modes_route_when_sufficient_context_exists(mode, context):
    payload = {
        "mode": mode,
        "inhibitor": "fixture-inhibitor",
        "inhibitor_concentration_uM": 2.0,
        "concentration_basis": ConcentrationBasis.UNBOUND,
        "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
        **context,
    }

    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={"benzene": payload},
    )
    benzene = result.substrates["benzene"]

    assert benzene.centralized_resolver_used is True
    assert benzene.kinetic_resolution_status == InhibitionResolutionStatus.RESOLVED_DIRECT.value
    assert benzene.inhibition_mode == mode.value
    assert benzene.kinetic_modifier is not None
    assert benzene.competitive_flux < benzene.single_flux


def test_absence_of_inhibition_evidence_leaves_legacy_single_substrate_behaviour_unchanged():
    result = interaction_engine.competitive_inhibition_flux("CYP2E1", {"benzene": 10.0})
    benzene = result.substrates["benzene"]

    assert benzene.centralized_resolver_used is False
    assert benzene.kinetic_resolution_status == "mechanism_absent"
    assert benzene.kinetic_modifier == 1.0
    assert benzene.single_flux == benzene.competitive_flux
    assert benzene.modified_flux == benzene.competitive_flux


def test_unknown_mode_does_not_default_to_competitive_in_live_path():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={
            "benzene": {
                "mode": InhibitionMode.UNKNOWN,
                "ic50_uM": 10.0,
                "inhibitor_concentration_uM": 2.0,
                "concentration_basis": ConcentrationBasis.UNBOUND,
                "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
            }
        },
    )
    benzene = result.substrates["benzene"]

    assert benzene.kinetic_modifier is None
    assert benzene.competitive_flux == benzene.single_flux
    assert "UNKNOWN_INHIBITION_MODE" in benzene.kinetic_warning_codes
    assert benzene.inhibition_mode == InhibitionMode.UNKNOWN.value


def test_missing_ki_or_context_does_not_fabricate_live_quantitative_modifier():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={
            "benzene": {
                "mode": InhibitionMode.COMPETITIVE,
                "inhibitor_concentration_uM": 2.0,
                "concentration_basis": ConcentrationBasis.UNBOUND,
                "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
            }
        },
    )
    benzene = result.substrates["benzene"]

    assert benzene.kinetic_modifier is None
    assert benzene.competitive_flux == benzene.single_flux
    assert "KI_MISSING" in benzene.kinetic_warning_codes


def test_unsupported_ic50_only_context_does_not_fabricate_live_modifier():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={
            "benzene": {
                "mode": InhibitionMode.MIXED,
                "ic50_uM": 10.0,
                "inhibitor_concentration_uM": 2.0,
                "concentration_basis": ConcentrationBasis.UNBOUND,
                "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
            }
        },
    )
    benzene = result.substrates["benzene"]

    assert benzene.kinetic_modifier is None
    assert benzene.competitive_flux == benzene.single_flux
    assert "MIXED_INHIBITION_REQUIRES_TWO_CONSTANTS" in benzene.kinetic_warning_codes


def test_nominal_concentration_basis_does_not_silently_produce_live_modifier():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0},
        inhibition_contexts={
            "benzene": {
                "mode": InhibitionMode.COMPETITIVE,
                "ki_free_enzyme_uM": 5.0,
                "inhibitor_concentration_uM": 2.0,
                "concentration_basis": ConcentrationBasis.NOMINAL,
                "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
            }
        },
    )
    benzene = result.substrates["benzene"]

    assert benzene.kinetic_modifier is None
    assert benzene.competitive_flux == benzene.single_flux
    assert "CONCENTRATION_BASIS_MISMATCH" in benzene.kinetic_warning_codes


def test_multiple_explicit_inhibition_contexts_are_deferred_without_quantitative_modifier():
    result = interaction_engine.competitive_inhibition_flux(
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
                    "inhibitor_concentration_uM": 1.5,
                    "concentration_basis": ConcentrationBasis.UNBOUND,
                    "parameter_concentration_basis": ConcentrationBasis.UNBOUND,
                },
            ]
        },
    )
    benzene = result.substrates["benzene"]

    assert benzene.kinetic_modifier is None
    assert benzene.modified_flux is None
    assert benzene.competitive_flux == benzene.single_flux
    assert benzene.kinetic_resolution_status == InhibitionResolutionStatus.REVIEW_REQUIRED.value
    assert "MULTIPLE_INHIBITORS_NOT_JOINTLY_RESOLVED" in benzene.kinetic_warning_codes
    assert benzene.discrepancy_classification == "multiple_inhibitors_not_jointly_resolved"
    _finite_walk(benzene.__dict__)


def test_implicit_mixed_competitor_resolution_does_not_quantify_partial_aggregate():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "ethanol": 2000.0, "unknown_competitor": 100.0},
    )
    benzene = result.substrates["benzene"]

    assert benzene.kinetic_modifier is None
    assert benzene.modified_flux is None
    assert benzene.competitive_flux == benzene.single_flux
    assert benzene.kinetic_resolution_status == InhibitionResolutionStatus.REVIEW_REQUIRED.value
    assert benzene.discrepancy_classification == "incomplete_competitor_resolution"
    assert "INCOMPLETE_COMPETITOR_KI_RESOLUTION" in benzene.kinetic_warning_codes
    assert "KI_MISSING" in benzene.kinetic_warning_codes
    assert benzene.aggregate_resolution["active_competitor_count"] == 2
    assert benzene.aggregate_resolution["resolved_competitor_count"] == 1
    assert benzene.aggregate_resolution["unresolved_competitor_count"] == 1
    assert benzene.aggregate_resolution["all_active_competitors_resolved"] is False
    _finite_walk(benzene.__dict__)


def test_implicit_aggregate_preserves_per_inhibitor_provenance_in_biological_output():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "ethanol": 2000.0, "unknown_competitor": 100.0},
    )
    benzene = result.substrates["benzene"]
    aggregate = benzene.biological_output["kinetic_effect"]["provenance"]["aggregate_resolution"]
    competitors = {item["inhibitor"]: item for item in aggregate["competitors"]}

    assert aggregate["all_active_competitors_resolved"] is False
    assert competitors["ethanol"]["resolved"] is True
    assert competitors["ethanol"]["ki_uM"] == pytest.approx(13000.0)
    assert competitors["ethanol"]["resolution_method"] == "measured_value"
    assert competitors["ethanol"]["source_kind"] == "curated"
    assert competitors["unknown_competitor"]["resolved"] is False
    assert "KI_MISSING" in competitors["unknown_competitor"]["warnings"]
    assert aggregate["aggregate_status"] == InhibitionResolutionStatus.REVIEW_REQUIRED.value
    _finite_walk(aggregate)


def test_live_internal_results_remain_finite():
    result = interaction_engine.competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "ethanol": 2000.0, "NDMA": 1.0},
    )

    for substrate in result.substrates.values():
        _finite_walk(substrate.__dict__)


def test_interaction_and_flux_engines_do_not_duplicate_reversible_equation_source():
    interaction_source = Path(interaction_engine.__file__).read_text(encoding="utf-8")
    flux_source = Path(flux_engine.__file__).read_text(encoding="utf-8")

    assert "Km_power * (1.0 + inhibition_term) + substrate_power" not in interaction_source
    assert "alpha_E" not in interaction_source
    assert "alpha_ES" not in interaction_source
    assert "alpha_E" not in flux_source
    assert "alpha_ES" not in flux_source


def test_public_interaction_compat_payload_keeps_legacy_fields_and_adds_optional_biological_block():
    result = interaction_engine.compute_interaction_matrix(
        {"benzene": 1.0, "ethanol": 1.0},
        enable_induction=False,
        enable_gsh_depletion=False,
    )
    payload = interaction_engine._interaction_matrix_to_compat_dict(result)
    benzene = payload["competitive_effects"]["CYP2E1"]["benzene"]

    assert {
        "single_flux",
        "competitive_flux",
        "flux_change_fraction",
        "inhibition_term",
        "activated_product_flux",
        "Km_uM",
        "concentration_uM",
        "product",
        "product_carcinogenic",
    }.issubset(benzene)
    assert "biological_output" in benzene
    assert benzene["biological_output"]["kinetic_effect"]["status"]


def test_transparency_module_does_not_call_live_inhibition_machinery():
    transparency_path = Path(interaction_engine.__file__).with_name("model_transparency.py")
    tree = ast.parse(transparency_path.read_text(encoding="utf-8"))
    forbidden_modules = {"interaction_engine"}
    forbidden_calls = {
        "compute_interaction_matrix",
        "competitive_inhibition_flux",
        "get_ki",
        "resolve_reversible_inhibition",
    }
    imported_live_symbols = set()
    called_live_symbols = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_live_symbols.update(
                alias.name.rsplit(".", 1)[-1]
                for alias in node.names
                if alias.name.rsplit(".", 1)[-1] in forbidden_modules
            )
        elif isinstance(node, ast.ImportFrom):
            module_name = (node.module or "").rsplit(".", 1)[-1]
            if module_name in forbidden_modules:
                imported_live_symbols.add(module_name)
            imported_live_symbols.update(
                alias.name for alias in node.names if alias.name in forbidden_modules
            )
            imported_live_symbols.update(
                alias.name for alias in node.names if alias.name in forbidden_calls
            )
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id in forbidden_calls:
                called_live_symbols.add(node.func.id)
            elif isinstance(node.func, ast.Attribute) and node.func.attr in forbidden_calls:
                called_live_symbols.add(node.func.attr)

    assert imported_live_symbols == set()
    assert called_live_symbols == set()
