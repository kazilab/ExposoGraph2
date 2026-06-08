import json
from math import isclose

import pytest

from ExposoGraph.effective_burden import compute_effective_carcinogenic_burden
from ExposoGraph.endpoint_toxic_flux import interpret_competitive_endpoint_flux
from ExposoGraph.gsh_redox_capacity import (
    compute_gsh_redox_capacity,
    compute_quasi_steady_gsh_fraction,
)
from ExposoGraph.interaction_schema import (
    AssumptionWarning,
    CompetitiveInteraction,
    EvidenceGrade,
    EvidenceRecord,
    KineticParameterSet,
    MetabolicReaction,
    ParameterUncertainty,
    ReactionRole,
    RiskDirectionIfFluxDecreases,
    SMEReviewStatus,
)
from ExposoGraph.kinetic_resolver import KiResolutionContext, KineticParameterResolver
from ExposoGraph.mechanism_attribution import (
    DEFAULT_MECHANISMS,
    MechanismName,
    MechanismState,
    compute_mechanism_attribution,
    evaluate_mechanism_states,
    validate_eight_state_values,
)
from ExposoGraph.model_transparency import (
    AssumptionCategory,
    build_transparency_report,
)
from ExposoGraph.parameter_resolution import (
    AffinityFallbackStatus,
    ParameterResolutionMethod,
    ParameterSourceKind,
)
from ExposoGraph.reaction_role_semantics import ReactionRoleAnnotation


def _warning_codes(warnings):
    return {warning.code for warning in warnings or []}


def _normalize(value):
    return "".join(character.lower() for character in str(value) if character.isalnum())


class _LocalInteractionProvider:
    def __init__(self, interactions):
        self._interactions = tuple(interactions)

    def get_competitive_interactions(self, enzyme):
        enzyme_key = _normalize(enzyme)
        return [
            interaction
            for interaction in self._interactions
            if _normalize(interaction.enzyme) == enzyme_key
        ]


@pytest.mark.parametrize(
    ("synthesis", "load"),
    [
        (0.0, 0.0),
        (0.0, 2.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (2.5, 7.5),
        (-1.0, 3.0),
        (4.0, -2.0),
    ],
)
def test_gsh_fraction_is_bounded_for_table_driven_synthesis_load_cases(synthesis, load):
    fraction = compute_quasi_steady_gsh_fraction(synthesis, load)

    assert 0.0 <= fraction <= 1.0

    result = compute_gsh_redox_capacity(
        tissue="liver",
        synthesis_capacity=synthesis,
        consumption_load=load,
        baseline_capacity=1.0,
    )
    assert 0.0 <= result.gsh_fraction <= 1.0
    assert 0.0 <= result.redox_capacity_ratio <= 1.0
    assert result.detox_penalty_multiplier >= 1.0


def test_gsh_detox_penalty_increases_as_redox_capacity_decreases():
    loads = [0.0, 0.5, 1.0, 3.0, 9.0]
    results = [
        compute_gsh_redox_capacity(
            tissue="default",
            synthesis_capacity=1.0,
            consumption_load=load,
            baseline_capacity=1.0,
        )
        for load in loads
    ]

    fractions = [result.gsh_fraction for result in results]
    penalties = [result.detox_penalty_multiplier for result in results]

    assert fractions == sorted(fractions, reverse=True)
    assert penalties == sorted(penalties)
    assert penalties[-1] > penalties[0]


def test_effective_burden_neutral_nonnegative_and_monotonic_in_positive_factor():
    neutral = compute_effective_carcinogenic_burden()
    assert neutral.effective_carcinogenic_burden_ratio == 1.0
    assert neutral.activation_burden_ratio == 1.0
    assert neutral.detox_failure_ratio == 1.0

    negative = compute_effective_carcinogenic_burden(
        activation_burden_ratio=-3.0,
        detox_failure_ratio=2.0,
        susceptibility_modifier=1.5,
    )
    assert negative.effective_carcinogenic_burden_ratio == 0.0
    assert "activation_burden_ratio_negative_clamped" in _warning_codes(negative.warnings)

    activation_values = [0.25, 1.0, 2.0, 4.0]
    burdens = [
        compute_effective_carcinogenic_burden(
            activation_burden_ratio=activation,
            detox_failure_ratio=1.25,
            susceptibility_modifier=1.1,
            gsh_relevant=True,
            gsh_detox_penalty_ratio=1.5,
        ).effective_carcinogenic_burden_ratio
        for activation in activation_values
    ]

    assert burdens == sorted(burdens)
    assert all(burden >= 0.0 for burden in burdens)
    assert burdens[-1] > burdens[0]


@pytest.mark.parametrize(
    ("annotation", "expected_warning"),
    [
        (
            ReactionRoleAnnotation(
                enzyme="CYPX",
                substrate="unknown_substrate",
                reaction_role=ReactionRole.UNKNOWN,
                risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.UNKNOWN,
                metadata={"product_carcinogenic": True},
            ),
            "endpoint_role_unknown",
        ),
        (
            ReactionRoleAnnotation(
                enzyme="CYPX",
                substrate="probe",
                reaction_role=ReactionRole.PROBE_ONLY,
                risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.NEUTRAL,
            ),
            "endpoint_probe_only_neutral",
        ),
        (
            ReactionRoleAnnotation(
                enzyme="CYPX",
                substrate="pending",
                reaction_role=ReactionRole.BIOACTIVATION,
                risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.DECREASE,
                review_status=SMEReviewStatus.PENDING_TEAM_AGREEMENT,
            ),
            "endpoint_role_pending",
        ),
        (
            ReactionRoleAnnotation(
                enzyme="CYPX",
                substrate="inactive",
                reaction_role=ReactionRole.DETOXIFICATION,
                risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.INCREASE,
                active=False,
            ),
            "endpoint_role_inactive",
        ),
    ],
)
def test_endpoint_flux_neutral_roles_do_not_infer_burden_from_metadata(annotation, expected_warning):
    result = interpret_competitive_endpoint_flux(0.25, annotation)

    assert result.endpoint_toxic_flux_ratio == 1.0
    assert result.activation_burden_ratio == 1.0
    assert result.detox_failure_ratio == 1.0
    assert expected_warning in _warning_codes(result.warnings)
    assert result.metadata["reaction_role_metadata"] == dict(annotation.metadata or {})

    for invalid_ratio in (0.0, -1.0):
        with pytest.raises(ValueError, match="flux_ratio must be greater than zero"):
            interpret_competitive_endpoint_flux(invalid_ratio, annotation)


@pytest.mark.parametrize(
    "evaluator",
    [
        lambda state: (
            1.0
            + (2.0 if state.induction else 0.0)
            + (3.0 if state.competition else 0.0)
            + (4.0 if state.gsh else 0.0)
        ),
        lambda state: (
            1.0
            + (2.0 if state.induction else 0.0)
            + (3.0 if state.competition else 0.0)
            + (4.0 if state.gsh else 0.0)
            + (0.5 if state.induction and state.competition else 0.0)
            + (0.25 if state.induction and state.gsh else 0.0)
            + (0.75 if state.competition and state.gsh else 0.0)
            + (1.25 if state.induction and state.competition and state.gsh else 0.0)
        ),
    ],
)
def test_mechanism_attribution_residuals_and_reconstruction_invariants(evaluator):
    state_values = evaluate_mechanism_states(evaluator)
    result = compute_mechanism_attribution(state_values, tolerance=1e-9)

    shapley_total = sum(effect.effect for effect in result.shapley_main_effects)
    explicit_total = (
        sum(term.effect for term in result.singleton_effects)
        + sum(term.effect for term in result.pairwise_interactions)
        + result.three_way_interaction.effect
    )

    assert len(result.state_values) == 8
    assert [state_value.state.key for state_value in result.state_values] == [
        "none",
        "induction",
        "competition",
        "gsh",
        "induction+competition",
        "induction+gsh",
        "competition+gsh",
        "induction+competition+gsh",
    ]
    assert isclose(shapley_total, result.total_effect, abs_tol=1e-9)
    assert isclose(explicit_total, result.total_effect, abs_tol=1e-9)
    assert result.shapley_residual == 0.0
    assert result.interaction_reconstruction_residual == 0.0
    assert result.residuals_are_zero_within_tolerance is True
    assert result.metadata["not_biological_causality"] is True


def test_mechanism_state_validation_rejects_missing_and_invalid_states():
    with pytest.raises(ValueError, match="Missing Phase 9 mechanism state values"):
        validate_eight_state_values({"none": 1.0})

    complete_values = {
        state.key: float(index)
        for index, state in enumerate(
            [MechanismState.from_active(active) for active in [(), *[(m,) for m in DEFAULT_MECHANISMS]]]
        )
    }
    complete_values["unknown"] = 9.0
    with pytest.raises(ValueError, match="Unknown Phase 9 mechanism"):
        validate_eight_state_values(complete_values)

    with pytest.raises(ValueError, match="tolerance must be a non-negative finite value"):
        compute_mechanism_attribution(
            {
                state.key: float(index)
                for index, state in enumerate(
                    [
                        MechanismState(),
                        MechanismState.from_active((MechanismName.INDUCTION,)),
                        MechanismState.from_active((MechanismName.COMPETITION,)),
                        MechanismState.from_active((MechanismName.GSH,)),
                        MechanismState.from_active((MechanismName.INDUCTION, MechanismName.COMPETITION)),
                        MechanismState.from_active((MechanismName.INDUCTION, MechanismName.GSH)),
                        MechanismState.from_active((MechanismName.COMPETITION, MechanismName.GSH)),
                        MechanismState.from_active(DEFAULT_MECHANISMS),
                    ]
                )
            },
            tolerance=-1.0,
        )


def test_kinetic_resolver_guards_ic50_affinity_and_warns_on_km_proxy():
    km_proxy_interaction = CompetitiveInteraction(
        enzyme="CYP1A1",
        substrate="BaP",
        kinetic_parameters=KineticParameterSet(
            km_uM=7.5,
            ki_uM=None,
            uncertainty=ParameterUncertainty(confidence="local"),
        ),
        evidence=EvidenceRecord(
            source="local fixture",
            grade=EvidenceGrade.PLACEHOLDER,
            confidence="local_metadata",
        ),
    )
    unresolved_interaction = CompetitiveInteraction(
        enzyme="CYP1A2",
        substrate="NoKi",
        kinetic_parameters=KineticParameterSet(km_uM=None, ki_uM=None),
    )
    resolver = KineticParameterResolver(
        provider=_LocalInteractionProvider([km_proxy_interaction, unresolved_interaction])
    )

    proxy = resolver.get_ki(
        "CYP1A1",
        "BaP",
        context=KiResolutionContext(allow_ic50_conversion=True, allow_km_proxy=True),
    )
    assert proxy.value == 7.5
    assert proxy.source_kind is ParameterSourceKind.ASSUMED
    assert proxy.resolution_method is ParameterResolutionMethod.ASSUMED_EQUAL_KM
    assert proxy.uncertainty.confidence == "low"
    assert proxy.metadata["is_curated_ki"] is False
    assert proxy.metadata["proxy_source_field"] == "Km_uM"
    assert {"ki_missing", "ic50_conversion_unavailable", "km_used_as_ki_proxy"} <= _warning_codes(
        proxy.warnings
    )

    unresolved = resolver.get_ki(
        "CYP1A2",
        "NoKi",
        context=KiResolutionContext(
            allow_ic50_conversion=True,
            allow_km_proxy=False,
            allow_affinity_fallback=True,
        ),
    )
    assert unresolved.value is None
    assert unresolved.resolution_method is ParameterResolutionMethod.UNRESOLVED
    assert unresolved.fallback_status is AffinityFallbackStatus.UNAVAILABLE
    assert {
        "ki_missing",
        "ic50_conversion_unavailable",
        "affinity_fallback_unavailable",
        "no_parameter_resolved",
    } <= _warning_codes(unresolved.warnings)


def test_transparency_report_ordering_caveats_and_deferrals_are_invariant():
    first_warning = AssumptionWarning(
        code="shapley_not_biological_causality",
        message="Shapley attribution metadata states model-output attribution only.",
    )
    second_warning = AssumptionWarning(
        code="gsh_not_pbpk_ode_validated",
        message="GSH redox-capacity model is not PBPK or ODE validated.",
    )

    left = build_transparency_report([first_warning, second_warning])
    right = build_transparency_report([second_warning, first_warning])

    assert [(record.category, record.code) for record in left.records] == [
        (record.category, record.code) for record in right.records
    ]
    assert left.accepted_non_blocking_caveats == []
    assert left.unresolved_blockers == []

    boundary_statements = {boundary.statement for boundary in left.model_boundaries}
    assert any("not PBPK/ODE/GSH-GSSG/Nrf2 validated" in statement for statement in boundary_statements)
    assert any("not proof of biological causality" in statement for statement in boundary_statements)

    queue_codes = {item.code for item in left.sme_review_queue}
    assert {"gsh_not_pbpk_ode_validated", "shapley_not_causality"} <= queue_codes

    deferral_subjects = {item.subject for item in left.deferrals}
    assert "RDKit/ECFP4/Tanimoto affinity fallback" in deferral_subjects
    assert any(item.category is AssumptionCategory.RELEASE_DEFERRAL for item in left.deferrals)
    assert all(item.review_status is SMEReviewStatus.DEFERRED_3_0 for item in left.deferrals)


def test_serializable_records_are_json_friendly_and_preserve_local_utf8_io(tmp_path):
    reaction = MetabolicReaction(
        enzyme="CYP2E1",
        substrate="beta-HCH",
        product="epoxide",
        reaction_role=ReactionRole.UNKNOWN,
        risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.UNKNOWN,
        evidence=EvidenceRecord(
            source="local fixture",
            grade=EvidenceGrade.PLACEHOLDER,
            metadata={"encoding_note": "UTF-8 local fixture"},
        ),
        warnings=[
            AssumptionWarning(
                code="reaction_role_unknown",
                message="Unknown role remains warned and neutral.",
            )
        ],
        metadata={"nested": {"roles": (ReactionRole.UNKNOWN, RiskDirectionIfFluxDecreases.UNKNOWN)}},
    )

    payload = reaction.to_dict()
    assert payload["reaction_role"] == "unknown"
    assert payload["risk_direction_if_flux_decreases"] == "unknown"
    assert payload["evidence"]["grade"] == "placeholder"
    assert payload["warnings"][0]["review_status"] == "unknown"
    assert payload["metadata"]["nested"]["roles"] == ["unknown", "unknown"]

    output_path = tmp_path / "reaction.json"
    output_path.write_text(json.dumps(payload, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    loaded = json.loads(output_path.read_text(encoding="utf-8"))

    assert loaded == payload
