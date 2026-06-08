"""Standalone Shapley and mechanism-interaction attribution for Phase 9.

This module attributes caller-provided scalar model outputs across induction,
competition, and GSH mechanism states. It does not import engines, compute risk
outputs, or make biological causality claims.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import factorial, isfinite
from typing import Callable, Iterable, Mapping

from .interaction_schema import (
    AssumptionWarning,
    EvidenceGrade,
    EvidenceRecord,
    JsonDict,
    SMEReviewStatus,
    SerializableRecord,
    ValueEnum,
)


class MechanismName(ValueEnum):
    """Phase 9 mechanism names for model-output attribution."""

    INDUCTION = "induction"
    COMPETITION = "competition"
    GSH = "gsh"


DEFAULT_MECHANISMS: tuple[MechanismName, MechanismName, MechanismName] = (
    MechanismName.INDUCTION,
    MechanismName.COMPETITION,
    MechanismName.GSH,
)


@dataclass(unsafe_hash=True)
class MechanismState(SerializableRecord):
    """Deterministic binary state for the three Phase 9 mechanisms."""

    induction: bool = False
    competition: bool = False
    gsh: bool = False

    @classmethod
    def from_active(cls, active_mechanisms: Iterable[MechanismName | str]) -> "MechanismState":
        active = {_coerce_mechanism(mechanism) for mechanism in active_mechanisms}
        return cls(
            induction=MechanismName.INDUCTION in active,
            competition=MechanismName.COMPETITION in active,
            gsh=MechanismName.GSH in active,
        )

    @property
    def active_mechanisms(self) -> tuple[MechanismName, ...]:
        return tuple(mechanism for mechanism in DEFAULT_MECHANISMS if self.is_active(mechanism))

    @property
    def key(self) -> str:
        active = self.active_mechanisms
        if not active:
            return "none"
        return "+".join(mechanism.value for mechanism in active)

    def is_active(self, mechanism: MechanismName | str) -> bool:
        mechanism_name = _coerce_mechanism(mechanism)
        if mechanism_name is MechanismName.INDUCTION:
            return self.induction
        if mechanism_name is MechanismName.COMPETITION:
            return self.competition
        return self.gsh


@dataclass
class MechanismStateValue(SerializableRecord):
    """Caller-provided model output for one mechanism state."""

    state: MechanismState
    value: float
    state_key: str | None = None
    metadata: JsonDict | None = None


@dataclass
class ShapleyEffect(SerializableRecord):
    """Shapley main effect for one mechanism."""

    mechanism: MechanismName
    effect: float
    metadata: JsonDict | None = None


@dataclass
class MechanismInteractionTerm(SerializableRecord):
    """Finite-difference singleton, pairwise, or three-way term."""

    mechanisms: tuple[MechanismName, ...]
    order: int
    effect: float
    term_type: str
    metadata: JsonDict | None = None


@dataclass
class MechanismAttributionResult(SerializableRecord):
    """Complete Phase 9 model-output attribution result."""

    mechanisms: list[MechanismName]
    baseline_value: float
    full_value: float
    total_effect: float
    state_values: list[MechanismStateValue]
    singleton_effects: list[MechanismInteractionTerm]
    shapley_main_effects: list[ShapleyEffect]
    pairwise_interactions: list[MechanismInteractionTerm]
    three_way_interaction: MechanismInteractionTerm
    shapley_residual: float
    interaction_reconstruction_residual: float
    tolerance: float
    residuals_are_zero_within_tolerance: bool
    warnings: list[AssumptionWarning]
    evidence: EvidenceRecord | None = None
    metadata: JsonDict | None = None


AttributionWarning = AssumptionWarning
StateValueInput = Mapping[MechanismState | str | tuple[str, ...] | frozenset[str], float] | Iterable[MechanismStateValue]


def generate_mechanism_states() -> list[MechanismState]:
    """Generate the required eight mechanism states in deterministic order."""

    states = [MechanismState()]
    for subset_size in range(1, len(DEFAULT_MECHANISMS) + 1):
        for subset in combinations(DEFAULT_MECHANISMS, subset_size):
            states.append(MechanismState.from_active(subset))
    return states


def validate_eight_state_values(state_values: StateValueInput) -> dict[MechanismState, float]:
    """Validate complete finite scalar outputs for all eight mechanism states."""

    values = _coerce_state_value_map(state_values)
    expected_states = generate_mechanism_states()
    expected = set(expected_states)
    provided = set(values)
    missing = [state.key for state in expected_states if state not in provided]
    if missing:
        raise ValueError(f"Missing Phase 9 mechanism state values: {', '.join(missing)}")
    extras = [state.key for state in provided if state not in expected]
    if extras:
        raise ValueError(f"Unexpected Phase 9 mechanism state values: {', '.join(sorted(extras))}")
    return {state: values[state] for state in expected_states}


def evaluate_mechanism_states(
    evaluator: Callable[[MechanismState], float],
) -> list[MechanismStateValue]:
    """Evaluate all eight states with a caller-provided scalar-output function."""

    evaluated: list[MechanismStateValue] = []
    for state in generate_mechanism_states():
        value = _validate_finite_value(evaluator(state), state.key)
        evaluated.append(MechanismStateValue(state=state, value=value, state_key=state.key))
    return evaluated


def compute_shapley_main_effects(state_values: StateValueInput) -> list[ShapleyEffect]:
    """Compute standard Shapley main effects for the three mechanisms."""

    values = validate_eight_state_values(state_values)
    n_mechanisms = len(DEFAULT_MECHANISMS)
    effects: list[ShapleyEffect] = []
    for mechanism in DEFAULT_MECHANISMS:
        effect = 0.0
        other_mechanisms = [candidate for candidate in DEFAULT_MECHANISMS if candidate is not mechanism]
        for subset_size in range(0, len(other_mechanisms) + 1):
            for subset in combinations(other_mechanisms, subset_size):
                subset_state = MechanismState.from_active(subset)
                with_mechanism_state = MechanismState.from_active((*subset, mechanism))
                weight = (
                    factorial(subset_size)
                    * factorial(n_mechanisms - subset_size - 1)
                    / factorial(n_mechanisms)
                )
                effect += weight * (values[with_mechanism_state] - values[subset_state])
        effects.append(
            ShapleyEffect(
                mechanism=mechanism,
                effect=effect,
                metadata={"formula": "weighted_average_marginal_contribution", "n_mechanisms": n_mechanisms},
            )
        )
    return effects


def compute_explicit_mechanism_interactions(
    state_values: StateValueInput,
) -> tuple[list[MechanismInteractionTerm], list[MechanismInteractionTerm], MechanismInteractionTerm]:
    """Compute singleton, pairwise, and three-way finite-difference terms."""

    values = validate_eight_state_values(state_values)
    empty = MechanismState()
    baseline = values[empty]

    singleton_terms: list[MechanismInteractionTerm] = []
    for mechanism in DEFAULT_MECHANISMS:
        state = MechanismState.from_active((mechanism,))
        singleton_terms.append(
            MechanismInteractionTerm(
                mechanisms=(mechanism,),
                order=1,
                effect=values[state] - baseline,
                term_type="singleton_effect",
                metadata={"state_key": state.key},
            )
        )

    pairwise_terms: list[MechanismInteractionTerm] = []
    for left, right in combinations(DEFAULT_MECHANISMS, 2):
        pair_state = MechanismState.from_active((left, right))
        left_state = MechanismState.from_active((left,))
        right_state = MechanismState.from_active((right,))
        pairwise_terms.append(
            MechanismInteractionTerm(
                mechanisms=(left, right),
                order=2,
                effect=values[pair_state] - values[left_state] - values[right_state] + baseline,
                term_type="pairwise_interaction",
                metadata={"state_key": pair_state.key, "formula": "mobius_pairwise_finite_difference"},
            )
        )

    induction = MechanismState.from_active((MechanismName.INDUCTION,))
    competition = MechanismState.from_active((MechanismName.COMPETITION,))
    gsh = MechanismState.from_active((MechanismName.GSH,))
    induction_competition = MechanismState.from_active((MechanismName.INDUCTION, MechanismName.COMPETITION))
    induction_gsh = MechanismState.from_active((MechanismName.INDUCTION, MechanismName.GSH))
    competition_gsh = MechanismState.from_active((MechanismName.COMPETITION, MechanismName.GSH))
    full = MechanismState.from_active(DEFAULT_MECHANISMS)
    three_way = MechanismInteractionTerm(
        mechanisms=DEFAULT_MECHANISMS,
        order=3,
        effect=(
            values[full]
            - values[induction_competition]
            - values[induction_gsh]
            - values[competition_gsh]
            + values[induction]
            + values[competition]
            + values[gsh]
            - baseline
        ),
        term_type="three_way_interaction",
        metadata={"state_key": full.key, "formula": "mobius_three_way_finite_difference"},
    )
    return singleton_terms, pairwise_terms, three_way


def compute_mechanism_attribution(
    state_values: StateValueInput,
    *,
    tolerance: float = 1e-9,
    evidence: EvidenceRecord | None = None,
    metadata: JsonDict | None = None,
) -> MechanismAttributionResult:
    """Compute Phase 9 Shapley and explicit interaction attribution."""

    if tolerance < 0.0 or not isfinite(float(tolerance)):
        raise ValueError("tolerance must be a non-negative finite value")

    values = validate_eight_state_values(state_values)
    state_value_records = [
        MechanismStateValue(state=state, value=values[state], state_key=state.key)
        for state in generate_mechanism_states()
    ]
    baseline_state = MechanismState()
    full_state = MechanismState.from_active(DEFAULT_MECHANISMS)
    baseline_value = values[baseline_state]
    full_value = values[full_state]
    total_effect = full_value - baseline_value

    shapley_effects = compute_shapley_main_effects(values)
    singleton_terms, pairwise_terms, three_way_term = compute_explicit_mechanism_interactions(values)

    shapley_residual = total_effect - sum(effect.effect for effect in shapley_effects)
    interaction_reconstruction = (
        sum(term.effect for term in singleton_terms)
        + sum(term.effect for term in pairwise_terms)
        + three_way_term.effect
    )
    interaction_residual = total_effect - interaction_reconstruction
    shapley_residual = _zero_if_close(shapley_residual, tolerance)
    interaction_residual = _zero_if_close(interaction_residual, tolerance)
    residuals_ok = abs(shapley_residual) <= tolerance and abs(interaction_residual) <= tolerance

    warnings: list[AssumptionWarning] = []
    if not residuals_ok:
        warnings.append(
            _warning(
                "mechanism_attribution_residual_nonzero",
                "Attribution residual exceeded tolerance; complete state values may be numerically inconsistent.",
                field="residuals_are_zero_within_tolerance",
            )
        )

    result_metadata: JsonDict = {
        "phase": "phase_9_shapley_interactions",
        "attribution_type": "model_output_attribution",
        "not_biological_causality": True,
        "mechanism_count": len(DEFAULT_MECHANISMS),
        "state_count": len(state_value_records),
        "state_generation": "all_binary_states_for_induction_competition_gsh",
        "residual_policy": "no_unexplained_residual_term",
        "engine_integration": False,
        "public_risk_output": "not_produced_or_modified",
        "phase10_behavior": False,
        "uncertainty_propagation": False,
    }
    if metadata:
        result_metadata["caller_metadata"] = dict(metadata)

    return MechanismAttributionResult(
        mechanisms=list(DEFAULT_MECHANISMS),
        baseline_value=baseline_value,
        full_value=full_value,
        total_effect=total_effect,
        state_values=state_value_records,
        singleton_effects=singleton_terms,
        shapley_main_effects=shapley_effects,
        pairwise_interactions=pairwise_terms,
        three_way_interaction=three_way_term,
        shapley_residual=shapley_residual,
        interaction_reconstruction_residual=interaction_residual,
        tolerance=float(tolerance),
        residuals_are_zero_within_tolerance=residuals_ok,
        warnings=warnings,
        evidence=evidence or _default_evidence(),
        metadata=result_metadata,
    )


def _coerce_state_value_map(state_values: StateValueInput) -> dict[MechanismState, float]:
    values: dict[MechanismState, float] = {}
    if isinstance(state_values, Mapping):
        iterator = state_values.items()
        for key, value in iterator:
            state = _coerce_state(key)
            if state in values:
                raise ValueError(f"Duplicate Phase 9 mechanism state value: {state.key}")
            values[state] = _validate_finite_value(value, state.key)
        return values

    for item in state_values:
        if not isinstance(item, MechanismStateValue):
            raise ValueError("State-value iterables must contain MechanismStateValue records.")
        state = _coerce_state(item.state)
        if state in values:
            raise ValueError(f"Duplicate Phase 9 mechanism state value: {state.key}")
        values[state] = _validate_finite_value(item.value, state.key)
    return values


def _coerce_state(value: MechanismState | str | tuple[str, ...] | frozenset[str]) -> MechanismState:
    if isinstance(value, MechanismState):
        return value
    if isinstance(value, str):
        if value == "none" or value.strip() == "":
            return MechanismState()
        return MechanismState.from_active(part.strip() for part in value.split("+") if part.strip())
    return MechanismState.from_active(value)


def _coerce_mechanism(value: MechanismName | str) -> MechanismName:
    if isinstance(value, MechanismName):
        return value
    try:
        return MechanismName(str(value))
    except ValueError as exc:
        raise ValueError(f"Unknown Phase 9 mechanism: {value}") from exc


def _validate_finite_value(value: float, state_key: str) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"State value for {state_key} must be numeric.") from exc
    if not isfinite(numeric):
        raise ValueError(f"State value for {state_key} must be finite.")
    return numeric


def _zero_if_close(value: float, tolerance: float) -> float:
    if abs(value) <= tolerance:
        return 0.0
    return value


def _default_evidence() -> EvidenceRecord:
    return EvidenceRecord(
        source="ExposoGraph Phase 9 local attribution layer",
        grade=EvidenceGrade.PLACEHOLDER,
        confidence="deterministic_model_output_attribution",
        notes="Shapley and finite-difference attribution over caller-provided model outputs; not biological causality.",
        metadata={"phase": "phase_9"},
    )


def _warning(code: str, message: str, *, field: str | None = None) -> AssumptionWarning:
    return AssumptionWarning(
        code=code,
        message=message,
        field=field,
        review_status=SMEReviewStatus.UNKNOWN,
    )
