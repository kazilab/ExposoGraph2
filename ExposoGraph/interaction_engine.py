"""Multi-carcinogen interaction modeling.

Models competitive enzyme inhibition, glutathione depletion, and
lifestyle-driven enzyme induction for combined carcinogen exposure profiles.

Data source: ``interaction_parameters.json`` in ``ExposoGraph/data/``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping

from .effective_burden import (
    EffectiveBurdenInput,
    EffectiveBurdenResult,
    GSHBurdenCouplingInput,
    compute_effective_carcinogenic_burden,
)
from .endpoint_toxic_flux import (
    EndpointFluxInput,
    EndpointToxicFluxResult,
    interpret_endpoint_toxic_flux,
)
from .flux_engine import apply_kinetic_modifier_once
from .gsh_redox_capacity import (
    GSHModelVersion,
    GSHRedoxCapacityInput,
    compute_gsh_redox_capacity,
)
from .interaction_schema import (
    ApplicabilityDomain,
    AssumptionWarning,
    ConcentrationBasis,
    InhibitionMode,
    ReactionRole,
    RiskDirectionIfFluxDecreases,
    SMEReviewStatus,
    enum_from_value,
    to_jsonable,
)
from .kinetic_resolver import get_ki, resolve_reversible_inhibition
from .mechanism_attribution import compute_mechanism_attribution, generate_mechanism_states
from .model_transparency import build_transparency_report
from .parameter_resolution import InhibitionResolutionStatus
from .reaction_role_semantics import get_default_reaction_role_registry

# ── Dataclasses ────────────────────────────────────────────────────────────


@dataclass
class EnzymeInductionProfile:
    """Fold-induction multipliers for CYP enzymes based on lifestyle."""

    enzyme_folds: dict[str, float]
    active_inducers: list[str] = field(default_factory=list)


@dataclass
class SubstrateFluxChange:
    """Flux change for a single substrate under competitive inhibition."""

    single_flux: float
    competitive_flux: float
    flux_change_fraction: float
    inhibition_term: float
    activated_product_flux: float
    Km_uM: float
    concentration_uM: float
    product: str
    product_carcinogenic: bool
    baseline_flux: float | None = None
    kinetic_modifier: float | None = None
    modified_flux: float | None = None
    kinetic_resolution_status: str = "not_assessable"
    kinetic_warning_codes: tuple[str, ...] = field(default_factory=tuple)
    inhibition_mode: str = "competitive"
    centralized_resolver_used: bool = False
    modifier_applied_once: bool = False
    discrepancy_classification: str = "not_compared"
    biological_output: dict[str, Any] | None = None
    aggregate_resolution: dict[str, Any] | None = None


@dataclass
class CompetitiveInhibitionResult:
    """Per-substrate flux changes under competitive inhibition for one enzyme."""

    enzyme: str
    substrates: dict[str, SubstrateFluxChange]


@dataclass
class GSHStatus:
    """Glutathione pool status under multi-carcinogen exposure."""

    baseline_gsh_mM: float
    steady_state_gsh_mM: float
    fraction_normal: float
    consumption_rate_umol_h_g: float
    synthesis_rate_umol_h_g: float
    net_rate_umol_h_g: float
    consumption_exceeds_synthesis: bool
    tipping_point_reached: bool
    tipping_point_multiplier: float
    impaired_pathways: list[str]
    individual_contributions: dict[str, Any]
    time_to_depletion_h: float | None
    tissue: str
    model_version: str = GSHModelVersion.LEGACY_DETACHED_GSH_PENALTY.value
    redox_capacity_ratio: float | None = None
    detox_penalty_multiplier: float | None = None
    warnings: list[dict[str, Any]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class MechanismResolvedRisk:
    """Per-carcinogen mechanism calculation used for adjusted relative risk."""

    carcinogen: str
    baseline_relative_risk: float
    induction_multiplier: float
    inhibition_burden_multiplier: float
    activation_burden_ratio: float
    detox_failure_ratio: float
    matrix_gsh_penalty: float
    gsh_pool_penalty: float
    susceptibility_modifier: float
    susceptibility_applied_in: str | None
    final_mechanism_multiplier: float
    adjusted_relative_risk: float
    inhibition_status: str
    review_required: bool
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return _json_sanitize(
            {
                "carcinogen": self.carcinogen,
                "baseline_relative_risk": self.baseline_relative_risk,
                "induction_multiplier": self.induction_multiplier,
                "inhibition_burden_multiplier": self.inhibition_burden_multiplier,
                "activation_burden_ratio": self.activation_burden_ratio,
                "detox_failure_ratio": self.detox_failure_ratio,
                "matrix_gsh_penalty": self.matrix_gsh_penalty,
                "gsh_pool_penalty": self.gsh_pool_penalty,
                "susceptibility_modifier": self.susceptibility_modifier,
                "susceptibility_applied_in": self.susceptibility_applied_in,
                "final_mechanism_multiplier": self.final_mechanism_multiplier,
                "adjusted_relative_risk": self.adjusted_relative_risk,
                "inhibition_status": self.inhibition_status,
                "review_required": self.review_required,
                "warnings": list(self.warnings),
                "provenance": deepcopy(self.provenance),
            }
        )


@dataclass
class _InhibitionBurdenResolution:
    burden_multiplier: float
    activation_burden_ratio: float
    detox_failure_ratio: float
    endpoint_toxic_flux_ratio: float | None
    status: str
    review_required: bool
    enzyme: str | None = None
    flux_substrate: str | None = None
    interpretation_substrate: str | None = None
    tissue: str | None = None
    flux_ratio: float | None = None
    warnings: list[str] = field(default_factory=list)
    reaction_role_annotation: Any | None = None
    reaction_role_block: dict[str, Any] | None = None
    endpoint_toxic_flux_result: EndpointToxicFluxResult | None = None
    effective_burden_result: EffectiveBurdenResult | None = None
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def reaction_role_interpretation(self) -> dict[str, Any]:
        if self.reaction_role_block is not None:
            return deepcopy(self.reaction_role_block)
        if self.status == "mechanism_absent":
            directional = "not_applicable_mechanism_absent"
        elif self.review_required:
            directional = "withheld_review_required"
        else:
            directional = "no_directional_change"
        return {
            "role": str(ReactionRole.UNKNOWN),
            "directional_interpretation": directional,
            "risk_direction_if_flux_decreases": str(RiskDirectionIfFluxDecreases.UNKNOWN),
            "annotation_record_id": self.provenance.get("annotation_record_id"),
            "review_required": self.review_required,
            "warnings": list(self.warnings),
            "sme_notes": [],
            "status": self.status,
        }

    @property
    def endpoint_toxic_flux(self) -> dict[str, Any] | None:
        if self.endpoint_toxic_flux_result is None:
            if self.enzyme is None or self.interpretation_substrate is None:
                return None
            reaction_role = self.reaction_role_interpretation
            endpoint_block = {
                "enzyme": self.enzyme,
                "substrate": self.interpretation_substrate,
                "tissue": self.tissue,
                "endpoint": None,
                "flux_ratio": self.flux_ratio,
                "endpoint_toxic_flux_ratio": self.endpoint_toxic_flux_ratio,
                "burden_multiplier": self.burden_multiplier,
                "activation_burden_ratio": self.activation_burden_ratio,
                "detox_failure_ratio": self.detox_failure_ratio,
                "reaction_role": reaction_role.get("role"),
                "risk_direction_if_flux_decreases": reaction_role.get(
                    "risk_direction_if_flux_decreases"
                ),
                "annotation_record_id": reaction_role.get("annotation_record_id"),
                "warnings": list(self.warnings),
                "role_warnings": deepcopy(reaction_role.get("warnings", [])),
                "sme_notes": deepcopy(reaction_role.get("sme_notes", [])),
                "metadata": {
                    "source": "interaction_engine.selected_inhibition_resolution",
                    "kinetic_mechanism_state": self.status,
                    "flux_substrate": self.flux_substrate,
                    "neutral_authoritative_resolution": True,
                },
                "status": self.status,
            }
        else:
            endpoint_block = self.endpoint_toxic_flux_result.to_dict()
        endpoint_block["review_required"] = self.review_required
        return endpoint_block

    @property
    def effective_burden(self) -> dict[str, Any] | None:
        if self.effective_burden_result is None:
            if self.enzyme is None or self.interpretation_substrate is None:
                return None
            reaction_role = self.reaction_role_interpretation
            effective_block = {
                "activation_burden_ratio": self.activation_burden_ratio,
                "detox_failure_ratio": self.detox_failure_ratio,
                "gsh_detox_penalty_ratio": 1.0,
                "susceptibility_modifier": 1.0,
                "effective_carcinogenic_burden_ratio": self.burden_multiplier,
                "gsh_relevant": False,
                "gsh_consumption_load": None,
                "gsh_consumption_load_scaled": None,
                "gsh_fraction": None,
                "redox_capacity_ratio": None,
                "model_boundary": "interaction_matrix_selected_inhibition_resolution",
                "warnings": list(self.warnings),
                "evidence": None,
                "gsh_coupling_result": None,
                "reaction_role": reaction_role.get("role"),
                "risk_direction_if_flux_decreases": reaction_role.get(
                    "risk_direction_if_flux_decreases"
                ),
                "annotation_record_id": reaction_role.get("annotation_record_id"),
                "role_warnings": deepcopy(reaction_role.get("warnings", [])),
                "sme_notes": deepcopy(reaction_role.get("sme_notes", [])),
                "metadata": {
                    "source": "interaction_engine.selected_inhibition_resolution",
                    "kinetic_mechanism_state": self.status,
                    "flux_substrate": self.flux_substrate,
                    "neutral_authoritative_resolution": True,
                },
                "status": self.status,
            }
        else:
            effective_block = self.effective_burden_result.to_dict()
        effective_block["review_required"] = self.review_required
        return effective_block

    def to_provenance(self) -> dict[str, Any]:
        return _json_sanitize(
            {
                "inhibition": {
                    "enzyme": self.enzyme,
                    "flux_substrate": self.flux_substrate,
                    "interpretation_substrate": self.interpretation_substrate,
                    "flux_ratio": self.flux_ratio,
                    "endpoint_toxic_flux_ratio": self.endpoint_toxic_flux_ratio,
                    "status": self.status,
                    "review_required": self.review_required,
                    "warnings": list(self.warnings),
                    "reaction_role_interpretation": deepcopy(self.reaction_role_interpretation),
                    "reaction_role_annotation": (
                        self.reaction_role_annotation.to_dict()
                        if self.reaction_role_annotation is not None
                        and hasattr(self.reaction_role_annotation, "to_dict")
                        else deepcopy(self.reaction_role_annotation)
                    ),
                    "endpoint_toxic_flux": deepcopy(self.endpoint_toxic_flux),
                    "effective_burden": deepcopy(self.effective_burden),
                    **deepcopy(self.provenance),
                }
            }
        )


@dataclass
class CriticalInteraction:
    """A genotype-specific critical interaction warning."""

    interaction: str
    severity: str
    mechanism: str
    affected_carcinogens: list[str]
    genotype_amplification: float
    clinical_note: str


@dataclass
class InteractionMatrixResult:
    """Full multi-carcinogen interaction analysis."""

    individual_risks: dict[str, float]
    interaction_adjusted_risks: dict[str, float]
    synergy_matrix: dict[str, float]
    gsh_status: GSHStatus
    induction_effects: EnzymeInductionProfile
    competitive_effects: dict[str, CompetitiveInhibitionResult]
    total_independent_risk: float
    total_interaction_risk: float
    interaction_factor: float
    tissue: str
    genotypes: dict[str, str]
    lifestyle: dict[str, Any]
    summary: str
    mechanism_attribution: dict[str, Any] | None = None
    mechanism_resolved_risks: dict[str, MechanismResolvedRisk] = field(default_factory=dict)


@dataclass
class SynergyDecomposition:
    """Pairwise synergy attribution with compatibility singleton fields."""

    pair: str
    composite: float
    delta_comp: float
    delta_gsh: float
    delta_ind: float
    additive_estimate: float
    residual: float
    shapley_decomposition: dict[str, Any] = field(default_factory=dict)
    main_effects: dict[str, float] = field(default_factory=dict)
    pairwise_interactions: dict[str, float] = field(default_factory=dict)
    three_way_interaction: float = 0.0
    reconstruction_residual: float = 0.0
    shapley_residual: float = 0.0
    residuals_are_zero_within_tolerance: bool = True
    state_values: dict[str, float] = field(default_factory=dict)
    decomposition_basis: str = "eight_state_shapley"
    residual_policy: str = "numerical_reconstruction_check_only"
    compatibility_fields: dict[str, Any] = field(default_factory=dict)
    dominant_mechanism: str = "near-additive"


@dataclass
class SynergyConfidenceInterval:
    """Bootstrap confidence intervals for a carcinogen-pair synergy."""

    pair: str
    n_iterations: int
    composite_mean: float
    composite_ci95: tuple[float, float]
    delta_comp_mean: float
    delta_comp_ci95: tuple[float, float]
    delta_gsh_mean: float
    delta_gsh_ci95: tuple[float, float]
    delta_ind_mean: float
    delta_ind_ci95: tuple[float, float]


# ── Data loading ───────────────────────────────────────────────────────────

_INTERACTION_PARAMS_FILE = Path(__file__).parent / "data" / "interaction_parameters.json"
_PROVENANCE_FILE = Path(__file__).parent / "data" / "parameter_provenance.json"
_INTERACTION_CACHE: dict[str, Any] | None = None
_PROVENANCE_CACHE: dict[str, Any] | None = None


def _load_interaction_params() -> dict[str, Any]:
    """Lazy-load and cache the interaction parameters JSON."""
    global _INTERACTION_CACHE
    if _INTERACTION_CACHE is None:
        if not _INTERACTION_PARAMS_FILE.exists():
            raise FileNotFoundError(
                f"Interaction parameters not found at {_INTERACTION_PARAMS_FILE}. "
                "Ensure interaction_parameters.json is in the data/ directory."
            )
        with open(_INTERACTION_PARAMS_FILE, "r") as fh:
            _INTERACTION_CACHE = json.load(fh)
    return _INTERACTION_CACHE


def _get_interaction_params() -> dict[str, Any]:
    """Return the cached interaction parameters."""
    return _load_interaction_params()


def get_parameter_provenance() -> dict[str, Any]:
    """Return a defensive copy of the kinetic-parameter provenance catalog.

    Includes Km/Vmax source references, confidence grades, and an explicit
    ``ki_status`` ("curated" vs "assumed_equal_km") for each enzyme-substrate
    pair in the interaction engine.
    """
    global _PROVENANCE_CACHE
    if _PROVENANCE_CACHE is None:
        if not _PROVENANCE_FILE.exists():
            raise FileNotFoundError(
                f"Provenance catalog not found at {_PROVENANCE_FILE}."
            )
        with open(_PROVENANCE_FILE, "r") as fh:
            _PROVENANCE_CACHE = json.load(fh)
    return deepcopy(_PROVENANCE_CACHE)


def get_interaction_source_catalog() -> list[dict[str, Any]]:
    """Return the prioritized source catalog for interaction-model curation."""
    metadata = _get_interaction_params().get("_metadata", {})
    return deepcopy(metadata.get("primary_data_sources", []))


def get_interaction_expansion_backlog() -> dict[str, Any]:
    """Return the structured interaction-expansion backlog from JSON metadata."""
    metadata = _get_interaction_params().get("_metadata", {})
    return deepcopy(metadata.get("expansion_backlog", {}))


def assumed_ki_pairs() -> list[tuple[str, str]]:
    """List enzyme-substrate pairs whose Ki is assumed equal to Km."""
    provenance = get_parameter_provenance().get("pairs", {})
    pairs: list[tuple[str, str]] = []
    for enzyme, substrates in provenance.items():
        for sub_name, entry in substrates.items():
            if entry.get("ki_status") == "assumed_equal_km":
                pairs.append((enzyme, sub_name))
    return pairs


# ── Public constants ───────────────────────────────────────────────────────

BASELINE_RISK_SCORES: dict[str, float] = {
    "PAH": 25.0,
    "HCA": 18.0,
    "NNK": 30.0,
    "benzene": 20.0,
    "NDMA": 22.0,
    "formaldehyde": 10.0,
    "chromium_VI": 28.0,
    "arsenic": 22.0,
    "AFB1": 35.0,
    "acetaldehyde": 15.0,
    "cadmium": 18.0,
    "vinyl_chloride": 20.0,
    "acrolein": 8.0,
}

CARCINOGEN_ENZYME_MAP: dict[str, list[str]] = {
    "PAH": ["CYP1A1", "CYP1B1"],
    "HCA": ["CYP1A2"],
    "NNK": ["CYP1A2", "CYP2E1"],
    "benzene": ["CYP2E1"],
    "NDMA": ["CYP2E1"],
    "vinyl_chloride": ["CYP2E1"],
    "AFB1": ["CYP3A4"],
    "acetaldehyde": [],
    "formaldehyde": [],
    "chromium_VI": [],
    "arsenic": [],
    "cadmium": [],
    "acrolein": [],
}

CARCINOGEN_GSH_DETOX: dict[str, str | None] = {
    "PAH": "PAH_GSTM1",
    "HCA": None,
    "NNK": None,
    "benzene": None,
    "NDMA": None,
    "vinyl_chloride": None,
    "AFB1": "BPDE_conjugation",
    "acetaldehyde": None,
    "formaldehyde": None,
    "chromium_VI": "chromium_VI",
    "arsenic": "arsenic_methylation",
    "cadmium": "cadmium",
    "acrolein": "acrolein",
}


# ── Private constants ──────────────────────────────────────────────────────

_DEFAULT_ENZYMES = ["CYP1A1", "CYP1A2", "CYP1B1", "CYP2E1", "CYP3A4"]

_CARCINOGEN_ALIASES: dict[str, str] = {
    "pah": "PAH",
    "bap": "PAH",
    "benzoapyrene": "PAH",
    "hca": "HCA",
    "phip": "HCA",
    "nnk": "NNK",
    "tobacconitrosamines": "NNK",
    "benzene": "benzene",
    "ndma": "NDMA",
    "formaldehyde": "formaldehyde",
    "chromium": "chromium_VI",
    "chromiumvi": "chromium_VI",
    "chromium_vi": "chromium_VI",
    "cr": "chromium_VI",
    "arsenic": "arsenic",
    "afb1": "AFB1",
    "aflatoxin": "AFB1",
    "aflatoxinb1": "AFB1",
    "acetaldehyde": "acetaldehyde",
    "aldehyde": "acetaldehyde",
    "cadmium": "cadmium",
    "vinylchloride": "vinyl_chloride",
    "vinyl_chloride": "vinyl_chloride",
    "acrolein": "acrolein",
    "ethanol": "ethanol",
}

_REFERENCE_CONCENTRATIONS_UM: dict[str, float] = {
    "benzene": 10.0,
    "NDMA": 0.5,
    "vinyl_chloride": 5.0,
    "ethanol": 500.0,
    "PAH": 0.1,
    "HCA": 1.0,
    "AFB1": 1.0,
    "NNK": 0.5,
}

_GSH_REFERENCE_RATES: dict[str, float] = {
    "PAH": 0.5,
    "chromium_VI": 0.3,
    "arsenic": 0.2,
    "cadmium": 0.1,
    "acrolein": 0.4,
    "ethanol": 0.8,
}

_SEVERITY_RANK = {"CRITICAL": 3, "HIGH": 2, "MODERATE": 1}


EXPOSURE_PROFILES: dict[str, dict[str, Any]] = {
    "smoker": {
        "_description": (
            "Active cigarette smoker (1 pack/day). Includes PAH from tobacco smoke, "
            "tobacco-specific nitrosamines (NNK), benzene, formaldehyde, cadmium, and "
            "acrolein. CYP1A2 induction amplifies HCA activation."
        ),
        "exposure": {
            "PAH": 3.0,
            "NNK": 4.0,
            "HCA": 1.5,
            "benzene": 6.0,
            "formaldehyde": 2.0,
            "cadmium": 2.0,
            "acrolein": 5.0,
        },
        "lifestyle": {"smoking": True, "pack_years": 20},
    },
    "heavy_drinker": {
        "_description": (
            "Chronic heavy drinker (>4 drinks/day). Elevated acetaldehyde, NDMA activation "
            "via CYP2E1, and GSH depletion from ethanol-induced ROS."
        ),
        "exposure": {
            "acetaldehyde": 3.0,
            "NDMA": 1.5,
            "benzene": 1.5,
            "ethanol": 8.0,
            "acrolein": 1.5,
        },
        "lifestyle": {"alcohol_heavy": True, "chronic_alcohol": True},
    },
    "smoker_heavy_drinker": {
        "_description": (
            "Critical interaction scenario: combined smoking + heavy drinking. Smoking "
            "induces CYP1A2 while alcohol induces CYP2E1 and depletes GSH."
        ),
        "exposure": {
            "PAH": 3.0,
            "NNK": 4.0,
            "HCA": 2.0,
            "benzene": 6.0,
            "NDMA": 2.0,
            "formaldehyde": 2.0,
            "cadmium": 2.0,
            "acrolein": 6.0,
            "acetaldehyde": 3.0,
            "ethanol": 8.0,
        },
        "lifestyle": {
            "smoking": True,
            "alcohol_heavy": True,
            "chronic_alcohol": True,
            "pack_years": 25,
        },
    },
    "industrial_worker": {
        "_description": (
            "Industrial worker with occupational chromium(VI), benzene, and formaldehyde "
            "exposure. Very high GSH consumption from Cr(VI)."
        ),
        "exposure": {
            "chromium_VI": 10.0,
            "benzene": 10.0,
            "formaldehyde": 5.0,
            "PAH": 2.0,
            "cadmium": 3.0,
        },
        "lifestyle": {"smoking": False, "alcohol_heavy": False},
    },
    "smoker_industrial_worker": {
        "_description": (
            "Worst case for GSH depletion: industrial worker who smokes. Smoking PAH plus "
            "industrial Cr(VI) and benzene creates three simultaneous GSH consumers."
        ),
        "exposure": {
            "PAH": 5.0,
            "NNK": 4.0,
            "benzene": 15.0,
            "chromium_VI": 10.0,
            "formaldehyde": 6.0,
            "cadmium": 4.0,
            "acrolein": 5.0,
            "HCA": 1.5,
            "ethanol": 0.5,
        },
        "lifestyle": {"smoking": True, "pack_years": 30},
    },
    "moderate_drinker": {
        "_description": (
            "Moderate drinker (1-2 drinks/day). Some CYP2E1 induction (~2x), mild GSH "
            "depletion, and mildly elevated NDMA risk."
        ),
        "exposure": {
            "acetaldehyde": 1.5,
            "NDMA": 1.2,
            "benzene": 1.0,
            "ethanol": 3.0,
        },
        "lifestyle": {"alcohol_moderate": True},
    },
    "smoker_moderate_drinker": {
        "_description": (
            "Smoker + moderate drinker — profile for patient JHBUI-10030. Moderate CYP2E1 "
            "induction plus full CYP1A2 induction from smoking."
        ),
        "exposure": {
            "PAH": 3.0,
            "NNK": 4.0,
            "HCA": 1.8,
            "benzene": 6.0,
            "NDMA": 1.5,
            "formaldehyde": 2.0,
            "cadmium": 2.0,
            "acrolein": 5.0,
            "acetaldehyde": 1.5,
            "ethanol": 3.0,
        },
        "lifestyle": {"smoking": True, "alcohol_moderate": True, "pack_years": 20},
    },
    "JHBUI_10030": {
        "_description": (
            "Patient JHBUI-10030: smoker + moderate drinker profile from the clinical "
            "ExposoGraph reference set."
        ),
        "exposure": {
            "PAH": 3.0,
            "NNK": 4.0,
            "HCA": 1.8,
            "benzene": 6.0,
            "NDMA": 1.5,
            "formaldehyde": 2.0,
            "cadmium": 2.0,
            "acrolein": 5.0,
            "acetaldehyde": 1.5,
            "ethanol": 3.0,
        },
        "lifestyle": {"smoking": True, "alcohol_moderate": True, "pack_years": 20},
        "genotypes": {
            "GSTM1": "null",
            "CYP2E1": "NM",
            "CYP1A2": "NM",
            "NAT2": "intermediate",
        },
    },
}


# ── Private helpers ────────────────────────────────────────────────────────


def _round(value: float, digits: int) -> float:
    """Round floats consistently for stable public results."""
    return round(float(value), digits)


def _canonical_carcinogen_key(name: str) -> str:
    """Normalize supported carcinogen aliases to canonical keys."""
    if name in BASELINE_RISK_SCORES or name == "ethanol":
        return name
    cleaned = name.strip().replace(" ", "").replace("-", "").replace("[", "").replace("]", "")
    if cleaned in _CARCINOGEN_ALIASES:
        return _CARCINOGEN_ALIASES[cleaned]
    lowered = cleaned.lower()
    return _CARCINOGEN_ALIASES.get(lowered, name)


def _normalize_exposure_profile(exposure_profile: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize carcinogen keys while preserving explicit rate keys."""
    normalized: dict[str, Any] = {}
    for key, value in exposure_profile.items():
        if key.endswith("_umol_h_g"):
            normalized[key] = value
            continue
        normalized[_canonical_carcinogen_key(key)] = value
    return normalized


def _extract_exposure_multiplier(value: float | dict[str, Any]) -> float:
    """Return an exposure multiplier from a scalar or mapping payload."""
    if isinstance(value, dict):
        if "exposure_multiplier" in value:
            return float(value["exposure_multiplier"])
        if "multiplier" in value:
            return float(value["multiplier"])
        return 1.0
    return float(value)


def _extract_concentration_uM(carcinogen: str, value: float | dict[str, Any]) -> float:
    """Return a substrate concentration for competition modeling."""
    if isinstance(value, dict) and "concentration_uM" in value:
        return float(value["concentration_uM"])
    ref = _REFERENCE_CONCENTRATIONS_UM.get(carcinogen, 1.0)
    return _extract_exposure_multiplier(value) * ref


def _collect_known_enzymes(params: dict[str, Any]) -> list[str]:
    """Collect the enzyme names that should always receive baseline folds."""
    enzymes = set(_DEFAULT_ENZYMES)
    for section in ("smoking", "chronic_alcohol", "TCDD_dioxin", "obesity_insulin_resistance"):
        entries = params["enzyme_induction"].get(section, {})
        for enzyme in entries:
            if not enzyme.startswith("_"):
                enzymes.add(enzyme)
    return sorted(enzymes)


def _apply_enzyme_folds(
    enzyme_folds: dict[str, float],
    config: dict[str, Any],
) -> None:
    """Apply fold-induction values from a parameter section."""
    for enzyme, values in config.items():
        if enzyme.startswith("_"):
            continue
        fold = float(values.get("fold_induction", 1.0))
        current = enzyme_folds.get(enzyme, 1.0)
        enzyme_folds[enzyme] = min(current, fold) if fold < 1.0 else max(current, fold)


def _genotype_activity_multiplier(gene: str, phenotype: str | None) -> float:
    """Return the activity multiplier defined in the interaction JSON."""
    if not phenotype:
        return 1.0
    params = _get_interaction_params()["genotype_modifiers"].get(gene, {})
    if phenotype in params:
        return float(params[phenotype].get("activity_multiplier", 1.0))
    lowered = phenotype.lower()
    for key, values in params.items():
        if key.lower() == lowered:
            return float(values.get("activity_multiplier", 1.0))
    return 1.0


_PULMONARY_TISSUES: frozenset[str] = frozenset(
    {"lung", "bronchus", "nasal", "nasopharynx", "trachea"}
)


def _build_competitive_substrates(
    exposure_profile: dict[str, Any],
    *,
    tissue: str = "Liver",
) -> dict[str, dict[str, float]]:
    """Build the substrate concentration maps used by the interaction model.

    When ``tissue`` is pulmonary, CYP2A13 and CYP2F1 are added as additional
    benzene / NNK activators to ambient (non-occupational) benzene metabolism
    may be dominated by pulmonary
    CYP2A13 / CYP2F1 rather than hepatic CYP2E1.
    """
    substrates: dict[str, dict[str, float]] = {}
    tissue_lower = str(tissue).lower() if tissue else "liver"
    is_pulmonary = tissue_lower in _PULMONARY_TISSUES

    cyp2e1: dict[str, float] = {}
    for key in ("benzene", "NDMA", "vinyl_chloride", "ethanol"):
        if key in exposure_profile:
            cyp2e1[key] = _extract_concentration_uM(key, exposure_profile[key])
    if len(cyp2e1) > 1:
        substrates["CYP2E1"] = cyp2e1

    cyp1a1: dict[str, float] = {}
    if "PAH" in exposure_profile:
        cyp1a1["BaP"] = _extract_concentration_uM("PAH", exposure_profile["PAH"])
    if "HCA" in exposure_profile:
        cyp1a1["PhIP"] = _extract_concentration_uM("HCA", exposure_profile["HCA"])
    if len(cyp1a1) > 1:
        substrates["CYP1A1"] = cyp1a1

    if is_pulmonary:
        cyp2a13: dict[str, float] = {}
        if "benzene" in exposure_profile:
            cyp2a13["benzene"] = _extract_concentration_uM(
                "benzene", exposure_profile["benzene"]
            )
        if "NNK" in exposure_profile:
            cyp2a13["NNK"] = _extract_concentration_uM(
                "NNK", exposure_profile["NNK"]
            )
        if len(cyp2a13) >= 1:
            substrates["CYP2A13"] = cyp2a13

        cyp2f1: dict[str, float] = {}
        if "benzene" in exposure_profile:
            cyp2f1["benzene"] = _extract_concentration_uM(
                "benzene", exposure_profile["benzene"]
            )
        if len(cyp2f1) >= 1:
            substrates["CYP2F1"] = cyp2f1

    return substrates


def _to_gsh_rate_map(exposure_profile: dict[str, Any]) -> dict[str, float]:
    """Normalize exposure input into rates consumable by the GSH model."""
    rate_map: dict[str, float] = {}
    for key, value in exposure_profile.items():
        if key.endswith("_umol_h_g"):
            scalar = float(value["flux_umol_h_g"] if isinstance(value, dict) else value)
            rate_map[key] = scalar
            continue

        canonical = _canonical_carcinogen_key(key)
        if canonical not in _GSH_REFERENCE_RATES:
            continue
        if isinstance(value, dict) and "flux_umol_h_g" in value:
            rate_map[f"{canonical}_umol_h_g"] = float(value["flux_umol_h_g"])
        else:
            rate_map[f"{canonical}_umol_h_g"] = _extract_exposure_multiplier(value) * _GSH_REFERENCE_RATES[canonical]
    return rate_map


def _gsh_consumer_key(base_key: str) -> str:
    return {
        "PAH": "PAH_GSTM1",
        "chromium_VI": "chromium_VI",
        "arsenic": "arsenic_methylation",
        "cadmium": "cadmium",
        "acrolein": "acrolein",
        "ethanol": "ethanol_ROS",
        "BPDE": "BPDE_conjugation",
        "acetaminophen": "acetaminophen_NAPQI",
    }.get(base_key, base_key)


def _warning_dict(code: str, message: str, *, field: str | None = None) -> dict[str, Any]:
    payload: dict[str, Any] = {"code": code, "message": message}
    if field is not None:
        payload["field"] = field
    return payload


def _warning_to_dict(warning: Any) -> dict[str, Any]:
    if hasattr(warning, "to_dict"):
        return warning.to_dict()
    if isinstance(warning, Mapping):
        return dict(warning)
    return _warning_dict(str(warning), str(warning))


def _explicit_dk_activation_scale(
    base_key: str,
    rate_key: str,
    exposure_profile: Mapping[str, Any],
) -> tuple[float | None, dict[str, float] | None]:
    for key in (base_key, rate_key):
        value = exposure_profile.get(key)
        if not isinstance(value, Mapping):
            continue
        if "d_factor" not in value or "k_factor" not in value:
            continue
        try:
            d_factor = max(0.0, float(value["d_factor"]))
            k_factor = max(0.0, float(value["k_factor"]))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(d_factor) or not math.isfinite(k_factor):
            continue
        return d_factor * k_factor, {"d_factor": d_factor, "k_factor": k_factor}
    return None, None


def _matrix_gsh_activation_scale(
    base_key: str,
    rate_key: str,
    exposure_profile: Mapping[str, Any],
    inhibition_burdens: Mapping[str, _InhibitionBurdenResolution],
    warnings: list[dict[str, Any]],
) -> tuple[float, str, dict[str, Any]]:
    direct = inhibition_burdens.get(base_key)
    direct_fallback_details: dict[str, Any] = {}
    if direct is not None:
        try:
            direct_activation_ratio = float(direct.activation_burden_ratio)
        except (TypeError, ValueError):
            direct_activation_ratio = math.nan
        direct_fallback_details = {
            "direct_activation_available": False,
            "inhibition_status": direct.status,
            "review_required": direct.review_required,
            "activation_burden_ratio": _round(direct_activation_ratio, 6)
            if math.isfinite(direct_activation_ratio)
            else None,
        }
        has_quantified_direct_activation = (
            direct.status == "mechanism_resolved"
            and not direct.review_required
            and math.isfinite(direct_activation_ratio)
            and direct_activation_ratio >= 0.0
            and not math.isclose(direct_activation_ratio, 1.0)
        )
    else:
        direct_activation_ratio = math.nan
        has_quantified_direct_activation = False

    if direct is not None and has_quantified_direct_activation:
        return (
            max(0.0, direct_activation_ratio),
            "direct_activation_burden_ratio",
            {
                "activation_burden_ratio": _round(direct_activation_ratio, 6),
                "inhibition_status": direct.status,
            },
        )

    dk_scale, dk_details = _explicit_dk_activation_scale(base_key, rate_key, exposure_profile)
    if dk_scale is not None:
        return (
            dk_scale,
            "d_times_k_approximation",
            {
                **direct_fallback_details,
                **(dk_details or {}),
                "internal_d_or_k_computation": False,
            },
        )

    warnings.append(
        _warning_dict(
            "gsh_upstream_activation_missing_neutral",
            "GSH-relevant matrix contribution lacked direct activation burden or explicit D/K factors; using neutral scaling.",
            field="upstream_activation_burden_ratio",
        )
    )
    return (
        1.0,
        "neutral_missing_upstream_activation",
        {**direct_fallback_details, "upstream_activation_burden_ratio": 1.0},
    )


def _compute_matrix_gsh_redox_status(
    exposure_profile: Mapping[str, Any],
    *,
    genotypes: Mapping[str, str],
    tissue: str,
    inhibition_burdens: Mapping[str, _InhibitionBurdenResolution],
) -> GSHStatus:
    params = _get_interaction_params()["gsh_depletion"]
    consumers = params["consumers"]

    baseline_mM = float(params["baseline_gsh_mM"])
    critical_fraction = float(params["critical_threshold_fraction"])
    synthesis_capacity = float(params["synthesis_rate_umol_h_g"])
    half_life_h = float(params.get("half_life_h", 2.5))
    liver_water_fraction = 0.70
    baseline_umol_g = baseline_mM * liver_water_fraction
    turnover_capacity = baseline_umol_g * math.log(2) / half_life_h

    rate_map = _to_gsh_rate_map(dict(exposure_profile))
    if str(genotypes.get("GSTM1", "active")).lower() in {"null", "null/null", "deletion", "0"}:
        if "PAH_umol_h_g" in rate_map:
            rate_map["PAH_umol_h_g"] *= 0.1

    warnings: list[dict[str, Any]] = []
    total_base_consumption = 0.0
    total_scaled_consumption = 0.0
    contributions: dict[str, dict[str, Any]] = {}

    for rate_key, flux_umol_h_g in rate_map.items():
        base_key = rate_key.removesuffix("_umol_h_g")
        consumer_key = _gsh_consumer_key(base_key)
        consumer = consumers.get(consumer_key)
        gsh_ratio = float(consumer.get("gsh_per_umol_substrate", 1.0)) if consumer else 1.0
        substrate_flux = float(flux_umol_h_g)
        base_gsh_drain = substrate_flux * gsh_ratio
        upstream_scale, scaling_source, scaling_details = _matrix_gsh_activation_scale(
            base_key,
            rate_key,
            exposure_profile,
            inhibition_burdens,
            warnings,
        )
        scaled_gsh_drain = base_gsh_drain * upstream_scale
        total_base_consumption += base_gsh_drain
        total_scaled_consumption += scaled_gsh_drain
        contributions[consumer_key] = {
            "substrate_flux_umol_h_g": _round(substrate_flux, 4),
            "base_gsh_consumption_umol_h_g": _round(base_gsh_drain, 4),
            "gsh_consumption_umol_h_g": _round(scaled_gsh_drain, 4),
            "stoichiometry": gsh_ratio,
            "upstream_activation_scale": _round(upstream_scale, 6),
            "upstream_scaling_source": scaling_source,
            "upstream_scaling_provenance": scaling_details,
            "fraction_of_total": 0.0,
        }

    for contribution in contributions.values():
        if total_scaled_consumption > 0:
            contribution["fraction_of_total"] = _round(
                contribution["gsh_consumption_umol_h_g"] / total_scaled_consumption,
                3,
            )

    gsh_result = compute_gsh_redox_capacity(
        GSHRedoxCapacityInput(
            tissue=tissue,
            consumption_load=total_scaled_consumption,
            synthesis_capacity=synthesis_capacity,
            turnover_capacity=turnover_capacity,
            baseline_capacity=1.0,
            model_version=GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY,
            metadata={
                "source": "interaction_engine.compute_interaction_matrix",
                "matrix_level_shared_pool": True,
                "upstream_activation_scaled": True,
                "base_gsh_consumption_load": total_base_consumption,
                "scaled_gsh_consumption_load": total_scaled_consumption,
                "upstream_activation_source_preference": "direct_activation_burden_ratio",
                "d_times_k_fallback_only_when_explicit": True,
            },
        )
    )
    warnings.extend(_warning_to_dict(warning) for warning in gsh_result.warnings)

    fraction_normal = gsh_result.gsh_fraction
    steady_state_mM = baseline_mM * fraction_normal
    net_rate = synthesis_capacity - total_scaled_consumption
    consumption_exceeds_synthesis = bool(total_scaled_consumption >= synthesis_capacity)
    tipping_point_reached = consumption_exceeds_synthesis

    time_to_depletion_h: float | None = None
    if total_scaled_consumption > synthesis_capacity:
        net_drain = total_scaled_consumption - synthesis_capacity
        gsh_to_lose = baseline_umol_g * (1.0 - critical_fraction)
        time_to_depletion_h = _round(gsh_to_lose / net_drain, 2)

    impaired_pathways: list[str] = []
    if fraction_normal < 0.30:
        impaired_pathways.append("Phase II GST conjugation (general)")
    if fraction_normal < critical_fraction:
        impaired_pathways.extend(
            [
                "BPDE-GSH conjugation (PAH detox)",
                "Arsenic methylation (GSTO1/GSTO2)",
                "Chromium(VI) reduction",
                "ROS scavenging (GPx)",
            ]
        )
    if fraction_normal < 0.10:
        impaired_pathways.append("Thioredoxin/GSSG reductase backup overwhelmed")

    return GSHStatus(
        baseline_gsh_mM=baseline_mM,
        steady_state_gsh_mM=_round(steady_state_mM, 3),
        fraction_normal=_round(fraction_normal, 6),
        consumption_rate_umol_h_g=_round(total_scaled_consumption, 4),
        synthesis_rate_umol_h_g=_round(synthesis_capacity, 4),
        net_rate_umol_h_g=_round(net_rate, 4),
        consumption_exceeds_synthesis=consumption_exceeds_synthesis,
        tipping_point_reached=tipping_point_reached,
        tipping_point_multiplier=_round(total_scaled_consumption / synthesis_capacity, 3)
        if synthesis_capacity > 0
        else 0.0,
        impaired_pathways=impaired_pathways,
        individual_contributions=contributions,
        time_to_depletion_h=time_to_depletion_h,
        tissue=tissue,
        model_version=str(gsh_result.model_version),
        redox_capacity_ratio=gsh_result.redox_capacity_ratio,
        detox_penalty_multiplier=gsh_result.detox_penalty_multiplier,
        warnings=warnings,
        metadata={
            "source": "gsh_redox_capacity.compute_gsh_redox_capacity",
            "legacy_matrix_gsh_behavior": False,
            "base_gsh_consumption_load": _round(total_base_consumption, 6),
            "scaled_gsh_consumption_load": _round(total_scaled_consumption, 6),
            "turnover_at_current_fraction": gsh_result.metadata.get("turnover_at_current_fraction")
            if gsh_result.metadata
            else None,
            "redox_capacity_result": gsh_result.to_dict(),
        },
    )


def _severity_sorted(interactions: list[CriticalInteraction]) -> list[CriticalInteraction]:
    """Sort critical interactions by severity then amplification."""
    return sorted(
        interactions,
        key=lambda item: (_SEVERITY_RANK.get(item.severity, 0), item.genotype_amplification),
        reverse=True,
    )


def _make_inert_gsh_status(tissue: str) -> GSHStatus:
    """Return a GSHStatus representing an unperturbed GSH pool."""
    params = _get_interaction_params().get("gsh_depletion", {})
    baseline = float(params.get("baseline_gsh_mM", 7.0))
    synthesis = float(params.get("synthesis_rate_umol_h_g", 1.0))
    return GSHStatus(
        baseline_gsh_mM=baseline,
        steady_state_gsh_mM=baseline,
        fraction_normal=1.0,
        consumption_rate_umol_h_g=0.0,
        synthesis_rate_umol_h_g=synthesis,
        net_rate_umol_h_g=synthesis,
        consumption_exceeds_synthesis=False,
        tipping_point_reached=False,
        tipping_point_multiplier=0.0,
        impaired_pathways=[],
        individual_contributions={},
        time_to_depletion_h=None,
        tissue=tissue,
        model_version=str(GSHModelVersion.PHASE7_QUASI_STEADY_RELATIVE_CAPACITY),
        redox_capacity_ratio=1.0,
        detox_penalty_multiplier=1.0,
        metadata={
            "source": "interaction_engine._make_inert_gsh_status",
            "gsh_depletion_enabled": False,
            "legacy_matrix_gsh_behavior": False,
        },
    )


def _percentile(values: list[float], pct: float) -> float:
    """Return a simple linear-interpolation percentile of ``values``."""
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (pct / 100.0) * (len(ordered) - 1)
    lower = int(math.floor(rank))
    upper = int(math.ceil(rank))
    if lower == upper:
        return ordered[lower]
    fraction = rank - lower
    return ordered[lower] + (ordered[upper] - ordered[lower]) * fraction


# ── Public API ────────────────────────────────────────────────────────────


def enzyme_induction_modifier(
    lifestyle: Mapping[str, bool | int | float],
) -> EnzymeInductionProfile:
    """Return fold-induction multipliers for CYP enzymes based on lifestyle.

    Supported keys include ``smoking``, ``heavy_smoking``, ``pack_years``,
    ``alcohol_moderate``, ``alcohol_heavy``, ``chronic_alcohol``,
    ``dioxin_exposed``, ``TCDD_exposed``, ``obesity``, and ``NAFLD``.
    """
    params = _get_interaction_params()
    induction = params["enzyme_induction"]
    enzyme_folds = {enzyme: 1.0 for enzyme in _collect_known_enzymes(params)}
    active_inducers: list[str] = []

    if lifestyle.get("smoking") or lifestyle.get("heavy_smoking"):
        _apply_enzyme_folds(enzyme_folds, induction["smoking"])
        active_inducers.append("smoking")

    if lifestyle.get("alcohol_moderate"):
        enzyme_folds["CYP2E1"] = max(enzyme_folds.get("CYP2E1", 1.0), 2.0)
        active_inducers.append("alcohol_moderate")

    if lifestyle.get("alcohol_heavy") or lifestyle.get("chronic_alcohol"):
        _apply_enzyme_folds(enzyme_folds, induction["chronic_alcohol"])
        active_inducers.append("chronic_alcohol")

    if (lifestyle.get("heavy_smoking") or lifestyle.get("pack_years", 0) >= 30) and (
        lifestyle.get("alcohol_heavy") or lifestyle.get("chronic_alcohol")
    ):
        _apply_enzyme_folds(enzyme_folds, induction["heavy_smoking_chronic_alcohol_combined"])
        active_inducers.append("heavy_smoking_chronic_alcohol_combined")

    if lifestyle.get("dioxin_exposed") or lifestyle.get("TCDD_exposed"):
        _apply_enzyme_folds(enzyme_folds, induction["TCDD_dioxin"])
        active_inducers.append("TCDD_dioxin")

    if lifestyle.get("obesity") or lifestyle.get("NAFLD") or lifestyle.get("insulin_resistance"):
        _apply_enzyme_folds(enzyme_folds, induction["obesity_insulin_resistance"])
        active_inducers.append("obesity_insulin_resistance")

    return EnzymeInductionProfile(
        enzyme_folds={key: _round(value, 3) for key, value in sorted(enzyme_folds.items())},
        active_inducers=sorted(dict.fromkeys(active_inducers)),
    )


def competitive_inhibition_flux(
    enzyme: str,
    substrates: dict[str, float],
    *,
    genotype_modifiers: dict[str, float] | None = None,
    tissue: str = "Liver",
    param_perturbations: dict[str, dict[str, float]] | None = None,
    inhibition_contexts: Mapping[str, Any] | None = None,
    include_biological_outputs: bool = True,
) -> CompetitiveInhibitionResult:
    """Compute adjusted flux using the centralized reversible-inhibition resolver.

    param_perturbations optionally rescales Km and Vmax for Monte Carlo
    uncertainty propagation. inhibition_contexts is an internal live interaction integration seam for
    typed resolver-compatible contexts and is not exposed through public API
    output.
    """
    params = _get_interaction_params()["competitive_inhibition"]
    if enzyme not in params or enzyme.startswith("_"):
        raise ValueError(f"Unknown enzyme for competitive inhibition: {enzyme}")

    genotype_modifiers = genotype_modifiers or {}
    perturbations = param_perturbations or {}
    enzyme_data = params[enzyme]
    param_substrates = enzyme_data["substrates"]
    genotype_multiplier = float(genotype_modifiers.get(enzyme, 1.0))
    results: dict[str, SubstrateFluxChange] = {}

    def _scale_for(name: str, field_name: str) -> float:
        entry = perturbations.get(name) or {}
        return float(entry.get(field_name, 1.0))

    for sub_name, conc_uM in substrates.items():
        concentration = float(conc_uM)
        if concentration <= 0:
            continue

        sub_params = param_substrates.get(
            sub_name,
            {
                "Km_uM": 50.0,
                "Vmax_relative": 0.5,
                "product": "unknown",
                "product_carcinogenic": False,
            },
        )
        Km_A = float(sub_params["Km_uM"]) * _scale_for(sub_name, "Km")
        Vmax_A = (
            float(sub_params.get("Vmax_relative", 1.0))
            * genotype_multiplier
            * _scale_for(sub_name, "Vmax")
        )
        hill_n = float(sub_params.get("hill_coefficient", 1.0))
        substrate_power = concentration ** hill_n
        Km_power = Km_A ** hill_n
        single_flux = (Vmax_A * substrate_power) / (Km_power + substrate_power)

        resolution = _resolve_live_inhibition_modifier(
            enzyme=enzyme,
            target_substrate=sub_name,
            target_concentration=concentration,
            target_substrate_power=substrate_power,
            target_km_effective=Km_power,
            substrates=substrates,
            param_substrates=param_substrates,
            scale_for=_scale_for,
            inhibition_contexts=inhibition_contexts,
        )
        if resolution["modifier"] is None:
            competitive_flux = single_flux
            kinetic_modifier = None
            modified_flux = None
            modifier_applied_once = False
        else:
            application = apply_kinetic_modifier_once(single_flux, resolution["modifier"])
            competitive_flux = application.modified_flux
            kinetic_modifier = application.kinetic_modifier
            modified_flux = application.modified_flux
            modifier_applied_once = application.applied_once

        flux_change_fraction = (
            (competitive_flux - single_flux) / single_flux if single_flux > 0 else 0.0
        )

        product_carcinogenic = bool(sub_params.get("product_carcinogenic", False))
        flux_change = SubstrateFluxChange(
            single_flux=_round(single_flux, 6),
            competitive_flux=_round(competitive_flux, 6),
            flux_change_fraction=_round(flux_change_fraction, 4),
            inhibition_term=_round(float(resolution.get("inhibition_term", 0.0)), 4),
            activated_product_flux=_round(competitive_flux if product_carcinogenic else 0.0, 6),
            Km_uM=Km_A,
            concentration_uM=concentration,
            product=str(sub_params.get("product", "unknown")),
            product_carcinogenic=product_carcinogenic,
            baseline_flux=_round(single_flux, 6),
            kinetic_modifier=_round(kinetic_modifier, 8) if kinetic_modifier is not None else None,
            modified_flux=_round(modified_flux, 6) if modified_flux is not None else None,
            kinetic_resolution_status=str(resolution["status"]),
            kinetic_warning_codes=tuple(sorted(resolution["warning_codes"])),
            inhibition_mode=str(resolution["mode"]),
            centralized_resolver_used=bool(resolution["centralized_resolver_used"]),
            modifier_applied_once=modifier_applied_once,
            discrepancy_classification=str(resolution["discrepancy_classification"]),
            aggregate_resolution=(
                deepcopy(resolution["aggregate_resolution"])
                if resolution.get("aggregate_resolution") is not None
                else None
            ),
        )
        if include_biological_outputs:
            flux_change.biological_output = _build_live_biological_output(
                enzyme=enzyme,
                substrate=sub_name,
                tissue=tissue,
                flux=flux_change,
                substrate_parameters=sub_params,
            )
        results[sub_name] = flux_change

    return CompetitiveInhibitionResult(enzyme=enzyme, substrates=results)


def _resolve_live_inhibition_modifier(
    *,
    enzyme: str,
    target_substrate: str,
    target_concentration: float,
    target_substrate_power: float,
    target_km_effective: float,
    substrates: Mapping[str, float],
    param_substrates: Mapping[str, Mapping[str, Any]],
    scale_for,
    inhibition_contexts: Mapping[str, Any] | None,
) -> dict[str, Any]:
    explicit_entries = _inhibition_context_entries(inhibition_contexts, target_substrate)
    if explicit_entries is not None:
        return _resolve_explicit_live_inhibition_contexts(
            entries=explicit_entries,
            target_substrate=target_substrate,
            target_concentration=target_concentration,
            target_km_effective=target_km_effective,
        )

    positive_competitors = [
        (name, float(conc))
        for name, conc in substrates.items()
        if name != target_substrate and float(conc) > 0.0
    ]
    if not positive_competitors:
        return _modifier_record(
            modifier=1.0,
            status="mechanism_absent",
            mode=InhibitionMode.COMPETITIVE,
            applied_resolver=False,
            discrepancy_classification="not_applicable_no_competitor",
        )

    normalized_inhibitor_load = 0.0
    warning_codes: set[str] = set()
    resolved_any = False
    competitor_details: list[dict[str, Any]] = []
    for other_name, other_conc in positive_competitors:
        ki = get_ki(enzyme, other_name, target_substrate=target_substrate)
        competitor_warnings = [
            warning.code.upper()
            for warning in (ki.warnings or [])
            if getattr(warning, "code", None)
        ]
        if ki.warnings:
            warning_codes.update(competitor_warnings)
        detail: dict[str, Any] = {
            "inhibitor": other_name,
            "concentration_uM": other_conc,
            "ki_uM": None,
            "normalized_load_contribution": None,
            "resolution_method": str(ki.resolution_method),
            "source_kind": str(ki.source_kind),
            "confidence": ki.uncertainty.confidence if ki.uncertainty else None,
            "warnings": competitor_warnings,
            "resolved": False,
            "metadata": deepcopy(ki.metadata) if ki.metadata is not None else None,
        }
        if ki.value is None:
            warning_codes.add("KI_MISSING")
            if "KI_MISSING" not in detail["warnings"]:
                detail["warnings"].append("KI_MISSING")
            competitor_details.append(detail)
            continue
        ki_value = float(ki.value)
        if ki.resolution_method.value == "assumed_equal_km":
            ki_value *= float(scale_for(other_name, "Km"))
        if ki_value <= 0.0 or not math.isfinite(ki_value):
            warning_codes.add("OUTSIDE_REVERSIBLE_INHIBITION_DOMAIN")
            detail["warnings"].append("OUTSIDE_REVERSIBLE_INHIBITION_DOMAIN")
            competitor_details.append(detail)
            continue
        contribution = other_conc / ki_value
        normalized_inhibitor_load += contribution
        resolved_any = True
        detail.update(
            {
                "ki_uM": _round(ki_value, 6),
                "normalized_load_contribution": _round(contribution, 8),
                "resolved": True,
            }
        )
        competitor_details.append(detail)

    unresolved_count = sum(1 for detail in competitor_details if not detail["resolved"])
    aggregate_resolution = {
        "mode": "implicit_competitive_load",
        "target_substrate": target_substrate,
        "active_competitor_count": len(positive_competitors),
        "resolved_competitor_count": len(positive_competitors) - unresolved_count,
        "unresolved_competitor_count": unresolved_count,
        "all_active_competitors_resolved": unresolved_count == 0,
        "normalized_inhibitor_load": _round(normalized_inhibitor_load, 8),
        "competitors": competitor_details,
    }

    if not resolved_any:
        warning_codes.add("INCOMPLETE_COMPETITOR_KI_RESOLUTION")
        aggregate_resolution["aggregate_status"] = str(InhibitionResolutionStatus.REVIEW_REQUIRED)
        aggregate_resolution["aggregate_warnings"] = sorted(warning_codes or {"KI_MISSING"})
        return _modifier_record(
            modifier=None,
            status=InhibitionResolutionStatus.REVIEW_REQUIRED,
            mode=InhibitionMode.COMPETITIVE,
            warning_codes=warning_codes or {"KI_MISSING"},
            applied_resolver=False,
            inhibition_term=0.0,
            discrepancy_classification="changed_fallback_behavior",
            aggregate_resolution=aggregate_resolution,
        )

    if unresolved_count:
        warning_codes.add("INCOMPLETE_COMPETITOR_KI_RESOLUTION")
        aggregate_resolution["aggregate_status"] = str(InhibitionResolutionStatus.REVIEW_REQUIRED)
        aggregate_resolution["aggregate_warnings"] = sorted(warning_codes)
        return _modifier_record(
            modifier=None,
            status=InhibitionResolutionStatus.REVIEW_REQUIRED,
            mode=InhibitionMode.COMPETITIVE,
            warning_codes=warning_codes,
            applied_resolver=False,
            inhibition_term=normalized_inhibitor_load,
            discrepancy_classification="incomplete_competitor_resolution",
            aggregate_resolution=aggregate_resolution,
        )

    resolved = resolve_reversible_inhibition(
        {
            "mode": InhibitionMode.COMPETITIVE,
            "enzyme": enzyme,
            "inhibitor": "aggregate_competitive_load",
            "target_substrate": target_substrate,
            "km_uM": target_km_effective,
            "ki_free_enzyme_uM": 1.0,
            "inhibitor_concentration_uM": normalized_inhibitor_load,
            "substrate_concentration_uM": target_substrate_power,
            "vmax": 1.0,
            "concentration_basis": ConcentrationBasis.MODEL_DERIVED,
            "parameter_concentration_basis": ConcentrationBasis.MODEL_DERIVED,
            "applicability_domain": ApplicabilityDomain.IN_DOMAIN,
            "metadata": {"live_adapter": "aggregate_competitive_load_v1"},
        }
    )
    warning_codes.update(warning.code for warning in (resolved.warnings or []))
    aggregate_resolution["aggregate_status"] = str(resolved.status)
    aggregate_resolution["aggregate_warnings"] = sorted(warning_codes)
    if resolved.status not in {
        InhibitionResolutionStatus.RESOLVED_DIRECT,
        InhibitionResolutionStatus.RESOLVED_DERIVED,
    } or resolved.kernel_result is None:
        return _modifier_record(
            modifier=None,
            status=resolved.status,
            mode=resolved.mode,
            warning_codes=warning_codes,
            applied_resolver=True,
            inhibition_term=normalized_inhibitor_load,
            discrepancy_classification="actual_defect",
            aggregate_resolution=aggregate_resolution,
        )

    return _modifier_record(
        modifier=float(resolved.kernel_result.flux_modifier),
        status=resolved.status,
        mode=resolved.mode,
        warning_codes=warning_codes,
        applied_resolver=True,
        inhibition_term=normalized_inhibitor_load,
        discrepancy_classification="numerical_precision",
        aggregate_resolution=aggregate_resolution,
    )


def _resolve_explicit_live_inhibition_contexts(
    *,
    entries: list[Any],
    target_substrate: str,
    target_concentration: float,
    target_km_effective: float,
) -> dict[str, Any]:
    warning_codes: set[str] = set()
    statuses: list[str] = []
    modes: list[str] = []

    if len(entries) > 1:
        for entry in entries:
            payload = _context_payload(entry)
            mode = enum_from_value(InhibitionMode, payload.get("mode"), InhibitionMode.UNKNOWN)
            modes.append(str(mode))
        warning_codes.add("MULTIPLE_INHIBITORS_NOT_JOINTLY_RESOLVED")
        return _modifier_record(
            modifier=None,
            status=InhibitionResolutionStatus.REVIEW_REQUIRED,
            mode=modes[0] if len(set(modes)) == 1 else InhibitionMode.UNKNOWN,
            warning_codes=warning_codes,
            applied_resolver=True,
            discrepancy_classification="multiple_inhibitors_not_jointly_resolved",
        )

    for entry in entries:
        payload = _context_payload(entry)
        payload.setdefault("target_substrate", target_substrate)
        payload.setdefault("km_uM", target_km_effective)
        payload.setdefault("substrate_concentration_uM", target_concentration)
        payload.setdefault("vmax", 1.0)
        resolved = resolve_reversible_inhibition(payload)
        statuses.append(str(resolved.status))
        modes.append(str(resolved.mode))
        warning_codes.update(warning.code for warning in (resolved.warnings or []))
        if resolved.status in {
            InhibitionResolutionStatus.RESOLVED_DIRECT,
            InhibitionResolutionStatus.RESOLVED_DERIVED,
        } and resolved.kernel_result is not None:
            return _modifier_record(
                modifier=float(resolved.kernel_result.flux_modifier),
                status=resolved.status,
                mode=resolved.mode,
                warning_codes=warning_codes,
                applied_resolver=True,
                discrepancy_classification="not_applicable_explicit_resolver_context",
            )

    return _modifier_record(
        modifier=None,
        status=statuses[0] if statuses else InhibitionResolutionStatus.REVIEW_REQUIRED,
        mode=modes[0] if modes else InhibitionMode.UNKNOWN,
        warning_codes=warning_codes or {"KI_MISSING"},
        applied_resolver=True,
        discrepancy_classification="prior_unsupported_ic50_or_missing_context_behavior",
    )


def _inhibition_context_entries(
    inhibition_contexts: Mapping[str, Any] | None,
    target_substrate: str,
) -> list[Any] | None:
    if inhibition_contexts is None or target_substrate not in inhibition_contexts:
        return None
    entries = inhibition_contexts[target_substrate]
    if isinstance(entries, list):
        return list(entries)
    return [entries]


def _context_payload(entry: Any) -> dict[str, Any]:
    if hasattr(entry, "to_dict"):
        return dict(entry.to_dict())
    return dict(entry)


def _modifier_record(
    *,
    modifier: float | None,
    status: Any,
    mode: Any,
    warning_codes: set[str] | None = None,
    applied_resolver: bool,
    inhibition_term: float = 0.0,
    discrepancy_classification: str,
    aggregate_resolution: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "modifier": modifier,
        "status": str(status),
        "mode": str(enum_from_value(InhibitionMode, mode, InhibitionMode.UNKNOWN)),
        "warning_codes": set(warning_codes or set()),
        "centralized_resolver_used": applied_resolver,
        "inhibition_term": inhibition_term,
        "discrepancy_classification": discrepancy_classification,
        "aggregate_resolution": aggregate_resolution,
    }

def _build_live_biological_output(
    *,
    enzyme: str,
    substrate: str,
    tissue: str,
    flux: SubstrateFluxChange,
    substrate_parameters: Mapping[str, Any],
    selected_resolution: _InhibitionBurdenResolution | None = None,
) -> dict[str, Any]:
    mechanism_state = (
        selected_resolution.status
        if selected_resolution is not None
        else _kinetic_mechanism_state(flux.kinetic_resolution_status)
    )
    review_required = (
        selected_resolution.review_required
        if selected_resolution is not None
        else mechanism_state == "mechanism_unresolved"
    )
    flux_ratio = _safe_flux_ratio(flux.competitive_flux, flux.single_flux)
    warnings = list(flux.kinetic_warning_codes)
    if selected_resolution is not None:
        warnings = list(dict.fromkeys([*warnings, *selected_resolution.warnings]))
    interpretation_substrate = (
        selected_resolution.interpretation_substrate
        if selected_resolution is not None and selected_resolution.interpretation_substrate
        else substrate
    )
    output_role = (
        "selected_authoritative_inhibition_effect"
        if selected_resolution is not None
        else "per_substrate_diagnostic"
    )

    kinetic_effect = {
        "status": flux.kinetic_resolution_status,
        "mechanism_state": mechanism_state,
        "applicability": _kinetic_applicability(flux, review_required),
        "mode": flux.inhibition_mode,
        "equation_id": "reversible_inhibition.mixed.v1" if flux.centralized_resolver_used else None,
        "baseline_flux": flux.baseline_flux,
        "flux_modifier": flux.kinetic_modifier,
        "modified_flux": flux.modified_flux,
        "warnings": warnings,
        "review_required": review_required,
        "centralized_resolver_used": flux.centralized_resolver_used,
        "modifier_applied_once": flux.modifier_applied_once,
        "discrepancy_classification": flux.discrepancy_classification,
        "provenance": {
            "source": "interaction_engine.competitive_inhibition_flux",
            "enzyme": enzyme,
            "substrate": substrate,
            "interpretation_substrate": interpretation_substrate,
            "output_role": output_role,
        },
    }
    if flux.aggregate_resolution is not None:
        kinetic_effect["provenance"]["aggregate_resolution"] = deepcopy(flux.aggregate_resolution)

    annotation = None
    if selected_resolution is not None:
        reaction_role = selected_resolution.reaction_role_interpretation
    else:
        registry = get_default_reaction_role_registry()
        annotation = registry.lookup(enzyme, interpretation_substrate, tissue=tissue)
        reaction_role = _reaction_role_interpretation(annotation, flux_ratio, review_required)
    endpoint_result = selected_resolution.endpoint_toxic_flux_result if selected_resolution else None
    endpoint_block: dict[str, Any]
    if selected_resolution is not None:
        endpoint_block = selected_resolution.endpoint_toxic_flux or {}
        endpoint_block["live_engine_integration"] = True
        endpoint_block["selected_authoritative_effect"] = True
        endpoint_block["diagnostic_only"] = False
    elif flux_ratio is None or review_required:
        endpoint_block = _unresolved_endpoint_block(
            enzyme,
            interpretation_substrate,
            tissue,
            flux_ratio,
            annotation,
            review_required=True,
        )
    else:
        endpoint_result = interpret_endpoint_toxic_flux(
            EndpointFluxInput(
                enzyme=enzyme,
                substrate=interpretation_substrate,
                flux_ratio=flux_ratio,
                tissue=tissue,
                annotation=annotation,
                metadata={
                    "source": "live_interaction_engine",
                    "kinetic_resolution_status": flux.kinetic_resolution_status,
                },
            )
        )
        endpoint_block = endpoint_result.to_dict()
        endpoint_block["live_engine_integration"] = True
        endpoint_block["selected_authoritative_effect"] = False
        endpoint_block["review_required"] = reaction_role["review_required"]

    gsh_relevance = _resolve_live_gsh_relevance(enzyme, substrate, substrate_parameters)
    gsh_result = None
    gsh_block: dict[str, Any]
    if not gsh_relevance["gsh_relevant"]:
        gsh_block = {
            **gsh_relevance,
            "review_required": False,
            "gsh_fraction": None,
            "redox_capacity_ratio": None,
            "detox_penalty_multiplier": 1.0,
            "live_engine_integration": True,
            "diagnostic_only": True,
            "included_in_authoritative_adjusted_risk": False,
        }
    elif endpoint_result is None or review_required:
        gsh_block = {
            **gsh_relevance,
            "review_required": True,
            "gsh_fraction": None,
            "redox_capacity_ratio": None,
            "detox_penalty_multiplier": 1.0,
            "warnings": ["gsh_not_quantified_without_resolved_endpoint_flux"],
            "live_engine_integration": True,
            "diagnostic_only": True,
            "included_in_authoritative_adjusted_risk": False,
        }
    else:
        base_load = max(0.0, float(flux.single_flux)) * float(gsh_relevance["gsh_per_umol_substrate"])
        scaled_load = base_load * float(endpoint_result.activation_burden_ratio)
        gsh_result = compute_gsh_redox_capacity(
            GSHRedoxCapacityInput(
                tissue=tissue,
                consumption_load=scaled_load,
                metadata={
                    "source": "live_interaction_engine",
                    "gsh_relevance_reason": gsh_relevance["gsh_relevance_reason"],
                    "base_gsh_consumption_load": base_load,
                    "endpoint_activation_burden_ratio": endpoint_result.activation_burden_ratio,
                },
            )
        )
        gsh_block = gsh_result.to_dict()
        gsh_block.update(gsh_relevance)
        gsh_block["base_gsh_consumption_load"] = round(base_load, 6)
        gsh_block["scaled_gsh_consumption_load"] = round(scaled_load, 6)
        gsh_block["review_required"] = False
        gsh_block["live_engine_integration"] = True
        gsh_block["diagnostic_only"] = True
        gsh_block["included_in_authoritative_adjusted_risk"] = False

    effective_result = selected_resolution.effective_burden_result if selected_resolution else None
    if selected_resolution is not None:
        effective_block = selected_resolution.effective_burden or {}
        effective_block["live_engine_integration"] = True
        effective_block["selected_authoritative_effect"] = True
        effective_block["diagnostic_only"] = False
        effective_block["includes_diagnostic_gsh_capacity"] = False
    elif endpoint_result is not None and not review_required:
        effective_result = compute_effective_carcinogenic_burden(
            EffectiveBurdenInput(
                endpoint_toxic_flux_result=endpoint_result,
                gsh_redox_capacity_result=gsh_result,
                gsh_relevant=bool(gsh_relevance["gsh_relevant"]),
                gsh_coupling=(
                    GSHBurdenCouplingInput(
                        gsh_relevant=True,
                        base_gsh_consumption_load=max(0.0, float(flux.single_flux))
                        * float(gsh_relevance["gsh_per_umol_substrate"]),
                        upstream_activation_burden_ratio=endpoint_result.activation_burden_ratio,
                        tissue=tissue,
                        metadata={"source": "live_interaction_engine"},
                    )
                    if gsh_relevance["gsh_relevant"]
                    else None
                ),
                metadata={"source": "live_interaction_engine"},
            )
        )
        effective_block = effective_result.to_dict()
        effective_block["live_engine_integration"] = True
        effective_block["review_required"] = False
        effective_block["selected_authoritative_effect"] = False
        effective_block["diagnostic_only"] = True
    else:
        effective_block = {
            "review_required": True,
            "effective_carcinogenic_burden_ratio": None,
            "warnings": ["effective_burden_not_quantified_without_resolved_endpoint_flux"],
            "live_engine_integration": True,
            "selected_authoritative_effect": selected_resolution is not None,
            "diagnostic_only": selected_resolution is None,
        }

    transparency_inputs: list[Any] = []
    if annotation is not None:
        transparency_inputs.append(annotation)
    elif selected_resolution is not None and selected_resolution.reaction_role_annotation is not None:
        transparency_inputs.append(selected_resolution.reaction_role_annotation)
    if endpoint_result is not None:
        transparency_inputs.append(endpoint_result)
    if gsh_result is not None:
        transparency_inputs.append(gsh_result)
    if effective_result is not None:
        transparency_inputs.append(effective_result)
    transparency_inputs.extend(_kinetic_warning_records(flux))
    transparency = build_transparency_report(
        *transparency_inputs,
        validation_summary={
            "live_biological_output_integration": True,
            "kinetic_mechanism_state": mechanism_state,
            "review_required": review_required,
            "gsh_relevant": gsh_relevance["gsh_relevant"],
        },
    ).to_dict()
    transparency["live_engine_integration"] = True

    return _json_sanitize(
        {
            "kinetic_effect": kinetic_effect,
            "reaction_role_interpretation": reaction_role,
            "endpoint_toxic_flux": endpoint_block,
            "gsh_capacity_effect": gsh_block,
            "effective_burden": effective_block,
            "model_transparency": transparency,
            "selected_authoritative_effect": selected_resolution is not None,
            "diagnostic_role": output_role,
            "interpretation_substrate": interpretation_substrate,
        }
    )


def _competitive_substrate_parameters(enzyme: str, substrate: str) -> Mapping[str, Any]:
    enzyme_data = _get_interaction_params()["competitive_inhibition"].get(enzyme, {})
    param_substrates = enzyme_data.get("substrates", {})
    return param_substrates.get(
        substrate,
        {
            "Km_uM": 50.0,
            "Vmax_relative": 0.5,
            "product": "unknown",
            "product_carcinogenic": False,
        },
    )


def _attach_live_biological_outputs(
    competitive_effects: dict[str, CompetitiveInhibitionResult],
    selected_resolutions: Mapping[tuple[str, str], _InhibitionBurdenResolution],
    *,
    tissue: str,
) -> None:
    for enzyme, enzyme_result in competitive_effects.items():
        for substrate, flux in enzyme_result.substrates.items():
            selected_resolution = selected_resolutions.get((enzyme, substrate))
            flux.biological_output = _build_live_biological_output(
                enzyme=enzyme,
                substrate=substrate,
                tissue=tissue,
                flux=flux,
                substrate_parameters=_competitive_substrate_parameters(enzyme, substrate),
                selected_resolution=selected_resolution,
            )


def _kinetic_mechanism_state(status: str) -> str:
    if status == "mechanism_disabled_for_attribution":
        return "mechanism_disabled_for_attribution"
    if status == "mechanism_absent":
        return "mechanism_absent"
    if status in {
        InhibitionResolutionStatus.RESOLVED_DIRECT.value,
        InhibitionResolutionStatus.RESOLVED_DERIVED.value,
    }:
        return "mechanism_resolved"
    return "mechanism_unresolved"


def _kinetic_applicability(flux: SubstrateFluxChange, review_required: bool) -> str:
    if flux.kinetic_resolution_status == "mechanism_absent":
        return "not_applicable_mechanism_absent"
    if review_required:
        return ApplicabilityDomain.NOT_ASSESSABLE.value
    if "CONCENTRATION_BASIS_MISMATCH" in flux.kinetic_warning_codes:
        return ApplicabilityDomain.OUTSIDE_DOMAIN.value
    return ApplicabilityDomain.IN_DOMAIN.value


def _safe_flux_ratio(numerator: float, denominator: float) -> float | None:
    try:
        ratio = float(numerator) / float(denominator)
    except (TypeError, ValueError, ZeroDivisionError):
        return None
    if not math.isfinite(ratio) or ratio <= 0.0:
        return None
    return round(ratio, 8)


def _reaction_role_interpretation(annotation: Any, flux_ratio: float | None, review_required: bool) -> dict[str, Any]:
    role_review_required = False
    if review_required or flux_ratio is None:
        directional = "withheld_review_required"
        role_review_required = True
    elif annotation.reaction_role is ReactionRole.UNKNOWN:
        directional = "withheld_unknown_role"
        role_review_required = True
    elif abs(flux_ratio - 1.0) <= 1e-12:
        directional = "no_directional_change"
        role_review_required = False
    elif flux_ratio < 1.0:
        if annotation.reaction_role is ReactionRole.BIOACTIVATION:
            directional = "toxic_product_formation_may_decrease"
        elif annotation.reaction_role is ReactionRole.DETOXIFICATION:
            directional = "burden_may_increase"
        elif annotation.reaction_role is ReactionRole.CLEARANCE:
            directional = "internal_exposure_may_increase"
        elif annotation.reaction_role is ReactionRole.PROTECTIVE_REPAIR:
            directional = "vulnerability_may_increase"
        else:
            directional = "withheld_unknown_role"
            role_review_required = True
    else:
        if annotation.reaction_role is ReactionRole.BIOACTIVATION:
            directional = "toxic_product_formation_may_increase"
        elif annotation.reaction_role is ReactionRole.DETOXIFICATION:
            directional = "burden_may_decrease"
        elif annotation.reaction_role is ReactionRole.CLEARANCE:
            directional = "internal_exposure_may_decrease"
        elif annotation.reaction_role is ReactionRole.PROTECTIVE_REPAIR:
            directional = "vulnerability_may_decrease"
        else:
            directional = "withheld_unknown_role"
            role_review_required = True

    return {
        "role": str(annotation.reaction_role),
        "directional_interpretation": directional,
        "risk_direction_if_flux_decreases": str(annotation.risk_direction_if_flux_decreases),
        "annotation_record_id": annotation.record_id,
        "review_required": role_review_required or getattr(annotation, "review_status", SMEReviewStatus.UNKNOWN)
        in {SMEReviewStatus.UNKNOWN, SMEReviewStatus.PENDING_TEAM_AGREEMENT, SMEReviewStatus.CANDIDATE},
        "warnings": [warning.to_dict() for warning in (annotation.warnings or [])],
        "sme_notes": [note.to_dict() for note in (annotation.sme_notes or [])],
    }


def _unresolved_endpoint_block(
    enzyme: str,
    substrate: str,
    tissue: str,
    flux_ratio: float | None,
    annotation: Any,
    *,
    review_required: bool,
) -> dict[str, Any]:
    return {
        "enzyme": enzyme,
        "substrate": substrate,
        "tissue": tissue,
        "flux_ratio": flux_ratio,
        "endpoint_toxic_flux_ratio": None,
        "reaction_role": str(annotation.reaction_role),
        "annotation_record_id": annotation.record_id,
        "review_required": review_required,
        "warnings": ["endpoint_toxic_flux_not_quantified_without_resolved_kinetic_effect"],
        "live_engine_integration": True,
    }


def _resolve_live_gsh_relevance(
    enzyme: str,
    substrate: str,
    substrate_parameters: Mapping[str, Any],
) -> dict[str, Any]:
    product = str(substrate_parameters.get("product", ""))
    notes = str(substrate_parameters.get("notes", ""))
    combined = " ".join([enzyme, substrate, product, notes]).lower()
    consumers = _get_interaction_params().get("gsh_depletion", {}).get("consumers", {})

    ignored_consumer_tokens = {"and", "from", "gsh", "gst", "cyp2e1", "cyp1a1", "cyp1a2"}
    for consumer_name, consumer in consumers.items():
        consumer_enzyme = str(consumer.get("enzyme", "")).lower()
        consumer_class = str(consumer.get("substrate_class", "")).lower()
        consumer_notes = str(consumer.get("notes", "")).lower()
        class_tokens = [
            token
            for token in consumer_class.replace("-", "_").split("_")
            if len(token) > 2 and token not in ignored_consumer_tokens
        ]
        enzyme_match = bool(consumer_enzyme) and (
            enzyme.lower() in consumer_enzyme or consumer_enzyme in enzyme.lower()
        )
        substrate_match = substrate.lower() in consumer_name.lower() or substrate.lower() in consumer_notes
        class_match = any(token in combined for token in class_tokens)
        if substrate_match or (enzyme_match and class_match):
            return {
                "gsh_relevant": True,
                "gsh_relevance_reason": "matched_gsh_consumer_annotation",
                "gsh_annotation_source": consumer_name,
                "gsh_per_umol_substrate": float(consumer.get("gsh_per_umol_substrate", 1.0)),
            }

    if "gsh" in combined or "glutathione" in combined:
        return {
            "gsh_relevant": True,
            "gsh_relevance_reason": "explicit_gsh_text_annotation",
            "gsh_annotation_source": "substrate_parameters",
            "gsh_per_umol_substrate": 1.0,
        }
    if any(token in combined for token in ("epoxide", "quinone", "napqi", "ros", "reactive metabolite")):
        return {
            "gsh_relevant": True,
            "gsh_relevance_reason": "reactive_metabolite_annotation",
            "gsh_annotation_source": "substrate_parameters",
            "gsh_per_umol_substrate": 1.0,
        }
    return {
        "gsh_relevant": False,
        "gsh_relevance_reason": "no_biological_gsh_annotation",
        "gsh_annotation_source": None,
        "gsh_per_umol_substrate": 0.0,
    }


def _kinetic_warning_records(flux: SubstrateFluxChange) -> list[AssumptionWarning]:
    return [
        AssumptionWarning(
            code=code,
            message=f"Live kinetic resolver warning: {code}.",
            field="kinetic_effect",
            review_status=SMEReviewStatus.UNKNOWN,
        )
        for code in flux.kinetic_warning_codes
    ]


def _compute_live_mechanism_attribution(
    exposure_profile: Mapping[str, float | dict[str, Any]],
    *,
    genotypes: dict[str, str],
    tissue: str,
    lifestyle: Mapping[str, bool | int | float],
    param_perturbations: dict[str, dict[str, float]] | None,
    expression_perturbations: dict[str, float] | None,
) -> dict[str, Any]:
    state_values = []
    unresolved_states: list[str] = []
    for state in generate_mechanism_states():
        state_result = compute_interaction_matrix(
            exposure_profile,
            genotypes=genotypes,
            tissue=tissue,
            lifestyle=lifestyle,
            enable_induction=state.induction,
            enable_competition=state.competition,
            enable_gsh_depletion=state.gsh,
            param_perturbations=param_perturbations,
            expression_perturbations=expression_perturbations,
            include_biological_outputs=False,
        )
        state_values.append((state, state_result.total_interaction_risk))
        if state.competition and _result_has_unresolved_inhibition(state_result):
            unresolved_states.append(state.key)

    attribution = compute_mechanism_attribution(
        dict(state_values),
        metadata={
            "state_source": "interaction_engine.compute_interaction_matrix",
            "engine_generated_states": True,
            "inhibition_disabled_state": "mechanism_disabled_for_attribution",
            "disabled_inhibition_modifier": 1.0,
            "unresolved_inhibition_not_treated_as_absent_or_disabled": True,
            "unresolved_inhibition_state_keys": unresolved_states,
        },
    ).to_dict()
    attribution["live_engine_integration"] = True
    attribution["state_calculation_source"] = "interaction_engine.compute_interaction_matrix"
    attribution["mechanism_state_distinctions"] = {
        "mechanism_absent": "no inhibition evidence on a live substrate",
        "mechanism_disabled_for_attribution": "competition toggle off; M_inh fixed at 1.0 for attribution only",
        "mechanism_resolved": "live kinetic modifier resolved and applied once",
        "mechanism_unresolved": "review-required inhibition, not absent or disabled",
    }
    return _json_sanitize(attribution)


def _result_has_unresolved_inhibition(result: InteractionMatrixResult) -> bool:
    for enzyme_result in result.competitive_effects.values():
        for flux in enzyme_result.substrates.values():
            if _kinetic_mechanism_state(flux.kinetic_resolution_status) == "mechanism_unresolved":
                return True
    return False


def _json_sanitize(value: Any) -> Any:
    converted = to_jsonable(value)
    if isinstance(converted, float):
        return converted if math.isfinite(converted) else None
    if isinstance(converted, list):
        return [_json_sanitize(item) for item in converted]
    if isinstance(converted, dict):
        return {str(key): _json_sanitize(item) for key, item in converted.items()}
    return converted


def gsh_depletion_model(
    exposure_profile: Mapping[str, float | dict[str, Any]],
    *,
    tissue: str = "Liver",
) -> GSHStatus:
    """Model glutathione pool dynamics under multi-carcinogen exposure.

    ``exposure_profile`` may use either exposure multipliers
    (for example ``{"PAH": 2.0}``) or direct rates with keys ending
    ``"_umol_h_g"``.
    """
    params = _get_interaction_params()["gsh_depletion"]
    consumers = params["consumers"]

    baseline_mM = float(params["baseline_gsh_mM"])
    critical_fraction = float(params["critical_threshold_fraction"])
    synthesis_rate = float(params["synthesis_rate_umol_h_g"])
    half_life_h = float(params.get("half_life_h", 2.5))
    liver_water_fraction = 0.70

    rate_map = _to_gsh_rate_map(_normalize_exposure_profile(exposure_profile))
    total_consumption = 0.0
    contributions: dict[str, dict[str, float]] = {}

    for rate_key, flux_umol_h_g in rate_map.items():
        base_key = rate_key.removesuffix("_umol_h_g")
        consumer_key = {
            "PAH": "PAH_GSTM1",
            "chromium_VI": "chromium_VI",
            "arsenic": "arsenic_methylation",
            "cadmium": "cadmium",
            "acrolein": "acrolein",
            "ethanol": "ethanol_ROS",
            "BPDE": "BPDE_conjugation",
            "acetaminophen": "acetaminophen_NAPQI",
        }.get(base_key, base_key)

        consumer = consumers.get(consumer_key)
        gsh_ratio = float(consumer.get("gsh_per_umol_substrate", 1.0)) if consumer else 1.0
        substrate_flux = float(flux_umol_h_g)
        gsh_drain = substrate_flux * gsh_ratio
        total_consumption += gsh_drain
        contributions[consumer_key] = {
            "substrate_flux_umol_h_g": _round(substrate_flux, 4),
            "gsh_consumption_umol_h_g": _round(gsh_drain, 4),
            "stoichiometry": gsh_ratio,
            "fraction_of_total": 0.0,
        }

    for contribution in contributions.values():
        if total_consumption > 0:
            contribution["fraction_of_total"] = _round(
                contribution["gsh_consumption_umol_h_g"] / total_consumption,
                3,
            )

    baseline_umol_g = baseline_mM * liver_water_fraction
    natural_turnover = baseline_umol_g * math.log(2) / half_life_h
    net_rate = synthesis_rate - total_consumption

    if total_consumption <= 0:
        steady_state_fraction = 1.0
    elif net_rate > 0:
        steady_state_fraction = max(
            0.0,
            1.0 - (total_consumption / (synthesis_rate + natural_turnover)),
        )
    else:
        steady_state_fraction = 0.0

    steady_state_umol_g = max(0.0, baseline_umol_g * steady_state_fraction)
    steady_state_mM = steady_state_umol_g / liver_water_fraction
    fraction_normal = steady_state_mM / baseline_mM if baseline_mM > 0 else 0.0

    time_to_depletion_h: float | None = None
    if total_consumption > synthesis_rate:
        net_drain = total_consumption - synthesis_rate
        gsh_to_lose = baseline_umol_g * (1.0 - critical_fraction)
        time_to_depletion_h = _round(gsh_to_lose / net_drain, 2)

    tipping_point_reached = bool(total_consumption >= synthesis_rate)

    impaired_pathways: list[str] = []
    if fraction_normal < 0.30:
        impaired_pathways.append("Phase II GST conjugation (general)")
    if fraction_normal < critical_fraction:
        impaired_pathways.extend(
            [
                "BPDE-GSH conjugation (PAH detox)",
                "Arsenic methylation (GSTO1/GSTO2)",
                "Chromium(VI) reduction",
                "ROS scavenging (GPx)",
            ]
        )
    if fraction_normal < 0.10:
        impaired_pathways.append("Thioredoxin/GSSG reductase backup overwhelmed")

    return GSHStatus(
        baseline_gsh_mM=baseline_mM,
        steady_state_gsh_mM=_round(steady_state_mM, 3),
        fraction_normal=_round(fraction_normal, 3),
        consumption_rate_umol_h_g=_round(total_consumption, 4),
        synthesis_rate_umol_h_g=synthesis_rate,
        net_rate_umol_h_g=_round(net_rate, 4),
        consumption_exceeds_synthesis=tipping_point_reached,
        tipping_point_reached=tipping_point_reached,
        tipping_point_multiplier=_round(total_consumption / synthesis_rate, 3) if synthesis_rate > 0 else 0.0,
        impaired_pathways=impaired_pathways,
        individual_contributions=contributions,
        time_to_depletion_h=time_to_depletion_h,
        tissue=tissue,
    )


def gsh_depletion_biology_model(
    exposure_profile: Mapping[str, float | dict[str, Any]],
    *,
    tissue: str = "Liver",
    synthesis_scale: float = 1.0,
) -> GSHStatus:
    """Sigmoidal GSH steady state with feedback synthesis + saturable consumption.

    Sibling to :func:`gsh_depletion_model`. Replaces the linear approximation's
    constant synthesis and zero-order consumption with:

    - Hill product-inhibition on the synthesis rate (gamma-GCS de-inhibition
      as GSH falls): ``V_synth(G) = V_max / (1 + (G/Ki)**n)``.
    - Michaelis-Menten saturable consumption in the GSH pool:
      ``V_cons(G) = total_consumption * G / (K_m + G)``.

    The transcendental steady state ``V_synth(G) = k*G + V_cons(G)`` is solved
    by bisection. ``synthesis_scale`` rescales V_max for reduced-capacity
    variants (e.g. 0.70 for GCLC, 0.50 for GCLM); the linear model has no
    direct equivalent.

    Output reuses :class:`GSHStatus` with biology-specific semantics:

    - ``synthesis_rate_umol_h_g``: V_synth at the steady state. Rises above
      the resting rate under exposure due to feedback de-repression.
    - ``consumption_rate_umol_h_g``: substrate-driven (G-saturated) input
      rate, identical to the linear model's value.
    - ``tipping_point_reached``: ``fraction_normal < critical_threshold_fraction``.
      There is no ``consumption >= synthesis`` cliff; the pool stabilises at
      a non-zero floor under feedback de-repression.
    - ``tipping_point_multiplier``: ``total_consumption / V_max_synth`` -
      how saturated the de-inhibited synthesis ceiling is.
    - ``time_to_depletion_h``: Euler-integrated time for the pool to fall
      below the critical threshold from baseline. ``None`` if the steady
      state is above the threshold.
    """
    params = _get_interaction_params()["gsh_depletion"]
    biology_params = params.get("biology_model")
    if biology_params is None:
        raise KeyError(
            "gsh_depletion.biology_model is missing from interaction_parameters.json"
        )
    consumers = params["consumers"]

    baseline_mM = float(params["baseline_gsh_mM"])
    critical_fraction = float(params["critical_threshold_fraction"])
    half_life_h = float(params.get("half_life_h", 2.5))
    liver_water_fraction = 0.70

    Ki_mM = float(biology_params["Ki_feedback_mM"])
    n_feedback = float(biology_params["n_feedback"])
    Km_mM = float(biology_params["Km_GST_mM"])

    G_baseline = baseline_mM * liver_water_fraction
    Ki_amount = Ki_mM * liver_water_fraction
    Km_amount = Km_mM * liver_water_fraction
    k_turnover = math.log(2) / half_life_h

    feedback_at_baseline = 1.0 + (G_baseline / Ki_amount) ** n_feedback
    v_max_synth = synthesis_scale * k_turnover * G_baseline * feedback_at_baseline

    rate_map = _to_gsh_rate_map(_normalize_exposure_profile(exposure_profile))
    total_consumption = 0.0
    contributions: dict[str, dict[str, float]] = {}

    for rate_key, flux_umol_h_g in rate_map.items():
        base_key = rate_key.removesuffix("_umol_h_g")
        consumer_key = {
            "PAH": "PAH_GSTM1",
            "chromium_VI": "chromium_VI",
            "arsenic": "arsenic_methylation",
            "cadmium": "cadmium",
            "acrolein": "acrolein",
            "ethanol": "ethanol_ROS",
            "BPDE": "BPDE_conjugation",
            "acetaminophen": "acetaminophen_NAPQI",
        }.get(base_key, base_key)

        consumer = consumers.get(consumer_key)
        gsh_ratio = float(consumer.get("gsh_per_umol_substrate", 1.0)) if consumer else 1.0
        substrate_flux = float(flux_umol_h_g)
        gsh_drain = substrate_flux * gsh_ratio
        total_consumption += gsh_drain
        contributions[consumer_key] = {
            "substrate_flux_umol_h_g": _round(substrate_flux, 4),
            "gsh_consumption_umol_h_g": _round(gsh_drain, 4),
            "stoichiometry": gsh_ratio,
            "fraction_of_total": 0.0,
        }

    for contribution in contributions.values():
        if total_consumption > 0:
            contribution["fraction_of_total"] = _round(
                contribution["gsh_consumption_umol_h_g"] / total_consumption,
                3,
            )

    def synth(G: float) -> float:
        return v_max_synth / (1.0 + (G / Ki_amount) ** n_feedback)

    def cons(G: float) -> float:
        if G <= 0 or total_consumption <= 0:
            return 0.0
        return total_consumption * G / (Km_amount + G)

    def f(G: float) -> float:
        return synth(G) - k_turnover * G - cons(G)

    G_lo, G_hi = 1e-12, G_baseline
    if f(G_hi) >= 0:
        G_ss = G_hi
    elif f(G_lo) <= 0:
        G_ss = 0.0
    else:
        for _ in range(200):
            G_mid = 0.5 * (G_lo + G_hi)
            if f(G_mid) > 0:
                G_lo = G_mid
            else:
                G_hi = G_mid
            if (G_hi - G_lo) < 1e-12:
                break
        G_ss = 0.5 * (G_lo + G_hi)

    fraction_normal = G_ss / G_baseline if G_baseline > 0 else 0.0
    steady_state_mM = G_ss / liver_water_fraction
    synth_at_ss = synth(G_ss)
    actual_cons_at_ss = cons(G_ss)
    net_rate = synth_at_ss - k_turnover * G_ss - actual_cons_at_ss

    tipping_point_reached = bool(fraction_normal < critical_fraction)
    tipping_point_multiplier = (
        _round(total_consumption / v_max_synth, 3) if v_max_synth > 0 else 0.0
    )

    impaired_pathways: list[str] = []
    if fraction_normal < 0.30:
        impaired_pathways.append("Phase II GST conjugation (general)")
    if fraction_normal < critical_fraction:
        impaired_pathways.extend(
            [
                "BPDE-GSH conjugation (PAH detox)",
                "Arsenic methylation (GSTO1/GSTO2)",
                "Chromium(VI) reduction",
                "ROS scavenging (GPx)",
            ]
        )
    if fraction_normal < 0.10:
        impaired_pathways.append("Thioredoxin/GSSG reductase backup overwhelmed")

    time_to_depletion_h: float | None = None
    if total_consumption > 0 and fraction_normal < critical_fraction:
        G = G_baseline
        G_critical = critical_fraction * G_baseline
        t = 0.0
        dt = 0.01
        max_t = 200.0
        while t < max_t:
            dG = synth(G) - k_turnover * G - cons(G)
            G_next = G + dG * dt
            if G_next <= G_critical:
                if G > G_next:
                    frac = (G - G_critical) / (G - G_next)
                    time_to_depletion_h = _round(t + frac * dt, 2)
                else:
                    time_to_depletion_h = _round(t, 2)
                break
            if abs(dG) < 1e-9:
                break
            G = G_next
            t += dt

    return GSHStatus(
        baseline_gsh_mM=baseline_mM,
        steady_state_gsh_mM=_round(steady_state_mM, 3),
        fraction_normal=_round(fraction_normal, 3),
        consumption_rate_umol_h_g=_round(total_consumption, 4),
        synthesis_rate_umol_h_g=_round(synth_at_ss, 4),
        net_rate_umol_h_g=_round(net_rate, 4),
        consumption_exceeds_synthesis=bool(total_consumption > v_max_synth),
        tipping_point_reached=tipping_point_reached,
        tipping_point_multiplier=tipping_point_multiplier,
        impaired_pathways=impaired_pathways,
        individual_contributions=contributions,
        time_to_depletion_h=time_to_depletion_h,
        tissue=tissue,
    )


def _compute_gsh_detox_penalty(
    carcinogen: str,
    gsh_fraction: float,
    genotypes: dict[str, str],
) -> float:
    """Return risk multiplier (>1) based on GSH depletion and genotype."""
    return _compute_gsh_detox_components(carcinogen, gsh_fraction, genotypes)[2]


def _compute_gsh_detox_components(
    carcinogen: str,
    gsh_fraction: float,
    genotypes: dict[str, str],
    *,
    redox_detox_penalty: float | None = None,
) -> tuple[float, float, float]:
    if CARCINOGEN_GSH_DETOX.get(carcinogen) is None:
        return 1.0, 1.0, 1.0

    genotype_factor = 1.0
    gstm1 = str(genotypes.get("GSTM1", "active")).lower()
    gstp1 = genotypes.get("GSTP1", "Ile105Ile")

    if carcinogen == "PAH" and gstm1 in {"null", "null/null", "deletion", "0"}:
        genotype_factor = 2.5
    elif carcinogen == "PAH" and gstp1 == "Val105Val":
        genotype_factor = 1.5

    if redox_detox_penalty is not None:
        gsh_penalty = max(1.0, float(redox_detox_penalty))
    else:
        if gsh_fraction >= 0.30:
            gsh_penalty = 1.0
        elif gsh_fraction >= 0.20:
            gsh_penalty = 1.0 + (0.30 - gsh_fraction) / 0.10 * 0.5
        else:
            gsh_penalty = 1.5 + (0.20 - gsh_fraction) / 0.20 * 3.0

    return genotype_factor, _round(gsh_penalty, 3), _round(genotype_factor * gsh_penalty, 3)


def _selected_competitive_effect(
    carcinogen: str,
    competitive_effects: dict[str, CompetitiveInhibitionResult],
) -> tuple[str, str, str, SubstrateFluxChange] | None:
    if carcinogen == "benzene":
        pulmonary_candidates: list[tuple[float, str, str, str, SubstrateFluxChange]] = []
        for pulmonary_enzyme in ("CYP2A13", "CYP2F1"):
            enzyme_result = competitive_effects.get(pulmonary_enzyme)
            if enzyme_result is None:
                continue
            sub_effect = enzyme_result.substrates.get("benzene")
            if sub_effect is not None:
                pulmonary_candidates.append(
                    (
                        sub_effect.flux_change_fraction,
                        pulmonary_enzyme,
                        "benzene",
                        "benzene",
                        sub_effect,
                    )
                )
        if pulmonary_candidates:
            resolved_candidates = [
                candidate
                for candidate in pulmonary_candidates
                if _kinetic_mechanism_state(candidate[4].kinetic_resolution_status)
                == "mechanism_resolved"
            ]
            selectable_candidates = resolved_candidates or pulmonary_candidates
            _, enzyme, flux_substrate, interpretation_substrate, sub_effect = max(
                selectable_candidates,
                key=lambda item: item[0],
            )
            return enzyme, flux_substrate, interpretation_substrate, sub_effect
        cyp2e1 = competitive_effects.get("CYP2E1")
        if cyp2e1 is not None:
            sub_effect = cyp2e1.substrates.get("benzene")
            if sub_effect is not None:
                return "CYP2E1", "benzene", "benzene", sub_effect
        return None

    if carcinogen in {"NDMA", "vinyl_chloride"}:
        cyp2e1 = competitive_effects.get("CYP2E1")
        if cyp2e1 is not None:
            sub_effect = cyp2e1.substrates.get(carcinogen)
            if sub_effect is not None:
                return "CYP2E1", carcinogen, carcinogen, sub_effect
        return None

    if carcinogen == "HCA":
        cyp1a1 = competitive_effects.get("CYP1A1")
        if cyp1a1 is not None:
            sub_effect = cyp1a1.substrates.get("PhIP")
            if sub_effect is not None:
                return "CYP1A1", "PhIP", "HCA", sub_effect
        return None

    return None


def _neutral_inhibition_resolution(
    status: str,
    *,
    warnings: list[str] | None = None,
    review_required: bool = False,
    enzyme: str | None = None,
    flux_substrate: str | None = None,
    interpretation_substrate: str | None = None,
    tissue: str | None = None,
    flux_ratio: float | None = None,
    reaction_role_annotation: Any | None = None,
    reaction_role_block: dict[str, Any] | None = None,
    provenance: dict[str, Any] | None = None,
) -> _InhibitionBurdenResolution:
    return _InhibitionBurdenResolution(
        burden_multiplier=1.0,
        activation_burden_ratio=1.0,
        detox_failure_ratio=1.0,
        endpoint_toxic_flux_ratio=1.0 if flux_ratio is not None else None,
        status=status,
        review_required=review_required,
        enzyme=enzyme,
        flux_substrate=flux_substrate,
        interpretation_substrate=interpretation_substrate,
        tissue=tissue,
        flux_ratio=flux_ratio,
        warnings=warnings or [],
        reaction_role_annotation=reaction_role_annotation,
        reaction_role_block=reaction_role_block,
        provenance=provenance or {},
    )


def _warning_codes_from_records(records: list[Any] | tuple[Any, ...]) -> list[str]:
    codes: list[str] = []
    for record in records:
        if hasattr(record, "code"):
            code = str(record.code)
        elif isinstance(record, Mapping):
            code = str(record.get("code", record))
        else:
            code = str(record)
        if code and code not in codes:
            codes.append(code)
    return codes


def _resolve_endpoint_inhibition_burden(
    carcinogen: str,
    competitive_effects: dict[str, CompetitiveInhibitionResult],
    *,
    tissue: str,
    enable_competition: bool,
) -> _InhibitionBurdenResolution:
    if not enable_competition:
        return _neutral_inhibition_resolution(
            "mechanism_disabled",
            provenance={"selection": "competition_toggle_off"},
        )

    selected = _selected_competitive_effect(carcinogen, competitive_effects)
    if selected is None:
        return _neutral_inhibition_resolution(
            "mechanism_absent",
            provenance={"selection": "no_supported_competitive_effect"},
        )

    enzyme, flux_substrate, interpretation_substrate, flux = selected
    mechanism_state = _kinetic_mechanism_state(flux.kinetic_resolution_status)
    flux_ratio = _safe_flux_ratio(flux.competitive_flux, flux.single_flux)
    warning_codes = _warning_codes_from_records(flux.kinetic_warning_codes)
    kinetic_review_required = mechanism_state not in {
        "mechanism_resolved",
        "mechanism_absent",
    } or flux_ratio is None
    registry = get_default_reaction_role_registry()
    annotation = registry.lookup(enzyme, interpretation_substrate, tissue=tissue)
    reaction_role = _reaction_role_interpretation(
        annotation,
        flux_ratio,
        review_required=kinetic_review_required,
    )
    warning_codes.extend(_warning_codes_from_records(reaction_role.get("warnings", [])))
    warning_codes = list(dict.fromkeys(warning_codes))
    base_provenance = {
        "selection": "existing_carcinogen_enzyme_substrate_precedence",
        "kinetic_resolution_status": flux.kinetic_resolution_status,
        "kinetic_mechanism_state": mechanism_state,
        "modifier_applied_once": flux.modifier_applied_once,
        "centralized_resolver_used": flux.centralized_resolver_used,
        "annotation_record_id": annotation.record_id,
        "reaction_role": str(annotation.reaction_role),
        "risk_direction_if_flux_decreases": str(annotation.risk_direction_if_flux_decreases),
    }

    if mechanism_state == "mechanism_absent":
        return _neutral_inhibition_resolution(
            mechanism_state,
            warnings=warning_codes,
            review_required=bool(reaction_role["review_required"]),
            enzyme=enzyme,
            flux_substrate=flux_substrate,
            interpretation_substrate=interpretation_substrate,
            tissue=tissue,
            flux_ratio=flux_ratio,
            reaction_role_annotation=annotation,
            reaction_role_block=reaction_role,
            provenance=base_provenance,
        )
    if mechanism_state != "mechanism_resolved" or flux_ratio is None:
        return _neutral_inhibition_resolution(
            mechanism_state,
            warnings=warning_codes or ["inhibition_burden_not_quantified_without_resolved_flux"],
            review_required=True,
            enzyme=enzyme,
            flux_substrate=flux_substrate,
            interpretation_substrate=interpretation_substrate,
            tissue=tissue,
            flux_ratio=flux_ratio,
            reaction_role_annotation=annotation,
            reaction_role_block=reaction_role,
            provenance=base_provenance,
        )

    endpoint_result = interpret_endpoint_toxic_flux(
        EndpointFluxInput(
            enzyme=enzyme,
            substrate=interpretation_substrate,
            flux_ratio=flux_ratio,
            tissue=tissue,
            annotation=annotation,
            metadata={
                "source": "interaction_engine",
                "flux_substrate": flux_substrate,
                "kinetic_resolution_status": flux.kinetic_resolution_status,
            },
        )
    )
    effective_result = compute_effective_carcinogenic_burden(
        EffectiveBurdenInput(
            endpoint_toxic_flux_result=endpoint_result,
            susceptibility_modifier=1.0,
            gsh_relevant=False,
            metadata={"source": "interaction_engine"},
        )
    )
    warning_codes.extend(_warning_codes_from_records(endpoint_result.warnings))
    warning_codes.extend(_warning_codes_from_records(effective_result.warnings))
    warning_codes = list(dict.fromkeys(warning_codes))

    return _InhibitionBurdenResolution(
        burden_multiplier=effective_result.effective_carcinogenic_burden_ratio,
        activation_burden_ratio=endpoint_result.activation_burden_ratio,
        detox_failure_ratio=endpoint_result.detox_failure_ratio,
        endpoint_toxic_flux_ratio=endpoint_result.endpoint_toxic_flux_ratio,
        status=mechanism_state,
        review_required=bool(reaction_role["review_required"]),
        enzyme=enzyme,
        flux_substrate=flux_substrate,
        interpretation_substrate=interpretation_substrate,
        tissue=tissue,
        flux_ratio=flux_ratio,
        warnings=warning_codes,
        reaction_role_annotation=annotation,
        reaction_role_block=reaction_role,
        endpoint_toxic_flux_result=endpoint_result,
        effective_burden_result=effective_result,
        provenance=base_provenance,
    )


def compute_interaction_matrix(
    exposure_profile: Mapping[str, float | dict[str, Any]],
    *,
    genotypes: dict[str, str] | None = None,
    tissue: str = "Liver",
    lifestyle: Mapping[str, bool | int | float] | None = None,
    enable_induction: bool = True,
    enable_competition: bool = True,
    enable_gsh_depletion: bool = True,
    param_perturbations: dict[str, dict[str, float]] | None = None,
    expression_perturbations: dict[str, float] | None = None,
    include_biological_outputs: bool = True,
) -> InteractionMatrixResult:
    """Combine induction, competition, and GSH depletion into risk scoring.

    The three mechanism toggles (``enable_induction``, ``enable_competition``,
    ``enable_gsh_depletion``) allow the caller to isolate each contribution for
    synergy decomposition. ``param_perturbations`` (substrate Km/Vmax scales)
    and ``expression_perturbations`` (per-enzyme activity scales) provide
    Monte Carlo uncertainty hooks.
    """
    genotypes = dict(genotypes or {})
    lifestyle = dict(lifestyle or {})
    expression_perturbations = expression_perturbations or {}
    normalized_exposure = _normalize_exposure_profile(exposure_profile)
    params = _get_interaction_params()
    rules = params["interaction_rules"]

    if enable_induction:
        induction_effects = enzyme_induction_modifier(lifestyle)
    else:
        induction_effects = EnzymeInductionProfile(enzyme_folds={}, active_inducers=[])

    genotype_activity = {
        enzyme: _genotype_activity_multiplier(enzyme, genotypes.get(enzyme))
        for enzyme in ("CYP1A2", "CYP2E1", "CYP3A4")
    }
    combined_enzyme_activity = {
        enzyme: (
            induction_effects.enzyme_folds.get(enzyme, 1.0)
            * genotype_activity.get(enzyme, 1.0)
            * float(expression_perturbations.get(enzyme, 1.0))
        )
        for enzyme in set(induction_effects.enzyme_folds) | set(genotype_activity)
    }

    individual_risks: dict[str, float] = {}
    present_carcinogens: list[str] = []
    for carcinogen, value in normalized_exposure.items():
        if carcinogen not in BASELINE_RISK_SCORES:
            continue
        present_carcinogens.append(carcinogen)
        individual_risks[carcinogen] = _round(
            BASELINE_RISK_SCORES[carcinogen] * _extract_exposure_multiplier(value),
            3,
        )

    competitive_effects: dict[str, CompetitiveInhibitionResult] = {}
    if enable_competition:
        for enzyme, substrates in _build_competitive_substrates(
            normalized_exposure, tissue=tissue
        ).items():
            competitive_effects[enzyme] = competitive_inhibition_flux(
                enzyme,
                substrates,
                genotype_modifiers=genotype_activity,
                tissue=tissue,
                param_perturbations=param_perturbations,
                include_biological_outputs=False,
            )

    induction_multipliers: dict[str, float] = {}
    inhibition_burdens: dict[str, _InhibitionBurdenResolution] = {}
    selected_inhibition_resolutions: dict[tuple[str, str], _InhibitionBurdenResolution] = {}
    for carcinogen in present_carcinogens:
        if enable_induction:
            induction_multiplier = 1.0
            for enzyme in CARCINOGEN_ENZYME_MAP.get(carcinogen, []):
                induction_multiplier = max(induction_multiplier, combined_enzyme_activity.get(enzyme, 1.0))
        else:
            induction_multiplier = 1.0
        induction_multipliers[carcinogen] = induction_multiplier

        inhibition_burden = _resolve_endpoint_inhibition_burden(
            carcinogen,
            competitive_effects,
            tissue=tissue,
            enable_competition=enable_competition,
        )
        inhibition_burdens[carcinogen] = inhibition_burden
        if inhibition_burden.enzyme and inhibition_burden.flux_substrate:
            selected_inhibition_resolutions[
                (inhibition_burden.enzyme, inhibition_burden.flux_substrate)
            ] = inhibition_burden

    if enable_gsh_depletion:
        gsh_status = _compute_matrix_gsh_redox_status(
            normalized_exposure,
            genotypes=genotypes,
            tissue=tissue,
            inhibition_burdens=inhibition_burdens,
        )
    else:
        gsh_status = _make_inert_gsh_status(tissue)

    mechanism_resolved_risks: dict[str, MechanismResolvedRisk] = {}
    interaction_adjusted_risks: dict[str, float] = {}
    gsh_warning_codes = _warning_codes_from_records(gsh_status.warnings)
    for carcinogen in present_carcinogens:
        base_risk = individual_risks[carcinogen]
        induction_multiplier = induction_multipliers[carcinogen]
        inhibition_burden = inhibition_burdens[carcinogen]

        if enable_gsh_depletion:
            susceptibility_modifier, gsh_pool_penalty, gsh_penalty = _compute_gsh_detox_components(
                carcinogen,
                gsh_status.fraction_normal,
                genotypes,
                redox_detox_penalty=gsh_status.detox_penalty_multiplier,
            )
        else:
            susceptibility_modifier = 1.0
            gsh_pool_penalty = 1.0
            gsh_penalty = 1.0

        final_multiplier = _round(
            induction_multiplier * inhibition_burden.burden_multiplier * gsh_penalty,
            6,
        )
        adjusted_risk = _round(base_risk * final_multiplier, 3)
        resolved = MechanismResolvedRisk(
            carcinogen=carcinogen,
            baseline_relative_risk=base_risk,
            induction_multiplier=_round(induction_multiplier, 6),
            inhibition_burden_multiplier=_round(inhibition_burden.burden_multiplier, 6),
            activation_burden_ratio=_round(inhibition_burden.activation_burden_ratio, 6),
            detox_failure_ratio=_round(inhibition_burden.detox_failure_ratio, 6),
            matrix_gsh_penalty=_round(gsh_penalty, 6),
            gsh_pool_penalty=_round(gsh_pool_penalty, 6),
            susceptibility_modifier=_round(susceptibility_modifier, 6),
            susceptibility_applied_in=(
                "matrix_gsh_penalty" if susceptibility_modifier != 1.0 else None
            ),
            final_mechanism_multiplier=final_multiplier,
            adjusted_relative_risk=adjusted_risk,
            inhibition_status=inhibition_burden.status,
            review_required=inhibition_burden.review_required,
            warnings=list(dict.fromkeys([*inhibition_burden.warnings, *gsh_warning_codes])),
            provenance={
                "baseline_risk_source": "BASELINE_RISK_SCORES",
                "induction_source": "enzyme_induction_modifier",
                "gsh_source": "gsh_redox_capacity.compute_gsh_redox_capacity",
                "gsh": {
                    "model_version": gsh_status.model_version,
                    "redox_capacity_ratio": gsh_status.redox_capacity_ratio,
                    "detox_penalty_multiplier": gsh_status.detox_penalty_multiplier,
                    "matrix_gsh_penalty_applied_once": True,
                    "diagnostic_gsh_capacity_included": False,
                    "legacy_matrix_gsh_behavior": False,
                    "warnings": gsh_warning_codes,
                },
                **inhibition_burden.to_provenance(),
            },
        )
        mechanism_resolved_risks[carcinogen] = resolved
        interaction_adjusted_risks[carcinogen] = resolved.adjusted_relative_risk

    synergy_matrix: dict[str, float] = {}
    for index, left in enumerate(present_carcinogens):
        for right in present_carcinogens[index + 1 :]:
            independent_total = individual_risks[left] + individual_risks[right]
            adjusted_total = interaction_adjusted_risks[left] + interaction_adjusted_risks[right]
            synergy_matrix[f"{left}_x_{right}"] = _round(
                adjusted_total / independent_total if independent_total > 0 else 1.0,
                3,
            )

    total_independent = _round(sum(individual_risks.values()), 3)
    total_interaction = _round(sum(interaction_adjusted_risks.values()), 3)
    interaction_factor = _round(
        total_interaction / total_independent if total_independent > 0 else 1.0,
        4,
    )

    summary_parts: list[str] = []
    induced = [name for name, fold in induction_effects.enzyme_folds.items() if fold > 1.5]
    if induced:
        summary_parts.append(
            f"Enzyme induction ({', '.join(induced)}) increases activation capacity."
        )
    if gsh_status.tipping_point_reached:
        summary_parts.append(
            f"GSH tipping point reached ({gsh_status.tipping_point_multiplier:.1f}x synthesis load)."
        )
    elif gsh_status.fraction_normal < 0.5:
        summary_parts.append(
            f"GSH depleted to {gsh_status.fraction_normal * 100:.0f}% of normal."
        )
    for carcinogen, resolved in mechanism_resolved_risks.items():
        inhibition = resolved.provenance.get("inhibition", {})
        if inhibition.get("enzyme") != "CYP2E1":
            continue
        burden_change = (resolved.inhibition_burden_multiplier - 1.0) * 100.0
        if abs(burden_change) > 10:
            direction = "increased" if burden_change > 0 else "decreased"
            summary_parts.append(
                f"CYP2E1 competition: {carcinogen} endpoint burden {direction} by {abs(burden_change):.0f}%."
            )

    synergy_threshold = float(rules.get("synergy_threshold", 1.2))
    antagonism_threshold = float(rules.get("antagonism_threshold", 0.8))
    if interaction_factor > synergy_threshold:
        classification = "synergistic"
    elif interaction_factor < antagonism_threshold:
        classification = "antagonistic"
    else:
        classification = "near-additive"
    summary_parts.append(f"Overall interaction factor: {interaction_factor:.2f}x ({classification}).")

    mechanism_attribution = None
    if include_biological_outputs:
        _attach_live_biological_outputs(
            competitive_effects,
            selected_inhibition_resolutions,
            tissue=tissue,
        )
        mechanism_attribution = _compute_live_mechanism_attribution(
            normalized_exposure,
            genotypes=genotypes,
            tissue=tissue,
            lifestyle=lifestyle,
            param_perturbations=param_perturbations,
            expression_perturbations=expression_perturbations,
        )

    return InteractionMatrixResult(
        individual_risks=individual_risks,
        interaction_adjusted_risks=interaction_adjusted_risks,
        synergy_matrix=synergy_matrix,
        gsh_status=gsh_status,
        induction_effects=induction_effects,
        competitive_effects=competitive_effects,
        total_independent_risk=total_independent,
        total_interaction_risk=total_interaction,
        interaction_factor=interaction_factor,
        tissue=tissue,
        genotypes=genotypes,
        lifestyle=lifestyle,
        summary=" ".join(summary_parts),
        mechanism_attribution=mechanism_attribution,
        mechanism_resolved_risks=(
            mechanism_resolved_risks if include_biological_outputs else {}
        ),
    )


def identify_critical_interactions(
    genotypes: dict[str, str],
) -> list[CriticalInteraction]:
    """Identify multi-carcinogen interaction scenarios most dangerous for a genotype."""
    gstm1 = str(genotypes.get("GSTM1", "active")).lower()
    gstp1 = genotypes.get("GSTP1", "Ile105Ile")
    cyp2e1 = genotypes.get("CYP2E1", "NM")
    cyp1a2 = genotypes.get("CYP1A2", "NM")
    nat2 = genotypes.get("NAT2", "intermediate")
    nqo1 = genotypes.get("NQO1", "WT")
    aldh2 = genotypes.get("ALDH2", "*1/*1")

    interactions: list[CriticalInteraction] = []

    if gstm1 in {"null", "null/null", "deletion", "0"}:
        interactions.append(
            CriticalInteraction(
                interaction="GSTM1-null × GSH depletion (Cr(VI), As, acrolein)",
                severity="CRITICAL",
                mechanism=(
                    "Loss of PAH-GSH conjugation removes a key detox pathway while other "
                    "GSH consumers can collapse the remaining shared glutathione pool."
                ),
                affected_carcinogens=["PAH", "chromium_VI", "arsenic", "acrolein", "AFB1"],
                genotype_amplification=2.5,
                clinical_note=(
                    "Elevated smoker and mixed occupational-exposure risk; avoid combined "
                    "PAH and heavy-metal exposure when possible."
                ),
            )
        )
        interactions.append(
            CriticalInteraction(
                interaction="GSTM1-null × Smoking (PAH) × Chromium(VI)",
                severity="CRITICAL",
                mechanism=(
                    "Smoking-derived PAH burden plus chromium-driven GSH depletion can trigger "
                    "cascading failure of GST-dependent detoxification."
                ),
                affected_carcinogens=["PAH", "chromium_VI"],
                genotype_amplification=4.0,
                clinical_note="Highest concern for mixed smoking and industrial inhalation profiles.",
            )
        )

    if aldh2 in {"*1/*2", "heterozygote"}:
        interactions.append(
            CriticalInteraction(
                interaction="ALDH2*2 heterozygote × Alcohol × Acetaldehyde",
                severity="HIGH",
                mechanism=(
                    "Reduced aldehyde clearance increases acetaldehyde persistence during alcohol "
                    "exposure, amplifying mucosal and esophageal carcinogenic stress."
                ),
                affected_carcinogens=["acetaldehyde"],
                genotype_amplification=3.0,
                clinical_note="Particularly relevant for regular drinkers and upper GI cancer risk.",
            )
        )
    elif aldh2 in {"*2/*2", "PM_ALDH2"}:
        interactions.append(
            CriticalInteraction(
                interaction="ALDH2*2 homozygous × Alcohol × Acetaldehyde",
                severity="CRITICAL",
                mechanism=(
                    "Near-complete aldehyde dehydrogenase loss leads to profound acetaldehyde "
                    "accumulation during alcohol exposure."
                ),
                affected_carcinogens=["acetaldehyde"],
                genotype_amplification=5.0,
                clinical_note="Very high esophageal cancer risk if alcohol exposure persists.",
            )
        )

    if cyp2e1 in {"NM", "UM_c1c1"}:
        interactions.append(
            CriticalInteraction(
                interaction="Chronic alcohol (CYP2E1 induction) × Benzene × NDMA",
                severity="HIGH",
                mechanism=(
                    "Chronic ethanol exposure induces CYP2E1, increasing low-concentration "
                    "activation of benzene and NDMA; high acute ethanol can still compete "
                    "for the same active site."
                ),
                affected_carcinogens=["benzene", "NDMA", "vinyl_chloride"],
                genotype_amplification=2.0 if cyp2e1 == "UM_c1c1" else 1.5,
                clinical_note="Relevant for moderate/heavy drinkers with solvent or nitrosamine exposure.",
            )
        )

    if cyp1a2 == "UM":
        interactions.append(
            CriticalInteraction(
                interaction="CYP1A2 ultra-rapid × Smoking × HCA",
                severity="HIGH",
                mechanism=(
                    "High constitutive CYP1A2 activity plus smoking-driven induction can sharply "
                    "increase heterocyclic amine activation."
                ),
                affected_carcinogens=["HCA", "NNK"],
                genotype_amplification=3.0,
                clinical_note="Highest concern in heavy smokers with frequent cooked-meat exposure.",
            )
        )

    if nat2 == "slow":
        interactions.append(
            CriticalInteraction(
                interaction="NAT2 slow × HCA × Smoking",
                severity="HIGH",
                mechanism=(
                    "Reduced acetylation impairs HCA detoxification while smoking-induced CYP1A2 "
                    "can increase upstream activation."
                ),
                affected_carcinogens=["HCA"],
                genotype_amplification=2.0,
                clinical_note="Relevant for colorectal and bladder-risk enrichment.",
            )
        )

    if nqo1 == "homozygous_609TT":
        interactions.append(
            CriticalInteraction(
                interaction="NQO1 null × Benzene × Smoking",
                severity="CRITICAL",
                mechanism=(
                    "Loss of quinone detoxification can amplify benzene-derived bone marrow stress "
                    "and hematologic carcinogenicity."
                ),
                affected_carcinogens=["benzene"],
                genotype_amplification=3.5,
                clinical_note="Strongest concern for mixed benzene and smoking exposure.",
            )
        )

    if gstp1 == "Val105Val":
        interactions.append(
            CriticalInteraction(
                interaction="GSTP1 Val105Val × AFB1 × Alcohol",
                severity="HIGH",
                mechanism=(
                    "Reduced GSTP1-mediated epoxide detoxification can elevate AFB1 DNA-adduct "
                    "burden, especially in alcohol-affected liver tissue."
                ),
                affected_carcinogens=["AFB1"],
                genotype_amplification=1.8,
                clinical_note="Most relevant for chronic liver injury or aflatoxin-endemic settings.",
            )
        )

    return _severity_sorted(interactions)


_SYNERGY_MECHANISM_LABELS = {
    "induction": "enzyme induction",
    "competition": "competitive inhibition",
    "gsh": "GSH depletion",
    "induction+competition": "induction x competition",
    "induction+gsh": "induction x GSH",
    "competition+gsh": "competition x GSH",
    "induction+competition+gsh": "induction x competition x GSH",
}


def _effect_by_mechanism(effects: list[dict[str, Any]]) -> dict[str, float]:
    return {
        str(effect["mechanism"]): float(effect["effect"])
        for effect in effects
    }


def _interaction_by_key(terms: list[dict[str, Any]]) -> dict[str, float]:
    keyed: dict[str, float] = {}
    for term in terms:
        key = "+".join(str(mechanism) for mechanism in term["mechanisms"])
        keyed[key] = float(term["effect"])
    return keyed


def _dominant_synergy_mechanism(
    main_effects: dict[str, float],
    pairwise_interactions: dict[str, float],
    three_way_interaction: float,
) -> str:
    terms = {
        _SYNERGY_MECHANISM_LABELS[key]: value
        for key, value in {
            **main_effects,
            **pairwise_interactions,
            "induction+competition+gsh": three_way_interaction,
        }.items()
    }
    positive_terms = {key: value for key, value in terms.items() if value > 0.0}
    if positive_terms:
        return max(positive_terms, key=positive_terms.get)
    if main_effects.get("competition", 0.0) < 0.0:
        return "competitive antagonism"
    return "near-additive"


def _compute_synergy_state_results(
    exposure_profile: Mapping[str, float | dict[str, Any]],
    *,
    genotypes: dict[str, str] | None,
    tissue: str,
    lifestyle: Mapping[str, bool | int | float] | None,
    param_perturbations: dict[str, dict[str, float]] | None = None,
    expression_perturbations: dict[str, float] | None = None,
) -> dict[str, InteractionMatrixResult]:
    state_results: dict[str, InteractionMatrixResult] = {}
    for state in generate_mechanism_states():
        state_results[state.key] = compute_interaction_matrix(
            exposure_profile,
            genotypes=genotypes,
            tissue=tissue,
            lifestyle=lifestyle,
            enable_induction=state.induction,
            enable_competition=state.competition,
            enable_gsh_depletion=state.gsh,
            param_perturbations=param_perturbations,
            expression_perturbations=expression_perturbations,
            include_biological_outputs=False,
        )
    return state_results


def _build_pair_synergy_decomposition(
    pair: str,
    state_results: dict[str, InteractionMatrixResult],
    *,
    tolerance: float,
) -> SynergyDecomposition:
    state_values = {
        state.key: state_results[state.key].synergy_matrix.get(pair, 1.0)
        for state in generate_mechanism_states()
    }
    attribution = compute_mechanism_attribution(state_values, tolerance=tolerance).to_dict()
    main_effects = _effect_by_mechanism(attribution["shapley_main_effects"])
    singleton_effects = _interaction_by_key(attribution["singleton_effects"])
    pairwise_interactions = _interaction_by_key(attribution["pairwise_interactions"])
    three_way_interaction = float(attribution["three_way_interaction"]["effect"])
    reconstruction_residual = float(attribution["interaction_reconstruction_residual"])
    shapley_residual = float(attribution["shapley_residual"])
    additive_estimate = (
        1.0
        + sum(singleton_effects.values())
        + sum(pairwise_interactions.values())
        + three_way_interaction
    )
    compatibility_fields = {
        "policy": "compatibility_only",
        "basis": "singleton_effects_from_eight_state_decomposition",
        "residual": "numerical_reconstruction_check_only",
        "primary_fields": [
            "main_effects",
            "pairwise_interactions",
            "three_way_interaction",
            "reconstruction_residual",
        ],
    }
    shapley_decomposition = {
        "decomposition_basis": "eight_state_shapley",
        "state_count": len(state_values),
        "state_values": {key: _round(value, 6) for key, value in state_values.items()},
        "main_effects": {key: _round(value, 6) for key, value in main_effects.items()},
        "singleton_effects": {key: _round(value, 6) for key, value in singleton_effects.items()},
        "pairwise_interactions": {
            key: _round(value, 6) for key, value in pairwise_interactions.items()
        },
        "three_way_interaction": _round(three_way_interaction, 6),
        "reconstruction_residual": _round(reconstruction_residual, 12),
        "shapley_residual": _round(shapley_residual, 12),
        "residuals_are_zero_within_tolerance": bool(
            attribution["residuals_are_zero_within_tolerance"]
        ),
        "tolerance": float(attribution["tolerance"]),
        "residual_policy": "numerical_reconstruction_check_only",
        "source": "compute_interaction_matrix",
    }
    return SynergyDecomposition(
        pair=pair,
        composite=_round(float(attribution["full_value"]), 4),
        delta_comp=_round(singleton_effects.get("competition", 0.0), 4),
        delta_gsh=_round(singleton_effects.get("gsh", 0.0), 4),
        delta_ind=_round(singleton_effects.get("induction", 0.0), 4),
        additive_estimate=_round(additive_estimate, 4),
        residual=_round(reconstruction_residual, 12),
        shapley_decomposition=shapley_decomposition,
        main_effects={key: _round(value, 6) for key, value in main_effects.items()},
        pairwise_interactions={
            key: _round(value, 6) for key, value in pairwise_interactions.items()
        },
        three_way_interaction=_round(three_way_interaction, 6),
        reconstruction_residual=_round(reconstruction_residual, 12),
        shapley_residual=_round(shapley_residual, 12),
        residuals_are_zero_within_tolerance=bool(
            attribution["residuals_are_zero_within_tolerance"]
        ),
        state_values={key: _round(value, 6) for key, value in state_values.items()},
        compatibility_fields=compatibility_fields,
        dominant_mechanism=_dominant_synergy_mechanism(
            main_effects,
            pairwise_interactions,
            three_way_interaction,
        ),
    )


def decompose_synergy(
    exposure_profile: dict[str, float | dict[str, Any]],
    *,
    genotypes: dict[str, str] | None = None,
    tissue: str = "Liver",
    lifestyle: dict[str, bool] | None = None,
    param_perturbations: dict[str, dict[str, float]] | None = None,
    expression_perturbations: dict[str, float] | None = None,
    include_biological_outputs: bool = False,
    tolerance: float = 1e-9,
) -> dict[str, SynergyDecomposition]:
    """Decompose pairwise synergy with the eight-state attribution layer."""
    _ = include_biological_outputs
    state_results = _compute_synergy_state_results(
        exposure_profile,
        genotypes=genotypes,
        tissue=tissue,
        lifestyle=lifestyle,
        param_perturbations=param_perturbations,
        expression_perturbations=expression_perturbations,
    )
    full_result = state_results["induction+competition+gsh"]
    return {
        pair: _build_pair_synergy_decomposition(pair, state_results, tolerance=tolerance)
        for pair in full_result.synergy_matrix
    }


def monte_carlo_synergy_ci(
    exposure_profile: dict[str, float | dict[str, Any]],
    *,
    genotypes: dict[str, str] | None = None,
    tissue: str = "Liver",
    lifestyle: dict[str, bool] | None = None,
    n_iterations: int = 200,
    km_sigma: float = 0.5,
    expression_sigma: float = 0.3,
    seed: int | None = None,
) -> dict[str, SynergyConfidenceInterval]:
    """Bootstrap 95% CIs for pairwise synergy factors and singleton effects.

    Km/Vmax parameters are perturbed by multiplicative lognormal noise with
    ``sigma=km_sigma`` (default ≈ ±50% one-sigma spread); enzyme expression
    weights are perturbed with ``sigma=expression_sigma`` (default ≈ ±30%).
    For each iteration the eight-state decomposition is summarized as mean and
    2.5/97.5 percentile bounds per carcinogen pair.
    """
    if n_iterations < 2:
        raise ValueError("n_iterations must be >= 2 for confidence intervals")

    import random

    rng = random.Random(seed)
    params = _get_interaction_params()
    competitive = params.get("competitive_inhibition", {})
    substrate_names: set[str] = set()
    for enzyme, entry in competitive.items():
        if enzyme.startswith("_"):
            continue
        for sub_name in entry.get("substrates", {}):
            substrate_names.add(sub_name)

    enzyme_names: set[str] = set(_DEFAULT_ENZYMES) | {"CYP1A2", "CYP2E1", "CYP3A4"}

    composite_draws: dict[str, list[float]] = {}
    delta_comp_draws: dict[str, list[float]] = {}
    delta_gsh_draws: dict[str, list[float]] = {}
    delta_ind_draws: dict[str, list[float]] = {}

    for _ in range(n_iterations):
        param_perturbations = {
            name: {
                "Km": math.exp(rng.gauss(0.0, km_sigma)),
                "Vmax": math.exp(rng.gauss(0.0, km_sigma)),
            }
            for name in substrate_names
        }
        expression_perturbations = {
            enzyme: math.exp(rng.gauss(0.0, expression_sigma)) for enzyme in enzyme_names
        }

        decomposed = decompose_synergy(
            exposure_profile,
            genotypes=genotypes,
            tissue=tissue,
            lifestyle=lifestyle,
            param_perturbations=param_perturbations,
            expression_perturbations=expression_perturbations,
            include_biological_outputs=False,
        )

        for pair, decomposition in decomposed.items():
            composite_draws.setdefault(pair, []).append(decomposition.composite)
            delta_comp_draws.setdefault(pair, []).append(decomposition.delta_comp)
            delta_gsh_draws.setdefault(pair, []).append(decomposition.delta_gsh)
            delta_ind_draws.setdefault(pair, []).append(decomposition.delta_ind)

    def _summary(values: list[float]) -> tuple[float, tuple[float, float]]:
        mean = sum(values) / len(values) if values else 0.0
        return _round(mean, 4), (_round(_percentile(values, 2.5), 4), _round(_percentile(values, 97.5), 4))

    intervals: dict[str, SynergyConfidenceInterval] = {}
    for pair in composite_draws:
        composite_mean, composite_ci = _summary(composite_draws[pair])
        comp_mean, comp_ci = _summary(delta_comp_draws[pair])
        gsh_mean, gsh_ci = _summary(delta_gsh_draws[pair])
        ind_mean, ind_ci = _summary(delta_ind_draws[pair])
        intervals[pair] = SynergyConfidenceInterval(
            pair=pair,
            n_iterations=n_iterations,
            composite_mean=composite_mean,
            composite_ci95=composite_ci,
            delta_comp_mean=comp_mean,
            delta_comp_ci95=comp_ci,
            delta_gsh_mean=gsh_mean,
            delta_gsh_ci95=gsh_ci,
            delta_ind_mean=ind_mean,
            delta_ind_ci95=ind_ci,
        )
    return intervals


def get_interaction_profiles() -> dict[str, dict[str, Any]]:
    """Return a defensive copy of the predefined interaction exposure profiles."""
    return deepcopy(EXPOSURE_PROFILES)


def _competitive_effects_to_compat_dict(
    competitive_effects: dict[str, CompetitiveInhibitionResult],
) -> dict[str, dict[str, Any]]:
    """Convert typed competitive-effect results into source-style dictionaries."""
    return {
        enzyme: {
            substrate: {
                "single_flux": flux.single_flux,
                "competitive_flux": flux.competitive_flux,
                "flux_change_fraction": flux.flux_change_fraction,
                "inhibition_term": flux.inhibition_term,
                "activated_product_flux": flux.activated_product_flux,
                "Km_uM": flux.Km_uM,
                "concentration_uM": flux.concentration_uM,
                "product": flux.product,
                "product_carcinogenic": flux.product_carcinogenic,
                "biological_output": deepcopy(flux.biological_output) if flux.biological_output is not None else None,
            }
            for substrate, flux in result.substrates.items()
        }
        for enzyme, result in competitive_effects.items()
    }


def _interaction_matrix_to_compat_dict(result: InteractionMatrixResult) -> dict[str, Any]:
    """Convert an interaction result into a source-style JSON-serializable dict."""
    payload = {
        "individual_risks": dict(result.individual_risks),
        "interaction_adjusted_risks": dict(result.interaction_adjusted_risks),
        "synergy_matrix": dict(result.synergy_matrix),
        "gsh_status": {
            "baseline_gsh_mM": result.gsh_status.baseline_gsh_mM,
            "steady_state_gsh_mM": result.gsh_status.steady_state_gsh_mM,
            "fraction_normal": result.gsh_status.fraction_normal,
            "gsh_consumption_rate_umol_h_g": result.gsh_status.consumption_rate_umol_h_g,
            "gsh_synthesis_rate_umol_h_g": result.gsh_status.synthesis_rate_umol_h_g,
            "net_gsh_rate_umol_h_g": result.gsh_status.net_rate_umol_h_g,
            "consumption_exceeds_synthesis": result.gsh_status.consumption_exceeds_synthesis,
            "tipping_point_reached": result.gsh_status.tipping_point_reached,
            "tipping_point_multiplier": result.gsh_status.tipping_point_multiplier,
            "impaired_pathways": list(result.gsh_status.impaired_pathways),
            "individual_contributions": deepcopy(result.gsh_status.individual_contributions),
            "time_to_depletion_h": result.gsh_status.time_to_depletion_h,
            "tissue": result.gsh_status.tissue,
            "model_version": result.gsh_status.model_version,
            "redox_capacity_ratio": result.gsh_status.redox_capacity_ratio,
            "detox_penalty_multiplier": result.gsh_status.detox_penalty_multiplier,
            "warnings": deepcopy(result.gsh_status.warnings),
            "metadata": deepcopy(result.gsh_status.metadata),
        },
        "induction_effects": dict(result.induction_effects.enzyme_folds),
        "competitive_effects": _competitive_effects_to_compat_dict(result.competitive_effects),
        "total_independent_risk": result.total_independent_risk,
        "total_interaction_risk": result.total_interaction_risk,
        "interaction_factor": result.interaction_factor,
        "tissue": result.tissue,
        "genotypes": deepcopy(result.genotypes),
        "lifestyle": deepcopy(result.lifestyle),
        "summary": result.summary,
        "mechanism_attribution": deepcopy(result.mechanism_attribution),
    }
    if result.mechanism_attribution is not None:
        payload["mechanism_resolved_risks"] = {
            carcinogen: resolved.to_dict()
            for carcinogen, resolved in result.mechanism_resolved_risks.items()
        }
    return payload


def _critical_interactions_to_compat_list(
    interactions: list[CriticalInteraction],
) -> list[dict[str, Any]]:
    """Convert typed critical-interaction warnings into source-style dicts."""
    return [asdict(item) for item in interactions]


def run_validation_case_1() -> tuple[InteractionMatrixResult, InteractionMatrixResult]:
    """Case 1: smoking-only vs smoking + heavy drinking."""
    print("\n" + "=" * 70)
    print("VALIDATION CASE 1: Smoking-only vs Smoking + Heavy Drinking")
    print("=" * 70)

    genotypes = {"GSTM1": "active", "CYP2E1": "NM", "CYP1A2": "NM"}
    smoker_profile = EXPOSURE_PROFILES["smoker"]
    smoker_drinker_profile = EXPOSURE_PROFILES["smoker_heavy_drinker"]

    result_smoker = compute_interaction_matrix(
        smoker_profile["exposure"],
        genotypes=genotypes,
        tissue="Liver",
        lifestyle=smoker_profile["lifestyle"],
    )
    result_smoker_drinker = compute_interaction_matrix(
        smoker_drinker_profile["exposure"],
        genotypes=genotypes,
        tissue="Liver",
        lifestyle=smoker_drinker_profile["lifestyle"],
    )

    print("\n--- Smoking Only ---")
    print(f"  CYP induction effects: {result_smoker.induction_effects.enzyme_folds}")
    print(f"  Individual risks: {result_smoker.individual_risks}")
    print(f"  Interaction-adjusted risks: {result_smoker.interaction_adjusted_risks}")
    print(f"  GSH fraction normal: {result_smoker.gsh_status.fraction_normal:.2%}")
    print(f"  Total independent risk: {result_smoker.total_independent_risk:.1f}")
    print(f"  Total interaction risk: {result_smoker.total_interaction_risk:.1f}")
    print(f"  Interaction factor: {result_smoker.interaction_factor:.3f}x")

    print("\n--- Smoking + Heavy Drinking ---")
    print(f"  CYP induction effects: {result_smoker_drinker.induction_effects.enzyme_folds}")
    print(
        "  CYP2E1 induction (from alcohol): "
        f"{result_smoker_drinker.induction_effects.enzyme_folds.get('CYP2E1', 1.0):.1f}x"
    )
    print(
        "  CYP1A2 induction (from smoking): "
        f"{result_smoker_drinker.induction_effects.enzyme_folds.get('CYP1A2', 1.0):.2f}x"
    )
    cyp2e1_effects = result_smoker_drinker.competitive_effects.get("CYP2E1")
    if cyp2e1_effects is not None:
        print("  CYP2E1 competitive effects:")
        for substrate, flux in cyp2e1_effects.substrates.items():
            print(f"    {substrate}: {flux.flux_change_fraction * 100:+.1f}% flux change vs single-substrate")
    print(f"  Individual risks: {result_smoker_drinker.individual_risks}")
    print(f"  Interaction-adjusted risks: {result_smoker_drinker.interaction_adjusted_risks}")
    print(
        "  GSH status: "
        f"{result_smoker_drinker.gsh_status.steady_state_gsh_mM:.2f} mM "
        f"({result_smoker_drinker.gsh_status.fraction_normal:.1%} of normal)"
    )
    if result_smoker_drinker.gsh_status.impaired_pathways:
        print(f"  Impaired pathways: {result_smoker_drinker.gsh_status.impaired_pathways}")
    print(f"  Total independent risk: {result_smoker_drinker.total_independent_risk:.1f}")
    print(f"  Total interaction risk: {result_smoker_drinker.total_interaction_risk:.1f}")
    print(f"  Interaction factor: {result_smoker_drinker.interaction_factor:.3f}x")

    synergy_ratio = (
        result_smoker_drinker.total_interaction_risk / result_smoker.total_interaction_risk
        if result_smoker.total_interaction_risk > 0
        else float("inf")
    )
    print(f"\n  -> Synergy: Smoking+Drinking risk is {synergy_ratio:.2f}x higher than smoking alone")
    print(
        "  -> CYP1A2-driven HCA risk delta: "
        f"{result_smoker_drinker.interaction_adjusted_risks.get('HCA', 0):.1f} vs "
        f"{result_smoker.interaction_adjusted_risks.get('HCA', 0):.1f}"
    )
    print(
        "  -> CYP2E1-driven NDMA risk delta: "
        f"{result_smoker_drinker.interaction_adjusted_risks.get('NDMA', 0):.1f} vs "
        f"{result_smoker.interaction_adjusted_risks.get('NDMA', 0):.1f}"
    )
    print(f"  -> Summary: {result_smoker_drinker.summary}")

    return result_smoker, result_smoker_drinker


def run_validation_case_2() -> InteractionMatrixResult:
    """Case 2: GSH depletion tipping point in the industrial-worker scenario."""
    print("\n" + "=" * 70)
    print("VALIDATION CASE 2: GSH Depletion Tipping Point (Industrial Worker)")
    print("=" * 70)

    print("\n--- GSH consumption vs synthesis rate at escalating Cr(VI) exposure ---")
    print("  Baseline: PAH (0.5 umol/h/g) + Arsenic (0.2 umol/h/g) + Acrolein (0.5 umol/h/g)")
    print(f"  GSH synthesis rate: {_get_interaction_params()['gsh_depletion']['synthesis_rate_umol_h_g']} umol/h/g")
    print()
    print(f"  {'Cr(VI) umol/h/g':<18} {'Total GSH drain':<18} {'GSH fraction':<15} {'Tipping Point':<15} {'Impaired?'}")
    print(f"  {'-' * 75}")

    tipping_point_announced = False
    for chromium_level in [0, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0]:
        status = gsh_depletion_model(
            {
                "PAH_umol_h_g": 0.5,
                "arsenic_umol_h_g": 0.2,
                "acrolein_umol_h_g": 0.5,
                "chromium_VI_umol_h_g": chromium_level,
            }
        )
        tipping_label = ""
        if status.tipping_point_reached and not tipping_point_announced:
            tipping_label = "*** TIPPING POINT ***"
            tipping_point_announced = True
        impaired = "YES" if status.impaired_pathways else "no"
        print(
            f"  {chromium_level:<18.1f} {status.consumption_rate_umol_h_g:<18.2f} "
            f"{status.fraction_normal:<15.1%} {str(status.tipping_point_reached):<15} {tipping_label}{impaired}"
        )

    print("\n  --- GSTM1-null comparison at Cr(VI) = 5.0 umol/h/g ---")
    for label, genotype in [("GSTM1-active", {"GSTM1": "active"}), ("GSTM1-null", {"GSTM1": "null"})]:
        rate_map = {
            "PAH_umol_h_g": 0.5 * (0.1 if genotype.get("GSTM1") == "null" else 1.0),
            "arsenic_umol_h_g": 0.2,
            "acrolein_umol_h_g": 0.5,
            "chromium_VI_umol_h_g": 5.0,
        }
        status = gsh_depletion_model(rate_map)
        print(
            f"  {label}: GSH = {status.steady_state_gsh_mM:.2f} mM ({status.fraction_normal:.1%} normal), "
            f"Tipping: {status.tipping_point_reached}, Time to critical depletion: {status.time_to_depletion_h} h"
        )
        if status.impaired_pathways:
            print(f"    Impaired pathways: {status.impaired_pathways}")

    industrial_profile = EXPOSURE_PROFILES["smoker_industrial_worker"]
    result = compute_interaction_matrix(
        industrial_profile["exposure"],
        genotypes={"GSTM1": "null", "CYP2E1": "NM"},
        tissue="Liver",
        lifestyle=industrial_profile["lifestyle"],
    )
    synth_rate = result.gsh_status.synthesis_rate_umol_h_g
    print("\n  Smoker + Industrial Worker (GSTM1-null) full profile:")
    print(
        f"    GSH: {result.gsh_status.steady_state_gsh_mM:.2f} mM "
        f"({result.gsh_status.fraction_normal:.1%} normal)"
    )
    print(
        f"    Consumption: {result.gsh_status.consumption_rate_umol_h_g:.2f} umol/h/g "
        f"vs synthesis {synth_rate} umol/h/g"
    )
    print(
        f"    Tipping point: {result.gsh_status.tipping_point_reached} "
        f"({result.gsh_status.tipping_point_multiplier:.2f}x synthesis rate)"
    )
    if result.gsh_status.time_to_depletion_h is not None:
        print(f"    Time to critical depletion: {result.gsh_status.time_to_depletion_h:.1f} hours")
    print(f"    Impaired: {result.gsh_status.impaired_pathways}")
    print(f"    Interaction factor: {result.interaction_factor:.2f}x")

    return result


def run_validation_case_3() -> tuple[CompetitiveInhibitionResult, CompetitiveInhibitionResult]:
    """Case 3: competitive inhibition and the ethanol paradox at CYP2E1."""
    print("\n" + "=" * 70)
    print("VALIDATION CASE 3: CYP2E1 Competitive Inhibition (Ethanol Paradox)")
    print("=" * 70)
    print("  Reference: PMID:11509752 — benzene cancer risk -62% in DBTEX mixture")

    result_no_ethanol = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "NDMA": 0.5},
        tissue="Liver",
    )
    result_low_ethanol = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "NDMA": 0.5, "ethanol": 200.0},
        tissue="Liver",
    )
    result_high_ethanol = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "NDMA": 0.5, "ethanol": 2000.0},
        tissue="Liver",
    )
    result_dbtex = competitive_inhibition_flux(
        "CYP2E1",
        {"benzene": 10.0, "NDMA": 0.5, "ethanol": 13000.0, "vinyl_chloride": 5.0},
        tissue="Liver",
    )

    print("\n  --- CYP2E1 Flux (benzene activation): No competition vs Ethanol present ---")
    print(f"  {'Scenario':<30} {'Benzene flux':<15} {'Change vs alone':<18} {'NDMA flux':<13} {'Inhibition term'}")
    print(f"  {'-' * 80}")
    for label, result in [
        ("Benzene+NDMA alone", result_no_ethanol),
        ("Low ethanol (200 uM)", result_low_ethanol),
        ("High ethanol (2000 uM)", result_high_ethanol),
        ("DBTEX (Km ethanol, 13000 uM)", result_dbtex),
    ]:
        benzene = result.substrates["benzene"]
        ndma = result.substrates.get("NDMA")
        print(
            f"  {label:<30} {benzene.competitive_flux:<15.5f} "
            f"{benzene.flux_change_fraction * 100:>+.1f}%{'':<11} "
            f"{(ndma.competitive_flux if ndma else 0):<13.5f} {benzene.inhibition_term:.3f}"
        )

    single_flux = result_no_ethanol.substrates["benzene"].single_flux
    dbtex_flux = result_dbtex.substrates["benzene"].competitive_flux
    reduction_pct = (1 - dbtex_flux / single_flux) * 100 if single_flux > 0 else 0.0
    print(f"\n  *** Benzene activation reduction at saturating ethanol: {reduction_pct:.0f}% ***")
    print(f"      (Haddad et al. 2001 reported 62% in DBTEX; model predicts {reduction_pct:.0f}%)")

    print("\n  --- Ethanol CYP2E1 kinetics detail ---")
    print("  Key insight: Ethanol (Km=13,000 uM) has lowest priority for CYP2E1 at physiological concentrations.")
    print("  NDMA (Km=15 uM) has highest affinity — it preferentially monopolizes CYP2E1.")
    print("  At high alcohol concentrations, competitive mass-action overwhelms CYP2E1.")
    print("  Paradox: Drinking reduces benzene activation acutely, but chronic drinking induces CYP2E1.")

    ethanol = result_high_ethanol.substrates.get("ethanol")
    if ethanol is not None:
        print(f"\n  Ethanol flux at 2000 uM: {ethanol.competitive_flux:.5f} (relative units)")
        print("  Acetaldehyde production is proportional to ethanol flux and peaks when ethanol dominates CYP2E1.")

    return result_no_ethanol, result_dbtex


def run_validation_case_4(patient_id: str = "JHBUI-10030") -> InteractionMatrixResult:
    """Case 4: full patient interaction profile for JHBUI-10030."""
    print("\n" + "=" * 70)
    print(f"VALIDATION CASE 4: Patient {patient_id} — Full Interaction Profile")
    print("=" * 70)

    profile_data = EXPOSURE_PROFILES.get("JHBUI_10030", EXPOSURE_PROFILES["smoker_moderate_drinker"])
    exposure = profile_data["exposure"]
    lifestyle = profile_data.get("lifestyle", {})
    genotypes = profile_data.get(
        "genotypes",
        {"GSTM1": "null", "CYP2E1": "NM", "CYP1A2": "NM", "NAT2": "intermediate"},
    )

    print(f"\n  Patient: {patient_id}")
    print("  Profile: Smoker (20 pack-years) + Moderate drinker")
    print(f"  Genotypes: {genotypes}")
    print(f"  Lifestyle: {lifestyle}")
    print("  Tissue: Liver (primary)")

    result = compute_interaction_matrix(
        exposure,
        genotypes=genotypes,
        tissue="Liver",
        lifestyle=lifestyle,
    )

    print("\n  --- Enzyme Induction Effects ---")
    for enzyme, fold in sorted(result.induction_effects.enzyme_folds.items()):
        if fold != 1.0:
            print(f"    {enzyme}: {fold:.2f}x")

    print("\n  --- Individual vs Interaction-Adjusted Risks ---")
    print(f"  {'Carcinogen':<18} {'Indiv. Risk':<14} {'Adj. Risk':<14} {'Change':<10} {'Driver'}")
    print(f"  {'-' * 70}")
    for carcinogen in sorted(result.individual_risks):
        individual = result.individual_risks[carcinogen]
        adjusted = result.interaction_adjusted_risks.get(carcinogen, individual)
        percent_change = (adjusted - individual) / individual * 100 if individual > 0 else 0.0
        enzymes = CARCINOGEN_ENZYME_MAP.get(carcinogen, [])
        driver = "+".join(
            f"{enzyme}({result.induction_effects.enzyme_folds.get(enzyme, 1.0):.1f}x)"
            for enzyme in enzymes
            if result.induction_effects.enzyme_folds.get(enzyme, 1.0) > 1.0
        ) or "GSH/direct"
        print(f"  {carcinogen:<18} {individual:<14.2f} {adjusted:<14.2f} {percent_change:>+.1f}%{'':<3} {driver}")

    print("\n  --- GSH Status ---")
    print(
        f"    Steady-state GSH: {result.gsh_status.steady_state_gsh_mM:.2f} mM "
        f"({result.gsh_status.fraction_normal:.1%} of normal)"
    )
    print(f"    GSH consumption: {result.gsh_status.consumption_rate_umol_h_g:.2f} umol/h/g")
    print(f"    GSH synthesis:   {result.gsh_status.synthesis_rate_umol_h_g:.2f} umol/h/g")
    print(
        f"    Tipping point:   {result.gsh_status.tipping_point_reached} "
        f"({result.gsh_status.tipping_point_multiplier:.2f}x synthesis rate)"
    )
    if result.gsh_status.impaired_pathways:
        print(f"    Impaired paths:  {result.gsh_status.impaired_pathways}")

    cyp2e1_effects = result.competitive_effects.get("CYP2E1")
    if cyp2e1_effects is not None:
        print("\n  --- CYP2E1 Competitive Effects ---")
        for substrate, flux in cyp2e1_effects.substrates.items():
            print(f"    {substrate}: {flux.flux_change_fraction * 100:+.1f}% flux change")

    print("\n  --- Summary ---")
    print(f"    Total independent risk:     {result.total_independent_risk:.2f}")
    print(f"    Total interaction-adj risk: {result.total_interaction_risk:.2f}")
    print(f"    Interaction factor:         {result.interaction_factor:.4f}x")
    classification = (
        "SYNERGISTIC"
        if result.interaction_factor > 1.2
        else "ANTAGONISTIC" if result.interaction_factor < 0.8 else "NEAR-ADDITIVE"
    )
    print(f"    Classification:             {classification}")
    print("\n  --- Critical Interactions for Genotype ---")
    for item in identify_critical_interactions(genotypes)[:4]:
        print(f"    [{item.severity}] {item.interaction}")
        print(f"      Amplification: {item.genotype_amplification:.1f}x")
        print(f"      Note: {item.clinical_note}")
    print(f"\n  Narrative: {result.summary}")

    return result


def run_interaction_validation_cases() -> None:
    """Run the full packaged interaction validation suite."""
    print("\n" + "=" * 70)
    print("ExposoGraph Multi-Carcinogen Interaction Engine — Validation Suite")
    print("=" * 70)
    run_validation_case_1()
    run_validation_case_2()
    run_validation_case_3()
    run_validation_case_4()
    print("\n" + "=" * 70)
    print("All validation cases complete.")
    print("=" * 70)


def run_validation_cases() -> None:
    """Compatibility alias for the interaction validation suite."""
    run_interaction_validation_cases()


def cli_main(argv: list[str] | None = None) -> int:
    """CLI entrypoint compatible with the original standalone interaction module."""
    parser = argparse.ArgumentParser(
        description="ExposoGraph multi-carcinogen interaction engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m ExposoGraph.interaction_cli --validate
  python -m ExposoGraph.interaction_cli --profile smoker_heavy_drinker --genotypes '{"GSTM1":"null","CYP2E1":"NM"}' --tissue Liver
  python -m ExposoGraph.interaction_cli --profile JHBUI_10030
  python -m ExposoGraph.interaction_cli --critical-interactions --genotypes '{"GSTM1":"null","NQO1":"homozygous_609TT"}'
        """,
    )
    parser.add_argument("--profile", type=str, default=None, help="Predefined exposure profile name")
    parser.add_argument("--genotypes", type=str, default="{}", help='JSON dict of genotypes, e.g. \'{"GSTM1":"null","CYP2E1":"NM"}\'')
    parser.add_argument("--tissue", type=str, default="Liver", help="Target tissue (default: Liver)")
    parser.add_argument("--validate", action="store_true", help="Run all 4 validation cases and exit")
    parser.add_argument("--list-profiles", action="store_true", help="List all available profiles")
    parser.add_argument("--critical-interactions", action="store_true", help="Identify critical interactions for given genotype")
    parser.add_argument("--output-json", type=str, default=None, help="Save results to JSON file")

    args = parser.parse_args(argv)

    if args.list_profiles:
        print("\nAvailable exposure profiles:")
        for name, profile in EXPOSURE_PROFILES.items():
            print(f"  {name:<30} {profile['_description'][:70]}...")
        return 0

    if args.validate:
        run_interaction_validation_cases()
        return 0

    try:
        genotypes = json.loads(args.genotypes)
    except json.JSONDecodeError as exc:
        print(f"Error: Could not parse genotypes JSON: {exc}", file=sys.stderr)
        return 1

    if args.critical_interactions:
        interactions = identify_critical_interactions(genotypes)
        if args.output_json:
            with open(args.output_json, "w") as handle:
                json.dump(_critical_interactions_to_compat_list(interactions), handle, indent=2)
            print(f"\nResults saved to: {args.output_json}")
            return 0
        print(f"\nCritical interactions for genotype: {genotypes}")
        for index, item in enumerate(interactions, 1):
            print(f"\n{index}. [{item.severity}] {item.interaction}")
            print(f"   Mechanism: {item.mechanism}")
            print(f"   Amplification: {item.genotype_amplification:.1f}x")
            print(f"   Clinical note: {item.clinical_note}")
        return 0

    if args.profile:
        if args.profile not in EXPOSURE_PROFILES:
            print(f"Error: Unknown profile '{args.profile}'. Use --list-profiles to see options.", file=sys.stderr)
            return 1

        profile = EXPOSURE_PROFILES[args.profile]
        exposure = deepcopy(profile["exposure"])
        lifestyle = deepcopy(profile.get("lifestyle", {}))
        if args.genotypes == "{}" and "genotypes" in profile:
            genotypes = deepcopy(profile["genotypes"])

        print(f"\nRunning profile: {args.profile}")
        print(f"Description: {profile['_description']}")
        result = compute_interaction_matrix(
            exposure,
            genotypes=genotypes,
            tissue=args.tissue,
            lifestyle=lifestyle,
        )

        if args.output_json:
            with open(args.output_json, "w") as handle:
                json.dump(_interaction_matrix_to_compat_dict(result), handle, indent=2)
            print(f"\nResults saved to: {args.output_json}")
            return 0

        print(f"\n--- Results for: {args.profile} ---")
        print(f"Tissue: {args.tissue} | Genotypes: {genotypes}")
        print(f"\nEnzyme Induction: {result.induction_effects.enzyme_folds}")
        print("\nIndividual risks:")
        for carcinogen, risk in sorted(result.individual_risks.items()):
            adjusted = result.interaction_adjusted_risks.get(carcinogen, risk)
            percent_change = (adjusted - risk) / risk * 100 if risk > 0 else 0.0
            print(f"  {carcinogen:<20} {risk:.2f} -> {adjusted:.2f} ({percent_change:+.1f}%)")
        print(
            f"\nGSH Status: {result.gsh_status.steady_state_gsh_mM:.2f} mM "
            f"({result.gsh_status.fraction_normal:.1%} normal)"
        )
        if result.gsh_status.impaired_pathways:
            print(f"Impaired pathways: {result.gsh_status.impaired_pathways}")
        print(f"\nTotal independent risk:    {result.total_independent_risk:.2f}")
        print(f"Total interaction risk:    {result.total_interaction_risk:.2f}")
        print(f"Interaction factor:        {result.interaction_factor:.4f}x")
        print(f"\nSummary: {result.summary}")
        return 0

    parser.print_help()
    return 0


def main(argv: list[str] | None = None) -> int:
    """Package CLI wrapper function."""
    return cli_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
