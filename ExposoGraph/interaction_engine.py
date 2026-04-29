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


@dataclass
class SynergyDecomposition:
    """Additive decomposition of a carcinogen-pair synergy factor.
    ``S_composite ≈ 1 + ΔS_comp + ΔS_gsh + ΔS_ind``.
    ``residual`` captures the departure from pure additivity (cross-mechanism
    coupling from multiplicative combination in the full model).
    """

    pair: str
    composite: float
    delta_comp: float
    delta_gsh: float
    delta_ind: float
    additive_estimate: float
    residual: float


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
) -> CompetitiveInhibitionResult:
    """Compute adjusted flux for substrates competing for the same CYP enzyme.

    ``param_perturbations`` is an optional mapping
    ``{substrate_name: {"Km": scale, "Vmax": scale}}`` that rescales Km / Vmax
    for Monte Carlo uncertainty propagation.
    """
    del tissue  # Reserved for future tissue-specific competition adjustments.

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

        inhibition_term = 0.0
        for other_name, other_conc in substrates.items():
            if other_name == sub_name or float(other_conc) <= 0:
                continue
            other_params = param_substrates.get(other_name, {"Km_uM": 50.0})
            other_Km = float(other_params["Km_uM"]) * _scale_for(other_name, "Km")
            inhibition_term += float(other_conc) / other_Km

        substrate_power = concentration ** hill_n
        Km_power = Km_A ** hill_n
        single_flux = (Vmax_A * substrate_power) / (Km_power + substrate_power)
        competitive_flux = (Vmax_A * substrate_power) / (Km_power * (1.0 + inhibition_term) + substrate_power)
        flux_change_fraction = (
            (competitive_flux - single_flux) / single_flux if single_flux > 0 else 0.0
        )

        product_carcinogenic = bool(sub_params.get("product_carcinogenic", False))
        results[sub_name] = SubstrateFluxChange(
            single_flux=_round(single_flux, 6),
            competitive_flux=_round(competitive_flux, 6),
            flux_change_fraction=_round(flux_change_fraction, 4),
            inhibition_term=_round(inhibition_term, 4),
            activated_product_flux=_round(competitive_flux if product_carcinogenic else 0.0, 6),
            Km_uM=Km_A,
            concentration_uM=concentration,
            product=str(sub_params.get("product", "unknown")),
            product_carcinogenic=product_carcinogenic,
        )

    return CompetitiveInhibitionResult(enzyme=enzyme, substrates=results)


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


def _compute_gsh_detox_penalty(
    carcinogen: str,
    gsh_fraction: float,
    genotypes: dict[str, str],
) -> float:
    """Return risk multiplier (>1) based on GSH depletion and genotype."""
    if CARCINOGEN_GSH_DETOX.get(carcinogen) is None:
        return 1.0

    genotype_factor = 1.0
    gstm1 = str(genotypes.get("GSTM1", "active")).lower()
    gstp1 = genotypes.get("GSTP1", "Ile105Ile")

    if carcinogen == "PAH" and gstm1 in {"null", "null/null", "deletion", "0"}:
        genotype_factor = 2.5
    elif carcinogen == "PAH" and gstp1 == "Val105Val":
        genotype_factor = 1.5

    if gsh_fraction >= 0.30:
        gsh_penalty = 1.0
    elif gsh_fraction >= 0.20:
        gsh_penalty = 1.0 + (0.30 - gsh_fraction) / 0.10 * 0.5
    else:
        gsh_penalty = 1.5 + (0.20 - gsh_fraction) / 0.20 * 3.0

    return _round(genotype_factor * gsh_penalty, 3)


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
            )

    if enable_gsh_depletion:
        gsh_exposure = _to_gsh_rate_map(normalized_exposure)
        if str(genotypes.get("GSTM1", "active")).lower() in {"null", "null/null", "deletion", "0"}:
            if "PAH_umol_h_g" in gsh_exposure:
                # Source-parity logic: absent GSTM1 greatly reduces PAH-GSH conjugation,
                # so PAH contributes far less to shared-pool GSH consumption.
                gsh_exposure["PAH_umol_h_g"] *= 0.1
        gsh_status = gsh_depletion_model(gsh_exposure, tissue=tissue)
    else:
        gsh_status = _make_inert_gsh_status(tissue)

    interaction_adjusted_risks: dict[str, float] = {}
    for carcinogen in present_carcinogens:
        base_risk = individual_risks[carcinogen]

        if enable_induction:
            induction_multiplier = 1.0
            for enzyme in CARCINOGEN_ENZYME_MAP.get(carcinogen, []):
                induction_multiplier = max(induction_multiplier, combined_enzyme_activity.get(enzyme, 1.0))
        else:
            induction_multiplier = 1.0

        competition_multiplier = 1.0
        if enable_competition:
            if carcinogen == "benzene":
                pulmonary_fractions: list[float] = []
                for pulmonary_enzyme in ("CYP2A13", "CYP2F1"):
                    if pulmonary_enzyme in competitive_effects:
                        sub_effect = competitive_effects[pulmonary_enzyme].substrates.get("benzene")
                        if sub_effect is not None:
                            pulmonary_fractions.append(sub_effect.flux_change_fraction)
                if pulmonary_fractions:
                    competition_multiplier = 1.0 + max(pulmonary_fractions)
                elif "CYP2E1" in competitive_effects:
                    sub_effect = competitive_effects["CYP2E1"].substrates.get("benzene")
                    if sub_effect is not None:
                        competition_multiplier = 1.0 + sub_effect.flux_change_fraction
            elif carcinogen in {"NDMA", "vinyl_chloride"} and "CYP2E1" in competitive_effects:
                sub_effect = competitive_effects["CYP2E1"].substrates.get(carcinogen)
                if sub_effect is not None:
                    competition_multiplier = 1.0 + sub_effect.flux_change_fraction
            elif carcinogen == "HCA" and "CYP1A1" in competitive_effects:
                sub_effect = competitive_effects["CYP1A1"].substrates.get("PhIP")
                if sub_effect is not None:
                    competition_multiplier = min(competition_multiplier, 1.0 + sub_effect.flux_change_fraction)

        if enable_gsh_depletion:
            gsh_penalty = _compute_gsh_detox_penalty(carcinogen, gsh_status.fraction_normal, genotypes)
        else:
            gsh_penalty = 1.0

        interaction_adjusted_risks[carcinogen] = _round(
            base_risk * induction_multiplier * competition_multiplier * gsh_penalty,
            3,
        )

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
    cyp2e1_effects = competitive_effects.get("CYP2E1")
    if cyp2e1_effects is not None:
        for substrate, result in cyp2e1_effects.substrates.items():
            percent_change = result.flux_change_fraction * 100
            if abs(percent_change) > 10:
                direction = "reduced" if percent_change < 0 else "increased"
                summary_parts.append(
                    f"CYP2E1 competition: {substrate} activation {direction} by {abs(percent_change):.0f}%."
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


def decompose_synergy(
    exposure_profile: dict[str, float | dict[str, Any]],
    *,
    genotypes: dict[str, str] | None = None,
    tissue: str = "Liver",
    lifestyle: dict[str, bool] | None = None,
) -> dict[str, SynergyDecomposition]:
    """Decompose pairwise synergy into competition, GSH, and induction deltas.

    Runs ``compute_interaction_matrix`` four times (full, competition-only,
    GSH-only, induction-only) and reports
    ``S = 1 + ΔS_comp + ΔS_gsh + ΔS_ind + residual`` for each carcinogen pair
    in the exposure profile.
    """
    full_result = compute_interaction_matrix(
        exposure_profile,
        genotypes=genotypes,
        tissue=tissue,
        lifestyle=lifestyle,
    )
    comp_only = compute_interaction_matrix(
        exposure_profile,
        genotypes=genotypes,
        tissue=tissue,
        lifestyle=lifestyle,
        enable_induction=False,
        enable_competition=True,
        enable_gsh_depletion=False,
    )
    gsh_only = compute_interaction_matrix(
        exposure_profile,
        genotypes=genotypes,
        tissue=tissue,
        lifestyle=lifestyle,
        enable_induction=False,
        enable_competition=False,
        enable_gsh_depletion=True,
    )
    ind_only = compute_interaction_matrix(
        exposure_profile,
        genotypes=genotypes,
        tissue=tissue,
        lifestyle=lifestyle,
        enable_induction=True,
        enable_competition=False,
        enable_gsh_depletion=False,
    )

    decomposed: dict[str, SynergyDecomposition] = {}
    for pair, composite in full_result.synergy_matrix.items():
        delta_comp = comp_only.synergy_matrix.get(pair, 1.0) - 1.0
        delta_gsh = gsh_only.synergy_matrix.get(pair, 1.0) - 1.0
        delta_ind = ind_only.synergy_matrix.get(pair, 1.0) - 1.0
        additive = 1.0 + delta_comp + delta_gsh + delta_ind
        decomposed[pair] = SynergyDecomposition(
            pair=pair,
            composite=_round(composite, 4),
            delta_comp=_round(delta_comp, 4),
            delta_gsh=_round(delta_gsh, 4),
            delta_ind=_round(delta_ind, 4),
            additive_estimate=_round(additive, 4),
            residual=_round(composite - additive, 4),
        )
    return decomposed


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
    """Bootstrap 95% CIs for pairwise synergy factors and mechanism deltas.

    Km/Vmax parameters are perturbed by multiplicative lognormal noise with
    ``sigma=km_sigma`` (default ≈ ±50% one-sigma spread); enzyme expression
    weights are perturbed with ``sigma=expression_sigma`` (default ≈ ±30%).
    For each iteration the synergy is decomposed and the aggregated draws are
    summarized as mean and 2.5/97.5 percentile bounds per carcinogen pair.
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

        def _run(**flags: bool) -> dict[str, float]:
            return compute_interaction_matrix(
                exposure_profile,
                genotypes=genotypes,
                tissue=tissue,
                lifestyle=lifestyle,
                param_perturbations=param_perturbations,
                expression_perturbations=expression_perturbations,
                **flags,
            ).synergy_matrix

        full = _run()
        comp = _run(enable_induction=False, enable_competition=True, enable_gsh_depletion=False)
        gsh = _run(enable_induction=False, enable_competition=False, enable_gsh_depletion=True)
        ind = _run(enable_induction=True, enable_competition=False, enable_gsh_depletion=False)

        for pair, composite in full.items():
            composite_draws.setdefault(pair, []).append(composite)
            delta_comp_draws.setdefault(pair, []).append(comp.get(pair, 1.0) - 1.0)
            delta_gsh_draws.setdefault(pair, []).append(gsh.get(pair, 1.0) - 1.0)
            delta_ind_draws.setdefault(pair, []).append(ind.get(pair, 1.0) - 1.0)

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
            }
            for substrate, flux in result.substrates.items()
        }
        for enzyme, result in competitive_effects.items()
    }


def _interaction_matrix_to_compat_dict(result: InteractionMatrixResult) -> dict[str, Any]:
    """Convert an interaction result into a source-style JSON-serializable dict."""
    return {
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
    }


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
