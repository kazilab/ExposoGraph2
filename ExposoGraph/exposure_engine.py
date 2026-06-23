"""Exposure integration engine for Gene x Environment risk scoring.

Combines genotype-based flux ratios with real-world exposure scenarios
across 14 carcinogen classes to produce exposure-weighted cancer risk
estimates.  Includes EPA lifetime excess cancer risk (LECR) calculations,
questionnaire-driven profiling, and scenario comparison utilities.

Parts of this code were adapted from ExposoGraph with assistance from
ChatGPT Codex and Claude Code.
"""

from __future__ import annotations

import argparse
import json
import sys
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Sequence, TypeAlias, cast

# ── Enums ─────────────────────────────────────────────────────────────────────


class ExposureCarcinogenClass(str, Enum):
    """All 14 carcinogen classes supported by exposure integration."""

    PAH = "PAH"
    HCA = "HCA"
    AROMATIC_AMINES = "AromaticAmines"
    TOBACCO_NITROSAMINES = "TobaccoNitrosamines"
    AFLATOXIN_B1 = "AflatoxinB1"
    ESTROGEN_METABOLITES = "EstrogenMetabolites"
    BENZENE = "Benzene"
    VINYL_CHLORIDE = "VinylChloride"
    HEAVY_METALS = "HeavyMetals"
    ALDEHYDE = "Aldehyde"
    DIOXINS_AHR = "Dioxins_AhR"
    DIETARY_NITROSO = "DietaryNitroso"
    CHLORINATED_SOLVENTS = "ChlorinatedSolvents"
    UV_RADIATION = "UV_Radiation"


class ExposureRiskCategory(str, Enum):
    """Risk categories for exposure-weighted risk scores.

    Thresholds (applied to combined risk score):
        score <= 0.8   -> LOW
        0.8 < score <= 2.0  -> AVERAGE
        2.0 < score <= 10.0 -> ELEVATED
        score > 10.0   -> HIGH
    """

    LOW = "Low Risk"
    AVERAGE = "Population Average"
    ELEVATED = "Elevated Risk"
    HIGH = "High Risk"


JsonDict: TypeAlias = dict[str, Any]


# ── Dataclasses ───────────────────────────────────────────────────────────────


@dataclass
class ExposureScenarioInfo:
    """A single exposure scenario for a carcinogen class."""

    scenario_id: str
    label: str
    multiplier: float
    tissue_conc_uM: float
    source: str
    tier: int = 1
    tier_label: str = "Tier 1 population baseline"


@dataclass
class ExposureWeightedRisk:
    """Result of compute_exposure_weighted_risk."""

    carcinogen_class: str
    class_label: str
    exposure_scenario: str
    scenario_label: str
    tissue: str
    flux_ratio: float
    exposure_multiplier: float
    exposure_tier: int
    exposure_tier_label: str
    tissue_factor: float
    combined_risk_score: float
    risk_category: ExposureRiskCategory
    tissue_conc_uM: float
    biomarker_dose_estimate: dict[str, Any] | None
    interpretation: str
    sources: list[str] = field(default_factory=list)


@dataclass
class LifetimeCancerRisk:
    """Result of compute_lifetime_cancer_risk (EPA LECR)."""

    carcinogen_class: str
    tissue: str
    daily_dose_mg_kg: float
    duration_years: float
    slope_factor: float
    genotype_modifier: float
    lecr: float
    exceeds_action_threshold: bool  # >= 1e-4
    exceeds_de_minimis: bool  # >= 1e-6


@dataclass
class ExposureProfileResult:
    """Result of apply_exposure_profile across all 14 classes."""

    tissue: str
    scenario_mapping: dict[str, str]
    risk_scores: list[ExposureWeightedRisk]
    n_classes_assessed: int
    max_risk_class: str | None
    max_risk_score: float | None
    high_risk_classes: list[str]
    elevated_risk_classes: list[str]
    mean_combined_risk: float | None
    category_counts: dict[str, int]
    exposure_answers_used: dict[str, Any] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)

    @property
    def patient_tissue(self) -> str:
        """Compatibility alias used by the original exposure module."""
        return self.tissue

    @property
    def summary(self) -> dict[str, Any]:
        """Expose source-style summary fields alongside the typed API."""
        return {
            "n_classes_assessed": self.n_classes_assessed,
            "max_risk_class": self.max_risk_class,
            "max_risk_score": self.max_risk_score,
            "high_risk_classes": list(self.high_risk_classes),
            "elevated_risk_classes": list(self.elevated_risk_classes),
            "mean_combined_risk": self.mean_combined_risk,
            "category_counts": dict(self.category_counts),
        }


@dataclass
class ScenarioComparison:
    """Result of compare_scenarios for a single carcinogen class."""

    carcinogen_class: str
    tissue: str
    genotypes: dict[str, str]
    scenario_a: ExposureWeightedRisk
    scenario_b: ExposureWeightedRisk
    fold_change_b_vs_a: float
    absolute_delta: float
    interpretation: str


@dataclass
class ProfileComparison:
    """Result of full_profile_comparison across all classes."""

    tissue: str
    genotypes: dict[str, str]
    profile_a_name: str
    profile_b_name: str
    comparison_table: list[dict[str, Any]]
    profile_a_total_score: float
    profile_b_total_score: float
    overall_fold_change: float | None
    profile_a_high_risk_classes: list[str]
    profile_b_high_risk_classes: list[str]

    @property
    def patient_tissue(self) -> str:
        """Compatibility alias used by the original exposure module."""
        return self.tissue

    @property
    def aggregate(self) -> dict[str, Any]:
        """Expose source-style aggregate fields alongside the typed API."""
        return {
            f"{self.profile_a_name}_total_score": self.profile_a_total_score,
            f"{self.profile_b_name}_total_score": self.profile_b_total_score,
            "overall_fold_change_B_vs_A": self.overall_fold_change,
            f"{self.profile_a_name}_high_risk_classes": list(
                self.profile_a_high_risk_classes
            ),
            f"{self.profile_b_name}_high_risk_classes": list(
                self.profile_b_high_risk_classes
            ),
        }


# ── Data loading (lazy cache) ────────────────────────────────────────────────

_EXPOSURE_DB_FILE = Path(__file__).parent / "data" / "exposure_database.json"
_EXPOSURE_CACHE: JsonDict | None = None


def _load_exposure_db() -> JsonDict:
    """Load the exposure database JSON from disk."""
    with open(_EXPOSURE_DB_FILE, encoding="utf-8") as f:
        return cast(JsonDict, json.load(f))


def _get_exposure_db() -> JsonDict:
    """Return the cached exposure database (lazy-loaded)."""
    global _EXPOSURE_CACHE
    if _EXPOSURE_CACHE is None:
        _EXPOSURE_CACHE = _load_exposure_db()
    return _EXPOSURE_CACHE


def get_database() -> JsonDict:
    """Return a defensive copy of the packaged exposure database.

    This mirrors the original ``05_exposure`` helper while protecting the
    package-level lazy cache from accidental mutation by callers.
    """
    return deepcopy(_get_exposure_db())


# ── Private helpers ───────────────────────────────────────────────────────────


def _resolve_class_key(query: str, classes: Mapping[str, Any]) -> str | None:
    """Case-insensitive / fuzzy key resolution for carcinogen class names."""
    if query in classes:
        return query
    query_lower = query.lower().replace(" ", "").replace("-", "")
    aliases: dict[str, str] = {
        "pah": "PAH",
        "polycyclicaromatichydrocarbon": "PAH",
        "benzoapyrene": "PAH",
        "bap": "PAH",
        "hca": "HCA",
        "heterocyclicamine": "HCA",
        "phip": "HCA",
        "aromaticamine": "AromaticAmines",
        "aromaticamines": "AromaticAmines",
        "4abp": "AromaticAmines",
        "tobacconitrosamine": "TobaccoNitrosamines",
        "tobacconitrosamines": "TobaccoNitrosamines",
        "nnk": "TobaccoNitrosamines",
        "tsna": "TobaccoNitrosamines",
        "aflatoxin": "AflatoxinB1",
        "aflatoxinb1": "AflatoxinB1",
        "afb1": "AflatoxinB1",
        "estrogen": "EstrogenMetabolites",
        "estrogenmetabolites": "EstrogenMetabolites",
        "catecholestrogen": "EstrogenMetabolites",
        "benzene": "Benzene",
        "vinylchloride": "VinylChloride",
        "vcm": "VinylChloride",
        "aldehyde": "Aldehyde",
        "aldehydes": "Aldehyde",
        "acetaldehyde": "Aldehyde",
        "formaldehyde": "Aldehyde",
        "alcohol": "Aldehyde",
        "dioxin": "Dioxins_AhR",
        "dioxins": "Dioxins_AhR",
        "ahr": "Dioxins_AhR",
        "tcdd": "Dioxins_AhR",
        "dietarynitroso": "DietaryNitroso",
        "ndma": "DietaryNitroso",
        "nitrosoamine": "DietaryNitroso",
        "chlorinatedsolvent": "ChlorinatedSolvents",
        "chlorinatedsolvents": "ChlorinatedSolvents",
        "tce": "ChlorinatedSolvents",
        "trichloroethylene": "ChlorinatedSolvents",
        "heavymetal": "HeavyMetals",
        "heavymetals": "HeavyMetals",
        "arsenic": "HeavyMetals",
        "cadmium": "HeavyMetals",
        "chromium": "HeavyMetals",
        "uv": "UV_Radiation",
        "uvradiation": "UV_Radiation",
        "ultraviolet": "UV_Radiation",
    }

    if query_lower in aliases:
        target = aliases[query_lower]
        if target in classes:
            return target

    # Partial match fallback
    for key in classes:
        if query_lower in key.lower().replace("_", "").replace(" ", ""):
            return key

    return None


def _resolve_scenario_key(query: str, scenarios: Mapping[str, Any]) -> str | None:
    """Resolve a scenario key with fuzzy matching."""
    if query in scenarios:
        return query
    q = query.lower().replace(" ", "_").replace("-", "_")
    for k in scenarios:
        if k.lower() == q:
            return k
    for k in scenarios:
        if q in k.lower() or k.lower() in q:
            return k
    # Universal fallback for 'general_population'
    universal_fallback_priority = [
        "general_population",
        "us_low_exposure",
        "premenopausal_no_exogenous",
        "light_drinker",
        "general_population_low_processed_meat",
        "minimal",
        "unexposed",
        "never_smoker",
        "nonsmoker_no_occupational",
    ]
    if q == "general_population":
        for fb in universal_fallback_priority:
            if fb in scenarios:
                return fb
    return None


_EXPOSURE_TIER_LABELS = {
    1: "Tier 1 population baseline",
    2: "Tier 2 lifestyle adjustment",
    3: "Tier 3 occupational/environmental or measured biomarker",
}


def _scenario_tier(
    scenario_id: str, scenario_data: Mapping[str, Any]
) -> tuple[int, str]:
    """Classify a scenario into the manuscript three-tier exposure framework."""
    explicit = scenario_data.get("exposure_tier")
    if explicit is not None:
        try:
            tier = int(explicit)
        except (TypeError, ValueError):
            tier = 1
        tier = min(3, max(1, tier))
        return tier, _EXPOSURE_TIER_LABELS[tier]

    text = f"{scenario_id} {scenario_data.get('label', '')}".lower()
    tier1_ids = {
        "general_population",
        "us_low_exposure",
        "premenopausal_no_exogenous",
        "minimal",
        "unexposed",
        "never_smoker",
        "nonsmoker_no_occupational",
        "nondrinker",
        "vegetarian_low_meat",
        "general_population_low_processed_meat",
    }
    if scenario_id in tier1_ids:
        return 1, _EXPOSURE_TIER_LABELS[1]

    tier3_tokens = (
        "occupational",
        "worker",
        "workplace",
        "industrial",
        "superfund",
        "contaminated",
        "well_water",
        "urban_high_traffic",
        "arsenic",
        "chromium",
        "cadmium",
        "petroleum",
        "foundry",
        "coke",
        "dry cleaning",
        "rubber",
        "tanning",
        "formaldehyde",
    )
    tier2_tokens = (
        "smoker",
        "smoke",
        "grilled",
        "meat",
        "drinker",
        "alcohol",
        "diet",
        "fish",
        "seafood",
        "tanning_bed",
        "outdoor",
        "processed",
    )
    tier1_tokens = (
        "general_population",
        "typical",
        "minimal",
        "unexposed",
        "never",
        "nonsmoker",
        "low",
        "vegetarian",
        "nondrinker",
    )

    if any(token in text for token in tier3_tokens):
        tier = 3
    elif any(token in text for token in tier2_tokens):
        tier = 2
    elif any(token in text for token in tier1_tokens):
        tier = 1
    else:
        tier = 1
    return tier, _EXPOSURE_TIER_LABELS[tier]


_BIOMARKER_CLASS_ALIASES: dict[str, set[str]] = {
    "PAH": {"PAH"},
    "HCA": {"HCA"},
    "AromaticAmines": {"AromaticAmines", "AromaticAmine"},
    "TobaccoNitrosamines": {"TobaccoNitrosamines", "Nitrosamine", "NNK", "NNN"},
    "AflatoxinB1": {"AflatoxinB1", "Aflatoxin"},
    "EstrogenMetabolites": {"EstrogenMetabolites", "Estrogen"},
    "Benzene": {"Benzene"},
    "VinylChloride": {"VinylChloride"},
    "HeavyMetals": {"HeavyMetals", "HeavyMetal"},
    "Aldehyde": {"Aldehyde", "Ethanol"},
    "Dioxins_AhR": {"Dioxins_AhR", "Dioxin"},
    "DietaryNitroso": {"DietaryNitroso", "Nitrosamine", "NDMA", "NDEA"},
    "ChlorinatedSolvents": {"ChlorinatedSolvents", "ChlorinatedSolvent"},
    "UV_Radiation": {"UV_Radiation"},
}


def _resolve_biomarker_dose_for_class(
    carcinogen_class: str,
    biomarker_measurements: Mapping[str, float] | None,
    *,
    tissue: str | None = None,
) -> dict[str, Any] | None:
    """Return the strongest matching biomarker-derived dose estimate."""
    if not biomarker_measurements:
        return None

    from .biomarker_mapping import (
        biomarker_dose_estimate_to_dict,
        estimate_tissue_dose_from_biomarker,
        get_biomarker_entries,
    )

    aliases = {
        carcinogen_class,
        *(_BIOMARKER_CLASS_ALIASES.get(carcinogen_class, set())),
    }
    normalized_measurements = {
        str(name).strip().lower(): float(value)
        for name, value in biomarker_measurements.items()
    }

    estimates: list[dict[str, Any]] = []
    for entry in get_biomarker_entries():
        if entry.carcinogen_class not in aliases:
            continue
        measured = normalized_measurements.get(entry.biomarker.lower())
        if measured is None:
            continue
        estimate = estimate_tissue_dose_from_biomarker(
            entry.biomarker,
            measured,
            entry=entry,
        )
        estimates.append(biomarker_dose_estimate_to_dict(estimate))

    if not estimates:
        return None

    if tissue:
        tissue_key = tissue.replace(" ", "").replace("-", "_").lower()
        tissue_matches = [
            est
            for est in estimates
            if str(est["target_tissue"]).replace(" ", "").replace("-", "_").lower()
            == tissue_key
        ]
        if tissue_matches:
            estimates = tissue_matches

    return max(estimates, key=lambda est: float(est["s_over_km"]))


# put into JSON
_TISSUE_FACTORS: dict[tuple[str, str], float] = {
    ("PAH", "Lung"): 1.5,
    ("PAH", "Liver"): 0.8,
    ("PAH", "Colon"): 1.2,
    ("PAH", "Skin"): 1.3,
    ("HCA", "Colon"): 1.4,
    ("HCA", "Liver"): 1.2,
    ("HCA", "Breast"): 1.1,
    ("AromaticAmines", "Bladder"): 2.0,
    ("AromaticAmines", "Liver"): 1.0,
    ("TobaccoNitrosamines", "Lung"): 2.0,
    ("TobaccoNitrosamines", "Liver"): 0.8,
    ("AflatoxinB1", "Liver"): 2.5,
    ("AflatoxinB1", "Lung"): 0.5,
    ("EstrogenMetabolites", "Breast"): 2.0,
    ("EstrogenMetabolites", "Endometrium"): 1.8,
    ("Benzene", "Liver"): 1.5,
    ("Benzene", "BoneMarrow"): 1.8,
    ("Benzene", "Lung"): 1.2,
    ("VinylChloride", "Liver"): 2.5,
    ("Aldehyde", "Liver"): 2.0,
    ("Aldehyde", "Esophagus"): 1.8,
    ("Aldehyde", "Colon"): 1.2,
    ("DietaryNitroso", "Liver"): 2.0,
    ("DietaryNitroso", "Esophagus"): 1.5,
    ("DietaryNitroso", "Stomach"): 1.5,
    ("ChlorinatedSolvents", "Kidney"): 2.5,
    ("ChlorinatedSolvents", "Liver"): 1.5,
    ("HeavyMetals", "Lung"): 1.5,
    ("HeavyMetals", "Kidney"): 1.8,
    ("HeavyMetals", "Bladder"): 1.3,
    ("UV_Radiation", "Skin"): 3.0,
    ("UV_Radiation", "Eye"): 1.5,
    ("Dioxins_AhR", "Liver"): 1.5,
    ("Dioxins_AhR", "Lung"): 1.2,
}


def _get_tissue_factor(carcinogen_class: str, tissue: str) -> float:
    """Return tissue-specific expression weight for a carcinogen class."""
    class_key = _resolve_class_key(
        carcinogen_class,
        {k[0]: None for k in _TISSUE_FACTORS},
    )
    tissue_normalized = tissue.replace(" ", "").replace("-", "").lower()

    for (cls, tis), factor in _TISSUE_FACTORS.items():
        if (
            cls.lower() == (class_key or "").lower()
            and tis.lower() == tissue_normalized
        ):
            return factor

    return 1.0


def _classify_risk(score: float, thresholds: Mapping[str, Any]) -> ExposureRiskCategory:
    """Classify a combined risk score using database thresholds."""
    if score <= float(thresholds["low"]["max"]):
        return ExposureRiskCategory.LOW
    elif score <= float(thresholds["average"]["max"]):
        return ExposureRiskCategory.AVERAGE
    elif score <= float(thresholds["elevated"]["max"]):
        return ExposureRiskCategory.ELEVATED
    else:
        return ExposureRiskCategory.HIGH


def _interpret_risk(
    carcinogen_class: str,
    genotypes: dict[str, str],
    tissue: str,
    flux_ratio: float,
    exposure_mult: float,
    combined: float,
    category: ExposureRiskCategory,
    sc_label: str,
) -> str:
    """Generate a plain-language interpretation of the combined risk score."""
    geno_parts = []
    for gene, status in genotypes.items():
        if status.lower() not in ("*1/*1", "wt", "wildtype", "normal", "present"):
            geno_parts.append(f"{gene} {status}")
    geno_str = (", ".join(geno_parts) or "wildtype") + " genotype"

    if combined < 0.5:
        action = "Minimal concern -- exposure level is very low regardless of genotype."
    elif category == ExposureRiskCategory.LOW:
        action = (
            "Genotype provides some protection or exposure is below population average."
        )
    elif category == ExposureRiskCategory.AVERAGE:
        action = "Standard prevention and surveillance recommendations apply."
    elif category == ExposureRiskCategory.ELEVATED:
        action = (
            "Consider counseling on exposure reduction. "
            f"Genotype ({geno_str}) amplifies susceptibility at this exposure level."
        )
    else:
        action = (
            f"HIGH PRIORITY: {geno_str} combined with {sc_label} exposure creates "
            f"substantially elevated risk in {tissue}. Recommend urgent exposure "
            "reduction and enhanced clinical surveillance."
        )

    return (
        f"{category.value} | {carcinogen_class} -> {tissue}. "
        f"Flux ratio (genotype): {flux_ratio:.2f}x baseline. "
        f"Exposure scenario: '{sc_label}' ({exposure_mult:.1f}x population average). "
        f"Combined score: {combined:.2f}. {action}"
    )


def _answers_to_scenario_keys(answers: dict[str, Any]) -> dict[str, str]:
    """Map questionnaire answers to scenario keys for each carcinogen class."""
    mapping: dict[str, str] = {}

    smoking = answers.get("smoking_status", "never")
    grilled = answers.get("grilled_meat_frequency", "never")
    occupation = answers.get("occupation_category", "office_indoor")
    occ_chemicals = answers.get("occupational_chemicals", [])

    # PAH
    if "PAHs_coke_oven" in occ_chemicals or occupation in ["mining_metals"]:
        mapping["PAH"] = "occupational_high"
    elif smoking == "current" and grilled in ["weekly", "daily"]:
        mapping["PAH"] = "smoker_heavy_grilled_meat"
    elif smoking == "current":
        mapping["PAH"] = "smoker"
    elif grilled in ["weekly", "daily"]:
        mapping["PAH"] = "heavy_grilled_meat_diet"
    else:
        mapping["PAH"] = "general_population"

    # HCA
    red_meat = answers.get("red_meat_servings", 3)
    try:
        red_meat = float(red_meat)
    except (TypeError, ValueError):
        red_meat = 3

    if red_meat >= 10 and grilled in ["daily", "weekly"]:
        mapping["HCA"] = "high_meat_well_done"
    elif red_meat >= 5 or grilled in ["weekly", "daily"]:
        mapping["HCA"] = "moderate_meat_diet"
    elif red_meat <= 1:
        mapping["HCA"] = "vegetarian_low_meat"
    else:
        mapping["HCA"] = "general_population"

    # Aromatic Amines
    if occupation in ["petroleum_chemical", "manufacturing_industrial"]:
        mapping["AromaticAmines"] = "occupational_dye_rubber"
    elif smoking == "current":
        mapping["AromaticAmines"] = "smoker"
    else:
        mapping["AromaticAmines"] = "general_population"

    # Tobacco Nitrosamines
    cpd = answers.get("cigarettes_per_day", 0)
    try:
        cpd = float(cpd)
    except (TypeError, ValueError):
        cpd = 0

    if smoking == "current" and cpd >= 40:
        mapping["TobaccoNitrosamines"] = "heavy_smoker"
    elif smoking == "current" and cpd >= 20:
        mapping["TobaccoNitrosamines"] = "smoker_one_pack"
    elif smoking == "current":
        mapping["TobaccoNitrosamines"] = "smoker_half_pack"
    elif answers.get("secondhand_smoke", False):
        mapping["TobaccoNitrosamines"] = "secondhand_smoke"
    else:
        mapping["TobaccoNitrosamines"] = "never_smoker"

    # Aflatoxin
    if answers.get("aflatoxin_exposure", False):
        mapping["AflatoxinB1"] = "developing_country_moderate"
    else:
        mapping["AflatoxinB1"] = "us_low_exposure"

    # Estrogen (endogenous only)
    mapping["EstrogenMetabolites"] = "premenopausal_no_exogenous"

    # Benzene
    if "benzene_petroleum" in occ_chemicals or occupation == "petroleum_chemical":
        mapping["Benzene"] = "occupational_petroleum"
    elif smoking == "current":
        mapping["Benzene"] = "smoker"
    else:
        urban = answers.get("urban_rural", "suburban")
        if urban == "urban_high_traffic":
            mapping["Benzene"] = "urban_high_traffic"
        else:
            mapping["Benzene"] = "general_population"

    # Aldehyde (alcohol)
    alcohol = answers.get("alcohol_status", "nondrinker")
    dpw = answers.get("drinks_per_week", 0)
    try:
        dpw = float(dpw)
    except (TypeError, ValueError):
        dpw = 0

    if alcohol == "nondrinker" or dpw == 0:
        mapping["Aldehyde"] = "nondrinker"
    elif dpw <= 7:
        mapping["Aldehyde"] = "light_drinker"
    elif dpw <= 14:
        mapping["Aldehyde"] = "moderate_drinker"
    else:
        mapping["Aldehyde"] = "heavy_drinker"

    if "formaldehyde" in occ_chemicals or occupation == "funeral_embalming":
        mapping["Aldehyde"] = "formaldehyde_occupational"

    # Heavy Metals
    if "arsenic" in occ_chemicals or answers.get("arsenic_endemic_area", False):
        mapping["HeavyMetals"] = "arsenic_endemic_well_water"
    elif "chromium_vi" in occ_chemicals:
        mapping["HeavyMetals"] = "occupational_chromium_vi"
    elif "cadmium_batteries" in occ_chemicals:
        mapping["HeavyMetals"] = "occupational_cadmium"
    elif smoking == "current":
        mapping["HeavyMetals"] = "smoker"
    else:
        mapping["HeavyMetals"] = "general_population"

    # Dietary Nitroso
    processed = answers.get("processed_meat_frequency", "monthly")
    water = answers.get("water_source", "municipal_treated")

    if processed in ["daily"] and occupation in [
        "manufacturing_industrial",
        "petroleum_chemical",
    ]:
        mapping["DietaryNitroso"] = "occupational_rubber_tanning"
    elif processed == "daily":
        mapping["DietaryNitroso"] = "high_processed_meat"
    elif water in ["well_water_untested"]:
        mapping["DietaryNitroso"] = "contaminated_water_above_advisory"
    else:
        mapping["DietaryNitroso"] = "general_population"

    # Chlorinated Solvents
    if "TCE_solvents" in occ_chemicals:
        mapping["ChlorinatedSolvents"] = "occupational_moderate"
    elif answers.get("proximity_industrial", False):
        mapping["ChlorinatedSolvents"] = "near_superfund_site"
    else:
        mapping["ChlorinatedSolvents"] = "general_population"

    # UV Radiation
    sun = answers.get("sun_exposure", "moderate")
    sun_map = {
        "minimal": "minimal",
        "moderate": "general_population",
        "high_outdoor_worker": "outdoor_worker",
        "tanning_bed_user": "tanning_bed_user",
    }
    mapping["UV_Radiation"] = sun_map.get(sun, "general_population")

    # Dioxins
    seafood = answers.get("seafood_frequency", 2)
    try:
        seafood = float(seafood)
    except (TypeError, ValueError):
        seafood = 2

    if "PCBs_dioxins" in occ_chemicals:
        mapping["Dioxins_AhR"] = "occupational_herbicide_agent_orange"
    elif seafood >= 5 or answers.get("proximity_industrial", False):
        mapping["Dioxins_AhR"] = "near_industrial_source"
    elif seafood >= 3:
        mapping["Dioxins_AhR"] = "high_fatty_fish_diet"
    else:
        mapping["Dioxins_AhR"] = "general_population"

    return mapping


# ── Public API ────────────────────────────────────────────────────────────────

# Ordered list used by apply_exposure_profile (VinylChloride excluded —
# no questionnaire mapping; falls back to general_population internally).
_ALL_ASSESSED_CLASSES = [
    "PAH",
    "HCA",
    "AromaticAmines",
    "TobaccoNitrosamines",
    "AflatoxinB1",
    "EstrogenMetabolites",
    "Benzene",
    "Aldehyde",
    "HeavyMetals",
    "DietaryNitroso",
    "ChlorinatedSolvents",
    "UV_Radiation",
    "Dioxins_AhR",
]


# ── Source-compatibility presets ─────────────────────────────────────────────

JHBUI10030_GENOTYPES: dict[str, str] = {
    "CYP2C19": "*1/*1",
    "CYP2C9": "*1/*1",
    "CYP2D6": "*2/*4",
    "CYP3A4": "*1/*1",
    "CYP3A5": "*3/*3",
    "CYP2B6": "*6/*6",
    "CYP4F2": "*1/*3",
    "GSTM1": "present",
    "GSTT1": "present",
    "GSTP1": "*A/*A",
    "NAT2": "*4/*4",
    "ALDH2": "*1/*1",
    "ADH1B": "*1/*1",
    "CYP1A1": "*1/*1",
    "CYP1A2": "*1/*1",
    "CYP2A6": "*1/*1",
    "CYP2E1": "*1/*1",
}

JHBUI10030_PROFILE_A: dict[str, Any] = {
    "smoking_status": "never",
    "cigarettes_per_day": 0,
    "secondhand_smoke": False,
    "alcohol_status": "light",
    "drinks_per_week": 4,
    "binge_episodes": 0,
    "grilled_meat_frequency": "monthly",
    "processed_meat_frequency": "monthly",
    "red_meat_servings": 3,
    "seafood_frequency": 2,
    "aflatoxin_exposure": False,
    "occupation_category": "office_indoor",
    "occupation_years": 10,
    "occupational_chemicals": ["none"],
    "ppe_use": "never",
    "urban_rural": "suburban",
    "proximity_industrial": False,
    "water_source": "municipal_treated",
    "arsenic_endemic_area": False,
    "indoor_air_quality": False,
    "sun_exposure": "moderate",
}

JHBUI10030_PROFILE_B: dict[str, Any] = {
    "smoking_status": "current",
    "cigarettes_per_day": 20,
    "secondhand_smoke": True,
    "alcohol_status": "moderate",
    "drinks_per_week": 10,
    "binge_episodes": 2,
    "grilled_meat_frequency": "daily",
    "processed_meat_frequency": "weekly",
    "red_meat_servings": 8,
    "seafood_frequency": 1,
    "aflatoxin_exposure": False,
    "occupation_category": "manufacturing_industrial",
    "occupation_years": 20,
    "occupational_chemicals": ["PAHs_coke_oven", "benzene_petroleum", "formaldehyde"],
    "ppe_use": "sometimes",
    "urban_rural": "urban_high_traffic",
    "proximity_industrial": True,
    "water_source": "well_water_untested",
    "arsenic_endemic_area": False,
    "indoor_air_quality": True,
    "sun_exposure": "high_outdoor_worker",
}


def _source_risk_result(
    result: ExposureWeightedRisk | dict[str, Any],
) -> dict[str, Any]:
    """Convert a typed risk result into the original dict-style payload."""
    if isinstance(result, dict):
        return result
    return {
        "carcinogen_class": result.carcinogen_class,
        "class_label": result.class_label,
        "exposure_scenario": result.exposure_scenario,
        "scenario_label": result.scenario_label,
        "tissue": result.tissue,
        "flux_ratio": result.flux_ratio,
        "exposure_multiplier": result.exposure_multiplier,
        "exposure_tier": result.exposure_tier,
        "exposure_tier_label": result.exposure_tier_label,
        "tissue_factor": result.tissue_factor,
        "combined_risk_score": result.combined_risk_score,
        "risk_category": result.risk_category.value,
        "tissue_conc_uM": result.tissue_conc_uM,
        "biomarker_dose_estimate": deepcopy(result.biomarker_dose_estimate),
        "interpretation": result.interpretation,
        "sources": list(result.sources),
    }


def _source_profile_result(
    result: ExposureProfileResult | dict[str, Any],
) -> dict[str, Any]:
    """Convert a typed profile result into the original dict-style payload."""
    if isinstance(result, dict):
        return result
    return {
        "patient_tissue": result.tissue,
        "exposure_answers_used": deepcopy(result.exposure_answers_used),
        "scenario_mapping": dict(result.scenario_mapping),
        "risk_scores": [_source_risk_result(r) for r in result.risk_scores],
        "summary": result.summary,
        "errors": dict(result.errors),
    }


def _source_profile_comparison(
    comparison: ProfileComparison | dict[str, Any],
) -> dict[str, Any]:
    """Convert a typed profile comparison into the original dict-style payload."""
    if isinstance(comparison, dict):
        return comparison
    return {
        "patient_tissue": comparison.tissue,
        "genotypes": deepcopy(comparison.genotypes),
        "profile_a_name": comparison.profile_a_name,
        "profile_b_name": comparison.profile_b_name,
        "comparison_table": deepcopy(comparison.comparison_table),
        "aggregate": comparison.aggregate,
    }


def get_exposure_scenarios(
    carcinogen_class: str,
) -> list[ExposureScenarioInfo]:
    """Return all available exposure scenarios for a carcinogen class."""
    db = _get_exposure_db()
    classes = db["carcinogen_classes"]

    key = _resolve_class_key(carcinogen_class, classes)
    if key is None:
        available = list(classes.keys())
        raise ValueError(
            f"Unknown carcinogen class: '{carcinogen_class}'. "
            f"Available classes: {available}"
        )

    cls_data = classes[key]
    scenarios: list[ExposureScenarioInfo] = []
    for sc_id, sc_data in cls_data.get("exposure_scenarios", {}).items():
        mult = (
            sc_data.get("multiplier_vs_baseline")
            or sc_data.get("multiplier_vs_baseline_chromium")
            or 1.0
        )
        conc = (
            sc_data.get("estimated_tissue_conc_uM")
            or sc_data.get("tissue_conc_uM")
            or 0.0
        )
        # Handle dict-typed concentrations (e.g. HeavyMetals multi-metal)
        if isinstance(conc, dict):
            conc = float(max(conc.values())) if conc else 0.0
        tier, tier_label = _scenario_tier(sc_id, sc_data)
        scenarios.append(
            ExposureScenarioInfo(
                scenario_id=sc_id,
                label=sc_data.get("label", sc_id),
                multiplier=float(mult),
                tissue_conc_uM=float(conc),
                source=sc_data.get("source", "See exposure_database.json"),
                tier=tier,
                tier_label=tier_label,
            )
        )
    return scenarios


def compute_flux_ratio(
    genotypes: dict[str, str],
    carcinogen_class: str,
    tissue: str,
) -> float:
    """Compute simplified genotype-based flux ratio (multiplier-based).

    Returns a float where 1.0 = population average; >1 = net activation bias.
    This is NOT Michaelis-Menten kinetics — use flux_engine for that.
    """
    db = _get_exposure_db()
    genotype_modifiers = db["flux_model_integration"]["genotype_modifiers"]

    flux = 1.0
    norm_geno = {k.upper(): v.lower() for k, v in genotypes.items()}

    for modifier_key, class_multipliers in genotype_modifiers.items():
        gene, _, variant_desc = modifier_key.partition("_")
        gene_upper = gene.upper()
        geno_val = norm_geno.get(gene_upper, "")

        match = False
        if "null" in modifier_key and "null" in geno_val:
            match = True
        elif "heterozygote" in modifier_key and any(
            x in geno_val for x in ["*1/*2", "hetero", "het"]
        ):
            match = True
        elif "homozygous_variant" in modifier_key and any(
            x in geno_val for x in ["*2/*2", "homo", "hom"]
        ):
            match = True
        elif "slow" in modifier_key and "slow" in geno_val:
            match = True
        elif "PM" in modifier_key and any(x in geno_val for x in ["pm", "poor"]):
            match = True
        elif "high_activity" in modifier_key and any(
            x in geno_val for x in ["high", "ultra", "rapid", "*1f", "*1/*1f"]
        ):
            match = True
        elif "low_activity" in modifier_key and any(
            x in geno_val for x in ["low", "y113h", "his113"]
        ):
            match = True

        if match:
            class_key = _resolve_class_key(
                carcinogen_class, {k: None for k in class_multipliers}
            )
            if class_key and class_key in class_multipliers:
                flux *= class_multipliers[class_key]

    tissue_factor = _get_tissue_factor(carcinogen_class, tissue)
    flux *= tissue_factor

    return round(flux, 4)


def classify_exposure_risk(score: float) -> ExposureRiskCategory:
    """Classify a combined risk score into an ExposureRiskCategory.

    Uses the standard thresholds from the exposure database.
    """
    db = _get_exposure_db()
    thresholds = db["flux_model_integration"]["risk_thresholds"]
    return _classify_risk(score, thresholds)


def compute_exposure_weighted_risk(
    carcinogen_class: str,
    genotypes: dict[str, str],
    tissue: str,
    *,
    exposure_scenario: str = "general_population",
    flux_ratio_override: float | None = None,
    biomarker_measurements: Mapping[str, float] | None = None,
) -> ExposureWeightedRisk:
    """Core GxE risk scoring: flux_ratio x exposure_multiplier x tissue_factor.

    Parameters
    ----------
    flux_ratio_override : float | None
        If provided, use this value instead of the internal simplified
        flux ratio.  Allows passing a flux_engine-derived ratio for the
        7 overlapping carcinogen classes.
    biomarker_measurements : Mapping[str, float] | None
        Optional measured biomarker values keyed by biomarker identifier.
        When a matching biomarker exists for the carcinogen class, the
        biomarker-derived [S]/Km value replaces the scenario multiplier and
        its derived tissue concentration replaces the scenario concentration.
    """
    db = _get_exposure_db()
    classes = db["carcinogen_classes"]
    thresholds = db["flux_model_integration"]["risk_thresholds"]

    class_key = _resolve_class_key(carcinogen_class, classes)
    if class_key is None:
        raise ValueError(f"Unknown carcinogen class: '{carcinogen_class}'")

    cls_data = classes[class_key]
    scenarios = cls_data.get("exposure_scenarios", {})

    sc_key = _resolve_scenario_key(exposure_scenario, scenarios)
    if sc_key is None:
        available = list(scenarios.keys())
        raise ValueError(
            f"Unknown scenario '{exposure_scenario}' for {class_key}. "
            f"Available: {available}"
        )
    sc_data = scenarios[sc_key]

    exposure_mult = float(
        sc_data.get("multiplier_vs_baseline")
        or sc_data.get("multiplier_vs_baseline_chromium")
        or 1.0
    )
    raw_conc = (
        sc_data.get("estimated_tissue_conc_uM") or sc_data.get("tissue_conc_uM") or 0.0
    )
    if isinstance(raw_conc, dict):
        tissue_conc = float(max(raw_conc.values())) if raw_conc else 0.0
    else:
        tissue_conc = float(raw_conc)

    tier, tier_label = _scenario_tier(sc_key, sc_data)
    biomarker_dose = _resolve_biomarker_dose_for_class(
        class_key,
        biomarker_measurements,
        tissue=tissue,
    )
    if biomarker_dose is not None:
        exposure_mult = float(biomarker_dose["s_over_km"])
        tissue_conc = float(biomarker_dose["tissue_conc_uM"])
        tier = 3
        tier_label = _EXPOSURE_TIER_LABELS[tier]

    tissue_factor = _get_tissue_factor(class_key, tissue)
    if flux_ratio_override is not None:
        flux_ratio = float(flux_ratio_override)
    else:
        tissue_adjusted_flux = compute_flux_ratio(genotypes, class_key, tissue)
        flux_ratio = (
            round(tissue_adjusted_flux / tissue_factor, 4)
            if tissue_factor != 0
            else tissue_adjusted_flux
        )

    combined_risk_score = round(flux_ratio * exposure_mult * tissue_factor, 4)

    category = _classify_risk(combined_risk_score, thresholds)

    interp = _interpret_risk(
        carcinogen_class=class_key,
        genotypes=genotypes,
        tissue=tissue,
        flux_ratio=flux_ratio,
        exposure_mult=exposure_mult,
        combined=combined_risk_score,
        category=category,
        sc_label=sc_data.get("label", sc_key),
    )

    return ExposureWeightedRisk(
        carcinogen_class=class_key,
        class_label=cls_data.get("class_label", class_key),
        exposure_scenario=sc_key,
        scenario_label=sc_data.get("label", sc_key),
        tissue=tissue,
        flux_ratio=flux_ratio,
        exposure_multiplier=exposure_mult,
        exposure_tier=tier,
        exposure_tier_label=tier_label,
        tissue_factor=tissue_factor,
        combined_risk_score=combined_risk_score,
        risk_category=category,
        tissue_conc_uM=tissue_conc,
        biomarker_dose_estimate=biomarker_dose,
        interpretation=interp,
        sources=[
            sc_data.get("source", ""),
            cls_data.get("epa_cancer_slope_factor", {}).get("source", ""),
            *(
                biomarker_dose.get("references", [])
                if biomarker_dose is not None
                else []
            ),
        ],
    )


_SLOPE_FACTOR_FALLBACKS: dict[str, float] = {
    "PAH": 1.0,
    "Benzene": 0.055,
    "VinylChloride": 0.72,
    "AflatoxinB1": 1.0,
    "DietaryNitroso": 51.0,
    "HeavyMetals": 1.5,
}


def _slope_factor_per_mg_kg_day(sf_data: Mapping[str, Any], fallback: float) -> float:
    """Return a slope factor normalized to the LECR dose unit: mg/kg-day."""
    raw_value = sf_data.get("value")
    if raw_value is None:
        return float(fallback)

    value = float(raw_value)
    unit = str(sf_data.get("unit", "")).lower().replace("µ", "μ")
    if "per μg/kg" in unit or "per ug/kg" in unit:
        return value * 1000.0
    if "per ng/kg" in unit:
        return value * 1_000_000.0
    return value


def compute_lifetime_cancer_risk(
    carcinogen_class: str,
    genotypes: dict[str, str],
    daily_dose_mg_kg: float,
    *,
    duration_years: float = 70.0,
    tissue: str = "Liver",
    flux_ratio_override: float | None = None,
) -> LifetimeCancerRisk:
    """Estimate lifetime excess cancer risk (LECR) using EPA slope factor.

    Formula: LECR = slope_factor x daily_dose x genotype_mod x (duration / 70)
    """
    db = _get_exposure_db()
    classes = db["carcinogen_classes"]
    class_key = _resolve_class_key(carcinogen_class, classes)
    if class_key is None:
        raise ValueError(f"Unknown carcinogen class: '{carcinogen_class}'")

    cls_data = classes[class_key]
    sf_data = cls_data.get("epa_cancer_slope_factor", {})
    slope_factor = _slope_factor_per_mg_kg_day(
        sf_data,
        _SLOPE_FACTOR_FALLBACKS.get(class_key, 0.1),
    )

    if flux_ratio_override is not None:
        genotype_mod = flux_ratio_override
    else:
        genotype_mod = compute_flux_ratio(genotypes, class_key, tissue)

    duration_fraction = min(duration_years / 70.0, 1.0)
    lecr = round(slope_factor * daily_dose_mg_kg * genotype_mod * duration_fraction, 8)

    return LifetimeCancerRisk(
        carcinogen_class=class_key,
        tissue=tissue,
        daily_dose_mg_kg=daily_dose_mg_kg,
        duration_years=duration_years,
        slope_factor=slope_factor,
        genotype_modifier=genotype_mod,
        lecr=lecr,
        exceeds_action_threshold=lecr >= 1e-4,
        exceeds_de_minimis=lecr >= 1e-6,
    )


def generate_exposure_questionnaire() -> dict[str, Any]:
    """Return the structured questionnaire schema from the database."""
    db = _get_exposure_db()
    return cast(JsonDict, db["exposure_questionnaire_schema"])


def apply_exposure_profile(
    patient_genotypes: dict[str, str],
    exposure_answers: dict[str, Any],
    *,
    tissue: str = "Liver",
    biomarker_measurements: Mapping[str, float] | None = None,
) -> ExposureProfileResult:
    """Map questionnaire answers to scenarios for all 14 classes, return risk scores.

    Optional biomarker measurements override scenario multipliers for matching
    classes using the curated biomarker -> tissue [S]/Km mapping.
    """
    scenario_mapping = _answers_to_scenario_keys(exposure_answers)

    risk_scores: list[ExposureWeightedRisk] = []
    errors: dict[str, str] = {}
    for cls in _ALL_ASSESSED_CLASSES:
        sc = scenario_mapping.get(cls, "general_population")
        try:
            result = compute_exposure_weighted_risk(
                carcinogen_class=cls,
                genotypes=patient_genotypes,
                tissue=tissue,
                exposure_scenario=sc,
                biomarker_measurements=biomarker_measurements,
            )
            risk_scores.append(result)
        except Exception as exc:
            errors[cls] = str(exc)

    scores = [r.combined_risk_score for r in risk_scores]
    categories = [r.risk_category for r in risk_scores]

    max_idx = scores.index(max(scores)) if scores else -1

    return ExposureProfileResult(
        tissue=tissue,
        scenario_mapping=scenario_mapping,
        risk_scores=risk_scores,
        n_classes_assessed=len(risk_scores),
        max_risk_class=risk_scores[max_idx].carcinogen_class if max_idx >= 0 else None,
        max_risk_score=max(scores) if scores else None,
        high_risk_classes=[
            r.carcinogen_class
            for r in risk_scores
            if r.risk_category == ExposureRiskCategory.HIGH
        ],
        elevated_risk_classes=[
            r.carcinogen_class
            for r in risk_scores
            if r.risk_category == ExposureRiskCategory.ELEVATED
        ],
        mean_combined_risk=round(sum(scores) / len(scores), 4) if scores else None,
        category_counts={
            "High Risk": sum(1 for c in categories if c == ExposureRiskCategory.HIGH),
            "Elevated Risk": sum(
                1 for c in categories if c == ExposureRiskCategory.ELEVATED
            ),
            "Population Average": sum(
                1 for c in categories if c == ExposureRiskCategory.AVERAGE
            ),
            "Low Risk": sum(1 for c in categories if c == ExposureRiskCategory.LOW),
        },
        exposure_answers_used=deepcopy(exposure_answers),
        errors=errors,
    )


def compare_scenarios(
    genotypes: dict[str, str],
    tissue: str,
    scenario_a: str,
    scenario_b: str,
    *,
    carcinogen_class: str = "PAH",
) -> ScenarioComparison:
    """Compare risk between two exposure scenarios for the same genotype."""
    result_a = compute_exposure_weighted_risk(
        carcinogen_class, genotypes, tissue, exposure_scenario=scenario_a
    )
    result_b = compute_exposure_weighted_risk(
        carcinogen_class, genotypes, tissue, exposure_scenario=scenario_b
    )

    score_a = result_a.combined_risk_score
    score_b = result_b.combined_risk_score

    if score_a == 0:
        ratio = float("inf") if score_b > 0 else 1.0
    else:
        ratio = score_b / score_a

    delta = score_b - score_a

    interpretation = (
        f"Scenario B ('{scenario_b}') has {ratio:.1f}x the risk of "
        f"Scenario A ('{scenario_a}') for {carcinogen_class} in {tissue}. "
        f"Risk score difference: {delta:+.3f}."
    )
    if ratio >= 10:
        interpretation += (
            " This is a clinically significant difference. "
            "Exposure reduction from A->B would substantially reduce risk."
        )
    elif ratio >= 2:
        interpretation += " Meaningful exposure-attributable risk difference."
    else:
        interpretation += " Modest difference; genotype dominates risk profile."

    return ScenarioComparison(
        carcinogen_class=carcinogen_class,
        tissue=tissue,
        genotypes=genotypes,
        scenario_a=result_a,
        scenario_b=result_b,
        fold_change_b_vs_a=round(ratio, 3),
        absolute_delta=round(delta, 4),
        interpretation=interpretation,
    )


def full_profile_comparison(
    genotypes: dict[str, str],
    tissue: str,
    profile_a_answers: dict[str, Any],
    profile_b_answers: dict[str, Any],
    *,
    profile_a_name: str = "Profile A",
    profile_b_name: str = "Profile B",
) -> ProfileComparison:
    """Compare two complete exposure profiles across all 14 classes."""
    profile_a = apply_exposure_profile(genotypes, profile_a_answers, tissue=tissue)
    profile_b = apply_exposure_profile(genotypes, profile_b_answers, tissue=tissue)

    scores_a = {r.carcinogen_class: r for r in profile_a.risk_scores}
    scores_b = {r.carcinogen_class: r for r in profile_b.risk_scores}

    all_classes = sorted(set(list(scores_a.keys()) + list(scores_b.keys())))

    comparison: list[dict[str, Any]] = []
    for cls in all_classes:
        ra = scores_a.get(cls)
        rb = scores_b.get(cls)
        sa = ra.combined_risk_score if ra else 0
        sb = rb.combined_risk_score if rb else 0
        fold = round(sb / sa, 2) if sa > 0 else None
        comparison.append(
            {
                "carcinogen_class": cls,
                f"{profile_a_name}_score": sa,
                f"{profile_a_name}_category": ra.risk_category.value if ra else "N/A",
                f"{profile_a_name}_scenario": ra.scenario_label if ra else "N/A",
                f"{profile_b_name}_score": sb,
                f"{profile_b_name}_category": rb.risk_category.value if rb else "N/A",
                f"{profile_b_name}_scenario": rb.scenario_label if rb else "N/A",
                "fold_change_B_vs_A": fold,
                "higher_risk_profile": (
                    profile_b_name
                    if sb > sa
                    else profile_a_name if sa > sb else "Equal"
                ),
            }
        )

    total_a = sum(r[f"{profile_a_name}_score"] for r in comparison)
    total_b = sum(r[f"{profile_b_name}_score"] for r in comparison)

    return ProfileComparison(
        tissue=tissue,
        genotypes=genotypes,
        profile_a_name=profile_a_name,
        profile_b_name=profile_b_name,
        comparison_table=comparison,
        profile_a_total_score=round(total_a, 3),
        profile_b_total_score=round(total_b, 3),
        overall_fold_change=round(total_b / total_a, 2) if total_a > 0 else None,
        profile_a_high_risk_classes=[
            r["carcinogen_class"]
            for r in comparison
            if r[f"{profile_a_name}_category"] == "High Risk"
        ],
        profile_b_high_risk_classes=[
            r["carcinogen_class"]
            for r in comparison
            if r[f"{profile_b_name}_category"] == "High Risk"
        ],
    )


# ── Formatting utilities ─────────────────────────────────────────────────────


def _risk_color_label(category: ExposureRiskCategory | str) -> str:
    """Source-compatible emoji risk label helper."""
    val = category.value if isinstance(category, ExposureRiskCategory) else category
    labels = {
        "High Risk": "🔴",
        "Elevated Risk": "🟠",
        "Population Average": "⚪",
        "Low Risk": "🟢",
    }
    return labels.get(val, "⚪")


def _risk_text_label(category: ExposureRiskCategory | str) -> str:
    """Map a risk category to a bracketed text label (no emoji)."""
    val = category.value if isinstance(category, ExposureRiskCategory) else category
    labels = {
        "High Risk": "[HIGH]",
        "Elevated Risk": "[ELEVATED]",
        "Population Average": "[AVERAGE]",
        "Low Risk": "[LOW]",
    }
    return labels.get(val, "[?]")


def format_exposure_results_table(
    results: list[ExposureWeightedRisk],
    *,
    title: str = "Exposure-Weighted Risk Assessment",
) -> str:
    """Format a list of risk results into a readable ASCII table."""
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append(f"  {title}")
    lines.append("=" * 78)
    lines.append(
        f"{'Class':<25} {'Scenario':<22} {'Flux':>6} {'Exp*':>6} "
        f"{'Score':>7} {'Category':<20}"
    )
    lines.append("-" * 78)

    for r in results:
        cls = r.carcinogen_class[:24]
        sc = r.scenario_label[:21]
        label = _risk_text_label(r.risk_category)
        lines.append(
            f"  {cls:<25} {sc:<22} {r.flux_ratio:>6.2f} {r.exposure_multiplier:>6.1f} "
            f"{r.combined_risk_score:>7.3f} {label} {r.risk_category.value}"
        )

    lines.append("=" * 78)
    return "\n".join(lines)


def format_exposure_comparison_table(comparison: ProfileComparison) -> str:
    """Format a profile comparison into a readable ASCII table."""
    pa = comparison.profile_a_name
    pb = comparison.profile_b_name
    lines: list[str] = []
    lines.append("=" * 90)
    lines.append(f"  EXPOSURE PROFILE COMPARISON -- {pa}  vs  {pb}")
    lines.append(f"  Tissue: {comparison.tissue}")
    lines.append("=" * 90)
    lines.append(
        f"{'Class':<25} {pa[:18]:<20} {'Cat-A':<18} "
        f"{pb[:18]:<20} {'Cat-B':<18} {'Fold':>5}"
    )
    lines.append("-" * 90)

    for row in comparison.comparison_table:
        cls = row["carcinogen_class"][:24]
        sa = row.get(f"{pa}_score", 0)
        ca = row.get(f"{pa}_category", "")
        sb = row.get(f"{pb}_score", 0)
        cb = row.get(f"{pb}_category", "")
        fold = row.get("fold_change_B_vs_A")
        fold_str = f"{fold:.1f}x" if fold is not None else " N/A"
        label_a = _risk_text_label(ca) if isinstance(ca, str) else _risk_text_label(ca)
        label_b = _risk_text_label(cb) if isinstance(cb, str) else _risk_text_label(cb)
        lines.append(
            f"  {cls:<25} {sa:<8.3f} {label_a} {ca[:16]:<16} "
            f"{sb:<8.3f} {label_b} {cb[:16]:<16} {fold_str:>5}"
        )

    lines.append("-" * 90)
    lines.append(
        f"  TOTAL SCORE:              {comparison.profile_a_total_score:>8.3f}"
        f"                      {comparison.profile_b_total_score:>8.3f}"
        f"   {comparison.overall_fold_change or 'N/A'}x"
    )
    lines.append("=" * 90)
    return "\n".join(lines)


def format_results_table(
    results: Sequence[ExposureWeightedRisk | dict[str, Any]],
    title: str = "Exposure-Weighted Risk Assessment",
) -> str:
    """Source-compatible formatter using the original helper name."""
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append(f"  {title}")
    lines.append("=" * 78)
    lines.append(
        f"{'Class':<25} {'Scenario':<22} {'Flux':>6} {'Exp×':>6} "
        f"{'Score':>7} {'Category':<20}"
    )
    lines.append("-" * 78)

    for result in results:
        row = _source_risk_result(result)
        if row.get("error"):
            lines.append(f"  {row['carcinogen_class']:<25} ERROR: {row['error'][:40]}")
            continue
        cat = row.get("risk_category", "")
        icon = _risk_color_label(cat)
        lines.append(
            f"  {row.get('carcinogen_class', '')[:24]:<25} "
            f"{row.get('scenario_label', '')[:21]:<22} "
            f"{row.get('flux_ratio', 0):>6.2f} "
            f"{row.get('exposure_multiplier', 0):>6.1f} "
            f"{row.get('combined_risk_score', 0):>7.3f} {icon} {cat}"
        )

    lines.append("=" * 78)
    return "\n".join(lines)


def format_comparison_table(comparison: ProfileComparison | dict[str, Any]) -> str:
    """Source-compatible formatter using the original helper name."""
    payload = _source_profile_comparison(comparison)
    pa = payload["profile_a_name"]
    pb = payload["profile_b_name"]
    aggregate = payload["aggregate"]

    lines: list[str] = []
    lines.append("=" * 90)
    lines.append(f"  EXPOSURE PROFILE COMPARISON — {pa}  vs  {pb}")
    lines.append(f"  Tissue: {payload['patient_tissue']}")
    lines.append("=" * 90)
    lines.append(
        f"{'Class':<25} {pa[:18]:<20} {'Cat-A':<18} "
        f"{pb[:18]:<20} {'Cat-B':<18} {'Fold':>5}"
    )
    lines.append("-" * 90)

    for row in payload["comparison_table"]:
        cat_a = row.get(f"{pa}_category", "")
        cat_b = row.get(f"{pb}_category", "")
        fold = row.get("fold_change_B_vs_A")
        lines.append(
            f"  {row['carcinogen_class'][:24]:<25} "
            f"{row.get(f'{pa}_score', 0):<8.3f} {_risk_color_label(cat_a)} {cat_a[:16]:<16} "
            f"{row.get(f'{pb}_score', 0):<8.3f} {_risk_color_label(cat_b)} {cat_b[:16]:<16} "
            f"{f'{fold:.1f}×' if fold is not None else ' N/A':>5}"
        )

    lines.append("-" * 90)
    lines.append(
        f"  TOTAL SCORE:              {aggregate[f'{pa}_total_score']:>8.3f}"
        f"                      {aggregate[f'{pb}_total_score']:>8.3f}"
        f"   {aggregate.get('overall_fold_change_B_vs_A') or 'N/A'}×"
    )
    lines.append("=" * 90)
    return "\n".join(lines)


def run_exposure_validation_cases() -> None:
    """Run the original exposure-module reference cases with packaged code."""
    print("=" * 78)
    print("  CASE 1: PAH Risk in Lung — GSTM1-null Nonsmoker vs Smoker")
    print("=" * 78)
    gstm1_null = {"GSTM1": "null"}
    case1_nonsmoker = _source_risk_result(
        compute_exposure_weighted_risk(
            "PAH", gstm1_null, "Lung", exposure_scenario="general_population"
        )
    )
    case1_smoker = _source_risk_result(
        compute_exposure_weighted_risk(
            "PAH", gstm1_null, "Lung", exposure_scenario="smoker"
        )
    )
    case1_smoker_grilled = _source_risk_result(
        compute_exposure_weighted_risk(
            "PAH", gstm1_null, "Lung", exposure_scenario="smoker_heavy_grilled_meat"
        )
    )
    print(
        format_results_table(
            [case1_nonsmoker, case1_smoker, case1_smoker_grilled],
            title="PAH Lung Validation",
        )
    )
    print(
        f"  Fold-change smoker vs baseline: {case1_smoker['combined_risk_score'] / case1_nonsmoker['combined_risk_score']:.1f}×"
    )
    print()
    print("=" * 78)
    print(
        "  CASE 2: Aldehyde Risk in Liver — ALDH2 *1/*2 Nondrinker vs Moderate Drinker"
    )
    print("=" * 78)
    aldh2_het = {"ALDH2": "*1/*2"}
    case2_nondrinker = _source_risk_result(
        compute_exposure_weighted_risk(
            "Aldehyde", aldh2_het, "Liver", exposure_scenario="nondrinker"
        )
    )
    case2_moderate = _source_risk_result(
        compute_exposure_weighted_risk(
            "Aldehyde", aldh2_het, "Liver", exposure_scenario="moderate_drinker"
        )
    )
    case2_wt = _source_risk_result(
        compute_exposure_weighted_risk(
            "Aldehyde",
            {"ALDH2": "*1/*1"},
            "Liver",
            exposure_scenario="moderate_drinker",
        )
    )
    print(
        format_results_table(
            [case2_nondrinker, case2_moderate, case2_wt],
            title="Aldehyde Liver Validation",
        )
    )
    print(
        "  Fold-change ALDH2*1/*2 moderate vs nondrinker: "
        f"{case2_moderate['combined_risk_score'] / max(case2_nondrinker['combined_risk_score'], 0.001):.1f}×"
    )
    print()
    print("=" * 78)
    print(
        "  CASE 3: Heavy Metals (Chromium VI) Risk in Lung — General Population vs Occupational Worker"
    )
    print("=" * 78)
    case3_general = _source_risk_result(
        compute_exposure_weighted_risk(
            "HeavyMetals", {}, "Lung", exposure_scenario="general_population"
        )
    )
    case3_occ = _source_risk_result(
        compute_exposure_weighted_risk(
            "HeavyMetals", {}, "Lung", exposure_scenario="occupational_chromium_vi"
        )
    )
    case3_occ_null = _source_risk_result(
        compute_exposure_weighted_risk(
            "HeavyMetals",
            {"GSTM1": "null"},
            "Lung",
            exposure_scenario="occupational_chromium_vi",
        )
    )
    print(
        format_results_table(
            [case3_general, case3_occ, case3_occ_null],
            title="Heavy Metals Lung Validation",
        )
    )
    print(
        f"  Fold-change occupational vs baseline: {case3_occ['combined_risk_score'] / case3_general['combined_risk_score']:.1f}×"
    )
    print()
    print("=" * 78)
    print("  CASE 4: Patient JHBUI-10030 — Full Profile A vs Profile B Comparison")
    print("=" * 78)
    comparison = full_profile_comparison(
        genotypes=JHBUI10030_GENOTYPES,
        tissue="Liver",
        profile_a_answers=JHBUI10030_PROFILE_A,
        profile_b_answers=JHBUI10030_PROFILE_B,
        profile_a_name="Profile_A",
        profile_b_name="Profile_B",
    )
    print(format_comparison_table(comparison))
    aggregate = comparison.aggregate
    print("  CASE 4 AGGREGATE STATISTICS:")
    print(
        f"  Profile A total risk score (all classes): {aggregate['Profile_A_total_score']:.3f}"
    )
    print(
        f"  Profile B total risk score (all classes): {aggregate['Profile_B_total_score']:.3f}"
    )
    print(
        f"  Overall fold-change Profile B vs A: {aggregate['overall_fold_change_B_vs_A']:.1f}×"
    )
    print()
    print("=" * 78)
    print("  SUPPLEMENTAL: Lifetime Excess Cancer Risk (LECR) Estimates")
    print("=" * 78)
    lecr_cases = [
        (
            "PAH",
            {"GSTM1": "null"},
            "Lung",
            0.000005,
            70,
            "GSTM1-null nonsmoker, PAH oral dose",
        ),
        (
            "PAH",
            {"GSTM1": "null"},
            "Lung",
            0.000015,
            40,
            "GSTM1-null smoker (40 yr exposure)",
        ),
        (
            "Aldehyde",
            {"ALDH2": "*1/*2"},
            "Liver",
            0.001,
            30,
            "ALDH2*1/*2 moderate drinker, 30 yr",
        ),
        (
            "DietaryNitroso",
            {},
            "Liver",
            0.0001,
            70,
            "High processed meat diet, wildtype",
        ),
        ("HeavyMetals", {}, "Lung", 0.01, 25, "Occupational Cr(VI) 25 yr, wildtype"),
    ]
    print(f"  {'Case Description':<45} {'LECR':>12}  {'vs Threshold':>15}")
    print(f"  {'-' * 45} {'-' * 12}  {'-' * 15}")
    for carcinogen_class, genotypes, tissue, dose, duration, description in lecr_cases:
        lecr = compute_lifetime_cancer_risk(
            carcinogen_class,
            genotypes,
            dose,
            duration_years=duration,
            tissue=tissue,
        )
        if lecr.lecr >= 1e-4:
            threshold_status = "⚠ EXCEEDS 1×10⁻⁴"
        elif lecr.lecr >= 1e-6:
            threshold_status = "⚠ above 1×10⁻⁶"
        else:
            threshold_status = "✓ below 1×10⁻⁶"
        print(f"  {description:<45} {lecr.lecr:>12.2e}  {threshold_status}")
    print()
    print("=" * 78)
    print("  VALIDATION COMPLETE — All 4 cases run successfully.")
    print("=" * 78)


def run_validation_cases() -> None:
    """Source-compatible alias for the exposure validation runner."""
    run_exposure_validation_cases()


def cli_main(argv: list[str] | None = None) -> int:
    """CLI entrypoint compatible with the original standalone exposure module."""
    parser = argparse.ArgumentParser(
        description="ExposoGraph exposure integration engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m ExposoGraph.exposure_cli --genotypes '{"GSTM1":"null","ALDH2":"*1/*2"}' --tissue Liver --class Aldehyde --scenario heavy_drinker
  python -m ExposoGraph.exposure_cli --genotypes '{"GSTM1":"null"}' --tissue Lung --scenario smoker --all-classes
  python -m ExposoGraph.exposure_cli --list-scenarios PAH
  python -m ExposoGraph.exposure_cli --lecr --class PAH --daily-dose 0.001 --duration 40
  python -m ExposoGraph.exposure_cli --patient JHBUI-10030 --full-comparison
  python -m ExposoGraph.exposure_cli --validate
        """,
    )
    parser.add_argument(
        "--genotypes",
        type=str,
        default="{}",
        help='JSON dict of genotypes, e.g. \'{"GSTM1":"null"}\'',
    )
    parser.add_argument(
        "--tissue", type=str, default="Liver", help="Target tissue (default: Liver)"
    )
    parser.add_argument(
        "--class",
        dest="carcinogen_class",
        type=str,
        default="PAH",
        help="Carcinogen class",
    )
    parser.add_argument(
        "--scenario",
        type=str,
        default="general_population",
        help="Exposure scenario key",
    )
    parser.add_argument(
        "--all-classes",
        action="store_true",
        help="Run all carcinogen classes with a given scenario",
    )
    parser.add_argument(
        "--list-scenarios",
        type=str,
        metavar="CLASS",
        help="List available scenarios for a class",
    )
    parser.add_argument(
        "--questionnaire",
        action="store_true",
        help="Print the exposure questionnaire schema",
    )
    parser.add_argument(
        "--lecr", action="store_true", help="Compute lifetime excess cancer risk"
    )
    parser.add_argument(
        "--daily-dose",
        type=float,
        default=0.001,
        help="Daily dose in mg/kg-day for --lecr",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=70.0,
        help="Exposure duration in years for --lecr",
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("SCENARIO_A", "SCENARIO_B"),
        help="Compare two exposure scenarios",
    )
    parser.add_argument(
        "--patient",
        type=str,
        help="Patient ID shortcut (currently supports JHBUI-10030)",
    )
    parser.add_argument(
        "--full-comparison",
        action="store_true",
        help="Run full profile A vs B comparison for patient",
    )
    parser.add_argument(
        "--validate", action="store_true", help="Run built-in exposure validation cases"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON instead of formatted tables",
    )

    args = parser.parse_args(argv)

    if args.validate:
        run_exposure_validation_cases()
        return 0

    try:
        genotypes = json.loads(args.genotypes)
    except json.JSONDecodeError as exc:
        print(f"ERROR: Invalid genotypes JSON: {exc}", file=sys.stderr)
        return 1

    if args.patient and args.patient.upper() == "JHBUI-10030":
        genotypes = deepcopy(JHBUI10030_GENOTYPES)

    if args.list_scenarios:
        scenarios = get_exposure_scenarios(args.list_scenarios)
        if args.json:
            print(json.dumps([asdict(s) for s in scenarios], indent=2, default=str))
        else:
            print(f"\nExposure scenarios for '{args.list_scenarios}':")
            print("-" * 60)
            for scenario in scenarios:
                print(
                    f"  {scenario.scenario_id:<35} ×{scenario.multiplier:.1f}  "
                    f"{scenario.label[:40]}"
                )
            print()
        return 0

    if args.questionnaire:
        print(json.dumps(generate_exposure_questionnaire(), indent=2))
        return 0

    if args.full_comparison and args.patient:
        profile_comparison = full_profile_comparison(
            genotypes=JHBUI10030_GENOTYPES,
            tissue=args.tissue,
            profile_a_answers=JHBUI10030_PROFILE_A,
            profile_b_answers=JHBUI10030_PROFILE_B,
            profile_a_name="Profile_A",
            profile_b_name="Profile_B",
        )
        if args.json:
            print(
                json.dumps(
                    _source_profile_comparison(profile_comparison),
                    indent=2,
                    default=str,
                )
            )
        else:
            print(format_comparison_table(profile_comparison))
        return 0

    if args.compare:
        scenario_a, scenario_b = args.compare
        scenario_comparison = compare_scenarios(
            genotypes=genotypes,
            tissue=args.tissue,
            scenario_a=scenario_a,
            scenario_b=scenario_b,
            carcinogen_class=args.carcinogen_class,
        )
        if args.json:
            print(json.dumps(asdict(scenario_comparison), indent=2, default=str))
        else:
            print(f"\n{'=' * 60}")
            print(f"  SCENARIO COMPARISON: {args.carcinogen_class} -> {args.tissue}")
            print(f"{'=' * 60}")
            print(f"  Scenario A: {scenario_comparison.scenario_a.scenario_label}")
            print(
                f"    Combined risk: {scenario_comparison.scenario_a.combined_risk_score:.4f}  "
                f"[{scenario_comparison.scenario_a.risk_category.value}]"
            )
            print(f"  Scenario B: {scenario_comparison.scenario_b.scenario_label}")
            print(
                f"    Combined risk: {scenario_comparison.scenario_b.combined_risk_score:.4f}  "
                f"[{scenario_comparison.scenario_b.risk_category.value}]"
            )
            print(
                f"\n  Fold-change B vs A: {scenario_comparison.fold_change_b_vs_a:.2f}×"
            )
            print(f"\n  Interpretation: {scenario_comparison.interpretation}\n")
        return 0

    if args.lecr:
        lecr = compute_lifetime_cancer_risk(
            carcinogen_class=args.carcinogen_class,
            genotypes=genotypes,
            daily_dose_mg_kg=args.daily_dose,
            duration_years=args.duration,
            tissue=args.tissue,
        )
        if args.json:
            print(json.dumps(asdict(lecr), indent=2, default=str))
        else:
            print("\nLifetime Excess Cancer Risk (LECR):")
            print(f"  Carcinogen class : {args.carcinogen_class}")
            print(f"  Tissue           : {args.tissue}")
            print(f"  Daily dose       : {args.daily_dose} mg/kg-day")
            print(f"  Duration         : {args.duration} years")
            print(f"  Genotypes        : {genotypes}")
            print(f"  LECR             : {lecr.lecr:.2e}  ({lecr.lecr:.8f})")
            print("  (EPA de minimis risk = 1×10⁻⁶; regulatory concern = 1×10⁻⁴)")
            if lecr.exceeds_action_threshold:
                print("  ⚠ EXCEEDS regulatory action threshold of 1 in 10,000")
            elif lecr.exceeds_de_minimis:
                print("  ⚠ Above de minimis threshold but below action level")
            else:
                print("  ✓ Below de minimis threshold (<1 in 1,000,000)")
            print()
        return 0

    if args.all_classes:
        results = [
            compute_exposure_weighted_risk(
                carcinogen_class,
                genotypes,
                args.tissue,
                exposure_scenario=args.scenario,
            )
            for carcinogen_class in _ALL_ASSESSED_CLASSES
        ]
        if args.json:
            print(json.dumps([asdict(r) for r in results], indent=2, default=str))
        else:
            print(
                format_results_table(
                    results,
                    title=f"All Classes — {args.tissue} | Scenario: {args.scenario}",
                )
            )
        return 0

    result = compute_exposure_weighted_risk(
        carcinogen_class=args.carcinogen_class,
        genotypes=genotypes,
        tissue=args.tissue,
        exposure_scenario=args.scenario,
    )
    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print(
            format_results_table(
                [result], title=f"{args.carcinogen_class} — {args.tissue}"
            )
        )
        print(f"\nInterpretation: {result.interpretation}\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Package CLI wrapper function."""
    return cli_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
