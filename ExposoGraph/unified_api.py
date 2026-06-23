"""Typed patient-level orchestration across the risk engines.

Provides one package-native entry point for patient risk assessment by
combining tissue context, genotype-specific flux modeling, exposure-weighted
risk, and multi-carcinogen interaction analysis into a single typed result.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from .chip_modifier import CHIPEffect, CHIPStatus
    from .epigenomic_modifier import EpigeneticEffect, MethylationStatus
    from .exposure_engine import ExposureProfileResult
    from .flux_engine import FullProfileResult
    from .interaction_engine import CriticalInteraction, InteractionMatrixResult
    from .models import KnowledgeGraph


_DEFAULT_LIFESTYLE = {
    "smoking": False,
    "alcohol_heavy": False,
    "alcohol_moderate": False,
    "occupational_exposure": False,
}

_FLUX_TO_EXPOSURE_CLASS = {
    "PAH": "PAH",
    "HCA": "HCA",
    "AromaticAmines": "AromaticAmines",
    "Aflatoxin": "AflatoxinB1",
    "EstrogenMetabolites": "EstrogenMetabolites",
    "Aldehyde": "Aldehyde",
    "Nitrosamine": "TobaccoNitrosamines",
    "NDMA": "DietaryNitroso",
    "NDEA": "DietaryNitroso",
    "Benzene": "Benzene",
    "VinylChloride": "VinylChloride",
    "UV_Radiation": "UV_Radiation",
    "Dioxin": "Dioxins_AhR",
    "ChlorinatedSolvent": "ChlorinatedSolvents",
    "HeavyMetal": "HeavyMetals",
}

_BIOMARKER_TO_FLUX_CLASS = {
    "PAH": "PAH",
    "HCA": "HCA",
    "AromaticAmines": "AromaticAmines",
    "AromaticAmine": "AromaticAmines",
    "TobaccoNitrosamines": "Nitrosamine",
    "NNK": "Nitrosamine",
    "NNN": "Nitrosamine",
    "Aflatoxin": "Aflatoxin",
    "AflatoxinB1": "Aflatoxin",
    "Ethanol": "Aldehyde",
    "Aldehyde": "Aldehyde",
    "Benzene": "Benzene",
    "VinylChloride": "VinylChloride",
    "DietaryNitroso": "NDMA",
    "NDMA": "NDMA",
    "NDEA": "NDEA",
    "ChlorinatedSolvent": "ChlorinatedSolvent",
    "ChlorinatedSolvents": "ChlorinatedSolvent",
    "HeavyMetal": "HeavyMetal",
    "HeavyMetals": "HeavyMetal",
}


@dataclass
class FluxClassEvidence:
    """Compact patient-facing summary of one flux class result."""

    carcinogen_class: str
    risk_classification: str
    net_ratio: float
    susceptibility_score_log2: float
    reactive_intermediate_uM: float | None
    time_to_steady_state_days: float | None
    model_kind: str
    parameter_source: str


@dataclass
class PatientRiskProfile:
    """Unified patient-level risk assessment result."""

    tissue: str
    genotypes: dict[str, str]
    lifestyle: dict[str, bool]
    exposure_scenario: str
    exposure_answers: dict[str, Any]
    tissue_report: str | None
    tissue_weight_count: int | None
    top_tissue_genes: list[tuple[str, float]]
    flux_profile: FullProfileResult | None
    exposure_profile: ExposureProfileResult | None
    interactions: InteractionMatrixResult | None
    biomarker_measurements: dict[str, float] = field(default_factory=dict)
    biomarker_dose_estimates: list[dict[str, Any]] = field(default_factory=list)
    flux_class_evidence: list[FluxClassEvidence] = field(default_factory=list)
    critical_warnings: list[CriticalInteraction] = field(default_factory=list)
    pipeline_warnings: list[str] = field(default_factory=list)
    pipeline_errors: dict[str, str] = field(default_factory=dict)
    chip_status: CHIPStatus | None = None
    chip_effect: CHIPEffect | None = None
    chip_adjusted_risks: dict[str, float] = field(default_factory=dict)
    methylation_status: MethylationStatus | None = None
    methylation_effect: EpigeneticEffect | None = None
    methylation_adjusted_risks: dict[str, float] = field(default_factory=dict)
    biological_output_integration: dict[str, Any] = field(default_factory=dict)
    summary: str = ""


@dataclass
class PatientProfileComparison:
    """Structured comparison between two patient risk profiles."""

    left_label: str
    right_label: str
    tissue: str
    left_total_interaction_risk: float | None
    right_total_interaction_risk: float | None
    total_interaction_risk_delta: float | None
    left_interaction_factor: float | None
    right_interaction_factor: float | None
    interaction_factor_delta: float | None
    left_high_risk_classes: list[str]
    right_high_risk_classes: list[str]
    shared_high_risk_classes: list[str]
    left_only_high_risk_classes: list[str]
    right_only_high_risk_classes: list[str]
    more_concerning_profile: str | None
    summary: str


_RISK_PRIORITY = {
    "HIGH": 0,
    "ELEVATED": 1,
    "MODERATE": 2,
    "LOW": 3,
    "PROTECTIVE": 4,
    "INSUFFICIENT_DATA": 5,
}


def _normalize_lifestyle(lifestyle: dict[str, Any] | None) -> dict[str, bool]:
    """Normalize lifestyle inputs to the expected boolean shape."""
    normalized = dict(_DEFAULT_LIFESTYLE)
    if not lifestyle:
        return normalized
    for key in normalized:
        normalized[key] = bool(lifestyle.get(key, False))
    return normalized


def _infer_lifestyle_from_answers(exposure_answers: dict[str, Any] | None) -> dict[str, bool]:
    """Infer coarse lifestyle flags from an exposure questionnaire payload."""
    if not exposure_answers:
        return dict(_DEFAULT_LIFESTYLE)

    smoking_status = str(exposure_answers.get("smoking_status", "")).lower()
    alcohol_status = str(exposure_answers.get("alcohol_status", "")).lower()
    occupation_category = str(exposure_answers.get("occupation_category", "")).lower()
    occupational_chemicals = exposure_answers.get("occupational_chemicals", [])

    occupational = bool(exposure_answers.get("proximity_industrial", False))
    occupational = occupational or "industrial" in occupation_category
    occupational = occupational or "manufacturing" in occupation_category
    occupational = occupational or bool(occupational_chemicals)

    return {
        "smoking": smoking_status in {"current", "daily", "some_day"},
        "alcohol_heavy": alcohol_status == "heavy",
        "alcohol_moderate": alcohol_status == "moderate",
        "occupational_exposure": occupational,
    }


def _derive_exposure_scenario(lifestyle: dict[str, bool]) -> str:
    """Map lifestyle factors to a compact named exposure scenario."""
    smoking = lifestyle.get("smoking", False)
    heavy_alcohol = lifestyle.get("alcohol_heavy", False)
    moderate_alcohol = lifestyle.get("alcohol_moderate", False)
    occupational = lifestyle.get("occupational_exposure", False)

    if smoking and heavy_alcohol and occupational:
        return "smoker_industrial_heavy_drinker"
    if smoking and heavy_alcohol:
        return "smoker_heavy_drinker"
    if smoking and occupational:
        return "smoker_industrial_worker"
    if smoking and moderate_alcohol:
        return "smoker_moderate_drinker"
    if smoking:
        return "smoker"
    if heavy_alcohol:
        return "heavy_drinker"
    if occupational:
        return "industrial_worker"
    if moderate_alcohol:
        return "moderate_drinker"
    return "general_population"


def _build_exposure_answers(
    exposure_scenario: str,
    lifestyle: dict[str, bool],
) -> dict[str, Any]:
    """Build a minimal questionnaire payload from a named scenario."""
    smoking = lifestyle.get("smoking", False)
    heavy_alcohol = lifestyle.get("alcohol_heavy", False)
    moderate_alcohol = lifestyle.get("alcohol_moderate", False)
    occupational = lifestyle.get("occupational_exposure", False)

    answers: dict[str, Any] = {
        "smoking_status": "current" if smoking else "never",
        "cigarettes_per_day": 25 if smoking else 0,
        "secondhand_smoke": False,
        "alcohol_status": "heavy" if heavy_alcohol else "moderate" if moderate_alcohol else "nondrinker",
        "drinks_per_week": 20 if heavy_alcohol else 8 if moderate_alcohol else 0,
        "grilled_meat_frequency": "weekly",
        "processed_meat_frequency": "monthly",
        "red_meat_servings": 4,
        "seafood_frequency": 1,
        "aflatoxin_exposure": False,
        "occupation_category": "manufacturing_industrial" if occupational else "office_indoor",
        "occupational_chemicals": [],
        "urban_rural": "suburban",
        "proximity_industrial": occupational,
        "water_source": "municipal_treated",
        "arsenic_endemic_area": False,
        "sun_exposure": "moderate",
    }

    if occupational:
        answers["occupational_chemicals"] = ["benzene_petroleum", "formaldehyde", "chromium_vi", "TCE_solvents"]
    if "heavy_drinker" in exposure_scenario:
        answers["processed_meat_frequency"] = "daily"
    if "smoker" in exposure_scenario:
        answers["secondhand_smoke"] = True
        answers["grilled_meat_frequency"] = "daily"

    return answers


def _build_interaction_exposure_profile(
    exposure_scenario: str,
    lifestyle: dict[str, bool],
    exposure_answers: dict[str, Any],
) -> dict[str, float]:
    """Build the compact interaction-engine exposure profile."""
    profile: dict[str, float] = {}

    if lifestyle.get("smoking", False) or "smoker" in exposure_scenario:
        profile.update(
            {
                "PAH": 3.0,
                "NNK": 4.0,
                "HCA": 1.5,
                "benzene": 6.0,
                "formaldehyde": 2.0,
                "cadmium": 2.0,
                "acrolein": 5.0,
            }
        )

    if lifestyle.get("alcohol_moderate", False) or "moderate_drinker" in exposure_scenario:
        profile["acetaldehyde"] = max(profile.get("acetaldehyde", 1.0), 1.5)
        profile["NDMA"] = max(profile.get("NDMA", 1.0), 1.2)
        profile["ethanol"] = max(profile.get("ethanol", 1.0), 3.0)

    if lifestyle.get("alcohol_heavy", False) or "heavy_drinker" in exposure_scenario:
        profile["acetaldehyde"] = max(profile.get("acetaldehyde", 1.0), 3.0)
        profile["NDMA"] = max(profile.get("NDMA", 1.0), 2.0)
        profile["benzene"] = max(profile.get("benzene", 1.0), 1.5)
        profile["ethanol"] = max(profile.get("ethanol", 1.0), 8.0)
        profile["acrolein"] = max(profile.get("acrolein", 1.0), 1.5)

    if lifestyle.get("occupational_exposure", False) or "industrial" in exposure_scenario:
        profile["chromium_VI"] = max(profile.get("chromium_VI", 1.0), 10.0)
        profile["benzene"] = max(profile.get("benzene", 1.0), 10.0)
        profile["formaldehyde"] = max(profile.get("formaldehyde", 1.0), 5.0)
        profile["PAH"] = max(profile.get("PAH", 1.0), 2.0)
        profile["cadmium"] = max(profile.get("cadmium", 1.0), 3.0)
        profile["vinyl_chloride"] = max(profile.get("vinyl_chloride", 1.0), 2.0)

    chemicals = [str(item).lower() for item in exposure_answers.get("occupational_chemicals", [])]
    if any("chromium" in item for item in chemicals):
        profile["chromium_VI"] = max(profile.get("chromium_VI", 1.0), 10.0)
    if any("tce" in item or "solvent" in item for item in chemicals):
        profile["vinyl_chloride"] = max(profile.get("vinyl_chloride", 1.0), 2.0)
        profile["benzene"] = max(profile.get("benzene", 1.0), 4.0)
    if any("benzene" in item for item in chemicals):
        profile["benzene"] = max(profile.get("benzene", 1.0), 10.0)
    if exposure_answers.get("arsenic_endemic_area", False):
        profile["arsenic"] = max(profile.get("arsenic", 1.0), 4.0)
    if exposure_answers.get("aflatoxin_exposure", False):
        profile["AFB1"] = max(profile.get("AFB1", 1.0), 4.0)

    processed_meat = str(exposure_answers.get("processed_meat_frequency", "")).lower()
    grilled_meat = str(exposure_answers.get("grilled_meat_frequency", "")).lower()
    if processed_meat in {"daily", "frequent", "high"}:
        profile["NDMA"] = max(profile.get("NDMA", 1.0), 2.0)
    if grilled_meat in {"daily", "frequent", "high"}:
        profile["HCA"] = max(profile.get("HCA", 1.0), 2.0)
        profile["PAH"] = max(profile.get("PAH", 1.0), 2.0)

    return profile or {"PAH": 1.0}


def _top_tissue_genes(weights: dict[str, float], n: int = 5) -> list[tuple[str, float]]:
    """Return the highest-weight tissue genes."""
    ranked = sorted(weights.items(), key=lambda item: (-item[1], item[0]))
    return [(gene, round(weight, 4)) for gene, weight in ranked[:n]]


def _rebuild_exposure_profile(
    tissue: str,
    scenario_mapping: dict[str, str],
    risk_scores: list[Any],
) -> ExposureProfileResult:
    """Rebuild an ExposureProfileResult after flux-based overrides."""
    from .exposure_engine import ExposureProfileResult, ExposureRiskCategory

    scores = [item.combined_risk_score for item in risk_scores]
    max_index = scores.index(max(scores)) if scores else -1

    return ExposureProfileResult(
        tissue=tissue,
        scenario_mapping=scenario_mapping,
        risk_scores=risk_scores,
        n_classes_assessed=len(risk_scores),
        max_risk_class=risk_scores[max_index].carcinogen_class if max_index >= 0 else None,
        max_risk_score=max(scores) if scores else None,
        high_risk_classes=[
            item.carcinogen_class
            for item in risk_scores
            if item.risk_category == ExposureRiskCategory.HIGH
        ],
        elevated_risk_classes=[
            item.carcinogen_class
            for item in risk_scores
            if item.risk_category == ExposureRiskCategory.ELEVATED
        ],
        mean_combined_risk=round(sum(scores) / len(scores), 4) if scores else None,
        category_counts={
            "High Risk": sum(1 for item in risk_scores if item.risk_category == ExposureRiskCategory.HIGH),
            "Elevated Risk": sum(1 for item in risk_scores if item.risk_category == ExposureRiskCategory.ELEVATED),
            "Population Average": sum(1 for item in risk_scores if item.risk_category == ExposureRiskCategory.AVERAGE),
            "Low Risk": sum(1 for item in risk_scores if item.risk_category == ExposureRiskCategory.LOW),
        },
    )


def _build_flux_exposure_profile_from_biomarkers(
    biomarker_measurements: dict[str, float] | None,
) -> tuple[dict[str, float], list[dict[str, Any]], list[str]]:
    """Convert biomarker measurements into flux-engine substrate overrides."""
    if not biomarker_measurements:
        return {}, [], []

    from .biomarker_mapping import (
        biomarker_dose_estimate_to_dict,
        estimate_tissue_dose_from_biomarker,
        get_entry_by_biomarker,
    )

    exposure_profile: dict[str, float] = {}
    estimates: list[dict[str, Any]] = []
    warnings: list[str] = []

    for biomarker, value in biomarker_measurements.items():
        entry = get_entry_by_biomarker(biomarker)
        if entry is None:
            warnings.append(f"Unknown biomarker measurement ignored: {biomarker}")
            continue

        estimate = estimate_tissue_dose_from_biomarker(
            entry.biomarker,
            float(value),
            entry=entry,
        )
        estimate_payload = biomarker_dose_estimate_to_dict(estimate)
        estimates.append(estimate_payload)

        flux_class = _BIOMARKER_TO_FLUX_CLASS.get(entry.carcinogen_class)
        if flux_class is None:
            warnings.append(
                f"No flux class mapping for biomarker {entry.biomarker} "
                f"({entry.carcinogen_class})"
            )
            continue

        existing = exposure_profile.get(flux_class)
        concentration = float(estimate.tissue_conc_uM)
        exposure_profile[flux_class] = (
            concentration if existing is None else max(existing, concentration)
        )

    return exposure_profile, estimates, warnings


def _build_exposure_profile(
    genotypes: dict[str, str],
    tissue: str,
    exposure_answers: dict[str, Any],
    flux_profile: FullProfileResult | None,
    biomarker_measurements: dict[str, float] | None = None,
) -> tuple[ExposureProfileResult, list[str]]:
    """Compute the exposure profile and override overlapping classes with flux results."""
    from .exposure_engine import apply_exposure_profile, compute_exposure_weighted_risk

    warnings: list[str] = []
    base_profile = apply_exposure_profile(
        patient_genotypes=genotypes,
        exposure_answers=exposure_answers,
        tissue=tissue,
        biomarker_measurements=biomarker_measurements,
    )
    if flux_profile is None:
        return base_profile, warnings

    flux_overrides: dict[str, float] = {}
    for flux_class, exposure_class in _FLUX_TO_EXPOSURE_CLASS.items():
        pathway_result = flux_profile.per_class_results.get(flux_class)
        if pathway_result is not None:
            prior = flux_overrides.get(exposure_class)
            flux_overrides[exposure_class] = (
                pathway_result.net_ratio
                if prior is None
                else max(prior, pathway_result.net_ratio)
            )

    if not flux_overrides:
        return base_profile, warnings

    overridden_scores: list[Any] = []
    for risk_score in base_profile.risk_scores:
        override = flux_overrides.get(risk_score.carcinogen_class)
        if override is None:
            overridden_scores.append(risk_score)
            continue
        try:
            overridden_scores.append(
                compute_exposure_weighted_risk(
                    carcinogen_class=risk_score.carcinogen_class,
                    genotypes=genotypes,
                    tissue=tissue,
                    exposure_scenario=risk_score.exposure_scenario,
                    flux_ratio_override=override,
                    biomarker_measurements=biomarker_measurements,
                )
            )
        except Exception as exc:  # pragma: no cover - defensive path
            warnings.append(
                f"Exposure override failed for {risk_score.carcinogen_class}: "
                f"{type(exc).__name__}: {exc}"
            )
            overridden_scores.append(risk_score)

    return _rebuild_exposure_profile(tissue, base_profile.scenario_mapping, overridden_scores), warnings


def _build_flux_class_evidence(flux_profile: FullProfileResult | None) -> list[FluxClassEvidence]:
    """Build a compact per-class metadata summary from the flux profile."""
    if flux_profile is None:
        return []

    entries = [
        FluxClassEvidence(
            carcinogen_class=class_name,
            risk_classification=result.risk_classification.value,
            net_ratio=result.net_ratio,
            susceptibility_score_log2=result.susceptibility_score_log2,
            reactive_intermediate_uM=result.steady_state_concentrations_uM.get(
                "reactive_intermediate_uM"
            ),
            time_to_steady_state_days=result.steady_state_model.get(
                "time_to_steady_state_days"
            ),
            model_kind=result.model_kind,
            parameter_source=result.parameter_source,
        )
        for class_name, result in flux_profile.per_class_results.items()
    ]
    return sorted(
        entries,
        key=lambda item: (
            _RISK_PRIORITY.get(item.risk_classification, 99),
            -item.net_ratio,
            item.carcinogen_class,
        ),
    )


def _collect_high_risk_classes(profile: PatientRiskProfile) -> list[str]:
    """Collect high-signal risk classes across the flux and exposure layers."""
    classes: set[str] = set()
    if profile.flux_profile is not None:
        classes.update(profile.flux_profile.elevated_or_high_risk_classes)
    if profile.exposure_profile is not None:
        classes.update(profile.exposure_profile.high_risk_classes)
        classes.update(profile.exposure_profile.elevated_risk_classes)
    return sorted(classes)


def _format_flux_model_kind(model_kind: str) -> str:
    """Render a flux model kind into compact patient-facing text."""
    mapping = {
        "measured_kinetics": "measured kinetics",
        "semi_quantitative_proxy": "semi-quantitative proxy",
        "damage_repair_proxy": "damage/repair proxy",
        "receptor_mediated_proxy": "receptor-mediated proxy",
        "unavailable": "no quantitative model",
    }
    return mapping.get(model_kind, model_kind.replace("_", " "))


def _normalize_graph_token(value: str | None) -> str:
    """Normalize graph labels/groups into a matching token."""
    if not value:
        return ""
    return "".join(char for char in value.lower() if char.isalnum())


def _graph_node_tokens(node: Any) -> list[str]:
    """Collect normalized matching tokens from a graph node."""
    return [
        token
        for token in {
            _normalize_graph_token(getattr(node, "id", None)),
            _normalize_graph_token(getattr(node, "label", None)),
            _normalize_graph_token(getattr(node, "group", None)),
            _normalize_graph_token(getattr(node, "canonical_id", None)),
            _normalize_graph_token(getattr(node, "canonical_label", None)),
        }
        if token
    ]


_EXPOSURE_NODE_ALIASES = {
    "bap": "PAH",
    "benzapyrene": "PAH",
    "benzoapyrene": "PAH",
    "dmba": "PAH",
    "phip": "HCA",
    "meiqx": "HCA",
    "4abp": "AromaticAmines",
    "4aminobiphenyl": "AromaticAmines",
    "benzidine": "AromaticAmines",
    "nnk": "TobaccoNitrosamines",
    "ndma": "DietaryNitroso",
    "ndea": "DietaryNitroso",
    "afb1": "AflatoxinB1",
    "aflatoxinb1": "AflatoxinB1",
    "e2": "EstrogenMetabolites",
    "estradiol": "EstrogenMetabolites",
    "17betaestradiol": "EstrogenMetabolites",
    "benzene": "Benzene",
    "vinylchloride": "VinylChloride",
    "formaldehyde": "Aldehyde",
    "acetaldehyde": "Aldehyde",
    "tcdd": "Dioxins_AhR",
    "pcb126": "Dioxins_AhR",
    "tce": "ChlorinatedSolvents",
    "trichloroethylene": "ChlorinatedSolvents",
    "pce": "ChlorinatedSolvents",
    "perchloroethylene": "ChlorinatedSolvents",
    "arsenic": "HeavyMetals",
    "cadmium": "HeavyMetals",
    "chromium": "HeavyMetals",
    "chromiumvi": "HeavyMetals",
    "nickel": "HeavyMetals",
    "beryllium": "HeavyMetals",
    "lead": "HeavyMetals",
    "pb": "HeavyMetals",
    "mercury": "HeavyMetals",
    "hg": "HeavyMetals",
    "cobalt": "HeavyMetals",
    "co": "HeavyMetals",
    "antimony": "HeavyMetals",
    "sb": "HeavyMetals",
}

_EXPOSURE_GROUP_ALIASES = {
    "pah": "PAH",
    "hca": "HCA",
    "aromaticamine": "AromaticAmines",
    "mycotoxin": "AflatoxinB1",
    "estrogen": "EstrogenMetabolites",
    "aldehyde": "Aldehyde",
    "dioxin": "Dioxins_AhR",
    "pcb": "Dioxins_AhR",
    "chlorinatedsolvent": "ChlorinatedSolvents",
    "heavymetal": "HeavyMetals",
}

_KINETIC_NODE_ALIASES = {
    "bap": "PAH",
    "benzapyrene": "PAH",
    "benzoapyrene": "PAH",
    "dmba": "PAH",
    "phip": "HCA",
    "meiqx": "HCA",
    "nnk": "Nitrosamine",
    "ndma": "NDMA",
    "afb1": "Aflatoxin",
    "aflatoxinb1": "Aflatoxin",
    "formaldehyde": "Aldehyde",
    "acetaldehyde": "Aldehyde",
    "benzene": "Benzene",
    "tcdd": "Dioxin",
    "pcb126": "Dioxin",
    "tce": "ChlorinatedSolvent",
    "trichloroethylene": "ChlorinatedSolvent",
    "pce": "ChlorinatedSolvent",
    "perchloroethylene": "ChlorinatedSolvent",
    "arsenic": "HeavyMetal",
    "cadmium": "HeavyMetal",
    "chromium": "HeavyMetal",
    "chromiumvi": "HeavyMetal",
    "nickel": "HeavyMetal",
    "beryllium": "HeavyMetal",
    "lead": "HeavyMetal",
    "pb": "HeavyMetal",
    "mercury": "HeavyMetal",
    "hg": "HeavyMetal",
    "cobalt": "HeavyMetal",
    "co": "HeavyMetal",
    "antimony": "HeavyMetal",
    "sb": "HeavyMetal",
}

_KINETIC_GROUP_ALIASES = {
    "pah": "PAH",
    "hca": "HCA",
    "mycotoxin": "Aflatoxin",
    "aldehyde": "Aldehyde",
    "dioxin": "Dioxin",
    "pcb": "Dioxin",
    "chlorinatedsolvent": "ChlorinatedSolvent",
    "heavymetal": "HeavyMetal",
}


def _resolve_exposure_class_for_node(node: Any) -> str | None:
    """Map a graph carcinogen node onto an exposure database class."""
    tokens = _graph_node_tokens(node)
    for token in tokens:
        if token in _EXPOSURE_NODE_ALIASES:
            return _EXPOSURE_NODE_ALIASES[token]
    group_token = _normalize_graph_token(getattr(node, "group", None))
    return _EXPOSURE_GROUP_ALIASES.get(group_token)


def _resolve_kinetic_class_for_node(node: Any) -> str | None:
    """Map a graph carcinogen node onto a flux-engine kinetic class."""
    tokens = _graph_node_tokens(node)
    for token in tokens:
        if token in _KINETIC_NODE_ALIASES:
            return _KINETIC_NODE_ALIASES[token]
    group_token = _normalize_graph_token(getattr(node, "group", None))
    return _KINETIC_GROUP_ALIASES.get(group_token)


def _match_kinetics_entry(
    pathway: dict[str, Any],
    enzyme_id: str,
) -> tuple[str, dict[str, Any]] | None:
    """Resolve the best kinetic-parameter entry for a graph enzyme id."""
    if enzyme_id in pathway:
        return enzyme_id, pathway[enzyme_id]

    normalized_enzyme = _normalize_graph_token(enzyme_id)
    for key, params in pathway.items():
        if _normalize_graph_token(key) == normalized_enzyme:
            return key, params

    for key, params in pathway.items():
        base = _normalize_graph_token(key.split("_")[0])
        if base == normalized_enzyme:
            return key, params

    return None


def _iter_kinetic_pathways_for_edge(
    pathways: dict[str, dict[str, Any]],
    *,
    preferred_section: str | None,
) -> list[tuple[str, dict[str, Any]]]:
    """Return kinetic pathway sections in edge-matching priority order."""
    ordered: list[tuple[str, dict[str, Any]]] = []
    if preferred_section is not None and preferred_section in pathways:
        ordered.append((preferred_section, pathways[preferred_section]))
    for section_name, section in pathways.items():
        if section_name != preferred_section:
            ordered.append((section_name, section))
    return ordered


def _match_kinetics_for_edge(
    kinetic_class: str,
    kinetic_payload: dict[str, Any],
    edge: Any,
    carcinogen_node: Any | None,
    *,
    preferred_section: str | None,
) -> tuple[str, str, dict[str, Any]] | None:
    """Resolve the best kinetic parameter block for an enriched graph edge."""
    pathways = kinetic_payload.get("pathways", {})
    carcinogen_id = getattr(carcinogen_node, "id", None)
    touches_carcinogen = carcinogen_id in {edge.source, edge.target}

    for section_name, pathway in _iter_kinetic_pathways_for_edge(
        pathways,
        preferred_section=preferred_section,
    ):
        if not isinstance(pathway, dict):
            continue
        # Dioxin kinetics are receptor-binding metadata, so only attach them
        # to the direct carcinogen -> receptor edge rather than every downstream
        # AHR signaling edge that shares the same protein symbol.
        if kinetic_class == "Dioxin" and section_name == "AhR_signaling" and not touches_carcinogen:
            continue

        enzyme_key = _match_kinetics_entry(pathway, edge.source)
        if enzyme_key is None:
            enzyme_key = _match_kinetics_entry(pathway, edge.target)
        if enzyme_key is not None:
            matched_key, params = enzyme_key
            return section_name, matched_key, params

    return None


def enrich_knowledge_graph(
    graph: Any,
    *,
    attach_tissue_weights: bool = True,
    attach_kinetic_params: bool = True,
    attach_exposure_data: bool = True,
) -> KnowledgeGraph:
    """Return a new KnowledgeGraph with quantitative annotation fields attached."""
    from .exposure_engine import _get_exposure_db
    from .flux_engine import _load_kinetic_params
    from .models import EdgeType, KnowledgeGraph, NodeType
    from .tissue_subgraphs import get_available_gtex_tissues, get_tissue_weights

    if isinstance(graph, KnowledgeGraph):
        base_graph = graph
    else:
        base_graph = KnowledgeGraph(**graph)

    node_lookup = {node.id: node for node in base_graph.nodes}

    tissue_weights_by_gene: dict[str, dict[str, float]] = {}
    if attach_tissue_weights:
        tissue_names = get_available_gtex_tissues()
        per_tissue = {tissue: get_tissue_weights(tissue) for tissue in tissue_names}
        for node in base_graph.nodes:
            if node.type != NodeType.ENZYME:
                continue
            weights: dict[str, float] = {}
            for tissue_name, values in per_tissue.items():
                if node.id in values:
                    weights[tissue_name] = round(float(values[node.id]), 6)
            if weights:
                tissue_weights_by_gene[node.id] = weights

    exposure_classes = _get_exposure_db()["carcinogen_classes"] if attach_exposure_data else {}
    kinetic_classes = _load_kinetic_params()["carcinogen_classes"] if attach_kinetic_params else {}

    enriched_nodes = []
    for node in base_graph.nodes:
        updates: dict[str, Any] = {}
        if attach_tissue_weights and node.id in tissue_weights_by_gene:
            merged_weights = dict(node.tissue_weights or {})
            merged_weights.update(tissue_weights_by_gene[node.id])
            updates["tissue_weights"] = merged_weights

        if attach_exposure_data and node.type == NodeType.CARCINOGEN:
            exposure_class = _resolve_exposure_class_for_node(node)
            if exposure_class and exposure_class in exposure_classes:
                merged_scenarios = dict(node.exposure_scenarios or {})
                merged_scenarios.update(deepcopy(exposure_classes[exposure_class].get("exposure_scenarios", {})))
                updates["exposure_scenarios"] = merged_scenarios

        enriched_nodes.append(node.model_copy(update=updates) if updates else node.model_copy(deep=True))

    enriched_edges = []
    for edge in base_graph.edges:
        edge_updates: dict[str, Any] = {}
        if attach_kinetic_params and edge.type in {EdgeType.ACTIVATES, EdgeType.DETOXIFIES}:
            carcinogen_node = None
            if edge.carcinogen:
                carcinogen_node = node_lookup.get(edge.carcinogen)
            if carcinogen_node is None:
                source_node = node_lookup.get(edge.source)
                target_node = node_lookup.get(edge.target)
                if source_node is not None and source_node.type == NodeType.CARCINOGEN:
                    carcinogen_node = source_node
                elif target_node is not None and target_node.type == NodeType.CARCINOGEN:
                    carcinogen_node = target_node

            kinetic_class = _resolve_kinetic_class_for_node(carcinogen_node) if carcinogen_node is not None else None
            if kinetic_class and kinetic_class in kinetic_classes:
                edge_role = "activation" if edge.type == EdgeType.ACTIVATES else "detoxification"
                kinetic_match = _match_kinetics_for_edge(
                    kinetic_class,
                    kinetic_classes[kinetic_class],
                    edge,
                    carcinogen_node,
                    preferred_section=edge_role,
                )
                if kinetic_match is not None:
                    section_name, matched_key, params = kinetic_match
                    merged_kinetics = dict(edge.kinetics or {})
                    merged_kinetics[kinetic_class] = {
                        "pathway_type": section_name,
                        "edge_role": edge_role,
                        "matched_enzyme_key": matched_key,
                        **deepcopy(params),
                    }
                    edge_updates["kinetics"] = merged_kinetics

        enriched_edges.append(
            edge.model_copy(update=edge_updates) if edge_updates else edge.model_copy(deep=True)
        )

    return KnowledgeGraph(nodes=enriched_nodes, edges=enriched_edges)


def patient_risk_query(
    genotypes: dict[str, str],
    tissue: str = "Liver",
    lifestyle: dict[str, Any] | None = None,
    exposure_scenario: str = "general_population",
    exposure_answers: dict[str, Any] | None = None,
    biomarker_measurements: dict[str, float] | None = None,
    *,
    include_interactions: bool = True,
    include_tissue_report: bool = True,
    chip_status: Any = None,
    methylation_status: Any = None,
) -> PatientRiskProfile:
    """Run the full patient-level risk pipeline with partial-failure handling.

    ``chip_status`` may be a :class:`~ExposoGraph.chip_modifier.CHIPStatus`
    instance, a list of CHIP-call dicts (``{"gene": "DNMT3A", "vaf": 0.08}``),
    or ``None`` to skip the somatic modifier layer.

    ``methylation_status`` may be a
    :class:`~ExposoGraph.epigenomic_modifier.MethylationStatus` instance, a
    list of methylation-call dicts
    (``{"gene": "MGMT", "beta_value": 0.72, "status": "hyper"}``), or ``None``
    to skip the epigenomic modifier layer.
    """
    from .chip_modifier import (
        CHIPStatus,
        apply_chip_to_risk_map,
        compute_chip_modifier,
        derive_chip_status,
        parse_chip_calls,
    )
    from .epigenomic_modifier import (
        MethylationStatus,
        apply_methylation_to_risk_map,
        compute_epigenetic_modifier,
        derive_methylation_status,
        parse_methylation_calls,
    )
    from .flux_engine import compute_full_profile
    from .interaction_engine import compute_interaction_matrix, identify_critical_interactions
    from .tissue_subgraphs import generate_tissue_report, get_tissue_weights

    pipeline_warnings: list[str] = []
    pipeline_errors: dict[str, str] = {}
    flux_profile: FullProfileResult | None = None
    exposure_profile: ExposureProfileResult | None = None
    interactions: InteractionMatrixResult | None = None
    critical_warnings: list[CriticalInteraction] = []
    tissue_report: str | None = None
    tissue_weight_count: int | None = None
    top_tissue_genes: list[tuple[str, float]] = []
    biomarker_exposure_profile: dict[str, float] = {}
    biomarker_dose_estimates: list[dict[str, Any]] = []

    normalized_lifestyle = _normalize_lifestyle(lifestyle)
    if exposure_answers:
        inferred = _infer_lifestyle_from_answers(exposure_answers)
        for key, value in inferred.items():
            if value:
                normalized_lifestyle[key] = True

    if exposure_answers is None:
        if exposure_scenario == "general_population" and normalized_lifestyle != _DEFAULT_LIFESTYLE:
            exposure_scenario = _derive_exposure_scenario(normalized_lifestyle)
        exposure_answers = _build_exposure_answers(exposure_scenario, normalized_lifestyle)
    elif exposure_scenario == "general_population" and normalized_lifestyle != _DEFAULT_LIFESTYLE:
        exposure_scenario = _derive_exposure_scenario(normalized_lifestyle)

    if include_tissue_report:
        try:
            tissue_weights = get_tissue_weights(tissue)
            tissue_weight_count = len(tissue_weights)
            top_tissue_genes = _top_tissue_genes(tissue_weights)
            tissue_report = generate_tissue_report(tissue)
        except Exception as exc:
            pipeline_errors["tissue"] = f"{type(exc).__name__}: {exc}"

    try:
        (
            biomarker_exposure_profile,
            biomarker_dose_estimates,
            biomarker_warnings,
        ) = _build_flux_exposure_profile_from_biomarkers(biomarker_measurements)
        pipeline_warnings.extend(biomarker_warnings)
    except Exception as exc:
        pipeline_errors["biomarkers"] = f"{type(exc).__name__}: {exc}"

    try:
        flux_profile = compute_full_profile(
            genotypes,
            tissue,
            exposure_profile=biomarker_exposure_profile or None,
            lifestyle=normalized_lifestyle,
        )
    except Exception as exc:
        pipeline_errors["flux"] = f"{type(exc).__name__}: {exc}"

    try:
        exposure_profile, exposure_warnings = _build_exposure_profile(
            genotypes,
            tissue,
            exposure_answers,
            flux_profile,
            biomarker_measurements,
        )
        pipeline_warnings.extend(exposure_warnings)
    except Exception as exc:
        pipeline_errors["exposure"] = f"{type(exc).__name__}: {exc}"

    if include_interactions:
        try:
            interactions = compute_interaction_matrix(
                _build_interaction_exposure_profile(
                    exposure_scenario,
                    normalized_lifestyle,
                    exposure_answers,
                ),
                genotypes=genotypes,
                tissue=tissue,
                lifestyle=normalized_lifestyle,
            )
        except Exception as exc:
            pipeline_errors["interactions"] = f"{type(exc).__name__}: {exc}"

    try:
        critical_warnings = identify_critical_interactions(genotypes)
    except Exception as exc:
        pipeline_errors["critical_warnings"] = f"{type(exc).__name__}: {exc}"

    base_interaction_risks: dict[str, float] | None = None
    if interactions is not None and interactions.interaction_adjusted_risks:
        base_interaction_risks = dict(interactions.interaction_adjusted_risks)

    resolved_chip_status = None
    chip_effect = None
    chip_adjusted_risks: dict[str, float] = {}
    if chip_status is not None:
        try:
            if isinstance(chip_status, CHIPStatus):
                resolved_chip_status = chip_status
            elif isinstance(chip_status, list):
                resolved_chip_status = parse_chip_calls(chip_status)
            elif isinstance(chip_status, dict):
                if "mutations" in chip_status or "is_present" in chip_status:
                    resolved_chip_status = derive_chip_status(
                        chip_status.get("mutations"),
                        is_present=chip_status.get("is_present"),
                    )
                else:
                    resolved_chip_status = parse_chip_calls([chip_status])
            else:
                raise TypeError(
                    "chip_status must be a CHIPStatus, a list of CHIP-call dicts, or None"
                )
            chip_effect = compute_chip_modifier(resolved_chip_status)
            if base_interaction_risks is not None:
                chip_adjusted_risks = apply_chip_to_risk_map(
                    dict(base_interaction_risks),
                    chip_effect,
                )
        except Exception as exc:
            pipeline_errors["chip"] = f"{type(exc).__name__}: {exc}"

    resolved_methylation_status = None
    methylation_effect = None
    methylation_adjusted_risks: dict[str, float] = {}
    if methylation_status is not None:
        try:
            if isinstance(methylation_status, MethylationStatus):
                resolved_methylation_status = methylation_status
            elif isinstance(methylation_status, list):
                resolved_methylation_status = parse_methylation_calls(methylation_status)
            elif isinstance(methylation_status, dict):
                if "marks" in methylation_status or "present" in methylation_status:
                    resolved_methylation_status = derive_methylation_status(
                        methylation_status.get("marks"),
                        is_present=methylation_status.get("present"),
                    )
                else:
                    resolved_methylation_status = parse_methylation_calls([methylation_status])
            else:
                raise TypeError(
                    "methylation_status must be a MethylationStatus, a list of "
                    "methylation-call dicts, or None"
                )
            methylation_effect = compute_epigenetic_modifier(resolved_methylation_status)
            if base_interaction_risks is not None:
                methylation_adjusted_risks = apply_methylation_to_risk_map(
                    dict(base_interaction_risks),
                    methylation_effect,
                )
        except Exception as exc:
            pipeline_errors["methylation"] = f"{type(exc).__name__}: {exc}"

    profile = PatientRiskProfile(
        tissue=tissue,
        genotypes=dict(genotypes),
        lifestyle=normalized_lifestyle,
        exposure_scenario=exposure_scenario,
        exposure_answers=dict(exposure_answers),
        biomarker_measurements=dict(biomarker_measurements or {}),
        biomarker_dose_estimates=biomarker_dose_estimates,
        tissue_report=tissue_report,
        tissue_weight_count=tissue_weight_count,
        top_tissue_genes=top_tissue_genes,
        flux_profile=flux_profile,
        flux_class_evidence=_build_flux_class_evidence(flux_profile),
        exposure_profile=exposure_profile,
        interactions=interactions,
        critical_warnings=critical_warnings,
        pipeline_warnings=pipeline_warnings,
        pipeline_errors=pipeline_errors,
        chip_status=resolved_chip_status,
        chip_effect=chip_effect,
        chip_adjusted_risks=chip_adjusted_risks,
        methylation_status=resolved_methylation_status,
        methylation_effect=methylation_effect,
        methylation_adjusted_risks=methylation_adjusted_risks,
        biological_output_integration=_build_biological_output_integration(interactions),
    )
    profile.summary = summarize_risk_profile(profile)
    return profile


def _build_biological_output_integration(interactions: InteractionMatrixResult | None) -> dict[str, Any]:
    if interactions is None:
        return {}
    substrate_outputs: dict[str, Any] = {}
    for enzyme, enzyme_result in interactions.competitive_effects.items():
        for substrate, flux in enzyme_result.substrates.items():
            if flux.biological_output is not None:
                substrate_outputs[f"{enzyme}:{substrate}"] = deepcopy(flux.biological_output)
    return {
        "substrate_outputs": substrate_outputs,
        "mechanism_attribution": deepcopy(interactions.mechanism_attribution),
        "source": "interaction_engine.live_biological_output",
    }



def summarize_risk_profile(profile: PatientRiskProfile) -> str:
    """Generate a plain-language summary from a patient risk profile."""
    lines = [
        f"Patient risk assessment for {profile.tissue}",
        f"Exposure scenario: {profile.exposure_scenario}",
    ]

    if profile.top_tissue_genes:
        gene_text = ", ".join(f"{gene} {weight:.2f}" for gene, weight in profile.top_tissue_genes[:3])
        lines.append(f"Top tissue-weighted genes: {gene_text}")

    if profile.flux_profile is not None:
        high_flux = profile.flux_profile.elevated_or_high_risk_classes
        moderate_flux = profile.flux_profile.moderate_risk_classes
        if high_flux:
            lines.append(f"Elevated/high flux classes: {', '.join(high_flux)}")
        elif moderate_flux:
            lines.append(f"Moderate flux classes: {', '.join(moderate_flux)}")

        flagged_flux = set(high_flux or moderate_flux)
        if flagged_flux and profile.flux_class_evidence:
            evidence = [
                item
                for item in profile.flux_class_evidence
                if item.carcinogen_class in flagged_flux
            ]
            measured = [
                item.carcinogen_class
                for item in evidence
                if item.model_kind == "measured_kinetics"
            ]
            proxy = [
                f"{item.carcinogen_class} ({_format_flux_model_kind(item.model_kind)})"
                for item in evidence
                if item.model_kind != "measured_kinetics"
            ]
            if measured:
                lines.append(f"Measured-kinetics flux support: {', '.join(measured)}")
            if proxy:
                lines.append(f"Proxy-backed flux support: {', '.join(proxy)}")

    if profile.biomarker_dose_estimates:
        biomarkers = ", ".join(
            f"{item['biomarker']} -> {item['tissue_conc_uM']:.3g} uM"
            for item in profile.biomarker_dose_estimates[:3]
        )
        lines.append(f"Biomarker-derived substrate inputs: {biomarkers}")

    if profile.exposure_profile is not None:
        exposure_hits = profile.exposure_profile.high_risk_classes or profile.exposure_profile.elevated_risk_classes
        if exposure_hits:
            lines.append(f"Exposure-prioritized classes: {', '.join(exposure_hits)}")
        if profile.exposure_profile.max_risk_class and profile.exposure_profile.max_risk_score is not None:
            lines.append(
                f"Highest exposure-weighted class: {profile.exposure_profile.max_risk_class} "
                f"({profile.exposure_profile.max_risk_score:.2f})"
            )

    if profile.interactions is not None:
        lines.append(
            f"Interaction factor: {profile.interactions.interaction_factor:.2f}x; "
            f"total interaction risk {profile.interactions.total_interaction_risk:.2f}"
        )
        gsh = profile.interactions.gsh_status
        if gsh.tipping_point_reached:
            lines.append(
                f"GSH tipping point reached at {gsh.fraction_normal * 100:.0f}% of normal."
            )
        elif gsh.fraction_normal < 0.5:
            lines.append(f"GSH depleted to {gsh.fraction_normal * 100:.0f}% of normal.")

    if profile.critical_warnings:
        lines.append(
            "Critical interactions: "
            + "; ".join(f"[{item.severity}] {item.interaction}" for item in profile.critical_warnings[:3])
        )

    if profile.pipeline_warnings:
        lines.append("Pipeline warnings: " + "; ".join(profile.pipeline_warnings))
    if profile.pipeline_errors:
        lines.append(
            "Pipeline errors: "
            + "; ".join(f"{step}={message}" for step, message in sorted(profile.pipeline_errors.items()))
        )

    return "\n".join(lines)


def compare_patient_profiles(
    left: PatientRiskProfile,
    right: PatientRiskProfile,
    *,
    left_label: str = "Profile A",
    right_label: str = "Profile B",
) -> PatientProfileComparison:
    """Compare two patient risk profiles across interaction burden and high-risk classes."""
    left_risk = left.interactions.total_interaction_risk if left.interactions is not None else None
    right_risk = right.interactions.total_interaction_risk if right.interactions is not None else None
    left_factor = left.interactions.interaction_factor if left.interactions is not None else None
    right_factor = right.interactions.interaction_factor if right.interactions is not None else None

    risk_delta = None
    if left_risk is not None and right_risk is not None:
        risk_delta = round(left_risk - right_risk, 4)

    factor_delta = None
    if left_factor is not None and right_factor is not None:
        factor_delta = round(left_factor - right_factor, 4)

    left_classes = _collect_high_risk_classes(left)
    right_classes = _collect_high_risk_classes(right)
    shared = sorted(set(left_classes) & set(right_classes))
    left_only = sorted(set(left_classes) - set(right_classes))
    right_only = sorted(set(right_classes) - set(left_classes))

    more_concerning: str | None = None
    if left_risk is not None and right_risk is not None and abs(left_risk - right_risk) > 1e-9:
        more_concerning = left_label if left_risk > right_risk else right_label
    elif left_factor is not None and right_factor is not None and abs(left_factor - right_factor) > 1e-9:
        more_concerning = left_label if left_factor > right_factor else right_label
    elif len(left_classes) != len(right_classes):
        more_concerning = left_label if len(left_classes) > len(right_classes) else right_label

    summary_parts = [
        f"{left_label} high-risk classes: {', '.join(left_classes) if left_classes else 'none'}",
        f"{right_label} high-risk classes: {', '.join(right_classes) if right_classes else 'none'}",
    ]
    if more_concerning is not None:
        summary_parts.append(f"Overall more concerning profile: {more_concerning}")

    return PatientProfileComparison(
        left_label=left_label,
        right_label=right_label,
        tissue=left.tissue if left.tissue == right.tissue else f"{left.tissue} vs {right.tissue}",
        left_total_interaction_risk=left_risk,
        right_total_interaction_risk=right_risk,
        total_interaction_risk_delta=risk_delta,
        left_interaction_factor=left_factor,
        right_interaction_factor=right_factor,
        interaction_factor_delta=factor_delta,
        left_high_risk_classes=left_classes,
        right_high_risk_classes=right_classes,
        shared_high_risk_classes=shared,
        left_only_high_risk_classes=left_only,
        right_only_high_risk_classes=right_only,
        more_concerning_profile=more_concerning,
        summary=" | ".join(summary_parts),
    )
