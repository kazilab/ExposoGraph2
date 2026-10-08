"""Exposure integration engine for Gene x Environment risk scoring.

Scenarios are read at call time from exposure_database_revised.json.
Callers pass a carcinogen_group_id and a base_scenario. A duplicated
base_scenario returns a list. The file is not loaded onto graph edges.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, TypeAlias, cast


class ExposureRiskCategory(str, Enum):
    """Risk categories for exposure-weighted risk scores."""

    LOW = "Low Risk"
    AVERAGE = "Population Average"
    ELEVATED = "Elevated Risk"
    HIGH = "High Risk"


JsonDict: TypeAlias = dict[str, Any]


@dataclass
class ExposureScenarioInfo:
    """One exposure scenario for a carcinogen group."""

    scenario_id: str
    base_scenario: str
    carcinogen_group_id: str
    label: str
    multiplier: float
    tissue_conc_uM: float
    source: str
    tier: int
    tier_label: str
    exposure_source: str | None = None
    member_carcinogen_ids: list[str] = field(default_factory=list)


@dataclass
class ExposureWeightedRisk:
    """One result of compute_exposure_weighted_risk."""

    carcinogen_group_id: str
    scenario_id: str
    base_scenario: str
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
    interpretation: str
    exposure_source: str | None = None
    member_carcinogen_ids: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)


@dataclass
class LifetimeCancerRisk:
    """One member result of compute_lifetime_cancer_risk."""

    carcinogen_group_id: str
    member_carcinogen_id: str
    tissue: str
    daily_dose_mg_kg: float
    duration_years: float
    slope_factor: float
    genotype_modifier: float
    lecr: float
    exceeds_action_threshold: bool
    exceeds_de_minimis: bool


@dataclass
class ScenarioComparison:
    """One paired comparison of two scenario results."""

    carcinogen_group_id: str
    tissue: str
    genotypes: dict[str, str]
    scenario_a: ExposureWeightedRisk
    scenario_b: ExposureWeightedRisk
    fold_change_b_vs_a: float
    absolute_delta: float
    interpretation: str


_EXPOSURE_DB_FILE = Path(__file__).parent / "data" / "exposure_database_revised.json"
_EXPOSURE_CACHE: JsonDict | None = None


def _load_exposure_db() -> JsonDict:
    with open(_EXPOSURE_DB_FILE, encoding="utf-8") as handle:
        return cast(JsonDict, json.load(handle))


def _get_exposure_db() -> JsonDict:
    global _EXPOSURE_CACHE
    if _EXPOSURE_CACHE is None:
        _EXPOSURE_CACHE = _load_exposure_db()
    return _EXPOSURE_CACHE


def get_database() -> JsonDict:
    """Return a defensive copy of the runtime exposure database."""
    from copy import deepcopy

    return deepcopy(_get_exposure_db())


def _scenarios_for_group(carcinogen_group_id: str) -> list[tuple[str, JsonDict]]:
    matches = [
        (scenario_id, scenario)
        for scenario_id, scenario in _get_exposure_db()["exposure_scenarios"].items()
        if scenario.get("carcinogen_group_id") == carcinogen_group_id
    ]
    if not matches:
        known = sorted({
            scenario.get("carcinogen_group_id")
            for scenario in _get_exposure_db()["exposure_scenarios"].values()
        })
        raise ValueError(
            f"Unknown carcinogen_group_id '{carcinogen_group_id}'. "
            f"Available: {known}"
        )
    return matches


def _scenarios_for_base(carcinogen_group_id: str, base_scenario: str) -> list[tuple[str, JsonDict]]:
    matches = [
        item
        for item in _scenarios_for_group(carcinogen_group_id)
        if item[1].get("base_scenario") == base_scenario
    ]
    if not matches:
        available = sorted({item[1].get("base_scenario") for item in _scenarios_for_group(carcinogen_group_id)})
        raise ValueError(
            f"Unknown base_scenario '{base_scenario}' for {carcinogen_group_id}. "
            f"Available: {available}"
        )
    return matches


def _numeric(scenario: Mapping[str, Any], *keys: str, default: float = 0.0) -> float:
    for key in keys:
        value = scenario.get(key)
        if value is None:
            continue
        if isinstance(value, dict):
            numbers = [float(item) for item in value.values() if isinstance(item, (int, float))]
            return max(numbers) if numbers else default
        if isinstance(value, (int, float)):
            return float(value)
    return default


def _tissue_factor(group_id: str, exposure_source: str | None, tissue: str) -> float:
    row = _get_exposure_db().get("tissue_factors", {}).get(group_id, {})
    if exposure_source and isinstance(row.get(exposure_source), dict):
        row = row[exposure_source]
    query = tissue.replace(" ", "").replace("-", "").lower()
    for key, factor in row.items():
        if isinstance(factor, (int, float)) and key.lower() == query:
            return float(factor)
    return 1.0


def _genotype_matches(modifier_key: str, geno_val: str) -> bool:
    if "null" in modifier_key and "null" in geno_val:
        return True
    if "heterozygote" in modifier_key and any(token in geno_val for token in ("*1/*2", "hetero", "het")):
        return True
    if "homozygous_variant" in modifier_key and any(token in geno_val for token in ("*2/*2", "homo", "hom")):
        return True
    if "slow" in modifier_key and "slow" in geno_val:
        return True
    if "PM" in modifier_key and any(token in geno_val for token in ("pm", "poor")):
        return True
    if "high_activity" in modifier_key and any(token in geno_val for token in ("high", "ultra", "rapid", "*1f", "*1/*1f")):
        return True
    if "low_activity" in modifier_key and any(token in geno_val for token in ("low", "y113h", "his113")):
        return True
    return False


def compute_flux_ratio(genotypes: dict[str, str], member_carcinogen_ids: list[str]) -> float:
    """Genotype multiplier for one scenario. Applied once per matching modifier."""
    modifiers = _get_exposure_db()["flux_model_integration"]["genotype_modifiers"]
    members = set(member_carcinogen_ids)
    norm = {key.upper(): value.lower() for key, value in genotypes.items()}
    flux = 1.0
    for modifier_key, by_member in modifiers.items():
        gene = modifier_key.partition("_")[0].upper()
        if not _genotype_matches(modifier_key, norm.get(gene, "")):
            continue
        matched = [float(value) for mid, value in by_member.items() if mid in members]
        if matched:
            flux *= matched[0]
    return round(flux, 4)


def _classify_risk(score: float, thresholds: Mapping[str, Any]) -> ExposureRiskCategory:
    if score <= float(thresholds["low"]["max"]):
        return ExposureRiskCategory.LOW
    if score <= float(thresholds["average"]["max"]):
        return ExposureRiskCategory.AVERAGE
    if score <= float(thresholds["elevated"]["max"]):
        return ExposureRiskCategory.ELEVATED
    return ExposureRiskCategory.HIGH


def classify_exposure_risk(score: float) -> ExposureRiskCategory:
    thresholds = _get_exposure_db()["flux_model_integration"]["risk_thresholds"]
    return _classify_risk(score, thresholds)


def _interpret_risk(
    group_id: str,
    genotypes: dict[str, str],
    tissue: str,
    flux_ratio: float,
    exposure_mult: float,
    combined: float,
    category: ExposureRiskCategory,
    sc_label: str,
) -> str:
    geno_parts = [
        f"{gene} {status}"
        for gene, status in genotypes.items()
        if status.lower() not in ("*1/*1", "wt", "wildtype", "normal", "present")
    ]
    geno_str = (", ".join(geno_parts) or "wildtype") + " genotype"
    if combined < 0.5:
        action = "Minimal concern -- exposure level is very low regardless of genotype."
    elif category == ExposureRiskCategory.LOW:
        action = "Genotype provides some protection or exposure is below population average."
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
        f"{category.value} | {group_id} -> {tissue}. "
        f"Flux ratio (genotype): {flux_ratio:.2f}x baseline. "
        f"Exposure scenario: '{sc_label}' ({exposure_mult:.1f}x population average). "
        f"Combined score: {combined:.2f}. {action}"
    )


def _scenario_info(scenario_id: str, scenario: Mapping[str, Any]) -> ExposureScenarioInfo:
    return ExposureScenarioInfo(
        scenario_id=scenario_id,
        base_scenario=str(scenario.get("base_scenario", "")),
        carcinogen_group_id=str(scenario.get("carcinogen_group_id", "")),
        label=str(scenario.get("label", scenario_id)),
        multiplier=_numeric(scenario, "multiplier_vs_baseline", "multiplier_vs_baseline_chromium", default=1.0),
        tissue_conc_uM=_numeric(scenario, "estimated_tissue_conc_uM", "tissue_conc_uM"),
        source=str(scenario.get("source", "")),
        tier=int(scenario.get("exposure_tier", 1)),
        tier_label=str(scenario.get("tier_label", "")),
        exposure_source=scenario.get("exposure_source"),
        member_carcinogen_ids=list(scenario.get("member_carcinogen_ids", [])),
    )


def get_exposure_scenarios(carcinogen_group_id: str) -> list[ExposureScenarioInfo]:
    """Return every scenario for a carcinogen group."""
    return [_scenario_info(sid, scenario) for sid, scenario in _scenarios_for_group(carcinogen_group_id)]


def compute_exposure_weighted_risk(
    carcinogen_group_id: str,
    base_scenario: str,
    genotypes: dict[str, str],
    tissue: str,
    *,
    flux_ratio_override: float | None = None,
) -> list[ExposureWeightedRisk]:
    """Score every scenario with this group id and base scenario name.

    A duplicated base_scenario, such as nitrosamine general_population,
    returns one result per scenario. Callers do not pass an exposure source.
    """
    thresholds = _get_exposure_db()["flux_model_integration"]["risk_thresholds"]
    results: list[ExposureWeightedRisk] = []
    for scenario_id, scenario in _scenarios_for_base(carcinogen_group_id, base_scenario):
        members = list(scenario.get("member_carcinogen_ids", []))
        source = scenario.get("exposure_source")
        tissue_factor = _tissue_factor(carcinogen_group_id, source, tissue)
        flux_ratio = (
            float(flux_ratio_override)
            if flux_ratio_override is not None
            else compute_flux_ratio(genotypes, members)
        )
        exposure_mult = _numeric(
            scenario, "multiplier_vs_baseline", "multiplier_vs_baseline_chromium", default=1.0
        )
        tissue_conc = _numeric(scenario, "estimated_tissue_conc_uM", "tissue_conc_uM")
        combined = round(flux_ratio * exposure_mult * tissue_factor, 4)
        category = _classify_risk(combined, thresholds)
        label = str(scenario.get("label", scenario_id))
        slope_source = ""
        for member in members:
            record = _get_exposure_db().get("risk_coefficients", {}).get(member, {})
            slope = record.get("epa_cancer_slope_factor", {})
            if isinstance(slope, dict) and slope.get("source"):
                slope_source = str(slope["source"])
                break
        results.append(
            ExposureWeightedRisk(
                carcinogen_group_id=carcinogen_group_id,
                scenario_id=scenario_id,
                base_scenario=str(scenario.get("base_scenario", base_scenario)),
                scenario_label=label,
                tissue=tissue,
                flux_ratio=flux_ratio,
                exposure_multiplier=exposure_mult,
                exposure_tier=int(scenario["exposure_tier"]),
                exposure_tier_label=str(scenario.get("tier_label", "")),
                tissue_factor=tissue_factor,
                combined_risk_score=combined,
                risk_category=category,
                tissue_conc_uM=tissue_conc,
                interpretation=_interpret_risk(
                    carcinogen_group_id, genotypes, tissue, flux_ratio,
                    exposure_mult, combined, category, label,
                ),
                exposure_source=source,
                member_carcinogen_ids=members,
                sources=[str(scenario.get("source", "")), slope_source],
            )
        )
    return results


def _slope_factor_per_mg_kg_day(sf_data: Mapping[str, Any]) -> float | None:
    raw_value = sf_data.get("value")
    if not isinstance(raw_value, (int, float)):
        return None
    value = float(raw_value)
    unit = str(sf_data.get("unit", "")).lower().replace("µ", "μ")
    if "per μg/kg" in unit or "per ug/kg" in unit:
        return value * 1000.0
    if "per ng/kg" in unit:
        return value * 1_000_000.0
    return value


def compute_lifetime_cancer_risk(
    carcinogen_group_id: str,
    genotypes: dict[str, str],
    daily_dose_mg_kg: float,
    *,
    duration_years: float = 70.0,
    tissue: str = "Liver",
    flux_ratio_override: float | None = None,
) -> list[LifetimeCancerRisk]:
    """One result per member of the group that has a top-level slope value.

    Nested oral or inhalation objects are not coerced into a slope factor.
    """
    member_ids: list[str] = []
    seen: set[str] = set()
    for _, scenario in _scenarios_for_group(carcinogen_group_id):
        for member in scenario.get("member_carcinogen_ids", []):
            if member not in seen:
                seen.add(member)
                member_ids.append(member)
    duration_fraction = min(duration_years / 70.0, 1.0)
    results: list[LifetimeCancerRisk] = []
    coefficients = _get_exposure_db().get("risk_coefficients", {})
    for member in member_ids:
        record = coefficients.get(member, {})
        slope_data = record.get("epa_cancer_slope_factor", {})
        if not isinstance(slope_data, dict):
            continue
        slope_factor = _slope_factor_per_mg_kg_day(slope_data)
        if slope_factor is None:
            continue
        source = None
        for _, scenario in _scenarios_for_group(carcinogen_group_id):
            if member in scenario.get("member_carcinogen_ids", []):
                source = scenario.get("exposure_source")
                break
        if flux_ratio_override is not None:
            genotype_mod = float(flux_ratio_override)
        else:
            genotype_mod = compute_flux_ratio(genotypes, [member])
            genotype_mod = round(genotype_mod * _tissue_factor(carcinogen_group_id, source, tissue), 4)
        lecr = round(slope_factor * daily_dose_mg_kg * genotype_mod * duration_fraction, 8)
        results.append(
            LifetimeCancerRisk(
                carcinogen_group_id=carcinogen_group_id,
                member_carcinogen_id=member,
                tissue=tissue,
                daily_dose_mg_kg=daily_dose_mg_kg,
                duration_years=duration_years,
                slope_factor=slope_factor,
                genotype_modifier=genotype_mod,
                lecr=lecr,
                exceeds_action_threshold=lecr >= 1e-4,
                exceeds_de_minimis=lecr >= 1e-6,
            )
        )
    return results


def compare_scenarios(
    carcinogen_group_id: str,
    base_scenario_a: str,
    base_scenario_b: str,
    genotypes: dict[str, str],
    tissue: str,
) -> list[ScenarioComparison]:
    """Compare two base scenarios. Pair duplicate results by exposure source."""
    left = compute_exposure_weighted_risk(carcinogen_group_id, base_scenario_a, genotypes, tissue)
    right = compute_exposure_weighted_risk(carcinogen_group_id, base_scenario_b, genotypes, tissue)
    pairs: list[tuple[ExposureWeightedRisk, ExposureWeightedRisk]] = []
    if len(left) == 1 and len(right) == 1:
        pairs = [(left[0], right[0])]
    else:
        right_by_source = {item.exposure_source: item for item in right}
        for item in left:
            match = right_by_source.get(item.exposure_source)
            if match is not None:
                pairs.append((item, match))
        if not pairs:
            raise ValueError(
                f"Could not pair {base_scenario_a} with {base_scenario_b} for {carcinogen_group_id}."
            )
    comparisons: list[ScenarioComparison] = []
    for result_a, result_b in pairs:
        score_a = result_a.combined_risk_score
        score_b = result_b.combined_risk_score
        ratio = (score_b / score_a) if score_a else (float("inf") if score_b else 1.0)
        delta = score_b - score_a
        interpretation = (
            f"Scenario B ('{result_b.scenario_id}') has {ratio:.1f}x the risk of "
            f"Scenario A ('{result_a.scenario_id}') for {carcinogen_group_id} in {tissue}. "
            f"Risk score difference: {delta:+.3f}."
        )
        comparisons.append(
            ScenarioComparison(
                carcinogen_group_id=carcinogen_group_id,
                tissue=tissue,
                genotypes=genotypes,
                scenario_a=result_a,
                scenario_b=result_b,
                fold_change_b_vs_a=round(ratio, 3),
                absolute_delta=round(delta, 4),
                interpretation=interpretation,
            )
        )
    return comparisons


def cli_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="ExposoGraph exposure integration engine")
    parser.add_argument("--group", dest="carcinogen_group_id", default="group_pahs")
    parser.add_argument("--scenario", dest="base_scenario", default="general_population")
    parser.add_argument("--tissue", default="Liver")
    parser.add_argument("--genotypes", default="{}")
    parser.add_argument("--list-scenarios", metavar="GROUP")
    parser.add_argument("--lecr", action="store_true")
    parser.add_argument("--daily-dose", type=float, default=0.001)
    parser.add_argument("--duration", type=float, default=70.0)
    parser.add_argument("--compare", nargs=2, metavar=("BASE_A", "BASE_B"))
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    try:
        genotypes = json.loads(args.genotypes)
    except json.JSONDecodeError as exc:
        print(f"ERROR: Invalid genotypes JSON: {exc}", file=sys.stderr)
        return 1
    if args.list_scenarios:
        payload = [asdict(item) for item in get_exposure_scenarios(args.list_scenarios)]
        print(json.dumps(payload, indent=2))
        return 0
    if args.lecr:
        payload = [
            asdict(item)
            for item in compute_lifetime_cancer_risk(
                args.carcinogen_group_id,
                genotypes,
                args.daily_dose,
                duration_years=args.duration,
                tissue=args.tissue,
            )
        ]
        print(json.dumps(payload, indent=2))
        return 0
    if args.compare:
        payload = [
            asdict(item)
            for item in compare_scenarios(
                args.carcinogen_group_id, args.compare[0], args.compare[1], genotypes, args.tissue
            )
        ]
        print(json.dumps(payload, indent=2, default=str))
        return 0
    payload = [
        asdict(item)
        for item in compute_exposure_weighted_risk(
            args.carcinogen_group_id, args.base_scenario, genotypes, args.tissue
        )
    ]
    print(json.dumps(payload, indent=2, default=str))
    return 0


def main(argv: list[str] | None = None) -> int:
    return cli_main(argv)
