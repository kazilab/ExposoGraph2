"""Toxicokinetic-Toxicodynamic (TK-TD) modeling for metal carcinogens.

This module provides physiologically-based pharmacokinetic (PBPK) and
empirical toxicokinetic models for predicting metal tissue concentrations
and dose-response relationships.

Features:
- 1-compartment and 2-compartment TK models
- Route-specific absorption (inhalation, oral, dermal)
- Target tissue dosimetry for key organs
- Dose-response modeling for oxidative stress endpoints
- Integration with biomarker prediction
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np
from scipy.integrate import odeint


class ExposureRoute(str, Enum):
    """Exposure routes for metal intake."""

    INHALATION = "inhalation"
    ORAL = "oral"
    DERMAL = "dermal"
    INJECTION = "injection"


class Species(str, Enum):
    """Species for TK modeling."""

    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"
    RABBIT = "rabbit"
    DOG = "dog"
    MONKEY = "monkey"


@dataclass
class TKParameters:
    """Toxicokinetic parameters for a metal."""

    metal: str
    species: Species

    # Absorption parameters
    oral_bioavailability: float  # Fraction absorbed (0-1)
    inhalation_absorption: float  # Fraction inhaled that reaches blood
    dermal_absorption: float  # Fraction absorbed through skin

    # Distribution parameters
    volume_of_distribution: float  # L/kg body weight
    tissue_partition_coeffs: dict[str, float]  # Tissue:Blood ratios

    # Elimination parameters
    elimination_rate_constant: float  # 1/day
    clearance: float | None = None  # L/day/kg
    half_life_days: float | None = None

    # Route-specific modifiers
    exposure_specific_params: dict[str, float] = field(default_factory=dict)


@dataclass
class ExposureScenario:
    """Exposure scenario definition."""

    route: ExposureRoute
    dose_per_day: float  # mg/day
    duration_days: int
    body_weight_kg: float

    # Route-specific parameters
    exposure_frequency: str = "daily"  # daily, intermittent, chronic
    particle_size: str | None = None  # For inhalation: PM10, PM2.5, ultrafine
    solubility: str | None = None  # For oral: soluble, insoluble


@dataclass
class TKResults:
    """Results from TK simulation."""

    time_days: np.ndarray
    blood_concentration: np.ndarray  # mg/L
    tissue_concentrations: dict[str, np.ndarray]  # mg/L or mg/kg
    cumulative_dose: np.ndarray  # mg
    auc: float  # Area under curve (mg·day/L)
    steady_state_concentration: float | None = None
    time_to_steady_state: float | None = None


@dataclass
class TDParameters:
    """Toxicodynamic parameters for dose-response modeling."""

    endpoint: str  # e.g., "MDA", "8-OHdG", "cell_death"

    # Hill equation parameters
    ec50: float  # Concentration producing 50% response
    emax: float  # Maximum response
    hill_coefficient: float = 1.0

    # Background/baseline
    baseline_response: float = 0.0

    # Time-dependent parameters
    time_to_response_days: float = 7.0
    recovery_rate: float | None = None  # 1/day


# ── Metal-Specific TK Parameters Database ─────────────────────────────────

# Parameters based on literature data (simplified for modeling)
METAL_TK_DATABASE: dict[tuple[str, Species], TKParameters] = {
    # Arsenic - Human
    ("Arsenic", Species.HUMAN): TKParameters(
        metal="Arsenic",
        species=Species.HUMAN,
        oral_bioavailability=0.95,  # Highly bioavailable
        inhalation_absorption=0.30,
        dermal_absorption=0.03,
        volume_of_distribution=0.12,  # L/kg in this simplified blood-focused model
        tissue_partition_coeffs={
            "liver": 3.0,
            "kidney": 5.0,
            "lung": 2.0,
            "skin": 1.5,
            "hair": 10.0,
        },
        elimination_rate_constant=0.11,  # 1/day (t1/2 ~6.3 days for As)
        clearance=0.0132,  # L/day/kg
        half_life_days=6.3,
    ),
    # Arsenic - Rat
    ("Arsenic", Species.RAT): TKParameters(
        metal="Arsenic",
        species=Species.RAT,
        oral_bioavailability=0.90,
        inhalation_absorption=0.25,
        dermal_absorption=0.02,
        volume_of_distribution=0.30,
        tissue_partition_coeffs={
            "liver": 4.0,
            "kidney": 6.0,
            "lung": 2.5,
        },
        elimination_rate_constant=0.17,
        clearance=0.051,
        half_life_days=4.1,
    ),
    # Arsenic - Mouse
    ("Arsenic", Species.MOUSE): TKParameters(
        metal="Arsenic",
        species=Species.MOUSE,
        oral_bioavailability=0.92,
        inhalation_absorption=0.28,
        dermal_absorption=0.025,
        volume_of_distribution=0.35,
        tissue_partition_coeffs={
            "liver": 4.5,
            "kidney": 7.0,
        },
        elimination_rate_constant=0.23,
        clearance=0.0805,
        half_life_days=3.0,
    ),
    # Cadmium - Human
    ("Cadmium", Species.HUMAN): TKParameters(
        metal="Cadmium",
        species=Species.HUMAN,
        oral_bioavailability=0.05,  # Low bioavailability
        inhalation_absorption=0.50,  # Higher for respirable particles
        dermal_absorption=0.01,
        volume_of_distribution=0.50,
        tissue_partition_coeffs={
            "kidney": 15.0,  # Highly concentrated
            "liver": 8.0,
            "lung": 5.0,
            "bone": 3.0,
        },
        elimination_rate_constant=0.00014,  # Very slow (t1/2 ~17 years)
        clearance=0.00007,
        half_life_days=17 * 365,
    ),
    # Cadmium - Rat
    ("Cadmium", Species.RAT): TKParameters(
        metal="Cadmium",
        species=Species.RAT,
        oral_bioavailability=0.02,
        inhalation_absorption=0.40,
        dermal_absorption=0.005,
        volume_of_distribution=0.40,
        tissue_partition_coeffs={
            "kidney": 20.0,
            "liver": 10.0,
        },
        elimination_rate_constant=0.0004,
        clearance=0.00016,
        half_life_days=4.75 * 365,
    ),
    # Lead - Human
    ("Lead", Species.HUMAN): TKParameters(
        metal="Lead",
        species=Species.HUMAN,
        oral_bioavailability=0.15,  # Adults (children ~50%)
        inhalation_absorption=0.40,
        dermal_absorption=0.01,
        volume_of_distribution=0.30,
        tissue_partition_coeffs={
            "bone": 50.0,  # Major storage site
            "kidney": 5.0,
            "liver": 3.0,
            "brain": 1.0,
        },
        elimination_rate_constant=0.000063,  # t1/2 ~30 years in bone
        clearance=0.000019,
        half_life_days=30 * 365,
    ),
    # Lead - Rat
    ("Lead", Species.RAT): TKParameters(
        metal="Lead",
        species=Species.RAT,
        oral_bioavailability=0.10,
        inhalation_absorption=0.35,
        dermal_absorption=0.008,
        volume_of_distribution=0.25,
        tissue_partition_coeffs={
            "bone": 40.0,
            "kidney": 4.0,
        },
        elimination_rate_constant=0.0005,
        clearance=0.000125,
        half_life_days=3.8 * 365,
    ),
    # Mercury - Human (methylmercury)
    ("Mercury", Species.HUMAN): TKParameters(
        metal="Mercury",
        species=Species.HUMAN,
        oral_bioavailability=0.95,  # Nearly complete absorption
        inhalation_absorption=0.80,  # Elemental Hg vapor
        dermal_absorption=0.05,
        volume_of_distribution=0.50,
        tissue_partition_coeffs={
            "brain": 5.0,  # Crosses BBB
            "kidney": 6.0,
            "liver": 4.0,
            "hair": 100.0,  # Highly accumulated
        },
        elimination_rate_constant=0.014,  # t1/2 ~50-70 days
        clearance=0.007,
        half_life_days=70,
    ),
    # Mercury - Mouse
    ("Mercury", Species.MOUSE): TKParameters(
        metal="Mercury",
        species=Species.MOUSE,
        oral_bioavailability=0.93,
        inhalation_absorption=0.75,
        dermal_absorption=0.04,
        volume_of_distribution=0.60,
        tissue_partition_coeffs={
            "brain": 6.0,
            "kidney": 8.0,
        },
        elimination_rate_constant=0.023,
        clearance=0.0138,
        half_life_days=30,
    ),
    # Chromium - Human
    ("Chromium", Species.HUMAN): TKParameters(
        metal="Chromium",
        species=Species.HUMAN,
        oral_bioavailability=0.02,  # Cr(III) poorly absorbed
        inhalation_absorption=0.15,
        dermal_absorption=0.01,
        volume_of_distribution=0.15,
        tissue_partition_coeffs={
            "liver": 2.0,
            "kidney": 3.0,
            "lung": 2.0,
        },
        elimination_rate_constant=0.08,
        clearance=0.012,
        half_life_days=8.7,
    ),
    # Nickel - Human
    ("Nickel", Species.HUMAN): TKParameters(
        metal="Nickel",
        species=Species.HUMAN,
        oral_bioavailability=0.05,
        inhalation_absorption=0.30,
        dermal_absorption=0.02,
        volume_of_distribution=0.20,
        tissue_partition_coeffs={
            "kidney": 4.0,
            "lung": 3.0,
            "liver": 2.0,
        },
        elimination_rate_constant=0.35,
        clearance=0.07,
        half_life_days=2.0,
    ),
    # Cobalt - Human
    ("Cobalt", Species.HUMAN): TKParameters(
        metal="Cobalt",
        species=Species.HUMAN,
        oral_bioavailability=0.30,
        inhalation_absorption=0.40,
        dermal_absorption=0.01,
        volume_of_distribution=0.25,
        tissue_partition_coeffs={
            "liver": 3.0,
            "kidney": 2.5,
            "heart": 2.0,
            "lung": 2.0,
        },
        elimination_rate_constant=0.35,  # t1/2 ~2 days
        clearance=0.0875,
        half_life_days=2.0,
    ),
    # Antimony - Human
    ("Antimony", Species.HUMAN): TKParameters(
        metal="Antimony",
        species=Species.HUMAN,
        oral_bioavailability=0.15,
        inhalation_absorption=0.20,
        dermal_absorption=0.01,
        volume_of_distribution=0.20,
        tissue_partition_coeffs={
            "liver": 2.0,
            "lung": 3.0,
        },
        elimination_rate_constant=0.35,
        clearance=0.07,
        half_life_days=2.0,
    ),
}


# TD parameters for oxidative stress endpoints
TD_DATABASE: dict[tuple[str, str], TDParameters] = {
    # MDA (lipid peroxidation)
    ("Arsenic", "MDA"): TDParameters(
        endpoint="MDA",
        ec50=0.5,  # mg/L blood concentration
        emax=3.0,  # 3x baseline
        hill_coefficient=1.5,
        baseline_response=1.5,  # μmol/L baseline
        time_to_response_days=14.0,
    ),
    ("Cadmium", "MDA"): TDParameters(
        endpoint="MDA",
        ec50=0.2,
        emax=2.5,
        hill_coefficient=1.2,
        baseline_response=1.5,
        time_to_response_days=30.0,
    ),
    ("Lead", "MDA"): TDParameters(
        endpoint="MDA",
        ec50=0.3,
        emax=2.0,
        hill_coefficient=1.0,
        baseline_response=1.5,
        time_to_response_days=60.0,
    ),
    # 8-OHdG (DNA oxidation)
    ("Arsenic", "8-OHdG"): TDParameters(
        endpoint="8-OHdG",
        ec50=0.3,
        emax=2.5,
        hill_coefficient=2.0,
        baseline_response=25.0,  # ng/mg creatinine
        time_to_response_days=21.0,
    ),
    ("Cadmium", "8-OHdG"): TDParameters(
        endpoint="8-OHdG",
        ec50=0.15,
        emax=2.0,
        hill_coefficient=1.8,
        baseline_response=25.0,
        time_to_response_days=45.0,
    ),
    ("Chromium", "8-OHdG"): TDParameters(
        endpoint="8-OHdG",
        ec50=0.1,
        emax=3.0,
        hill_coefficient=2.5,
        baseline_response=25.0,
        time_to_response_days=7.0,
    ),
    # SOD (antioxidant enzyme - typically decreases)
    ("Arsenic", "SOD"): TDParameters(
        endpoint="SOD",
        ec50=0.4,
        emax=0.5,  # 50% of baseline (decrease)
        hill_coefficient=1.5,
        baseline_response=100.0,  # U/mg protein
        time_to_response_days=30.0,
    ),
    ("Cadmium", "SOD"): TDParameters(
        endpoint="SOD",
        ec50=0.3,
        emax=0.6,
        hill_coefficient=1.0,
        baseline_response=100.0,
        time_to_response_days=45.0,
    ),
}


def _trapezoid_integral(y: np.ndarray, x: np.ndarray) -> float:
    """Return the trapezoidal integral across NumPy versions."""
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is not None:
        return float(trapezoid(y, x))
    return float(np.trapz(y, x))


# ── TK Model Functions ────────────────────────────────────────────────────


def get_tk_params(metal: str, species: Species) -> TKParameters | None:
    """Get TK parameters for a metal-species combination.

    Args:
        metal: Metal name
        species: Species enum

    Returns:
        TKParameters if available, None otherwise
    """
    return METAL_TK_DATABASE.get((metal, species))


def one_compartment_model(
    params: TKParameters,
    scenario: ExposureScenario,
    time_points: np.ndarray | None = None,
) -> TKResults:
    """Simulate 1-compartment TK model.

    Args:
        params: TK parameters
        scenario: Exposure scenario
        time_points: Optional custom time points (days)

    Returns:
        TKResults with concentration-time profile
    """
    if time_points is None:
        time_points = np.linspace(0, scenario.duration_days, 1000)

    # Get absorption fraction based on route
    if scenario.route == ExposureRoute.ORAL:
        absorption = params.oral_bioavailability
    elif scenario.route == ExposureRoute.INHALATION:
        absorption = params.inhalation_absorption
    elif scenario.route == ExposureRoute.DERMAL:
        absorption = params.dermal_absorption
    else:
        absorption = 1.0

    # Calculate dose rate (mg/day)
    dose_rate = scenario.dose_per_day * absorption

    # Differential equation: dC/dt = (Dose/Vd) - k*C
    # where C = concentration, k = elimination rate constant
    Vd = params.volume_of_distribution * scenario.body_weight_kg  # Total Vd in L
    k = params.elimination_rate_constant

    def model(y: list[float], t: float) -> list[float]:
        # y is [C] - concentration as list
        C = y[0]
        dCdt = dose_rate / Vd - k * C
        return [dCdt]

    # Solve ODE
    C0 = [0.0]  # Initial concentration as list
    solution = odeint(model, C0, time_points)
    blood_conc = solution[:, 0]

    # Calculate tissue concentrations
    tissue_concs = {}
    for tissue, partition in params.tissue_partition_coeffs.items():
        tissue_concs[tissue] = blood_conc * partition

    # Calculate cumulative dose
    cumulative = np.cumsum(np.full_like(time_points, dose_rate)) * np.diff(
        np.concatenate([[0], time_points])
    )

    # Calculate AUC (trapezoidal rule)
    auc = _trapezoid_integral(blood_conc, time_points)

    # Estimate steady state
    ss_conc = dose_rate / (k * Vd) if k > 0 else np.inf

    return TKResults(
        time_days=time_points,
        blood_concentration=blood_conc,
        tissue_concentrations=tissue_concs,
        cumulative_dose=cumulative,
        auc=auc,
        steady_state_concentration=ss_conc,
        time_to_steady_state=4.0 / k if k > 0 else None,  # ~4 half-lives
    )


def two_compartment_model(
    params: TKParameters,
    scenario: ExposureScenario,
    k12: float = 0.5,  # Central to peripheral
    k21: float = 0.3,  # Peripheral to central
    time_points: np.ndarray | None = None,
) -> TKResults:
    """Simulate 2-compartment TK model.

    Args:
        params: TK parameters
        scenario: Exposure scenario
        k12: Transfer rate central to peripheral (1/day)
        k21: Transfer rate peripheral to central (1/day)
        time_points: Optional custom time points (days)

    Returns:
        TKResults with concentration-time profile
    """
    if time_points is None:
        time_points = np.linspace(0, scenario.duration_days, 1000)

    # Get absorption fraction
    if scenario.route == ExposureRoute.ORAL:
        absorption = params.oral_bioavailability
    elif scenario.route == ExposureRoute.INHALATION:
        absorption = params.inhalation_absorption
    elif scenario.route == ExposureRoute.DERMAL:
        absorption = params.dermal_absorption
    else:
        absorption = 1.0

    dose_rate = scenario.dose_per_day * absorption

    # Compartment volumes (simplified)
    Vc = params.volume_of_distribution * scenario.body_weight_kg * 0.7  # Central
    Vp = params.volume_of_distribution * scenario.body_weight_kg * 0.3  # Peripheral

    k_el = params.elimination_rate_constant

    # Differential equations
    def model(y: list[float], t: float) -> list[float]:
        Cc, Cp = y[0], y[1]  # Central and peripheral concentrations
        dCc_dt = dose_rate / Vc - k12 * Cc + k21 * Cp * (Vp / Vc) - k_el * Cc
        dCp_dt = k12 * Cc * (Vc / Vp) - k21 * Cp
        return [dCc_dt, dCp_dt]

    # Solve
    y0 = [0.0, 0.0]
    solution = odeint(model, y0, time_points)
    blood_conc = solution[:, 0]  # Central compartment = blood

    # Tissue concentrations based on central compartment
    tissue_concs = {}
    for tissue, partition in params.tissue_partition_coeffs.items():
        tissue_concs[tissue] = blood_conc * partition

    cumulative = np.cumsum(np.full_like(time_points, dose_rate)) * np.diff(
        np.concatenate([[0], time_points])
    )

    auc = _trapezoid_integral(blood_conc, time_points)

    return TKResults(
        time_days=time_points,
        blood_concentration=blood_conc,
        tissue_concentrations=tissue_concs,
        cumulative_dose=cumulative,
        auc=auc,
    )


def predict_biomarker_response(
    metal: str,
    endpoint: str,
    blood_concentration: float,
    exposure_duration_days: float,
) -> dict[str, Any]:
    """Predict biomarker response based on TK-TD model.

    Args:
        metal: Metal name
        endpoint: Biomarker endpoint (e.g., "MDA", "8-OHdG")
        blood_concentration: Current blood concentration (mg/L)
        exposure_duration_days: Duration of exposure

    Returns:
        Dictionary with predicted response
    """
    td_params = TD_DATABASE.get((metal, endpoint))
    if not td_params:
        return {"error": f"No TD parameters for {metal}-{endpoint}"}

    # Hill equation for concentration-response
    # E = Ebaseline + (Emax * C^n) / (EC50^n + C^n)
    C = blood_concentration
    n = td_params.hill_coefficient

    if endpoint in ["SOD", "CAT", "GPx", "GSH"]:
        # For decreasing biomarkers, emax is a fraction
        fold_change = td_params.emax * (C ** n) / (td_params.ec50 ** n + C ** n)
        predicted = td_params.baseline_response * (1 - fold_change)
    else:
        # For increasing biomarkers
        fold_change = 1 + (td_params.emax - 1) * (C ** n) / (td_params.ec50 ** n + C ** n)
        predicted = td_params.baseline_response * fold_change

    # Time factor (response develops over time)
    time_fraction = min(exposure_duration_days / td_params.time_to_response_days, 1.0)
    time_adjusted_prediction = td_params.baseline_response + (
        predicted - td_params.baseline_response
    ) * time_fraction

    return {
        "metal": metal,
        "endpoint": endpoint,
        "blood_concentration_mg_L": blood_concentration,
        "baseline": td_params.baseline_response,
        "predicted_response": round(time_adjusted_prediction, 2),
        "fold_change": round(time_adjusted_prediction / td_params.baseline_response, 2),
        "time_to_full_response_days": td_params.time_to_response_days,
        "response_development": f"{time_fraction * 100:.0f}%",
    }


def calculate_daily_intake(
    metal: str,
    environmental_concentration: float,
    exposure_route: ExposureRoute,
    exposure_duration_hours: float = 24,
    inhalation_rate: float = 20,  # m3/day
    water_intake: float = 2,  # L/day
    body_weight: float = 70,  # kg
) -> float:
    """Calculate daily intake from environmental concentration.

    Args:
        metal: Metal name (for absorption factor lookup)
        environmental_concentration: mg/m3 (air), mg/L (water), mg/cm2 (skin)
        exposure_route: Route of exposure
        exposure_duration_hours: Hours per day exposed
        inhalation_rate: m3/day
        water_intake: L/day
        body_weight: kg

    Returns:
        Dose in mg/day
    """
    if exposure_route == ExposureRoute.INHALATION:
        # mg/m3 * m3/day = mg/day
        dose = environmental_concentration * inhalation_rate * (exposure_duration_hours / 24)
    elif exposure_route == ExposureRoute.ORAL:
        # mg/L * L/day = mg/day
        dose = environmental_concentration * water_intake
    elif exposure_route == ExposureRoute.DERMAL:
        # Simplified dermal calculation
        # Assume 2000 cm2 exposed skin, absorption
        exposed_area_cm2 = 2000
        dose = environmental_concentration * exposed_area_cm2 * 0.1  # mg/day
    else:
        dose = environmental_concentration

    return dose


def estimate_time_to_target_concentration(
    metal: str,
    species: Species,
    exposure_scenario: ExposureScenario,
    target_concentration: float,
    tolerance: float = 0.05,
) -> float | None:
    """Estimate time to reach target blood concentration.

    Args:
        metal: Metal name
        species: Species
        exposure_scenario: Exposure scenario
        target_concentration: Target concentration (mg/L)
        tolerance: Fractional tolerance (default 5%)

    Returns:
        Time in days to reach target, or None if not achievable
    """
    params = get_tk_params(metal, species)
    if not params:
        return None

    # Run simulation for extended period
    extended_duration = max(exposure_scenario.duration_days, 365 * 5)  # Up to 5 years
    extended_scenario = ExposureScenario(
        route=exposure_scenario.route,
        dose_per_day=exposure_scenario.dose_per_day,
        duration_days=extended_duration,
        body_weight_kg=exposure_scenario.body_weight_kg,
    )

    results = one_compartment_model(params, extended_scenario)

    # Find when within tolerance of target.
    # When the trajectory crosses the target between sampled time points,
    # interpolate linearly so the estimate is not dependent on grid spacing.
    lower_bound = target_concentration * (1 - tolerance)
    upper_bound = target_concentration * (1 + tolerance)

    peak_concentration = float(np.max(results.blood_concentration))
    if lower_bound > peak_concentration:
        return None

    for i, conc in enumerate(results.blood_concentration):
        if lower_bound <= conc <= upper_bound:
            return float(results.time_days[i])
        if i == 0:
            continue
        previous_conc = float(results.blood_concentration[i - 1])
        if previous_conc == conc:
            continue
        crossed_target = (
            previous_conc < target_concentration <= conc
            or previous_conc > target_concentration >= conc
        )
        if crossed_target:
            previous_time = float(results.time_days[i - 1])
            current_time = float(results.time_days[i])
            fraction = (target_concentration - previous_conc) / (conc - previous_conc)
            return float(previous_time + fraction * (current_time - previous_time))

    # Check if target is above steady state
    if target_concentration > (results.steady_state_concentration or float("inf")):
        return None

    return None


def get_tissue_dose_metrics(
    results: TKResults,
    tissue: str,
) -> dict[str, Any]:
    """Get dose metrics for a specific tissue.

    Args:
        results: TK simulation results
        tissue: Target tissue name

    Returns:
        Dictionary with tissue-specific metrics
    """
    if tissue not in results.tissue_concentrations:
        return {"error": f"Tissue {tissue} not in results"}

    tissue_conc = results.tissue_concentrations[tissue]

    return {
        "tissue": tissue,
        "cmax": round(float(np.max(tissue_conc)), 4),
        "cavg": round(float(np.mean(tissue_conc)), 4),
        "auc_tissue": round(_trapezoid_integral(tissue_conc, results.time_days), 4),
        "time_at_cmax": round(float(results.time_days[np.argmax(tissue_conc)]), 1),
    }


def compare_exposure_scenarios(
    metal: str,
    species: Species,
    scenarios: list[ExposureScenario],
    metrics: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Compare multiple exposure scenarios for a metal.

    Args:
        metal: Metal name
        species: Species
        scenarios: List of exposure scenarios to compare
        metrics: List of metrics to include (default: all)

    Returns:
        List of comparison results
    """
    params = get_tk_params(metal, species)
    if not params:
        return [{"error": f"No TK parameters for {metal}-{species.value}"}]

    results = []
    for i, scenario in enumerate(scenarios):
        tk_results = one_compartment_model(params, scenario)

        result = {
            "scenario_id": i + 1,
            "route": scenario.route.value,
            "dose_per_day_mg": scenario.dose_per_day,
            "duration_days": scenario.duration_days,
            "auc": round(tk_results.auc, 2),
            "cmax_mg_L": round(float(np.max(tk_results.blood_concentration)), 4),
            "steady_state_mg_L": round(tk_results.steady_state_concentration or 0, 4),
        }

        results.append(result)

    return results
