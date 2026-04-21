"""Cross-species extrapolation for metal carcinogen risk assessment.

This module provides tools for extrapolating toxicological data from
animal species to humans, including:
- Allometric scaling for dose/concentration conversions
- Body weight and lifespan normalization
- Metabolic rate corrections
- Uncertainty factor application
- Life-stage considerations (juvenile vs adult)
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

import numpy as np


class LifeStage(str, Enum):
    """Life stages for risk assessment."""

    FETAL = "fetal"
    NEONATAL = "neonatal"
    JUVENILE = "juvenile"
    ADULT = "adult"
    ELDERLY = "elderly"


class ExtrapolationType(str, Enum):
    """Types of cross-species extrapolation."""

    DOSE = "dose"  # mg/day
    DOSE_NORMALIZED = "dose_normalized"  # mg/kg/day
    CONCENTRATION = "concentration"  # mg/L
    AUC = "auc"  # mg·day/L
    CLEARANCE = "clearance"  # L/day
    VOLUME_DISTRIBUTION = "vd"  # L/kg


@dataclass
class SpeciesCharacteristics:
    """Physiological characteristics of a species."""

    species: str
    body_weight_kg: float  # Typical adult weight
    lifespan_years: float
    lifespan_days: float
    heart_rate_bpm: float
    respiratory_rate: float  # breaths/min
    cardiac_output_L_min: float
    hepatic_blood_flow_pct: float  # % of cardiac output
    renal_blood_flow_pct: float  # % of cardiac output
    glomerular_filtration_rate_ml_min: float
    basal_metabolic_rate_kcal_day: float
    surface_area_m2: float
    liver_weight_pct: float  # % of body weight
    kidney_weight_pct: float  # % of body weight
    brain_weight_pct: float  # % of body weight


@dataclass
class UncertaintyFactors:
    """Uncertainty factors for risk assessment."""

    interspecies_uf: float = 10.0  # Animal to human (default 10)
    intraspecies_uf: float = 10.0  # Human variability (default 10)
    subchronic_to_chronic: float = 1.0  # If using subchronic data
    LOAEL_to_NOAEL: float = 1.0  # If extrapolating from LOAEL
    database_incomplete: float = 1.0  # If database gaps
    modifying_factor: float = 1.0  # Professional judgment

    @property
    def total_uf(self) -> float:
        """Calculate total uncertainty factor."""
        return (
            self.interspecies_uf
            * self.intraspecies_uf
            * self.subchronic_to_chronic
            * self.LOAEL_to_NOAEL
            * self.database_incomplete
            * self.modifying_factor
        )


@dataclass
class AllometricCoefficients:
    """Allometric scaling coefficients (Y = a * BW^b)."""

    parameter: str
    coefficient_a: float
    exponent_b: float
    units: str


# ── Species Characteristics Database ──────────────────────────────────────

SPECIES_DATA: dict[str, SpeciesCharacteristics] = {
    "human": SpeciesCharacteristics(
        species="human",
        body_weight_kg=70.0,
        lifespan_years=80.0,
        lifespan_days=80 * 365,
        heart_rate_bpm=70.0,
        respiratory_rate=15.0,
        cardiac_output_L_min=5.0,
        hepatic_blood_flow_pct=25.0,
        renal_blood_flow_pct=20.0,
        glomerular_filtration_rate_ml_min=125.0,
        basal_metabolic_rate_kcal_day=1500.0,
        surface_area_m2=1.8,
        liver_weight_pct=2.5,
        kidney_weight_pct=0.4,
        brain_weight_pct=2.0,
    ),
    "mouse": SpeciesCharacteristics(
        species="mouse",
        body_weight_kg=0.025,  # 25g
        lifespan_years=2.0,
        lifespan_days=2 * 365,
        heart_rate_bpm=600.0,
        respiratory_rate=120.0,
        cardiac_output_L_min=0.015,
        hepatic_blood_flow_pct=25.0,
        renal_blood_flow_pct=15.0,
        glomerular_filtration_rate_ml_min=0.28,
        basal_metabolic_rate_kcal_day=3.5,
        surface_area_m2=0.0045,
        liver_weight_pct=5.0,
        kidney_weight_pct=1.5,
        brain_weight_pct=1.7,
    ),
    "rat": SpeciesCharacteristics(
        species="rat",
        body_weight_kg=0.25,  # 250g
        lifespan_years=2.5,
        lifespan_days=2.5 * 365,
        heart_rate_bpm=350.0,
        respiratory_rate=80.0,
        cardiac_output_L_min=0.07,
        hepatic_blood_flow_pct=25.0,
        renal_blood_flow_pct=15.0,
        glomerular_filtration_rate_ml_min=1.5,
        basal_metabolic_rate_kcal_day=20.0,
        surface_area_m2=0.025,
        liver_weight_pct=4.0,
        kidney_weight_pct=0.7,
        brain_weight_pct=0.6,
    ),
    "rabbit": SpeciesCharacteristics(
        species="rabbit",
        body_weight_kg=2.5,
        lifespan_years=8.0,
        lifespan_days=8 * 365,
        heart_rate_bpm=200.0,
        respiratory_rate=40.0,
        cardiac_output_L_min=0.35,
        hepatic_blood_flow_pct=25.0,
        renal_blood_flow_pct=18.0,
        glomerular_filtration_rate_ml_min=8.0,
        basal_metabolic_rate_kcal_day=150.0,
        surface_area_m2=0.15,
        liver_weight_pct=4.0,
        kidney_weight_pct=0.6,
        brain_weight_pct=0.4,
    ),
    "dog": SpeciesCharacteristics(
        species="dog",
        body_weight_kg=10.0,
        lifespan_years=12.0,
        lifespan_days=12 * 365,
        heart_rate_bpm=90.0,
        respiratory_rate=20.0,
        cardiac_output_L_min=1.2,
        hepatic_blood_flow_pct=25.0,
        renal_blood_flow_pct=20.0,
        glomerular_filtration_rate_ml_min=50.0,
        basal_metabolic_rate_kcal_day=800.0,
        surface_area_m2=0.5,
        liver_weight_pct=3.5,
        kidney_weight_pct=0.5,
        brain_weight_pct=0.8,
    ),
    "monkey": SpeciesCharacteristics(
        species="monkey",
        body_weight_kg=5.0,
        lifespan_years=25.0,
        lifespan_days=25 * 365,
        heart_rate_bpm=120.0,
        respiratory_rate=30.0,
        cardiac_output_L_min=0.6,
        hepatic_blood_flow_pct=25.0,
        renal_blood_flow_pct=20.0,
        glomerular_filtration_rate_ml_min=30.0,
        basal_metabolic_rate_kcal_day=400.0,
        surface_area_m2=0.3,
        liver_weight_pct=2.8,
        kidney_weight_pct=0.5,
        brain_weight_pct=1.2,
    ),
}


# Allometric scaling coefficients for physiological parameters
# Y = a * BW^b where BW is in kg
ALLOMETRIC_COEFFICIENTS: dict[str, AllometricCoefficients] = {
    "basal_metabolic_rate": AllometricCoefficients(
        parameter="BMR",
        coefficient_a=70.0,  # kcal/day
        exponent_b=0.75,  # Kleiber's law
        units="kcal/day",
    ),
    "cardiac_output": AllometricCoefficients(
        parameter="Cardiac Output",
        coefficient_a=0.166,
        exponent_b=0.75,
        units="L/min",
    ),
    "glomerular_filtration_rate": AllometricCoefficients(
        parameter="GFR",
        coefficient_a=4.2,
        exponent_b=0.77,
        units="mL/min",
    ),
    "liver_blood_flow": AllometricCoefficients(
        parameter="Liver Blood Flow",
        coefficient_a=0.042,
        exponent_b=0.75,
        units="L/min",
    ),
    "kidney_blood_flow": AllometricCoefficients(
        parameter="Kidney Blood Flow",
        coefficient_a=0.033,
        exponent_b=0.75,
        units="L/min",
    ),
    "tidal_volume": AllometricCoefficients(
        parameter="Tidal Volume",
        coefficient_a=7.0,
        exponent_b=1.0,
        units="mL",
    ),
    "minute_volume": AllometricCoefficients(
        parameter="Minute Ventilation",
        coefficient_a=0.21,
        exponent_b=0.75,
        units="L/min",
    ),
}


# Life-stage specific adjustment factors (relative to adult)
LIFE_STAGE_FACTORS: dict[str, dict[str, float]] = {
    "absorption": {
        LifeStage.FETAL.value: 1.0,
        LifeStage.NEONATAL.value: 1.2,  # Often higher
        LifeStage.JUVENILE.value: 1.1,
        LifeStage.ADULT.value: 1.0,
        LifeStage.ELDERLY.value: 0.9,  # Often reduced
    },
    "metabolism": {
        LifeStage.FETAL.value: 0.3,
        LifeStage.NEONATAL.value: 0.5,
        LifeStage.JUVENILE.value: 0.8,
        LifeStage.ADULT.value: 1.0,
        LifeStage.ELDERLY.value: 0.7,
    },
    "elimination": {
        LifeStage.FETAL.value: 0.2,
        LifeStage.NEONATAL.value: 0.4,
        LifeStage.JUVENILE.value: 0.8,
        LifeStage.ADULT.value: 1.0,
        LifeStage.ELDERLY.value: 0.7,
    },
    "body_water": {
        LifeStage.FETAL.value: 0.8,
        LifeStage.NEONATAL.value: 0.75,
        LifeStage.JUVENILE.value: 0.65,
        LifeStage.ADULT.value: 0.60,  # Adult male
        LifeStage.ELDERLY.value: 0.55,
    },
}


# ── Allometric Scaling Functions ─────────────────────────────────────────


def get_species_data(species: str) -> SpeciesCharacteristics | None:
    """Get physiological data for a species.

    Args:
        species: Species name (human, mouse, rat, etc.)

    Returns:
        SpeciesCharacteristics if available
    """
    return SPECIES_DATA.get(species.lower())


def allometric_scale(
    value: float,
    from_species_bw: float,
    to_species_bw: float,
    exponent: float = 0.75,
) -> float:
    """Apply allometric scaling between species.

    Args:
        value: Parameter value in from_species
        from_species_bw: Body weight of source species (kg)
        to_species_bw: Body weight of target species (kg)
        exponent: Allometric exponent (default 0.75 for most parameters)

    Returns:
        Scaled value for target species
    """
    return float(value * (to_species_bw / from_species_bw) ** exponent)


def convert_dose_mg_per_kg(
    dose: float,
    from_species: str,
    to_species: str,
    extrapolation_type: ExtrapolationType = ExtrapolationType.DOSE_NORMALIZED,
) -> dict[str, Any]:
    """Convert dose from one species to another.

    Supports multiple extrapolation approaches:
    - Body weight scaling (1.0 exponent)
    - Allometric scaling (0.75 exponent - metabolic rate)
    - Brain weight scaling (0.66 exponent - CNS effects)
    - Life span scaling (0.25 exponent - chronic effects)

    Args:
        dose: Dose in mg/kg/day (from_species)
        from_species: Source species name
        to_species: Target species name
        extrapolation_type: Type of extrapolation

    Returns:
        Dictionary with different extrapolation results
    """
    from_data = get_species_data(from_species)
    to_data = get_species_data(to_species)

    if not from_data or not to_data:
        return {"error": "Species data not available"}

    # Different scaling approaches
    extrapolations: dict[str, float] = {}
    results: dict[str, Any] = {
        "from_species": from_species,
        "to_species": to_species,
        "input_dose_mg_kg_day": dose,
        "from_bw_kg": from_data.body_weight_kg,
        "to_bw_kg": to_data.body_weight_kg,
        "extrapolations": extrapolations,
    }

    # 1. Simple body weight scaling (dose in mg/day, then normalize)
    dose_mg_day = dose * from_data.body_weight_kg
    bw_scaled = dose_mg_day / to_data.body_weight_kg
    extrapolations["body_weight_scaled"] = round(bw_scaled, 4)

    # 2. Allometric scaling (0.75 exponent - metabolic rate basis)
    # HED = Animal Dose × (Animal BW / Human BW)^0.33
    allometric_factor = (from_data.body_weight_kg / to_data.body_weight_kg) ** 0.33
    allometric_scaled = dose * allometric_factor
    extrapolations["allometric_scaled_0.75"] = round(allometric_scaled, 4)

    # 3. Surface area scaling (0.67 exponent - EPA approach)
    surface_factor = (from_data.body_weight_kg / to_data.body_weight_kg) ** (1/3)
    surface_scaled = dose * surface_factor
    extrapolations["surface_area_scaled"] = round(surface_scaled, 4)

    # 4. Brain weight scaling (0.66 exponent - for CNS effects)
    # Brain weight ratio approximated from literature
    brain_factor = (from_data.body_weight_kg / to_data.body_weight_kg) ** 0.66
    brain_scaled = dose * brain_factor
    extrapolations["brain_weight_scaled"] = round(brain_scaled, 4)

    # 5. Lifespan adjustment (for chronic effects)
    lifespan_ratio = from_data.lifespan_days / to_data.lifespan_days
    lifespan_adjusted = allometric_scaled / lifespan_ratio
    extrapolations["lifespan_adjusted"] = round(lifespan_adjusted, 4)

    # Recommended approach for metals (tends to be allometric or surface area)
    results["recommended_dose_mg_kg_day"] = round(allometric_scaled, 4)
    results["recommended_approach"] = "allometric_scaled_0.75"

    return results


def convert_clearance(
    clearance: float,
    from_species: str,
    to_species: str,
) -> dict[str, Any]:
    """Convert clearance between species using allometric scaling.

    Args:
        clearance: Clearance value (L/hr or L/day)
        from_species: Source species name
        to_species: Target species name

    Returns:
        Dictionary with scaled clearance values
    """
    from_data = get_species_data(from_species)
    to_data = get_species_data(to_species)

    if not from_data or not to_data:
        return {"error": "Species data not available"}

    # Clearance scales with BMR exponent (0.75)
    scaled_clearance = allometric_scale(
        clearance,
        from_data.body_weight_kg,
        to_data.body_weight_kg,
        exponent=0.75,
    )

    # Alternative: Linear scaling (for glomerular filtration)
    linear_scaled = allometric_scale(
        clearance,
        from_data.body_weight_kg,
        to_data.body_weight_kg,
        exponent=1.0,
    )

    return {
        "from_species": from_species,
        "to_species": to_species,
        "input_clearance": clearance,
        "allometric_scaled_0.75": round(scaled_clearance, 4),
        "linear_scaled": round(linear_scaled, 4),
        "recommended": round(scaled_clearance, 4),
    }


def convert_volume_distribution(
    vd: float,
    from_species: str,
    to_species: str,
) -> dict[str, Any]:
    """Convert volume of distribution between species.

    Args:
        vd: Volume of distribution (L or L/kg)
        from_species: Source species name
        to_species: Target species name

    Returns:
        Dictionary with scaled Vd values
    """
    from_data = get_species_data(from_species)
    to_data = get_species_data(to_species)

    if not from_data or not to_data:
        return {"error": "Species data not available"}

    # Vd as total volume scales with body weight (1.0 exponent)
    total_vd_scaled = allometric_scale(
        vd,
        from_data.body_weight_kg,
        to_data.body_weight_kg,
        exponent=1.0,
    )

    # Vd per kg body weight stays approximately constant
    vd_per_kg = vd / from_data.body_weight_kg
    normalized_vd = vd_per_kg * to_data.body_weight_kg

    return {
        "from_species": from_species,
        "to_species": to_species,
        "input_vd": vd,
        "total_volume_scaled": round(total_vd_scaled, 4),
        "normalized_per_kg": round(vd_per_kg, 4),
        "recommended_total": round(normalized_vd, 4),
    }


def convert_auc(
    auc: float,
    from_species: str,
    to_species: str,
    dose_mg: float,
    clearance_method: str = "allometric",
) -> dict[str, Any]:
    """Convert AUC between species.

    AUC = Dose / Clearance, so scaling depends on clearance scaling.

    Args:
        auc: Area under curve (mg·day/L or mg·hr/L)
        from_species: Source species name
        to_species: Target species name
        dose_mg: Dose used to generate AUC
        clearance_method: Method for clearance scaling

    Returns:
        Dictionary with scaled AUC
    """
    from_data = get_species_data(from_species)
    to_data = get_species_data(to_species)

    if not from_data or not to_data:
        return {"error": "Species data not available"}

    # First scale the dose using allometric approach
    dose_ratio = (to_data.body_weight_kg / from_data.body_weight_kg) ** 0.75
    scaled_dose = dose_mg * dose_ratio

    # Scale clearance
    if clearance_method == "allometric":
        clearance_ratio = (to_data.body_weight_kg / from_data.body_weight_kg) ** 0.75
    else:  # linear
        clearance_ratio = to_data.body_weight_kg / from_data.body_weight_kg

    # AUC scales with (dose ratio) / (clearance ratio)
    auc_scale_factor = dose_ratio / clearance_ratio
    scaled_auc = auc * auc_scale_factor

    return {
        "from_species": from_species,
        "to_species": to_species,
        "input_auc": auc,
        "dose_mg": dose_mg,
        "scaled_dose_mg": round(scaled_dose, 4),
        "scaled_auc": round(scaled_auc, 4),
        "scale_factor": round(auc_scale_factor, 4),
    }


# ── Uncertainty Factor Functions ─────────────────────────────────────────


def apply_uncertainty_factors(
    dose: float,
    ufs: UncertaintyFactors,
    apply_interspecies: bool = True,
    apply_intraspecies: bool = True,
) -> dict[str, Any]:
    """Apply uncertainty factors to a dose for risk assessment.

    Args:
        dose: Dose value (mg/kg/day or other unit)
        ufs: UncertaintyFactors object
        apply_interspecies: Whether to apply interspecies UF
        apply_intraspecies: Whether to apply intraspecies UF

    Returns:
        Dictionary with factored doses and explanation
    """
    # Calculate individual factors to apply
    effective_interspecies = ufs.interspecies_uf if apply_interspecies else 1.0
    effective_intraspecies = ufs.intraspecies_uf if apply_intraspecies else 1.0

    other_factors = (
        ufs.subchronic_to_chronic
        * ufs.LOAEL_to_NOAEL
        * ufs.database_incomplete
        * ufs.modifying_factor
    )

    total_applied = effective_interspecies * effective_intraspecies * other_factors

    factored_dose = dose / total_applied

    return {
        "original_dose": dose,
        "adjusted_dose": round(factored_dose, 6),
        "total_uf_applied": round(total_applied, 1),
        "factors": {
            "interspecies": effective_interspecies if apply_interspecies else "N/A",
            "intraspecies": effective_intraspecies if apply_intraspecies else "N/A",
            "subchronic_to_chronic": ufs.subchronic_to_chronic,
            "LOAEL_to_NOAEL": ufs.LOAEL_to_NOAEL,
            "database_incomplete": ufs.database_incomplete,
            "modifying_factor": ufs.modifying_factor,
        },
    }


def calculate_reference_dose(
    noael_or_loael: float,
    is_loael: bool = False,
    critical_effect: str = "",
    database_completeness: str = "complete",
    use_default_ufs: bool = True,
) -> dict[str, Any]:
    """Calculate Reference Dose (RfD) or Tolerable Daily Intake.

    Args:
        noael_or_loael: NOAEL or LOAEL from animal/human study
        is_loael: Whether the value is LOAEL (not NOAEL)
        critical_effect: Description of critical effect
        database_completeness: "complete", "partial", or "minimal"
        use_default_ufs: Use default 10x interspecies and 10x intraspecies

    Returns:
        Dictionary with RfD calculation
    """
    ufs = UncertaintyFactors()

    if use_default_ufs:
        ufs.interspecies_uf = 10.0
        ufs.intraspecies_uf = 10.0

    if is_loael:
        ufs.LOAEL_to_NOAEL = 3.0  # Default LOAEL to NOAEL factor

    if database_completeness == "partial":
        ufs.database_incomplete = 3.0
    elif database_completeness == "minimal":
        ufs.database_incomplete = 10.0

    # Apply uncertainty factors
    rfd = noael_or_loael / ufs.total_uf

    return {
        "starting_point": "LOAEL" if is_loael else "NOAEL",
        "starting_value": noael_or_loael,
        "critical_effect": critical_effect,
        "reference_dose_mg_kg_day": round(rfd, 6),
        "total_uf": ufs.total_uf,
        "uf_breakdown": {
            "interspecies": ufs.interspecies_uf,
            "intraspecies": ufs.intraspecies_uf,
            "LOAEL_to_NOAEL": ufs.LOAEL_to_NOAEL,
            "database_incomplete": ufs.database_incomplete,
        },
    }


# ── Life-Stage Functions ──────────────────────────────────────────────────


def get_life_stage_factor(
    parameter: str,
    life_stage: LifeStage | str,
) -> float:
    """Get adjustment factor for a life stage.

    Args:
        parameter: Parameter name (absorption, metabolism, elimination, body_water)
        life_stage: Life stage enum or string

    Returns:
        Adjustment factor (multiplier relative to adult)
    """
    stage_key = life_stage.value if isinstance(life_stage, LifeStage) else life_stage

    if parameter not in LIFE_STAGE_FACTORS:
        return 1.0

    return LIFE_STAGE_FACTORS[parameter].get(stage_key, 1.0)


def adjust_for_life_stage(
    dose: float,
    life_stage: LifeStage | str,
    adjustment_parameters: list[str] | None = None,
) -> dict[str, Any]:
    """Adjust dose for life stage considerations.

    Args:
        dose: Adult-equivalent dose (mg/kg/day)
        life_stage: Target life stage
        adjustment_parameters: Parameters to adjust (default: all)

    Returns:
        Dictionary with adjusted dose and factors
    """
    if adjustment_parameters is None:
        adjustment_parameters = ["absorption", "metabolism", "elimination", "body_water"]

    stage_key = life_stage.value if isinstance(life_stage, LifeStage) else life_stage

    factors = {}
    for param in adjustment_parameters:
        factors[param] = get_life_stage_factor(param, life_stage)

    # Combined adjustment (geometric mean of individual factors)
    # This is a simplified approach
    combined_factor = np.prod(list(factors.values())) ** (1 / len(factors))

    adjusted_dose = dose * combined_factor

    return {
        "adult_dose": dose,
        "life_stage": stage_key,
        "adjustment_factors": factors,
        "combined_factor": round(combined_factor, 3),
        "adjusted_dose": round(adjusted_dose, 6),
    }


def pediatric_adjustment(
    adult_dose: float,
    child_age_months: int,
    metal: str | None = None,
) -> dict[str, Any]:
    """Specialized pediatric dose adjustment.

    Args:
        adult_dose: Adult reference dose (mg/kg/day)
        child_age_months: Child age in months
        metal: Optional metal name for metal-specific adjustments

    Returns:
        Dictionary with pediatric dose recommendations
    """
    # Determine life stage
    if child_age_months < 12:
        stage = LifeStage.NEONATAL
    elif child_age_months < 36:
        stage = LifeStage.JUVENILE
    else:
        stage = LifeStage.ADULT

    # Base adjustment
    base_adjustment = adjust_for_life_stage(adult_dose, stage)

    # Additional considerations for children
    considerations: dict[str, Any] = {
        "higher_absorption": child_age_months < 12,  # Gut more permeable
        "developing_organs": child_age_months < 72,  # Brain, kidney developing
        "weight_normalized": True,
    }

    # Lead-specific: children absorb ~50% vs adults ~15%
    if metal and metal.lower() == "lead":
        if child_age_months < 72:
            # Higher GI absorption in children
            absorption_boost = 3.0  # 50% / 15% ≈ 3.3
            base_adjustment["adjusted_dose"] /= absorption_boost
            considerations["lead_child_adjustment"] = f"divided by {absorption_boost:.1f}"

    return {
        "adult_dose": adult_dose,
        "child_age_months": child_age_months,
        "life_stage": stage.value,
        "recommended_dose_mg_kg_day": round(base_adjustment["adjusted_dose"], 6),
        "adjustment_factor": base_adjustment["combined_factor"],
        "considerations": considerations,
    }


# ── Species Comparison Functions ─────────────────────────────────────────


def compare_species_physiology(
    species_list: list[str],
) -> dict[str, Any]:
    """Compare physiological parameters across species.

    Args:
        species_list: List of species names to compare

    Returns:
        Dictionary with comparison table
    """
    data = {}
    for species in species_list:
        char = get_species_data(species)
        if char:
            data[species] = {
                "body_weight_kg": char.body_weight_kg,
                "lifespan_years": char.lifespan_years,
                "heart_rate_bpm": char.heart_rate_bpm,
                "gfr_ml_min": char.glomerular_filtration_rate_ml_min,
                "bmr_kcal_day": char.basal_metabolic_rate_kcal_day,
                "surface_area_m2": char.surface_area_m2,
            }

    return {
        "comparison": data,
        "body_weight_ratios": {
            f"{s1}:{s2}": round(
                data[s1]["body_weight_kg"] / data[s2]["body_weight_kg"], 2
            )
            for s1 in data
            for s2 in data
            if s1 != s2
        },
    }


def extrapolate_study_duration(
        duration_days_source: int,
        from_species: str,
        to_species: str,
        method: str = "lifespan_fraction",
    ) -> dict[str, Any]:
        """Extrapolate equivalent study duration between species.

        Args:
            duration_days_source: Study duration in source species
            from_species: Source species name
            to_species: Target species name
            method: "lifespan_fraction", "body_weight", or "chronic_ratio"

        Returns:
            Dictionary with extrapolated duration
        """
        from_data = get_species_data(from_species)
        to_data = get_species_data(to_species)

        if not from_data or not to_data:
            return {"error": "Species data not available"}

        if method == "lifespan_fraction":
            # Fraction of lifespan approach
            lifespan_fraction = duration_days_source / from_data.lifespan_days
            equivalent_days = lifespan_fraction * to_data.lifespan_days
        elif method == "body_weight":
            # Allometric scaling of time (exponent -0.25)
            time_ratio = (to_data.body_weight_kg / from_data.body_weight_kg) ** (-0.25)
            equivalent_days = duration_days_source * time_ratio
        elif method == "chronic_ratio":
            # Standard chronic study ratios (2 years rat = lifetime human equivalent)
            if from_species == "rat" and duration_days_source >= 2 * 365 - 30:
                # Chronic rat study ≈ 70-year human equivalent
                equivalent_days = 70 * 365
            else:
                lifespan_fraction = duration_days_source / from_data.lifespan_days
                equivalent_days = lifespan_fraction * to_data.lifespan_days
        else:
            return {"error": f"Unknown method: {method}"}

        return {
            "from_species": from_species,
            "to_species": to_species,
            "input_duration_days": duration_days_source,
            "method": method,
            "equivalent_duration_days": round(equivalent_days, 1),
            "equivalent_years": round(equivalent_days / 365, 2),
            "lifespan_fraction_source": round(
                duration_days_source / from_data.lifespan_days * 100, 2
            ),
        }
