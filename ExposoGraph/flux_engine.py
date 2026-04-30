"""Quantitative metabolic flux modeling engine.

Computes carcinogen activation and detoxification fluxes using published
enzyme kinetic parameters (Michaelis-Menten / Hill equation) with genotype
and tissue expression modifiers.  Produces net flux ratios interpretable
as individualized cancer-risk indicators.

Primary measured kinetics live in ``kinetic_parameters.json``.
Classes that currently require receptor-mediated or semi-quantitative proxy
models load their coefficients from ``proxy_flux_parameters.json`` and
supporting exposure defaults from ``exposure_database.json``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Mapping, TypeAlias, cast

# ── Enums ──────────────────────────────────────────────────────────────────


class CarcinogenClass(str, Enum):
    """Supported carcinogen classes for quantitative flux modeling."""

    PAH = "PAH"
    AFLATOXIN = "Aflatoxin"
    ALDEHYDE = "Aldehyde"
    NITROSAMINE = "Nitrosamine"
    NDMA = "NDMA"
    NDEA = "NDEA"
    HCA = "HCA"
    AROMATIC_AMINES = "AromaticAmines"
    ESTROGEN_METABOLITES = "EstrogenMetabolites"
    BENZENE = "Benzene"
    VINYL_CHLORIDE = "VinylChloride"
    CHLORINATED_SOLVENT = "ChlorinatedSolvent"
    UV_RADIATION = "UV_Radiation"
    DIOXIN = "Dioxin"
    HEAVY_METAL = "HeavyMetal"


class RiskClassification(str, Enum):
    """Risk tier derived from activation / detoxification net ratio."""

    PROTECTIVE = "PROTECTIVE"
    LOW = "LOW"
    MODERATE = "MODERATE"
    ELEVATED = "ELEVATED"
    HIGH = "HIGH"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


class FluxTissueWeightSource(str, Enum):
    """Source used for tissue expression weights in flux calculations."""

    CURATED = "curated"
    GTEX = "gtex"


JsonDict: TypeAlias = dict[str, Any]
GenotypeMap: TypeAlias = dict[str, str]
FluxResultDict: TypeAlias = dict[str, Any]
FluxCalculator: TypeAlias = Callable[
    [GenotypeMap, str, float, FluxTissueWeightSource],
    FluxResultDict,
]
LifestyleMap: TypeAlias = Mapping[str, bool | int | float]


# ── Dataclasses ────────────────────────────────────────────────────────────


@dataclass
class EnzymeFlux:
    """Flux result for a single enzyme in a pathway.

    ``model_kind`` / ``parameter_source`` identify whether a term came from
    measured kinetics or a proxy block. Proxy terms also expose provenance
    fields so downstream reporting can cite bundled evidence without reopening
    the JSON files.
    """

    enzyme: str
    flux: float
    genotype_modifier: float
    tissue_weight: float
    confidence: str
    induction_modifier: float = 1.0
    qivive_scale: float = 1.0
    fraction: float = 0.0
    kinetics: str = "michaelis_menten"
    note: str = ""
    model_kind: str = "measured_kinetics"
    parameter_source: str = "kinetic_parameters.json"
    provenance_ref: str = ""
    provenance_sources: list[str] = field(default_factory=list)
    parameter_basis: str = ""


@dataclass
class PathwayFluxResult:
    """Result of pathway flux computation for one carcinogen class."""

    carcinogen_class: str
    tissue: str
    substrate_concentration_uM: float
    genotypes_used: dict[str, str]
    activation_enzymes: list[EnzymeFlux]
    detox_enzymes: list[EnzymeFlux]
    total_activation: float
    total_detox: float
    net_ratio: float
    susceptibility_score_log2: float
    risk_classification: RiskClassification
    tissue_weight_source: FluxTissueWeightSource
    model_kind: str = "measured_kinetics"
    parameter_source: str = "kinetic_parameters.json"
    unit_note: str = ""
    warnings: list[str] = field(default_factory=list)
    induction_factors_used: dict[str, float] = field(default_factory=dict)
    qivive_applied: bool = False
    qivive_context: dict[str, float] = field(default_factory=dict)
    steady_state_concentrations_uM: dict[str, float] = field(default_factory=dict)
    steady_state_model: dict[str, Any] = field(default_factory=dict)
    steady_state_concentration_proxy_uM: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class FluxSteadyStateResult:
    """Steady-state concentrations from the flux-coupled PBPK compartment model."""

    concentrations_uM: dict[str, float]
    model: dict[str, Any]


@dataclass
class FullProfileResult:
    """Result of computing flux across all carcinogen classes."""

    tissue: str
    genotypes: dict[str, str]
    per_class_results: dict[str, PathwayFluxResult]
    elevated_or_high_risk_classes: list[str]
    moderate_risk_classes: list[str]
    total_classes_modeled: int
    tissue_weight_source: FluxTissueWeightSource


@dataclass
class SensitivityResult:
    """Result of single-gene sensitivity analysis."""

    carcinogen_class: str
    gene_varied: str
    tissue: str
    baseline_ratio: float
    results_by_phenotype: dict[str, dict[str, Any]]
    max_fold_change: float | None
    tissue_weight_source: FluxTissueWeightSource


# ── Data loading (lazy cache) ──────────────────────────────────────────────

_KINETIC_PARAMS_FILE = Path(__file__).parent / "data" / "kinetic_parameters.json"
_KINETIC_CACHE: JsonDict | None = None
_EXPOSURE_DB_FILE = Path(__file__).parent / "data" / "exposure_database.json"
_EXPOSURE_DB_CACHE: JsonDict | None = None
_INTERACTION_PARAMS_FILE = Path(__file__).parent / "data" / "interaction_parameters.json"
_INTERACTION_CACHE: JsonDict | None = None
_PROXY_FLUX_PARAMS_FILE = Path(__file__).parent / "data" / "proxy_flux_parameters.json"
_PROXY_FLUX_CACHE: JsonDict | None = None
_PROXY_FLUX_PROVENANCE_FILE = Path(__file__).parent / "data" / "proxy_flux_provenance.json"
_PROXY_FLUX_PROVENANCE_CACHE: JsonDict | None = None


def _load_kinetic_params() -> JsonDict:
    """Lazy-load and cache the kinetic parameters JSON."""
    global _KINETIC_CACHE
    if _KINETIC_CACHE is None:
        if not _KINETIC_PARAMS_FILE.exists():
            raise FileNotFoundError(
                f"Kinetic parameters not found at {_KINETIC_PARAMS_FILE}. "
                "Ensure kinetic_parameters.json is in the data/ directory."
            )
        with open(_KINETIC_PARAMS_FILE, "r", encoding="utf-8") as fh:
            _KINETIC_CACHE = cast(JsonDict, json.load(fh))
    return _KINETIC_CACHE


def _load_exposure_db() -> JsonDict:
    """Lazy-load and cache the exposure database JSON."""
    global _EXPOSURE_DB_CACHE
    if _EXPOSURE_DB_CACHE is None:
        with open(_EXPOSURE_DB_FILE, "r", encoding="utf-8") as fh:
            _EXPOSURE_DB_CACHE = cast(JsonDict, json.load(fh))
    return _EXPOSURE_DB_CACHE


def _load_interaction_params() -> JsonDict:
    """Lazy-load and cache the interaction-parameter JSON."""
    global _INTERACTION_CACHE
    if _INTERACTION_CACHE is None:
        with open(_INTERACTION_PARAMS_FILE, "r", encoding="utf-8") as fh:
            _INTERACTION_CACHE = cast(JsonDict, json.load(fh))
    return _INTERACTION_CACHE


def _load_proxy_flux_params() -> JsonDict:
    """Lazy-load and cache the proxy flux parameter JSON."""
    global _PROXY_FLUX_CACHE
    if _PROXY_FLUX_CACHE is None:
        with open(_PROXY_FLUX_PARAMS_FILE, "r", encoding="utf-8") as fh:
            _PROXY_FLUX_CACHE = cast(JsonDict, json.load(fh))
    return _PROXY_FLUX_CACHE


def _load_proxy_flux_provenance() -> JsonDict:
    """Lazy-load and cache the proxy flux provenance JSON."""
    global _PROXY_FLUX_PROVENANCE_CACHE
    if _PROXY_FLUX_PROVENANCE_CACHE is None:
        with open(_PROXY_FLUX_PROVENANCE_FILE, "r", encoding="utf-8") as fh:
            _PROXY_FLUX_PROVENANCE_CACHE = cast(JsonDict, json.load(fh))
    return _PROXY_FLUX_PROVENANCE_CACHE


# ── Core kinetic equations ─────────────────────────────────────────────────


def michaelis_menten(S: float, Vmax: float, Km: float) -> float:
    """Calculate reaction velocity using Michaelis-Menten kinetics.

    v = Vmax * S / (Km + S)

    Args:
        S: Substrate concentration (uM).
        Vmax: Maximum velocity.
        Km: Michaelis constant (uM).

    Returns:
        Reaction velocity in same units as *Vmax*.
    """
    if Km <= 0:
        raise ValueError(f"Km must be positive, got {Km}")
    if S < 0:
        raise ValueError(f"Substrate concentration cannot be negative, got {S}")
    return Vmax * S / (Km + S)


def hill_equation(S: float, Vmax: float, K50: float, n: float) -> float:
    """Calculate reaction velocity using Hill (cooperative) kinetics.

    v = Vmax * S^n / (K50^n + S^n)

    Args:
        S: Substrate concentration (uM).
        Vmax: Maximum velocity.
        K50: Half-maximal concentration (uM).
        n: Hill coefficient (n=1 reduces to Michaelis-Menten).

    Returns:
        Reaction velocity in same units as *Vmax*.
    """
    if K50 <= 0:
        raise ValueError(f"K50 must be positive, got {K50}")
    if n <= 0:
        raise ValueError(f"Hill coefficient n must be positive, got {n}")
    if S < 0:
        raise ValueError(f"Substrate concentration cannot be negative, got {S}")
    Sn = S ** n
    K50n = K50 ** n
    return float(Vmax * Sn / (K50n + Sn))


# ── Modifier functions ─────────────────────────────────────────────────────


_GST_NULL_RESIDUAL_ACTIVITY = 0.05


def genotype_modifier(diplotype: str, gene: str) -> float:
    """Return Vmax scaling factor (0.0-2.0) based on metaboliser phenotype.

    Standard scale: PM=0.0, IM=0.5, NM=1.0, RM=1.5, UM=2.0.
    Special cases for ALDH2 heterozygotes, GSTM1/GSTT1 copy-number states,
    and cohort-specific CYP aliases used by the extension modules.

    Args:
        diplotype: Phenotype label (for example ``"NM"``, ``"null"``,
            or ``"*1/*2"``).
        gene: Gene name (for example ``"GSTM1"`` or ``"ALDH2"``).

    Returns:
        Scaling factor between 0.0 and 2.0.
    """
    params = _load_kinetic_params()
    special = params["genotype_modifiers"]["special_cases"]
    std = params["genotype_modifiers"]["standard_scale"]

    diplotype_lower = diplotype.lower().strip()
    gene_upper = gene.upper().strip()

    # Gene-specific manuscript/reference aliases that are more precise than
    # the generic PM/IM/NM/RM/UM scale.
    if gene_upper == "CYP1A2":
        if diplotype_lower in ("*1f/*1f", "1f/1f", "cyp1a2*1f/*1f", "um_1f_1f"):
            return 1.5
        if diplotype_lower in ("*1a/*1f", "1a/1f", "*1f/*1a", "1f/1a"):
            # Heterozygous *1F: intermediate inducibility (Sachse et al. 1999;
            # Ghotbi et al. 2007). Splits the difference between *1A/*1A (1.0)
            # and *1F/*1F (1.5).
            return 1.25
        if diplotype_lower in ("*1a/*1a", "1a/1a"):
            return 1.0
        if diplotype_lower in ("*1k", "*1k/*1k", "1k/1k"):
            return 0.5
        if diplotype_lower in ("pm", "poor", "poor metabolizer", "poor_metabolizer"):
            return 0.3

    if gene_upper == "NAT2":
        if diplotype_lower in ("slow", "slow acetylator", "slow_acetylator", "sa"):
            return 0.2
        if diplotype_lower in ("intermediate", "intermediate acetylator", "intermediate_acetylator"):
            return 0.5
        if diplotype_lower in ("rapid", "rapid acetylator", "rapid_acetylator", "ra"):
            return 1.0

    if gene_upper == "NAT1":
        # NAT1 phenotype assignments: *4 is the reference (rapid), *10 has
        # been associated with modestly increased acetylation activity in
        # several reports (Bell et al. 1995; Lin et al. 1998), and *14 is a
        # well-established slow allele (Hughes et al. 1998).
        if diplotype_lower in ("*4/*4", "4/4", "rapid", "ra"):
            return 1.0
        if diplotype_lower in ("*4/*10", "*10/*4", "4/10", "10/4"):
            return 1.05
        if diplotype_lower in ("*10/*10", "10/10"):
            return 1.1
        if diplotype_lower in ("*4/*14", "*14/*4", "4/14", "14/4"):
            return 0.75
        if diplotype_lower in ("*10/*14", "*14/*10", "10/14", "14/10"):
            return 0.8
        if diplotype_lower in ("*14/*14", "14/14", "slow", "sa"):
            return 0.5

    if gene_upper == "CYP2D6":
        if diplotype_lower in ("*1/*1", "*1/*2", "*2/*2"):
            return 1.0
        if diplotype_lower in ("*1/*4", "*1/*5", "*2/*4", "*10/*10", "im"):
            return 0.5
        if diplotype_lower in ("*4/*4", "*5/*5", "*4/*5", "pm", "poor"):
            return 0.0
        if "x2" in diplotype_lower or diplotype_lower in ("um", "ultrarapid"):
            return 2.0

    # Special cases
    if gene_upper == "ALDH2" and diplotype in ("*1/*2", "heterozygote"):
        return float(special["ALDH2_star1_star2"]["activity_fraction"])

    if gene_upper == "ALDH2" and diplotype in ("*2/*2", "PM_ALDH2"):
        return float(
            special.get("ALDH2_star2_homozygous", {}).get("activity_fraction", 0.001)
        )

    if gene_upper == "GSTM1" and diplotype_lower in (
        "null",
        "null/null",
        "deletion",
        "deleted",
        "0",
        "0/0",
    ):
        return _GST_NULL_RESIDUAL_ACTIVITY

    if gene_upper == "GSTT1" and diplotype_lower in (
        "null",
        "null/null",
        "deletion",
        "deleted",
        "0",
        "0/0",
    ):
        return _GST_NULL_RESIDUAL_ACTIVITY

    if gene_upper in {"GSTM1", "GSTT1"} and diplotype_lower in (
        "present",
        "active",
        "wt",
        "wildtype",
        "*1/*1",
        "1/1",
    ):
        return 1.0

    # Extension cohort aliases.
    if gene_upper == "CYP2E1" and diplotype_lower in ("um_c1c1", "*1c/*1c", "c1/c1"):
        # Source interaction model defines CYP2E1*1C/*1C as a 140% activity state.
        return 1.4

    if gene_upper == "CYP1A1":
        if diplotype_lower in ("*1/*2a", "*1/2a", "wt/*2a", "*2a carrier"):
            # CYP1A1*2A is an inducibility/risk allele rather than a well-calibrated
            # kinetic phenotype, so use a conservative step-up above NM.
            return 1.25
        if diplotype_lower in ("*2a/*2a", "2a/2a"):
            return 1.5

    # Gene-specific phenotype scales (override the generic PM=0/IM=0.5/NM=1.0
    # standard for genes whose null/slow phenotype retains substantial residual
    # activity in vivo).
    if gene_upper == "EPHX1":
        # Hassett 1994; Smith 1997: Y113H/Y113H ("slow") retains ~30-50% epoxide
        # hydrolase activity, not zero.
        if diplotype_lower in ("pm", "slow"):
            return 0.4
        if diplotype_lower in ("im", "intermediate"):
            return 0.7
        if diplotype_lower in ("rm", "rapid", "fast"):
            return 1.3

    if gene_upper == "NQO1":
        # Siegel 1999; Ross 2004: NQO1*2 (Pro187Ser) homozygotes retain ~3-5%
        # activity due to ubiquitin-mediated degradation; heterozygotes ~50%.
        if diplotype_lower in ("pm", "*2/*2"):
            return 0.05
        if diplotype_lower in ("im", "*1/*2"):
            return 0.5

    if gene_upper == "GSTP1":
        # Watson 1998; Hu 1997: Ile105Val (Val/Val) retains ~30-50% activity
        # for many PAH-diol epoxide substrates (reduced thermal stability and
        # affinity, not loss of function).
        if diplotype_lower in ("pm", "val/val"):
            return 0.4
        if diplotype_lower in ("im", "ile/val"):
            return 0.7

    if gene_upper == "CYP1B1":
        # Bailey 1998; Shimada 1999: CYP1B1*3 (L432V) carriers have modestly
        # elevated catalysis (~25-50%), not the standard 2x ultrarapid scale.
        if diplotype_lower in ("rm", "*1/*3", "leu/val"):
            return 1.25
        if diplotype_lower in ("um", "*3/*3", "val/val"):
            return 1.5

    # Standard phenotype scale
    phenotype_map = {
        "pm": std["PM"],
        "poor": std["PM"],
        "im": std["IM"],
        "intermediate": std["IM"],
        "nm": std["NM"],
        "normal": std["NM"],
        "wt": std["NM"],
        "wildtype": std["NM"],
        "*1/*1": std["NM"],
        "rm": std["RM"],
        "rapid": std["RM"],
        "um": std["UM"],
        "ultrarapid": std["UM"],
        "null": std["PM"],
        "0/0": 0.0,
    }

    key = diplotype_lower
    if key in phenotype_map:
        return float(phenotype_map[key])

    # Try numeric
    try:
        val = float(diplotype)
        if 0.0 <= val <= 2.0:
            return val
    except (ValueError, TypeError):
        pass

    warnings.warn(
        f"Unrecognized diplotype '{diplotype}' for gene '{gene}'; defaulting to NM (1.0)",
        stacklevel=2,
    )
    return 1.0


def _proxy_genotype_modifier(diplotype: str, gene: str | None) -> float:
    """Return a silent genotype scaling factor for proxy-model terms."""
    if not gene:
        return 1.0

    gene_upper = gene.upper().strip()
    diplotype_lower = str(diplotype).lower().strip()
    params = _load_kinetic_params()
    special = params["genotype_modifiers"]["special_cases"]
    std = params["genotype_modifiers"]["standard_scale"]

    if gene_upper == "CYP1A2":
        if diplotype_lower in ("*1f/*1f", "1f/1f", "cyp1a2*1f/*1f", "um_1f_1f"):
            return 1.5
        if diplotype_lower in ("*1a/*1f", "1a/1f", "*1f/*1a", "1f/1a"):
            return 1.25
        if diplotype_lower in ("*1a/*1a", "1a/1a"):
            return 1.0
        if diplotype_lower in ("*1k", "*1k/*1k", "1k/1k"):
            return 0.5
        if diplotype_lower in ("pm", "poor", "poor metabolizer", "poor_metabolizer"):
            return 0.3

    if gene_upper == "NAT2":
        if diplotype_lower in ("slow", "slow acetylator", "slow_acetylator", "sa"):
            return 0.2
        if diplotype_lower in ("intermediate", "intermediate acetylator", "intermediate_acetylator"):
            return 0.5
        if diplotype_lower in ("rapid", "rapid acetylator", "rapid_acetylator", "ra"):
            return 1.0

    if gene_upper == "NAT1":
        if diplotype_lower in ("*4/*4", "4/4", "rapid", "ra"):
            return 1.0
        if diplotype_lower in ("*4/*10", "*10/*4", "4/10", "10/4"):
            return 1.05
        if diplotype_lower in ("*10/*10", "10/10"):
            return 1.1
        if diplotype_lower in ("*4/*14", "*14/*4", "4/14", "14/4"):
            return 0.75
        if diplotype_lower in ("*10/*14", "*14/*10", "10/14", "14/10"):
            return 0.8
        if diplotype_lower in ("*14/*14", "14/14", "slow", "sa"):
            return 0.5

    if gene_upper == "CYP2D6":
        if diplotype_lower in ("*1/*1", "*1/*2", "*2/*2"):
            return 1.0
        if diplotype_lower in ("*1/*4", "*1/*5", "*2/*4", "*10/*10", "im"):
            return 0.5
        if diplotype_lower in ("*4/*4", "*5/*5", "*4/*5", "pm", "poor"):
            return 0.0
        if "x2" in diplotype_lower or diplotype_lower in ("um", "ultrarapid"):
            return 2.0

    if gene_upper == "ALDH2" and diplotype in ("*1/*2", "heterozygote"):
        return float(special["ALDH2_star1_star2"]["activity_fraction"])
    if gene_upper == "ALDH2" and diplotype in ("*2/*2", "PM_ALDH2"):
        return float(
            special.get("ALDH2_star2_homozygous", {}).get("activity_fraction", 0.001)
        )
    if gene_upper in {"GSTM1", "GSTT1"} and diplotype_lower in (
        "null",
        "null/null",
        "deletion",
        "deleted",
        "0",
        "0/0",
    ):
        return _GST_NULL_RESIDUAL_ACTIVITY
    if gene_upper in {"GSTM1", "GSTT1"} and diplotype_lower in (
        "present",
        "active",
        "wt",
        "wildtype",
        "*1/*1",
        "1/1",
    ):
        return 1.0
    if gene_upper == "CYP2E1" and diplotype_lower in ("um_c1c1", "*1c/*1c", "c1/c1"):
        return 1.4
    if gene_upper == "CYP1A1":
        if diplotype_lower in ("*1/*2a", "*1/2a", "wt/*2a", "*2a carrier"):
            return 1.25
        if diplotype_lower in ("*2a/*2a", "2a/2a"):
            return 1.5

    if gene_upper == "EPHX1":
        if diplotype_lower in ("pm", "slow"):
            return 0.4
        if diplotype_lower in ("im", "intermediate"):
            return 0.7
        if diplotype_lower in ("rm", "rapid", "fast"):
            return 1.3
    if gene_upper == "NQO1":
        if diplotype_lower in ("pm", "*2/*2"):
            return 0.05
        if diplotype_lower in ("im", "*1/*2"):
            return 0.5
    if gene_upper == "GSTP1":
        if diplotype_lower in ("pm", "val/val"):
            return 0.4
        if diplotype_lower in ("im", "ile/val"):
            return 0.7
    if gene_upper == "CYP1B1":
        if diplotype_lower in ("rm", "*1/*3", "leu/val"):
            return 1.25
        if diplotype_lower in ("um", "*3/*3", "val/val"):
            return 1.5

    phenotype_map = {
        "pm": std["PM"],
        "poor": std["PM"],
        "im": std["IM"],
        "intermediate": std["IM"],
        "nm": std["NM"],
        "normal": std["NM"],
        "wt": std["NM"],
        "wildtype": std["NM"],
        "*1/*1": std["NM"],
        "rm": std["RM"],
        "rapid": std["RM"],
        "um": std["UM"],
        "ultrarapid": std["UM"],
        "null": std["PM"],
        "0/0": 0.0,
    }
    if diplotype_lower in phenotype_map:
        return float(phenotype_map[diplotype_lower])
    try:
        val = float(diplotype)
        if 0.0 <= val <= 40.0:
            return val
    except (ValueError, TypeError):
        pass
    return 1.0


# ── Tissue weight (GTEx integration) ──────────────────────────────────────

_FLUX_TISSUE_TO_GTEX: dict[str, str] = {
    "liver": "Liver",
    "lung": "Lung",
    "prostate": "Prostate",
    "bladder": "Bladder",
    "colon": "Colon",
    "breast": "Breast",
    "kidney": "Kidney",
    "esophagus": "Esophagus",
}

_TISSUE_ALIASES: dict[str, str] = {
    "liver": "liver",
    "hepatic": "liver",
    "lung": "lung",
    "pulmonary": "lung",
    "prostate": "prostate",
    "breast": "breast",
    "colon": "colon",
    "colorectal": "colon",
    "kidney": "kidney",
    "renal": "kidney",
    "bladder": "bladder",
    "lymphocyte": "lymphocyte",
    "blood": "lymphocyte",
    "esophagus": "esophagus",
    "esophageal": "esophagus",
    "stomach": "stomach",
    "gastric": "stomach",
    "intestine": "intestine",
    "nasal_mucosa": "nasal_mucosa",
    "brain": "brain",
    "heart": "heart",
    "muscle": "muscle",
    "adipose": "adipose",
    "skin": "skin",
    "placenta": "placenta",
}


def _normalize_tissue(tissue: str) -> str:
    """Normalize a tissue name string to a canonical key."""
    return _TISSUE_ALIASES.get(tissue.lower().strip(), tissue.lower().strip())


def _normalize_tissue_weight_source(
    tissue_weight_source: FluxTissueWeightSource | str,
) -> FluxTissueWeightSource:
    """Normalize tissue-weight source labels to a supported enum."""
    if isinstance(tissue_weight_source, FluxTissueWeightSource):
        return tissue_weight_source

    normalized = str(tissue_weight_source).strip().lower()
    if normalized == FluxTissueWeightSource.CURATED.value:
        return FluxTissueWeightSource.CURATED
    if normalized == FluxTissueWeightSource.GTEX.value:
        return FluxTissueWeightSource.GTEX
    raise ValueError(
        f"Unknown tissue_weight_source '{tissue_weight_source}'. "
        "Expected 'curated' or 'gtex'."
    )


def get_flux_tissue_weight(
    gene: str,
    tissue: str,
    tissue_weight_source: FluxTissueWeightSource | str = FluxTissueWeightSource.CURATED,
) -> float:
    """Return tissue expression weight for an enzyme in a given tissue.

    The default ``curated`` mode preserves the original ``04_flux_model``
    scientific calibration, which was developed and validated against the
    curated weights stored in ``kinetic_parameters.json``. The optional
    ``gtex`` mode uses quantitative GTEx weights where available, and falls
    back to the curated table for non-GTEx tissues or genes missing from GTEx.

    Args:
        gene: Gene/enzyme symbol (e.g. "CYP1A1").
        tissue: Tissue name (e.g. "Lung", "liver").
        tissue_weight_source: ``"curated"`` (default) or ``"gtex"``.

    Returns:
        Expression weight between 0.0 and 1.0.  Returns 0.5 if the gene
        is not found in any source (moderate default).
    """
    source = _normalize_tissue_weight_source(tissue_weight_source)
    tissue_key = _normalize_tissue(tissue)

    # Explicit GTEx mode: use quantitative GTEx data when available.
    if source == FluxTissueWeightSource.GTEX:
        gtex_name = _FLUX_TISSUE_TO_GTEX.get(tissue_key)
        if gtex_name is not None:
            try:
                from .tissue_subgraphs import get_tissue_weights

                weights = get_tissue_weights(gtex_name)
                if gene in weights:
                    return float(weights[gene])
            except (ImportError, FileNotFoundError, ValueError):
                pass

    # Curated source model, used by the original standalone flux extension.
    params = _load_kinetic_params()
    tw = cast(JsonDict, params["tissue_expression_weights"])

    if gene in tw:
        gene_weights = cast(JsonDict, tw[gene])
        if tissue_key in gene_weights:
            return float(gene_weights[tissue_key])
        # Partial match
        for k, v in gene_weights.items():
            if tissue_key in k or k in tissue_key:
                return float(v)
        # Gene exists but tissue not listed
        return 0.2
    else:
        # Gene not in table — moderate default
        return 0.5


def tissue_weight(
    gene: str,
    tissue: str,
    tissue_weight_source: FluxTissueWeightSource | str = FluxTissueWeightSource.CURATED,
) -> float:
    """Source-compatible alias for :func:`get_flux_tissue_weight`."""
    return get_flux_tissue_weight(gene, tissue, tissue_weight_source=tissue_weight_source)


# ── Risk classification ────────────────────────────────────────────────────


def classify_risk(net_ratio: float) -> RiskClassification:
    """Classify an activation/detoxification net ratio into a risk tier.

    Thresholds:
        <0.5 PROTECTIVE, <1.0 LOW, <2.0 MODERATE, <5.0 ELEVATED, >=5.0 HIGH
    """
    if net_ratio < 0.5:
        return RiskClassification.PROTECTIVE
    elif net_ratio < 1.0:
        return RiskClassification.LOW
    elif net_ratio < 2.0:
        return RiskClassification.MODERATE
    elif net_ratio < 5.0:
        return RiskClassification.ELEVATED
    else:
        return RiskClassification.HIGH


def _classify_risk(net_ratio: float) -> str:
    """Source-compatible string wrapper around :func:`classify_risk`."""
    return classify_risk(net_ratio).value


# ── Helper ─────────────────────────────────────────────────────────────────


def _get_default_concentration(carcinogen_class: str) -> float:
    """Return default environmental exposure concentration in uM."""
    params = _load_kinetic_params()
    defaults = params["metadata"]["exposure_defaults_uM"]
    class_map = {
        "PAH": defaults["BaP"],
        "Aflatoxin": defaults["AFB1"],
        "Aldehyde": defaults["acetaldehyde"],
        "Nitrosamine": defaults["NNK"],
        "NDMA": defaults["NDMA"],
        "HCA": defaults["PhIP"],
        "Benzene": defaults["benzene"],
        "ChlorinatedSolvent": defaults["TCE"],
        "Dioxin": defaults["TCDD"],
        "HeavyMetal": defaults["arsenic"],
    }
    if carcinogen_class in class_map:
        return float(class_map[carcinogen_class])

    proxy_classes = _load_proxy_flux_params()["classes"]
    proxy_cfg = proxy_classes.get(carcinogen_class)
    if proxy_cfg is not None:
        source = proxy_cfg["exposure_default"]["source"]
        if source == "kinetic_parameters":
            field = proxy_cfg["exposure_default"]["field"]
            value = defaults.get(field)
            if value is not None:
                return float(value)
        elif source == "exposure_database":
            exposure_db = _load_exposure_db()["carcinogen_classes"]
            class_name = proxy_cfg["exposure_default"]["class"]
            scenario_name = proxy_cfg["exposure_default"]["scenario"]
            field_name = proxy_cfg["exposure_default"]["field"]
            scenario = exposure_db.get(class_name, {}).get("exposure_scenarios", {}).get(scenario_name, {})
            value = scenario.get(field_name)
            if value is not None:
                return float(value)

    return 0.1


_KNOWN_FLUX_GENE_PREFIXES = tuple(
    sorted(
        {
            "ABCB1", "ABCC2", "ABCG2", "ADH1B", "ADH5", "AHRR", "ALDH1A1", "ALDH2",
            "AS3MT", "COMT", "CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2A13",
            "CYP2D6", "CYP2E1", "CYP3A4", "EPHX1", "ERCC2", "GSTA1", "GSTM1",
            "GSTP1", "GSTT1", "MGMT", "NAT1", "NAT2", "NQO1", "OGG1", "POLH",
            "SULT1E1", "UGT2B7", "XPC", "XRCC1",
        },
        key=len,
        reverse=True,
    )
)


def _term_gene_name(term_name: str) -> str | None:
    """Resolve a result term name such as ``CYP3A4_AFQ1`` to a gene symbol."""
    for gene in _KNOWN_FLUX_GENE_PREFIXES:
        if term_name == gene or term_name.startswith(f"{gene}_"):
            return gene
    return None


def _resolve_induction_factors(
    lifestyle: LifestyleMap | None = None,
    induction_factors: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Resolve optional co-exposure induction inputs into per-enzyme Vmax folds."""
    resolved: dict[str, float] = {}

    if lifestyle:
        try:
            from .interaction_engine import enzyme_induction_modifier

            resolved.update(enzyme_induction_modifier(lifestyle).enzyme_folds)
        except Exception as exc:
            warnings.warn(
                f"Could not resolve lifestyle induction factors; using explicit/default factors only: {exc}",
                stacklevel=2,
            )

    if induction_factors:
        for gene, factor in induction_factors.items():
            try:
                numeric = float(factor)
            except (TypeError, ValueError):
                continue
            if numeric > 0:
                resolved[str(gene).upper()] = numeric

    return {
        gene: round(factor, 6)
        for gene, factor in sorted(resolved.items())
        if math.isfinite(factor) and factor > 0 and not math.isclose(factor, 1.0)
    }


def _rescale_flux_section_for_induction(
    enzymes: dict[str, Any],
    induction_factors: Mapping[str, float],
) -> tuple[float, float]:
    """Apply induction folds to enzyme-term fluxes and return old/new sums."""
    old_sum = 0.0
    new_sum = 0.0
    for term_name, edata in enzymes.items():
        if not isinstance(edata, dict):
            continue
        try:
            old_flux = float(edata.get("flux", 0.0))
        except (TypeError, ValueError):
            old_flux = 0.0
        gene = _term_gene_name(term_name)
        factor = float(induction_factors.get(gene or "", 1.0))
        new_flux = old_flux * factor
        edata["induction_modifier"] = round(factor, 6)
        if not math.isclose(factor, 1.0):
            edata["flux"] = round(new_flux, 6)
        old_sum += old_flux
        new_sum += new_flux
    return old_sum, new_sum


def _apply_induction_modifiers(
    result: FluxResultDict,
    induction_factors: Mapping[str, float],
) -> FluxResultDict:
    """Apply resolved Vmax induction folds to an internal flux-result payload."""
    if not induction_factors:
        return result

    for section_name, total_name in (
        ("activation_enzymes", "total_activation"),
        ("detox_enzymes", "total_detox"),
    ):
        enzymes = result.get(section_name, {})
        if not isinstance(enzymes, dict):
            continue
        old_sum, new_sum = _rescale_flux_section_for_induction(enzymes, induction_factors)
        if old_sum > 0 and total_name in result:
            result[total_name] = round(float(result[total_name]) * new_sum / old_sum, 6)

    return result


_FALLBACK_QIVIVE_TISSUES: dict[str, dict[str, float]] = {
    "liver": {"mppgl_mg_per_g": 40.0, "organ_weight_g": 1500.0},
    "lung": {"mppgl_mg_per_g": 20.0, "organ_weight_g": 1000.0},
    "kidney": {"mppgl_mg_per_g": 12.0, "organ_weight_g": 300.0},
    "intestine": {"mppgl_mg_per_g": 35.0, "organ_weight_g": 900.0},
}

_STEADY_STATE_DEFAULTS: dict[str, float] = {
    "body_weight_kg": 70.0,
    "volume_l_per_kg": 0.7,
    "absorption_fraction": 1.0,
    "exposure_frequency_per_day": 1.0,
    "cardiac_output_l_per_day": 7200.0,
    "background_clearance_rate_per_day": 0.05,
    "reactive_intermediate_loss_rate_per_day": 1.0,
    "detoxified_metabolite_loss_rate_per_day": 1.0,
    "flux_rate_scale_per_day": 1.0,
    "rate_reference_concentration_uM": 1.0,
}

_FALLBACK_STEADY_STATE_TISSUES: dict[str, dict[str, float]] = {
    "liver": {
        "organ_weight_g": 1500.0,
        "tissue_partition_coefficient": 1.0,
        "tissue_blood_flow_fraction": 0.25,
    },
    "lung": {
        "organ_weight_g": 1000.0,
        "tissue_partition_coefficient": 0.8,
        "tissue_blood_flow_fraction": 1.0,
    },
    "kidney": {
        "organ_weight_g": 300.0,
        "tissue_partition_coefficient": 1.1,
        "tissue_blood_flow_fraction": 0.2,
    },
    "intestine": {
        "organ_weight_g": 900.0,
        "tissue_partition_coefficient": 0.9,
        "tissue_blood_flow_fraction": 0.12,
    },
    "bladder": {
        "organ_weight_g": 150.0,
        "tissue_partition_coefficient": 0.7,
        "tissue_blood_flow_fraction": 0.02,
    },
    "breast": {
        "organ_weight_g": 500.0,
        "tissue_partition_coefficient": 1.4,
        "tissue_blood_flow_fraction": 0.03,
    },
    "colon": {
        "organ_weight_g": 600.0,
        "tissue_partition_coefficient": 0.9,
        "tissue_blood_flow_fraction": 0.08,
    },
    "prostate": {
        "organ_weight_g": 30.0,
        "tissue_partition_coefficient": 0.8,
        "tissue_blood_flow_fraction": 0.01,
    },
    "esophagus": {
        "organ_weight_g": 40.0,
        "tissue_partition_coefficient": 0.8,
        "tissue_blood_flow_fraction": 0.01,
    },
    "skin": {
        "organ_weight_g": 3300.0,
        "tissue_partition_coefficient": 1.2,
        "tissue_blood_flow_fraction": 0.05,
    },
}


def qivive_intrinsic_clearance(
    vmax: float,
    km: float,
    *,
    microsomal_protein_mg_per_g_tissue: float,
    organ_weight_g: float,
) -> float:
    """Upscale in vitro intrinsic clearance using MPPGL and organ weight.

    The returned value preserves the caller's Vmax/Km unit family, multiplied by
    mg microsomal protein per gram tissue and organ mass.
    """
    if km <= 0:
        raise ValueError(f"Km must be positive, got {km}")
    if microsomal_protein_mg_per_g_tissue <= 0:
        raise ValueError("microsomal_protein_mg_per_g_tissue must be positive")
    if organ_weight_g <= 0:
        raise ValueError("organ_weight_g must be positive")
    return (vmax / km) * microsomal_protein_mg_per_g_tissue * organ_weight_g


def _qivive_context_for_tissue(
    tissue: str,
    overrides: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Return MPPGL/organ-weight context for optional QIVIVE flux scaling."""
    params = _load_kinetic_params()
    metadata = cast(JsonDict, params.get("metadata", {}))
    qivive_defaults = cast(JsonDict, metadata.get("qivive_defaults", {}))
    tissue_defaults = cast(JsonDict, qivive_defaults.get("tissues", {}))
    tissue_key = _normalize_tissue(tissue)
    source = cast(dict[str, float], tissue_defaults.get(tissue_key, _FALLBACK_QIVIVE_TISSUES.get(tissue_key, _FALLBACK_QIVIVE_TISSUES["liver"])))
    context = {
        "mppgl_mg_per_g": float(source["mppgl_mg_per_g"]),
        "organ_weight_g": float(source["organ_weight_g"]),
    }
    if overrides:
        if "mppgl_mg_per_g" in overrides:
            context["mppgl_mg_per_g"] = float(overrides["mppgl_mg_per_g"])
        if "organ_weight_g" in overrides:
            context["organ_weight_g"] = float(overrides["organ_weight_g"])
    context["scale"] = round(context["mppgl_mg_per_g"] * context["organ_weight_g"], 6)
    return context


def _apply_qivive_scale(result: FluxResultDict, qivive_context: Mapping[str, float]) -> FluxResultDict:
    """Apply a common tissue-level QIVIVE scale to reported flux magnitudes."""
    scale = float(qivive_context.get("scale", 1.0))
    if math.isclose(scale, 1.0):
        return result

    for section_name in ("activation_enzymes", "detox_enzymes"):
        enzymes = result.get(section_name, {})
        if not isinstance(enzymes, dict):
            continue
        for edata in enzymes.values():
            if not isinstance(edata, dict):
                continue
            edata["qivive_scale"] = round(scale, 6)
            if "flux" in edata:
                edata["flux"] = round(float(edata["flux"]) * scale, 6)

    for key in ("total_activation", "total_detox"):
        if key in result:
            result[key] = round(float(result[key]) * scale, 6)

    note = str(result.get("unit_note", "")).strip()
    qivive_note = (
        "QIVIVE common tissue scale applied using MPPGL "
        f"{qivive_context['mppgl_mg_per_g']} mg/g and organ weight "
        f"{qivive_context['organ_weight_g']} g."
    )
    result["unit_note"] = f"{note}; {qivive_note}" if note else qivive_note
    return result


def _susceptibility_score_log2(net_ratio: float) -> float:
    """Return log2 activation/detoxification susceptibility score."""
    if net_ratio <= 0 or not math.isfinite(net_ratio):
        return 0.0
    return round(math.log2(net_ratio), 4)


def _positive_context_float(context: Mapping[str, Any], key: str, fallback: float) -> float:
    """Read a positive numeric context value with a conservative fallback."""
    try:
        value = float(context.get(key, fallback))
    except (TypeError, ValueError):
        return fallback
    if value <= 0 or not math.isfinite(value):
        return fallback
    return value


def _bounded_fraction_context(context: Mapping[str, Any], key: str, fallback: float) -> float:
    """Read a fraction constrained to the open interval used by PBPK rates."""
    try:
        value = float(context.get(key, fallback))
    except (TypeError, ValueError):
        return fallback
    if not math.isfinite(value):
        return fallback
    return min(max(value, 1e-6), 1.0)


def _round_steady_state_value(value: float) -> float:
    """Round steady-state outputs without hiding very small non-zero values."""
    if not math.isfinite(value) or value < 0:
        return 0.0
    if value == 0:
        return 0.0
    return round(value, 8 if abs(value) < 1e-4 else 6)


def _steady_state_context_for_tissue(
    tissue: str,
    overrides: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Return validated defaults for the flux-coupled steady-state solver."""
    params = _load_kinetic_params()
    metadata = cast(JsonDict, params.get("metadata", {}))
    configured = cast(JsonDict, metadata.get("steady_state_defaults", {}))
    tissue_key = _normalize_tissue(tissue)

    qivive_context = _qivive_context_for_tissue(tissue)
    fallback_tissue = _FALLBACK_STEADY_STATE_TISSUES.get(
        tissue_key,
        _FALLBACK_STEADY_STATE_TISSUES["liver"],
    )
    configured_tissues = cast(JsonDict, configured.get("tissues", {}))
    configured_tissue = cast(JsonDict, configured_tissues.get(tissue_key, {}))
    tissue_source: dict[str, Any] = {
        **fallback_tissue,
        **configured_tissue,
        "organ_weight_g": configured_tissue.get(
            "organ_weight_g",
            qivive_context.get("organ_weight_g", fallback_tissue["organ_weight_g"]),
        ),
    }

    context: dict[str, float] = {}
    for key, fallback in _STEADY_STATE_DEFAULTS.items():
        context[key] = _positive_context_float(configured, key, fallback)
    context["absorption_fraction"] = _bounded_fraction_context(
        configured,
        "absorption_fraction",
        _STEADY_STATE_DEFAULTS["absorption_fraction"],
    )
    context["organ_weight_g"] = _positive_context_float(
        tissue_source,
        "organ_weight_g",
        fallback_tissue["organ_weight_g"],
    )
    context["tissue_partition_coefficient"] = _positive_context_float(
        tissue_source,
        "tissue_partition_coefficient",
        fallback_tissue["tissue_partition_coefficient"],
    )
    context["tissue_blood_flow_fraction"] = _bounded_fraction_context(
        tissue_source,
        "tissue_blood_flow_fraction",
        fallback_tissue["tissue_blood_flow_fraction"],
    )

    if overrides:
        for key, value in overrides.items():
            if key in {"absorption_fraction", "tissue_blood_flow_fraction"}:
                context[key] = _bounded_fraction_context(overrides, key, context[key])
            else:
                context[key] = _positive_context_float(overrides, key, context.get(key, 1.0))

    context["central_volume_l"] = round(
        context["body_weight_kg"] * context["volume_l_per_kg"],
        6,
    )
    context["tissue_volume_l"] = round(context["organ_weight_g"] / 1000.0, 6)
    context["tissue_blood_flow_l_per_day"] = round(
        context["cardiac_output_l_per_day"] * context["tissue_blood_flow_fraction"],
        6,
    )
    return context


def solve_flux_steady_state(
    substrate_conc_uM: float,
    activation_flux: float,
    detox_flux: float,
    tissue: str,
    *,
    context: Mapping[str, float] | None = None,
) -> FluxSteadyStateResult:
    """Solve a one-tissue PBPK steady-state model coupled to pathway flux.

    The solver treats activation and detoxification fluxes as concentration-
    normalized first-order metabolic rate constants, then solves central and
    tissue steady state with perfusion-limited tissue extraction. This is more
    explicit than the former proportional proxy: every reported concentration
    is derived from volume, organ mass, blood flow, partitioning, and clearance
    rates carried in the returned model payload.
    """
    if substrate_conc_uM < 0:
        raise ValueError("substrate_conc_uM cannot be negative")
    if activation_flux < 0:
        raise ValueError("activation_flux cannot be negative")
    if detox_flux < 0:
        raise ValueError("detox_flux cannot be negative")

    solver_context = _steady_state_context_for_tissue(tissue, context)
    central_volume_l = solver_context["central_volume_l"]
    tissue_volume_l = solver_context["tissue_volume_l"]
    tissue_flow_l_per_day = solver_context["tissue_blood_flow_l_per_day"]
    partition = solver_context["tissue_partition_coefficient"]
    background_clearance_rate = solver_context["background_clearance_rate_per_day"]
    flux_rate_scale = solver_context["flux_rate_scale_per_day"]
    reference_conc = max(
        substrate_conc_uM,
        solver_context["rate_reference_concentration_uM"],
        1e-12,
    )

    activation_rate = max(activation_flux, 0.0) / reference_conc * flux_rate_scale
    detox_rate = max(detox_flux, 0.0) / reference_conc * flux_rate_scale
    metabolic_rate = activation_rate + detox_rate
    intrinsic_clearance_l_per_day = metabolic_rate * tissue_volume_l
    extraction_ratio = (
        intrinsic_clearance_l_per_day / (tissue_flow_l_per_day + intrinsic_clearance_l_per_day)
        if tissue_flow_l_per_day + intrinsic_clearance_l_per_day > 0
        else 0.0
    )
    tissue_clearance_l_per_day = tissue_flow_l_per_day * extraction_ratio
    background_clearance_l_per_day = background_clearance_rate * central_volume_l
    total_clearance_l_per_day = background_clearance_l_per_day + tissue_clearance_l_per_day

    input_rate_umol_per_day = (
        substrate_conc_uM
        * central_volume_l
        * solver_context["absorption_fraction"]
        * solver_context["exposure_frequency_per_day"]
    )
    central_conc = (
        input_rate_umol_per_day / total_clearance_l_per_day
        if total_clearance_l_per_day > 0
        else 0.0
    )
    tissue_conc = (
        partition * central_conc * tissue_flow_l_per_day
        / (tissue_flow_l_per_day + intrinsic_clearance_l_per_day)
        if tissue_flow_l_per_day + intrinsic_clearance_l_per_day > 0
        else partition * central_conc
    )
    reactive_loss_rate = solver_context["reactive_intermediate_loss_rate_per_day"] + detox_rate
    detoxified_loss_rate = solver_context["detoxified_metabolite_loss_rate_per_day"]
    reactive_conc = (
        tissue_conc * activation_rate / reactive_loss_rate
        if reactive_loss_rate > 0
        else 0.0
    )
    detoxified_conc = (
        tissue_conc * detox_rate / detoxified_loss_rate
        if detoxified_loss_rate > 0
        else 0.0
    )

    central_rate = total_clearance_l_per_day / central_volume_l if central_volume_l > 0 else 0.0
    tissue_exchange_rate = (
        tissue_flow_l_per_day / (tissue_volume_l * partition) + metabolic_rate
        if tissue_volume_l > 0 and partition > 0
        else metabolic_rate
    )
    steady_rates = [
        rate
        for rate in (
            central_rate,
            tissue_exchange_rate,
            reactive_loss_rate,
            detoxified_loss_rate,
        )
        if rate > 0 and math.isfinite(rate)
    ]
    time_to_steady_state_days = 4.0 / min(steady_rates) if steady_rates else 0.0

    concentrations = {
        "central_substrate_uM": _round_steady_state_value(central_conc),
        "tissue_substrate_uM": _round_steady_state_value(tissue_conc),
        "reactive_intermediate_uM": _round_steady_state_value(reactive_conc),
        "detoxified_metabolite_uM": _round_steady_state_value(detoxified_conc),
    }
    model = {
        "model": "one_tissue_perfusion_limited_pbpk_steady_state",
        "input_rate_umol_per_day": _round_steady_state_value(input_rate_umol_per_day),
        "activation_rate_per_day": _round_steady_state_value(activation_rate),
        "detox_rate_per_day": _round_steady_state_value(detox_rate),
        "metabolic_rate_per_day": _round_steady_state_value(metabolic_rate),
        "background_clearance_l_per_day": _round_steady_state_value(
            background_clearance_l_per_day
        ),
        "tissue_clearance_l_per_day": _round_steady_state_value(tissue_clearance_l_per_day),
        "total_clearance_l_per_day": _round_steady_state_value(total_clearance_l_per_day),
        "extraction_ratio": _round_steady_state_value(extraction_ratio),
        "time_to_steady_state_days": _round_steady_state_value(time_to_steady_state_days),
        **{
            key: _round_steady_state_value(value)
            for key, value in solver_context.items()
        },
    }
    return FluxSteadyStateResult(concentrations_uM=concentrations, model=model)


def _steady_state_concentration_proxy(
    substrate_conc_uM: float,
    act: float,
    det: float,
    tissue: str = "Liver",
    context: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Deprecated compatibility alias for historical proxy-shaped payloads."""
    steady_state = solve_flux_steady_state(substrate_conc_uM, act, det, tissue, context=context)
    return {
        "reactive_intermediate_proxy_uM": steady_state.concentrations_uM[
            "reactive_intermediate_uM"
        ],
        "detoxified_metabolite_proxy_uM": steady_state.concentrations_uM[
            "detoxified_metabolite_uM"
        ],
    }


def _relative_capacity_scale(label: str | None) -> float:
    """Map qualitative intrinsic-clearance labels onto coarse numeric scales."""
    if label is None or str(label).strip() == "":
        return 1.0

    mapping = {
        "very_low": 0.1,
        "low": 0.2,
        "low_to_moderate": 0.35,
        "moderate": 0.5,
        "moderate_to_high": 0.75,
        "high": 1.0,
    }
    return mapping.get(str(label or "").strip().lower(), 0.35)


def _pathway_tissue_weight(
    gene: str | None,
    tissue: str,
    tissue_weight_source: FluxTissueWeightSource,
    supported_tissues: list[str] | None = None,
) -> float:
    """Return a tissue weight, falling back to pathway-level tissue support."""
    params = _load_kinetic_params()
    tissue_weights = cast(JsonDict, params["tissue_expression_weights"])
    if gene and gene in tissue_weights:
        return get_flux_tissue_weight(gene, tissue, tissue_weight_source)

    if supported_tissues:
        tissue_key = _normalize_tissue(tissue)
        supported_keys = {_normalize_tissue(name) for name in supported_tissues}
        return 1.0 if tissue_key in supported_keys else 0.2

    if gene:
        return get_flux_tissue_weight(gene, tissue, tissue_weight_source)
    return 0.5


def _proxy_diplotype_for_gene(genotypes: GenotypeMap, gene: str | None) -> str:
    """Return a safe proxy-model diplotype label for an optional gene key."""
    if not gene:
        return "NM"
    return genotypes.get(gene, "NM")


def _compute_proxy_mm_term(
    term: JsonDict,
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> tuple[float, float, float]:
    """Evaluate a Michaelis-Menten proxy term."""
    gene = cast(str | None, term.get("gene"))
    gm = _proxy_genotype_modifier(_proxy_diplotype_for_gene(genotypes, gene), gene)
    tw = _pathway_tissue_weight(
        gene,
        tissue,
        tissue_weight_source,
        supported_tissues=term.get("supported_tissues"),
    )
    vmax = (
        float(term["vmax"])
        * float(term.get("vmax_relative", 1.0))
        * _relative_capacity_scale(term.get("relative_capacity"))
    )
    flux = michaelis_menten(S, vmax * gm * tw, float(term["km"]))
    return flux, gm, tw


def _compute_proxy_saturating_term(
    term: JsonDict,
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> tuple[float, float, float]:
    """Evaluate a simple saturating proxy term."""
    gene = cast(str | None, term.get("gene"))
    gm = _proxy_genotype_modifier(_proxy_diplotype_for_gene(genotypes, gene), gene)
    tw = _pathway_tissue_weight(
        gene,
        tissue,
        tissue_weight_source,
        supported_tissues=term.get("supported_tissues"),
    )
    flux = float(term["scale"]) * gm * tw * S / (float(term["km"]) + S)
    return flux, gm, tw


def _compute_proxy_hill_term(
    term: JsonDict,
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> tuple[float, float, float]:
    """Evaluate a Hill-style damage or signaling proxy term."""
    gene = cast(str | None, term.get("gene"))
    gm = _proxy_genotype_modifier(_proxy_diplotype_for_gene(genotypes, gene), gene)
    tw = _pathway_tissue_weight(
        gene,
        tissue,
        tissue_weight_source,
        supported_tissues=term.get("supported_tissues"),
    )
    substrate = S * float(term.get("substrate_scale", 1.0))
    flux = hill_equation(
        substrate,
        float(term["vmax"]) * gm * tw,
        float(term["k50"]),
        float(term.get("hill_n", 1.0)),
    )
    return flux, gm, tw


def _compute_proxy_repair_term(
    activation_flux: float,
    term: JsonDict,
    genotypes: GenotypeMap,
    tissue: str,
    tissue_weight_source: FluxTissueWeightSource,
) -> tuple[float, float, float]:
    """Evaluate a repair-capacity proxy term derived from activation burden."""
    gene = cast(str | None, term.get("gene"))
    gm = _proxy_genotype_modifier(_proxy_diplotype_for_gene(genotypes, gene), gene)
    tw = _pathway_tissue_weight(
        gene,
        tissue,
        tissue_weight_source,
        supported_tissues=term.get("supported_tissues"),
    )
    return activation_flux * float(term["scale"]) * gm * tw, gm, tw


def _get_proxy_class_params(class_name: str) -> JsonDict:
    """Return the proxy flux config for a class."""
    proxy_classes = cast(JsonDict, _load_proxy_flux_params()["classes"])
    return cast(JsonDict, proxy_classes[class_name])


def _resolve_nested_ref(doc: JsonDict, ref: str) -> JsonDict:
    """Resolve a dotted reference inside a nested mapping."""
    node: JsonDict = doc
    for part in ref.split("."):
        node = cast(JsonDict, node[part])
    return node


def _class_parameter_metadata(carcinogen_class: str) -> dict[str, str]:
    """Return class-level parameter metadata for measured or proxy models."""
    proxy_cfg = _load_proxy_flux_params()["classes"].get(carcinogen_class)
    if proxy_cfg is None:
        return {
            "model_kind": "measured_kinetics",
            "parameter_source": _KINETIC_PARAMS_FILE.name,
        }
    return {
        "model_kind": proxy_cfg["model_kind"],
        "parameter_source": _PROXY_FLUX_PARAMS_FILE.name,
    }


def _proxy_term_metadata(
    carcinogen_class: str,
    term_name: str,
    *,
    activation_term: bool,
) -> JsonDict | None:
    """Return provenance metadata for a proxy-model term."""
    cfg = _get_proxy_class_params(carcinogen_class)
    sections = ("activation_terms",) if activation_term else ("detox_terms", "repair_terms")

    for section in sections:
        term_cfg = cfg.get(section, {}).get(term_name)
        if term_cfg is None:
            continue

        ref = term_cfg.get("provenance_ref", "")
        sources: list[str] = []
        basis = ""
        if ref:
            entry = _resolve_nested_ref(_load_proxy_flux_provenance(), ref)
            sources = list(entry.get("sources", []))
            basis = entry.get("parameter_basis", "")

        return {
            "model_kind": cfg["model_kind"],
            "parameter_source": _PROXY_FLUX_PARAMS_FILE.name,
            "provenance_ref": ref,
            "provenance_sources": sources,
            "parameter_basis": basis,
        }

    return None


def _annotate_flux_result_metadata(
    carcinogen_class: str,
    result: FluxResultDict,
) -> FluxResultDict:
    """Attach class- and enzyme-level parameter metadata to a flux result."""
    class_meta = _class_parameter_metadata(carcinogen_class)
    result.setdefault("model_kind", class_meta["model_kind"])
    result.setdefault("parameter_source", class_meta["parameter_source"])

    proxy_mode = class_meta["parameter_source"] == _PROXY_FLUX_PARAMS_FILE.name
    for activation_term, enzymes in (
        (True, result.get("activation_enzymes", {})),
        (False, result.get("detox_enzymes", {})),
    ):
        for name, edata in enzymes.items():
            if not isinstance(edata, dict):
                continue

            edata.setdefault("model_kind", class_meta["model_kind"])
            edata.setdefault("parameter_source", class_meta["parameter_source"])
            edata.setdefault("provenance_ref", "")
            edata.setdefault("provenance_sources", [])
            edata.setdefault("parameter_basis", "")

            if not proxy_mode:
                continue

            proxy_meta = _proxy_term_metadata(
                carcinogen_class,
                name,
                activation_term=activation_term,
            )
            if proxy_meta is None:
                continue
            edata.update(proxy_meta)

    return result


# ── Pathway-specific flux calculators ──────────────────────────────────────


def _compute_pah_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute PAH (BaP) activation and detoxification fluxes."""
    params = _load_kinetic_params()
    p = params["carcinogen_classes"]["PAH"]["pathways"]
    activation_enzymes: dict[str, Any] = {}
    detox_enzymes: dict[str, Any] = {}

    # --- ACTIVATION ---
    # CYP1A1
    cyp1a1_p = p["activation"]["CYP1A1"]
    gm = genotype_modifier(genotypes.get("CYP1A1", "NM"), "CYP1A1")
    tw = get_flux_tissue_weight("CYP1A1", tissue, tissue_weight_source)
    v_cyp1a1 = michaelis_menten(
        S, cyp1a1_p["Vmax_pmol_min_pmolP450"] * gm * tw, cyp1a1_p["Km_uM"]
    )
    activation_enzymes["CYP1A1"] = {
        "flux": round(v_cyp1a1, 4),
        "genotype_modifier": gm,
        "tissue_weight": tw,
        "confidence": cyp1a1_p["confidence"],
    }

    # CYP1B1
    cyp1b1_p = p["activation"]["CYP1B1"]
    gm1b1 = genotype_modifier(genotypes.get("CYP1B1", "NM"), "CYP1B1")
    tw1b1 = get_flux_tissue_weight("CYP1B1", tissue, tissue_weight_source)
    v_cyp1b1 = michaelis_menten(
        S, cyp1b1_p["Vmax_pmol_min_pmolP450"] * gm1b1 * tw1b1, cyp1b1_p["Km_uM"]
    )
    activation_enzymes["CYP1B1"] = {
        "flux": round(v_cyp1b1, 4),
        "genotype_modifier": gm1b1,
        "tissue_weight": tw1b1,
        "confidence": cyp1b1_p["confidence"],
    }

    # EPHX1 (intermediate efficiency factor)
    ephx1_p = p["activation"]["EPHX1"]
    ephx1_gm = genotype_modifier(genotypes.get("EPHX1", "NM"), "EPHX1")
    ephx1_tw = get_flux_tissue_weight("EPHX1", tissue, tissue_weight_source)
    ephx1_efficiency = min(
        1.0, (ephx1_p["CLint"] / 40.0) * ephx1_gm * ephx1_tw
    )

    total_activation_raw = v_cyp1a1 + v_cyp1b1
    total_activation = total_activation_raw * ephx1_efficiency

    # --- DETOXIFICATION ---
    # GSTM1
    gstm1_p = p["detoxification"]["GSTM1"]
    gstm1_gm = genotype_modifier(genotypes.get("GSTM1", "NM"), "GSTM1")
    gstm1_tw = get_flux_tissue_weight("GSTM1", tissue, tissue_weight_source)
    v_gstm1 = gstm1_p["kcat_s"] * gstm1_gm * gstm1_tw * S / (gstm1_p["Km_uM"] + S)
    detox_enzymes["GSTM1"] = {
        "flux": round(v_gstm1, 6),
        "genotype_modifier": gstm1_gm,
        "tissue_weight": gstm1_tw,
        "confidence": gstm1_p["confidence"],
    }

    # GSTP1
    gstp1_p = p["detoxification"]["GSTP1"]
    gstp1_gm = genotype_modifier(genotypes.get("GSTP1", "NM"), "GSTP1")
    gstp1_tw = get_flux_tissue_weight("GSTP1", tissue, tissue_weight_source)
    v_gstp1 = gstp1_p["kcat_s"] * gstp1_gm * gstp1_tw * S / (gstp1_p["Km_uM"] + S)
    detox_enzymes["GSTP1"] = {
        "flux": round(v_gstp1, 6),
        "genotype_modifier": gstp1_gm,
        "tissue_weight": gstp1_tw,
        "confidence": gstp1_p["confidence"],
    }

    total_detox_raw = v_gstm1 + v_gstp1
    # CLint normalization: CYP1A1 CLint 1100 / (GSTM1 CLint 0.024 * 60) = 764
    GST_SCALE = 764.0
    total_detox = total_detox_raw * GST_SCALE

    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": detox_enzymes,
        "total_activation": round(total_activation, 4),
        "total_detox": round(total_detox, 4),
        "ephx1_efficiency": round(ephx1_efficiency, 3),
        "unit_note": "Activation in pmol/min/pmolP450; detox scaled via CLint normalization",
    }


def _compute_aflatoxin_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute AFB1 activation and detoxification fluxes."""
    params = _load_kinetic_params()
    p = params["carcinogen_classes"]["Aflatoxin"]["pathways"]
    activation_enzymes: dict[str, Any] = {}
    detox_enzymes: dict[str, Any] = {}

    # CYP3A4 — Hill kinetics
    cyp3a4_p = p["activation"]["CYP3A4"]
    gm3a4 = genotype_modifier(genotypes.get("CYP3A4", "NM"), "CYP3A4")
    tw3a4 = get_flux_tissue_weight("CYP3A4", tissue, tissue_weight_source)
    v_cyp3a4 = hill_equation(
        S,
        cyp3a4_p["Vmax_pmol_min_pmolP450"] * gm3a4 * tw3a4,
        cyp3a4_p["Km_uM"],
        cyp3a4_p["hill_n"],
    )
    activation_enzymes["CYP3A4"] = {
        "flux": round(v_cyp3a4, 4),
        "kinetics": "hill",
        "n": cyp3a4_p["hill_n"],
        "genotype_modifier": gm3a4,
        "tissue_weight": tw3a4,
        "confidence": cyp3a4_p["confidence"],
        "fraction_contribution": 0.45,
    }

    # CYP1A2 — Michaelis-Menten
    cyp1a2_p = p["activation"]["CYP1A2"]
    gm1a2 = genotype_modifier(genotypes.get("CYP1A2", "NM"), "CYP1A2")
    tw1a2 = get_flux_tissue_weight("CYP1A2", tissue, tissue_weight_source)
    v_cyp1a2 = michaelis_menten(
        S,
        cyp1a2_p["Vmax_pmol_min_pmolP450"] * gm1a2 * tw1a2,
        cyp1a2_p["Km_uM"],
    )
    activation_enzymes["CYP1A2"] = {
        "flux": round(v_cyp1a2, 4),
        "kinetics": "michaelis_menten",
        "genotype_modifier": gm1a2,
        "tissue_weight": tw1a2,
        "confidence": cyp1a2_p["confidence"],
        "fraction_contribution": 0.49,
    }

    total_activation = v_cyp3a4 + v_cyp1a2

    # CYP3A4 AFQ1 (detox)
    afq1_p = p["detoxification"]["CYP3A4_AFQ1"]
    v_afq1 = hill_equation(
        S,
        afq1_p["Vmax_pmol_min_pmolP450"] * gm3a4 * tw3a4,
        afq1_p["Km_uM"],
        n=2.0,
    )
    detox_enzymes["CYP3A4_AFQ1"] = {
        "flux": round(v_afq1, 4),
        "genotype_modifier": gm3a4,
        "tissue_weight": tw3a4,
        "confidence": afq1_p["confidence"],
    }

    # GSTA1 (estimated)
    gsta1_gm = genotype_modifier(genotypes.get("GSTA1", "NM"), "GSTA1")
    gsta1_tw = get_flux_tissue_weight("GSTM1", tissue, tissue_weight_source)  # proxy
    GSTA1_VMAX_EST = 0.05
    GSTA1_KM_EST = 10.0
    v_gsta1 = michaelis_menten(S, GSTA1_VMAX_EST * gsta1_gm * gsta1_tw, GSTA1_KM_EST)
    detox_enzymes["GSTA1_conjugation"] = {
        "flux": round(v_gsta1, 6),
        "genotype_modifier": gsta1_gm,
        "tissue_weight": gsta1_tw,
        "confidence": "low",
    }

    total_detox = v_afq1 + v_gsta1

    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": detox_enzymes,
        "total_activation": round(total_activation, 4),
        "total_detox": round(total_detox, 6),
        "unit_note": "All fluxes in pmol/min/pmolP450; GSTA1 Vmax is estimated",
    }


def _compute_aldehyde_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute aldehyde (acetaldehyde) clearance flux."""
    params = _load_kinetic_params()
    p = params["carcinogen_classes"]["Aldehyde"]["pathways"]["acetaldehyde_clearance"]
    detox_enzymes: dict[str, Any] = {}

    # Determine ALDH2 genotype
    aldh2_gt = genotypes.get("ALDH2", "*1/*1")
    if aldh2_gt in ("*1/*1", "NM", "WT", "wildtype"):
        aldh2_p = p["ALDH2_star1"]
        aldh2_gm = 1.0
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]
    elif aldh2_gt in ("*1/*2", "heterozygote"):
        aldh2_p = p["ALDH2_star1"]
        aldh2_gm = genotype_modifier("*1/*2", "ALDH2")  # 0.25
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]
    elif aldh2_gt in ("*2/*2", "PM", "PM_ALDH2"):
        aldh2_p = p["ALDH2_star2_homozygous"]
        aldh2_gm = 1.0  # parameters already reflect variant
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]
    else:
        aldh2_p = p["ALDH2_star1"]
        aldh2_gm = genotype_modifier(aldh2_gt, "ALDH2")
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]

    tw_aldh2 = get_flux_tissue_weight("ALDH2", tissue, tissue_weight_source)
    v_aldh2 = michaelis_menten(S, aldh2_vmax * aldh2_gm * tw_aldh2, aldh2_km)
    detox_enzymes["ALDH2"] = {
        "flux": round(v_aldh2, 6),
        "genotype": aldh2_gt,
        "genotype_modifier": aldh2_gm,
        "tissue_weight": tw_aldh2,
        "CLint": round((aldh2_vmax * aldh2_gm * tw_aldh2) / aldh2_km, 4),
        "confidence": "high",
    }

    # ALDH1A1 (backup)
    aldh1a1_p = p["ALDH1A1"]
    aldh1a1_gm = genotype_modifier(genotypes.get("ALDH1A1", "NM"), "ALDH1A1")
    tw_aldh1a1 = get_flux_tissue_weight("ALDH1A1", tissue, tissue_weight_source)
    v_aldh1a1 = michaelis_menten(
        S, aldh1a1_p["Vmax_U_per_mg"] * aldh1a1_gm * tw_aldh1a1, aldh1a1_p["Km_uM"]
    )
    detox_enzymes["ALDH1A1"] = {
        "flux": round(v_aldh1a1, 6),
        "genotype_modifier": aldh1a1_gm,
        "tissue_weight": tw_aldh1a1,
        "confidence": "high",
    }

    total_detox = v_aldh2 + v_aldh1a1

    # Ethanol -> Acetaldehyde production
    adh_p = params["carcinogen_classes"]["Aldehyde"]["pathways"]["ethanol_oxidation"]
    adh_gt = genotypes.get("ADH1B", "*1/*1")
    if adh_gt in ("*2/*2", "fast", "RM"):
        adh_params = adh_p["ADH1B_star2"]
    else:
        adh_params = adh_p["ADH1B_star1"]

    eth_conc = params["metadata"]["exposure_defaults_uM"]["ethanol"]
    v_adh = michaelis_menten(
        eth_conc, adh_params["Vmax_U_per_mg"], adh_params["Km_uM"]
    )

    return {
        "activation_enzymes": {
            "ADH1B": {
                "reaction": "Ethanol -> Acetaldehyde",
                "flux": round(v_adh, 6),
                "genotype": adh_gt,
                "genotype_modifier": 1.0,
                "tissue_weight": 1.0,
                "confidence": "high",
            }
        },
        "detox_enzymes": detox_enzymes,
        "total_activation": round(v_adh, 6),
        "total_detox": round(total_detox, 6),
        "unit_note": "Flux in U/mg; ratio dimensionally consistent",
    }


def _compute_nitrosamine_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute NNK activation flux."""
    params = _load_kinetic_params()
    p = params["carcinogen_classes"]["Nitrosamine"]["pathways"]
    activation_enzymes: dict[str, Any] = {}

    # CYP2A13
    cyp2a13_p = p["activation"]["CYP2A13"]
    gm2a13 = genotype_modifier(genotypes.get("CYP2A13", "NM"), "CYP2A13")
    tw2a13 = get_flux_tissue_weight("CYP2A13", tissue, tissue_weight_source)
    v_cyp2a13 = michaelis_menten(
        S,
        cyp2a13_p["Vmax_pmol_min_pmolP450"] * gm2a13 * tw2a13,
        cyp2a13_p["Km_uM"],
    )
    activation_enzymes["CYP2A13"] = {
        "flux": round(v_cyp2a13, 6),
        "genotype_modifier": gm2a13,
        "tissue_weight": tw2a13,
        "confidence": cyp2a13_p["confidence"],
    }

    # CYP2A6 (~200-fold lower CLint)
    gm2a6 = genotype_modifier(genotypes.get("CYP2A6", "NM"), "CYP2A6")
    tw2a6 = get_flux_tissue_weight("CYP2A6", tissue, tissue_weight_source)
    divisor = gm2a13 * tw2a13 + 1e-9
    v_cyp2a6 = v_cyp2a13 * 0.005 * gm2a6 * tw2a6 / divisor
    activation_enzymes["CYP2A6"] = {
        "flux": round(v_cyp2a6, 6),
        "genotype_modifier": gm2a6,
        "tissue_weight": tw2a6,
        "confidence": "moderate",
    }

    total_activation = v_cyp2a13 + v_cyp2a6

    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": {
            "NNAL_reduction": {
                "flux": total_activation * 0.3,
                "genotype_modifier": 1.0,
                "tissue_weight": 1.0,
                "confidence": "estimated",
            }
        },
        "total_activation": round(total_activation, 6),
        "total_detox": round(total_activation * 0.3, 6),
        "unit_note": "Flux in pmol/min/pmolP450",
    }


def _compute_ndma_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute NDMA activation flux via CYP2E1."""
    params = _load_kinetic_params()
    p = params["carcinogen_classes"]["NDMA"]["pathways"]["activation"]["CYP2E1"]
    gm2e1 = genotype_modifier(genotypes.get("CYP2E1", "NM"), "CYP2E1")
    tw2e1 = get_flux_tissue_weight("CYP2E1", tissue, tissue_weight_source)
    v_cyp2e1 = michaelis_menten(
        S, p["Vmax_pmol_min_mg"] * gm2e1 * tw2e1, p["Km_uM"]
    )

    return {
        "activation_enzymes": {
            "CYP2E1": {
                "flux": round(v_cyp2e1, 4),
                "genotype_modifier": gm2e1,
                "tissue_weight": tw2e1,
                "confidence": p["confidence"],
            }
        },
        "detox_enzymes": {},
        "total_activation": round(v_cyp2e1, 4),
        "total_detox": v_cyp2e1 * 0.1,
        "unit_note": "Flux in pmol/min/mg microsomal protein",
    }


def _compute_hca_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute PhIP (HCA) activation and detoxification fluxes."""
    params = _load_kinetic_params()
    p = params["carcinogen_classes"]["HCA"]["pathways"]
    activation_enzymes: dict[str, Any] = {}

    # CYP1A2
    cyp1a2_p = p["activation"]["CYP1A2"]
    gm1a2 = genotype_modifier(genotypes.get("CYP1A2", "NM"), "CYP1A2")
    tw1a2 = get_flux_tissue_weight("CYP1A2", tissue, tissue_weight_source)
    v_1a2 = michaelis_menten(
        S,
        cyp1a2_p["Vmax_nmol_min_nmolP450"] * gm1a2 * tw1a2,
        cyp1a2_p["Km_uM"],
    )
    activation_enzymes["CYP1A2"] = {
        "flux": round(v_1a2, 6),
        "genotype_modifier": gm1a2,
        "tissue_weight": tw1a2,
        "confidence": cyp1a2_p["confidence"],
    }

    # CYP1A1
    cyp1a1_p = p["activation"]["CYP1A1"]
    gm1a1 = genotype_modifier(genotypes.get("CYP1A1", "NM"), "CYP1A1")
    tw1a1 = get_flux_tissue_weight("CYP1A1", tissue, tissue_weight_source)
    v_1a1 = michaelis_menten(
        S,
        cyp1a1_p["Vmax_nmol_min_nmolP450"] * gm1a1 * tw1a1,
        cyp1a1_p["Km_uM"],
    )
    activation_enzymes["CYP1A1"] = {
        "flux": round(v_1a1, 6),
        "genotype_modifier": gm1a1,
        "tissue_weight": tw1a1,
        "confidence": cyp1a1_p["confidence"],
    }

    # CYP1B1
    cyp1b1_p = p["activation"]["CYP1B1"]
    gm1b1 = genotype_modifier(genotypes.get("CYP1B1", "NM"), "CYP1B1")
    tw1b1 = get_flux_tissue_weight("CYP1B1", tissue, tissue_weight_source)
    v_1b1 = michaelis_menten(
        S,
        cyp1b1_p["Vmax_nmol_min_nmolP450"] * gm1b1 * tw1b1,
        cyp1b1_p["Km_uM"],
    )
    activation_enzymes["CYP1B1"] = {
        "flux": round(v_1b1, 6),
        "genotype_modifier": gm1b1,
        "tissue_weight": tw1b1,
        "confidence": cyp1b1_p["confidence"],
    }

    total_activation = v_1a2 + v_1a1 + v_1b1

    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": {
            "NAT2_acetylation": {
                "flux": total_activation * 0.2,
                "genotype_modifier": 1.0,
                "tissue_weight": 1.0,
                "confidence": "estimated",
            }
        },
        "total_activation": round(total_activation, 6),
        "total_detox": round(total_activation * 0.2, 6),
        "unit_note": "Flux in nmol/min/nmolP450",
    }


def _compute_benzene_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute benzene activation flux via CYP2E1."""
    params = _load_kinetic_params()
    p = params["carcinogen_classes"]["Benzene"]["pathways"]["activation"]["CYP2E1_liver"]
    gm2e1 = genotype_modifier(genotypes.get("CYP2E1", "NM"), "CYP2E1")
    tw2e1 = get_flux_tissue_weight("CYP2E1", tissue, tissue_weight_source)
    km = p["Km_uM"]
    vmax_estimate = 100.0
    v_cyp2e1 = michaelis_menten(S, vmax_estimate * gm2e1 * tw2e1, km)

    # NQO1 detox (benzoquinone reduction)
    nqo1_gm = genotype_modifier(genotypes.get("NQO1", "NM"), "NQO1")
    v_nqo1 = v_cyp2e1 * 0.4 * nqo1_gm

    return {
        "activation_enzymes": {
            "CYP2E1": {
                "flux": round(v_cyp2e1, 4),
                "genotype_modifier": gm2e1,
                "tissue_weight": tw2e1,
                "confidence": "moderate",
            }
        },
        "detox_enzymes": {
            "NQO1": {
                "flux": round(v_nqo1, 4),
                "genotype_modifier": nqo1_gm,
                "tissue_weight": 1.0,
                "confidence": "low",
            }
        },
        "total_activation": round(v_cyp2e1, 4),
        "total_detox": round(v_nqo1, 4),
        "unit_note": "Flux in pmol/min/mg (estimated Vmax)",
    }


def _compute_aromatic_amine_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute a semi-quantitative 4-ABP / aromatic-amine activation proxy."""
    cfg = _get_proxy_class_params("AromaticAmines")
    act1_p = cfg["activation_terms"]["CYP1A2"]
    v_1a2, gm1a2, tw1a2 = _compute_proxy_mm_term(
        act1_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    act2_p = cfg["activation_terms"]["NAT1"]
    v_nat1, gm_nat1, tw_nat1 = _compute_proxy_mm_term(
        act2_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    total_activation = v_1a2 + v_nat1

    detox_nat2 = cfg["detox_terms"]["NAT2"]
    v_nat2, gm_nat2, tw_nat2 = _compute_proxy_mm_term(
        detox_nat2,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_gstm1 = cfg["detox_terms"]["GSTM1"]
    v_gstm1, gm_gstm1, tw_gstm1 = _compute_proxy_mm_term(
        detox_gstm1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_gstp1 = cfg["detox_terms"]["GSTP1"]
    v_gstp1, gm_gstp1, tw_gstp1 = _compute_proxy_mm_term(
        detox_gstp1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    repair_xpc_p = cfg["repair_terms"]["XPC"]
    repair_xpc, gm_xpc, tw_xpc = _compute_proxy_repair_term(
        total_activation,
        repair_xpc_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    repair_ercc2_p = cfg["repair_terms"]["ERCC2"]
    repair_ercc2, gm_ercc2, tw_ercc2 = _compute_proxy_repair_term(
        total_activation,
        repair_ercc2_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    total_detox = v_nat2 + v_gstm1 + v_gstp1 + repair_xpc + repair_ercc2

    return {
        "activation_enzymes": {
            "CYP1A2": {
                "flux": round(v_1a2, 6),
                "genotype_modifier": gm1a2,
                "tissue_weight": tw1a2,
                "confidence": act1_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": act1_p["note"],
            },
            "NAT1": {
                "flux": round(v_nat1, 6),
                "genotype_modifier": gm_nat1,
                "tissue_weight": tw_nat1,
                "confidence": act2_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": act2_p["note"],
            },
        },
        "detox_enzymes": {
            "NAT2": {
                "flux": round(v_nat2, 6),
                "genotype_modifier": gm_nat2,
                "tissue_weight": tw_nat2,
                "confidence": detox_nat2["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_nat2["note"],
            },
            "GSTM1": {
                "flux": round(v_gstm1, 6),
                "genotype_modifier": gm_gstm1,
                "tissue_weight": tw_gstm1,
                "confidence": detox_gstm1["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_gstm1["note"],
            },
            "GSTP1": {
                "flux": round(v_gstp1, 6),
                "genotype_modifier": gm_gstp1,
                "tissue_weight": tw_gstp1,
                "confidence": detox_gstp1["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_gstp1["note"],
            },
            "XPC": {
                "flux": round(repair_xpc, 6),
                "genotype_modifier": gm_xpc,
                "tissue_weight": tw_xpc,
                "confidence": repair_xpc_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_xpc_p["note"],
            },
            "ERCC2": {
                "flux": round(repair_ercc2, 6),
                "genotype_modifier": gm_ercc2,
                "tissue_weight": tw_ercc2,
                "confidence": repair_ercc2_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_ercc2_p["note"],
            },
        },
        "total_activation": round(total_activation, 6),
        "total_detox": round(total_detox, 6),
        "unit_note": cfg["unit_note"],
    }


def _compute_estrogen_metabolite_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute catechol-estrogen / estradiol-quinone burden with detox proxies."""
    cfg = _get_proxy_class_params("EstrogenMetabolites")

    act1_p = cfg["activation_terms"]["CYP1B1"]
    v_1b1, gm1b1, tw1b1 = _compute_proxy_mm_term(
        act1_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    act2_p = cfg["activation_terms"]["CYP1A1"]
    v_1a1, gm1a1, tw1a1 = _compute_proxy_mm_term(
        act2_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    act3_p = cfg["activation_terms"]["CYP1A2"]
    v_1a2, gm1a2, tw1a2 = _compute_proxy_mm_term(
        act3_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    total_activation = v_1b1 + v_1a1 + v_1a2

    detox_comt = cfg["detox_terms"]["COMT"]
    v_comt, gm_comt, tw_comt = _compute_proxy_mm_term(
        detox_comt,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_sult = cfg["detox_terms"]["SULT1E1"]
    v_sult1e1, gm_sult1e1, tw_sult1e1 = _compute_proxy_mm_term(
        detox_sult,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_ugt = cfg["detox_terms"]["UGT2B7"]
    v_ugt2b7, gm_ugt2b7, tw_ugt2b7 = _compute_proxy_mm_term(
        detox_ugt,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_gstp1 = cfg["detox_terms"]["GSTP1"]
    v_gstp1, gm_gstp1, tw_gstp1 = _compute_proxy_mm_term(
        detox_gstp1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    repair_xrcc1_p = cfg["repair_terms"]["XRCC1"]
    repair_xrcc1, gm_xrcc1, tw_xrcc1 = _compute_proxy_repair_term(
        total_activation,
        repair_xrcc1_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    total_detox = v_comt + v_sult1e1 + v_ugt2b7 + v_gstp1 + repair_xrcc1

    return {
        "activation_enzymes": {
            "CYP1B1": {
                "flux": round(v_1b1, 6),
                "genotype_modifier": gm1b1,
                "tissue_weight": tw1b1,
                "confidence": act1_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": act1_p["note"],
            },
            "CYP1A1": {
                "flux": round(v_1a1, 6),
                "genotype_modifier": gm1a1,
                "tissue_weight": tw1a1,
                "confidence": act2_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": act2_p["note"],
            },
            "CYP1A2": {
                "flux": round(v_1a2, 6),
                "genotype_modifier": gm1a2,
                "tissue_weight": tw1a2,
                "confidence": act3_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": act3_p["note"],
            },
        },
        "detox_enzymes": {
            "COMT": {
                "flux": round(v_comt, 6),
                "genotype_modifier": gm_comt,
                "tissue_weight": tw_comt,
                "confidence": detox_comt["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_comt["note"],
            },
            "SULT1E1": {
                "flux": round(v_sult1e1, 6),
                "genotype_modifier": gm_sult1e1,
                "tissue_weight": tw_sult1e1,
                "confidence": detox_sult["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_sult["note"],
            },
            "UGT2B7": {
                "flux": round(v_ugt2b7, 6),
                "genotype_modifier": gm_ugt2b7,
                "tissue_weight": tw_ugt2b7,
                "confidence": detox_ugt["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_ugt["note"],
            },
            "GSTP1": {
                "flux": round(v_gstp1, 6),
                "genotype_modifier": gm_gstp1,
                "tissue_weight": tw_gstp1,
                "confidence": detox_gstp1["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_gstp1["note"],
            },
            "XRCC1": {
                "flux": round(repair_xrcc1, 6),
                "genotype_modifier": gm_xrcc1,
                "tissue_weight": tw_xrcc1,
                "confidence": repair_xrcc1_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_xrcc1_p["note"],
            },
        },
        "total_activation": round(total_activation, 6),
        "total_detox": round(total_detox, 6),
        "unit_note": cfg["unit_note"],
    }


def _compute_ndea_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute NDEA-specific dietary nitroso activation with pulmonary/liver proxies."""
    cfg = _get_proxy_class_params("NDEA")

    act_2a13 = cfg["activation_terms"]["CYP2A13"]
    v_2a13, gm2a13, tw2a13 = _compute_proxy_mm_term(
        act_2a13,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    act_2e1 = cfg["activation_terms"]["CYP2E1"]
    v_2e1, gm2e1, tw2e1 = _compute_proxy_mm_term(
        act_2e1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    act_2a6 = cfg["activation_terms"]["CYP2A6"]
    v_2a6, gm2a6, tw2a6 = _compute_proxy_mm_term(
        act_2a6,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    total_activation = v_2a13 + v_2e1 + v_2a6

    repair_mgmt_p = cfg["repair_terms"]["MGMT"]
    repair_mgmt, gm_mgmt, tw_mgmt = _compute_proxy_repair_term(
        total_activation,
        repair_mgmt_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    detox_gstp1 = cfg["detox_terms"]["GSTP1"]
    v_gstp1, gm_gstp1, tw_gstp1 = _compute_proxy_mm_term(
        detox_gstp1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )
    total_detox = repair_mgmt + v_gstp1

    return {
        "activation_enzymes": {
            "CYP2A13": {
                "flux": round(v_2a13, 6),
                "genotype_modifier": gm2a13,
                "tissue_weight": tw2a13,
                "confidence": act_2a13["confidence"],
                "note": act_2a13["note"],
            },
            "CYP2E1": {
                "flux": round(v_2e1, 6),
                "genotype_modifier": gm2e1,
                "tissue_weight": tw2e1,
                "confidence": act_2e1["confidence"],
                "note": act_2e1["note"],
            },
            "CYP2A6": {
                "flux": round(v_2a6, 6),
                "genotype_modifier": gm2a6,
                "tissue_weight": tw2a6,
                "confidence": act_2a6["confidence"],
                "note": act_2a6["note"],
            },
        },
        "detox_enzymes": {
            "MGMT": {
                "flux": round(repair_mgmt, 6),
                "genotype_modifier": gm_mgmt,
                "tissue_weight": tw_mgmt,
                "confidence": repair_mgmt_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_mgmt_p["note"],
            },
            "GSTP1": {
                "flux": round(v_gstp1, 6),
                "genotype_modifier": gm_gstp1,
                "tissue_weight": tw_gstp1,
                "confidence": detox_gstp1["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_gstp1["note"],
            },
        },
        "total_activation": round(total_activation, 6),
        "total_detox": round(total_detox, 6),
        "unit_note": cfg["unit_note"],
    }


def _compute_vinyl_chloride_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute vinyl-chloride activation with GST/repair attenuation proxies."""
    cfg = _get_proxy_class_params("VinylChloride")

    act_2e1 = cfg["activation_terms"]["CYP2E1"]
    v_2e1, gm2e1, tw2e1 = _compute_proxy_mm_term(
        act_2e1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_gstt1 = cfg["detox_terms"]["GSTT1"]
    v_gstt1, gm_gstt1, tw_gstt1 = _compute_proxy_mm_term(
        detox_gstt1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_gstm1 = cfg["detox_terms"]["GSTM1"]
    v_gstm1, gm_gstm1, tw_gstm1 = _compute_proxy_mm_term(
        detox_gstm1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_ephx1 = cfg["detox_terms"]["EPHX1"]
    v_ephx1, gm_ephx1, tw_ephx1 = _compute_proxy_mm_term(
        detox_ephx1,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    detox_aldh2 = cfg["detox_terms"]["ALDH2"]
    v_aldh2, gm_aldh2, tw_aldh2 = _compute_proxy_mm_term(
        detox_aldh2,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    repair_ercc2_p = cfg["repair_terms"]["ERCC2"]
    repair_ercc2, gm_ercc2, tw_ercc2 = _compute_proxy_repair_term(
        v_2e1,
        repair_ercc2_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    total_detox = v_gstt1 + v_gstm1 + v_ephx1 + v_aldh2 + repair_ercc2

    return {
        "activation_enzymes": {
            "CYP2E1": {
                "flux": round(v_2e1, 6),
                "genotype_modifier": gm2e1,
                "tissue_weight": tw2e1,
                "confidence": act_2e1["confidence"],
                "note": act_2e1["note"],
            }
        },
        "detox_enzymes": {
            "GSTT1": {
                "flux": round(v_gstt1, 6),
                "genotype_modifier": gm_gstt1,
                "tissue_weight": tw_gstt1,
                "confidence": detox_gstt1["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_gstt1["note"],
            },
            "GSTM1": {
                "flux": round(v_gstm1, 6),
                "genotype_modifier": gm_gstm1,
                "tissue_weight": tw_gstm1,
                "confidence": detox_gstm1["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_gstm1["note"],
            },
            "EPHX1": {
                "flux": round(v_ephx1, 6),
                "genotype_modifier": gm_ephx1,
                "tissue_weight": tw_ephx1,
                "confidence": detox_ephx1["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_ephx1["note"],
            },
            "ALDH2": {
                "flux": round(v_aldh2, 6),
                "genotype_modifier": gm_aldh2,
                "tissue_weight": tw_aldh2,
                "confidence": detox_aldh2["confidence"],
                "kinetics": "semi_quantitative",
                "note": detox_aldh2["note"],
            },
            "ERCC2": {
                "flux": round(repair_ercc2, 6),
                "genotype_modifier": gm_ercc2,
                "tissue_weight": tw_ercc2,
                "confidence": repair_ercc2_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_ercc2_p["note"],
            },
        },
        "total_activation": round(v_2e1, 6),
        "total_detox": round(total_detox, 6),
        "unit_note": cfg["unit_note"],
    }


def _compute_uv_radiation_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute UV-driven DNA-damage burden as photoproducts versus repair capacity."""
    cfg = _get_proxy_class_params("UV_Radiation")
    uvb_p = cfg["activation_terms"]["UVB_photoproduct_burden"]
    v_uvb, _, tissue_signal = _compute_proxy_hill_term(
        uvb_p,
        {},
        tissue,
        S,
        tissue_weight_source,
    )
    uva_p = cfg["activation_terms"]["UVA_oxidative_tail"]
    v_uva, _, _ = _compute_proxy_hill_term(
        uva_p,
        {},
        tissue,
        S,
        tissue_weight_source,
    )
    total_activation = v_uvb + v_uva

    repair_xpc_p = cfg["repair_terms"]["XPC"]
    repair_xpc, gm_xpc, tw_xpc = _compute_proxy_repair_term(
        total_activation,
        repair_xpc_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    repair_ercc2_p = cfg["repair_terms"]["ERCC2"]
    repair_ercc2, gm_ercc2, tw_ercc2 = _compute_proxy_repair_term(
        total_activation,
        repair_ercc2_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    repair_ogg1_p = cfg["repair_terms"]["OGG1"]
    repair_ogg1, gm_ogg1, tw_ogg1 = _compute_proxy_repair_term(
        total_activation,
        repair_ogg1_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    repair_polh_p = cfg["repair_terms"]["POLH"]
    repair_polh, gm_polh, tw_polh = _compute_proxy_repair_term(
        total_activation,
        repair_polh_p,
        genotypes,
        tissue,
        tissue_weight_source,
    )
    total_detox = repair_xpc + repair_ercc2 + repair_ogg1 + repair_polh

    return {
        "activation_enzymes": {
            "UVB_photoproduct_burden": {
                "flux": round(v_uvb, 6),
                "genotype_modifier": 1.0,
                "tissue_weight": tissue_signal,
                "confidence": uvb_p["confidence"],
                "kinetics": "damage_proxy",
                "note": uvb_p["note"],
            },
            "UVA_oxidative_tail": {
                "flux": round(v_uva, 6),
                "genotype_modifier": 1.0,
                "tissue_weight": tissue_signal,
                "confidence": uva_p["confidence"],
                "kinetics": "damage_proxy",
                "note": uva_p["note"],
            },
        },
        "detox_enzymes": {
            "XPC": {
                "flux": round(repair_xpc, 6),
                "genotype_modifier": gm_xpc,
                "tissue_weight": tw_xpc,
                "confidence": repair_xpc_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_xpc_p["note"],
            },
            "ERCC2": {
                "flux": round(repair_ercc2, 6),
                "genotype_modifier": gm_ercc2,
                "tissue_weight": tw_ercc2,
                "confidence": repair_ercc2_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_ercc2_p["note"],
            },
            "OGG1": {
                "flux": round(repair_ogg1, 6),
                "genotype_modifier": gm_ogg1,
                "tissue_weight": tw_ogg1,
                "confidence": repair_ogg1_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_ogg1_p["note"],
            },
            "POLH": {
                "flux": round(repair_polh, 6),
                "genotype_modifier": gm_polh,
                "tissue_weight": tw_polh,
                "confidence": repair_polh_p["confidence"],
                "kinetics": "repair_proxy",
                "note": repair_polh_p["note"],
            },
        },
        "total_activation": round(total_activation, 6),
        "total_detox": round(total_detox, 6),
        "unit_note": cfg["unit_note"],
    }


def _compute_chlorinated_solvent_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute TCE-centered chlorinated-solvent bioactivation with proxy clearance."""
    cfg = _get_proxy_class_params("ChlorinatedSolvent")
    activation_enzymes: dict[str, Any] = {}

    oxidation_p = cfg["activation_terms"]["CYP2E1"]
    v_oxidation, gm2e1, tw2e1 = _compute_proxy_mm_term(
        oxidation_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )
    activation_enzymes["CYP2E1"] = {
        "flux": round(v_oxidation, 6),
        "genotype_modifier": gm2e1,
        "tissue_weight": tw2e1,
        "confidence": oxidation_p["confidence"],
        "note": oxidation_p["note"],
    }

    gsh_p = cfg["activation_terms"]["GSTT1"]
    v_gsh, gm_gstt1, tw_gstt1 = _compute_proxy_mm_term(
        gsh_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )
    activation_enzymes["GSTT1"] = {
        "flux": round(v_gsh, 6),
        "genotype_modifier": gm_gstt1,
        "tissue_weight": tw_gstt1,
        "confidence": gsh_p["confidence"],
        "note": gsh_p["note"],
    }

    detox_p = cfg["detox_terms"]["non_genotoxic_clearance_proxy"]
    v_clearance = v_oxidation * float(detox_p["scale"])
    detox_enzymes = {
        "non_genotoxic_clearance_proxy": {
            "flux": round(v_clearance, 6),
            "genotype_modifier": 1.0,
            "tissue_weight": max(tw2e1, 0.2),
            "confidence": detox_p["confidence"],
            "note": detox_p["note"],
        }
    }

    total_activation = v_oxidation + v_gsh

    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": detox_enzymes,
        "total_activation": round(total_activation, 6),
        "total_detox": round(v_clearance, 6),
        "unit_note": cfg["unit_note"],
    }


def _compute_dioxin_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute receptor-mediated dioxin signaling as an induction burden score."""
    cfg = _get_proxy_class_params("Dioxin")
    signal_p = cfg["signal"]
    signal_strength, _, tissue_signal = _compute_proxy_hill_term(
        signal_p,
        {},
        tissue,
        S,
        tissue_weight_source,
    )

    act1_p = cfg["activation_terms"]["CYP1A1"]
    gm1a1 = _proxy_genotype_modifier(genotypes.get("CYP1A1", "NM"), "CYP1A1")
    tw1a1 = get_flux_tissue_weight("CYP1A1", tissue, tissue_weight_source)
    v1a1 = signal_strength * float(act1_p["induction_factor"]) * gm1a1 * tw1a1

    act2_p = cfg["activation_terms"]["CYP1B1"]
    gm1b1 = _proxy_genotype_modifier(genotypes.get("CYP1B1", "NM"), "CYP1B1")
    tw1b1 = get_flux_tissue_weight("CYP1B1", tissue, tissue_weight_source)
    v1b1 = signal_strength * float(act2_p["induction_factor"]) * gm1b1 * tw1b1

    feedback_p = cfg["detox_terms"]["AHRR_feedback"]
    v_feedback = signal_strength * float(feedback_p["scale_of_signal"])

    return {
        "activation_enzymes": {
            "CYP1A1": {
                "flux": round(v1a1, 6),
                "genotype_modifier": gm1a1,
                "tissue_weight": tw1a1,
                "confidence": act1_p["confidence"],
                "kinetics": "receptor_mediated",
                "note": act1_p["note"],
            },
            "CYP1B1": {
                "flux": round(v1b1, 6),
                "genotype_modifier": gm1b1,
                "tissue_weight": tw1b1,
                "confidence": act2_p["confidence"],
                "kinetics": "receptor_mediated",
                "note": act2_p["note"],
            },
        },
        "detox_enzymes": {
            "AHRR_feedback": {
                "flux": round(v_feedback, 6),
                "genotype_modifier": 1.0,
                "tissue_weight": tissue_signal,
                "confidence": feedback_p["confidence"],
                "kinetics": "receptor_mediated",
                "note": feedback_p["note"],
            }
        },
        "total_activation": round(v1a1 + v1b1, 6),
        "total_detox": round(v_feedback, 6),
        "unit_note": cfg["unit_note"],
    }


def _compute_heavy_metal_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> FluxResultDict:
    """Compute semi-quantitative heavy-metal risk as ROS burden vs methylation clearance."""
    cfg = _get_proxy_class_params("HeavyMetal")

    as3mt_p = cfg["detox_terms"]["AS3MT"]
    v_as3mt, gm_as3mt, tw_as3mt = _compute_proxy_mm_term(
        as3mt_p,
        genotypes,
        tissue,
        S,
        tissue_weight_source,
    )

    ros_p = cfg["activation_terms"]["general_ROS"]
    v_ros, _, tw_ros = _compute_proxy_saturating_term(
        ros_p,
        {},
        tissue,
        S,
        tissue_weight_source,
    )

    cadmium_p = cfg["activation_terms"]["cadmium_stress_proxy"]
    v_cadmium, _, tw_cadmium = _compute_proxy_saturating_term(
        cadmium_p,
        {},
        tissue,
        S,
        tissue_weight_source,
    )

    total_activation = v_ros + v_cadmium

    return {
        "activation_enzymes": {
            "general_ROS": {
                "flux": round(v_ros, 6),
                "genotype_modifier": 1.0,
                "tissue_weight": tw_ros,
                "confidence": ros_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": ros_p["note"],
            },
            "cadmium_stress_proxy": {
                "flux": round(v_cadmium, 6),
                "genotype_modifier": 1.0,
                "tissue_weight": tw_cadmium,
                "confidence": cadmium_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": cadmium_p["note"],
            },
        },
        "detox_enzymes": {
            "AS3MT": {
                "flux": round(v_as3mt, 6),
                "genotype_modifier": gm_as3mt,
                "tissue_weight": tw_as3mt,
                "confidence": as3mt_p["confidence"],
                "kinetics": "semi_quantitative",
                "note": as3mt_p["note"],
            }
        },
        "total_activation": round(total_activation, 6),
        "total_detox": round(v_as3mt, 6),
        "unit_note": cfg["unit_note"],
    }


# ── Helpers for dataclass conversion ───────────────────────────────────────


def _enzyme_flux_from_dict(name: str, d: JsonDict) -> EnzymeFlux:
    """Convert an internal enzyme dict to an EnzymeFlux dataclass."""
    return EnzymeFlux(
        enzyme=name,
        flux=float(d.get("flux", 0.0)),
        genotype_modifier=float(d.get("genotype_modifier", 1.0)),
        tissue_weight=float(d.get("tissue_weight", 1.0)),
        confidence=str(d.get("confidence", "unknown")),
        induction_modifier=float(d.get("induction_modifier", 1.0)),
        qivive_scale=float(d.get("qivive_scale", 1.0)),
        fraction=float(d.get("fraction", 0.0)),
        kinetics=str(d.get("kinetics", "michaelis_menten")),
        note=str(d.get("note", "")),
        model_kind=str(d.get("model_kind", "measured_kinetics")),
        parameter_source=str(d.get("parameter_source", _KINETIC_PARAMS_FILE.name)),
        provenance_ref=str(d.get("provenance_ref", "")),
        provenance_sources=[str(source) for source in d.get("provenance_sources", [])],
        parameter_basis=str(d.get("parameter_basis", "")),
    )


# ── Public API ─────────────────────────────────────────────────────────────

_DISPATCH: dict[str, FluxCalculator] = {
    "PAH": _compute_pah_flux,
    "Aflatoxin": _compute_aflatoxin_flux,
    "Aldehyde": _compute_aldehyde_flux,
    "Nitrosamine": _compute_nitrosamine_flux,
    "NDMA": _compute_ndma_flux,
    "NDEA": _compute_ndea_flux,
    "HCA": _compute_hca_flux,
    "AromaticAmines": _compute_aromatic_amine_flux,
    "EstrogenMetabolites": _compute_estrogen_metabolite_flux,
    "Benzene": _compute_benzene_flux,
    "VinylChloride": _compute_vinyl_chloride_flux,
    "ChlorinatedSolvent": _compute_chlorinated_solvent_flux,
    "UV_Radiation": _compute_uv_radiation_flux,
    "Dioxin": _compute_dioxin_flux,
    "HeavyMetal": _compute_heavy_metal_flux,
}


def compute_pathway_flux(
    carcinogen_class: CarcinogenClass | str,
    genotypes: dict[str, str],
    tissue: str,
    substrate_conc_uM: float | None = None,
    tissue_weight_source: FluxTissueWeightSource | str = FluxTissueWeightSource.CURATED,
    *,
    lifestyle: LifestyleMap | None = None,
    induction_factors: Mapping[str, float] | None = None,
    qivive: bool = False,
    qivive_context: Mapping[str, float] | None = None,
    steady_state_context: Mapping[str, float] | None = None,
) -> PathwayFluxResult:
    """Compute activation/detoxification flux for a carcinogen class.

    Args:
        carcinogen_class: Carcinogen class enum or string.
        genotypes: Gene-to-phenotype mapping (e.g. ``{"CYP1A1": "NM"}``).
        tissue: Tissue name (e.g. "Lung", "Liver").
        substrate_conc_uM: Substrate concentration; defaults to
            environmental exposure default.
        tissue_weight_source: ``"curated"`` (source-parity default) or
            ``"gtex"`` for quantitative GTEx weighting.
        lifestyle: Optional co-exposure state. When provided, lifestyle-driven
            induction folds from :mod:`ExposoGraph.interaction_engine` are
            applied as Vmax multipliers.
        induction_factors: Optional explicit per-gene Vmax folds. These
            override lifestyle-derived values for matching genes.
        qivive: If ``True``, apply a common tissue-level QIVIVE scale based on
            MPPGL and organ weight to reported flux magnitudes.
        qivive_context: Optional ``{"mppgl_mg_per_g": x, "organ_weight_g": y}``
            override for QIVIVE scaling.
        steady_state_context: Optional PBPK steady-state context override
            (body weight, central volume, tissue partition, blood-flow
            fraction, background clearance, and related first-order rates).

    Returns:
        :class:`PathwayFluxResult` with activation, detoxification,
        net ratio and risk classification. Returned results also expose
        ``model_kind`` / ``parameter_source`` at the class level, while
        proxy-backed enzyme terms include resolved provenance metadata.
    """
    cls_str = carcinogen_class.value if isinstance(carcinogen_class, CarcinogenClass) else carcinogen_class
    weight_source = _normalize_tissue_weight_source(tissue_weight_source)

    if cls_str not in _DISPATCH:
        return PathwayFluxResult(
            carcinogen_class=cls_str,
            tissue=tissue,
            substrate_concentration_uM=0.0,
            genotypes_used=genotypes,
            activation_enzymes=[],
            detox_enzymes=[],
            total_activation=0.0,
            total_detox=0.0,
            net_ratio=0.0,
            susceptibility_score_log2=0.0,
            risk_classification=RiskClassification.INSUFFICIENT_DATA,
            tissue_weight_source=weight_source,
            model_kind="unavailable",
            parameter_source="",
            warnings=[f"No quantitative model for '{cls_str}'"],
        )

    if substrate_conc_uM is None:
        substrate_conc_uM = _get_default_concentration(cls_str)

    resolved_induction = _resolve_induction_factors(lifestyle, induction_factors)
    result = _DISPATCH[cls_str](genotypes, tissue, substrate_conc_uM, weight_source)
    result = _annotate_flux_result_metadata(cls_str, result)
    result = _apply_induction_modifiers(result, resolved_induction)

    qivive_used_context: dict[str, float] = {}
    if qivive:
        qivive_used_context = _qivive_context_for_tissue(tissue, qivive_context)
        result = _apply_qivive_scale(result, qivive_used_context)

    act = float(result["total_activation"])
    det = float(result["total_detox"])

    if det > 0:
        net_ratio = act / det
    elif act > 0:
        net_ratio = 999.0
    else:
        net_ratio = 1.0

    # Fractional contributions
    for edata in result.get("activation_enzymes", {}).values():
        if isinstance(edata, dict) and "flux" in edata and act > 0:
            edata["fraction"] = round(edata["flux"] / act, 4)
    for edata in result.get("detox_enzymes", {}).values():
        if isinstance(edata, dict) and "flux" in edata and det > 0:
            edata["fraction"] = round(edata["flux"] / det, 4)

    risk = classify_risk(net_ratio)
    susceptibility_score = _susceptibility_score_log2(net_ratio)
    steady_state_input = dict(steady_state_context or {})
    if qivive_used_context and "organ_weight_g" not in steady_state_input:
        steady_state_input["organ_weight_g"] = qivive_used_context["organ_weight_g"]
    steady_state = solve_flux_steady_state(
        substrate_conc_uM,
        act,
        det,
        tissue,
        context=steady_state_input or None,
    )
    steady_state_proxy = {
        "reactive_intermediate_proxy_uM": steady_state.concentrations_uM[
            "reactive_intermediate_uM"
        ],
        "detoxified_metabolite_proxy_uM": steady_state.concentrations_uM[
            "detoxified_metabolite_uM"
        ],
    }

    # Convert to dataclasses
    act_enzymes = [
        _enzyme_flux_from_dict(name, data)
        for name, data in result.get("activation_enzymes", {}).items()
        if isinstance(data, dict)
    ]
    det_enzymes = [
        _enzyme_flux_from_dict(name, data)
        for name, data in result.get("detox_enzymes", {}).items()
        if isinstance(data, dict)
    ]

    # Warnings
    warn_list: list[str] = []
    all_enzymes = {
        **result.get("activation_enzymes", {}),
        **result.get("detox_enzymes", {}),
    }
    if any(
        isinstance(v, dict) and v.get("confidence") in ("estimated", "low")
        for v in all_enzymes.values()
    ):
        warn_list.append("ESTIMATED_PARAMS")
    if resolved_induction:
        warn_list.append("INDUCTION_FACTORS_APPLIED")
    if qivive:
        warn_list.append("QIVIVE_SCALE_APPLIED")

    return PathwayFluxResult(
        carcinogen_class=cls_str,
        tissue=tissue,
        substrate_concentration_uM=substrate_conc_uM,
        genotypes_used=genotypes,
        activation_enzymes=act_enzymes,
        detox_enzymes=det_enzymes,
        total_activation=round(act, 6),
        total_detox=round(det, 6),
        net_ratio=round(net_ratio, 4),
        susceptibility_score_log2=susceptibility_score,
        risk_classification=risk,
        tissue_weight_source=weight_source,
        model_kind=result.get("model_kind", "measured_kinetics"),
        parameter_source=result.get("parameter_source", _KINETIC_PARAMS_FILE.name),
        unit_note=result.get("unit_note", ""),
        warnings=warn_list,
        induction_factors_used=resolved_induction,
        qivive_applied=qivive,
        qivive_context=qivive_used_context,
        steady_state_concentrations_uM=steady_state.concentrations_uM,
        steady_state_model=steady_state.model,
        steady_state_concentration_proxy_uM=steady_state_proxy,
    )


def compute_full_profile(
    genotypes: dict[str, str],
    tissue: str,
    exposure_profile: dict[str, float] | None = None,
    tissue_weight_source: FluxTissueWeightSource | str = FluxTissueWeightSource.CURATED,
    *,
    lifestyle: LifestyleMap | None = None,
    induction_factors: Mapping[str, float] | None = None,
    qivive: bool = False,
    qivive_context: Mapping[str, float] | None = None,
    steady_state_context: Mapping[str, float] | None = None,
) -> FullProfileResult:
    """Compute flux across all supported carcinogen classes.

    Args:
        genotypes: Gene-to-phenotype mapping.
        tissue: Tissue name.
        exposure_profile: Optional overrides ``{class: concentration_uM}``.
        tissue_weight_source: ``"curated"`` (default) or ``"gtex"``.
        lifestyle: Optional co-exposure state used to resolve enzyme induction.
        induction_factors: Optional explicit per-gene Vmax induction folds.
        qivive: Apply optional MPPGL/organ-weight QIVIVE scaling.
        qivive_context: Optional QIVIVE context override.
        steady_state_context: Optional PBPK steady-state context override.

    Returns:
        :class:`FullProfileResult` with per-class results and summary.
    """
    classes = list(_DISPATCH)
    exposure_profile = exposure_profile or {}
    weight_source = _normalize_tissue_weight_source(tissue_weight_source)

    results: dict[str, PathwayFluxResult] = {}
    for cls in classes:
        conc = exposure_profile.get(cls)
        results[cls] = compute_pathway_flux(
            cls,
            genotypes,
            tissue,
            conc,
            tissue_weight_source=weight_source,
            lifestyle=lifestyle,
            induction_factors=induction_factors,
            qivive=qivive,
            qivive_context=qivive_context,
            steady_state_context=steady_state_context,
        )

    elevated = [
        c
        for c, r in results.items()
        if r.risk_classification in (RiskClassification.ELEVATED, RiskClassification.HIGH)
    ]
    moderate = [
        c for c, r in results.items() if r.risk_classification == RiskClassification.MODERATE
    ]

    return FullProfileResult(
        tissue=tissue,
        genotypes=genotypes,
        per_class_results=results,
        elevated_or_high_risk_classes=elevated,
        moderate_risk_classes=moderate,
        total_classes_modeled=len(classes),
        tissue_weight_source=weight_source,
    )


def sensitivity_analysis(
    carcinogen_class: CarcinogenClass | str,
    gene: str,
    tissue: str,
    substrate_conc_uM: float | None = None,
    tissue_weight_source: FluxTissueWeightSource | str = FluxTissueWeightSource.CURATED,
) -> SensitivityResult:
    """Assess how changing one gene's genotype shifts the net ratio.

    Tests PM, IM, NM, RM, UM (or null for GSTM1/GSTT1, or star alleles
    for ALDH2) while holding all other genes at NM.

    Args:
        carcinogen_class: Carcinogen class.
        gene: Gene to vary.
        tissue: Tissue name.
        substrate_conc_uM: Optional substrate concentration.
        tissue_weight_source: ``"curated"`` (default) or ``"gtex"``.

    Returns:
        :class:`SensitivityResult` with per-phenotype ratios.
    """
    cls_str = carcinogen_class.value if isinstance(carcinogen_class, CarcinogenClass) else carcinogen_class
    weight_source = _normalize_tissue_weight_source(tissue_weight_source)

    base_genotypes = {
        g: "NM"
        for g in [
            "CYP1A1", "CYP1B1", "CYP1A2", "CYP3A4", "CYP3A5",
            "CYP2A13", "CYP2A6", "CYP2E1", "EPHX1",
            "GSTM1", "GSTT1", "GSTP1", "ALDH2", "ALDH1A1", "ADH1B",
            "NQO1", "CYP2D6", "CYP2B6", "NAT1", "NAT2",
            "COMT", "SULT1E1", "UGT2B7", "AS3MT",
            "XPC", "ERCC2", "XRCC1", "OGG1", "MGMT", "POLH",
        ]
    }

    test_phenotypes = ["PM", "IM", "NM", "RM", "UM"]
    if gene in ("GSTM1", "GSTT1"):
        test_phenotypes = ["null", "NM", "RM"]
    elif gene == "ALDH2":
        test_phenotypes = ["*1/*1", "*1/*2", "*2/*2"]

    baseline = compute_pathway_flux(
        cls_str,
        base_genotypes,
        tissue,
        substrate_conc_uM,
        tissue_weight_source=weight_source,
    )
    baseline_ratio = baseline.net_ratio

    sensitivity_results: dict[str, dict[str, Any]] = {}
    for pt in test_phenotypes:
        gt = dict(base_genotypes)
        gt[gene] = pt
        r = compute_pathway_flux(
            cls_str,
            gt,
            tissue,
            substrate_conc_uM,
            tissue_weight_source=weight_source,
        )
        net = r.net_ratio

        if isinstance(baseline_ratio, (int, float)) and baseline_ratio > 0:
            delta = round(net - baseline_ratio, 4)
            fold_change = round(net / baseline_ratio, 4)
        else:
            delta = None
            fold_change = None

        sensitivity_results[pt] = {
            "net_ratio": net,
            "risk_classification": r.risk_classification.value,
            "delta_from_NM_baseline": delta,
            "fold_change_from_NM": fold_change,
        }

    max_fc = max(
        (
            v["fold_change_from_NM"]
            for v in sensitivity_results.values()
            if v["fold_change_from_NM"] is not None
        ),
        default=None,
    )

    return SensitivityResult(
        carcinogen_class=cls_str,
        gene_varied=gene,
        tissue=tissue,
        baseline_ratio=baseline_ratio,
        results_by_phenotype=sensitivity_results,
        max_fold_change=max_fc,
        tissue_weight_source=weight_source,
    )


def run_validation_cases(
    tissue_weight_source: FluxTissueWeightSource | str = FluxTissueWeightSource.CURATED,
) -> None:
    """Run built-in validation cases and print a human-readable summary."""
    weight_source = _normalize_tissue_weight_source(tissue_weight_source)
    print("=" * 72)
    print("ExposoGraph Flux Engine — Validation Cases")
    print(f"Tissue weights: {weight_source.value}")
    print("=" * 72)

    baseline_pah = compute_pathway_flux(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "NM", "CYP1B1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        "Lung",
        0.1,
        tissue_weight_source=weight_source,
    )
    gstm1_null = compute_pathway_flux(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "null", "CYP1B1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        "Lung",
        0.1,
        tissue_weight_source=weight_source,
    )
    print("\n[CASE 1] PAH in Lung: CYP1A1 NM + GSTM1 null")
    print(f"  Baseline net ratio : {baseline_pah.net_ratio}")
    print(f"  GSTM1-null ratio   : {gstm1_null.net_ratio}")
    print(f"  Risk classification: {gstm1_null.risk_classification.value}")

    aldh2_wt = compute_pathway_flux(
        "Aldehyde",
        {"ALDH2": "*1/*1", "ALDH1A1": "NM", "ADH1B": "*1/*1"},
        "Liver",
        10.0,
        tissue_weight_source=weight_source,
    )
    aldh2_het = compute_pathway_flux(
        "Aldehyde",
        {"ALDH2": "*1/*2", "ALDH1A1": "NM", "ADH1B": "*1/*1"},
        "Liver",
        10.0,
        tissue_weight_source=weight_source,
    )
    print("\n[CASE 2] Aldehyde clearance: ALDH2 *1/*2 heterozygote")
    print(f"  Wildtype detox flux : {aldh2_wt.total_detox:.6f}")
    print(f"  Heterozygote detox  : {aldh2_het.total_detox:.6f}")
    print(f"  Risk classification : {aldh2_het.risk_classification.value}")

    um_null = compute_pathway_flux(
        "PAH",
        {"CYP1A1": "UM", "GSTM1": "null", "CYP1B1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        "Lung",
        0.1,
        tissue_weight_source=weight_source,
    )
    print("\n[CASE 3] PAH in Lung: CYP1A1 UM + GSTM1 null")
    print(f"  Net ratio           : {um_null.net_ratio}")
    print(f"  Risk classification : {um_null.risk_classification.value}")

    print("\n" + "=" * 72)
    print("Validation complete.")
    print("=" * 72)


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint compatible with the original standalone flux module."""
    parser = argparse.ArgumentParser(
        description="ExposoGraph metabolic flux engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m ExposoGraph.flux_engine --genotypes '{"CYP1A1":"NM","GSTM1":"null"}' --tissue Lung --carcinogen PAH
  python -m ExposoGraph.flux_engine --full-profile --genotypes '{"ALDH2":"*1/*2"}' --tissue Liver
  python -m ExposoGraph.flux_engine --sensitivity --gene GSTM1 --carcinogen PAH --tissue Lung
  python -m ExposoGraph.flux_engine --validate
        """,
    )
    parser.add_argument("--genotypes", type=str, default="{}", help="JSON string of gene-to-phenotype mappings")
    parser.add_argument("--tissue", type=str, default="Liver", help="Target tissue")
    parser.add_argument("--carcinogen", type=str, default=None, help="Carcinogen class to model")
    parser.add_argument("--concentration", type=float, default=None, help="Substrate concentration in uM")
    parser.add_argument("--full-profile", action="store_true", help="Compute all carcinogen classes")
    parser.add_argument("--validate", action="store_true", help="Run built-in validation cases")
    parser.add_argument("--sensitivity", action="store_true", help="Run sensitivity analysis for one gene")
    parser.add_argument("--gene", type=str, default=None, help="Gene to vary for sensitivity analysis")
    parser.add_argument("--output-json", action="store_true", help="Output JSON")
    parser.add_argument("--lifestyle", type=str, default="{}", help="JSON lifestyle/co-exposure flags for induction modeling")
    parser.add_argument("--induction-factors", type=str, default="{}", help="JSON per-gene explicit Vmax induction factors")
    parser.add_argument("--qivive", action="store_true", help="Apply MPPGL/organ-weight QIVIVE scaling to flux magnitudes")
    parser.add_argument("--mppgl", type=float, default=None, help="QIVIVE override: microsomal protein mg/g tissue")
    parser.add_argument("--organ-weight-g", type=float, default=None, help="QIVIVE override: organ weight in grams")
    parser.add_argument(
        "--tissue-weight-source",
        choices=[FluxTissueWeightSource.CURATED.value, FluxTissueWeightSource.GTEX.value],
        default=FluxTissueWeightSource.CURATED.value,
        help="Use curated source-parity tissue weights (default) or GTEx quantitative weights.",
    )

    args = parser.parse_args(argv)
    weight_source = _normalize_tissue_weight_source(args.tissue_weight_source)

    if args.validate:
        run_validation_cases(weight_source)
        return 0

    try:
        genotypes = json.loads(args.genotypes)
    except json.JSONDecodeError as exc:
        print(f"ERROR: Could not parse --genotypes JSON: {exc}", file=sys.stderr)
        return 1
    try:
        lifestyle = json.loads(args.lifestyle)
    except json.JSONDecodeError as exc:
        print(f"ERROR: Could not parse --lifestyle JSON: {exc}", file=sys.stderr)
        return 1
    try:
        induction_factors = json.loads(args.induction_factors)
    except json.JSONDecodeError as exc:
        print(f"ERROR: Could not parse --induction-factors JSON: {exc}", file=sys.stderr)
        return 1

    qivive_context = {}
    if args.mppgl is not None:
        qivive_context["mppgl_mg_per_g"] = args.mppgl
    if args.organ_weight_g is not None:
        qivive_context["organ_weight_g"] = args.organ_weight_g

    if args.sensitivity:
        if not args.gene:
            print("ERROR: --sensitivity requires --gene", file=sys.stderr)
            return 1
        result_obj: Any = sensitivity_analysis(
            args.carcinogen or "PAH",
            args.gene,
            args.tissue,
            args.concentration,
            tissue_weight_source=weight_source,
        )
    elif args.full_profile or not args.carcinogen:
        result_obj = compute_full_profile(
            genotypes,
            args.tissue,
            tissue_weight_source=weight_source,
            lifestyle=lifestyle,
            induction_factors=induction_factors,
            qivive=args.qivive,
            qivive_context=qivive_context or None,
        )
    else:
        result_obj = compute_pathway_flux(
            args.carcinogen,
            genotypes,
            args.tissue,
            args.concentration,
            tissue_weight_source=weight_source,
            lifestyle=lifestyle,
            induction_factors=induction_factors,
            qivive=args.qivive,
            qivive_context=qivive_context or None,
        )

    print(json.dumps(asdict(result_obj), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
