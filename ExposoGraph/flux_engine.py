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
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Mapping, TypeAlias, cast

from .flux_equations import (
    activation_detox_ratio,
    finite_nonnegative as _finite_nonnegative,
    hill_equation,
    michaelis_menten,
    saturating_flux,
    scaled_vmax,
    susceptibility_score_log2 as _equation_susceptibility_score_log2,
)

if TYPE_CHECKING:
    from .engine import FluxReaction, GraphEngine

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


@dataclass(frozen=True)
class KineticFluxApplication:
    """Internal record for applying one resolved kinetic modifier to a flux."""

    baseline_flux: float
    kinetic_modifier: float
    modified_flux: float
    applied_once: bool = True


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


# ── Parameter-source labels (output vocabulary; data served by GraphEngine) ──

_KINETIC_PARAMETER_SOURCE = "kinetic_parameters.json"
_PROXY_PARAMETER_SOURCE = "proxy_flux_parameters.json"
# Global output precision (significant figures) for flux entries, derived
# terms, and totals -- replaces the former per-class entry_round/total_round
# aggregation syntax. Significant figures (rather than decimal places) keep
# uniform relative precision across flux magnitudes spanning ~1e-5 to ~100.
_FLUX_SIGNIFICANT_FIGURES: int = 4


def _round_flux(value: float) -> float:
    """Round a flux output to ``_FLUX_SIGNIFICANT_FIGURES`` significant figures."""
    if value == 0 or not math.isfinite(value):
        return value
    exponent = math.floor(math.log10(abs(value)))
    return round(value, _FLUX_SIGNIFICANT_FIGURES - 1 - exponent)


# ── Core kinetic equations ─────────────────────────────────────────────────


# Michaelis-Menten, Hill, and finite-input equations live in flux_equations.py.


def apply_kinetic_modifier_once(baseline_flux: float, kinetic_modifier: float) -> KineticFluxApplication:
    """Apply one resolved kinetic modifier to a non-negative baseline flux."""

    baseline = _finite_nonnegative(baseline_flux, "baseline_flux")
    modifier = _finite_nonnegative(kinetic_modifier, "kinetic_modifier")
    return KineticFluxApplication(
        baseline_flux=baseline,
        kinetic_modifier=modifier,
        modified_flux=baseline * modifier,
    )


# ── Modifier functions ─────────────────────────────────────────────────────


def genotype_modifier(diplotype: str, gene: str) -> float:
    """Return Vmax scaling factor (0.0-2.0) based on metaboliser phenotype.

    Standard scale: PM=0.0, IM=0.5, NM=1.0, RM=1.5, UM=2.0.
    Special cases for ALDH2 heterozygotes, GSTM1/GSTT1 copy-number states,
    and cohort-specific CYP aliases used by the extension modules.

    Delegates to ``GraphEngine.get_genotype_modifier`` over the retained
    ``genotype_modifiers`` tables.

    Args:
        diplotype: Phenotype label (for example ``"NM"``, ``"null"``,
            or ``"*1/*2"``).
        gene: Gene name (for example ``"GSTM1"`` or ``"ALDH2"``).

    Returns:
        Scaling factor between 0.0 and 2.0.
    """
    return _get_flux_contract_engine().get_genotype_modifier(diplotype, gene)


def _proxy_genotype_modifier(diplotype: str, gene: str | None) -> float:
    """Return a silent genotype scaling factor for proxy-model terms."""
    return _get_flux_contract_engine().get_proxy_genotype_modifier(diplotype, gene)


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
    """Validate the requested source label and report the source actually used.

    Tissue weights now always come from the engine's GTEx expression
    table, so every valid request resolves to ``GTEX``; invalid labels
    still raise. The request parameter remains part of the public API
    for compatibility.
    """
    if isinstance(tissue_weight_source, FluxTissueWeightSource):
        return FluxTissueWeightSource.GTEX

    normalized = str(tissue_weight_source).strip().lower()
    if normalized in (FluxTissueWeightSource.CURATED.value, FluxTissueWeightSource.GTEX.value):
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

    Weights come from the engine's GTEx expression table
    (``tissue_expression_data_raw.json``, normalized per enzyme by its
    most-expressing tissue) -- the same source ``GraphEngine`` bakes onto
    enzyme nodes at ``load_reference_graph`` time. The hand-curated
    ``tissue_expression_weights`` table formerly shipped in
    ``kinetic_parameters.json`` was removed; ``tissue_weight_source`` is
    kept for API compatibility and no longer selects a different source.

    Args:
        gene: Gene/enzyme symbol (e.g. "CYP1A1").
        tissue: Tissue name (e.g. "Lung", "liver").
        tissue_weight_source: Ignored (retained for API compatibility).

    Returns:
        Expression weight between 0.0 and 1.0. Returns 0.2 when the gene
        is known but the tissue is not covered by the expression source,
        and 0.5 when the gene has no entry in the source (moderate
        default).
    """
    _normalize_tissue_weight_source(tissue_weight_source)
    tissue_key = _normalize_tissue(tissue)

    engine = _get_flux_contract_engine()
    weights = engine.get_tissue_expression(gene)
    if weights is None:
        # Gene not in the expression source -- moderate default
        return 0.5
    expression_tissue = _FLUX_TISSUE_TO_GTEX.get(tissue_key)
    if expression_tissue is None:
        # Tissue not covered by the expression source
        return 0.2
    weight = weights.get(expression_tissue)
    if weight is None:
        return 0.2
    return float(weight)


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
    return float(_get_flux_contract_engine().get_default_concentration(carcinogen_class))


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
    interaction_params: Mapping[str, Any] | None = None,
) -> dict[str, float]:
    """Resolve optional co-exposure induction inputs into per-enzyme Vmax folds."""
    resolved: dict[str, float] = {}

    if lifestyle:
        try:
            from .interaction_engine import enzyme_induction_modifier

            resolved.update(
                enzyme_induction_modifier(lifestyle, interaction_params=interaction_params).enzyme_folds
            )
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
        gene: _round_flux(factor)
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
        edata["induction_modifier"] = _round_flux(factor)
        if not math.isclose(factor, 1.0):
            edata["flux"] = _round_flux(new_flux)
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
            result[total_name] = float(result[total_name]) * new_sum / old_sum

    return result


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
    return _get_flux_contract_engine().get_qivive_context(_normalize_tissue(tissue), overrides)


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
            edata["qivive_scale"] = _round_flux(scale)
            if "flux" in edata:
                edata["flux"] = _round_flux(float(edata["flux"]) * scale)

    for key in ("total_activation", "total_detox"):
        if key in result:
            result[key] = float(result[key]) * scale

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
    return _equation_susceptibility_score_log2(net_ratio)


def _round_steady_state_value(value: float) -> float:
    """Clamp invalid steady-state outputs; round valid ones to flux precision."""
    if not math.isfinite(value) or value < 0:
        return 0.0
    return _round_flux(value)


def _steady_state_context_for_tissue(
    tissue: str,
    overrides: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Return validated defaults for the flux-coupled steady-state solver."""
    return _get_flux_contract_engine().get_steady_state_context(_normalize_tissue(tissue), overrides)


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
    engine = _get_flux_contract_engine()
    if gene and engine.get_tissue_expression(gene) is not None:
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
    vmax = scaled_vmax(
        float(term["vmax"]),
        gm,
        tw,
        vmax_relative=float(term.get("vmax_relative", 1.0)),
        relative_capacity_scale=_relative_capacity_scale(term.get("relative_capacity")),
    )
    flux = michaelis_menten(S, vmax, float(term["km"]))
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
    flux = saturating_flux(S, float(term["scale"]) * gm * tw, float(term["km"]))
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
    """Return the proxy flux config for a class (served by GraphEngine)."""
    cfg = _get_flux_contract_engine().get_flux_class_config(class_name)
    if not cfg:
        raise KeyError(class_name)
    return cast(JsonDict, cfg)


_FLUX_CONTRACT_ENGINE: "GraphEngine | None" = None


def _get_flux_contract_engine() -> "GraphEngine":
    """Return the module-level bare engine serving the flux-reaction contract.

    A bare GraphEngine is enough: ``get_flux_reactions`` and friends lazily
    build the side index from the bundled parameter files without loading
    the full graph, and ``get_edge_flux_reactions`` degrades to the index
    when no graph is loaded. Callers that already hold a reference engine
    (with the graph loaded and edge kinetics baked) pass it explicitly
    through ``compute_pathway_flux(engine=...)`` so the mechanistic reads
    walk ``Edge.kinetics`` instead.
    """
    global _FLUX_CONTRACT_ENGINE
    if _FLUX_CONTRACT_ENGINE is None:
        from .engine import GraphEngine

        _FLUX_CONTRACT_ENGINE = GraphEngine()
    return _FLUX_CONTRACT_ENGINE


# Output labels for the proxy kinetics field, keyed by rate law.
_PROXY_TERM_KINETICS_LABELS: dict[str, str] = {
    "michaelis_menten": "semi_quantitative",
    "hill": "damage_proxy",
    "saturating": "semi_quantitative",
    "repair": "repair_proxy",
}

# Legacy output-shape compatibility: the hand-written NDEA and VinylChloride
# functions omitted the "kinetics" label on activation entries, so
# _enzyme_flux_from_dict defaulted those to "michaelis_menten". Preserved so
# the generic loop is output-identical; drop when output normalization is
# accepted as a deliberate change.
_PROXY_ACTIVATION_WITHOUT_KINETICS = frozenset({"NDEA", "VinylChloride"})


def _compute_generic_proxy_flux(
    carcinogen_class: str,
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
    *,
    engine: "GraphEngine | None" = None,
) -> FluxResultDict:
    """Compute a semi-quantitative proxy flux from engine FluxReaction records.

    Replaces the six hand-written proxy-class functions (AromaticAmines,
    EstrogenMetabolites, NDEA, VinylChloride, UV_Radiation, HeavyMetal).
    Enzyme scope, role, and rate law come from
    ``GraphEngine.get_flux_reactions`` -- the flux contract -- and the
    per-term parameters ride on the reaction records. Dioxin
    (receptor/signaling model) and ChlorinatedSolvent (derived clearance)
    keep their dedicated functions: their model structures are not
    term-sum proxies.
    """
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    reactions = active_engine.get_flux_reactions(carcinogen_class)
    cfg = _get_proxy_class_params(carcinogen_class)

    def _entry(reaction: "FluxReaction", value: float, gm: float, tw: float) -> JsonDict:
        entry: JsonDict = {
            "flux": _round_flux(value),
            "genotype_modifier": gm,
            "tissue_weight": tw,
            "confidence": reaction.confidence,
        }
        if reaction.role != "activation" or carcinogen_class not in _PROXY_ACTIVATION_WITHOUT_KINETICS:
            entry["kinetics"] = _PROXY_TERM_KINETICS_LABELS.get(reaction.rate_law, "semi_quantitative")
        entry["note"] = str(reaction.params.get("note", ""))
        return entry

    def _evaluate(reaction: "FluxReaction", activation_base: float | None = None) -> tuple[float, float, float]:
        term = reaction.params
        if reaction.rate_law == "repair":
            if activation_base is None:
                raise ValueError(f"Repair term {reaction.term_key} evaluated without an activation base")
            return _compute_proxy_repair_term(activation_base, term, genotypes, tissue, tissue_weight_source)
        if reaction.rate_law == "michaelis_menten":
            return _compute_proxy_mm_term(term, genotypes, tissue, S, tissue_weight_source)
        if reaction.rate_law == "hill":
            return _compute_proxy_hill_term(term, genotypes, tissue, S, tissue_weight_source)
        if reaction.rate_law == "saturating":
            return _compute_proxy_saturating_term(term, genotypes, tissue, S, tissue_weight_source)
        raise ValueError(
            f"Unsupported proxy rate law {reaction.rate_law!r} for "
            f"{carcinogen_class}/{reaction.term_key}"
        )

    activation_enzymes: dict[str, Any] = {}
    total_activation = 0.0
    for reaction in reactions:
        if reaction.role != "activation":
            continue
        value, gm, tw = _evaluate(reaction)
        total_activation += value
        activation_enzymes[reaction.term_key] = _entry(reaction, value, gm, tw)

    detox_enzymes: dict[str, Any] = {}
    total_detox = 0.0
    for reaction in reactions:
        if reaction.role == "activation":
            continue
        value, gm, tw = _evaluate(reaction, activation_base=total_activation)
        total_detox += value
        detox_enzymes[reaction.term_key] = _entry(reaction, value, gm, tw)

    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": detox_enzymes,
        "total_activation": total_activation,
        "total_detox": total_detox,
        "unit_note": cfg["unit_note"],
    }


def _compute_generic_mechanistic_flux(
    carcinogen_class: str,
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
    *,
    engine: "GraphEngine | None" = None,
) -> FluxResultDict:
    """Compute a measured-kinetics flux from engine records and an aggregation spec.

    Replaces the hand-written PAH, Nitrosamine, NDMA, HCA, and Benzene
    functions. Direct terms (roster entries without a ``dormant`` flag)
    come from ``GraphEngine.get_edge_flux_reactions`` -- term parameters
    are read off the substrate→enzyme edges' ``kinetics["flux_terms"]``
    payloads, with the side index supplying only the roster order and the
    handful of terms not yet edge-anchored; how they combine --
    per-role vmax fields, derived terms, efficiency scaling, detox
    fractions, total scaling, rounding, and unit notes -- comes from the
    class's ``aggregation`` block in kinetic_parameters.json via
    ``GraphEngine.get_flux_aggregation``. Aflatoxin and Aldehyde keep
    dedicated functions (their aggregation blocks carry a
    ``dedicated_function`` status note).
    """
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    reactions = active_engine.get_edge_flux_reactions(carcinogen_class)
    agg = active_engine.get_flux_aggregation(carcinogen_class)
    vmax_fields: dict[str, Any] = dict(agg.get("vmax_field", {}))
    km_field = str(agg.get("km_field", "Km_uM"))
    derived_specs: dict[str, Any] = dict(agg.get("derived_terms", {}))
    entry_names: dict[str, str] = dict(agg.get("entry_names", {}))
    efficiency_spec: JsonDict | None = agg.get("activation_efficiency_term")

    def _gene(reaction: "FluxReaction") -> str:
        gene = reaction.params.get("gene")
        if isinstance(gene, str) and gene:
            return gene
        return str(reaction.enzyme_id or reaction.term_key)

    def _entry(value: float, gm: float, tw: float, role: str, confidence: str) -> JsonDict:
        return {
            "flux": _round_flux(value),
            "genotype_modifier": gm,
            "tissue_weight": tw,
            "confidence": confidence,
        }

    def _evaluate(reaction: "FluxReaction", role: str) -> tuple[float, float, float]:
        term = reaction.params
        gene = _gene(reaction)
        gm = genotype_modifier(genotypes.get(gene, "NM"), gene)
        tw = get_flux_tissue_weight(gene, tissue, tissue_weight_source)
        vmax_field = vmax_fields.get(role) or vmax_fields.get("activation")
        if vmax_field is None:
            raise ValueError(
                f"No vmax field for {carcinogen_class}/{reaction.term_key} role {role!r}"
            )
        vmax = float(term[vmax_field]) * gm * tw
        km = float(term[km_field])
        if reaction.rate_law == "hill":
            return hill_equation(S, vmax, km, float(term.get("hill_n", 1.0))), gm, tw
        # michaelis_menten and clint_normalized (kcat-form GST) terms share
        # the MM expression Vmax * S / (Km + S).
        return michaelis_menten(S, vmax, km), gm, tw

    activation_entries: dict[str, Any] = {}
    detox_entries: dict[str, Any] = {}
    computed: dict[str, tuple[float, float, float]] = {}
    term_roles: dict[str, str] = {}

    # Pass 1: direct terms (roster entries the model evaluates as-is).
    for reaction in reactions:
        if reaction.term_key in derived_specs or reaction.params.get("dormant"):
            continue
        if reaction.rate_law == "clint_ceiling":
            if efficiency_spec is None or efficiency_spec.get("term") != reaction.term_key:
                raise ValueError(
                    f"clint_ceiling term {carcinogen_class}/{reaction.term_key} "
                    "has no activation_efficiency_term spec"
                )
            continue
        role = reaction.role
        if role not in ("activation", "detoxification"):
            continue
        value, gm, tw = _evaluate(reaction, role)
        computed[reaction.term_key] = (value, gm, tw)
        term_roles[reaction.term_key] = role
        entries = activation_entries if role == "activation" else detox_entries
        entries[entry_names.get(reaction.term_key, reaction.term_key)] = _entry(
            value, gm, tw, role, reaction.confidence
        )

    # Pass 2: derived terms (flux expressed from an already-evaluated base).
    for name, spec in derived_specs.items():
        base_key = str(spec["base_term"])
        if base_key not in computed:
            raise ValueError(f"Derived term {name!r} has unevaluated base {base_key!r}")
        base_value, base_gm, base_tw = computed[base_key]
        fraction = float(spec["fraction"])
        gm = genotype_modifier(genotypes.get(name, "NM"), name)
        form = str(spec.get("form", "fraction_of_base"))
        if form == "fraction_of_base":
            tw = 1.0
            value = base_value * fraction * gm
        elif form == "relative_to_base":
            tw = get_flux_tissue_weight(name, tissue, tissue_weight_source)
            value = base_value * fraction * gm * tw / (base_gm * base_tw + 1e-9)
        else:
            raise ValueError(f"Unsupported derived form {form!r} for {name}")
        role = next((r.role for r in reactions if r.term_key == name), None)
        if role not in ("activation", "detoxification"):
            raise ValueError(f"Derived term {name!r} has no activation/detoxification role")
        computed[name] = (value, gm, tw)
        term_roles[name] = role
        entries = activation_entries if role == "activation" else detox_entries
        entries[entry_names.get(name, name)] = _entry(
            value, gm, tw, role, str(spec.get("confidence", ""))
        )

    # Totals: raw sums in entry order, then class-level aggregation.
    total_activation = 0.0
    total_detox = 0.0
    for key, (value, _, _) in computed.items():
        if term_roles[key] == "activation":
            total_activation += value
        else:
            total_detox += value

    if efficiency_spec is not None:
        eff_reaction = next(r for r in reactions if r.term_key == efficiency_spec["term"])
        gene = _gene(eff_reaction)
        eff_gm = genotype_modifier(genotypes.get(gene, "NM"), gene)
        eff_tw = get_flux_tissue_weight(gene, tissue, tissue_weight_source)
        efficiency = min(
            1.0,
            (float(eff_reaction.params["CLint"]) / float(efficiency_spec["reference_CLint"]))
            * eff_gm
            * eff_tw,
        )
        total_activation = total_activation * efficiency

    detox_fraction = agg.get("detox_fraction_of_activation")
    if detox_fraction is not None:
        total_detox = total_activation * float(detox_fraction)
        detox_entry_spec = agg.get("detox_entry")
        if detox_entry_spec:
            detox_entries[str(detox_entry_spec["name"])] = {
                "flux": total_detox,
                "genotype_modifier": 1.0,
                "tissue_weight": 1.0,
                "confidence": str(detox_entry_spec.get("confidence", "estimated")),
            }

    detox_scale = agg.get("detox_total_scale")
    if detox_scale is not None:
        total_detox = total_detox * float(detox_scale)

    result: FluxResultDict = {
        "activation_enzymes": activation_entries,
        "detox_enzymes": detox_entries,
        "total_activation": total_activation,
        "total_detox": total_detox,
    }
    if efficiency_spec is not None:
        result[str(efficiency_spec["output_field"])] = _round_flux(efficiency)
    result["unit_note"] = str(agg.get("unit_note", ""))
    return result


def _class_parameter_metadata(carcinogen_class: str) -> dict[str, str]:
    """Return class-level parameter metadata for measured or proxy models."""
    proxy_cfg = _get_flux_contract_engine().get_flux_class_config(carcinogen_class)
    if not proxy_cfg:
        return {
            "model_kind": "measured_kinetics",
            "parameter_source": _KINETIC_PARAMETER_SOURCE,
        }
    return {
        "model_kind": proxy_cfg["model_kind"],
        "parameter_source": _PROXY_PARAMETER_SOURCE,
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
            entry = _get_flux_contract_engine().get_flux_provenance_entry(ref)
            sources = list(entry.get("sources", []))
            basis = entry.get("parameter_basis", "")

        return {
            "model_kind": cfg["model_kind"],
            "parameter_source": _PROXY_PARAMETER_SOURCE,
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

    proxy_mode = class_meta["parameter_source"] == _PROXY_PARAMETER_SOURCE
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


def _compute_aflatoxin_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
    *,
    engine: "GraphEngine | None" = None,
) -> FluxResultDict:
    """Compute AFB1 activation and detoxification fluxes."""
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    terms = {r.term_key: r for r in active_engine.get_edge_flux_reactions("Aflatoxin")}
    agg = active_engine.get_flux_aggregation("Aflatoxin")
    activation_enzymes: dict[str, Any] = {}
    detox_enzymes: dict[str, Any] = {}

    # CYP3A4 — Hill kinetics
    cyp3a4_p = terms["CYP3A4"].params
    gm3a4 = genotype_modifier(genotypes.get("CYP3A4", "NM"), "CYP3A4")
    tw3a4 = get_flux_tissue_weight("CYP3A4", tissue, tissue_weight_source)
    v_cyp3a4 = hill_equation(
        S,
        cyp3a4_p["Vmax_pmol_min_pmolP450"] * gm3a4 * tw3a4,
        cyp3a4_p["Km_uM"],
        cyp3a4_p["hill_n"],
    )
    activation_enzymes["CYP3A4"] = {
        "flux": _round_flux(v_cyp3a4),
        "kinetics": "hill",
        "n": cyp3a4_p["hill_n"],
        "genotype_modifier": gm3a4,
        "tissue_weight": tw3a4,
        "confidence": terms["CYP3A4"].confidence,
        "fraction_contribution": cyp3a4_p["fraction_contribution"],
    }

    # CYP1A2 — Michaelis-Menten
    cyp1a2_p = terms["CYP1A2"].params
    gm1a2 = genotype_modifier(genotypes.get("CYP1A2", "NM"), "CYP1A2")
    tw1a2 = get_flux_tissue_weight("CYP1A2", tissue, tissue_weight_source)
    v_cyp1a2 = michaelis_menten(
        S,
        cyp1a2_p["Vmax_pmol_min_pmolP450"] * gm1a2 * tw1a2,
        cyp1a2_p["Km_uM"],
    )
    activation_enzymes["CYP1A2"] = {
        "flux": _round_flux(v_cyp1a2),
        "kinetics": "michaelis_menten",
        "genotype_modifier": gm1a2,
        "tissue_weight": tw1a2,
        "confidence": terms["CYP1A2"].confidence,
        "fraction_contribution": cyp1a2_p["fraction_contribution"],
    }

    total_activation = v_cyp3a4 + v_cyp1a2

    # CYP3A4 AFQ1 (detox)
    afq1_p = terms["CYP3A4_AFQ1"].params
    v_afq1 = hill_equation(
        S,
        afq1_p["Vmax_pmol_min_pmolP450"] * gm3a4 * tw3a4,
        afq1_p["Km_uM"],
        afq1_p["hill_n"],
    )
    detox_enzymes["CYP3A4_AFQ1"] = {
        "flux": _round_flux(v_afq1),
        "genotype_modifier": gm3a4,
        "tissue_weight": tw3a4,
        "confidence": terms["CYP3A4_AFQ1"].confidence,
    }

    # GSTA1 (estimated)
    gsta1_p = terms["GSTA1"].params
    gsta1_gm = genotype_modifier(genotypes.get("GSTA1", "NM"), "GSTA1")
    gsta1_tw = get_flux_tissue_weight(
        gsta1_p.get("tissue_weight_gene", "GSTA1"), tissue, tissue_weight_source
    )
    v_gsta1 = michaelis_menten(
        S,
        gsta1_p["Vmax_pmol_min_pmolP450_estimated"] * gsta1_gm * gsta1_tw,
        gsta1_p["Km_uM"],
    )
    detox_enzymes["GSTA1_conjugation"] = {
        "flux": _round_flux(v_gsta1),
        "genotype_modifier": gsta1_gm,
        "tissue_weight": gsta1_tw,
        "confidence": terms["GSTA1"].confidence,
    }

    total_detox = v_afq1 + v_gsta1

    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": detox_enzymes,
        "total_activation": total_activation,
        "total_detox": total_detox,
        "unit_note": agg["unit_note"],
    }


def _compute_aldehyde_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
    *,
    engine: "GraphEngine | None" = None,
) -> FluxResultDict:
    """Compute aldehyde (acetaldehyde) clearance flux."""
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    terms = {r.term_key: r for r in active_engine.get_edge_flux_reactions("Aldehyde")}
    agg = active_engine.get_flux_aggregation("Aldehyde")
    detox_enzymes: dict[str, Any] = {}

    # Determine ALDH2 genotype
    aldh2_gt = genotypes.get("ALDH2", "*1/*1")
    if aldh2_gt in ("*1/*1", "NM", "WT", "wildtype"):
        aldh2_term = terms["ALDH2_star1"]
        aldh2_p = aldh2_term.params
        aldh2_gm = 1.0
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]
    elif aldh2_gt in ("*1/*2", "heterozygote"):
        aldh2_term = terms["ALDH2_star1"]
        aldh2_p = aldh2_term.params
        aldh2_gm = genotype_modifier("*1/*2", "ALDH2")  # 0.25
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]
    elif aldh2_gt in ("*2/*2", "PM", "PM_ALDH2"):
        aldh2_term = terms["ALDH2_star2_homozygous"]
        aldh2_p = aldh2_term.params
        aldh2_gm = 1.0  # parameters already reflect variant
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]
    else:
        aldh2_term = terms["ALDH2_star1"]
        aldh2_p = aldh2_term.params
        aldh2_gm = genotype_modifier(aldh2_gt, "ALDH2")
        aldh2_km = aldh2_p["Km_uM"]
        aldh2_vmax = aldh2_p["Vmax_U_per_mg"]

    tw_aldh2 = get_flux_tissue_weight("ALDH2", tissue, tissue_weight_source)
    v_aldh2 = michaelis_menten(S, aldh2_vmax * aldh2_gm * tw_aldh2, aldh2_km)
    detox_enzymes["ALDH2"] = {
        "flux": _round_flux(v_aldh2),
        "genotype": aldh2_gt,
        "genotype_modifier": aldh2_gm,
        "tissue_weight": tw_aldh2,
        "CLint": _round_flux((aldh2_vmax * aldh2_gm * tw_aldh2) / aldh2_km),
        "confidence": aldh2_term.confidence,
    }

    # ALDH1A1 (backup)
    aldh1a1_p = terms["ALDH1A1"].params
    aldh1a1_gm = genotype_modifier(genotypes.get("ALDH1A1", "NM"), "ALDH1A1")
    tw_aldh1a1 = get_flux_tissue_weight("ALDH1A1", tissue, tissue_weight_source)
    v_aldh1a1 = michaelis_menten(
        S, aldh1a1_p["Vmax_U_per_mg"] * aldh1a1_gm * tw_aldh1a1, aldh1a1_p["Km_uM"]
    )
    detox_enzymes["ALDH1A1"] = {
        "flux": _round_flux(v_aldh1a1),
        "genotype_modifier": aldh1a1_gm,
        "tissue_weight": tw_aldh1a1,
        "confidence": terms["ALDH1A1"].confidence,
    }

    total_detox = v_aldh2 + v_aldh1a1

    # Ethanol -> Acetaldehyde production
    adh_gt = genotypes.get("ADH1B", "*1/*1")
    if adh_gt in ("*2/*2", "fast", "RM"):
        adh_params = terms["ADH1B_star2"].params
        adh_confidence = terms["ADH1B_star2"].confidence
    else:
        adh_params = terms["ADH1B_star1"].params
        adh_confidence = terms["ADH1B_star1"].confidence

    eth_conc = active_engine.get_flux_metadata()["exposure_defaults_uM"]["ethanol"]
    v_adh = michaelis_menten(
        eth_conc, adh_params["Vmax_U_per_mg"], adh_params["Km_uM"]
    )

    return {
        "activation_enzymes": {
            "ADH1B": {
                "reaction": "Ethanol -> Acetaldehyde",
                "flux": _round_flux(v_adh),
                "genotype": adh_gt,
                "genotype_modifier": 1.0,
                "tissue_weight": 1.0,
                "confidence": adh_confidence,
            }
        },
        "detox_enzymes": detox_enzymes,
        "total_activation": v_adh,
        "total_detox": total_detox,
        "unit_note": agg["unit_note"],
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
        "flux": _round_flux(v_oxidation),
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
        "flux": _round_flux(v_gsh),
        "genotype_modifier": gm_gstt1,
        "tissue_weight": tw_gstt1,
        "confidence": gsh_p["confidence"],
        "note": gsh_p["note"],
    }

    detox_p = cfg["detox_terms"]["non_genotoxic_clearance_proxy"]
    v_clearance = v_oxidation * float(detox_p["scale"])
    detox_enzymes = {
        "non_genotoxic_clearance_proxy": {
            "flux": _round_flux(v_clearance),
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
        "total_activation": total_activation,
        "total_detox": v_clearance,
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
                "flux": _round_flux(v1a1),
                "genotype_modifier": gm1a1,
                "tissue_weight": tw1a1,
                "confidence": act1_p["confidence"],
                "kinetics": "receptor_mediated",
                "note": act1_p["note"],
            },
            "CYP1B1": {
                "flux": _round_flux(v1b1),
                "genotype_modifier": gm1b1,
                "tissue_weight": tw1b1,
                "confidence": act2_p["confidence"],
                "kinetics": "receptor_mediated",
                "note": act2_p["note"],
            },
        },
        "detox_enzymes": {
            "AHRR_feedback": {
                "flux": _round_flux(v_feedback),
                "genotype_modifier": 1.0,
                "tissue_weight": tissue_signal,
                "confidence": feedback_p["confidence"],
                "kinetics": "receptor_mediated",
                "note": feedback_p["note"],
            }
        },
        "total_activation": v1a1 + v1b1,
        "total_detox": v_feedback,
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
        parameter_source=str(d.get("parameter_source", _KINETIC_PARAMETER_SOURCE)),
        provenance_ref=str(d.get("provenance_ref", "")),
        provenance_sources=[str(source) for source in d.get("provenance_sources", [])],
        parameter_basis=str(d.get("parameter_basis", "")),
    )


# ── Public API ─────────────────────────────────────────────────────────────

# Proxy classes routed through the generic engine-contract loop.
# Measured-kinetics classes routed through the generic mechanistic loop are
# listed in _GENERIC_MECHANISTIC_FLUX_CLASSES below. Aflatoxin and Aldehyde
# (in _DEDICATED_ENGINE_FLUX_CLASSES) keep hand-written control flow but
# source their term parameters and aggregation blocks through the engine
# contract. The remaining entries in _DISPATCH -- Dioxin and
# ChlorinatedSolvent -- keep fully hand-written functions whose model
# structures are not term-sum proxies.
_GENERIC_PROXY_FLUX_CLASSES = frozenset(
    {"AromaticAmines", "EstrogenMetabolites", "NDEA", "VinylChloride", "UV_Radiation", "HeavyMetal"}
)

# Measured-kinetics classes whose term evaluation and aggregation are fully
# described by the per-class "aggregation" block in kinetic_parameters.json.
_GENERIC_MECHANISTIC_FLUX_CLASSES = frozenset({"PAH", "Nitrosamine", "NDMA", "HCA", "Benzene"})

# Dedicated functions that still take the engine kwarg for their parameters.
_DEDICATED_ENGINE_FLUX_CLASSES = frozenset({"Aflatoxin", "Aldehyde"})

_DISPATCH: dict[str, FluxCalculator] = {
    "PAH": partial(_compute_generic_mechanistic_flux, "PAH"),
    "Aflatoxin": _compute_aflatoxin_flux,
    "Aldehyde": _compute_aldehyde_flux,
    "Nitrosamine": partial(_compute_generic_mechanistic_flux, "Nitrosamine"),
    "NDMA": partial(_compute_generic_mechanistic_flux, "NDMA"),
    "NDEA": partial(_compute_generic_proxy_flux, "NDEA"),
    "HCA": partial(_compute_generic_mechanistic_flux, "HCA"),
    "AromaticAmines": partial(_compute_generic_proxy_flux, "AromaticAmines"),
    "EstrogenMetabolites": partial(_compute_generic_proxy_flux, "EstrogenMetabolites"),
    "Benzene": partial(_compute_generic_mechanistic_flux, "Benzene"),
    "VinylChloride": partial(_compute_generic_proxy_flux, "VinylChloride"),
    "ChlorinatedSolvent": _compute_chlorinated_solvent_flux,
    "UV_Radiation": partial(_compute_generic_proxy_flux, "UV_Radiation"),
    "Dioxin": _compute_dioxin_flux,
    "HeavyMetal": partial(_compute_generic_proxy_flux, "HeavyMetal"),
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
    engine: "GraphEngine | None" = None,
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
        engine: Optional :class:`~ExposoGraph.engine.GraphEngine`. When
            given, lifestyle-driven induction folds are resolved against
            the engine's interaction-parameter document instead of the
            interaction engine's own file loader.

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

    resolved_induction = _resolve_induction_factors(
        lifestyle,
        induction_factors,
        interaction_params=engine.get_interaction_parameters() if engine is not None else None,
    )
    if (
        cls_str in _GENERIC_PROXY_FLUX_CLASSES
        or cls_str in _GENERIC_MECHANISTIC_FLUX_CLASSES
        or cls_str in _DEDICATED_ENGINE_FLUX_CLASSES
    ):
        result = _DISPATCH[cls_str](genotypes, tissue, substrate_conc_uM, weight_source, engine=engine)
    else:
        result = _DISPATCH[cls_str](genotypes, tissue, substrate_conc_uM, weight_source)
    result = _annotate_flux_result_metadata(cls_str, result)
    result = _apply_induction_modifiers(result, resolved_induction)

    qivive_used_context: dict[str, float] = {}
    if qivive:
        qivive_used_context = _qivive_context_for_tissue(tissue, qivive_context)
        result = _apply_qivive_scale(result, qivive_used_context)

    act = float(result["total_activation"])
    det = float(result["total_detox"])
    net_ratio = activation_detox_ratio(act, det)

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
        total_activation=_round_flux(act),
        total_detox=_round_flux(det),
        net_ratio=_round_flux(net_ratio),
        susceptibility_score_log2=susceptibility_score,
        risk_classification=risk,
        tissue_weight_source=weight_source,
        model_kind=result.get("model_kind", "measured_kinetics"),
        parameter_source=result.get("parameter_source", _KINETIC_PARAMETER_SOURCE),
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
