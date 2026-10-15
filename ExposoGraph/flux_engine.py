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

    PAH = "group_pahs"
    AFLATOXIN = "group_mycotoxins"
    ALDEHYDE = "group_aldehydes"
    NITROSAMINE = "group_nitrosamines"
    NDMA = "NDMA"
    NDEA = "NDEA"
    HCA = "group_hcas"
    AROMATIC_AMINES = "group_aromatic_amines"
    ESTROGEN_METABOLITES = "group_estrogen"
    BENZENE = "Benzene"
    VINYL_CHLORIDE = "VinylChloride"
    CHLORINATED_SOLVENT = "group_chlorinated_solvents"
    UV_RADIATION = "group_uv_radiation"
    DIOXIN = "group_dioxins"
    HEAVY_METAL = "group_heavy_metals"


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


def apply_kinetic_modifier_once(
    baseline_flux: float, kinetic_modifier: float
) -> KineticFluxApplication:
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

    Delegates to ``GraphEngine.get_genotype_modifier``, which reads
    ``genotypes.json``.

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


def _normalize_tissue(tissue: str) -> str:
    """Lowercase a tissue string. Vocabulary keys are already lowercase."""
    return tissue.lower().strip()


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
    if normalized in (
        FluxTissueWeightSource.CURATED.value,
        FluxTissueWeightSource.GTEX.value,
    ):
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
    weight = weights.get(tissue_key)
    if weight is None:
        # Tissue not covered by the expression source
        return 0.2
    return float(weight)


def tissue_weight(
    gene: str,
    tissue: str,
    tissue_weight_source: FluxTissueWeightSource | str = FluxTissueWeightSource.CURATED,
) -> float:
    """Source-compatible alias for :func:`get_flux_tissue_weight`."""
    return get_flux_tissue_weight(
        gene, tissue, tissue_weight_source=tissue_weight_source
    )


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
    elif net_ratio < 4.0:
        return RiskClassification.ELEVATED
    else:
        return RiskClassification.HIGH


def _classify_risk(net_ratio: float) -> str:
    """Source-compatible string wrapper around :func:`classify_risk`."""
    return classify_risk(net_ratio).value


# ── Helper ─────────────────────────────────────────────────────────────────


def _get_default_concentration(carcinogen_class: str) -> float:
    """Return default environmental exposure concentration in uM."""
    return float(
        _get_flux_contract_engine().get_default_concentration(carcinogen_class)
    )


_KNOWN_FLUX_GENE_PREFIXES = tuple(
    sorted(
        {
            "ABCB1",
            "ABCC2",
            "ABCG2",
            "ADH1B",
            "ADH5",
            "AHRR",
            "ALDH1A1",
            "ALDH2",
            "AS3MT",
            "COMT",
            "CYP1A1",
            "CYP1A2",
            "CYP1B1",
            "CYP2A6",
            "CYP2A13",
            "CYP2D6",
            "CYP2E1",
            "CYP3A4",
            "EPHX1",
            "ERCC2",
            "GSTA1",
            "GSTM1",
            "GSTP1",
            "GSTT1",
            "MGMT",
            "NAT1",
            "NAT2",
            "NQO1",
            "OGG1",
            "POLH",
            "SULT1E1",
            "UGT2B7",
            "XPC",
            "XRCC1",
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
                enzyme_induction_modifier(
                    lifestyle, interaction_params=interaction_params
                ).enzyme_folds
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
        gene: factor
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
        old_sum, new_sum = _rescale_flux_section_for_induction(
            enzymes, induction_factors
        )
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
    return _get_flux_contract_engine().get_qivive_context(
        _normalize_tissue(tissue), overrides
    )


def _apply_qivive_scale(
    result: FluxResultDict, qivive_context: Mapping[str, float]
) -> FluxResultDict:
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
    return _get_flux_contract_engine().get_steady_state_context(
        _normalize_tissue(tissue), overrides
    )


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
        intrinsic_clearance_l_per_day
        / (tissue_flow_l_per_day + intrinsic_clearance_l_per_day)
        if tissue_flow_l_per_day + intrinsic_clearance_l_per_day > 0
        else 0.0
    )
    tissue_clearance_l_per_day = tissue_flow_l_per_day * extraction_ratio
    background_clearance_l_per_day = background_clearance_rate * central_volume_l
    total_clearance_l_per_day = (
        background_clearance_l_per_day + tissue_clearance_l_per_day
    )

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
        partition
        * central_conc
        * tissue_flow_l_per_day
        / (tissue_flow_l_per_day + intrinsic_clearance_l_per_day)
        if tissue_flow_l_per_day + intrinsic_clearance_l_per_day > 0
        else partition * central_conc
    )
    reactive_loss_rate = (
        solver_context["reactive_intermediate_loss_rate_per_day"] + detox_rate
    )
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

    central_rate = (
        total_clearance_l_per_day / central_volume_l if central_volume_l > 0 else 0.0
    )
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
        "tissue_clearance_l_per_day": _round_steady_state_value(
            tissue_clearance_l_per_day
        ),
        "total_clearance_l_per_day": _round_steady_state_value(
            total_clearance_l_per_day
        ),
        "extraction_ratio": _round_steady_state_value(extraction_ratio),
        "time_to_steady_state_days": _round_steady_state_value(
            time_to_steady_state_days
        ),
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
    steady_state = solve_flux_steady_state(
        substrate_conc_uM, act, det, tissue, context=context
    )
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


def _active_proxy_class_params(
    class_name: str, engine: "GraphEngine | None"
) -> JsonDict:
    """Return the proxy flux config for a class via the active engine.

    A loaded engine serves the class block graph-first (CarcinogenGroup
    ``flux_class_metadata`` carrier, JSON fallback per field); a bare
    engine (``None``) uses the module-level JSON-only contract engine.
    """
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    cfg = active_engine.get_flux_class_config(class_name)
    if not cfg:
        raise KeyError(class_name)
    return cast(JsonDict, cfg)


def _get_proxy_class_params(class_name: str) -> JsonDict:
    """Return the proxy flux config for a class (bare JSON-served engine)."""
    return _active_proxy_class_params(class_name, None)


_FLUX_CONTRACT_ENGINE: "GraphEngine | None" = None


def _induction_params_from_rules(rules) -> dict[str, Any]:
    """Rebuild the induction block the lifestyle modifier still reads.

    The modifier wants ``enzyme_induction`` keyed by context, then enzyme.
    The rules already carry that context, the enzyme, and the original
    record in ``metadata``. This does not add a selection rule.
    """
    sections: dict[str, dict] = {}
    for rule in rules:
        sections.setdefault(rule.exposure_context, {})[rule.enzyme] = dict(rule.metadata or {})
    return {"enzyme_induction": sections}


def _get_flux_contract_engine() -> "GraphEngine":
    """Return the loaded revised engine that serves the flux contract."""
    global _FLUX_CONTRACT_ENGINE
    if _FLUX_CONTRACT_ENGINE is None:
        from .engine_revised import GraphEngine

        engine = GraphEngine()
        engine.load_reference_graph()
        _FLUX_CONTRACT_ENGINE = engine
    return _FLUX_CONTRACT_ENGINE


# Output labels for the proxy kinetics field, keyed by rate law.
_PROXY_TERM_KINETICS_LABELS: dict[str, str] = {
    "michaelis_menten": "semi_quantitative",
    "hill": "damage_proxy",
    "saturating": "semi_quantitative",
    "repair": "repair_proxy",
    "derived_scale": "semi_quantitative",
}

# Legacy output-shape compatibility: the hand-written NDEA, VinylChloride,
# and ChlorinatedSolvent functions omitted the "kinetics" label on
# activation entries, so _enzyme_flux_from_dict defaulted those to
# "michaelis_menten". Preserved so the generic loop is output-identical;
# drop when output normalization is accepted as a deliberate change.
_PROXY_ACTIVATION_WITHOUT_KINETICS = frozenset(
    {"NDEA", "VinylChloride", "group_chlorinated_solvents"}
)


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

    Replaces the seven hand-written proxy-class functions (AromaticAmines,
    EstrogenMetabolites, NDEA, VinylChloride, UV_Radiation, HeavyMetal, and
    ChlorinatedSolvent). Enzyme scope, role, and rate law come from
    ``GraphEngine.get_edge_flux_reactions`` -- the graph-walk flux
    contract, reading term parameters off baked edge payloads and
    CarcinogenGroup class_level_terms carriers -- and the
    per-term parameters ride on the reaction records. ``derived_scale``
    detox terms scale a named activation term's computed value (the
    ChlorinatedSolvent non-genotoxic clearance: ``v_oxidation x scale``
    with a tissue-weight floor), the detox-direction analogue of how
    repair terms consume the activation total. Dioxin
    (receptor/signaling model) keeps its dedicated function: its model
    structure is not a term-sum proxy.
    """
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    reactions = active_engine.get_edge_flux_reactions(carcinogen_class)
    cfg = _active_proxy_class_params(carcinogen_class, engine)

    def _entry(
        reaction: "FluxReaction", value: float, gm: float, tw: float
    ) -> JsonDict:
        entry: JsonDict = {
            "flux": _round_flux(value),
            "genotype_modifier": gm,
            "tissue_weight": tw,
            "confidence": reaction.confidence,
        }
        if (
            reaction.role != "activation"
            or carcinogen_class not in _PROXY_ACTIVATION_WITHOUT_KINETICS
        ):
            entry["kinetics"] = _PROXY_TERM_KINETICS_LABELS.get(
                reaction.rate_law, "semi_quantitative"
            )
        entry["note"] = str(reaction.params.get("note", ""))
        return entry

    def _evaluate(
        reaction: "FluxReaction", activation_base: float | None = None
    ) -> tuple[float, float, float]:
        term = reaction.params
        if reaction.rate_law == "repair":
            if activation_base is None:
                raise ValueError(
                    f"Repair term {reaction.term_key} evaluated without an activation base"
                )
            return _compute_proxy_repair_term(
                activation_base, term, genotypes, tissue, tissue_weight_source
            )
        if reaction.rate_law == "michaelis_menten":
            return _compute_proxy_mm_term(
                term, genotypes, tissue, S, tissue_weight_source
            )
        if reaction.rate_law == "hill":
            return _compute_proxy_hill_term(
                term, genotypes, tissue, S, tissue_weight_source
            )
        if reaction.rate_law == "saturating":
            return _compute_proxy_saturating_term(
                term, genotypes, tissue, S, tissue_weight_source
            )
        if reaction.rate_law == "derived_scale":
            # Detox term derived from a named activation term's computed
            # value (e.g. ChlorinatedSolvent's non-genotoxic clearance,
            # v_oxidation x scale with a tissue-weight floor): the
            # detox-direction analogue of the repair term's coupling to
            # the activation total. The activation pass runs first, so the
            # base term's value and tissue weight are already recorded.
            base_key = str(term.get("scale_of", ""))
            base = computed_activation.get(base_key)
            if base is None:
                raise ValueError(
                    f"derived_scale term {reaction.term_key} references unknown "
                    f"activation term {base_key!r} for {carcinogen_class}"
                )
            base_value, base_tw = base
            floor = float(term.get("tissue_weight_floor", 0.0))
            return base_value * float(term["scale"]), 1.0, max(base_tw, floor)
        raise ValueError(
            f"Unsupported proxy rate law {reaction.rate_law!r} for "
            f"{carcinogen_class}/{reaction.term_key}"
        )

    activation_enzymes: dict[str, Any] = {}
    total_activation = 0.0
    computed_activation: dict[str, tuple[float, float]] = {}
    for reaction in reactions:
        if reaction.role != "activation":
            continue
        value, gm, tw = _evaluate(reaction)
        total_activation += value
        computed_activation[reaction.term_key] = (value, tw)
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

    def _entry(
        value: float, gm: float, tw: float, role: str, confidence: str
    ) -> JsonDict:
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
            if (
                efficiency_spec is None
                or efficiency_spec.get("term") != reaction.term_key
            ):
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
            raise ValueError(
                f"Derived term {name!r} has no activation/detoxification role"
            )
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
        eff_reaction = next(
            r for r in reactions if r.term_key == efficiency_spec["term"]
        )
        gene = _gene(eff_reaction)
        eff_gm = genotype_modifier(genotypes.get(gene, "NM"), gene)
        eff_tw = get_flux_tissue_weight(gene, tissue, tissue_weight_source)
        efficiency = min(
            1.0,
            (
                float(eff_reaction.params["CLint"])
                / float(efficiency_spec["reference_CLint"])
            )
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


def _class_parameter_metadata(
    carcinogen_class: str,
    engine: "GraphEngine | None" = None,
) -> dict[str, str]:
    """Return class-level parameter metadata for measured or proxy models."""
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    proxy_cfg = active_engine.get_flux_class_config(carcinogen_class)
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
    engine: "GraphEngine | None" = None,
) -> JsonDict | None:
    """Return provenance metadata for a proxy-model term."""
    cfg = _active_proxy_class_params(carcinogen_class, engine)
    sections = (
        ("activation",) if activation_term else ("detoxification", "repair")
    )
    pathways = cfg.get("pathways", {})

    for section in sections:
        term_cfg = pathways.get(section, {}).get(term_name)
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
    engine: "GraphEngine | None" = None,
) -> FluxResultDict:
    """Attach class- and enzyme-level parameter metadata to a flux result."""
    class_meta = _class_parameter_metadata(carcinogen_class, engine)
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
                engine=engine,
            )
            if proxy_meta is None:
                continue
            edata.update(proxy_meta)

    return result


# ── Pathway-specific flux calculators ──────────────────────────────────────


def _measured_constants(params: Mapping[str, Any]) -> tuple[float, float] | None:
    """Km and the numeric Vmax on a term, if both are present."""
    km = params.get("Km_uM")
    vmax = next(
        (value for key, value in params.items()
         if key.startswith("Vmax_") and isinstance(value, (int, float))),
        None,
    )
    if km is None or vmax is None:
        return None
    return float(km), float(vmax)


def _phenotype_aliases(engine: "GraphEngine", gene: str) -> list[set[str]]:
    catalog = engine.G.graph.get("genotypes") or {}
    block = (catalog.get("genes") or {}).get(gene) or {}
    return [
        {str(alias).lower() for alias in phenotype.get("aliases", [])}
        for phenotype in (block.get("phenotypes") or {}).values()
    ]


def _label_matches_diplotype(engine: "GraphEngine", gene: str, label: str, diplotype: Any) -> bool:
    """True when the caller label is this row, or shares a phenotype alias with it."""
    if not label or not diplotype:
        return False
    if str(label).lower() == str(diplotype).lower():
        return True
    key, dip = str(label).lower(), str(diplotype).lower()
    return any(key in aliases and dip in aliases for aliases in _phenotype_aliases(engine, gene))


def _reference_row(engine: "GraphEngine", gene: str, rows: list) -> "FluxReaction":
    """Measured sibling that is not itself a phenotype-alias row."""
    measured = [row for row in rows if _measured_constants(row.params)]
    aliased = set().union(*_phenotype_aliases(engine, gene)) if rows else set()
    plain = [
        row for row in measured
        if str(row.params.get("diplotype") or "").lower() not in aliased
    ]
    return (plain or measured or rows)[0]


def _select_row(engine: "GraphEngine", gene: str, label: str | None, rows: list):
    """Measurement for this label, and a multiplier only when the match has no constants.

    A matched row with Km and Vmax is used as measured. Its activity_fraction
    is not applied again. A matched row with only activity_fraction scales the
    reference measurement. No match leaves the multiplier to the catalog.
    """
    matched = [
        row for row in rows
        if _label_matches_diplotype(engine, gene, label or "", row.params.get("diplotype"))
    ]
    measured = [row for row in matched if _measured_constants(row.params)]
    if measured:
        return measured[0], 1.0
    if matched and matched[0].params.get("activity_fraction") is not None:
        return _reference_row(engine, gene, rows), float(matched[0].params["activity_fraction"])
    return _reference_row(engine, gene, rows), None


def _edge_groups(reactions: list) -> dict[tuple, list]:
    grouped: dict[tuple, list] = {}
    for reaction in reactions:
        key = (
            reaction.params.get("source_node_id"),
            reaction.params.get("target_node_id"),
            reaction.role,
        )
        grouped.setdefault(key, []).append(reaction)
    return grouped


def _evaluate_selected(
    engine: "GraphEngine",
    rows: list,
    genotypes: GenotypeMap,
    tissue: str,
    substrate: float,
    tissue_weight_source: FluxTissueWeightSource,
) -> tuple["FluxReaction", float, float, float] | None:
    gene = str(rows[0].enzyme_id or rows[0].params.get("target_node_id") or rows[0].term_key)
    label = genotypes.get(gene)
    row, fraction = _select_row(engine, gene, label, rows)
    constants = _measured_constants(row.params)
    if constants is None:
        return None
    km, vmax = constants
    modifier = fraction if fraction is not None else genotype_modifier(label or "NM", gene)
    tissue_gene = str(row.params.get("tissue_weight_gene") or gene)
    weight = get_flux_tissue_weight(tissue_gene, tissue, tissue_weight_source)
    scaled = vmax * modifier * weight
    if row.rate_law == "hill":
        value = hill_equation(substrate, scaled, km, float(row.params.get("hill_n", 1.0)))
    else:
        value = michaelis_menten(substrate, scaled, km)
    return row, value, modifier, weight


def _flux_entry(row: "FluxReaction", value: float, modifier: float, weight: float, label: str | None) -> JsonDict:
    entry: JsonDict = {
        "flux": _round_flux(value),
        "genotype": label or "",
        "genotype_modifier": modifier,
        "tissue_weight": weight,
        "confidence": row.confidence,
        "kinetics": row.rate_law,
    }
    if row.params.get("fraction_contribution") is not None:
        entry["fraction_contribution"] = row.params["fraction_contribution"]
    if row.rate_law == "hill" and row.params.get("hill_n") is not None:
        entry["n"] = row.params["hill_n"]
    constants = _measured_constants(row.params)
    if constants is not None:
        km, vmax = constants
        entry["CLint"] = _round_flux((vmax * modifier * weight) / km)
    if row.params.get("reaction"):
        entry["reaction"] = row.params["reaction"]
    return entry


def _compute_aflatoxin_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
    *,
    engine: "GraphEngine | None" = None,
) -> FluxResultDict:
    """AFB1 fluxes. Each edge is one measurement. Modifiers come from that edge."""
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    reactions = active_engine.get_edge_flux_reactions("group_mycotoxins")
    agg = active_engine.get_flux_aggregation("group_mycotoxins")
    activation_enzymes: dict[str, Any] = {}
    detox_enzymes: dict[str, Any] = {}
    total_activation = 0.0
    total_detox = 0.0
    for rows in _edge_groups(reactions).values():
        evaluated = _evaluate_selected(
            active_engine, rows, genotypes, tissue, S, tissue_weight_source
        )
        if evaluated is None:
            continue
        row, value, modifier, weight = evaluated
        gene = str(row.enzyme_id or row.params.get("target_node_id") or row.term_key)
        entry = _flux_entry(row, value, modifier, weight, genotypes.get(gene))
        if row.role == "activation":
            activation_enzymes[row.term_key] = entry
            total_activation += value
        else:
            detox_enzymes[row.term_key] = entry
            total_detox += value
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
    """Aldehyde clearance. Diplotype rows are alternatives, not extra enzymes.

    Copied rows on a second substrate for the same enzyme are not a second
    measurement. Which start owns ethanol oxidation is unchanged here.
    """
    active_engine = engine if engine is not None else _get_flux_contract_engine()
    reactions = active_engine.get_edge_flux_reactions("group_aldehydes")
    agg = active_engine.get_flux_aggregation("group_aldehydes")
    grouped: dict[tuple, tuple] = {}
    for (source, target, role), rows in _edge_groups(reactions).items():
        if (target, role) in grouped:
            continue
        if not any(_measured_constants(row.params) for row in rows):
            continue
        grouped[(target, role)] = (source, rows)

    activation_enzymes: dict[str, Any] = {}
    detox_enzymes: dict[str, Any] = {}
    total_activation = 0.0
    total_detox = 0.0
    ethanol = active_engine.get_flux_metadata()["exposure_defaults_uM"]["ethanol"]
    for (target, _role), (source, rows) in grouped.items():
        substrate = ethanol if source == "Ethanol" else S
        evaluated = _evaluate_selected(
            active_engine, rows, genotypes, tissue, substrate, tissue_weight_source
        )
        if evaluated is None:
            continue
        row, value, modifier, weight = evaluated
        gene = str(row.enzyme_id or target)
        entry = _flux_entry(row, value, modifier, weight, genotypes.get(gene))
        if source == "Ethanol":
            activation_enzymes[gene] = entry
            total_activation += value
        else:
            detox_enzymes[gene] = entry
            total_detox += value
    return {
        "activation_enzymes": activation_enzymes,
        "detox_enzymes": detox_enzymes,
        "total_activation": total_activation,
        "total_detox": total_detox,
        "unit_note": agg["unit_note"],
    }


def _compute_dioxin_flux(
    genotypes: GenotypeMap,
    tissue: str,
    S: float,
    tissue_weight_source: FluxTissueWeightSource,
    *,
    engine: "GraphEngine | None" = None,
) -> FluxResultDict:
    """Compute receptor-mediated dioxin signaling as an induction burden score."""
    cfg = _active_proxy_class_params("group_dioxins", engine)
    signal_p = cfg["signal"]
    signal_strength, _, tissue_signal = _compute_proxy_hill_term(
        signal_p,
        {},
        tissue,
        S,
        tissue_weight_source,
    )

    act1_p = cfg["pathways"]["activation"]["CYP1A1"]
    gm1a1 = _proxy_genotype_modifier(genotypes.get("CYP1A1", "NM"), "CYP1A1")
    tw1a1 = get_flux_tissue_weight("CYP1A1", tissue, tissue_weight_source)
    v1a1 = signal_strength * float(act1_p["induction_factor"]) * gm1a1 * tw1a1

    act2_p = cfg["pathways"]["activation"]["CYP1B1"]
    gm1b1 = _proxy_genotype_modifier(genotypes.get("CYP1B1", "NM"), "CYP1B1")
    tw1b1 = get_flux_tissue_weight("CYP1B1", tissue, tissue_weight_source)
    v1b1 = signal_strength * float(act2_p["induction_factor"]) * gm1b1 * tw1b1

    feedback_p = cfg["pathways"]["detoxification"]["AHRR_feedback"]
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
# contract. Dioxin keeps a fully hand-written function: its
# receptor-mediated model (latent AHR signal gating induced terms) is not
# a term-sum proxy. ChlorinatedSolvent collapsed into the generic loop when
# it gained the ``derived_scale`` rate law (its only structural wrinkle,
# the GSTT1 bioactivation inversion being a role annotation, not code).
_GENERIC_PROXY_FLUX_CLASSES = frozenset(
    {
        "group_aromatic_amines",
        "group_estrogen",
        "NDEA",
        "VinylChloride",
        "group_uv_radiation",
        "group_heavy_metals",
        "group_chlorinated_solvents",
    }
)

# Measured-kinetics classes whose term evaluation and aggregation are fully
# described by the per-class "aggregation" block in kinetic_parameters.json.
_GENERIC_MECHANISTIC_FLUX_CLASSES = frozenset(
    {"group_pahs", "group_nitrosamines", "NDMA", "group_hcas", "Benzene"}
)

# Dedicated proxy functions that now take the engine kwarg so their
# class-config reads are served graph-first when a loaded engine is
# passed (bare JSON fallback otherwise).
_DEDICATED_PROXY_ENGINE_FLUX_CLASSES = frozenset({"group_dioxins"})

# Dedicated functions that still choose which edges belong in the result.
# Row selection inside them follows the graph: a measured diplotype row is
# used as-is, and activity_fraction applies only when that row has no constants.
_DEDICATED_ENGINE_FLUX_CLASSES = frozenset({"group_mycotoxins", "group_aldehydes"})

_DISPATCH: dict[str, FluxCalculator] = {
    "group_pahs": partial(_compute_generic_mechanistic_flux, "group_pahs"),
    "group_mycotoxins": _compute_aflatoxin_flux,
    "group_aldehydes": _compute_aldehyde_flux,
    "group_nitrosamines": partial(
        _compute_generic_mechanistic_flux, "group_nitrosamines"
    ),
    "NDMA": partial(_compute_generic_mechanistic_flux, "NDMA"),
    "NDEA": partial(_compute_generic_proxy_flux, "NDEA"),
    "group_hcas": partial(_compute_generic_mechanistic_flux, "group_hcas"),
    "group_aromatic_amines": partial(
        _compute_generic_proxy_flux, "group_aromatic_amines"
    ),
    "group_estrogen": partial(_compute_generic_proxy_flux, "group_estrogen"),
    "Benzene": partial(_compute_generic_mechanistic_flux, "Benzene"),
    "VinylChloride": partial(_compute_generic_proxy_flux, "VinylChloride"),
    "group_chlorinated_solvents": partial(
        _compute_generic_proxy_flux, "group_chlorinated_solvents"
    ),
    "group_uv_radiation": partial(_compute_generic_proxy_flux, "group_uv_radiation"),
    "group_dioxins": _compute_dioxin_flux,
    "group_heavy_metals": partial(_compute_generic_proxy_flux, "group_heavy_metals"),
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
        engine: Optional loaded graph engine. When given, lifestyle-driven
            induction folds are resolved from ``get_induction_rules()``
            instead of the interaction engine's own file loader.

    Returns:
        :class:`PathwayFluxResult` with activation, detoxification,
        net ratio and risk classification. Returned results also expose
        ``model_kind`` / ``parameter_source`` at the class level, while
        proxy-backed enzyme terms include resolved provenance metadata.
    """
    cls_str = (
        carcinogen_class.value
        if isinstance(carcinogen_class, CarcinogenClass)
        else carcinogen_class
    )
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
        interaction_params=(
            _induction_params_from_rules(engine.get_induction_rules())
            if engine is not None else None
        ),
    )
    if (
        cls_str in _GENERIC_PROXY_FLUX_CLASSES
        or cls_str in _GENERIC_MECHANISTIC_FLUX_CLASSES
        or cls_str in _DEDICATED_ENGINE_FLUX_CLASSES
        or cls_str in _DEDICATED_PROXY_ENGINE_FLUX_CLASSES
    ):
        result = _DISPATCH[cls_str](
            genotypes, tissue, substrate_conc_uM, weight_source, engine=engine
        )
    else:
        result = _DISPATCH[cls_str](genotypes, tissue, substrate_conc_uM, weight_source)
    result = _annotate_flux_result_metadata(cls_str, result, engine=engine)
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
        induction_factors_used={
            gene: _round_flux(factor) for gene, factor in resolved_induction.items()
        },
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
        if r.risk_classification
        in (RiskClassification.ELEVATED, RiskClassification.HIGH)
    ]
    moderate = [
        c
        for c, r in results.items()
        if r.risk_classification == RiskClassification.MODERATE
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
        "group_pahs",
        {"CYP1A1": "NM", "GSTM1": "NM", "CYP1B1": "NM", "GSTP1": "NM", "EPHX1": "NM"},
        "Lung",
        0.1,
        tissue_weight_source=weight_source,
    )
    gstm1_null = compute_pathway_flux(
        "group_pahs",
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
        "group_aldehydes",
        {"ALDH2": "*1/*1", "ALDH1A1": "NM", "ADH1B": "*1/*1"},
        "Liver",
        10.0,
        tissue_weight_source=weight_source,
    )
    aldh2_het = compute_pathway_flux(
        "group_aldehydes",
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
        "group_pahs",
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
  python -m ExposoGraph.flux_engine --genotypes '{"CYP1A1":"NM","GSTM1":"null"}' --tissue Lung --carcinogen group_pahs
  python -m ExposoGraph.flux_engine --full-profile --genotypes '{"ALDH2":"*1/*2"}' --tissue Liver
  python -m ExposoGraph.flux_engine --sensitivity --gene GSTM1 --carcinogen group_pahs --tissue Lung
  python -m ExposoGraph.flux_engine --validate
        """,
    )
    parser.add_argument(
        "--genotypes",
        type=str,
        default="{}",
        help="JSON string of gene-to-phenotype mappings",
    )
    parser.add_argument("--tissue", type=str, default="Liver", help="Target tissue")
    parser.add_argument(
        "--carcinogen", type=str, default=None, help="Carcinogen class to model"
    )
    parser.add_argument(
        "--concentration",
        type=float,
        default=None,
        help="Substrate concentration in uM",
    )
    parser.add_argument(
        "--full-profile", action="store_true", help="Compute all carcinogen classes"
    )
    parser.add_argument(
        "--validate", action="store_true", help="Run built-in validation cases"
    )
    parser.add_argument(
        "--gene", type=str, default=None, help="Gene to vary for sensitivity analysis"
    )
    parser.add_argument("--output-json", action="store_true", help="Output JSON")
    parser.add_argument(
        "--lifestyle",
        type=str,
        default="{}",
        help="JSON lifestyle/co-exposure flags for induction modeling",
    )
    parser.add_argument(
        "--induction-factors",
        type=str,
        default="{}",
        help="JSON per-gene explicit Vmax induction factors",
    )
    parser.add_argument(
        "--qivive",
        action="store_true",
        help="Apply MPPGL/organ-weight QIVIVE scaling to flux magnitudes",
    )
    parser.add_argument(
        "--mppgl",
        type=float,
        default=None,
        help="QIVIVE override: microsomal protein mg/g tissue",
    )
    parser.add_argument(
        "--organ-weight-g",
        type=float,
        default=None,
        help="QIVIVE override: organ weight in grams",
    )
    parser.add_argument(
        "--tissue-weight-source",
        choices=[
            FluxTissueWeightSource.CURATED.value,
            FluxTissueWeightSource.GTEX.value,
        ],
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
        print(
            f"ERROR: Could not parse --induction-factors JSON: {exc}", file=sys.stderr
        )
        return 1

    qivive_context = {}
    if args.mppgl is not None:
        qivive_context["mppgl_mg_per_g"] = args.mppgl
    if args.organ_weight_g is not None:
        qivive_context["organ_weight_g"] = args.organ_weight_g

    if args.full_profile or not args.carcinogen:
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
