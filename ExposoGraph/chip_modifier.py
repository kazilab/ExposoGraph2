"""Clonal hematopoiesis of indeterminate potential (CHIP) somatic modifier layer.

Implements the CHIP module from ExposoGraph 2.0 (CHIP somatic
modifier layer) as a post-hoc somatic adjustment applied on top of germline
flux and interaction results.

The module exposes two integration mechanisms:

1. Aggregate risk modifier (default 1.4x) applied to CHIP-positive individuals
   to reflect the 1.3-1.6x solid-tumor risk elevation observed in published
   CHIP epidemiology (Jaiswal 2017, Bolton 2020, Niroula 2021).
2. Gene-specific pathway modifiers. CHIP driver mutations (DNMT3A, TET2, ASXL1,
   TP53) alter enzyme promoter methylation and damage response, producing
   carcinogen-class-specific multipliers (PAH, nitrosamines, alkylating agents,
   estrogens) when the relevant driver is present.

Germline x CHIP interaction is modelled multiplicatively: the CHIP-adjusted
risk is (germline risk) * (aggregate modifier) * (gene-specific pathway
multiplier). This matches with CHIP-positive individuals in the top germline 
quintile show 2.3x composite scores vs 1.8x for CHIP-negative top-quintile 
individuals.

All parameters are model-derived conservative estimates grounded in published
epidemiology; they are expected to be refined as CHIP-stratified expression
data for carcinogen-metabolizing enzymes become available.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable

# ── Constants ──────────────────────────────────────────────────────────────

# Aggregate CHIP risk modifier: conservative estimate based on
# 1.3-1.6x solid-tumor risk elevations reported for CHIP-positive individuals
# (Jaiswal 2017, Bolton 2020, Niroula 2021).
AGGREGATE_CHIP_MODIFIER: float = 1.4

# VAF threshold for calling CHIP. somatic mutations in
# hematologic malignancy-associated genes at variant allele frequency >=2%".
VAF_THRESHOLD: float = 0.02

# Known CHIP driver genes (Genovese 2014, Jaiswal 2014, Steensma 2015).
CHIP_DRIVER_GENES: frozenset[str] = frozenset(
    {
        "DNMT3A",
        "TET2",
        "ASXL1",
        "TP53",
        "JAK2",
        "SF3B1",
        "SRSF2",
        "PPM1D",
        "IDH1",
        "IDH2",
        "GNB1",
        "GNAS",
        "CBL",
        "BCOR",
    }
)

# Gene-specific pathway modifiers. Keys are CHIP driver genes; values describe
# the carcinogen classes whose risk is elevated when the driver is present,
# with a multiplier and biological justification.
#
# "DNMT3A and TET2 mutations ... may affect
# expression of carcinogen-metabolizing enzymes whose promoters are subject to
# methylation-dependent silencing (CYP1A1, GSTP1, MGMT)"; "TP53 mutations in
# CHIP clones may impair DNA damage response to carcinogen-induced adducts";
# "ASXL1 mutations affect chromatin remodeling and may alter inducibility of
# AhR-responsive genes".
CHIP_GENE_EFFECTS: dict[str, dict[str, Any]] = {
    "DNMT3A": {
        "description": (
            "DNMT3A loss-of-function alters de novo DNA methylation homeostasis. "
            "Expected effect: reduced promoter methylation of AhR-responsive "
            "enzymes (CYP1A1, CYP1B1) with possible over-expression, amplifying "
            "bioactivation of PAH, estrogens, and dioxins."
        ),
        "affected_enzymes": ("CYP1A1", "CYP1B1", "GSTP1"),
        "carcinogen_class_multipliers": {
            "PAH": 1.6,
            "Estrogen": 1.3,
            "Dioxin": 1.3,
            "HCA": 1.2,
        },
        "mechanism": "methylation_dysregulation",
    },
    "TET2": {
        "description": (
            "TET2 loss impairs 5mC->5hmC conversion, leading to promoter "
            "hypermethylation and silencing of GSTP1 and OGG1 (TET-mediated "
            "demethylation of detoxification/repair gene promoters, "
            ". Reduces Phase II conjugation of PAH and "
            "oxidative-damage repair."
        ),
        "affected_enzymes": ("GSTP1", "OGG1", "MGMT"),
        "carcinogen_class_multipliers": {
            "PAH": 1.5,
            "Nitrosamine": 1.5,
            "Aldehyde": 1.4,
            "NDMA": 1.5,
        },
        "mechanism": "repair_impairment",
    },
    "ASXL1": {
        "description": (
            "ASXL1 mutations affect polycomb-group chromatin remodeling and "
            "alter inducibility of AhR-responsive genes. Effect on carcinogen "
            "metabolism is intermediate between DNMT3A and TET2."
        ),
        "affected_enzymes": ("CYP1A1", "CYP1A2", "CYP1B1"),
        "carcinogen_class_multipliers": {
            "PAH": 1.3,
            "Dioxin": 1.5,
            "HCA": 1.3,
        },
        "mechanism": "chromatin_dysregulation",
    },
    "TP53": {
        "description": (
            "TP53 loss-of-function in CHIP clones impairs DNA damage response "
            "to carcinogen-induced adducts. Increases the downstream mutagenic "
            "burden of any DNA-reactive intermediate regardless of specific "
            "enzymatic activation pathway."
        ),
        "affected_enzymes": (),
        "carcinogen_class_multipliers": {
            "PAH": 1.8,
            "Aflatoxin": 1.8,
            "Nitrosamine": 1.8,
            "NDMA": 1.8,
            "AlkylatingAgent": 2.0,
            "HeavyMetal": 1.5,
            "Benzene": 1.6,
        },
        "mechanism": "damage_response_impairment",
    },
    "PPM1D": {
        "description": (
            "PPM1D gain-of-function dephosphorylates p53, further impairing "
            "DNA damage response in CHIP clones."
        ),
        "affected_enzymes": (),
        "carcinogen_class_multipliers": {
            "PAH": 1.4,
            "AlkylatingAgent": 1.6,
            "Nitrosamine": 1.4,
        },
        "mechanism": "damage_response_impairment",
    },
}


# ── Dataclasses ────────────────────────────────────────────────────────────


@dataclass
class CHIPMutation:
    """A single CHIP-associated somatic mutation call."""

    gene: str
    variant_allele_fraction: float
    variant: str | None = None
    consequence: str | None = None

    def is_driver(self) -> bool:
        """Return True when the mutation is in a recognized CHIP driver gene."""
        return self.gene.upper() in CHIP_DRIVER_GENES

    def passes_vaf_threshold(self, threshold: float = VAF_THRESHOLD) -> bool:
        """Return True when the VAF meets the CHIP-calling threshold."""
        return self.variant_allele_fraction >= threshold


@dataclass
class CHIPStatus:
    """CHIP status for a single individual."""

    present: bool
    mutations: list[CHIPMutation] = field(default_factory=list)
    max_vaf: float = 0.0

    @property
    def driver_genes(self) -> list[str]:
        """Return the set of CHIP driver genes carried by this individual."""
        seen: list[str] = []
        for mutation in self.mutations:
            gene = mutation.gene.upper()
            if mutation.is_driver() and mutation.passes_vaf_threshold() and gene not in seen:
                seen.append(gene)
        return seen


@dataclass
class CHIPEffect:
    """Result of applying the CHIP modifier layer to a germline risk profile."""

    chip_positive: bool
    aggregate_modifier: float
    per_class_multipliers: dict[str, float]
    driver_genes: list[str]
    affected_pathways: dict[str, list[str]]
    mechanism_summary: str
    notes: list[str] = field(default_factory=list)


# ── Parsing helpers ────────────────────────────────────────────────────────


def parse_chip_calls(calls: Iterable[dict[str, Any]] | None) -> CHIPStatus:
    """Build a :class:`CHIPStatus` from a list of variant-call records.

    Each call should be a mapping with at minimum ``gene`` and
    ``variant_allele_fraction`` (or ``vaf``). Calls below the VAF threshold or
    not in a recognized driver gene are retained for transparency but do not
    trigger a CHIP-positive classification.
    """
    if not calls:
        return CHIPStatus(present=False, mutations=[], max_vaf=0.0)

    mutations: list[CHIPMutation] = []
    for raw in calls:
        gene = str(raw.get("gene", "")).upper()
        if not gene:
            continue
        vaf_value = raw.get("variant_allele_fraction", raw.get("vaf", 0.0))
        try:
            vaf = float(vaf_value)
        except (TypeError, ValueError):
            vaf = 0.0
        mutations.append(
            CHIPMutation(
                gene=gene,
                variant_allele_fraction=round(vaf, 4),
                variant=str(raw["variant"]) if raw.get("variant") else None,
                consequence=str(raw["consequence"]) if raw.get("consequence") else None,
            )
        )

    qualifying = [m for m in mutations if m.is_driver() and m.passes_vaf_threshold()]
    max_vaf = max((m.variant_allele_fraction for m in qualifying), default=0.0)
    return CHIPStatus(
        present=bool(qualifying),
        mutations=mutations,
        max_vaf=round(max_vaf, 4),
    )


def derive_chip_status(
    mutations: Iterable[dict[str, Any]] | None = None,
    *,
    is_present: bool | None = None,
) -> CHIPStatus:
    """Convenience wrapper: build CHIPStatus from either calls or a binary flag.

    ``is_present`` takes precedence over empty call lists, enabling summary
    datasets that carry only a binary CHIP flag without variant-level detail.
    """
    status = parse_chip_calls(mutations)
    if is_present is not None and not status.mutations:
        status = CHIPStatus(present=bool(is_present), mutations=[], max_vaf=0.0)
    return status


# ── Core modifier computation ──────────────────────────────────────────────


def compute_chip_modifier(
    chip_status: CHIPStatus | None,
    *,
    aggregate_modifier: float = AGGREGATE_CHIP_MODIFIER,
) -> CHIPEffect:
    """Translate a CHIP status into per-carcinogen-class multipliers.

    Multipliers from multiple driver genes combine multiplicatively (consistent
    with the finding that CHIP amplifies, rather than overrides,
    germline susceptibility).
    """
    if chip_status is None or not chip_status.present:
        return CHIPEffect(
            chip_positive=False,
            aggregate_modifier=1.0,
            per_class_multipliers={},
            driver_genes=[],
            affected_pathways={},
            mechanism_summary="CHIP-negative: no somatic modifier applied.",
        )

    per_class: dict[str, float] = {}
    affected_pathways: dict[str, list[str]] = {}
    mechanisms: list[str] = []
    notes: list[str] = []

    for gene in chip_status.driver_genes:
        effect = CHIP_GENE_EFFECTS.get(gene)
        if effect is None:
            notes.append(
                f"CHIP driver {gene} lacks a curated pathway modifier; "
                "only the aggregate modifier is applied."
            )
            continue
        for class_name, multiplier in effect["carcinogen_class_multipliers"].items():
            per_class[class_name] = round(per_class.get(class_name, 1.0) * float(multiplier), 4)
            pathway_list = affected_pathways.setdefault(class_name, [])
            if gene not in pathway_list:
                pathway_list.append(gene)
        mechanisms.append(f"{gene} ({effect['mechanism']})")

    summary_bits: list[str] = []
    if mechanisms:
        summary_bits.append(f"Active CHIP drivers: {', '.join(mechanisms)}.")
    summary_bits.append(
        f"Aggregate modifier: {aggregate_modifier:.2f}x applied to germline risk."
    )
    if per_class:
        top_class = max(per_class.items(), key=lambda item: item[1])
        summary_bits.append(
            f"Largest pathway amplification: {top_class[0]} "
            f"({top_class[1]:.2f}x)."
        )

    return CHIPEffect(
        chip_positive=True,
        aggregate_modifier=round(float(aggregate_modifier), 4),
        per_class_multipliers=per_class,
        driver_genes=list(chip_status.driver_genes),
        affected_pathways=affected_pathways,
        mechanism_summary=" ".join(summary_bits),
        notes=notes,
    )


# ── Applying the modifier to existing risk outputs ─────────────────────────


def apply_chip_to_risk_score(
    risk_score: float,
    effect: CHIPEffect,
    carcinogen_class: str | None = None,
) -> float:
    """Return a CHIP-adjusted risk score for one carcinogen class.

    ``carcinogen_class`` may be ``None`` to apply only the aggregate modifier
    (useful for composite "any cancer" scores).
    """
    if not effect.chip_positive:
        return round(float(risk_score), 4)

    adjusted = float(risk_score) * float(effect.aggregate_modifier)
    if carcinogen_class is not None:
        class_multiplier = effect.per_class_multipliers.get(carcinogen_class, 1.0)
        adjusted *= float(class_multiplier)
    return round(adjusted, 4)


def apply_chip_to_risk_map(
    risks: dict[str, float],
    effect: CHIPEffect,
) -> dict[str, float]:
    """Apply :func:`apply_chip_to_risk_score` across a {class: score} mapping."""
    adjusted: dict[str, float] = {}
    for class_name, score in risks.items():
        adjusted[class_name] = apply_chip_to_risk_score(score, effect, class_name)
    return adjusted


def combine_with_germline_profile(
    germline_risk_scores: dict[str, float],
    effect: CHIPEffect,
    *,
    germline_quintile: int | None = None,
) -> dict[str, Any]:
    """Bundle germline-only and CHIP-adjusted risk scores for downstream use.

    When ``germline_quintile`` is supplied (1-5), the helper also returns a
    ``quintile_interaction_ratio`` that estimates the multiplicative
    germline x CHIP amplification for the top-quintile vs bottom-quintile
    contrast.
    """
    adjusted = apply_chip_to_risk_map(germline_risk_scores, effect)
    aggregate_adjusted_total = sum(adjusted.values())
    germline_total = sum(germline_risk_scores.values())
    fold_change = (
        round(aggregate_adjusted_total / germline_total, 4)
        if germline_total > 0
        else None
    )

    result: dict[str, Any] = {
        "germline_scores": dict(germline_risk_scores),
        "chip_adjusted_scores": adjusted,
        "aggregate_fold_change": fold_change,
        "chip_positive": effect.chip_positive,
        "driver_genes": list(effect.driver_genes),
    }
    if germline_quintile is not None and effect.chip_positive:
        # Approximate the top-vs-bottom quintile amplification reported in the
        # 2.3x top-quintile CHIP+ vs 1.8x top-quintile CHIP-.
        result["quintile_interaction_ratio"] = round(
            2.3 / 1.8 if germline_quintile >= 4 else 1.0,
            4,
        )
    return result


# ── JSON-friendly serialization ────────────────────────────────────────────


def chip_effect_to_dict(effect: CHIPEffect) -> dict[str, Any]:
    """Convert a :class:`CHIPEffect` into a JSON-serializable dictionary."""
    payload = asdict(effect)
    payload["affected_pathways"] = {
        class_name: list(genes) for class_name, genes in effect.affected_pathways.items()
    }
    return payload


def chip_status_to_dict(status: CHIPStatus) -> dict[str, Any]:
    """Convert a :class:`CHIPStatus` into a JSON-serializable dictionary."""
    return {
        "present": status.present,
        "max_vaf": status.max_vaf,
        "driver_genes": status.driver_genes,
        "mutations": [asdict(mutation) for mutation in status.mutations],
    }


def get_chip_gene_effects() -> dict[str, dict[str, Any]]:
    """Return a defensive copy of the CHIP gene-effect table."""
    return deepcopy(CHIP_GENE_EFFECTS)


__all__ = [
    "AGGREGATE_CHIP_MODIFIER",
    "CHIPEffect",
    "CHIPMutation",
    "CHIPStatus",
    "CHIP_DRIVER_GENES",
    "CHIP_GENE_EFFECTS",
    "VAF_THRESHOLD",
    "apply_chip_to_risk_map",
    "apply_chip_to_risk_score",
    "chip_effect_to_dict",
    "chip_status_to_dict",
    "combine_with_germline_profile",
    "compute_chip_modifier",
    "derive_chip_status",
    "get_chip_gene_effects",
    "parse_chip_calls",
]
