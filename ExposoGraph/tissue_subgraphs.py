"""Tissue-specific subgraph extraction and analysis.

Provides functions to extract tissue-specific subgraphs from the knowledge graph,
analyze tissue-specific metabolism, and identify cancer site-specific pathways.
"""

from __future__ import annotations

import argparse
import copy
import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from .engine import GraphEngine
from .graph_analysis import MetabolismChain, metabolism_chain
from .models import KnowledgeGraph, Node, NodeType


class TissueType(str, Enum):
    """Major tissues/organs relevant to metal carcinogenesis."""

    LIVER = "liver"
    LUNG = "lung"
    KIDNEY = "kidney"
    BLADDER = "bladder"
    SKIN = "skin"
    PROSTATE = "prostate"
    BREAST = "breast"
    COLON = "colon"
    ESOPHAGUS = "esophagus"
    STOMACH = "stomach"
    PANCREAS = "pancreas"
    BRAIN = "brain"
    BONE_MARROW = "bone_marrow"
    LYMPHOCYTES = "lymphocytes"
    NASAL_MUCOSA = "nasal_mucosa"
    INTESTINE = "intestine"
    PLEURA = "pleura"
    PLACENTA = "placenta"
    UTERUS = "uterus"


class CancerSite(str, Enum):
    """Cancer sites with established metal carcinogen associations."""

    LUNG_ADENOCARCINOMA = "lung_adenocarcinoma"
    LUNG_SQUAMOUS = "lung_squamous_cell_carcinoma"
    BLADDER_UROTHELIAL = "bladder_urothelial_carcinoma"
    SKIN_SQUAMOUS = "skin_squamous_cell_carcinoma"
    PROSTATE_ADENOCARCINOMA = "prostate_adenocarcinoma"
    BREAST_CARCINOMA = "breast_carcinoma"
    KIDNEY_RENAL_CELL = "kidney_renal_cell_carcinoma"
    LIVER_HEPATOCELLULAR = "liver_hepatocellular_carcinoma"
    COLON_ADENOCARCINOMA = "colon_adenocarcinoma"
    ESOPHAGEAL_SQUAMOUS = "esophageal_squamous_cell_carcinoma"
    NASAL_SINONASAL = "nasal_sinonasal_carcinoma"
    PLEURAL_MESOTHELIOMA = "pleural_mesothelioma"  # for beryllium


# Tissue expression patterns for metal metabolism genes
# Based on GTEx and tissue-specific expression databases
TISSUE_EXPRESSION_PATTERNS: dict[str, dict[str, Any]] = {
    # Phase I enzymes
    "CYP1A1": {
        "high": [TissueType.LUNG, TissueType.SKIN, TissueType.LYMPHOCYTES, TissueType.PLACENTA],
        "medium": [TissueType.LIVER],
        "low": [TissueType.KIDNEY, TissueType.COLON],
        "inducible": True,  # By AhR ligands (TCDD, PAHs)
    },
    "CYP1A2": {
        "high": [TissueType.LIVER],
        "medium": [],
        "low": [TissueType.LUNG],
        "inducible": True,
    },
    "CYP1B1": {
        "high": [TissueType.BREAST, TissueType.PROSTATE, TissueType.UTERUS, TissueType.LUNG],
        "medium": [TissueType.KIDNEY],
        "low": [TissueType.LIVER],
        "inducible": False,
    },
    "CYP2A6": {
        "high": [TissueType.LIVER, TissueType.NASAL_MUCOSA],
        "medium": [],
        "low": [TissueType.LUNG],
        "inducible": False,
    },
    "CYP2A13": {
        "high": [TissueType.LUNG, TissueType.NASAL_MUCOSA],
        "medium": [],
        "low": [],
        "inducible": False,
    },
    "CYP2E1": {
        "high": [TissueType.LIVER, TissueType.LUNG],
        "medium": [TissueType.KIDNEY, TissueType.BRAIN],
        "low": [TissueType.COLON],
        "inducible": True,  # By ethanol, acetone, isoniazid
    },
    "CYP3A4": {
        "high": [TissueType.LIVER, TissueType.INTESTINE],
        "medium": [TissueType.KIDNEY, TissueType.LUNG],
        "low": [TissueType.BRAIN],
        "inducible": True,
    },

    # Phase II enzymes
    "GSTM1": {
        "high": [TissueType.LIVER, TissueType.LUNG, TissueType.BRAIN],
        "medium": [TissueType.KIDNEY],
        "low": [],
        "inducible": True,
    },
    "GSTT1": {
        "high": [TissueType.LIVER, TissueType.KIDNEY],
        "medium": [TissueType.BONE_MARROW],
        "low": [],
        "inducible": False,
    },
    "GSTP1": {
        "high": [TissueType.LUNG, TissueType.BRAIN, TissueType.PLACENTA],
        "medium": [TissueType.LIVER, TissueType.KIDNEY],
        "low": [],
        "inducible": True,
    },
    "NQO1": {
        "high": [TissueType.LIVER, TissueType.LUNG, TissueType.COLON],
        "medium": [TissueType.KIDNEY],
        "low": [],
        "inducible": True,  # By Nrf2 activators
    },
    "NAT2": {
        "high": [TissueType.LIVER, TissueType.INTESTINE],
        "medium": [],
        "low": [TissueType.BLADDER, TissueType.COLON],
        "inducible": False,
    },
    "EPHX1": {
        "high": [TissueType.LIVER, TissueType.LUNG],
        "medium": [TissueType.KIDNEY],
        "low": [],
        "inducible": False,
    },
    "UGT1A1": {
        "high": [TissueType.LIVER, TissueType.INTESTINE],
        "medium": [TissueType.KIDNEY],
        "low": [],
        "inducible": True,
    },

    # Arsenic metabolism
    "AS3MT": {
        "high": [TissueType.LIVER, TissueType.KIDNEY],
        "medium": [TissueType.LUNG, TissueType.SKIN],
        "low": [TissueType.BRAIN],
        "inducible": False,
    },
    "GSTO1": {
        "high": [TissueType.LIVER, TissueType.KIDNEY],
        "medium": [TissueType.LUNG],
        "low": [TissueType.BRAIN],
        "inducible": False,
    },
    "GSTO2": {
        "high": [TissueType.LIVER, TissueType.KIDNEY],
        "medium": [],
        "low": [],
        "inducible": False,
    },

    # Cadmium metabolism
    "MT1A": {
        "high": [TissueType.LIVER, TissueType.KIDNEY],
        "medium": [TissueType.PANCREAS, TissueType.INTESTINE],
        "low": [TissueType.BRAIN],
        "inducible": True,  # By metal exposure
    },
    "MT2A": {
        "high": [TissueType.LIVER, TissueType.KIDNEY],
        "medium": [TissueType.PANCREAS],
        "low": [TissueType.BRAIN],
        "inducible": True,
    },

    # DNA repair - ubiquitous but with tissue differences
    "XRCC1": {
        "high": [TissueType.LIVER, TissueType.LUNG, TissueType.KIDNEY],
        "medium": [],
        "low": [],
        "ubiquitous": True,
    },
    "OGG1": {
        "high": [TissueType.LIVER, TissueType.LUNG],
        "medium": [TissueType.KIDNEY],
        "low": [],
        "ubiquitous": True,
        "mitochondrial": True,
    },
    "XPC": {
        "high": [TissueType.SKIN, TissueType.LIVER],
        "medium": [TissueType.LUNG],
        "low": [],
        "ubiquitous": True,
    },
    "ERCC2": {
        "high": [TissueType.LIVER, TissueType.LUNG],
        "medium": [],
        "low": [],
        "ubiquitous": True,
    },
    "MGMT": {
        "high": [TissueType.LIVER, TissueType.COLON, TissueType.BRAIN],
        "medium": [TissueType.LUNG],
        "low": [],
        "ubiquitous": True,
    },

    # Immune/HLA
    "HLA-DPB1": {
        "high": [TissueType.LYMPHOCYTES, TissueType.LUNG],
        "medium": [],
        "low": [],
        "immune_restricted": True,
    },
}

# ── GTEx v8 quantitative expression data ──────────────────────────────────

_GTEX_DATA_FILE = Path(__file__).parent / "data" / "tissue_expression_data.json"
_GTEX_CACHE: dict[str, Any] | None = None

# Mapping from TissueType enum values to GTEx tissue name strings
_TISSUE_TO_GTEX: dict[TissueType, str] = {
    TissueType.LIVER: "Liver",
    TissueType.LUNG: "Lung",
    TissueType.PROSTATE: "Prostate",
    TissueType.BLADDER: "Bladder",
    TissueType.COLON: "Colon",
    TissueType.BREAST: "Breast",
    TissueType.KIDNEY: "Kidney",
    TissueType.ESOPHAGUS: "Esophagus",
}

_PHASE1_GENES = [
    "CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2A13", "CYP2B6",
    "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP2F1", "CYP3A4",
    "CYP3A5", "CYP17A1", "CYP19A1", "EPHX1", "AHR", "ARNT", "AHRR",
]
_PHASE2_CONJUGATION_GENES = [
    "GSTM1", "GSTP1", "GSTT1", "GSTO1", "GSTO2",
    "NAT1", "NAT2", "SULT1A1",
    "UGT1A1", "UGT2B7", "UGT2B15", "UGT2B17",
    "COMT", "AS3MT",
]
_PHASE2_OTHER_GENES = [
    "NQO1", "ALDH2", "ADH1B", "ADH1C", "ADH5",
    "AKR1C2", "AKR1C3", "SRD5A1", "SRD5A2", "HSD3B2",
    "CCBL1", "MT1A", "MT2A",
]
_TRANSPORT_GENES = ["ABCB1", "ABCC2", "ABCG2", "SLC11A1"]
_DNA_REPAIR_GENES = [
    "XRCC1",
    "MGMT",
    "MLH1",
    "MSH2",
    "OGG1",
    "ERCC2",
    "XPC",
    "FANCD2",
    "PARP1",
]
_IMMUNE_GENES = ["HLA_DPB1"]
_TISSUE_REPORT_CATEGORIES: dict[str, list[str]] = {
    "Phase I (Bioactivation/Oxidation)": _PHASE1_GENES,
    "Phase II (Conjugation)": _PHASE2_CONJUGATION_GENES,
    "Phase II (Reduction/Other)": _PHASE2_OTHER_GENES,
    "Transporters": _TRANSPORT_GENES,
    "DNA Repair": _DNA_REPAIR_GENES,
    "Immune": _IMMUNE_GENES,
}
_TISSUE_CARCINOGEN_CLASS_RANKINGS: dict[str, list[tuple[str, str]]] = {
    "Liver": [
        (
            "PAHs (Polycyclic Aromatic Hydrocarbons)",
            "CYP1A1/1A2/3A4/EPHX1 all highly expressed; GSTM1 detox capacity present",
        ),
        (
            "Nitrosamines / N-Nitroso Compounds",
            "CYP2A6 (2nd-highest hepatic CYP) and CYP2E1 (dominant activator) both liver-dominant",
        ),
        (
            "Aflatoxins",
            "CYP3A4 + EPHX1 activation pathway fully expressed; UGT1A1/SULT1A1 detox present",
        ),
        (
            "Alcohols and Aldehydes",
            "ADH1B, ADH1C (highest overall expression), ALDH2 all liver-dominant",
        ),
        (
            "Aromatic Amines / HCAs",
            "NAT2 liver-dominant; SULT1A1 high; GSTM1 moderate; hepatic arylamine activation complete",
        ),
        (
            "Chlorinated Solvents",
            "CYP2E1 is the dominant activator for TCE, vinyl chloride; liver is primary exposure organ",
        ),
        (
            "Estrogens / Endocrine Disruptors",
            "CYP3A4 (estrogen hydroxylation), UGT1A1 conjugation; moderate CYP19A1",
        ),
        (
            "Heavy Metals (Arsenic, Cadmium)",
            "AS3MT arsenic methylation liver-dominant; MT1A/MT2A sequestration present",
        ),
    ],
    "Lung": [
        (
            "PAHs (Polycyclic Aromatic Hydrocarbons)",
            "CYP1A1 expressed (lung-inducible); AHR signaling (AHRR lung-enriched); GSTP1 very high detox",
        ),
        (
            "Tobacco-Specific Nitrosamines (TSNAs / NNK)",
            "CYP2A13 lung-specific CYP (NNK primary activator); CYP2E1 present",
        ),
        (
            "Dioxins / Furans (AhR Ligands)",
            "AHRR specifically enriched in lung; AHR/ARNT expressed; CYP1A1 inducible",
        ),
        (
            "Aldehydes (Acrolein, Formaldehyde)",
            "ADH1B and ALDH2 moderate lung expression; direct mucosal exposure",
        ),
        (
            "Heavy Metals (Cadmium, Hexavalent Chromium, Arsenic)",
            "AS3MT present; MT2A expressed; lung is primary inhalation route",
        ),
        (
            "HCAs (Heterocyclic Amines)",
            "NAT1 present; NQO1 moderate; lung exposure via secondhand smoke",
        ),
        (
            "Asbestos / Mineral Fibers",
            "Oxidative stress pathway; NQO1 present; PARP1 DNA repair active",
        ),
    ],
    "Prostate": [
        (
            "Androgens and Androgen-Activating Compounds",
            "SRD5A2 (DHT synthesis), SRD5A1, HSD3B2, CYP17A1 expressed; androgen axis complete",
        ),
        (
            "Estrogens / Estrogenic Compounds",
            "CYP1B1 is dominant Phase I enzyme in prostate; COMT catechol-estrogen detox present",
        ),
        (
            "PAHs (Polycyclic Aromatic Hydrocarbons)",
            "CYP1B1 high (PAH activator); GSTM1 expressed; epidemiological links to grilled-meat exposure",
        ),
        (
            "Cadmium and Heavy Metals",
            "Cadmium epidemiologically linked to prostate cancer; CCBL1, MT2A expressed",
        ),
        (
            "HCAs / Dietary Carcinogens",
            "NAT1 present; AKR1C3 expressed; well-done red meat associations",
        ),
        (
            "Pesticides / Organochlorines",
            "Lipophilic accumulation in prostate stroma; CYP3A5 expressed",
        ),
    ],
    "Bladder": [
        (
            "Aromatic Amines (Benzidine, 4-ABP, 2-Naphthylamine)",
            "CYP1A1 highly expressed in bladder (top IARC Group 1 mechanism); GSTP1 very high; NAT2 low -> reduced detox",
        ),
        (
            "PAHs (Polycyclic Aromatic Hydrocarbons)",
            "CYP1A1 high (bladder nTPM ~= liver levels); GSTP1 abundant detox",
        ),
        (
            "Chlorinated Solvents (TCE, Perchloroethylene)",
            "CYP2E1 present; GSTP1 detox; occupational/environmental exposure route",
        ),
        (
            "Nitrosamines",
            "Urinary concentration effect; CYP2E1 present in bladder epithelium",
        ),
        (
            "Acrolein (Cyclophosphamide metabolite)",
            "Direct urothelial toxin; GSTP1 conjugation; ALDH2 moderate",
        ),
        (
            "Arsenic",
            "Drinking water exposure; AS3MT low -> limited methylation detox",
        ),
    ],
    "Colon": [
        (
            "HCAs (Heterocyclic Amines from cooked meat)",
            "NAT1/NAT2 both expressed; NQO1 high; strong epidemiological association",
        ),
        (
            "Dietary N-Nitroso Compounds (DNOCs)",
            "CYP2E1 expressed; direct colonic exposure from processed meat",
        ),
        (
            "Secondary Bile Acids",
            "UGT2B7 present; AKR1C3 expressed; GSTP1 very high detox",
        ),
        (
            "PAHs (from dietary char/smoke)",
            "CYP1B1 and GSTM1 expressed; dietary PAH delivery to mucosa",
        ),
        (
            "Fiber-fermentation Metabolites / Butyrate",
            "NQO1 very high in colon; protective and carcinogenic metabolic balance",
        ),
        (
            "Alcohol (ethanol metabolites)",
            "ADH1B/ADH1C/ALDH2 all expressed; acetaldehyde mucosal exposure",
        ),
    ],
    "Breast": [
        (
            "Estrogens / Estrogenic Endocrine Disruptors",
            "CYP19A1 (aromatase) present; CYP1B1 high (4-OH-E2 formation); COMT catechol detox",
        ),
        (
            "PAHs (Polycyclic Aromatic Hydrocarbons)",
            "CYP1B1 dominant activator; NQO1 very high protective; GSTM1 expressed",
        ),
        (
            "Aldehydes (Acetaldehyde from alcohol)",
            "ADH1B highest breast expression of all tissues after liver; ALDH2 moderate",
        ),
        (
            "Organochlorine Pesticides / PCBs",
            "Lipophilic accumulation in breast adipose; CYP1B1 activation; AhR-mediated",
        ),
        (
            "HCAs",
            "AKR1C3 breast-enriched; NAT1 expressed; grilled-meat dietary exposure",
        ),
        (
            "Radiation / Oxidative DNA Damage",
            "BRCA pathway context; PARP1, XRCC1, OGG1 all expressed; NQO1 very high antioxidant",
        ),
    ],
    "Kidney": [
        (
            "Heavy Metals (Cadmium, Mercury, Lead, Arsenic)",
            "CCBL1 kidney-enriched (cysteine-metal conjugate processing); MT1A/MT2A expressed; primary excretory organ",
        ),
        (
            "Chlorinated Solvents (Trichloroethylene / TCE)",
            "GSTP1 extremely high (renal GSH-conjugate processing); CCBL1 renal-specific beta-lyase activation",
        ),
        (
            "Aflatoxins",
            "Renal proximal tubule is secondary target; EPHX1 expressed; GSTP1 detox",
        ),
        (
            "Aristolochic Acid",
            "NQO1 high (reductive activation); direct tubular cell exposure; endemic nephropathy",
        ),
        (
            "PAHs",
            "GSTP1 very high; CYP1B1 expressed; occupational/environmental exposure",
        ),
        (
            "Nitrosamines",
            "CYP2E1 expressed; UGT2B7 present; renal tubule concentration effect",
        ),
        (
            "Ochratoxin A (mycotoxin)",
            "Direct proximal tubule toxin; ABCC2/ABCG2 transporter expression determines accumulation",
        ),
    ],
    "Esophagus": [
        (
            "Acetaldehyde / Alcohol",
            "ADH1B/ADH1C/ALDH2 expressed; direct mucosal contact; ALDH2 polymorphism major risk factor",
        ),
        (
            "PAHs (from tobacco smoke, smoked foods)",
            "GSTP1 extremely high (highest of all tissues); CYP1A1 expressed; direct mucosal exposure",
        ),
        (
            "Nitrosamines (tobacco-specific)",
            "CYP2E1 present; direct contact from tobacco smoke and pickled foods",
        ),
        (
            "Thermal Injury / Hot Beverages",
            "Physical carcinogenesis + chemical co-exposure; NQO1/COMT expressed",
        ),
        (
            "Aflatoxins (dietary)",
            "EPHX1 very high in esophagus; direct dietary delivery; GSTM1 detox",
        ),
        (
            "Aromatic Amines",
            "GSTP1 very high provides detox capacity; CYP1A1 activation possible",
        ),
    ],
}


def _load_gtex_data() -> dict[str, Any]:
    """Lazy-load and cache the GTEx tissue expression JSON data."""
    global _GTEX_CACHE
    if _GTEX_CACHE is None:
        if not _GTEX_DATA_FILE.exists():
            raise FileNotFoundError(
                f"GTEx expression data not found at {_GTEX_DATA_FILE}. "
                "Ensure tissue_expression_data.json is in the data/ directory."
            )
        with open(_GTEX_DATA_FILE, "r") as fh:
            _GTEX_CACHE = json.load(fh)
    return _GTEX_CACHE


def get_available_gtex_tissues() -> list[str]:
    """Return the list of tissue names with GTEx expression data."""
    data = _load_gtex_data()
    return list(data["metadata"]["tissues"])


def get_available_gtex_genes() -> list[str]:
    """Return sorted list of gene symbols in the GTEx expression dataset."""
    data = _load_gtex_data()
    return sorted(data["expression"].keys())


def get_available_tissues() -> list[str]:
    """Compatibility alias returning all GTEx-backed tissue names."""
    return get_available_gtex_tissues()


def get_available_genes() -> list[str]:
    """Compatibility alias returning all GTEx-backed gene symbols."""
    return get_available_gtex_genes()


def get_tissue_weights(tissue: TissueType | str) -> dict[str, float]:
    """Return GTEx normalized expression weights for all genes in a tissue.

    Weights are normalized to each gene's maximum-expressing tissue (max = 1.0).

    Args:
        tissue: TissueType enum value or GTEx tissue name string.

    Returns:
        {gene_symbol: weight_float} dictionary.
    """
    gtex_name = _resolve_gtex_tissue(tissue)
    data = _load_gtex_data()
    return {gene: vals[gtex_name] for gene, vals in data["weights"].items()}


def get_tissue_expression(tissue: TissueType | str) -> dict[str, float]:
    """Return raw nTPM expression values for all genes in a tissue.

    Args:
        tissue: TissueType enum value or GTEx tissue name string.

    Returns:
        {gene_symbol: nTPM_float} dictionary.
    """
    gtex_name = _resolve_gtex_tissue(tissue)
    data = _load_gtex_data()
    return {gene: vals[gtex_name] for gene, vals in data["expression"].items()}


def weight_activity_scores(
    activity_scores: dict[str, float],
    tissue: TissueType | str,
) -> dict[str, float]:
    """Multiply per-gene activity scores by tissue-specific expression weights.

    Genes absent from the GTEx dataset are returned unchanged (weight 1.0).

    Args:
        activity_scores: {gene_symbol: score} mapping.
        tissue: TissueType enum value or GTEx tissue name string.

    Returns:
        {gene_symbol: weighted_score} dictionary.
    """
    weights = get_tissue_weights(tissue)
    return {
        gene: round(score * weights.get(gene, 1.0), 6)
        for gene, score in activity_scores.items()
    }


def generate_tissue_report(
    tissue: TissueType | str,
    threshold: float = 0.05,
) -> str:
    """Generate a human-readable summary of tissue-specific enzyme expression.

    Args:
        tissue: TissueType enum value or GTEx tissue name string.
        threshold: Minimum weight to include a gene (default 0.05).

    Returns:
        Multi-line report string.
    """
    gtex_name = _resolve_gtex_tissue(tissue)
    weights = get_tissue_weights(tissue)
    expression = get_tissue_expression(tissue)

    included = {g: (weights[g], expression[g]) for g in weights if weights[g] >= threshold}
    excluded = {g: (weights[g], expression[g]) for g in weights if weights[g] < threshold}

    sorted_included = sorted(included.items(), key=lambda x: x[1][0], reverse=True)

    lines = [
        "=" * 70,
        f"  ExposoGraph Tissue Report: {gtex_name.upper()}",
        f"  Threshold: weight >= {threshold} (>={threshold*100:.0f}% of maximum-expressing tissue)",
        "  Data: GTEx v8 via Human Protein Atlas v25 detail table | Units: nTPM weights",
        "=" * 70,
        "",
        f"INCLUDED ENZYMES ({len(included)} of {len(weights)} total):",
        "-" * 70,
    ]

    for category_name, category_genes in _TISSUE_REPORT_CATEGORIES.items():
        category_rows = [(g, w, e) for g, (w, e) in sorted_included if g in category_genes]
        if not category_rows:
            continue
        lines.append("")
        lines.append(f"  [{category_name}]")
        lines.append(f"  {'Gene':<12} {'Weight':>8}  {'nTPM':>10}  Bar")
        lines.append(f"  {'-'*12} {'-'*8}  {'-'*10}  {'-'*20}")
        for gene, w, e in category_rows:
            bar_len = max(1, int(w * 20))
            bar = "#" * bar_len
            lines.append(f"  {gene:<12} {w:>8.4f}  {e:>10.1f}  {bar}")

    lines.append("")
    lines.append(f"EXCLUDED ENZYMES ({len(excluded)}) [weight < {threshold}]:")
    lines.append("-" * 70)
    exc_str = ", ".join(
        f"{g}({weights[g]:.3f})"
        for g in sorted(excluded.keys())
        if weights[g] > 0
    )
    zero_str = ", ".join(g for g in sorted(excluded.keys()) if weights[g] == 0.0)
    if exc_str:
        lines.append(f"  Low expression:  {exc_str}")
    if zero_str:
        lines.append(f"  Zero expression: {zero_str}")

    lines.append("")
    lines.append("CARCINOGEN CLASS RELEVANCE:")
    lines.append("-" * 70)
    for rank, (class_name, rationale) in enumerate(get_top_carcinogen_classes_for_tissue(tissue), start=1):
        lines.append(f"  {rank}. {class_name}")
        lines.append(f"     {rationale}")

    lines.append("")
    lines.append("=" * 70)

    return "\n".join(lines)


def _resolve_gtex_tissue(tissue: TissueType | str) -> str:
    """Resolve a TissueType or string to a GTEx tissue name."""
    if isinstance(tissue, TissueType):
        gtex_name = _TISSUE_TO_GTEX.get(tissue)
        if gtex_name is None:
            raise ValueError(
                f"No GTEx data available for {tissue.value}. "
                f"Available: {list(_TISSUE_TO_GTEX.keys())}"
            )
        return gtex_name

    # String input — validate against available tissues
    available = get_available_gtex_tissues()
    if tissue not in available:
        raise ValueError(f"Unknown tissue '{tissue}'. Available: {available}")
    return tissue


def get_top_carcinogen_classes_for_tissue(tissue: TissueType | str) -> list[tuple[str, str]]:
    """Return ranked carcinogen classes for a tissue using the source module logic."""
    gtex_name = _resolve_gtex_tissue(tissue)
    return list(_TISSUE_CARCINOGEN_CLASS_RANKINGS.get(gtex_name, []))


def filter_graph_by_tissue(
    graph_data: dict[str, Any],
    tissue: TissueType | str,
    threshold: float = 0.05,
) -> dict[str, Any]:
    """Filter a raw graph-data dictionary to the enzymes expressed in *tissue*.

    This is a source-compatible JSON helper for users migrating from the
    standalone ``03_tissue_views`` module. Enzymes below the threshold are
    removed while non-enzyme nodes remain available for pathway context.
    """
    gtex_name = _resolve_gtex_tissue(tissue)
    weights = get_tissue_weights(gtex_name)
    filtered = copy.deepcopy(graph_data)

    included_ids: set[str] = set()
    excluded_ids: set[str] = set()
    filtered_nodes: list[dict[str, Any]] = []

    for node in filtered.get("nodes", []):
        node_id = node.get("id") or node.get("gene") or ""
        node_type_value = getattr(node.get("type", ""), "value", node.get("type", ""))
        node_type = str(node_type_value)
        node_type_norm = node_type.lower()

        if node_type_norm not in {"enzyme", "gene"}:
            filtered_nodes.append(node)
            if node_id:
                included_ids.add(node_id)
            continue

        weight = weights.get(node_id)
        if weight is None:
            node["tissue_weight"] = 1.0
            node["tissue_weight_source"] = "not_in_dataset"
            filtered_nodes.append(node)
            if node_id:
                included_ids.add(node_id)
            continue

        if weight >= threshold:
            node["tissue_weight"] = round(float(weight), 6)
            node["tissue_weight_source"] = "gtex_v8_hpa"
            filtered_nodes.append(node)
            if node_id:
                included_ids.add(node_id)
        elif node_id:
            excluded_ids.add(node_id)

    filtered["nodes"] = filtered_nodes

    filtered_edges: list[dict[str, Any]] = []
    for edge in filtered.get("edges", []):
        source = edge.get("source", "")
        target = edge.get("target", "")
        if source in included_ids and target in included_ids:
            filtered_edges.append(edge)
    filtered["edges"] = filtered_edges

    filtered.setdefault("tissue_metadata", {})
    filtered["tissue_metadata"].update(
        {
            "tissue": gtex_name,
            "threshold": threshold,
            "included_enzymes": sorted(node_id for node_id in included_ids if node_id in weights),
            "excluded_enzymes": sorted(excluded_ids),
            "data_source": "GTEx v8 via Human Protein Atlas v25 detail table (proteinatlas.org)",
            "units": "nTPM weights normalized to tissue maximum",
        }
    )
    return filtered


#   Sensitivity analysis across thresholds (0.10, 0.25, 0.50)
#    tissue edge retention ranges shift by 8-15 percentage points".
DEFAULT_THRESHOLD_SWEEP: tuple[float, ...] = (0.10, 0.25, 0.50)


@dataclass
class TissueThresholdRetention:
    """Edge/node retention statistics for one tissue at one threshold."""

    tissue: str
    threshold: float
    nodes_retained: int
    nodes_total: int
    edges_retained: int
    edges_total: int
    enzyme_coverage: float
    edge_retention: float


@dataclass
class ThresholdSweepResult:
    """Threshold sweep across tissues for a single graph."""

    thresholds: list[float]
    entries: list[TissueThresholdRetention]
    edge_retention_matrix: dict[str, dict[float, float]]


def tissue_threshold_sweep(
    graph_data: dict[str, Any],
    tissues: list[TissueType | str] | None = None,
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLD_SWEEP,
) -> ThresholdSweepResult:
    """Sweep GTEx-expression thresholds across tissues and return retention stats.

    for each (tissue, threshold) pair, report the fraction of enzyme nodes and
    edges retained once ``filter_graph_by_tissue`` is applied. This lets
    downstream callers reproduce 0.10 / 0.25 / 0.50 threshold
    comparison and identify threshold-sensitive tissues such as mammary.
    """
    if tissues is None:
        tissues = get_available_gtex_tissues()

    base_node_count = len(graph_data.get("nodes", []))
    base_edge_count = len(graph_data.get("edges", []))

    entries: list[TissueThresholdRetention] = []
    matrix: dict[str, dict[float, float]] = {}

    for tissue in tissues:
        try:
            tissue_label = _resolve_gtex_tissue(tissue)
        except Exception:
            tissue_label = str(tissue)
        matrix.setdefault(tissue_label, {})
        for threshold in thresholds:
            filtered = filter_graph_by_tissue(graph_data, tissue_label, threshold=threshold)
            nodes_retained = len(filtered.get("nodes", []))
            edges_retained = len(filtered.get("edges", []))
            included_enzymes = (
                filtered.get("tissue_metadata", {}).get("included_enzymes") or []
            )
            excluded_enzymes = (
                filtered.get("tissue_metadata", {}).get("excluded_enzymes") or []
            )
            enzyme_total = len(included_enzymes) + len(excluded_enzymes)
            enzyme_coverage = (
                len(included_enzymes) / enzyme_total if enzyme_total > 0 else 1.0
            )
            edge_retention = (
                edges_retained / base_edge_count if base_edge_count > 0 else 1.0
            )
            entries.append(
                TissueThresholdRetention(
                    tissue=tissue_label,
                    threshold=threshold,
                    nodes_retained=nodes_retained,
                    nodes_total=base_node_count,
                    edges_retained=edges_retained,
                    edges_total=base_edge_count,
                    enzyme_coverage=round(enzyme_coverage, 4),
                    edge_retention=round(edge_retention, 4),
                )
            )
            matrix[tissue_label][threshold] = round(edge_retention, 4)

    return ThresholdSweepResult(
        thresholds=list(thresholds),
        entries=entries,
        edge_retention_matrix=matrix,
    )


def annotate_graph_with_tissue_weights(
    graph: KnowledgeGraph | GraphEngine | dict[str, Any],
) -> KnowledgeGraph:
    """Attach GTEx ``tissue_weights`` annotations to enzyme and gene nodes.

    The result is a package-native :class:`KnowledgeGraph` that exports cleanly
    to the D3 reference map, so tissue-aware filters can operate on the bundled
    graph-data.js asset.
    """
    if isinstance(graph, GraphEngine):
        base_graph = graph.to_knowledge_graph()
    elif isinstance(graph, KnowledgeGraph):
        base_graph = graph
    else:
        base_graph = KnowledgeGraph(**graph)

    tissue_names = get_available_tissues()
    weights_by_tissue = {tissue_name: get_tissue_weights(tissue_name) for tissue_name in tissue_names}

    enriched_nodes: list[Node] = []
    for node in base_graph.nodes:
        if node.type not in {NodeType.ENZYME, NodeType.GENE}:
            enriched_nodes.append(node.model_copy(deep=True))
            continue

        tissue_weights: dict[str, float] = {}
        for tissue_name, tissue_map in weights_by_tissue.items():
            if node.id in tissue_map:
                tissue_weights[tissue_name] = round(float(tissue_map[node.id]), 6)

        if tissue_weights:
            merged_weights = dict(node.tissue_weights or {})
            merged_weights.update(tissue_weights)
            enriched_nodes.append(node.model_copy(update={"tissue_weights": merged_weights}))
        else:
            enriched_nodes.append(node.model_copy(deep=True))

    return KnowledgeGraph(
        nodes=enriched_nodes,
        edges=[edge.model_copy(deep=True) for edge in base_graph.edges],
    )


# Cancer site-specific metal associations
CANCER_SITE_METAL_PROFILES: dict[CancerSite, dict[str, Any]] = {
    CancerSite.LUNG_ADENOCARCINOMA: {
        "primary_metals": ["Arsenic", "Chromium(VI)", "Nickel"],
        "key_genes": ["CYP1A1", "CYP2A13", "GSTM1", "GSTT1", "OGG1", "XRCC1"],
        "mechanisms": ["DNA_adducts", "ROS", "chromosomal_instability"],
        "tissues": [TissueType.LUNG, TissueType.NASAL_MUCOSA],
    },
    CancerSite.LUNG_SQUAMOUS: {
        "primary_metals": ["Arsenic", "Nickel", "Chromium(VI)"],
        "key_genes": ["CYP1A1", "CYP2A6", "GSTM1", "GSTP1", "OGG1", "XPC"],
        "mechanisms": ["DNA_adducts", "inflammation", "TP53_mutations"],
        "tissues": [TissueType.LUNG, TissueType.SKIN],
    },
    CancerSite.BLADDER_UROTHELIAL: {
        "primary_metals": ["Arsenic"],
        "key_genes": ["AS3MT", "GSTO1", "GSTM1", "GSTT1", "NAT2", "XRCC1"],
        "mechanisms": ["DNA_adducts", "oxidative_stress", "chromosomal_instability"],
        "tissues": [TissueType.BLADDER, TissueType.KIDNEY],
    },
    CancerSite.SKIN_SQUAMOUS: {
        "primary_metals": ["Arsenic"],
        "key_genes": ["CYP1A1", "GSTM1", "GSTP1", "OGG1", "XPC"],
        "mechanisms": ["DNA_adducts", "ROS", "immune_suppression"],
        "tissues": [TissueType.SKIN],
    },
    CancerSite.PROSTATE_ADENOCARCINOMA: {
        "primary_metals": ["Cadmium"],
        "key_genes": ["MT1A", "MT2A", "CYP1B1", "SRD5A2", "GSTP1"],
        "mechanisms": ["androgen_signaling", "DNA_repair_inhibition", "epigenetic_changes"],
        "tissues": [TissueType.PROSTATE],
    },
    CancerSite.KIDNEY_RENAL_CELL: {
        "primary_metals": ["Cadmium", "Lead"],
        "key_genes": ["MT1A", "MT2A", "CYP2E1", "GSTT1", "NQO1"],
        "mechanisms": ["oxidative_stress", "mitochondrial_dysfunction", "VHL_pathway"],
        "tissues": [TissueType.KIDNEY],
    },
    CancerSite.LIVER_HEPATOCELLULAR: {
        "primary_metals": ["Arsenic"],
        "key_genes": ["CYP1A2", "CYP2E1", "CYP3A4", "GSTM1", "GSTT1", "NAT2"],
        "mechanisms": ["metabolic_activation", "hepatotoxicity", "cirrhosis"],
        "tissues": [TissueType.LIVER],
    },
    CancerSite.NASAL_SINONASAL: {
        "primary_metals": ["Nickel", "Chromium(VI)"],
        "key_genes": ["CYP2A6", "CYP2A13", "GSTM1", "GSTP1"],
        "mechanisms": ["DNA_protein_crosslinks", "epigenetic_silencing"],
        "tissues": [TissueType.NASAL_MUCOSA],
    },
    CancerSite.PLEURAL_MESOTHELIOMA: {
        "primary_metals": ["Beryllium"],
        "key_genes": ["HLA-DPB1"],
        "mechanisms": ["immune_mediated", "fibrosis", "genomic_instability"],
        "tissues": [TissueType.LUNG, TissueType.PLEURA],
    },
}


@dataclass
class TissueMetabolismProfile:
    """Metabolic capacity profile for a tissue."""

    tissue: TissueType
    phase_i_capacity: dict[str, float] = field(default_factory=dict)
    phase_ii_capacity: dict[str, float] = field(default_factory=dict)
    dna_repair_capacity: dict[str, float] = field(default_factory=dict)
    metal_specific_pathways: dict[str, float] = field(default_factory=dict)
    overall_risk_score: float = 0.0


# ── Core Functions ──────────────────────────────────────────────────────────


def get_tissue_expression_level(gene: str, tissue: TissueType) -> str:
    """Get expression level of a gene in a specific tissue.

    Checks the qualitative TISSUE_EXPRESSION_PATTERNS first. For genes not
    in that set, falls back to GTEx quantitative weights when the tissue has
    GTEx data available.

    Args:
        gene: Gene symbol
        tissue: Tissue type

    Returns:
        Expression level: "high", "medium", "low", or "unknown"
    """
    if gene in TISSUE_EXPRESSION_PATTERNS:
        pattern = TISSUE_EXPRESSION_PATTERNS[gene]

        if tissue in pattern.get("high", []):
            return "high"
        elif tissue in pattern.get("medium", []):
            return "medium"
        elif tissue in pattern.get("low", []):
            return "low"

        # Check if ubiquitous
        if pattern.get("ubiquitous", False):
            return "medium"  # Baseline expression

        return "low"

    # Fallback to GTEx quantitative data
    gtex_name = _TISSUE_TO_GTEX.get(tissue)
    if gtex_name is not None:
        try:
            weights = get_tissue_weights(tissue)
        except (FileNotFoundError, ValueError):
            return "unknown"
        w = weights.get(gene)
        if w is not None:
            if w >= 0.5:
                return "high"
            elif w >= 0.2:
                return "medium"
            elif w > 0:
                return "low"
            else:
                return "low"

    return "unknown"


def get_genes_by_tissue(
    tissue: TissueType,
    min_level: str = "medium",
    *,
    threshold: float | None = None,
) -> list[str]:
    """Get genes expressed in a tissue above threshold.

    When *threshold* is ``None`` (default) the function uses the qualitative
    ``min_level`` filter against ``TISSUE_EXPRESSION_PATTERNS`` (plus GTEx
    fallback). When *threshold* is a float (0-1), it switches to purely
    quantitative GTEx weights for the 8 tissues that have data.

    Args:
        tissue: Target tissue
        min_level: Minimum expression level ("high", "medium", or "low").
            Ignored when *threshold* is set.
        threshold: Optional GTEx weight threshold (0-1). When set, returns
            genes whose GTEx weight >= this value.

    Returns:
        Sorted list of gene symbols
    """
    if threshold is not None:
        # Quantitative GTEx mode
        try:
            weights = get_tissue_weights(tissue)
        except (FileNotFoundError, ValueError):
            return []
        return sorted(g for g, w in weights.items() if w >= threshold)

    # Qualitative mode (original behaviour)
    level_priority = {"high": 3, "medium": 2, "low": 1, "unknown": 0}
    min_priority = level_priority.get(min_level, 2)

    genes = []
    for gene in TISSUE_EXPRESSION_PATTERNS:
        level = get_tissue_expression_level(gene, tissue)
        if level_priority.get(level, 0) >= min_priority:
            genes.append(gene)

    return sorted(genes)


def extract_tissue_subgraph(
    engine: GraphEngine,
    tissue: TissueType,
    include_carcinogens: list[str] | None = None,
) -> KnowledgeGraph:
    """Extract a tissue-specific subgraph from the full knowledge graph.

    Args:
        engine: Full graph engine
        tissue: Target tissue
        include_carcinogens: Optional list of carcinogen IDs to focus on

    Returns:
        Tissue-specific KnowledgeGraph

    Example:
        >>> engine = GraphEngine()
        >>> engine.load(full_graph)
        >>> lung_graph = extract_tissue_subgraph(engine, TissueType.LUNG)
    """
    # Get genes expressed in this tissue
    tissue_genes = set(get_genes_by_tissue(tissue, min_level="low"))

    # Add tissue node itself if it exists
    tissue_nodes = []
    tissue_node_id = tissue.value

    # Find all relevant nodes
    for node_id in engine.G.nodes:
        node_data = engine.get_data(node_id)
        if node_data is None:
            continue

        # Include tissue node
        if node_id == tissue_node_id or node_data.get("type") == "Tissue":
            if tissue_node_id in node_id or node_data.get("id") == tissue_node_id:
                tissue_nodes.append(node_id)

        # Include genes expressed in tissue
        if node_id in tissue_genes:
            tissue_nodes.append(node_id)

        # Include specified carcinogens
        if include_carcinogens and node_id in include_carcinogens:
            tissue_nodes.append(node_id)

        # Include carcinogens with tissue annotation
        if node_data.get("type") == "Carcinogen":
            exposure = node_data.get("exposure", "")
            if tissue.value in exposure.lower():
                tissue_nodes.append(node_id)

        # Include metabolites and adducts connected to tissue genes
        if node_data.get("type") in ("Metabolite", "DNA_Adduct"):
            # Will be filtered by edge connections below
            pass

    # Deduplicate
    tissue_nodes = list(set(tissue_nodes))

    # Find edges connecting tissue nodes
    tissue_edges = []
    for u, v, data in engine.G.edges(data=True):
        # Include edges between tissue-expressed nodes
        if u in tissue_nodes or v in tissue_nodes:
            tissue_edges.append((u, v, dict(data)))

    # Also include metabolites, adducts, and referenced carcinogens
    connected_nodes = set()
    carcinogen_refs = set()
    for u, v, data in tissue_edges:
        etype = data.get("type", "")
        if etype in ("ACTIVATES", "DETOXIFIES", "FORMS_ADDUCT"):
            if u not in tissue_nodes:
                connected_nodes.add(u)
            if v not in tissue_nodes:
                connected_nodes.add(v)
        # Track carcinogen references from edges
        if data.get("carcinogen"):
            carcinogen_refs.add(data["carcinogen"])

    # Include carcinogen nodes that are referenced by edges
    for ref in carcinogen_refs:
        if ref in engine.G:
            connected_nodes.add(ref)

    # Final node list
    final_nodes = list(set(tissue_nodes) | connected_nodes)

    # Build KnowledgeGraph
    nodes = []
    for node_id in final_nodes:
        data = engine.get_data(node_id)
        if data:
            # Add tissue expression annotation to detail field
            if data.get("type") in ("Enzyme", "Gene"):
                expr_level = get_tissue_expression_level(node_id, tissue)
                existing_detail = data.get("detail", "")
                tissue_info = f"[{tissue.value} expression: {expr_level}]"
                if existing_detail:
                    data["detail"] = f"{existing_detail} | {tissue_info}"
                else:
                    data["detail"] = tissue_info
            nodes.append(Node(**data))

    edges = []
    for u, v, data in tissue_edges:
        if u in final_nodes and v in final_nodes:
            from .models import Edge, EdgeType
            # Remove source/target from data to avoid duplication
            edge_data = {k: v for k, v in data.items() if k not in ("type", "source", "target")}
            edges.append(Edge(
                source=u,
                target=v,
                type=EdgeType(data.get("type", "CUSTOM")),
                **edge_data
            ))

    return KnowledgeGraph(nodes=nodes, edges=edges)


def build_cancer_site_subgraph(
    engine: GraphEngine,
    cancer_site: CancerSite,
) -> KnowledgeGraph:
    """Build a cancer site-specific subgraph with relevant metals and pathways.

    Args:
        engine: Full graph engine
        cancer_site: Target cancer site

    Returns:
        Cancer site-specific KnowledgeGraph
    """
    profile = CANCER_SITE_METAL_PROFILES.get(cancer_site)
    if not profile:
        return KnowledgeGraph(nodes=[], edges=[])

    # Get primary metals and key genes for this cancer site
    primary_metals = profile["primary_metals"]
    key_genes = profile["key_genes"]
    tissues = profile["tissues"]

    # Find nodes
    nodes_to_include = set()

    for node_id in engine.G.nodes:
        node_data = engine.get_data(node_id)
        if node_data is None:
            continue

        # Include primary carcinogens (by label matching)
        if node_data.get("type") == "Carcinogen":
            label = node_data.get("label", "").upper()
            for metal in primary_metals:
                if metal.upper().replace("(", "").replace(")", "") in label:
                    nodes_to_include.add(node_id)

        # Include key genes
        if node_id in key_genes:
            nodes_to_include.add(node_id)

        # Include tissue nodes
        for tissue in tissues:
            if tissue.value in node_id.lower():
                nodes_to_include.add(node_id)

    # Find connected nodes (metabolites, adducts)
    connected = set()
    for u, v, data in engine.G.edges(data=True):
        etype = data.get("type", "")
        if u in nodes_to_include or v in nodes_to_include:
            if etype in ("ACTIVATES", "DETOXIFIES", "FORMS_ADDUCT", "REPAIRS", "PATHWAY"):
                connected.add(u)
                connected.add(v)

    final_nodes = nodes_to_include | connected

    # Build graph
    nodes = []
    for node_id in final_nodes:
        data = engine.get_data(node_id)
        if data:
            # Add cancer site annotation to detail field
            is_primary = node_id in key_genes or node_id in nodes_to_include
            existing_detail = data.get("detail", "")
            site_info = f"[Cancer site: {cancer_site.value}"
            if is_primary:
                site_info += " (primary)]"
            else:
                site_info += "]"
            if existing_detail:
                data["detail"] = f"{existing_detail} | {site_info}"
            else:
                data["detail"] = site_info
            nodes.append(Node(**data))

    edges = []
    for u, v, data in engine.G.edges(data=True):
        if u in final_nodes and v in final_nodes:
            from .models import Edge, EdgeType
            # Remove source/target from data to avoid duplication
            edge_data = {k: v for k, v in data.items() if k not in ("type", "source", "target")}
            edges.append(Edge(
                source=u,
                target=v,
                type=EdgeType(data.get("type", "CUSTOM")),
                **edge_data
            ))

    return KnowledgeGraph(nodes=nodes, edges=edges)


def tissue_metabolism_chain(
    engine: GraphEngine,
    carcinogen_id: str,
    tissue: TissueType,
) -> MetabolismChain:
    """Extract metabolism chain for a carcinogen in a specific tissue.

    Only includes enzymes and pathways relevant to the specified tissue.

    Args:
        engine: Graph engine
        carcinogen_id: Carcinogen to analyze
        tissue: Target tissue

    Returns:
        Tissue-filtered MetabolismChain
    """
    # Get full metabolism chain
    full_chain = metabolism_chain(engine, carcinogen_id)

    # Get tissue-expressed genes
    tissue_genes = set(get_genes_by_tissue(tissue, min_level="low"))

    # Filter chain edges to only those involving tissue-expressed enzymes
    filtered_edges = []
    for edge in full_chain.edges:
        source = edge.get("source", "")
        target = edge.get("target", "")

        # Include if source or target is tissue-expressed or not a gene
        source_in_tissue = source in tissue_genes or source == carcinogen_id
        target_node = engine.get_data(target)
        target_is_not_gene = target_node is not None and target_node.get("type") != "Enzyme"

        if source_in_tissue or target_is_not_gene:
            # Add tissue expression annotation
            edge_copy = dict(edge)
            if source in tissue_genes:
                expr = get_tissue_expression_level(source, tissue)
                edge_copy["source_tissue_expression"] = expr
            filtered_edges.append(edge_copy)

    # Create filtered chain
    chain = MetabolismChain(carcinogen_id=carcinogen_id)
    chain.edges = filtered_edges

    # Get unique nodes from filtered edges
    node_ids = set()
    for edge in filtered_edges:
        node_ids.add(edge.get("source"))
        node_ids.add(edge.get("target"))
    chain.node_ids = sorted([n for n in node_ids if n])

    return chain


def calculate_tissue_metabolism_capacity(
    engine: GraphEngine,
    tissue: TissueType,
    carcinogen_id: str | None = None,
) -> TissueMetabolismProfile:
    """Calculate the metabolic capacity profile for a tissue.

    Args:
        engine: Graph engine
        tissue: Target tissue
        carcinogen_id: Optional carcinogen to focus on

    Returns:
        TissueMetabolismProfile with capacity scores
    """
    profile = TissueMetabolismProfile(tissue=tissue)

    # Get tissue-expressed genes with expression weights
    tissue_genes = get_genes_by_tissue(tissue, min_level="low")
    expr_weights = {"high": 1.0, "medium": 0.7, "low": 0.4, "unknown": 0.2}

    # Try to load GTEx quantitative weights for more precise scoring
    gtex_weights: dict[str, float] | None = None
    if tissue in _TISSUE_TO_GTEX:
        try:
            gtex_weights = get_tissue_weights(tissue)
        except (FileNotFoundError, ValueError):
            pass

    for gene in tissue_genes:
        # Prefer GTEx quantitative weight when available
        if gtex_weights is not None and gene in gtex_weights:
            weight = gtex_weights[gene]
        else:
            level = get_tissue_expression_level(gene, tissue)
            weight = expr_weights.get(level, 0.2)

        node_data = engine.get_data(gene)
        if not node_data:
            continue

        phase = node_data.get("phase")
        role = node_data.get("role", "")

        # Categorize by phase
        if phase == "I":
            profile.phase_i_capacity[gene] = weight
        elif phase == "II":
            profile.phase_ii_capacity[gene] = weight
        elif role == "Repair":
            profile.dna_repair_capacity[gene] = weight

        # Metal-specific pathways
        if gene in ("AS3MT", "GSTO1", "GSTO2"):
            profile.metal_specific_pathways["arsenic_metabolism"] = max(
                profile.metal_specific_pathways.get("arsenic_metabolism", 0),
                weight
            )
        elif gene in ("MT1A", "MT2A"):
            profile.metal_specific_pathways["cadmium_detoxification"] = max(
                profile.metal_specific_pathways.get("cadmium_detoxification", 0),
                weight
            )
        elif gene == "HLA-DPB1":
            profile.metal_specific_pathways["beryllium_sensitivity"] = weight

    # Calculate overall risk score
    # Higher activation capacity + lower repair capacity = higher risk
    activation_score = sum(profile.phase_i_capacity.values())
    detox_score = sum(profile.phase_ii_capacity.values())
    repair_score = sum(profile.dna_repair_capacity.values())

    # Risk score: activation / (detox + repair + 1) normalized
    denominator = detox_score + repair_score + 1
    raw_score = activation_score / denominator
    profile.overall_risk_score = min(raw_score, 10.0)  # Cap at 10

    return profile


def compare_tissue_metabolism(
    engine: GraphEngine,
    carcinogen_id: str,
    tissues: list[TissueType],
) -> dict[str, Any]:
    """Compare metabolism across multiple tissues for a carcinogen.

    Args:
        engine: Graph engine
        carcinogen_id: Carcinogen to analyze
        tissues: List of tissues to compare

    Returns:
        Comparison dictionary with tissue rankings
    """
    comparisons = []

    for tissue in tissues:
        chain = tissue_metabolism_chain(engine, carcinogen_id, tissue)
        capacity = calculate_tissue_metabolism_capacity(engine, tissue, carcinogen_id)

        comparisons.append({
            "tissue": tissue.value,
            "activation_edges": len(chain.activation_edges),
            "detox_edges": len(chain.detox_edges),
            "adduct_edges": len(chain.adduct_edges),
            "repair_edges": len(chain.repair_edges),
            "phase_i_capacity": sum(capacity.phase_i_capacity.values()),
            "phase_ii_capacity": sum(capacity.phase_ii_capacity.values()),
            "dna_repair_capacity": sum(capacity.dna_repair_capacity.values()),
            "risk_score": capacity.overall_risk_score,
        })

    # Rank by risk score
    comparisons.sort(key=lambda x: x["risk_score"], reverse=True)

    return {
        "carcinogen": carcinogen_id,
        "tissue_rankings": comparisons,
        "highest_risk_tissue": comparisons[0]["tissue"] if comparisons else None,
        "lowest_risk_tissue": comparisons[-1]["tissue"] if comparisons else None,
    }


def get_metal_cancer_site_associations(metal: str) -> list[CancerSite]:
    """Get cancer sites associated with a specific metal.

    Args:
        metal: Metal name (e.g., "Arsenic", "Cadmium")

    Returns:
        List of associated CancerSite values
    """
    sites = []
    for site, profile in CANCER_SITE_METAL_PROFILES.items():
        for primary_metal in profile["primary_metals"]:
            if metal.lower() in primary_metal.lower():
                sites.append(site)
                break
    return sites


def build_multi_tissue_pathway(
    engine: GraphEngine,
    carcinogen_id: str,
    source_tissue: TissueType,
    target_tissue: TissueType,
) -> dict[str, Any]:
    """Trace a carcinogen pathway from absorption to target organ.

    Args:
        engine: Graph engine
        carcinogen_id: Carcinogen to trace
        source_tissue: Entry tissue (e.g., LUNG for inhalation)
        target_tissue: Target tissue for carcinogenesis

    Returns:
        Multi-tissue pathway analysis
    """
    # Get metabolism chains for each tissue
    source_chain = tissue_metabolism_chain(engine, carcinogen_id, source_tissue)
    target_chain = tissue_metabolism_chain(engine, carcinogen_id, target_tissue)

    # Identify transport routes (if any transporters in graph)
    transporters = []
    for node_id in engine.G.nodes:
        node_data = engine.get_data(node_id)
        if node_data and node_data.get("role") == "Transport":
            transporters.append(node_id)

    return {
        "carcinogen": carcinogen_id,
        "source_tissue": source_tissue.value,
        "target_tissue": target_tissue.value,
        "source_metabolism": {
            "activation_count": len(source_chain.activation_edges),
            "detox_count": len(source_chain.detox_edges),
        },
        "target_metabolism": {
            "activation_count": len(target_chain.activation_edges),
            "adduct_count": len(target_chain.adduct_edges),
            "repair_count": len(target_chain.repair_edges),
        },
        "transporters_available": transporters,
        "notes": (
            f"{carcinogen_id} enters via {source_tissue.value}, "
            f"metabolized ({len(source_chain.activation_edges)} activation steps), "
            f"causes damage in {target_tissue.value} "
            f"({len(target_chain.adduct_edges)} adduct-forming pathways)"
        ),
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m ExposoGraph.tissue_subgraphs",
        description="GTEx-backed tissue filtering and reporting for ExposoGraph graphs.",
    )
    parser.add_argument("--tissue", help="GTEx tissue name used for JSON graph filtering.")
    parser.add_argument("--threshold", type=float, default=0.05, help="Minimum tissue weight for enzyme inclusion.")
    parser.add_argument("--input", help="Input graph JSON path for --tissue filtering.")
    parser.add_argument("--output", help="Output path for tissue-filtered graph JSON.")
    parser.add_argument("--report", help="Print the formatted tissue report for a GTEx tissue.")
    parser.add_argument("--list-tissues", action="store_true", help="List available GTEx tissues.")
    parser.add_argument("--weights-only", help="Print GTEx tissue weights as JSON for a GTEx tissue.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Small CLI compatible with the original 03_tissue_views workflow."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.list_tissues:
        print("\n".join(get_available_tissues()))
        return 0

    if args.weights_only:
        print(json.dumps(get_tissue_weights(args.weights_only), indent=2, sort_keys=True))
        return 0

    if args.report:
        print(generate_tissue_report(args.report, threshold=args.threshold))
        return 0

    if args.tissue:
        if not args.input or not args.output:
            parser.error("--tissue filtering requires both --input and --output.")

        graph_data = json.loads(Path(args.input).read_text(encoding="utf-8"))
        filtered = filter_graph_by_tissue(graph_data, args.tissue, threshold=args.threshold)
        Path(args.output).write_text(
            json.dumps(filtered, indent=2),
            encoding="utf-8",
        )
        print(f"Wrote tissue-filtered graph for {args.tissue} to {args.output}")
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
