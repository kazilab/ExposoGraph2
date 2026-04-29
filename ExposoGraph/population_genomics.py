"""Population pharmacogenomics module for ancestry-aware risk assessment.

Provides allele frequency data across major population groups and utilities
for population-stratified pharmacogenomic analysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from .models import KnowledgeGraph, Node, NodeType


class Population(str, Enum):
    """Major ancestry groups from 1000 Genomes Project superpopulations."""

    AFR = "AFR"  # African
    AMR = "AMR"  # Admixed American
    EAS = "EAS"  # East Asian
    EUR = "EUR"  # European
    SAS = "SAS"  # South Asian
    GLOBAL = "GLOBAL"  # All populations combined


@dataclass
class AlleleFrequency:
    """Allele frequency data for a specific population."""

    allele: str
    population: Population
    frequency: float  # 0.0 to 1.0
    sample_size: int
    source: str = "1000G"  # 1000 Genomes, gnomAD, etc.
    source_version: str = "Phase 3"
    confidence: str = "High"  # High, Medium, Low
    pmid_references: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not 0.0 <= self.frequency <= 1.0:
            raise ValueError(f"Frequency must be between 0 and 1, got {self.frequency}")


@dataclass
class PopulationRiskProfile:
    """Risk profile for a specific gene in a population."""

    gene: str
    population: Population
    risk_allele: str
    risk_allele_frequency: float
    protective_allele: str | None = None
    protective_allele_frequency: float | None = None
    odds_ratio: float | None = None  # From meta-analyses
    cancer_sites: list[str] = field(default_factory=list)
    metal_exposure: str | None = None  # Associated metal
    pmid_references: list[str] = field(default_factory=list)


# ── Population Allele Frequency Database ─────────────────────────────────────

# Curated from 1000 Genomes Phase 3, gnomAD v3, and PharmGKB population summaries
# Key pharmacogenes with metal metabolism relevance

POPULATION_ALLELE_FREQUENCIES: dict[str, list[AlleleFrequency]] = {
    # AS3MT - Arsenic metabolism (critical for arsenic methylation efficiency)
    "AS3MT": [
        # M287T (rs11191439) - affects MMA/DMA ratios
        AlleleFrequency("T (287Thr)", Population.EUR, 0.095, 503, pmid_references=["28537708", "19371612"]),
        AlleleFrequency("T (287Thr)", Population.EAS, 0.299, 504, pmid_references=["28537708", "19371612"]),
        AlleleFrequency("T (287Thr)", Population.SAS, 0.182, 489, pmid_references=["28537708", "19371612"]),
        AlleleFrequency("T (287Thr)", Population.AFR, 0.045, 661, pmid_references=["28537708", "19371612"]),
        AlleleFrequency("T (287Thr)", Population.AMR, 0.112, 347, pmid_references=["28537708", "19371612"]),
        # Q144R (rs11191442)
        AlleleFrequency("R (144Arg)", Population.EUR, 0.001, 503),
        AlleleFrequency("R (144Arg)", Population.EAS, 0.156, 504),
        AlleleFrequency("R (144Arg)", Population.SAS, 0.002, 489),
    ],

    # GSTO1 - Arsenic metabolism (reduces MMA(V) to MMA(III))
    "GSTO1": [
        # A140D (rs4925) - alters arsenic metabolism
        AlleleFrequency("D (140Asp)", Population.EUR, 0.317, 503, pmid_references=["21798077"]),
        AlleleFrequency("D (140Asp)", Population.EAS, 0.091, 504, pmid_references=["21798077"]),
        AlleleFrequency("D (140Asp)", Population.SAS, 0.174, 489, pmid_references=["21798077"]),
        AlleleFrequency("D (140Asp)", Population.AFR, 0.243, 661, pmid_references=["21798077"]),
        AlleleFrequency("D (140Asp)", Population.AMR, 0.294, 347, pmid_references=["21798077"]),
    ],

    # GSTO2 - Arsenic metabolism
    "GSTO2": [
        # N142D (rs156697)
        AlleleFrequency("D (142Asp)", Population.EUR, 0.232, 503),
        AlleleFrequency("D (142Asp)", Population.EAS, 0.299, 504),
        AlleleFrequency("D (142Asp)", Population.SAS, 0.311, 489),
        AlleleFrequency("D (142Asp)", Population.AFR, 0.425, 661),
    ],

    # MT1A - Cadmium and zinc homeostasis
    "MT1A": [
        # rs11076161 promoter variant
        AlleleFrequency("A (promoter)", Population.EUR, 0.452, 503),
        AlleleFrequency("A (promoter)", Population.EAS, 0.342, 504),
        AlleleFrequency("A (promoter)", Population.SAS, 0.378, 489),
        AlleleFrequency("A (promoter)", Population.AFR, 0.512, 661),
    ],

    # MT2A - Cadmium detoxification
    "MT2A": [
        # rs1610216 - represented in MT haplotypes linked to OSCC susceptibility
        AlleleFrequency("T (rs1610216)", Population.EUR, 0.387, 503, pmid_references=["33340085"]),
        AlleleFrequency("T (rs1610216)", Population.EAS, 0.291, 504, pmid_references=["33340085"]),
        AlleleFrequency("T (rs1610216)", Population.SAS, 0.356, 489, pmid_references=["33340085"]),
    ],

    # HLA-DPB1 - Beryllium sensitivity
    "HLA-DPB1": [
        # Glu69 (rs9276293) - beryllium-specific risk allele
        AlleleFrequency("Glu69+", Population.EUR, 0.384, 503, pmid_references=["10558913"]),
        AlleleFrequency("Glu69+", Population.EAS, 0.315, 504, pmid_references=["10558913"]),
        AlleleFrequency("Glu69+", Population.AFR, 0.428, 661, pmid_references=["10558913"]),
        AlleleFrequency("Glu69+", Population.AMR, 0.401, 347, pmid_references=["10558913"]),
        AlleleFrequency("Glu69+", Population.SAS, 0.312, 489, pmid_references=["10558913"]),
    ],

    # GSTM1 - Deletion polymorphism (null genotype)
    "GSTM1": [
        AlleleFrequency("Null (del)", Population.EUR, 0.500, 503, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.EAS, 0.512, 504, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.SAS, 0.334, 489, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.AFR, 0.367, 661, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.AMR, 0.406, 347, pmid_references=["21798077", "32650499"]),
    ],

    # GSTT1 - Deletion polymorphism
    "GSTT1": [
        AlleleFrequency("Null (del)", Population.EUR, 0.176, 503, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.EAS, 0.571, 504, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.SAS, 0.216, 489, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.AFR, 0.224, 661, pmid_references=["21798077", "32650499"]),
        AlleleFrequency("Null (del)", Population.AMR, 0.207, 347, pmid_references=["21798077", "32650499"]),
    ],

    # CYP1A1 - PAH and estrogen metabolism
    "CYP1A1": [
        # Ile462Val (rs1048943) - high activity variant
        AlleleFrequency("Val (462Val)", Population.EUR, 0.045, 503, pmid_references=["15647817"]),
        AlleleFrequency("Val (462Val)", Population.EAS, 0.082, 504, pmid_references=["15647817"]),
        AlleleFrequency("Val (462Val)", Population.SAS, 0.068, 489, pmid_references=["15647817"]),
        AlleleFrequency("Val (462Val)", Population.AFR, 0.013, 661, pmid_references=["15647817"]),
        # MspI (rs4646903)
        AlleleFrequency("T (MspI)", Population.EUR, 0.289, 503),
        AlleleFrequency("T (MspI)", Population.EAS, 0.363, 504),
        AlleleFrequency("T (MspI)", Population.SAS, 0.421, 489),
    ],

    # NQO1 - Quinone reduction (redox cycling protection)
    "NQO1": [
        # Pro187Ser (rs1800566) - reduced activity
        AlleleFrequency("Ser (187Ser)", Population.EUR, 0.208, 503, pmid_references=["23860519"]),
        AlleleFrequency("Ser (187Ser)", Population.EAS, 0.156, 504, pmid_references=["23860519"]),
        AlleleFrequency("Ser (187Ser)", Population.SAS, 0.192, 489, pmid_references=["23860519"]),
        AlleleFrequency("Ser (187Ser)", Population.AFR, 0.076, 661, pmid_references=["23860519"]),
        AlleleFrequency("Ser (187Ser)", Population.AMR, 0.167, 347, pmid_references=["23860519"]),
    ],

    # OGG1 - Oxidative DNA damage repair
    "OGG1": [
        # Ser326Cys (rs1052133) - reduced repair capacity
        AlleleFrequency("Cys (326Cys)", Population.EUR, 0.233, 503, pmid_references=["25588927"]),
        AlleleFrequency("Cys (326Cys)", Population.EAS, 0.442, 504, pmid_references=["25588927"]),
        AlleleFrequency("Cys (326Cys)", Population.SAS, 0.378, 489, pmid_references=["25588927"]),
        AlleleFrequency("Cys (326Cys)", Population.AFR, 0.156, 661, pmid_references=["25588927"]),
        AlleleFrequency("Cys (326Cys)", Population.AMR, 0.244, 347, pmid_references=["25588927"]),
    ],

    # XRCC1 - DNA repair
    "XRCC1": [
        # Arg399Gln (rs25487)
        AlleleFrequency("Gln (399Gln)", Population.EUR, 0.364, 503, pmid_references=["20486219"]),
        AlleleFrequency("Gln (399Gln)", Population.EAS, 0.356, 504, pmid_references=["20486219"]),
        AlleleFrequency("Gln (399Gln)", Population.SAS, 0.421, 489, pmid_references=["20486219"]),
        AlleleFrequency("Gln (399Gln)", Population.AFR, 0.267, 661, pmid_references=["20486219"]),
    ],

    # XPC - NER capacity
    "XPC": [
        # Lys939Gln (rs2228001) - reduced NER
        AlleleFrequency("Gln (939Gln)", Population.EUR, 0.381, 503, pmid_references=["22592359"]),
        AlleleFrequency("Gln (939Gln)", Population.EAS, 0.294, 504, pmid_references=["22592359"]),
        AlleleFrequency("Gln (939Gln)", Population.SAS, 0.412, 489, pmid_references=["22592359"]),
        AlleleFrequency("Gln (939Gln)", Population.AFR, 0.298, 661, pmid_references=["22592359"]),
    ],
}


# Population-specific risk associations for metal exposure
POPULATION_RISK_PROFILES: list[PopulationRiskProfile] = [
    # Arsenic-related risks
    PopulationRiskProfile(
        gene="AS3MT",
        population=Population.EAS,
        risk_allele="T (287Thr)",
        risk_allele_frequency=0.299,
        odds_ratio=1.52,
        cancer_sites=["bladder", "lung", "skin"],
        metal_exposure="Arsenic",
        pmid_references=["28537708", "19680750"],
    ),
    PopulationRiskProfile(
        gene="GSTO1",
        population=Population.EUR,
        risk_allele="D (140Asp)",
        risk_allele_frequency=0.317,
        odds_ratio=1.38,
        cancer_sites=["lung", "bladder"],
        metal_exposure="Arsenic",
        pmid_references=["21798077"],
    ),
    PopulationRiskProfile(
        gene="GSTM1",
        population=Population.GLOBAL,
        risk_allele="Null (del)",
        risk_allele_frequency=0.50,
        odds_ratio=1.42,
        cancer_sites=["lung", "bladder"],
        metal_exposure="Arsenic",
        pmid_references=["21798077", "32650499"],
    ),
    PopulationRiskProfile(
        gene="GSTT1",
        population=Population.EAS,
        risk_allele="Null (del)",
        risk_allele_frequency=0.571,
        odds_ratio=1.68,
        cancer_sites=["bladder"],
        metal_exposure="Arsenic",
        pmid_references=["21798077", "32650499"],
    ),
    # Cadmium-related
    PopulationRiskProfile(
        gene="MT1A",
        population=Population.EUR,
        risk_allele="A (promoter)",
        risk_allele_frequency=0.452,
        odds_ratio=1.25,
        cancer_sites=["kidney", "prostate"],
        metal_exposure="Cadmium",
        pmid_references=["33340085"],
    ),
    # Beryllium-related
    PopulationRiskProfile(
        gene="HLA-DPB1",
        population=Population.EUR,
        risk_allele="Glu69+",
        risk_allele_frequency=0.384,
        odds_ratio=3.2,
        cancer_sites=["lung"],
        metal_exposure="Beryllium",
        pmid_references=["10558913", "17038120"],
    ),
    # Oxidative stress/DNA repair
    PopulationRiskProfile(
        gene="OGG1",
        population=Population.EAS,
        risk_allele="Cys (326Cys)",
        risk_allele_frequency=0.442,
        odds_ratio=1.35,
        cancer_sites=["lung", "esophagus"],
        metal_exposure="Multiple (As, Cr, Ni)",
        pmid_references=["25588927"],
    ),
    PopulationRiskProfile(
        gene="NQO1",
        population=Population.EUR,
        risk_allele="Ser (187Ser)",
        risk_allele_frequency=0.208,
        odds_ratio=1.28,
        cancer_sites=["lung", "colorectal"],
        metal_exposure="Multiple",
        pmid_references=["23860519"],
    ),
]


# ── Functions ───────────────────────────────────────────────────────────────


def get_allele_frequencies(gene: str) -> list[AlleleFrequency]:
    """Get all allele frequency data for a gene.

    Args:
        gene: Gene symbol (e.g., "AS3MT", "GSTM1")

    Returns:
        List of AlleleFrequency records across populations

    Example:
        >>> freqs = get_allele_frequencies("AS3MT")
        >>> eas_freq = [f for f in freqs if f.population == Population.EAS]
    """
    return POPULATION_ALLELE_FREQUENCIES.get(gene.upper(), [])


def get_population_frequency(
    gene: str, allele: str, population: Population
) -> AlleleFrequency | None:
    """Get frequency for a specific gene/allele/population combination.

    Args:
        gene: Gene symbol
        allele: Allele identifier
        population: Population group

    Returns:
        AlleleFrequency if found, None otherwise
    """
    for freq in POPULATION_ALLELE_FREQUENCIES.get(gene.upper(), []):
        if freq.allele == allele and freq.population == population:
            return freq
    return None


def get_population_risk_profiles(
    gene: str | None = None,
    population: Population | None = None,
    metal: str | None = None,
) -> list[PopulationRiskProfile]:
    """Get population-specific risk profiles with optional filtering.

    Args:
        gene: Filter by gene symbol
        population: Filter by population
        metal: Filter by metal exposure (e.g., "Arsenic", "Cadmium")

    Returns:
        List of matching PopulationRiskProfile records
    """
    results = POPULATION_RISK_PROFILES

    if gene:
        results = [r for r in results if r.gene.upper() == gene.upper()]
    if population:
        results = [r for r in results if r.population == population]
    if metal:
        results = [
            r for r in results if r.metal_exposure and metal.lower() in r.metal_exposure.lower()
        ]

    return results


def calculate_population_risk_score(
    gene: str,
    population: Population,
    diplotype: str | None = None,
) -> dict[str, Any]:
    """Calculate a population-weighted risk score for a gene.

    This combines allele frequency with odds ratio to estimate
    population-level attributable risk.

    Args:
        gene: Gene symbol
        population: Population group
        diplotype: Optional specific diplotype (for example, ``*1/*4``)

    Returns:
        Dictionary with risk metrics

    Example:
        >>> score = calculate_population_risk_score("AS3MT", Population.EAS)
        >>> print(f"Population AF: {score['risk_allele_frequency']:.3f}")
    """
    profiles = get_population_risk_profiles(gene=gene, population=population)

    if not profiles:
        return {
            "gene": gene,
            "population": population.value,
            "risk_score": None,
            "message": "No risk profile available for this gene/population combination",
        }

    # Use the highest OR profile if multiple exist
    profile = max(profiles, key=lambda p: p.odds_ratio or 1.0)

    # Calculate population attributable fraction (PAF)
    af = profile.risk_allele_frequency
    or_val = profile.odds_ratio or 1.0
    # PAF = AF * (OR - 1) / (1 + AF * (OR - 1))
    paf = (af * (or_val - 1)) / (1 + af * (or_val - 1)) if or_val != 1 else 0

    return {
        "gene": gene,
        "population": population.value,
        "risk_allele": profile.risk_allele,
        "risk_allele_frequency": af,
        "odds_ratio": or_val,
        "population_attributable_fraction": round(paf, 4),
        "cancer_sites": profile.cancer_sites,
        "metal_exposure": profile.metal_exposure,
        "references": profile.pmid_references,
    }


def compare_population_risks(gene: str) -> dict[str, Any]:
    """Compare risk profiles across all populations for a gene.

    Args:
        gene: Gene symbol

    Returns:
        Dictionary with cross-population comparison
    """
    populations = [p for p in Population if p != Population.GLOBAL]
    comparisons = []

    for pop in populations:
        score = calculate_population_risk_score(gene, pop)
        if score.get("risk_score") is not None or score.get("risk_allele_frequency"):
            comparisons.append({
                "population": pop.value,
                "frequency": score.get("risk_allele_frequency"),
                "odds_ratio": score.get("odds_ratio"),
                "paf": score.get("population_attributable_fraction"),
            })

    if not comparisons:
        return {"gene": gene, "comparisons": [], "max_risk_population": None}

    # Find population with highest risk allele frequency
    max_freq_pop = max(comparisons, key=lambda x: x.get("frequency") or 0)

    return {
        "gene": gene,
        "comparisons": comparisons,
        "max_risk_population": max_freq_pop["population"],
        "frequency_range": (
            min(c["frequency"] for c in comparisons if c["frequency"]),
            max(c["frequency"] for c in comparisons if c["frequency"]),
        ),
    }


def get_genes_by_population_prevalence(
    population: Population,
    min_frequency: float = 0.10,
) -> list[str]:
    """Get list of genes with risk alleles above frequency threshold in a population.

    Useful for prioritizing which pharmacogenes to test in a population.

    Args:
        population: Population group
        min_frequency: Minimum allele frequency (0.0-1.0)

    Returns:
        List of gene symbols
    """
    genes = set()
    for gene, freqs in POPULATION_ALLELE_FREQUENCIES.items():
        for freq in freqs:
            if freq.population == population and freq.frequency >= min_frequency:
                genes.add(gene)
                break
    return sorted(genes)


def annotate_nodes_with_population_data(
    nodes: list[Node],
    population: Population,
) -> list[Node]:
    """Annotate gene/enzyme nodes with population frequency data.

    Args:
        nodes: List of Node objects
        population: Population to annotate for

    Returns:
        New list of nodes with population data added to detail field
    """
    annotated = []
    for node in nodes:
        if (
            node.type in (NodeType.ENZYME, NodeType.GENE)
            and node.id in POPULATION_ALLELE_FREQUENCIES
        ):
            freqs = get_allele_frequencies(node.id)
            pop_freqs = [f for f in freqs if f.population == population]

            if pop_freqs:
                # Format annotation for detail field
                freq_strs = [f"{f.allele}: {f.frequency:.1%}" for f in pop_freqs[:3]]
                pop_detail = f"[{population.value}] " + "; ".join(freq_strs)

                # Create new node with updated detail
                node_data = node.model_dump()
                existing_detail = node_data.get("detail", "")
                if existing_detail:
                    node_data["detail"] = f"{existing_detail} | {pop_detail}"
                else:
                    node_data["detail"] = pop_detail
                annotated.append(Node(**node_data))
            else:
                annotated.append(node)
        else:
            annotated.append(node)

    return annotated


def build_population_aware_panel(population: Population) -> KnowledgeGraph:
    """Build a gene panel prioritized for a specific population.

    Prioritizes genes where:
    1. Risk alleles have frequency > 10%
    2. Population-specific risk associations exist
    3. Metal metabolism relevance is established

    Args:
        population: Target population

    Returns:
        KnowledgeGraph with prioritized gene panel
    """
    from .reference_data import build_full_panel

    # Get full panel
    full_panel = build_full_panel()

    # Get population-prevalent genes
    prevalent_genes = set(get_genes_by_population_prevalence(population, min_frequency=0.10))

    # Get genes with population risk profiles
    risk_genes = set()
    for profile in POPULATION_RISK_PROFILES:
        if profile.population == population or profile.population == Population.GLOBAL:
            risk_genes.add(profile.gene)

    # Priority genes (intersection of prevalence and risk)
    priority_genes = prevalent_genes & risk_genes

    # Build annotated nodes with population data and priority
    final_nodes = []
    for node in full_panel.nodes:
        if (
            node.type in (NodeType.ENZYME, NodeType.GENE)
            and node.id in POPULATION_ALLELE_FREQUENCIES
        ):
            freqs = get_allele_frequencies(node.id)
            pop_freqs = [f for f in freqs if f.population == population]

            node_data = node.model_dump()
            annotations = []

            # Add frequency info
            if pop_freqs:
                freq_strs = [f"{f.allele}: {f.frequency:.1%}" for f in pop_freqs[:3]]
                annotations.append(f"[{population.value}] " + "; ".join(freq_strs))

            # Add priority info
            if node.id in priority_genes:
                node_data["tier"] = 1  # Override as Tier 1 for this population
                annotations.append("Population Priority: HIGH")
            elif node.id in prevalent_genes or node.id in risk_genes:
                annotations.append("Population Priority: Moderate")

            # Update detail field
            existing_detail = node_data.get("detail", "")
            pop_info = " | ".join(annotations)
            if existing_detail:
                node_data["detail"] = f"{existing_detail} | {pop_info}"
            else:
                node_data["detail"] = pop_info

            final_nodes.append(Node(**node_data))
        else:
            final_nodes.append(node)

    return KnowledgeGraph(nodes=final_nodes, edges=[])
