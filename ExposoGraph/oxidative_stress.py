"""Oxidative stress biomarkers and metal-induced ROS mechanisms.

This module provides:
- Oxidative stress biomarker definitions
- Metal-specific ROS generation pathways
- Antioxidant enzyme response profiles
- Lipid peroxidation and DNA damage markers
- Integrated oxidative stress scoring
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from .models import Edge, EdgeType, KnowledgeGraph, Node, NodeType


class OSMarkerType(str, Enum):
    """Types of oxidative stress biomarkers."""

    ROS = "reactive_oxygen_species"
    RNS = "reactive_nitrogen_species"
    LIPID_PEROXIDATION = "lipid_peroxidation"
    PROTEIN_OXIDATION = "protein_oxidation"
    DNA_OXIDATION = "dna_oxidation"
    ANTIOXIDANT_ENZYME = "antioxidant_enzyme"
    GSH_METABOLISM = "glutathione_metabolism"
    IRON_SULFUR = "iron_sulfur_cluster"


class ReactiveSpecies(str, Enum):
    """Reactive oxygen and nitrogen species."""

    SUPEROXIDE = "superoxide_anion"  # O2•-
    HYDROGEN_PEROXIDE = "hydrogen_peroxide"  # H2O2
    HYDROXYL_RADICAL = "hydroxyl_radical"  # •OH
    SINGLET_OXYGEN = "singlet_oxygen"  # 1O2
    PEROXYNITRITE = "peroxynitrite"  # ONOO-
    NITRIC_OXIDE = "nitric_oxide"  # NO•
    HYPOCHLOROUS_ACID = "hypochlorous_acid"  # HOCl


@dataclass
class OxidativeStressMarker:
    """Oxidative stress biomarker with metadata."""

    id: str
    name: str
    abbreviation: str
    marker_type: OSMarkerType
    detection_method: str  # e.g., "ELISA", "LC-MS/MS", "spectrophotometric"
    sample_matrix: str  # e.g., "urine", "plasma", "tissue", "exhaled_breath"
    reference_range_low: float | None = None
    reference_range_high: float | None = None
    unit: str = ""
    metal_induced: bool = True
    associated_metals: list[str] = field(default_factory=list)
    pmid_references: list[str] = field(default_factory=list)


@dataclass
class ROSGenerationPathway:
    """Pathway for metal-induced ROS generation."""

    metal: str
    primary_ros: ReactiveSpecies
    generation_mechanism: str
    cellular_compartment: str  # mitochondria, cytosol, ER, nucleus
    enzyme_involved: str | None = None
    requires_fenton: bool = False
    requires_redox_cycling: bool = False
    metal_specific_notes: str = ""


# ── Oxidative Stress Biomarker Database ────────────────────────────────────

OXIDATIVE_STRESS_MARKERS: dict[str, OxidativeStressMarker] = {
    # Lipid peroxidation markers
    "MDA": OxidativeStressMarker(
        id="MDA",
        name="Malondialdehyde",
        abbreviation="MDA",
        marker_type=OSMarkerType.LIPID_PEROXIDATION,
        detection_method="TBARS assay, LC-MS/MS",
        sample_matrix="plasma, urine, tissue",
        unit="μmol/L",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury", "Nickel"],
        pmid_references=["19825587", "22023223"],
    ),
    "4_HNE": OxidativeStressMarker(
        id="4-HNE",
        name="4-Hydroxynonenal",
        abbreviation="4-HNE",
        marker_type=OSMarkerType.LIPID_PEROXIDATION,
        detection_method="ELISA, LC-MS/MS",
        sample_matrix="plasma, urine, tissue",
        unit="ng/mL",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium", "Lead", "Mercury"],
        pmid_references=["22023223", "24577282"],
    ),
    "8_isoprostane": OxidativeStressMarker(
        id="8-iso-PGF2α",
        name="8-Isoprostane",
        abbreviation="8-iso-PGF2α",
        marker_type=OSMarkerType.LIPID_PEROXIDATION,
        detection_method="ELISA, LC-MS/MS",
        sample_matrix="urine, plasma, exhaled breath condensate",
        unit="pg/mL",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Nickel"],
        pmid_references=["19825587", "15626646"],
    ),
    "TBARS": OxidativeStressMarker(
        id="TBARS",
        name="Thiobarbituric Acid Reactive Substances",
        abbreviation="TBARS",
        marker_type=OSMarkerType.LIPID_PEROXIDATION,
        detection_method="spectrophotometric",
        sample_matrix="plasma, serum",
        unit="μmol/L",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury", "Nickel"],
        pmid_references=["22023223"],
    ),

    # DNA oxidation markers
    "8_OHdG": OxidativeStressMarker(
        id="8-OHdG",
        name="8-Hydroxy-2'-deoxyguanosine",
        abbreviation="8-OHdG",
        marker_type=OSMarkerType.DNA_OXIDATION,
        detection_method="LC-MS/MS, ELISA",
        sample_matrix="urine, serum, tissue",
        unit="ng/mg creatinine",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Nickel", "Iron", "Copper"],
        pmid_references=["19825587", "22023223", "15626646"],
    ),
    "8_OHGua": OxidativeStressMarker(
        id="8-OH-Gua",
        name="8-Hydroxyguanine",
        abbreviation="8-OH-Gua",
        marker_type=OSMarkerType.DNA_OXIDATION,
        detection_method="LC-MS/MS, GC-MS",
        sample_matrix="urine, plasma",
        unit="ng/mg creatinine",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium", "Chromium"],
        pmid_references=["15626646"],
    ),
    "5_OHdC": OxidativeStressMarker(
        id="5-OH-dC",
        name="5-Hydroxy-2'-deoxycytidine",
        abbreviation="5-OH-dC",
        marker_type=OSMarkerType.DNA_OXIDATION,
        detection_method="LC-MS/MS",
        sample_matrix="urine",
        unit="ng/mg creatinine",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium"],
        pmid_references=["19825587"],
    ),

    # Protein oxidation markers
    "protein_carbonyls": OxidativeStressMarker(
        id="PC",
        name="Protein Carbonyls",
        abbreviation="PC",
        marker_type=OSMarkerType.PROTEIN_OXIDATION,
        detection_method="spectrophotometric (DNPH), ELISA",
        sample_matrix="plasma, serum, tissue",
        unit="nmol/mg protein",
        metal_induced=True,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury", "Nickel"],
        pmid_references=["22023223", "24577282"],
    ),
    "3_nitrotyrosine": OxidativeStressMarker(
        id="3-NT",
        name="3-Nitrotyrosine",
        abbreviation="3-NT",
        marker_type=OSMarkerType.PROTEIN_OXIDATION,
        detection_method="LC-MS/MS, ELISA",
        sample_matrix="plasma, urine",
        unit="ng/mL",
        metal_induced=True,
        associated_metals=["Arsenic", "Chromium", "Nickel"],
        pmid_references=["15626646"],
    ),
    "AOPP": OxidativeStressMarker(
        id="AOPP",
        name="Advanced Oxidation Protein Products",
        abbreviation="AOPP",
        marker_type=OSMarkerType.PROTEIN_OXIDATION,
        detection_method="spectrophotometric",
        sample_matrix="plasma, serum",
        unit="μmol/L",
        metal_induced=True,
        associated_metals=["Cadmium", "Lead", "Mercury"],
        pmid_references=["22023223"],
    ),

    # Antioxidant enzymes
    "SOD": OxidativeStressMarker(
        id="SOD",
        name="Superoxide Dismutase",
        abbreviation="SOD",
        marker_type=OSMarkerType.ANTIOXIDANT_ENZYME,
        detection_method="spectrophotometric, ELISA",
        sample_matrix="erythrocytes, plasma, tissue",
        unit="U/mg protein",
        metal_induced=False,  # Usually decreased under metal stress
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury", "Nickel"],
        pmid_references=["19825587", "22023223"],
    ),
    "catalase": OxidativeStressMarker(
        id="CAT",
        name="Catalase",
        abbreviation="CAT",
        marker_type=OSMarkerType.ANTIOXIDANT_ENZYME,
        detection_method="spectrophotometric",
        sample_matrix="erythrocytes, plasma, tissue",
        unit="U/mg protein",
        metal_induced=False,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Nickel"],
        pmid_references=["22023223"],
    ),
    "GPx": OxidativeStressMarker(
        id="GPx",
        name="Glutathione Peroxidase",
        abbreviation="GPx",
        marker_type=OSMarkerType.ANTIOXIDANT_ENZYME,
        detection_method="spectrophotometric",
        sample_matrix="erythrocytes, plasma",
        unit="U/mg protein",
        metal_induced=False,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury", "Selenium"],
        pmid_references=["19825587", "22023223"],
    ),
    "GR": OxidativeStressMarker(
        id="GR",
        name="Glutathione Reductase",
        abbreviation="GR",
        marker_type=OSMarkerType.ANTIOXIDANT_ENZYME,
        detection_method="spectrophotometric",
        sample_matrix="erythrocytes, plasma",
        unit="U/L",
        metal_induced=False,
        associated_metals=["Arsenic", "Cadmium", "Mercury"],
        pmid_references=["22023223"],
    ),

    # Glutathione metabolism
    "GSH": OxidativeStressMarker(
        id="GSH",
        name="Reduced Glutathione",
        abbreviation="GSH",
        marker_type=OSMarkerType.GSH_METABOLISM,
        detection_method="HPLC, spectrophotometric",
        sample_matrix="erythrocytes, plasma, tissue",
        unit="μmol/L",
        metal_induced=False,  # Depleted under metal stress
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury", "Nickel"],
        pmid_references=["19825587", "22023223"],
    ),
    "GSSG": OxidativeStressMarker(
        id="GSSG",
        name="Oxidized Glutathione",
        abbreviation="GSSG",
        marker_type=OSMarkerType.GSH_METABOLISM,
        detection_method="HPLC, spectrophotometric",
        sample_matrix="erythrocytes, plasma",
        unit="μmol/L",
        metal_induced=True,  # Increased under metal stress
        associated_metals=["Arsenic", "Cadmium", "Mercury", "Nickel"],
        pmid_references=["22023223"],
    ),
    "GSH_GSSG_ratio": OxidativeStressMarker(
        id="GSH/GSSG",
        name="GSH/GSSG Ratio",
        abbreviation="GSH/GSSG",
        marker_type=OSMarkerType.GSH_METABOLISM,
        detection_method="calculated",
        sample_matrix="erythrocytes, plasma",
        unit="ratio",
        metal_induced=False,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury"],
        pmid_references=["19825587"],
    ),

    # Other redox markers
    "TRX": OxidativeStressMarker(
        id="TRX",
        name="Thioredoxin",
        abbreviation="TRX",
        marker_type=OSMarkerType.ANTIOXIDANT_ENZYME,
        detection_method="ELISA, Western blot",
        sample_matrix="plasma, tissue",
        unit="ng/mL",
        metal_induced=False,
        associated_metals=["Arsenic", "Cadmium", "Nickel"],
        pmid_references=["15626646"],
    ),
    "total_antioxidant_capacity": OxidativeStressMarker(
        id="TAC",
        name="Total Antioxidant Capacity",
        abbreviation="TAC",
        marker_type=OSMarkerType.ANTIOXIDANT_ENZYME,
        detection_method="FRAP, ORAC, TEAC",
        sample_matrix="plasma, serum",
        unit="μmol Trolox eq/L",
        metal_induced=False,
        associated_metals=["Arsenic", "Cadmium", "Chromium", "Lead", "Mercury"],
        pmid_references=["22023223"],
    ),
}


# Metal-specific ROS generation pathways
METAL_ROS_PATHWAYS: dict[str, list[ROSGenerationPathway]] = {
    "Arsenic": [
        ROSGenerationPathway(
            metal="Arsenic",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="MMA(III) redox cycling via mitochondrial disruption",
            cellular_compartment="mitochondria",
            enzyme_involved="Complex I/III",
            requires_redox_cycling=True,
            metal_specific_notes="Trivalent metabolites MMA(III) and DMA(III) are most potent ROS inducers",
        ),
        ROSGenerationPathway(
            metal="Arsenic",
            primary_ros=ReactiveSpecies.HYDROGEN_PEROXIDE,
            generation_mechanism="Dismutation of superoxide, inhibition of antioxidant enzymes",
            cellular_compartment="cytosol",
            enzyme_involved="SOD",
            requires_redox_cycling=True,
            metal_specific_notes="Arsenic inhibits catalase and GPx, leading to H2O2 accumulation",
        ),
        ROSGenerationPathway(
            metal="Arsenic",
            primary_ros=ReactiveSpecies.HYDROXYL_RADICAL,
            generation_mechanism="Fenton-like reaction with released iron",
            cellular_compartment="mitochondria, cytosol",
            requires_fenton=True,
            metal_specific_notes="Arsenic releases iron from ferritin, amplifying Fenton chemistry",
        ),
    ],
    "Cadmium": [
        ROSGenerationPathway(
            metal="Cadmium",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="Mitochondrial dysfunction, Complex III inhibition",
            cellular_compartment="mitochondria",
            enzyme_involved="Complex III",
            requires_redox_cycling=False,
            metal_specific_notes="Cadmium displaces iron and copper from proteins, triggering redox activity",
        ),
        ROSGenerationPathway(
            metal="Cadmium",
            primary_ros=ReactiveSpecies.NITRIC_OXIDE,
            generation_mechanism="iNOS induction, inflammation",
            cellular_compartment="cytosol",
            enzyme_involved="iNOS",
            requires_redox_cycling=False,
            metal_specific_notes="Chronic cadmium exposure induces inflammatory NO production",
        ),
        ROSGenerationPathway(
            metal="Cadmium",
            primary_ros=ReactiveSpecies.PEROXYNITRITE,
            generation_mechanism="Reaction of NO with superoxide",
            cellular_compartment="mitochondria, cytosol",
            requires_redox_cycling=True,
            metal_specific_notes="ONOO- formation when both NO and O2•- are elevated",
        ),
    ],
    "Chromium": [
        ROSGenerationPathway(
            metal="Chromium",
            primary_ros=ReactiveSpecies.HYDROXYL_RADICAL,
            generation_mechanism="Cr(V)-mediated Fenton reaction",
            cellular_compartment="cytosol, nucleus",
            requires_fenton=True,
            metal_specific_notes="Cr(VI) reduction to Cr(V) generates radicals in redox cycle",
        ),
        ROSGenerationPathway(
            metal="Chromium",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="NADPH oxidase activation, mitochondrial dysfunction",
            cellular_compartment="cytosol, mitochondria",
            enzyme_involved="NADPH oxidase, Complex I",
            requires_redox_cycling=True,
            metal_specific_notes="Cr(V) can react with H2O2 to regenerate Cr(VI) and •OH",
        ),
        ROSGenerationPathway(
            metal="Chromium",
            primary_ros=ReactiveSpecies.SINGLET_OXYGEN,
            generation_mechanism="Photoactivation of Cr-DNA complexes",
            cellular_compartment="nucleus",
            requires_redox_cycling=False,
            metal_specific_notes="UV light can activate Cr-DNA adducts to generate 1O2",
        ),
    ],
    "Nickel": [
        ROSGenerationPathway(
            metal="Nickel",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="Hypoxia-mimic response, HIF-1α stabilization",
            cellular_compartment="mitochondria",
            enzyme_involved="Complex I/III",
            requires_redox_cycling=False,
            metal_specific_notes="Ni mimics hypoxia, disrupting oxygen sensing and increasing ROS",
        ),
        ROSGenerationPathway(
            metal="Nickel",
            primary_ros=ReactiveSpecies.HYDROXYL_RADICAL,
            generation_mechanism="Fenton reaction with Ni(II) bound to DNA",
            cellular_compartment="nucleus",
            requires_fenton=True,
            metal_specific_notes="Ni binds to DNA minor groove, positioning •OH generation at DNA",
        ),
    ],
    "Lead": [
        ROSGenerationPathway(
            metal="Lead",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="Mitochondrial dysfunction, delta-aminolevulinic acid oxidation",
            cellular_compartment="mitochondria, cytosol",
            enzyme_involved="ALAD (inhibited), delta-ALA",
            requires_redox_cycling=True,
            metal_specific_notes="Delta-ALA accumulates and undergoes auto-oxidation, generating ROS",
        ),
        ROSGenerationPathway(
            metal="Lead",
            primary_ros=ReactiveSpecies.HYDROGEN_PEROXIDE,
            generation_mechanism="Inhibition of antioxidant enzymes (SOD, catalase, GPx)",
            cellular_compartment="cytosol",
            requires_redox_cycling=False,
            metal_specific_notes="Lead directly binds to sulfhydryl groups in antioxidant enzymes",
        ),
    ],
    "Mercury": [
        ROSGenerationPathway(
            metal="Mercury",
            primary_ros=ReactiveSpecies.HYDROGEN_PEROXIDE,
            generation_mechanism="Thiol depletion, glutathione oxidation",
            cellular_compartment="cytosol, mitochondria",
            requires_redox_cycling=True,
            metal_specific_notes="Hg has high affinity for thiols, depleting GSH and protein thiols",
        ),
        ROSGenerationPathway(
            metal="Mercury",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="Mitochondrial dysfunction, thioredoxin system inhibition",
            cellular_compartment="mitochondria",
            enzyme_involved="Thioredoxin reductase",
            requires_redox_cycling=True,
            metal_specific_notes="Hg inhibits thioredoxin reductase, disrupting redox balance",
        ),
    ],
    "Cobalt": [
        ROSGenerationPathway(
            metal="Cobalt",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="Hypoxia mimicry, HIF-1α stabilization",
            cellular_compartment="mitochondria, cytosol",
            enzyme_involved="Prolyl hydroxylases (inhibited)",
            requires_redox_cycling=False,
            metal_specific_notes="Co(II) stabilizes HIF-1α by inhibiting PHDs, mimicking hypoxia",
        ),
        ROSGenerationPathway(
            metal="Cobalt",
            primary_ros=ReactiveSpecies.HYDROXYL_RADICAL,
            generation_mechanism="Fenton-like reaction with ascorbate",
            cellular_compartment="cytosol",
            requires_fenton=True,
            metal_specific_notes="Co can participate in ascorbate-driven Fenton chemistry",
        ),
    ],
    "Iron": [
        ROSGenerationPathway(
            metal="Iron",
            primary_ros=ReactiveSpecies.HYDROXYL_RADICAL,
            generation_mechanism="Classic Fenton reaction",
            cellular_compartment="cytosol, mitochondria, lysosomes",
            requires_fenton=True,
            metal_specific_notes="Fe(II) + H2O2 → Fe(III) + •OH + OH- is the canonical Fenton reaction",
        ),
        ROSGenerationPathway(
            metal="Iron",
            primary_ros=ReactiveSpecies.SUPEROXIDE,
            generation_mechanism="Iron-sulfur cluster destabilization",
            cellular_compartment="mitochondria, cytosol",
            enzyme_involved="Various Fe-S enzymes",
            requires_redox_cycling=True,
            metal_specific_notes="Excess iron destabilizes Fe-S clusters, releasing free iron and ROS",
        ),
    ],
    "Copper": [
        ROSGenerationPathway(
            metal="Copper",
            primary_ros=ReactiveSpecies.HYDROXYL_RADICAL,
            generation_mechanism="Fenton-like reaction",
            cellular_compartment="cytosol",
            requires_fenton=True,
            metal_specific_notes="Cu(I) can drive Fenton chemistry similarly to Fe(II)",
        ),
    ],
}


# Metal-induced antioxidant response patterns
ANTIOXIDANT_RESPONSE_PATTERNS: dict[str, dict[str, str]] = {
    "Arsenic": {
        "SOD": "decreased",
        "CAT": "decreased",
        "GPx": "decreased",
        "GSH": "depleted",
        "phase": "chronic",
    },
    "Cadmium": {
        "SOD": "variable",
        "CAT": "decreased",
        "GPx": "decreased",
        "GSH": "depleted",
        "phase": "chronic",
    },
    "Chromium": {
        "SOD": "increased_then_decreased",
        "CAT": "decreased",
        "GPx": "decreased",
        "GSH": "depleted",
        "phase": "acute_then_chronic",
    },
    "Nickel": {
        "SOD": "increased",
        "CAT": "variable",
        "GPx": "decreased",
        "GSH": "depleted",
        "phase": "variable",
    },
    "Lead": {
        "SOD": "decreased",
        "CAT": "decreased",
        "GPx": "decreased",
        "GSH": "depleted",
        "phase": "dose_dependent",
    },
    "Mercury": {
        "SOD": "variable",
        "CAT": "decreased",
        "GPx": "decreased",
        "GSH": "depleted",
        "phase": "acute",
    },
}


# ── Functions ─────────────────────────────────────────────────────────────


def get_os_marker(marker_id: str) -> OxidativeStressMarker | None:
    """Get oxidative stress marker by ID.

    Args:
        marker_id: Marker identifier (e.g., "MDA", "8-OHdG")

    Returns:
        OxidativeStressMarker if found, None otherwise
    """
    return OXIDATIVE_STRESS_MARKERS.get(marker_id)


def get_markers_by_metal(metal: str) -> list[OxidativeStressMarker]:
    """Get all oxidative stress markers associated with a metal.

    Args:
        metal: Metal name (e.g., "Arsenic", "Cadmium")

    Returns:
        List of associated markers
    """
    markers = []
    for marker in OXIDATIVE_STRESS_MARKERS.values():
        if any(metal.lower() in m.lower() for m in marker.associated_metals):
            markers.append(marker)
    return markers


def get_markers_by_type(marker_type: OSMarkerType) -> list[OxidativeStressMarker]:
    """Get markers by type category.

    Args:
        marker_type: Type of marker

    Returns:
        List of markers of that type
    """
    return [m for m in OXIDATIVE_STRESS_MARKERS.values() if m.marker_type == marker_type]


def get_ros_pathways(metal: str) -> list[ROSGenerationPathway]:
    """Get ROS generation pathways for a metal.

    Args:
        metal: Metal name

    Returns:
        List of ROS generation pathways
    """
    return METAL_ROS_PATHWAYS.get(metal, [])


def get_primary_ros_species(metal: str) -> list[ReactiveSpecies]:
    """Get primary ROS species generated by a metal.

    Args:
        metal: Metal name

    Returns:
        List of primary reactive species
    """
    pathways = get_ros_pathways(metal)
    species = {p.primary_ros for p in pathways}
    return list(species)


def get_antioxidant_response(metal: str) -> dict[str, str]:
    """Get expected antioxidant enzyme response pattern for a metal.

    Args:
        metal: Metal name

    Returns:
        Dictionary with expected changes ("increased", "decreased", "variable")
    """
    return ANTIOXIDANT_RESPONSE_PATTERNS.get(metal, {})


def requires_fenton_chemistry(metal: str) -> bool:
    """Check if metal generates ROS via Fenton chemistry.

    Args:
        metal: Metal name

    Returns:
        True if metal participates in Fenton reactions
    """
    pathways = get_ros_pathways(metal)
    return any(p.requires_fenton for p in pathways)


def calculate_os_risk_score(
    metal: str,
    exposure_duration_days: int,
    biomarker_levels: dict[str, float],
) -> dict[str, Any]:
    """Calculate an oxidative stress risk score.

    Combines metal-specific ROS generation potential with biomarker levels.

    Args:
        metal: Metal name
        exposure_duration_days: Duration of exposure
        biomarker_levels: Dict of {marker_id: level}

    Returns:
        Risk assessment dictionary
    """
    # Get metal ROS potential
    pathways = get_ros_pathways(metal)
    ros_potential = len(pathways)
    fenton_capable = requires_fenton_chemistry(metal)

    # Calculate biomarker severity
    severity_score = 0.0
    elevated_markers = []

    for marker_id, level in biomarker_levels.items():
        marker = get_os_marker(marker_id)
        if marker and marker.reference_range_high:
            if level > marker.reference_range_high:
                ratio = level / marker.reference_range_high
                severity_score += min(ratio - 1.0, 5.0)  # Cap at 5x
                elevated_markers.append(marker_id)

    # Calculate risk score
    # Base from metal potential (0-10), modified by biomarkers (0-10)
    metal_score = min(ros_potential * 2, 10)
    if fenton_capable:
        metal_score += 2

    biomarker_score = min(severity_score, 10)

    # Combined score weighted by exposure duration
    duration_factor = min(exposure_duration_days / 365, 2.0)  # Max 2x after 2 years
    total_score = (metal_score * 0.4 + biomarker_score * 0.6) * duration_factor

    return {
        "metal": metal,
        "ros_generation_pathways": ros_potential,
        "fenton_chemistry": fenton_capable,
        "biomarker_severity": biomarker_score,
        "metal_os_potential": metal_score,
        "elevated_markers": elevated_markers,
        "exposure_duration_factor": duration_factor,
        "total_os_risk_score": round(total_score, 2),
        "risk_category": (
            "high" if total_score > 7
            else "moderate" if total_score > 4
            else "low"
        ),
    }


def build_os_marker_graph(metal: str | None = None) -> KnowledgeGraph:
    """Build a knowledge graph of oxidative stress markers.

    Args:
        metal: Optional metal to filter by

    Returns:
        KnowledgeGraph with markers as nodes
    """
    if metal:
        markers = get_markers_by_metal(metal)
    else:
        markers = list(OXIDATIVE_STRESS_MARKERS.values())

    nodes = []
    edges = []

    # Create marker nodes
    for marker in markers:
        node = Node(
            id=marker.id,
            label=marker.name,
            type=NodeType.METABOLITE,  # Biochemical marker
            detail=f"{marker.marker_type.value} marker. Detection: {marker.detection_method}. "
                   f"Sample: {marker.sample_matrix}.",
            group="Oxidative_Stress_Marker",
            source_db="OS_Biomarker_Database",
            evidence=f"Associated with {', '.join(marker.associated_metals)}",
        )
        nodes.append(node)

    # Create metal nodes and connect to markers
    if metal:
        metal_node = Node(
            id=metal.replace(" ", "_").replace("(", "").replace(")", ""),
            label=metal,
            type=NodeType.CARCINOGEN,
            group="Heavy_Metal",
        )
        nodes.append(metal_node)

        for marker in markers:
            edges.append(Edge(
                source=metal_node.id,
                target=marker.id,
                type=EdgeType.INDUCES,
                evidence=f"Metal exposure induces {marker.marker_type.value}",
            ))

    # Add ROS nodes and connect pathways
    pathways = get_ros_pathways(metal) if metal else []
    for pathway in pathways:
        ros_id = f"ROS_{pathway.primary_ros.value}"
        ros_node = next((n for n in nodes if n.id == ros_id), None)

        if not ros_node:
            ros_node = Node(
                id=ros_id,
                label=pathway.primary_ros.value.replace("_", " ").title(),
                type=NodeType.METABOLITE,
                group="Reactive_Species",
                detail=f"Generated via {pathway.generation_mechanism}",
            )
            nodes.append(ros_node)

        # Connect metal to ROS
        if metal:
            metal_node_id = metal.replace(" ", "_").replace("(", "").replace(")", "")
            edges.append(Edge(
                source=metal_node_id,
                target=ros_id,
                type=EdgeType.INDUCES,
                evidence=pathway.generation_mechanism,
            ))

    return KnowledgeGraph(nodes=nodes, edges=edges)


def get_recommended_os_panel(metal: str, exposure_route: str | None = None) -> list[str]:
    """Get recommended oxidative stress biomarker panel for a metal.

    Args:
        metal: Metal name
        exposure_route: Optional route (inhalation, ingestion, dermal)

    Returns:
        List of recommended marker IDs
    """
    base_markers = []

    # All metal exposures should include core markers
    core = ["8-OHdG", "MDA", "SOD", "GSH", "GSSG"]

    metal_specific = {
        "Arsenic": ["4-HNE", "8-iso-PGF2α", "protein_carbonyls", "GPx", "GR"],
        "Cadmium": ["3-nitrotyrosine", "AOPP", "8-iso-PGF2α", "CAT"],
        "Chromium": ["4-HNE", "8-OH-Gua", "protein_carbonyls"],
        "Nickel": ["3-nitrotyrosine", "8-iso-PGF2α", "TRX"],
        "Lead": ["MDA", "protein_carbonyls", "AOPP", "CAT"],
        "Mercury": ["GSH", "GSSG", "GPx", "TRX"],
    }

    base_markers = core + metal_specific.get(metal, [])

    # Add route-specific markers
    if exposure_route == "inhalation":
        base_markers.append("8-iso-PGF2α")  # Good for lung oxidative stress

    return list(dict.fromkeys(base_markers))  # Remove duplicates, preserve order


def integrate_os_into_metabolism_chain(
    chain_nodes: list[Node],
    chain_edges: list[Edge],
    metal: str,
) -> tuple[list[Node], list[Edge]]:
    """Integrate oxidative stress nodes into a metabolism chain.

    Args:
        chain_nodes: Existing metabolism chain nodes
        chain_edges: Existing metabolism chain edges
        metal: Metal to integrate ROS pathways for

    Returns:
        Updated nodes and edges with OS integration
    """
    nodes = list(chain_nodes)
    edges = list(chain_edges)

    # Add ROS generation nodes
    pathways = get_ros_pathways(metal)
    added_ros = set()

    for pathway in pathways:
        ros_id = f"ROS_{pathway.primary_ros.value}"
        if ros_id not in added_ros:
            ros_node = Node(
                id=ros_id,
                label=pathway.primary_ros.value.replace("_", " ").title(),
                type=NodeType.METABOLITE,
                group="Reactive_Species",
                reactivity="Very High",
                detail=pathway.generation_mechanism,
            )
            nodes.append(ros_node)
            added_ros.add(ros_id)

        # Find enzyme node in chain that could generate this ROS
        for node in chain_nodes:
            if node.type == NodeType.ENZYME and pathway.enzyme_involved:
                if pathway.enzyme_involved.lower() in node.id.lower() or \
                   pathway.enzyme_involved.lower() in node.label.lower():
                    # Add edge from enzyme to ROS
                    edges.append(Edge(
                        source=node.id,
                        target=ros_id,
                        type=EdgeType.INDUCES,
                        carcinogen=metal.replace(" ", "_").replace("(", "").replace(")", ""),
                        evidence=pathway.generation_mechanism[:100],
                    ))

    # Add DNA damage node if 8-OHdG is relevant
    if any(metal.lower() in m.lower() for m in ["Arsenic", "Cadmium", "Chromium", "Nickel"]):
        # Check if 8-OHdG already exists
        ohDG_exists = any(n.id == "8_OHdG" for n in nodes)
        if not ohDG_exists:
            ohDG_node = Node(
                id="8_OHdG",
                label="8-OHdG",
                type=NodeType.DNA_ADDUCT,
                group="DNA_Oxidative_Damage",
            )
            nodes.append(ohDG_node)

            # Connect ROS to 8-OHdG
            for ros_id in added_ros:
                if "hydroxyl" in ros_id.lower() or "superoxide" in ros_id.lower():
                    edges.append(Edge(
                        source=ros_id,
                        target="8_OHdG",
                        type=EdgeType.FORMS_ADDUCT,
                        carcinogen=metal.replace(" ", "_").replace("(", "").replace(")", ""),
                        evidence="ROS-induced DNA oxidation",
                    ))

    return nodes, edges
