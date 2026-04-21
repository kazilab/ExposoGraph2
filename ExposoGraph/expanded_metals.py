"""Expanded metal coverage: Lead, Mercury, Cobalt, Antimony.

This module adds four additional heavy metals with established or probable
carcinogenicity (IARC Groups 2A and 2B) and their metabolic pathways.

Metals included:
- Lead (Pb) - IARC Group 2A (probable carcinogen)
- Mercury (Hg) - IARC Group 2B (possible carcinogen), especially methylmercury
- Cobalt (Co) - IARC Group 2B, includes hard metal disease mechanisms
- Antimony (Sb) - Often co-exposed with arsenic, respiratory toxicity
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .db_clients.iarc import IARCGroup
from .models import Edge, EdgeType, KnowledgeGraph, Node, NodeType


@dataclass
class MetalProfile:
    """Comprehensive profile for a metal carcinogen."""

    symbol: str
    name: str
    atomic_number: int
    iarc_group: IARCGroup
    primary_exposure_routes: list[str]
    target_organs: list[str]
    half_life_days: float
    key_health_effects: list[str]
    metabolism_genes: list[str]
    transport_genes: list[str]
    detoxification_pathway: str
    oxidative_stress_mechanism: str
    epigenetic_effects: list[str]
    pmid_references: list[str] = field(default_factory=list)


# ── Expanded Metal Database ────────────────────────────────────────────────

EXPANDED_METALS: dict[str, MetalProfile] = {
    "Lead": MetalProfile(
        symbol="Pb",
        name="Lead",
        atomic_number=82,
        iarc_group=IARCGroup.GROUP_2A,
        primary_exposure_routes=[
            "inhalation (occupational dust/fumes)",
            "ingestion (contaminated water, paint chips)",
            "dermal (limited)",
        ],
        target_organs=[
            "kidney",
            "brain (neurodevelopmental)",
            "blood (heme synthesis)",
            "bone (storage)",
            "cardiovascular",
        ],
        half_life_days=30 * 365,  # ~30 years in bone
        key_health_effects=[
            "Nephrotoxicity",
            "Neurotoxicity (especially in children)",
            "Hypertension",
            "Reproductive toxicity",
            "Renal cell carcinoma",
        ],
        metabolism_genes=[
            "ALAD",  # Delta-aminolevulinic acid dehydratase - major target
            "VDR",   # Vitamin D receptor - influences bone deposition
            "HFE",   # Hemochromatosis gene - iron-lead interactions
            "TF",    # Transferrin - transport
            "ALPL",  # Alkaline phosphatase
        ],
        transport_genes=[
            "SLC11A2",  # DMT1 - primary uptake transporter
            "SLC11A1",  # NRAMP1
            "SLC39A14", # ZIP14
            "ATP7B",    # Wilson disease protein
            "CFTR",     # Cystic fibrosis transmembrane conductance
        ],
        detoxification_pathway="Sequestration in bone as lead phosphate; chelation therapy for acute",
        oxidative_stress_mechanism="Inhibition of ALAD causes delta-ALA accumulation which auto-oxidizes generating ROS; disruption of heme synthesis",
        epigenetic_effects=[
            "DNA methylation changes in tumor suppressor genes",
            "Histone modification in kidney cells",
            "miRNA dysregulation (miR-21 upregulation)",
        ],
        pmid_references=["22935954", "21768835", "22023223"],
    ),

    "Mercury": MetalProfile(
        symbol="Hg",
        name="Mercury",
        atomic_number=80,
        iarc_group=IARCGroup.GROUP_2B,
        primary_exposure_routes=[
            "ingestion (methylmercury in fish)",
            "inhalation (elemental mercury vapor)",
            "dermal (elemental)",
        ],
        target_organs=[
            "brain (neurotoxicity)",
            "kidney",
            "immune system",
            "cardiovascular",
        ],
        half_life_days=70,  # 45-70 days for methylmercury
        key_health_effects=[
            "Neurotoxicity (Minamata disease)",
            "Nephrotic syndrome",
            "Autoimmune effects",
            "Cardiovascular disease",
        ],
        metabolism_genes=[
            "GSTM1",  # Glutathione conjugation
            "GSTT1",  # Glutathione conjugation
            "MT1A",   # Metallothionein sequestration
            "MT2A",   # Metallothionein sequestration
            "SEPP1",  # Selenoprotein P - protective
            "GPX1",   # Glutathione peroxidase 1
        ],
        transport_genes=[
            "SLC7A5",  # LAT1 - amino acid mimicry transport
            "SLC7A8",  # LAT2
            "SLC3A2",  # CD98
            "LRP2",    # Megalin - renal reabsorption
            "OAT1",    # Organic anion transporter
        ],
        detoxification_pathway="GSH conjugation and biliary excretion; MT sequestration; selenium antagonism",
        oxidative_stress_mechanism="High affinity for thiols depletes GSH; inhibits thioredoxin reductase and GPx",
        epigenetic_effects=[
            "Global DNA hypomethylation",
            "Altered histone acetylation in neuronal cells",
            "miRNA-146a upregulation in kidney",
        ],
        pmid_references=["22023223", "24577282"],
    ),

    "Cobalt": MetalProfile(
        symbol="Co",
        name="Cobalt",
        atomic_number=27,
        iarc_group=IARCGroup.GROUP_2B,
        primary_exposure_routes=[
            "inhalation (hard metal dust, carbide tools)",
            "ingestion (vitamin B12, contaminated water)",
            "implanted devices (hip prostheses)",
        ],
        target_organs=[
            "lung",
            "thyroid",
            "heart",
            "blood (polycythemia)",
            "bones (prosthetic implants)",
        ],
        half_life_days=10,  # 5-15 days
        key_health_effects=[
            "Hard metal disease (interstitial lung disease)",
            "Asthma",
            "Polycythemia vera",
            "Cardiomyopathy",
            "Thyroid dysfunction",
        ],
        metabolism_genes=[
            "HIF1A",   # Hypoxia inducible factor 1A - stabilized by Co
            "EGLN1",   # PHD2 - inhibited by Co
            "EPAS1",   # HIF-2α
            "VHL",     # Von Hippel-Lindau - HIF regulation
            "EPO",     # Erythropoietin - upregulated by Co
        ],
        transport_genes=[
            "SLC11A2",  # DMT1 - shared with iron
            "SLC40A1",  # Ferroportin - export
            "TF",       # Transferrin
            "TFR2",     # Transferrin receptor 2
            "FLVCR1",   # Heme export
        ],
        detoxification_pathway="Urinary excretion (primary); competes with iron for transporters",
        oxidative_stress_mechanism="Hypoxia mimicry stabilizes HIF-1α; ascorbate-driven Fenton chemistry",
        epigenetic_effects=[
            "HIF-mediated transcriptional reprogramming",
            "Histone demethylase inhibition",
            "Epigenetic silencing of tumor suppressors",
        ],
        pmid_references=["24577282", "22935954"],
    ),

    "Antimony": MetalProfile(
        symbol="Sb",
        name="Antimony",
        atomic_number=51,
        iarc_group=IARCGroup.GROUP_2B,
        primary_exposure_routes=[
            "inhalation (mining, smelting dust)",
            "ingestion (contaminated water)",
            "dermal (limited)",
        ],
        target_organs=[
            "lung",
            "liver",
            "cardiovascular",
        ],
        half_life_days=2,  # Rapid excretion
        key_health_effects=[
            "Pneumoconiosis (antimony lung disease)",
            "Cardiotoxicity",
            "Hepatotoxicity",
            "Often co-toxic with arsenic",
        ],
        metabolism_genes=[
            "GSTO1",   # Similar to arsenic metabolism
            "GSTM1",   # Glutathione conjugation
            "GSTT1",   # Glutathione conjugation
        ],
        transport_genes=[
            "AQP7",    # Aquaporin 7 - shared with arsenic
            "AQP9",    # Aquaporin 9
            "ABCC2",   # MRP2 - efflux
            "SLC11A2", # DMT1
        ],
        detoxification_pathway="Conjugation with GSH and biliary excretion (similar to arsenic)",
        oxidative_stress_mechanism="Pentavalent reduction generates ROS; thiol depletion",
        epigenetic_effects=[
            "Limited data on epigenetic effects",
            "May mimic arsenic in DNA methylation effects",
        ],
        pmid_references=["22935954"],
    ),
}


# Expanded metal node definitions for knowledge graph
EXPANDED_METAL_NODES: list[dict[str, Any]] = [
    {
        "id": "Pb",
        "label": "Lead",
        "type": "Carcinogen",
        "group": "Heavy_Metal",
        "iarc": "Group 2A",
        "exposure": "Contaminated water, occupational dust/fumes, legacy paint",
        "detail": "IARC Group 2A probable carcinogen. 30-year bone half-life. Target: kidney, brain, blood.",
        "atomic_number": 82,
    },
    {
        "id": "Hg",
        "label": "Mercury",
        "type": "Carcinogen",
        "group": "Heavy_Metal",
        "iarc": "Group 2B",
        "exposure": "Fish consumption (methylmercury), dental amalgams, occupational vapor",
        "detail": "IARC Group 2B possible carcinogen. 70-day half-life. Target: brain, kidney. High thiol affinity.",
        "atomic_number": 80,
    },
    {
        "id": "Co",
        "label": "Cobalt",
        "type": "Carcinogen",
        "group": "Heavy_Metal",
        "iarc": "Group 2B",
        "exposure": "Hard metal tools, carbide manufacturing, hip prostheses",
        "detail": "IARC Group 2B possible carcinogen. Hypoxia mimic. Target: lung, thyroid. 10-day half-life.",
        "atomic_number": 27,
    },
    {
        "id": "Sb",
        "label": "Antimony",
        "type": "Carcinogen",
        "group": "Heavy_Metal",
        "iarc": "Group 2B",
        "exposure": "Mining, smelting, often co-exposed with arsenic",
        "detail": "IARC Group 2B possible carcinogen. Often co-toxic with arsenic. Target: lung, heart, liver.",
        "atomic_number": 51,
    },
]


# Enzyme and transporter nodes for expanded metals
EXPANDED_METAL_ENZYMES: list[dict[str, Any]] = [
    # Lead-specific enzymes
    {
        "id": "ALAD",
        "label": "ALAD",
        "type": "Enzyme",
        "phase": "II",
        "role": "Detoxification",
        "detail": "Delta-aminolevulinic acid dehydratase. Primary lead target. Inhibition causes delta-ALA accumulation and ROS.",
        "tissue": "blood, liver, kidney",
        "group": "Heme_Synthesis",
        "metal_targets": ["Pb"],
    },
    {
        "id": "ALAS1",
        "label": "ALAS1",
        "type": "Enzyme",
        "phase": "I",
        "role": "Activation",
        "detail": "Aminolevulinic acid synthase 1. Rate-limiting in heme synthesis.",
        "tissue": "liver, erythroid",
        "group": "Heme_Synthesis",
    },

    # Mercury-specific
    {
        "id": "SEPP1",
        "label": "SEPP1",
        "type": "Enzyme",
        "phase": "II",
        "role": "Detoxification",
        "detail": "Selenoprotein P. Protective against mercury toxicity via selenium antagonism.",
        "tissue": "liver, plasma",
        "group": "Selenium_Protection",
        "metal_targets": ["Hg"],
    },

    # Cobalt-specific (HIF pathway)
    {
        "id": "EGLN1",
        "label": "EGLN1 (PHD2)",
        "type": "Enzyme",
        "phase": "I",
        "role": "Regulation",
        "detail": "Prolyl hydroxylase domain 2. Inhibited by cobalt, stabilizing HIF-1α (hypoxia mimicry).",
        "tissue": "ubiquitous",
        "group": "HIF_Pathway",
        "metal_targets": ["Co"],
    },
    {
        "id": "HIF1A",
        "label": "HIF-1α",
        "type": "Gene",  # Actually a transcription factor
        "detail": "Hypoxia inducible factor 1 alpha. Stabilized by cobalt, driving transcriptional reprogramming.",
        "tissue": "ubiquitous",
        "group": "HIF_Pathway",
    },
    {
        "id": "GPx",
        "label": "Glutathione peroxidase",
        "type": "Enzyme",
        "phase": "II",
        "role": "Detoxification",
        "detail": "Selenoprotein; mercury inhibits selenol-dependent activity.",
        "tissue": "ubiquitous",
        "group": "Antioxidant",
    },
    {
        "id": "VEGF",
        "label": "VEGFA",
        "type": "Gene",
        "detail": "Angiogenic factor upregulated under HIF-1α (cobalt hypoxia mimicry).",
        "tissue": "ubiquitous",
        "group": "HIF_Pathway",
    },
]


# Transporter nodes for expanded metals
EXPANDED_METAL_TRANSPORTERS: list[dict[str, Any]] = [
    {
        "id": "DMT1_Pb",
        "label": "DMT1 (SLC11A2)",
        "type": "Enzyme",
        "phase": "III",
        "role": "Transport",
        "detail": "Divalent metal transporter 1. Major entry route for lead (shared with iron).",
        "tissue": "intestine, kidney, brain",
        "group": "Metal_Transport",
        "metal_targets": ["Pb", "Cd", "Co"],
    },
    {
        "id": "LAT1_Hg",
        "label": "LAT1 (SLC7A5)",
        "type": "Enzyme",
        "phase": "III",
        "role": "Transport",
        "detail": "L-type amino acid transporter 1. Methylmercury-cysteine complex mimicry entry.",
        "tissue": "brain (BBB), placenta, intestine",
        "group": "Amino_Acid_Transport",
        "metal_targets": ["Hg"],
    },
    {
        "id": "ZIP14_Co",
        "label": "ZIP14 (SLC39A14)",
        "type": "Enzyme",
        "phase": "III",
        "role": "Transport",
        "detail": "Zinc/iron permease 14. Transports non-transferrin bound cobalt and iron.",
        "tissue": "liver, pancreas, heart",
        "group": "Metal_Transport",
        "metal_targets": ["Co", "Fe"],
    },
]


# Metabolite nodes for expanded metals
EXPANDED_METAL_METABOLITES: list[dict[str, Any]] = [
    {
        "id": "delta_ALA",
        "label": "Delta-ALA",
        "type": "Metabolite",
        "reactivity": "High",
        "detail": "Delta-aminolevulinic acid. Accumulates with lead exposure. Auto-oxidizes generating ROS.",
        "associated_metal": "Pb",
    },
    {
        "id": "ALA_portporphyria",
        "label": "ALA Porphyria",
        "type": "Metabolite",
        "reactivity": "High",
        "detail": "Acquired porphyria from lead-induced ALAD inhibition.",
        "associated_metal": "Pb",
    },
    {
        "id": "MeHg_GSH",
        "label": "Methylmercury-GSH",
        "type": "Metabolite",
        "reactivity": "Intermediate",
        "detail": "Methylmercury-glutathione complex. Transported and excreted.",
        "associated_metal": "Hg",
    },
    {
        "id": "Co_plus2",
        "label": "Co(II)",
        "type": "Metabolite",
        "reactivity": "High",
        "detail": "Cobalt(II) ion. Redox active, participates in Fenton-like reactions.",
        "associated_metal": "Co",
    },
    {
        "id": "GSH",
        "label": "Glutathione",
        "type": "Metabolite",
        "reactivity": "High",
        "detail": "Major cellular antioxidant; depleted by lead and mercury binding.",
    },
]


# Edge definitions for expanded metals
EXPANDED_METAL_EDGES: list[dict[str, Any]] = [
    # Lead pathways
    {"source": "Pb", "target": "ALAD", "type": "INHIBITS", "carcinogen": "Pb", "evidence": "Lead binds to sulfhydryl groups in ALAD"},
    {"source": "ALAD", "target": "delta_ALA", "type": "ACCUMULATES", "carcinogen": "Pb", "evidence": "Inhibition causes ALA accumulation"},
    {"source": "delta_ALA", "target": "ROS_metal", "type": "GENERATES", "carcinogen": "Pb", "evidence": "Auto-oxidation of ALA generates ROS"},
    {"source": "DMT1_Pb", "target": "Pb", "type": "TRANSPORTS", "carcinogen": "Pb", "evidence": "DMT1 mediates lead uptake"},
    {"source": "Pb", "target": "GSH", "type": "DEPLETES", "carcinogen": "Pb", "evidence": "Lead binds GSH"},

    # Mercury pathways
    {"source": "LAT1_Hg", "target": "Hg", "type": "TRANSPORTS", "carcinogen": "Hg", "evidence": "Methylmercury-cysteine mimicry"},
    {"source": "Hg", "target": "GSH", "type": "DEPLETES", "carcinogen": "Hg", "evidence": "High affinity thiol binding"},
    {"source": "Hg", "target": "GPx", "type": "INHIBITS", "carcinogen": "Hg", "evidence": "Selenol group binding"},
    {"source": "SEPP1", "target": "Hg", "type": "ANTAGONIZES", "carcinogen": "Hg", "evidence": "Selenium antagonism"},
    {"source": "Hg", "target": "MT1A", "type": "INDUCES", "carcinogen": "Hg", "evidence": "Mercury induces MT expression"},

    # Cobalt pathways
    {"source": "Co", "target": "EGLN1", "type": "INHIBITS", "carcinogen": "Co", "evidence": "PHD2 inhibition by Co(II)"},
    {"source": "EGLN1", "target": "HIF1A", "type": "REGULATES", "carcinogen": "Co", "evidence": "PHD2 normally degrades HIF-1α"},
    {"source": "Co", "target": "HIF1A", "type": "STABILIZES", "carcinogen": "Co", "evidence": "Hypoxia mimicry"},
    {"source": "HIF1A", "target": "VEGF", "type": "INDUCES", "carcinogen": "Co", "evidence": "HIF-1α induces angiogenic factors"},
    {"source": "ZIP14_Co", "target": "Co", "type": "TRANSPORTS", "carcinogen": "Co", "evidence": "Non-transferrin bound uptake"},
    {"source": "Co_plus2", "target": "ROS_metal", "type": "GENERATES", "carcinogen": "Co", "evidence": "Fenton-like chemistry"},

    # Antimony pathways
    {"source": "Sb", "target": "GSTO1", "type": "METABOLIZES", "carcinogen": "Sb", "evidence": "Similar to arsenic metabolism"},
    {"source": "Sb", "target": "As", "type": "CO_EXPOSED", "carcinogen": "Sb", "evidence": "Often co-exposed with arsenic in mining"},
]


# ── Functions ─────────────────────────────────────────────────────────────


def _edge_identity(edge: Edge) -> tuple[str, str, str, str | None, str | None]:
    return (
        edge.source,
        edge.target,
        edge.type.value,
        edge.custom_predicate,
        edge.carcinogen,
    )


def _merge_knowledge_graphs(*graphs: KnowledgeGraph) -> KnowledgeGraph:
    """Merge graphs without duplicating node ids or edge identities."""
    node_by_id: dict[str, Node] = {}
    for graph in graphs:
        for node in graph.nodes:
            node_by_id.setdefault(node.id, node)

    edge_by_identity: dict[tuple[str, str, str, str | None, str | None], Edge] = {}
    for graph in graphs:
        for edge in graph.edges:
            edge_by_identity.setdefault(_edge_identity(edge), edge)

    return KnowledgeGraph(
        nodes=list(node_by_id.values()),
        edges=list(edge_by_identity.values()),
    )


def load_core_metallo_reference_graph() -> KnowledgeGraph:
    """Load the bundled D3 reference map (133 nodes, 166 edges).

    Expanded-metal pathway edges reference nodes from this core graph
    (e.g. ``As``, ``GSTO1``, ``ROS_metal``); merge via
    :func:`merge_expanded_metals_into` or :func:`build_expanded_metal_graph`.
    """
    from .example_graphs import _normalize_reference_graph_ids
    from .exporter import parse_graph_data_js

    path = Path(__file__).resolve().parent / "map" / "graph-data.js"
    return _normalize_reference_graph_ids(parse_graph_data_js(path))


def load_heavy_metal_reference_graph() -> KnowledgeGraph:
    """Load the bundled module-01 heavy-metal reference graph.

    This is the package-native equivalent of the old ``01_heavy_metals``
    standalone JSON import. The returned graph already contains the five IARC
    Group 1 heavy metals together with the canonical shared panel entities and
    their integrated edges.
    """
    return load_core_metallo_reference_graph()


def merge_heavy_metal_reference_into(graph: KnowledgeGraph) -> KnowledgeGraph:
    """Merge the bundled module-01 heavy-metal reference graph into *graph*.

    This mirrors the old ``add_heavy_metals.py`` merge workflow, but operates
    directly on package-native ``KnowledgeGraph`` models.
    """
    return _merge_knowledge_graphs(graph, load_heavy_metal_reference_graph())


def _expanded_metal_node_list() -> list[Node]:
    """All nodes defined for the Pb/Hg/Co/Sb expansion layer."""
    nodes: list[Node] = []

    for metal_data in EXPANDED_METAL_NODES:
        nodes.append(
            Node(
                id=metal_data["id"],
                label=metal_data["label"],
                type=NodeType(metal_data["type"]),
                group=metal_data["group"],
                iarc=metal_data["iarc"],
                exposure=metal_data["exposure"],
                detail=metal_data["detail"],
            )
        )

    for enzyme_data in EXPANDED_METAL_ENZYMES:
        nodes.append(
            Node(
                id=enzyme_data["id"],
                label=enzyme_data["label"],
                type=NodeType(enzyme_data["type"]),
                phase=enzyme_data.get("phase"),
                role=enzyme_data.get("role"),
                detail=enzyme_data["detail"],
                tissue=enzyme_data.get("tissue"),
                group=enzyme_data.get("group"),
            )
        )

    for transporter_data in EXPANDED_METAL_TRANSPORTERS:
        nodes.append(
            Node(
                id=transporter_data["id"],
                label=transporter_data["label"],
                type=NodeType(transporter_data["type"]),
                phase=transporter_data.get("phase"),
                role=transporter_data.get("role"),
                detail=transporter_data["detail"],
                tissue=transporter_data.get("tissue"),
                group=transporter_data.get("group"),
            )
        )

    for metabolite_data in EXPANDED_METAL_METABOLITES:
        nodes.append(
            Node(
                id=metabolite_data["id"],
                label=metabolite_data["label"],
                type=NodeType(metabolite_data["type"]),
                reactivity=metabolite_data.get("reactivity"),
                detail=metabolite_data["detail"],
            )
        )

    return nodes


def _expanded_metal_edges_for_node_ids(node_ids: set[str]) -> list[Edge]:
    """Instantiate expansion edges whose endpoints exist in ``node_ids``."""
    edges: list[Edge] = []
    for edge_data in EXPANDED_METAL_EDGES:
        if edge_data["source"] not in node_ids or edge_data["target"] not in node_ids:
            continue

        etype_str = edge_data["type"]
        try:
            etype = EdgeType(etype_str)
        except ValueError:
            etype = EdgeType.CUSTOM

        edges.append(
            Edge(
                source=edge_data["source"],
                target=edge_data["target"],
                type=etype,
                carcinogen=edge_data.get("carcinogen"),
                evidence=edge_data.get("evidence"),
                custom_predicate=etype_str if etype == EdgeType.CUSTOM else None,
            )
        )
    return edges


def merge_expanded_metals_into(core_graph: KnowledgeGraph) -> KnowledgeGraph:
    """Merge Pb/Hg/Co/Sb nodes and edges into a core reference graph.

    All :data:`EXPANDED_METAL_EDGES` are retained when both endpoints exist
    in the union of ``core_graph`` and the expansion nodes (typically after
    loading :func:`load_core_metallo_reference_graph`).
    """
    expanded_nodes = _expanded_metal_node_list()
    node_by_id = {n.id: n for n in core_graph.nodes}
    for node in expanded_nodes:
        node_by_id.setdefault(node.id, node)

    combined_nodes = list(node_by_id.values())
    node_ids = set(node_by_id)

    expanded_edges = _expanded_metal_edges_for_node_ids(node_ids)
    existing_edges = {(e.source, e.target, e.type.value) for e in core_graph.edges}
    combined_edges = list(core_graph.edges)

    for edge in expanded_edges:
        edge_sig = (edge.source, edge.target, edge.type.value)
        if edge_sig not in existing_edges:
            combined_edges.append(edge)
            existing_edges.add(edge_sig)

    return KnowledgeGraph(nodes=combined_nodes, edges=combined_edges)


def get_metal_profile(metal: str) -> MetalProfile | None:
    """Get comprehensive profile for an expanded metal.

    Args:
        metal: Metal name (e.g., "Lead", "Mercury")

    Returns:
        MetalProfile if found, None otherwise
    """
    return EXPANDED_METALS.get(metal)


def get_all_expanded_metals() -> list[MetalProfile]:
    """Get all expanded metal profiles.

    Returns:
        List of all expanded metal profiles
    """
    return list(EXPANDED_METALS.values())


def get_metals_by_iarc_group(group: IARCGroup) -> list[MetalProfile]:
    """Get metals by IARC classification.

    Args:
        group: IARC classification group

    Returns:
        List of metals in that group
    """
    return [m for m in EXPANDED_METALS.values() if m.iarc_group == group]


def build_expanded_metal_graph(*, include_core: bool = True) -> KnowledgeGraph:
    """Build expanded-metal pathways merged into the core reference graph by default.

    The bundled reference map (133 nodes, 166 edges) includes the five IARC
    Group 1 heavy metals and shared entities such as ``ROS`` and ``GSTO1``.
    Merging ensures all :data:`EXPANDED_METAL_EDGES` resolve.

    Args:
        include_core: If ``True`` (default), merge into
            :func:`load_core_metallo_reference_graph`. If ``False``, return only
            the expansion nodes and edges that are internally consistent
            (several cross-core edges are omitted).
    """
    if include_core:
        return merge_expanded_metals_into(load_core_metallo_reference_graph())

    nodes = _expanded_metal_node_list()
    node_ids = {n.id for n in nodes}
    edges = _expanded_metal_edges_for_node_ids(node_ids)
    return KnowledgeGraph(nodes=nodes, edges=edges)


def get_metal_specific_nodes(metal: str) -> tuple[list[Node], list[Edge]]:
    """Get nodes and edges specific to a single expanded metal.

    Args:
        metal: Metal symbol or name

    Returns:
        Tuple of (nodes, edges) specific to that metal
    """
    full_graph = build_expanded_metal_graph(include_core=True)

    # Normalize metal identifier
    metal_id = metal.replace(" ", "_").replace("(", "").replace(")", "")

    # Find all connected nodes
    connected = {metal_id}

    # Add nodes referenced by edges
    for edge in full_graph.edges:
        if edge.carcinogen == metal_id:
            connected.add(edge.source)
            connected.add(edge.target)
        if edge.source == metal_id or edge.target == metal_id:
            connected.add(edge.source)
            connected.add(edge.target)

    # Filter nodes
    nodes = [n for n in full_graph.nodes if n.id in connected]

    # Filter edges
    edges = [e for e in full_graph.edges
             if e.source in connected and e.target in connected]

    return nodes, edges


def integrate_with_core_metals(core_graph: KnowledgeGraph) -> KnowledgeGraph:
    """Integrate expanded metals with the core ExposoGraph reference graph.

    Args:
        core_graph: Core graph (e.g. bundled 122/144 reference or a subgraph).

    Returns:
        Combined knowledge graph with Pb, Hg, Co, Sb nodes and full edge set.
    """
    return merge_expanded_metals_into(core_graph)


def get_co_exposure_profile(metals: list[str]) -> dict[str, Any]:
    """Analyze co-exposure scenarios for multiple metals.

    Args:
        metals: List of metal names

    Returns:
        Co-exposure analysis
    """
    # Known co-exposure scenarios
    co_exposure_pairs = {
        ("Arsenic", "Antimony"): "Mining/smelting co-exposure, shared metabolism (GSTO1)",
        ("Cadmium", "Lead"): "Tobacco smoke, occupational settings",
        ("Lead", "Mercury"): "Neurotoxic synergism, shared DMT1 transport",
        ("Cobalt", "Tungsten"): "Hard metal disease (CoWC)",
        ("Cadmium", "Nickel"): "Battery manufacturing, welding",
    }

    synergistic_effects = []
    for i, metal1 in enumerate(metals):
        for metal2 in metals[i+1:]:
            key = (metal1, metal2)
            reverse_key = (metal2, metal1)

            if key in co_exposure_pairs:
                synergistic_effects.append({
                    "metals": [metal1, metal2],
                    "scenario": co_exposure_pairs[key],
                })
            elif reverse_key in co_exposure_pairs:
                synergistic_effects.append({
                    "metals": [metal1, metal2],
                    "scenario": co_exposure_pairs[reverse_key],
                })

    # Check for shared pathways
    shared_transporters = set()
    for metal in metals:
        profile = get_metal_profile(metal)
        if profile:
            shared_transporters.update(profile.transport_genes)

    return {
        "metals": metals,
        "co_exposure_scenarios": synergistic_effects,
        "shared_transporters": list(shared_transporters),
        "risk_note": "Co-exposure may increase toxicity through shared pathways",
    }


def get_expanded_metal_os_integration(metal: str) -> dict[str, Any]:
    """Get oxidative stress integration data for an expanded metal.

    Args:
        metal: Metal name

    Returns:
        Integration data for oxidative stress module
    """
    profile = get_metal_profile(metal)
    if not profile:
        return {"error": f"Metal {metal} not found in expanded database"}

    # Cross-reference with oxidative stress module
    from .oxidative_stress import get_markers_by_metal, get_ros_pathways

    ros_pathways = get_ros_pathways(metal)
    os_markers = get_markers_by_metal(metal)

    return {
        "metal": metal,
        "primary_os_mechanism": profile.oxidative_stress_mechanism,
        "ros_pathways": len(ros_pathways),
        "recommended_markers": [m.id for m in os_markers[:5]],
        "epigenetic_effects": profile.epigenetic_effects,
        "detoxification_pathway": profile.detoxification_pathway,
    }
