"""Orchestrator for seeding the knowledge graph from public databases.

Converts KEGG, CTD, and IARC data into :class:`KnowledgeGraph` objects
that can be merged into an existing graph via the :class:`GraphEngine`.
"""

from __future__ import annotations

import logging
from typing import Optional

from .config import GraphMode
from .db_clients.ctd import ChemicalGeneInteraction, CTDClient
from .db_clients.iarc import IARCClassifier
from .db_clients.kegg import KEGGClient
from .grounding import prepare_knowledge_graph
from .models import Edge, EdgeType, KnowledgeGraph, Node, NodeType, ProvenanceRecord, RecordOrigin

logger = logging.getLogger(__name__)


def seed_from_kegg_pathway(
    pathway_id: str,
    *,
    client: Optional[KEGGClient] = None,
    mode: GraphMode | str = GraphMode.EXPLORATORY,
) -> KnowledgeGraph:
    """Build a KnowledgeGraph from a KEGG pathway.

    Creates a Pathway node and Gene nodes for all member genes,
    connected by PATHWAY edges.

    Parameters
    ----------
    pathway_id:
        KEGG pathway identifier, e.g. ``"hsa05204"``.
    client:
        Optional pre-configured :class:`KEGGClient`.
    """
    if client is None:
        client = KEGGClient()

    pathway = client.get_pathway(pathway_id)
    nodes: list[Node] = [
        Node(
            id=pathway.pathway_id,
            label=pathway.name or pathway.pathway_id,
            type=NodeType.PATHWAY,
            provenance=[ProvenanceRecord(source_db="KEGG", record_id=pathway.pathway_id)],
        )
    ]
    edges: list[Edge] = []

    for gene_symbol in pathway.genes:
        gene_id = gene_symbol.upper()
        nodes.append(
            Node(
                id=gene_id,
                label=gene_symbol,
                type=NodeType.GENE,
                provenance=[ProvenanceRecord(source_db="KEGG")],
            )
        )
        edges.append(
            Edge(
                source=gene_id,
                target=pathway.pathway_id,
                type=EdgeType.PATHWAY,
                source_db="KEGG",
            )
        )

    raw_graph = KnowledgeGraph(
        nodes=[node.model_copy(update={"origin": RecordOrigin.SEEDED}) for node in nodes],
        edges=[edge.model_copy(update={"origin": RecordOrigin.SEEDED}) for edge in edges],
    )
    prepared_graph, _warnings = prepare_knowledge_graph(
        raw_graph,
        mode=mode,
        reference_graphs=[("kegg", raw_graph)],
    )
    return prepared_graph


def seed_from_ctd(
    chemical_name: str,
    *,
    client: Optional[CTDClient] = None,
    organism: str = "Homo sapiens",
    mode: GraphMode | str = GraphMode.EXPLORATORY,
) -> KnowledgeGraph:
    """Build a KnowledgeGraph from CTD chemical-gene interactions.

    Creates a Carcinogen node for the chemical and Gene nodes for each
    interacting gene, connected by ACTIVATES or DETOXIFIES edges
    based on interaction text heuristics.

    Parameters
    ----------
    chemical_name:
        Chemical name to query (e.g. ``"Benzo(a)pyrene"``).
    client:
        Optional pre-configured :class:`CTDClient`.
    organism:
        Organism filter. Defaults to ``"Homo sapiens"``.
    """
    if client is None:
        client = CTDClient()

    interactions = client.get_chemical_gene_interactions(
        chemical_name, organism=organism,
    )

    chem_id = Node.generate_id(chemical_name)
    nodes: list[Node] = [
        Node(
            id=chem_id,
            label=chemical_name,
            type=NodeType.CARCINOGEN,
            provenance=[ProvenanceRecord(source_db="CTD")],
        )
    ]
    edges: list[Edge] = []
    seen_genes: set[str] = set()

    for ixn in interactions:
        gene_id = ixn.gene_symbol.upper()
        if gene_id not in seen_genes:
            seen_genes.add(gene_id)
            pmid = ixn.pubmed_ids[0] if ixn.pubmed_ids else None
            nodes.append(
                Node(
                    id=gene_id,
                    label=ixn.gene_symbol,
                    type=NodeType.GENE,
                    provenance=[
                        ProvenanceRecord(
                            source_db="CTD",
                            record_id=ixn.gene_id,
                            pmid=pmid,
                        )
                    ],
                )
            )

        edge_type = _infer_edge_type(ixn)
        edges.append(
            Edge(
                source=gene_id,
                target=chem_id,
                type=edge_type,
                label=ixn.interaction[:80] if ixn.interaction else None,
                carcinogen=chem_id,
                source_db="CTD",
                pmid=ixn.pubmed_ids[0] if ixn.pubmed_ids else None,
            )
        )

    raw_graph = KnowledgeGraph(
        nodes=[node.model_copy(update={"origin": RecordOrigin.SEEDED}) for node in nodes],
        edges=[edge.model_copy(update={"origin": RecordOrigin.SEEDED}) for edge in edges],
    )
    prepared_graph, _warnings = prepare_knowledge_graph(
        raw_graph,
        mode=mode,
        reference_graphs=[("ctd", raw_graph)],
    )
    return prepared_graph


def _infer_edge_type(ixn: ChemicalGeneInteraction) -> EdgeType:
    """Heuristically map a CTD interaction description to an EdgeType."""
    text = ixn.interaction.lower()
    if any(kw in text for kw in ("metabolis", "activat", "hydroxylat", "oxidat", "epoxid")):
        return EdgeType.ACTIVATES
    if any(kw in text for kw in ("conjugat", "detoxif", "glucuronid", "glutathione", "sulfat")):
        return EdgeType.DETOXIFIES
    if any(kw in text for kw in ("transport", "efflux", "export")):
        return EdgeType.TRANSPORTS
    if any(kw in text for kw in ("repair", "excision")):
        return EdgeType.REPAIRS
    if any(kw in text for kw in ("induc", "upregulat", "increas")):
        return EdgeType.INDUCES
    if any(kw in text for kw in ("inhibit", "downregulat", "decreas", "suppress")):
        return EdgeType.INHIBITS
    return EdgeType.ACTIVATES  # default


def seed_iarc_classification(
    chemical_name: str,
    *,
    classifier: Optional[IARCClassifier] = None,
) -> Optional[dict[str, str]]:
    """Look up IARC classification for a chemical.

    Returns a dict with ``group``, ``cas``, and ``category`` keys,
    or ``None`` if the chemical is not in the IARC dataset.

    This is a lightweight helper — it does not produce a full KnowledgeGraph
    but provides annotation data to enrich existing Carcinogen nodes.
    """
    if classifier is None:
        classifier = IARCClassifier()
    return classifier.get_entry(chemical_name)
