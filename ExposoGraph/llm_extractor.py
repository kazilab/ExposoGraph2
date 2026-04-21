"""LLM-powered entity / relation extraction for the knowledge graph.

Supports multiple LLM backends (OpenAI, Ollama) via a pluggable protocol.
Falls back to JSON-mode parsing when structured output is unavailable.
"""

from __future__ import annotations

import logging
from typing import Optional

from .config import GraphMode
from .grounding import prepare_knowledge_graph
from .llm_backend import LLMBackend, OpenAIBackend, UsageRecord
from .models import KnowledgeGraph, RecordOrigin

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """\
You are an expert biochemist and toxicologist. Your task is to extract a
structured knowledge graph from the user's natural-language description of
carcinogen metabolism, gene interactions, and DNA damage pathways.

Return **only** valid JSON matching the schema below. Do not include any
text outside the JSON block.

### Node types (use exactly these strings for "type"):
- Carcinogen  — chemical agents; include "group" (e.g. PAH, HCA, Aromatic_Amine,
  Nitrosamine, Mycotoxin, Estrogen, Androgen, Solvent, Alkylating, Heavy_Metal)
  and "iarc" classification (Group 1 / 2A / 2B / 3). For heavy metal carcinogens
  include "exposure" (primary human exposure routes).
- Enzyme      — include "phase" (I, II, III, when applicable) and "role"
  (Activation, Detoxification, Mixed, Transport, Repair). For DNA repair
  proteins, use "role": "Repair" and store the repair class in "group"
  (for example "DNA Repair (BER)" or "DNA Repair (NER)") instead of using
  "phase": "Repair".
- Gene        — a gene locus (use when distinguishing the gene from its encoded enzyme,
  e.g. for pharmacogenomic variants or tissue expression context).
- Metabolite  — include "reactivity" (High, Intermediate, Low).
- DNA_Adduct  — DNA lesion types.
- Pathway     — biological pathways; use KEGG IDs when possible.
- Tissue      — anatomical tissue or organ where expression/metabolism occurs.

### Edge types (use exactly these strings for "type"):
- ACTIVATES     — enzyme activates a procarcinogen → reactive metabolite
- DETOXIFIES    — enzyme conjugates / inactivates a metabolite
- TRANSPORTS    — efflux transporter moves a conjugate out of the cell
- FORMS_ADDUCT  — reactive metabolite covalently modifies DNA
- REPAIRS       — DNA repair enzyme removes a lesion
- PATHWAY       — node belongs to a biological pathway
- EXPRESSED_IN  — gene or enzyme is expressed in a tissue
- INDUCES       — substance or exposure induces enzyme expression/activity
- INHIBITS      — substance or exposure inhibits enzyme expression/activity
- ENCODES       — gene encodes an enzyme

### JSON Schema:
{
  "nodes": [
    {
      "id": "<short_unique_id>",
      "label": "<display name>",
      "type": "<NodeType>",
      "detail": "<one-line description>",
      "group": "<carcinogen class or repair class, or null>",
      "iarc": "<IARC group or null>",
      "phase": "<enzyme phase or null>",
      "role": "<enzyme role or null>",
      "reactivity": "<metabolite reactivity or null>",
      "source_db": "<supporting database(s) such as NCBI Gene, GTEx, ClinPGx, CTD, IARC, or KEGG, or null>",
      "evidence": "<brief evidence note or null>",
      "pmid": "<PubMed ID or null>",
      "tissue": "<relevant tissue context or null>",
      "variant": "<star allele or variant name or null>",
      "phenotype": "<functional phenotype such as poor metabolizer or null>",
      "activity_score": "<numeric activity score or null>",
      "tier": "<gene panel tier: 1, 2, or null>",
      "exposure": "<primary exposure routes for carcinogens, or null>"
    }
  ],
  "edges": [
    {
      "source": "<source node id>",
      "target": "<target node id>",
      "type": "<EdgeType>",
      "label": "<short description of the reaction>",
      "carcinogen": "<id of the parent carcinogen, if applicable, or null>",
      "source_db": "<supporting database(s) such as NCBI Gene, CTD, IARC, or KEGG, or null>",
      "evidence": "<brief evidence note or null>",
      "pmid": "<PubMed ID or null>",
      "tissue": "<relevant tissue context or null>"
    }
  ]
}

Guidelines:
- Generate concise, uppercase-safe IDs (e.g. "BaP", "CYP1A1", "BPDE_dG").
- Every edge's source and target MUST reference an id that exists in the nodes list.
- Include the full metabolic chain: activation → metabolite → adduct → repair.
- Also include detoxification / conjugation branches when mentioned.
- If the user mentions KEGG pathway IDs, include Pathway nodes.
- Add annotation fields only when supported by the text; otherwise return null.
- Use `source_db` to reflect database-style provenance such as NCBI Gene, GTEx, ClinPGx, CTD, IARC, and KEGG.
- Capture tissue specificity, pharmacogenomic variants, and metabolizer phenotype when the text provides them.
"""


def extract_graph(
    text: str,
    *,
    model: str = "gpt-4o",
    api_key: Optional[str] = None,
    base_url: Optional[str] = None,
    backend: Optional[LLMBackend] = None,
    mode: GraphMode | str = GraphMode.EXPLORATORY,
) -> KnowledgeGraph:
    """Send *text* to the LLM and return a validated KnowledgeGraph.

    If *backend* is provided it is used directly; otherwise an
    :class:`OpenAIBackend` is created from the given credentials.
    """
    result, _usage = extract_graph_with_usage(
        text, model=model, api_key=api_key, base_url=base_url, backend=backend, mode=mode,
    )
    return result


def extract_graph_with_usage(
    text: str,
    *,
    model: str = "gpt-4o",
    api_key: Optional[str] = None,
    base_url: Optional[str] = None,
    backend: Optional[LLMBackend] = None,
    mode: GraphMode | str = GraphMode.EXPLORATORY,
) -> tuple[KnowledgeGraph, UsageRecord]:
    """Like :func:`extract_graph` but also returns token usage metadata."""
    if backend is None:
        backend = OpenAIBackend(api_key=api_key, base_url=base_url)

    raw, usage = backend.extract_json(text, SYSTEM_PROMPT, model)
    kg = KnowledgeGraph(**raw)
    kg = KnowledgeGraph(
        nodes=[
            node.model_copy(update={"origin": RecordOrigin.LLM})
            for node in kg.nodes
        ],
        edges=[
            edge.model_copy(update={"origin": RecordOrigin.LLM})
            for edge in kg.edges
        ],
    )
    prepared_graph, warnings = prepare_knowledge_graph(kg, mode=mode)
    for warning in warnings:
        logger.warning(warning)
    return prepared_graph, usage


EXAMPLE_INPUT = """\
Inorganic arsenic (iAs) is an IARC Group 1 heavy metal carcinogen with
primary exposure from contaminated groundwater and rice. Arsenic
methyltransferase AS3MT methylates iAs(III) to MMA(V), then to DMA(V)
(the less toxic urinary excretion form). GSTO1 catalyzes the rate-limiting
reduction of MMA(V) to MMA(III), the most genotoxic arsenic metabolite.
MMA(III) generates reactive oxygen species (ROS) that cause 8-OHdG
oxidative DNA lesions, repaired by OGG1 via base excision repair. PARP1,
a zinc-finger BER enzyme, is inhibited by arsenic at low concentrations.
GSTT1 and GSTM1 detoxify arsenic via glutathione conjugation (GSTT1-null
genotype with high arsenic exposure gives OR 4.08 for urothelial carcinoma).
ABCC2 (MRP2) effluxes arsenic-glutathione conjugates. Key pharmacogenomic
variants include AS3MT haplotypes (rs3740393) that alter MMA/DMA ratios
and cancer risk, and GSTO1 Ala140Asp (rs4925).
"""
