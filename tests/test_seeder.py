"""Tests for ExposoGraph.seeder."""

from unittest.mock import MagicMock

import pytest

from ExposoGraph.config import GraphMode
from ExposoGraph.db_clients.ctd import ChemicalGeneInteraction, CTDClient
from ExposoGraph.db_clients.iarc import IARCClassifier
from ExposoGraph.db_clients.kegg import KEGGClient, KEGGPathway
from ExposoGraph.models import EdgeType, KnowledgeGraph, MatchStatus, NodeType, RecordOrigin
from ExposoGraph.seeder import (
    _infer_edge_type,
    seed_from_ctd,
    seed_from_kegg_pathway,
    seed_iarc_classification,
)

# ── KEGG seeder ───────────────────────────────────────────────────────────


class TestSeedFromKEGG:
    def test_creates_pathway_and_genes(self):
        mock_client = MagicMock(spec=KEGGClient)
        mock_client.get_pathway.return_value = KEGGPathway(
            pathway_id="hsa05204",
            name="Chemical Carcinogenesis",
            genes=["CYP1A1", "CYP1B1", "GSTM1"],
        )

        kg = seed_from_kegg_pathway("hsa05204", client=mock_client)

        assert isinstance(kg, KnowledgeGraph)
        node_ids = {n.id for n in kg.nodes}
        assert "hsa05204" in node_ids
        assert "CYP1A1" in node_ids
        assert "CYP1B1" in node_ids
        assert "GSTM1" in node_ids
        assert len(kg.nodes) == 4  # 1 pathway + 3 genes

    def test_creates_pathway_edges(self):
        mock_client = MagicMock(spec=KEGGClient)
        mock_client.get_pathway.return_value = KEGGPathway(
            pathway_id="hsa05204",
            name="Test",
            genes=["CYP1A1"],
        )

        kg = seed_from_kegg_pathway("hsa05204", client=mock_client)

        assert len(kg.edges) == 1
        assert kg.edges[0].type == EdgeType.PATHWAY
        assert kg.edges[0].source == "CYP1A1"
        assert kg.edges[0].target == "hsa05204"

    def test_provenance_set(self):
        mock_client = MagicMock(spec=KEGGClient)
        mock_client.get_pathway.return_value = KEGGPathway(
            pathway_id="hsa05204",
            name="Test",
            genes=[],
        )

        kg = seed_from_kegg_pathway("hsa05204", client=mock_client)

        assert kg.nodes[0].provenance[0].source_db == "KEGG"

    def test_empty_pathway(self):
        mock_client = MagicMock(spec=KEGGClient)
        mock_client.get_pathway.return_value = KEGGPathway(
            pathway_id="hsa00000", name="Empty", genes=[],
        )

        kg = seed_from_kegg_pathway("hsa00000", client=mock_client)
        assert len(kg.nodes) == 1  # just the pathway node
        assert len(kg.edges) == 0

    def test_strict_mode_keeps_seeded_pathway_graph(self):
        mock_client = MagicMock(spec=KEGGClient)
        mock_client.get_pathway.return_value = KEGGPathway(
            pathway_id="hsa05204",
            name="Chemical Carcinogenesis",
            genes=["CYP1A1"],
        )

        kg = seed_from_kegg_pathway("hsa05204", client=mock_client, mode=GraphMode.STRICT)

        assert len(kg.nodes) == 2
        assert len(kg.edges) == 1
        assert all(node.origin == RecordOrigin.SEEDED for node in kg.nodes)
        assert all(
            node.match_status in (MatchStatus.CANONICAL, MatchStatus.ALIAS) for node in kg.nodes
        )
        assert kg.edges[0].match_status == MatchStatus.CANONICAL


# ── CTD seeder ────────────────────────────────────────────────────────────


class TestSeedFromCTD:
    def test_creates_carcinogen_and_genes(self):
        mock_client = MagicMock(spec=CTDClient)
        mock_client.get_chemical_gene_interactions.return_value = [
            ChemicalGeneInteraction(
                chemical_name="BaP", chemical_id="D001564",
                gene_symbol="CYP1A1", gene_id="1543",
                interaction="metabolism", pubmed_ids=["12345678"],
            ),
            ChemicalGeneInteraction(
                chemical_name="BaP", chemical_id="D001564",
                gene_symbol="GSTM1", gene_id="2944",
                interaction="glutathione conjugation", pubmed_ids=["23456789"],
            ),
        ]

        kg = seed_from_ctd("BaP", client=mock_client)

        assert isinstance(kg, KnowledgeGraph)
        node_ids = {n.id for n in kg.nodes}
        assert "BaP" in node_ids
        assert "CYP1A1" in node_ids
        assert "GSTM1" in node_ids

    def test_edge_types_inferred(self):
        mock_client = MagicMock(spec=CTDClient)
        mock_client.get_chemical_gene_interactions.return_value = [
            ChemicalGeneInteraction(
                chemical_name="BaP", chemical_id="D001564",
                gene_symbol="CYP1A1", gene_id="1543",
                interaction="metabolic activation",
            ),
            ChemicalGeneInteraction(
                chemical_name="BaP", chemical_id="D001564",
                gene_symbol="GSTM1", gene_id="2944",
                interaction="glutathione conjugation detoxification",
            ),
        ]

        kg = seed_from_ctd("BaP", client=mock_client)

        edge_types = {e.source: e.type for e in kg.edges}
        assert edge_types["CYP1A1"] == EdgeType.ACTIVATES
        assert edge_types["GSTM1"] == EdgeType.DETOXIFIES

    def test_deduplicates_genes(self):
        mock_client = MagicMock(spec=CTDClient)
        mock_client.get_chemical_gene_interactions.return_value = [
            ChemicalGeneInteraction(
                chemical_name="BaP", chemical_id="D001564",
                gene_symbol="CYP1A1", gene_id="1543",
                interaction="activation",
            ),
            ChemicalGeneInteraction(
                chemical_name="BaP", chemical_id="D001564",
                gene_symbol="CYP1A1", gene_id="1543",
                interaction="hydroxylation",
            ),
        ]

        kg = seed_from_ctd("BaP", client=mock_client)

        gene_nodes = [n for n in kg.nodes if n.type == NodeType.GENE]
        assert len(gene_nodes) == 1  # CYP1A1 appears only once
        assert len(kg.edges) == 2  # but both interactions become edges

    def test_empty_interactions(self):
        mock_client = MagicMock(spec=CTDClient)
        mock_client.get_chemical_gene_interactions.return_value = []

        kg = seed_from_ctd("NonexistentChem", client=mock_client)

        assert len(kg.nodes) == 1  # just the carcinogen
        assert len(kg.edges) == 0

    def test_strict_mode_keeps_seeded_ctd_graph(self):
        mock_client = MagicMock(spec=CTDClient)
        mock_client.get_chemical_gene_interactions.return_value = [
            ChemicalGeneInteraction(
                chemical_name="Novel Chemical",
                chemical_id="D999999",
                gene_symbol="NOVEL1",
                gene_id="9999",
                interaction="metabolism",
            ),
        ]

        kg = seed_from_ctd("Novel Chemical", client=mock_client, mode=GraphMode.STRICT)

        assert len(kg.nodes) == 2
        assert len(kg.edges) == 1
        assert all(node.origin == RecordOrigin.SEEDED for node in kg.nodes)
        assert all(node.match_status == MatchStatus.CANONICAL for node in kg.nodes)
        assert kg.edges[0].match_status == MatchStatus.CANONICAL


# ── Edge type inference ───────────────────────────────────────────────────


class TestInferEdgeType:
    @pytest.mark.parametrize("text,expected", [
        ("metabolic activation", EdgeType.ACTIVATES),
        ("CYP1A1 hydroxylation of BaP", EdgeType.ACTIVATES),
        ("epoxidation reaction", EdgeType.ACTIVATES),
        ("glutathione conjugation", EdgeType.DETOXIFIES),
        ("detoxification via GSTM1", EdgeType.DETOXIFIES),
        ("efflux transport", EdgeType.TRANSPORTS),
        ("nucleotide excision repair", EdgeType.REPAIRS),
        ("induces expression", EdgeType.INDUCES),
        ("upregulates CYP1A1", EdgeType.INDUCES),
        ("inhibits enzyme activity", EdgeType.INHIBITS),
        ("downregulates expression", EdgeType.INHIBITS),
        ("unknown interaction", EdgeType.ACTIVATES),  # default
    ])
    def test_inference(self, text, expected):
        ixn = ChemicalGeneInteraction(
            chemical_name="X", chemical_id="Y",
            gene_symbol="Z", gene_id="0",
            interaction=text,
        )
        assert _infer_edge_type(ixn) == expected


# ── IARC seeder ───────────────────────────────────────────────────────────


class TestSeedIARC:
    def test_returns_entry(self):
        result = seed_iarc_classification("Benzo[a]pyrene")
        assert result is not None
        assert result["group"] == "Group 1"

    def test_returns_none_for_unknown(self):
        result = seed_iarc_classification("NotAChemical")
        assert result is None

    def test_with_custom_classifier(self):
        clf = IARCClassifier(
            extra={"TestChem": {"group": "Group 2A", "cas": "", "category": "Test"}}
        )
        result = seed_iarc_classification("TestChem", classifier=clf)
        assert result is not None
        assert result["group"] == "Group 2A"
