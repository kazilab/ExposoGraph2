"""Tests for ExposoGraph.reference_data."""

from ExposoGraph.engine import GraphEngine
from ExposoGraph.models import NodeType
from ExposoGraph.reference_data import (
    ACTIVITY_SCORE_METADATA,
    ACTIVITY_SCORES,
    CURATION_SOURCE_MANIFEST,
    REFERENCE_KEGG_PATHWAYS,
    TIER1_GENES,
    TIER2_GENES,
    build_full_panel,
    build_tier1_panel,
    build_tier2_panel,
    get_activity_score_metadata,
    get_activity_score_references,
    get_activity_scores,
)


class TestGenePanels:
    def test_curation_source_manifest_tracks_primary_and_supporting_sources(self):
        primary_sources = CURATION_SOURCE_MANIFEST["primary_sources"]
        supporting_sources = CURATION_SOURCE_MANIFEST["supporting_sources"]

        assert set(primary_sources) == {"IARC", "KEGG", "PharmVar", "CPIC", "CTD", "GTEx", "PubMed"}
        assert set(supporting_sources) == {"NCBI Gene", "ClinPGx"}
        assert primary_sources["KEGG"]["reference_pathways"] == [
            entry["pathway_id"] for entry in REFERENCE_KEGG_PATHWAYS
        ]

    def test_reference_kegg_pathways_cover_manuscript_catalog(self):
        pathway_ids = {entry["pathway_id"] for entry in REFERENCE_KEGG_PATHWAYS}
        assert pathway_ids == {"hsa00980", "hsa00140", "hsa05204", "hsa05208"}

    def test_tier1_has_13_genes(self):
        assert len(TIER1_GENES) == 13

    def test_tier2_has_23_genes(self):
        assert len(TIER2_GENES) == 23

    def test_build_tier1_panel(self):
        kg = build_tier1_panel()
        assert len(kg.nodes) == 13
        assert all(n.type == NodeType.ENZYME for n in kg.nodes)
        assert all(n.tier == 1 for n in kg.nodes)

    def test_build_tier2_panel(self):
        kg = build_tier2_panel()
        assert len(kg.nodes) == 23
        assert all(n.type == NodeType.ENZYME for n in kg.nodes)
        assert all(n.tier == 2 for n in kg.nodes)

    def test_build_full_panel(self):
        kg = build_full_panel()
        assert len(kg.nodes) == 36

    def test_no_duplicate_ids(self):
        kg = build_full_panel()
        ids = [n.id for n in kg.nodes]
        assert len(ids) == len(set(ids)), f"Duplicate IDs: {[x for x in ids if ids.count(x) > 1]}"

    def test_panel_loads_into_engine(self):
        engine = GraphEngine()
        kg = build_full_panel()
        warnings = engine.load(kg)
        assert warnings == []
        assert engine.node_count == 36

    def test_extended_panel_matches_manuscript_phase_distribution(self):
        panel = build_full_panel()
        phase_i = [node for node in panel.nodes if node.phase == "I"]
        phase_ii = [node for node in panel.nodes if node.phase == "II"]
        phase_iii = [node for node in panel.nodes if node.phase == "III"]
        repair = [node for node in panel.nodes if (node.group or "").startswith("DNA Repair")]

        assert len(phase_i) == 15
        assert len(phase_ii) == 13
        assert len(phase_iii) == 3
        assert len(repair) == 5

    def test_all_genes_have_required_fields(self):
        for gene_list in [TIER1_GENES, TIER2_GENES]:
            for g in gene_list:
                assert "id" in g
                assert "label" in g
                assert "role" in g
                assert "detail" in g
                assert "provenance" in g
                assert g["phase"] or g["group"]

    def test_all_genes_have_database_backed_provenance(self):
        for gene_list in [TIER1_GENES, TIER2_GENES]:
            for g in gene_list:
                provenance = g["provenance"]
                assert len(provenance) >= 2
                assert provenance[0]["source_db"] == "NCBI Gene"
                assert provenance[0]["record_id"]
                assert provenance[0]["url"].startswith("https://www.ncbi.nlm.nih.gov/gene/")

    def test_repair_genes_use_repair_class_not_phase(self):
        repair_ids = {"XRCC1", "XPC", "ERCC2", "OGG1", "MGMT"}
        repair_genes = [g for g in TIER2_GENES if g["id"] in repair_ids]
        assert len(repair_genes) == 5
        assert all(g["phase"] is None for g in repair_genes)
        assert all(g["group"].startswith("DNA Repair") for g in repair_genes)


class TestActivityScores:
    def test_known_gene(self):
        scores = get_activity_scores("CYP1A1")
        assert scores is not None
        assert len(scores) > 0

    def test_unknown_gene(self):
        assert get_activity_scores("FAKE_GENE") is None

    def test_all_entries_have_required_fields(self):
        for gene, entries in ACTIVITY_SCORES.items():
            for entry in entries:
                assert "allele" in entry, f"{gene}: missing allele"
                assert "value" in entry, f"{gene}: missing value"
                assert "phenotype" in entry, f"{gene}: missing phenotype"
                assert "confidence" in entry, f"{gene}: missing confidence"
                assert isinstance(entry["value"], (int, float)), f"{gene}: value not numeric"

    def test_activity_scores_cover_tier1_and_tier2(self):
        scored_genes = set(ACTIVITY_SCORES.keys())
        tier1_ids = {g["id"] for g in TIER1_GENES}
        tier2_ids = {g["id"] for g in TIER2_GENES}
        all_ids = tier1_ids | tier2_ids
        # At least the major genes should have scores
        assert len(scored_genes & all_ids) >= 15

    def test_activity_score_metadata_covers_all_scored_genes(self):
        assert set(ACTIVITY_SCORE_METADATA) == set(ACTIVITY_SCORES)

    def test_activity_score_references_have_source_and_url(self):
        for gene, metadata in ACTIVITY_SCORE_METADATA.items():
            assert metadata["evidence_basis"]
            assert metadata["note"]
            references = metadata["references"]
            assert len(references) >= 1
            for ref in references:
                assert ref["source_db"], f"{gene}: missing source_db"
                assert ref["citation"], f"{gene}: missing citation"
                assert ref["url"].startswith("https://"), f"{gene}: missing URL"

    def test_activity_score_metadata_lookup(self):
        metadata = get_activity_score_metadata("CYP2D6")
        assert metadata is not None
        assert metadata["evidence_basis"]
        references = get_activity_score_references("CYP2D6")
        assert references is not None
        assert len(references) >= 1

    def test_guideline_backed_activity_scores_have_explicit_cpic_references(self):
        for gene in ("CYP2D6", "CYP2C9", "CYP2C19", "UGT1A1"):
            references = get_activity_score_references(gene) or []
            assert "CPIC" in {ref["source_db"] for ref in references}

    def test_pharmvar_backed_pharmacogenes_keep_explicit_allele_definition_references(self):
        for gene in ("CYP2D6", "CYP2C9", "CYP2C19", "CYP2A6", "CYP3A4"):
            references = get_activity_score_references(gene) or []
            assert "PharmVar" in {ref["source_db"] for ref in references}


class TestNewNodeTypes:
    def test_gene_node_type(self):
        from ExposoGraph.models import Node
        node = Node(id="CYP1A1_gene", label="CYP1A1 gene", type=NodeType.GENE)
        assert node.type == NodeType.GENE

    def test_tissue_node_type(self):
        from ExposoGraph.models import Node
        node = Node(id="lung", label="Lung", type=NodeType.TISSUE)
        assert node.type == NodeType.TISSUE


class TestNewEdgeTypes:
    def test_expressed_in(self):
        from ExposoGraph.models import Edge, EdgeType
        edge = Edge(source="CYP1A1", target="lung", type=EdgeType.EXPRESSED_IN)
        assert edge.type == EdgeType.EXPRESSED_IN

    def test_induces(self):
        from ExposoGraph.models import Edge, EdgeType
        edge = Edge(source="BaP", target="CYP1A1", type=EdgeType.INDUCES)
        assert edge.type == EdgeType.INDUCES

    def test_inhibits(self):
        from ExposoGraph.models import Edge, EdgeType
        edge = Edge(source="drug", target="CYP1A1", type=EdgeType.INHIBITS)
        assert edge.type == EdgeType.INHIBITS

    def test_encodes(self):
        from ExposoGraph.models import Edge, EdgeType
        edge = Edge(source="CYP1A1_gene", target="CYP1A1", type=EdgeType.ENCODES)
        assert edge.type == EdgeType.ENCODES
