"""Tests for canonical reference and seeded showcase graphs."""

from collections import Counter
from inspect import signature

from ExposoGraph import (
    build_androgen_module_engine,
    build_androgen_module_graph,
    build_full_legends_architecture_summary,
    build_full_legends_engine,
    build_full_legends_graph,
    build_reference_architecture_summary,
    build_reference_engine,
    build_reference_graph,
    parse_graph_artifact,
    pathway_subgraph,
    write_reference_exports,
)


def test_reference_graph_matches_bundled_counts():
    kg = build_reference_graph()

    assert len(kg.nodes) == 214
    assert len(kg.edges) == 321


def test_reference_engine_loads_and_validates():
    engine = build_reference_engine()

    assert engine.node_count == 214
    assert engine.edge_count == 321
    assert engine.validate() == []


def test_reference_architecture_summary_matches_bundled_graph():
    summary = build_reference_architecture_summary()

    assert summary.node_count == 214
    assert summary.edge_count == 321
    assert summary.node_type_count == len(summary.node_type_counts)
    assert summary.edge_type_count == len(summary.edge_type_counts)
    assert sum(summary.node_type_counts.values()) == summary.node_count
    assert sum(summary.edge_type_counts.values()) == summary.edge_count


def test_reference_exports_write_consistent_artifacts(tmp_path):
    kwargs = {}
    if "bundle_dir" in signature(write_reference_exports).parameters:
        kwargs["bundle_dir"] = tmp_path / "reference_graph.js"

    artifacts = write_reference_exports(tmp_path, **kwargs)

    assert {"html", "json", "plotly_html", "graph_data_js"} <= set(artifacts)
    for path in artifacts.values():
        assert path.exists()

    for key in ("html", "json", "graph_data_js"):
        graph = parse_graph_artifact(artifacts[key])
        assert len(graph.nodes) == 214
        assert len(graph.edges) == 321


def test_full_legends_graph_matches_phase2_target_counts():
    kg = build_full_legends_graph()
    node_types = Counter(node.type.value for node in kg.nodes)
    edge_types = Counter(edge.type.value for edge in kg.edges)

    assert len(kg.nodes) == 214
    assert len(kg.edges) == 321
    assert node_types == {
        "Carcinogen": 56,
        "Enzyme": 60,
        "Metabolite": 59,
        "DNA_Adduct": 27,
        "Pathway": 12,
    }
    assert edge_types == {
        "ACTIVATES": 104,
        "DETOXIFIES": 47,
        "FORMS_ADDUCT": 37,
        "PATHWAY": 88,
        "REPAIRS": 33,
        "TRANSPORTS": 9,
        "INDUCES": 2,
        "INHIBITS": 1,
    }


def test_full_legends_engine_loads_and_validates():
    engine = build_full_legends_engine()

    assert engine.node_count == 214
    assert engine.edge_count == 321
    assert engine.validate() == []


def test_full_legends_graph_keeps_key_entities():
    engine = build_full_legends_engine()

    for node_id in (
        "DMBA",
        "MeIQx",
        "Benzidine",
        "NDMA",
        "E2",
        "DHT",
        "Benzene",
        "VinylChloride",
        "EthyleneOxide",
        "CYP17A1",
        "SRD5A1",
        "SRD5A2",
        "CYP19A1",
        "AKR1C3",
        "UGT2B17",
        "UGT2B15",
        "AKR1C2",
        "HSD3B2",
        "CYP3A5",
        "COMT",
        "CYP2B6",
        "CYP2C9",
        "CYP2C19",
        "CYP2D6",
        "CYP2F1",
    ):
        assert engine.get_node(node_id) is not None


def test_full_legends_graph_exposes_curated_kegg_pathways():
    engine = build_full_legends_engine()

    members = pathway_subgraph(engine, "hsa00140")
    assert "CYP17A1" in members
    assert "Testosterone" in members

    members = pathway_subgraph(engine, "hsa05204")
    assert "AFB1" in members
    assert "NNK" in members

    members = pathway_subgraph(engine, "hsa00980")
    assert "CYP2F1" in members

    members = pathway_subgraph(engine, "hsa00982")
    for gene in ("CYP2B6", "CYP2C9", "CYP2C19", "CYP2D6"):
        assert gene in members


def test_full_legends_architecture_summary_matches_seeded_graph():
    summary = build_full_legends_architecture_summary()

    assert summary.node_count == 214
    assert summary.edge_count == 321
    assert summary.node_type_count == 5
    assert summary.edge_type_count == 8
    assert summary.node_type_counts == {
        "Carcinogen": 56,
        "Enzyme": 60,
        "Metabolite": 59,
        "DNA_Adduct": 27,
        "Pathway": 12,
    }
    assert summary.edge_type_counts == {
        "ACTIVATES": 104,
        "DETOXIFIES": 47,
        "TRANSPORTS": 9,
        "FORMS_ADDUCT": 37,
        "REPAIRS": 33,
        "PATHWAY": 88,
        "INDUCES": 2,
        "INHIBITS": 1,
    }


def test_full_legends_architecture_summary_keeps_inventories():
    summary = build_full_legends_architecture_summary()
    carcinogen_classes = {group.name: group.count for group in summary.carcinogen_classes}
    enzyme_categories = {group.name: group.count for group in summary.enzyme_categories}

    assert carcinogen_classes == {
        "PAH": 2,
        "HCA": 2,
        "Aromatic Amines": 2,
        "Nitrosamines": 3,
        "Mycotoxins": 1,
        "Estrogens": 1,
        "Androgens": 2,
        "Solvents": 2,
        "Alkylating Agents": 9,
    }
    assert enzyme_categories == {
        "Phase I": 22,
        "Phase II": 21,
        "Phase III": 4,
        "DNA Repair": 8,
    }
    assert "DMBA" in summary.carcinogens
    assert "5a-DHT" in summary.carcinogens
    assert "CYP19A1" in summary.enzymes
    assert "UGT2B15" in summary.enzymes
    assert "Chemical carcinogenesis - DNA adducts" in summary.pathway_labels


def test_androgen_module_graph_has_receptor_variant_and_tissue_context():
    kg = build_androgen_module_graph()
    node_types = Counter(node.type.value for node in kg.nodes)
    edge_types = Counter(edge.type.value for edge in kg.edges)
    edge_index = {
        (edge.source, edge.target, edge.type.value, edge.custom_predicate)
        for edge in kg.edges
    }

    assert len(kg.nodes) == 31
    assert len(kg.edges) == 41
    assert node_types == {
        "Carcinogen": 3,
        "Enzyme": 9,
        "Gene": 5,
        "Metabolite": 6,
        "DNA_Adduct": 2,
        "Pathway": 3,
        "Tissue": 3,
    }
    assert edge_types == {
        "ACTIVATES": 3,
        "CUSTOM": 7,
        "DETOXIFIES": 3,
        "FORMS_ADDUCT": 2,
        "PATHWAY": 15,
        "EXPRESSED_IN": 7,
        "ENCODES": 4,
    }
    assert ("SRD5A2", "DHT", "CUSTOM", "CONVERTS_TO_DHT") in edge_index
    assert ("CYP19A1", "E2", "CUSTOM", "AROMATIZES_TO_ESTRADIOL") in edge_index
    assert ("DHT", "AR", "CUSTOM", "BINDS_RECEPTOR") in edge_index
    assert ("AR", "AR_signal_program", "CUSTOM", "ACTIVATES_TRANSCRIPTION") in edge_index
    assert ("AR", "Prostate", "EXPRESSED_IN", None) in edge_index
    assert ("SRD5A2_V89L", "SRD5A2", "ENCODES", None) in edge_index


def test_androgen_module_engine_loads_and_exposes_variant_annotations():
    engine = build_androgen_module_engine()
    ar_node = engine.get_node("AR")
    srd5a2_v89l = engine.get_node("SRD5A2_V89L")
    ugt2b17_deletion = engine.get_node("UGT2B17_copy_number_deletion")

    assert engine.node_count == 31
    assert engine.edge_count == 41
    assert engine.validate() == []
    assert engine.get_node("CYP3A5") is not None
    assert ar_node["type"] == "Gene"
    assert srd5a2_v89l["variant"] == "V89L"
    assert ugt2b17_deletion["phenotype"].startswith("Absent")


def test_full_legends_graph_can_merge_optional_androgen_module():
    kg = build_full_legends_graph(include_androgen_module=True)
    node_types = Counter(node.type.value for node in kg.nodes)
    edge_types = Counter(edge.type.value for edge in kg.edges)

    assert len(kg.nodes) == 225
    assert len(kg.edges) == 353
    assert node_types == {
        "Carcinogen": 56,
        "Enzyme": 60,
        "Gene": 5,
        "Metabolite": 59,
        "DNA_Adduct": 29,
        "Pathway": 13,
        "Tissue": 3,
    }
    assert edge_types == {
        "ACTIVATES": 106,
        "CUSTOM": 7,
        "DETOXIFIES": 48,
        "TRANSPORTS": 9,
        "FORMS_ADDUCT": 39,
        "REPAIRS": 33,
        "PATHWAY": 97,
        "EXPRESSED_IN": 7,
        "ENCODES": 4,
        "INDUCES": 2,
        "INHIBITS": 1,
    }


def test_full_legends_architecture_summary_can_include_androgen_module():
    summary = build_full_legends_architecture_summary(include_androgen_module=True)
    enzyme_categories = {group.name: group.count for group in summary.enzyme_categories}

    assert summary.node_count == 225
    assert summary.edge_count == 353
    assert summary.node_type_count == 7
    assert summary.edge_type_count == 11
    assert summary.node_type_counts == {
        "Carcinogen": 56,
        "Enzyme": 60,
        "Gene": 5,
        "Metabolite": 59,
        "DNA_Adduct": 29,
        "Pathway": 13,
        "Tissue": 3,
    }
    assert summary.edge_type_counts == {
        "ACTIVATES": 106,
        "DETOXIFIES": 48,
        "TRANSPORTS": 9,
        "FORMS_ADDUCT": 39,
        "REPAIRS": 33,
        "PATHWAY": 97,
        "EXPRESSED_IN": 7,
        "ENCODES": 4,
        "CUSTOM": 7,
        "INDUCES": 2,
        "INHIBITS": 1,
    }
    assert enzyme_categories == {
        "Phase I": 22,
        "Phase II": 21,
        "Phase III": 4,
        "DNA Repair": 8,
    }
    assert "AR proliferative transcriptional program" in summary.pathway_labels
    assert "CYP3A5" in summary.enzymes
