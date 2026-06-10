"""Integration test: build BaP graph → analyze → export → reimport → verify.

This test exercises the full pipeline end-to-end without any mocks,
using the pre-built example graph JSON.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ExposoGraph import (
    GraphEngine,
    KnowledgeGraph,
    all_shortest_paths,
    centrality,
    metabolism_chain,
    pathway_subgraph,
    shortest_path,
    variant_impact_score,
)
from ExposoGraph.db_clients.iarc import IARCClassifier, IARCGroup
from ExposoGraph.exporter import (
    parse_graph_data_text,
    to_gexf,
    to_graph_data_js,
    to_interactive_html_string,
    to_json,
)
from ExposoGraph.models import EdgeType, NodeType

EXAMPLE_JSON = Path(__file__).parent.parent / "examples" / "bap_graph.json"


@pytest.fixture(scope="module")
def bap_kg() -> KnowledgeGraph:
    """Load the pre-built BaP example graph."""
    with open(EXAMPLE_JSON, encoding="utf-8") as f:
        data = json.load(f)
    return KnowledgeGraph(**data)


@pytest.fixture(scope="module")
def engine(bap_kg: KnowledgeGraph) -> GraphEngine:
    eng = GraphEngine()
    eng.load(bap_kg)
    return eng


# ── Load & validate ──────────────────────────────────────────────────────


class TestLoadAndValidate:
    def test_node_count(self, engine: GraphEngine):
        assert engine.node_count == 20

    def test_edge_count(self, engine: GraphEngine):
        assert engine.edge_count == 20

    def test_all_node_types_present(self, engine: GraphEngine):
        kg = engine.to_knowledge_graph()
        node_types = {n.type for n in kg.nodes}
        assert NodeType.CARCINOGEN in node_types
        assert NodeType.ENZYME in node_types
        assert NodeType.METABOLITE in node_types
        assert NodeType.DNA_ADDUCT in node_types
        assert NodeType.PATHWAY in node_types
        assert NodeType.TISSUE in node_types

    def test_all_edge_types_present(self, engine: GraphEngine):
        kg = engine.to_knowledge_graph()
        edge_types = {e.type for e in kg.edges}
        assert EdgeType.ACTIVATES in edge_types
        assert EdgeType.DETOXIFIES in edge_types
        assert EdgeType.TRANSPORTS in edge_types
        assert EdgeType.FORMS_ADDUCT in edge_types
        assert EdgeType.REPAIRS in edge_types
        assert EdgeType.PATHWAY in edge_types
        assert EdgeType.EXPRESSED_IN in edge_types

    def test_validation_passes(self, engine: GraphEngine):
        errors = engine.validate()
        assert errors == []

    def test_key_nodes_present(self, engine: GraphEngine):
        for node_id in ["BaP", "CYP1A1", "CYP1B1", "GSTM1", "BPDE", "BPDE_dG", "XPC"]:
            assert engine.get_node(node_id) is not None, f"Missing node: {node_id}"


# ── Analysis ─────────────────────────────────────────────────────────────


class TestAnalysis:
    def test_shortest_path_activation_to_adduct(self, engine: GraphEngine):
        path = shortest_path(engine, "CYP1A1", "BPDE_dG")
        assert path is not None
        assert path[0] == "CYP1A1"
        assert path[-1] == "BPDE_dG"

    def test_all_shortest_paths(self, engine: GraphEngine):
        paths = all_shortest_paths(engine, "CYP1A1", "BPDE_dG")
        assert len(paths) >= 1

    def test_centrality_cyp1a1_high(self, engine: GraphEngine):
        scores = centrality(engine, method="degree")
        assert scores["CYP1A1"] > scores["Lung"]

    def test_betweenness_centrality(self, engine: GraphEngine):
        scores = centrality(engine, method="betweenness")
        assert len(scores) == engine.node_count

    def test_metabolism_chain(self, engine: GraphEngine):
        chain = metabolism_chain(engine, "BaP")
        assert len(chain.node_ids) >= 5
        assert len(chain.activation_edges) >= 2
        assert len(chain.detox_edges) >= 1
        assert len(chain.adduct_edges) >= 1
        assert len(chain.repair_edges) >= 1

    def test_pathway_subgraph(self, engine: GraphEngine):
        members = pathway_subgraph(engine, "hsa05204")
        assert "BaP" in members
        assert "CYP1A1" in members

    def test_variant_impact_cyp1a1(self, engine: GraphEngine):
        impact = variant_impact_score(engine, "CYP1A1")
        assert impact is not None
        assert impact.activity_score == 1.0
        assert impact.downstream_adduct_count >= 1
        assert impact.score >= 0

    def test_variant_impact_gstm1_null(self, engine: GraphEngine):
        impact = variant_impact_score(engine, "GSTM1")
        assert impact is not None
        assert impact.activity_score == 0.0


# ── Export & reimport ────────────────────────────────────────────────────


class TestExportReimport:
    def test_json_roundtrip(self, engine: GraphEngine, tmp_path):
        out = to_json(engine, tmp_path / "graph.json")
        data = json.loads(out.read_text(encoding="utf-8"))
        restored = KnowledgeGraph(**data)
        assert len(restored.nodes) == engine.node_count

    def test_json_roundtrip_via_engine(self, engine: GraphEngine):
        json_str = engine.to_json()
        data = json.loads(json_str)
        new_engine = GraphEngine()
        new_engine.load(KnowledgeGraph(**data))
        assert new_engine.node_count == engine.node_count
        assert new_engine.edge_count == engine.edge_count

    def test_graph_data_js_roundtrip(self, engine: GraphEngine, tmp_path):
        out = to_graph_data_js(engine, tmp_path / "graph-data.js")
        js_text = out.read_text(encoding="utf-8")
        assert "GRAPH_DATA" in js_text
        restored = parse_graph_data_text(js_text)
        assert len(restored.nodes) == engine.node_count

    def test_gexf_export(self, engine: GraphEngine, tmp_path):
        out = to_gexf(engine, tmp_path / "graph.gexf")
        gexf_str = out.read_text(encoding="utf-8")
        assert "<?xml" in gexf_str
        assert "BaP" in gexf_str
        assert "CYP1A1" in gexf_str

    def test_interactive_html_export(self, engine: GraphEngine):
        try:
            html = to_interactive_html_string(engine)
        except FileNotFoundError:
            pytest.skip("Viewer HTML template not found in this environment")
        assert "<html" in html.lower()
        assert "GRAPH_DATA" in html

    def test_reimported_graph_validates(self, engine: GraphEngine, tmp_path):
        out = to_json(engine, tmp_path / "graph.json")
        data = json.loads(out.read_text(encoding="utf-8"))
        restored = KnowledgeGraph(**data)
        new_engine = GraphEngine()
        new_engine.load(restored)
        assert new_engine.validate() == []


# ── IARC enrichment ──────────────────────────────────────────────────────


class TestIARCEnrichment:
    def test_bap_classified_group_1(self):
        clf = IARCClassifier()
        assert clf.classify("Benzo[a]pyrene") == IARCGroup.GROUP_1

    def test_enrichment_matches_graph(self, bap_kg: KnowledgeGraph):
        clf = IARCClassifier()
        bap_node = next(n for n in bap_kg.nodes if n.id == "BaP")
        iarc_entry = clf.get_entry("Benzo[a]pyrene")
        assert iarc_entry is not None
        assert bap_node.iarc == iarc_entry["group"]


class TestLocalBiomarkerScaffoldUtilities:
    def test_local_connector_scaffolds_use_curated_files_and_cache(self, tmp_path):
        from ExposoGraph._biomarker_scaffold.scripts.integrations.atsdR_connector import (
            CuratedToxicologyConnector as AtsdrConnector,
        )
        from ExposoGraph._biomarker_scaffold.scripts.integrations.brenda_connector import (
            BrendaConnector,
        )
        from ExposoGraph._biomarker_scaffold.scripts.integrations.comptox_connector import (
            CompToxConnector,
        )
        from ExposoGraph._biomarker_scaffold.scripts.integrations.iarc_connector import (
            CuratedToxicologyConnector as IarcConnector,
        )
        from ExposoGraph._biomarker_scaffold.scripts.integrations.iris_connector import (
            CuratedToxicologyConnector as IrisConnector,
        )
        from ExposoGraph._biomarker_scaffold.scripts.integrations.pubchem_connector import (
            PubChemConnector,
        )
        from ExposoGraph._biomarker_scaffold.scripts.integrations.source_base import (
            EvidenceSource,
        )

        source = EvidenceSource("fixture")
        with pytest.raises(NotImplementedError):
            source.search("benzene")
        with pytest.raises(NotImplementedError):
            source.fetch("benzene")
        assert source.normalize({"a": 1}) == {"a": 1}

        brenda_csv = tmp_path / "brenda.csv"
        brenda_csv.write_text(
            "enzyme,EC,substrate,Km_uM,confidence,organism,tissue,reference\n"
            "CYP2E1,1.14.14.1,benzene,42,0.8,human,liver,PMID:1\n",
            encoding="utf-8",
        )
        brenda = BrendaConnector(brenda_csv)
        assert brenda.search("benzene")[0]["Km_uM"] == 42.0
        assert brenda.fetch("CYP2E1")["confidence"] == 0.8
        assert brenda.fetch("missing")["source_status"] == "missing_curated_record"
        assert brenda.normalize({"enzyme": "CYP", "substrate": "x", "Km": "", "confidence": ""})[
            "confidence"
        ] == 0.5

        toxicology_yaml = tmp_path / "toxicology.yaml"
        toxicology_yaml.write_text(
            "records:\n"
            "  - chemical: benzene\n"
            "    identifier: benzene-id\n"
            "    cas: 71-43-2\n"
            "    carcinogen_class: known\n"
            "    reference: PMID:2\n",
            encoding="utf-8",
        )
        for connector_cls in (AtsdrConnector, IarcConnector, IrisConnector):
            connector = connector_cls(connector_cls.__module__.split(".")[-1], toxicology_yaml)
            assert connector.search("benzene")[0]["carcinogen_class"] == "known"
            assert connector.fetch("71-43-2")["source_reference"] == "PMID:2"
            assert connector.fetch("missing")["source_status"] == "missing_curated_record"

        toxicology_csv = tmp_path / "toxicology.csv"
        toxicology_csv.write_text(
            "chemical,identifier,cas,classification,source_status\n"
            "arsenic,arsenic-id,7440-38-2,Group 1,curated\n",
            encoding="utf-8",
        )
        assert AtsdrConnector("csv", toxicology_csv).fetch("arsenic-id")["classification"] == "Group 1"
        assert AtsdrConnector("missing").search("anything") == []

        comptox_cache = tmp_path / "comptox"
        comptox_cache.mkdir()
        (comptox_cache / "search_benzene.json").write_text(
            json.dumps({"results": [{"dtxsid": "DTXSID3039242"}]}),
            encoding="utf-8",
        )
        (comptox_cache / "DTXSID3039242.json").write_text(
            json.dumps(
                {
                    "dtxsid": "DTXSID3039242",
                    "dtxcid": "DTXCID30182",
                    "casrn": "71-43-2",
                    "preferredName": "Benzene",
                }
            ),
            encoding="utf-8",
        )
        comptox = CompToxConnector(comptox_cache)
        assert comptox.search("benzene")[0]["dtxsid"] == "DTXSID3039242"
        assert comptox.fetch("DTXSID3039242")["preferred_name"] == "Benzene"
        comptox._get_json = lambda _url: (_ for _ in ()).throw(RuntimeError("offline"))
        assert comptox.search("missing") == []
        assert comptox.fetch("DTXSID000")["source_status"] == "missing_cached_record"
        assert CompToxConnector(tmp_path / "dry-comptox", dry_run=True).fetch("DTXSID")[
            "source_status"
        ] == "dry_run"

        pubchem_cache = tmp_path / "pubchem"
        pubchem_cache.mkdir()
        (pubchem_cache / "search_benzene.json").write_text(
            json.dumps({"IdentifierList": {"CID": [241]}}),
            encoding="utf-8",
        )
        (pubchem_cache / "cid_241.json").write_text(
            json.dumps(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 241,
                                "MolecularFormula": "C6H6",
                                "MolecularWeight": 78.11,
                                "CanonicalSMILES": "C1=CC=CC=C1",
                                "InChIKey": "UHOVQNZJYSORNB-UHFFFAOYSA-N",
                            }
                        ]
                    }
                }
            ),
            encoding="utf-8",
        )
        pubchem = PubChemConnector(pubchem_cache)
        assert pubchem.search("benzene")[0]["cid"] == 241
        assert pubchem.fetch("241")["formula"] == "C6H6"
        assert PubChemConnector(tmp_path / "dry-pubchem", dry_run=True).search("benzene") == []
        assert PubChemConnector(tmp_path / "dry-pubchem-fetch", dry_run=True).fetch("241")[
            "source_status"
        ] == "dry_run"

    def test_nhanes_catalog_registry_metadata_is_local(self):
        from ExposoGraph._biomarker_scaffold.scripts.nhanes.biomarker_registry import (
            BIOMARKER_REGISTRY,
        )
        from ExposoGraph._biomarker_scaffold.scripts.nhanes.catalog import (
            available_cycles,
            class_available,
            get_cycle_files,
            get_file_url,
        )
        from ExposoGraph._biomarker_scaffold.scripts.nhanes.class_registry import (
            NHANES_CLASS_REGISTRY,
        )

        assert "2017-2018" in available_cycles()
        assert get_cycle_files("2017-2018")["PAH"] == "PAH_J"
        assert get_file_url("2017-2018", "PAH").endswith("/PAH_J.XPT")
        assert class_available("2017-2018", "PAH") is True
        assert NHANES_CLASS_REGISTRY["PAH"]["requires_creatinine"] is True
        assert BIOMARKER_REGISTRY["PAH"]["URXP10"]["parent_compound"] == "pyrene"
        with pytest.raises(KeyError):
            get_cycle_files("1900-1901")
        with pytest.raises(KeyError):
            get_file_url("2017-2018", "NOT_A_FILE")
        with pytest.raises(KeyError):
            class_available("2017-2018", "NOT_A_CLASS")

    def test_registry_resolution_validation_and_documents_use_fixture_records(
        self, tmp_path, monkeypatch, capsys
    ):
        from ExposoGraph._biomarker_scaffold.scripts.registries import check_mapping
        from ExposoGraph._biomarker_scaffold.scripts.registries.evidence import (
            classify_source_status,
            coverage_report,
        )
        from ExposoGraph._biomarker_scaffold.scripts.registries.loader import (
            load_json_mapping,
            load_registry_document,
            load_yaml_mapping,
            write_json_mapping,
        )
        from ExposoGraph._biomarker_scaffold.scripts.registries.mapping_document import (
            apply_update_list,
            build_mapping_document,
            normalize_mapping_document,
            validate_mapping_document,
        )
        import ExposoGraph._biomarker_scaffold.scripts.registries.loader as registry_loader

        monkeypatch.setattr(
            registry_loader,
            "load_yaml_registry",
            registry_loader.load_yaml_mapping,
            raising=False,
        )
        from ExposoGraph._biomarker_scaffold.scripts.registries.resolver import (
            _load_external_sources,
            _load_measurement_file,
            main as resolver_main,
            normalize_alias,
            resolve_all,
            resolve_biomarker,
        )
        from ExposoGraph._biomarker_scaffold.scripts.registries.schema import (
            BiomarkerRecord,
            EvidenceRecord,
            RegistryValidationReport,
            ValidationIssue,
        )
        from ExposoGraph._biomarker_scaffold.scripts.registries.validator import (
            validate_biomarker_record,
            validate_registry,
        )

        entries = [
            {
                "biomarker": "urinary_1_hydroxypyrene",
                "matrix": "urine",
                "reference_range": [0.1, 5.0],
                "reference_units": "ng/L",
                "source_status": "nhanes_native",
                "lifestyle_factor": "default",
            }
        ]
        document = build_mapping_document(entries, metadata={"schema_version": "1.2.0"})
        assert validate_mapping_document(document) == []
        normalized = normalize_mapping_document(
            {**document, "_update_list": "bad", "extra": {"kept": True}},
            metadata={"customizable": True},
        )
        assert normalized["extra"] == {"kept": True}
        updated = apply_update_list(
            normalized,
            [
                {
                    "entry_id": "urinary_1_hydroxypyrene::default",
                    "entry": {"reference_units": "ug/L"},
                },
                {
                    "entry_id": "blood_benzene::default",
                    "entry": {
                        "biomarker": "blood_benzene",
                        "matrix": "blood",
                        "reference_range": [0, 1],
                        "reference_units": "ng/mL",
                        "source_status": "literature_biomarker",
                    },
                },
                {"op": "remove", "entry_id": "blood_benzene::default"},
                {"op": "noop"},
            ],
        )
        assert updated["entries"][0]["reference_units"] == "ug/L"
        assert validate_mapping_document({"entries": "bad"}) == ["Document must include an 'entries' list."]
        assert validate_mapping_document({"_metadata": {}, "_update_list": [], "entries": [object()]})

        measurement = [
            {
                "biomarker": "urinary_1_hydroxypyrene",
                "nhanes_variable": "URXP10",
                "matrix": "urine",
                "chemical_class": "PAH",
            }
        ]
        mapping = {
            "entries": [
                {
                    "biomarker": "urinary_1_hydroxypyrene",
                    "reference_units": "ng/L",
                    "Km_uM": 12.5,
                },
                {"biomarker": "benzene exposure index"},
            ]
        }
        external_sources = {
            "pubchem": {
                "1_hydroxypyrene": {
                    "pubchem_cid": "123",
                    "references": ["PMID:3"],
                }
            }
        }
        resolved = resolve_biomarker("URXP10", measurement, mapping, external_sources)
        assert resolved["coverage_status"] == "nhanes_native"
        assert resolved["confidence_metadata"]["measurement_match"] == "direct_or_alias"
        assert normalize_alias("Urinary 1-Hydroxypyrene") == "1_hydroxypyrene"
        all_resolved = resolve_all(mapping, measurement, external_sources)
        assert len(all_resolved) == 2
        assert classify_source_status({"biomarker": "blood_ethanol"}) == "scenario_proxy"
        assert classify_source_status({"biomarker": "exposure index"}) == "model_proxy"
        assert classify_source_status({"references": ["PMID:1"]}) == "literature_biomarker"
        report = coverage_report(all_resolved)
        assert report["total_mapping_entries"] == 2
        assert report["entries_with_brenda_kinetic_support"] >= 1

        json_path = tmp_path / "mapping.json"
        yaml_path = tmp_path / "measurement.yaml"
        write_json_mapping(json_path, document)
        yaml_path.write_text(
            "records:\n"
            "  - biomarker: urinary_1_hydroxypyrene\n"
            "    nhanes_variable: URXP10\n",
            encoding="utf-8",
        )
        assert load_json_mapping(json_path)["entries"]
        assert load_yaml_mapping(yaml_path)["records"][0]["biomarker"] == "urinary_1_hydroxypyrene"
        assert load_registry_document(json_path)["entries"]
        assert _load_measurement_file(None) == []
        assert _load_external_sources([str(yaml_path)])["measurement"]["1_hydroxypyrene"]["nhanes_variable"] == "URXP10"
        with pytest.raises(ValueError):
            _load_measurement_file(str(tmp_path / "bad.txt"))

        record = BiomarkerRecord(
            biomarker_id="b1",
            canonical_name="Biomarker",
            matrix="urine",
            chemical_class="PAH",
            source_status="model_proxy",
            evidence=[EvidenceRecord(source="fixture")],
            pubchem_cid="123",
            Km_uM=1.0,
            reference_units="ng/L",
        )
        assert validate_biomarker_record(record)[0].level == "info"
        registry_report = validate_registry(
            [
                record,
                {
                    "biomarker_id": "b1",
                    "matrix": "",
                    "chemical_class": "PAH",
                    "source_status": "",
                },
            ]
        )
        assert registry_report.errors
        info_report = RegistryValidationReport([ValidationIssue("info", None, None, "ok")])
        assert info_report.ok is True
        assert info_report.infos[0].message == "ok"
        with pytest.raises(AttributeError):
            info_report.as_dict()

        invalid_document = {
            "_metadata": {},
            "_update_list": [],
            "entries": [{"biomarker": "x", "lifestyle_factor": "default"}],
        }
        fixed, fixed_errors = check_mapping.validate_biomarker_mapping_document(invalid_document, fix=True)
        assert fixed["entries"][0]["entry_id"] == "x::default"
        assert fixed_errors == ["x::default trace.source_registry is required"]
        _unchanged, errors = check_mapping.validate_biomarker_mapping_document(invalid_document, fix=False)
        assert errors

        out_path = tmp_path / "resolved.json"
        report_path = tmp_path / "coverage.json"
        measurement_path = tmp_path / "measurements.json"
        measurement_path.write_text(json.dumps({"records": measurement}), encoding="utf-8")
        resolver_main(
            [
                "--mapping",
                str(json_path),
                "--nhanes",
                str(measurement_path),
                "--external",
                str(yaml_path),
                "--out",
                str(out_path),
                "--report",
                str(report_path),
            ]
        )
        assert out_path.exists()
        assert report_path.exists()

        check_mapping_path = tmp_path / "check_mapping.json"
        write_json_mapping(
            check_mapping_path,
            {
                "_metadata": {"generated_at": "2026-01-01T00:00:00+00:00"},
                "_update_list": [],
                "entries": [
                    {
                        "biomarker": "urinary_1_hydroxypyrene",
                        "lifestyle_factor": "default",
                        "entry_id": "urinary_1_hydroxypyrene::default",
                        "trace": {"source_registry": "fixture"},
                    }
                ],
            },
        )
        monkeypatch.setattr(
            "sys.argv",
            ["check_mapping", "--mapping", str(check_mapping_path), "--fix"],
        )
        with pytest.raises(SystemExit) as excinfo:
            check_mapping.main()
        assert excinfo.value.code == 0
        assert json.loads(capsys.readouterr().out)["valid"] is True
