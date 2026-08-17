"""Focused WP06 checks for v2 local data/provider paths."""

from __future__ import annotations

import json
import socket
import subprocess
import sys
import urllib.request
from pathlib import Path

import pytest

import ExposoGraph
from ExposoGraph import (
    CarcinogenClass,
    InteractionMatrixResult,
    PathwayFluxResult,
    build_reference_engine,
    build_reference_graph,
    compute_interaction_matrix,
    compute_pathway_flux,
    exposure_engine,
    flux_engine,
    interaction_engine,
)
from ExposoGraph.exporter import parse_graph_data_js, to_graph_data_js
from ExposoGraph.models import NodeType
from ExposoGraph.parameter_provider import JSONInteractionParameterProvider

PACKAGE_ROOT = Path(ExposoGraph.__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent
DATA_DIR = PACKAGE_ROOT / "data"
MAP_DIR = PACKAGE_ROOT / "map"

RUNTIME_JSON_FILES = (
    DATA_DIR / "kinetic_parameters.json",
    DATA_DIR / "interaction_parameters.json",
    DATA_DIR / "exposure_database.json",
    DATA_DIR / "proxy_flux_parameters.json",
    DATA_DIR / "proxy_flux_provenance.json",
    DATA_DIR / "parameter_provenance.json",
    DATA_DIR / "tissue_expression_data.json",
    DATA_DIR / "biomarker_mapping.json",
    DATA_DIR / "mutational_signatures.json",
)

RUNTIME_SOURCE_FILES = (
    PACKAGE_ROOT / "__init__.py",
    PACKAGE_ROOT / "reference_data.py",
    PACKAGE_ROOT / "flux_engine.py",
    PACKAGE_ROOT / "exposure_engine.py",
    PACKAGE_ROOT / "interaction_engine.py",
    PACKAGE_ROOT / "parameter_provider.py",
    PACKAGE_ROOT / "parameter_resolution.py",
    PACKAGE_ROOT / "kinetic_resolver.py",
    PACKAGE_ROOT / "tissue_subgraphs.py",
    PACKAGE_ROOT / "expanded_metals.py",
    PACKAGE_ROOT / "unified_api.py",
    PACKAGE_ROOT / "exporter.py",
    PACKAGE_ROOT / "biomarker_mapping.py",
    PACKAGE_ROOT / "mutational_signatures.py",
)


def _assert_under(path: Path, parent: Path) -> None:
    path.resolve().relative_to(parent.resolve())


def test_packaged_data_files_exist() -> None:
    for path in RUNTIME_JSON_FILES:
        assert path.is_file(), f"Missing packaged data file: {path.name}"
        _assert_under(path, DATA_DIR)
        assert isinstance(json.loads(path.read_text(encoding="utf-8")), dict)

    assert (MAP_DIR / "graph-data.js").is_file()
    assert (MAP_DIR / "index.html").is_file()
    _assert_under(MAP_DIR / "graph-data.js", MAP_DIR)
    _assert_under(MAP_DIR / "index.html", MAP_DIR)


def test_module3_data_loads_from_packaged_paths() -> None:
    for path in (
        flux_engine._KINETIC_PARAMS_FILE,
        flux_engine._EXPOSURE_DB_FILE,
        flux_engine._INTERACTION_PARAMS_FILE,
        flux_engine._PROXY_FLUX_PARAMS_FILE,
        flux_engine._PROXY_FLUX_PROVENANCE_FILE,
        exposure_engine._EXPOSURE_DB_FILE,
    ):
        assert Path(path).is_file()
        _assert_under(Path(path), DATA_DIR)

    assert flux_engine._load_kinetic_params()["carcinogen_classes"]
    assert flux_engine._load_exposure_db()["carcinogen_classes"]
    assert flux_engine._load_proxy_flux_params()["classes"]
    assert flux_engine._load_proxy_flux_provenance()

    result = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )
    assert isinstance(result, PathwayFluxResult)
    assert result.total_activation > 0.0
    assert result.total_detox > 0.0


def test_module5_data_loads_from_packaged_paths() -> None:
    for path in (
        interaction_engine._INTERACTION_PARAMS_FILE,
        interaction_engine._PROVENANCE_FILE,
    ):
        assert Path(path).is_file()
        _assert_under(Path(path), DATA_DIR)

    params = interaction_engine._load_interaction_params()
    provenance = interaction_engine.get_parameter_provenance()
    assert params["competitive_inhibition"]
    assert params["gsh_depletion"]
    assert provenance["pairs"]

    provider = JSONInteractionParameterProvider()
    _assert_under(provider.interaction_path, DATA_DIR)
    _assert_under(provider.provenance_path, DATA_DIR)
    assert provider.get_competitive_interactions("CYP2E1")
    assert provider.get_gsh_consumers()

    result = compute_interaction_matrix({"PAH": 1.0, "HCA": 0.5}, tissue="Liver")
    assert isinstance(result, InteractionMatrixResult)
    assert result.total_independent_risk > 0.0
    assert result.total_interaction_risk > 0.0


def test_reference_graph_builds_without_remote_access(tmp_path: Path) -> None:
    graph_path = MAP_DIR / "graph-data.js"
    graph = parse_graph_data_js(graph_path)
    assert len(graph.nodes) > 0
    assert len(graph.edges) > 0

    reference_graph = build_reference_graph()
    # graph-data.json (the Python-side source of truth) now intentionally
    # carries additional NodeType.SUBSTRATE nodes that graph-data.js (the
    # bundled Streamlit/D3 viewer asset) does not -- these are queryable via
    # the engine but are not meant to render in the static map bundle. Node
    # counts diverge by exactly that many; edges are untouched and still
    # match exactly.
    graph_node_ids = {node.id for node in graph.nodes}
    reference_node_ids = {node.id for node in reference_graph.nodes}
    assert graph_node_ids.issubset(reference_node_ids)
    substrate_only_ids = {
        node.id for node in reference_graph.nodes if node.type == NodeType.SUBSTRATE
    }
    assert reference_node_ids - graph_node_ids == substrate_only_ids
    # graph-data.json also now carries the 58 enzyme->substrate topology
    # edges added for the competitive_inhibition pairs from
    # interaction_parameters.json that had no prior qualifying edge --
    # these are queryable via the engine but, like the Substrate nodes
    # above, are not yet rendered in the static map bundle (a render-
    # exclusion pass for both is a separate, later step). Edge counts
    # diverge by exactly that many; all pre-existing edges still match.
    graph_edge_pairs = {(edge.source, edge.target) for edge in graph.edges}
    reference_edge_pairs = {(edge.source, edge.target) for edge in reference_graph.edges}
    assert graph_edge_pairs.issubset(reference_edge_pairs)
    assert len(reference_edge_pairs) - len(graph_edge_pairs) == 58

    engine = build_reference_engine()
    assert engine.node_count == len(reference_graph.nodes)
    assert engine.edge_count == len(reference_graph.edges)
    assert engine.validate() == []

    out = to_graph_data_js(engine, tmp_path / "graph-data.js")
    assert out.is_file()
    restored = parse_graph_data_js(out)
    assert len(restored.nodes) == engine.node_count
    assert len(restored.edges) == engine.edge_count


def test_no_absolute_protected_repo_paths() -> None:
    forbidden_fragments = (
        "ExposoGraph2-main",
        "Repo\\ExposoGraph2-main",
        "Repo/ExposoGraph2-main",
    )
    for path in RUNTIME_SOURCE_FILES:
        text = path.read_text(encoding="utf-8")
        for fragment in forbidden_fragments:
            assert fragment not in text, f"Forbidden protected-path fragment in {path.name}"

    for path in (*RUNTIME_JSON_FILES, MAP_DIR / "graph-data.js", MAP_DIR / "index.html"):
        text = path.read_text(encoding="utf-8")
        for fragment in forbidden_fragments:
            assert fragment not in text, f"Forbidden protected-path fragment in {path.name}"


def test_no_missing_json_runtime_dependency() -> None:
    graph = build_reference_graph()
    assert graph.nodes
    assert graph.edges

    assert flux_engine._load_kinetic_params()
    assert flux_engine._load_exposure_db()
    assert flux_engine._load_proxy_flux_params()
    assert interaction_engine._load_interaction_params()
    assert interaction_engine.get_parameter_provenance()


def test_no_remote_dependency_for_v2_runtime(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_network(*args: object, **kwargs: object) -> None:
        raise AssertionError("Network access is not required by v2 runtime data paths")

    monkeypatch.setattr(socket, "create_connection", fail_network)
    monkeypatch.setattr(urllib.request, "urlopen", fail_network)

    try:
        import requests
    except ImportError:
        requests = None
    if requests is not None:
        monkeypatch.setattr(requests.sessions.Session, "request", fail_network)

    assert build_reference_graph().nodes
    assert build_reference_engine().node_count > 0
    assert compute_pathway_flux("PAH", {"CYP1A1": "NM", "GSTM1": "NM"}, tissue="Lung")
    assert compute_interaction_matrix({"PAH": 1.0, "HCA": 0.5}, tissue="Liver")


def test_platform_integration_harness_passes() -> None:
    completed = subprocess.run(
        [sys.executable, "tools/local_platform_integration_check.py"],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        capture_output=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Full-platform integration check: PASS" in completed.stdout
    assert "Repository metadata check" in completed.stdout
    assert "Runtime path disclosure" in completed.stdout
