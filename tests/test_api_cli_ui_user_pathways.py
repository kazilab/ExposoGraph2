from __future__ import annotations

import importlib
import json
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import pytest


def test_module3_api_example_runs():
    from ExposoGraph.flux_engine import CarcinogenClass, compute_pathway_flux

    result = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )
    payload = asdict(result)

    assert payload["carcinogen_class"] == "PAH"
    assert payload["tissue"] == "Lung"
    assert payload["total_activation"] >= 0
    assert payload["total_detox"] >= 0
    assert "synergy_matrix" not in payload
    assert "interaction_factor" not in payload
    json.dumps(payload, default=str, allow_nan=False)


def test_module5_api_example_runs():
    from ExposoGraph.interaction_engine import compute_interaction_matrix
    from ExposoGraph.interaction_engine import _interaction_matrix_to_compat_dict

    result = compute_interaction_matrix(
        {"PAH": 3.0, "NNK": 4.0, "benzene": 5.0, "ethanol": 3.0},
        genotypes={"GSTM1": "null", "CYP2E1": "NM"},
        lifestyle={"smoking": True, "alcohol_moderate": True},
        tissue="Liver",
    )
    payload = _interaction_matrix_to_compat_dict(result)

    assert payload["synergy_matrix"]
    assert payload["interaction_factor"] > 0
    assert payload["mechanism_attribution"]["decomposition_basis"] == "eight_state_shapley"
    assert payload["module5_model_card"]["diagnostic_output_policy"]
    json.dumps(payload, default=str, allow_nan=False)


def test_cli_imports_or_runs_if_advertised():
    for module_name in ("flux_cli", "exposure_cli", "interaction_cli"):
        completed = subprocess.run(
            [sys.executable, "-m", f"ExposoGraph.{module_name}", "--help"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        output = f"{completed.stdout}\n{completed.stderr}".lower()
        assert completed.returncode == 0, output
        assert "usage" in output


def test_streamlit_imports_without_starting_server():
    pytest.importorskip("streamlit")

    for module_name in (
        "ExposoGraph.ui_flux",
        "ExposoGraph.ui_map_viewer",
        "ExposoGraph.ui_d3_viewer",
        "ExposoGraph.ui_data",
        "ExposoGraph.ui_manual",
        "ExposoGraph.ui_sidebar",
    ):
        importlib.import_module(module_name)


def test_ui_labels_do_not_deprecate_module3():
    repo_root = Path(__file__).resolve().parents[1]
    label_files = [
        repo_root / "README.md",
        repo_root / "ExposoGraph" / "app.py",
        repo_root / "ExposoGraph" / "ui_flux.py",
        repo_root / "ExposoGraph" / "ui_map_viewer.py",
        repo_root / "ExposoGraph" / "ui_d3_viewer.py",
    ]
    label_text = "\n".join(path.read_text(encoding="utf-8") for path in label_files)
    lower = label_text.lower()

    assert "flux engine" in lower
    assert "compute_pathway_flux" in lower
    assert "compute_interaction_matrix" in lower
    assert "module 3 is deprecated" not in lower
    assert "module 3 deprecated" not in lower
    assert "module 5 replaces module 3" not in lower
    assert "module 5 supersedes module 3" not in lower


def test_graph_data_validation(tmp_path):
    from ExposoGraph.reference_data import build_reference_engine
    from ExposoGraph.config import GraphVisibility
    from ExposoGraph.exporter import parse_graph_data_js, to_graph_data_js
    from ExposoGraph.graph_filters import graph_visibility_label

    output_path = tmp_path / "graph-data.js"
    to_graph_data_js(
        build_reference_engine(),
        output_path,
        visibility=GraphVisibility.VALIDATED_ONLY,
    )
    graph = parse_graph_data_js(output_path)

    assert graph.nodes
    assert graph.edges
    assert graph_visibility_label(GraphVisibility.VALIDATED_ONLY) == "Validated Only"
    node_ids = {node.id for node in graph.nodes}
    assert all(edge.source in node_ids and edge.target in node_ids for edge in graph.edges)
