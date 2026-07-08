from __future__ import annotations

import json
from pathlib import Path

import ExposoGraph
from ExposoGraph import (
    InteractionMatrixResult,
    LocalV2DataProvider,
    PathwayFluxResult,
    get_default_v2_provider,
)


PACKAGE_ROOT = Path(ExposoGraph.__file__).resolve().parent


def _assert_under(path: Path, parent: Path) -> None:
    path.resolve().relative_to(parent.resolve())


def test_provider_facade_exposes_required_v2_accessors() -> None:
    provider = get_default_v2_provider()

    assert isinstance(provider, LocalV2DataProvider)
    assert provider.package_root == PACKAGE_ROOT
    assert provider.data_dir == PACKAGE_ROOT / "data"
    assert provider.map_dir == PACKAGE_ROOT / "map"
    assert provider.reference_graph_path == PACKAGE_ROOT / "map" / "graph-data.js"
    assert {
        "kinetic_parameters",
        "interaction_parameters",
        "exposure_database",
        "parameter_provenance",
        "reference_graph",
    } <= set(provider.runtime_paths())


def test_module3_uses_provider_or_provider_wrapped_path() -> None:
    provider = get_default_v2_provider()
    result = provider.compute_module3_pathway_flux(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )

    assert isinstance(result, PathwayFluxResult)
    assert result.total_activation > 0.0
    assert result.total_detox > 0.0


def test_module5_uses_provider_or_provider_wrapped_path() -> None:
    provider = get_default_v2_provider()
    result = provider.compute_module5_interaction_matrix(
        {"PAH": 1.0, "HCA": 0.5},
        tissue="Liver",
        include_biological_outputs=True,
    )

    assert isinstance(result, InteractionMatrixResult)
    assert result.mechanism_resolved_risks
    assert result.total_interaction_risk > 0.0


def test_reference_graph_provider_path_local_only() -> None:
    provider = get_default_v2_provider()
    paths = provider.runtime_paths()

    for path in paths.values():
        _assert_under(path, PACKAGE_ROOT)
        assert path.exists()

    graph = provider.build_reference_graph()
    engine = provider.build_reference_engine()

    assert graph.nodes
    assert graph.edges
    assert engine.node_count == len(graph.nodes)


def test_provider_outputs_are_json_safe_where_applicable() -> None:
    provider = get_default_v2_provider()
    module3 = provider.compute_module3_pathway_flux(
        "PAH",
        {"CYP1A1": "NM", "GSTM1": "NM"},
        tissue="Lung",
    )
    module5 = provider.compute_module5_interaction_matrix(
        {"benzene": 1.0, "ethanol": 10.0},
        include_biological_outputs=True,
    )

    json.dumps(provider.json_safe(module3), sort_keys=True, allow_nan=False)
    json.dumps(provider.json_safe(module5), sort_keys=True, allow_nan=False)

