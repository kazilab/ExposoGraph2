"""Smoke tests for tissue-specific subgraph extraction."""

from ExposoGraph import (
    DEFAULT_THRESHOLD_SWEEP,
    ThresholdSweepResult,
    TissueThresholdRetention,
    TissueType,
    build_reference_graph,
    get_tissue_expression,
    get_top_carcinogen_classes_for_tissue,
    tissue_threshold_sweep,
)


def test_liver_top_carcinogen_classes_list_is_populated():
    top = get_top_carcinogen_classes_for_tissue("Liver")

    assert len(top) >= 3
    assert all(isinstance(item, tuple) and len(item) == 2 for item in top)
    # Liver handles PAHs prominently via CYP1A1/3A4.
    assert any("PAH" in name for name, _ in top)


def test_tissue_expression_lookup_returns_mapping_for_liver():
    expr = get_tissue_expression(TissueType.LIVER)

    assert isinstance(expr, dict)
    assert all(isinstance(v, (int, float)) for v in expr.values())


def test_default_threshold_sweep_is_set():
    assert DEFAULT_THRESHOLD_SWEEP == (0.10, 0.25, 0.50)


def test_tissue_threshold_sweep_returns_retention_entries():
    graph = build_reference_graph()
    graph_data = {
        "nodes": [node.model_dump() for node in graph.nodes],
        "edges": [edge.model_dump() for edge in graph.edges],
    }
    result = tissue_threshold_sweep(graph_data, tissues=["Liver", "Lung"])
    assert isinstance(result, ThresholdSweepResult)
    assert set(result.thresholds) == set(DEFAULT_THRESHOLD_SWEEP)
    assert result.entries
    assert all(isinstance(e, TissueThresholdRetention) for e in result.entries)
    # Every tissue has one entry per threshold
    assert len(result.entries) == 2 * len(DEFAULT_THRESHOLD_SWEEP)


def test_tissue_threshold_sweep_retention_decreases_with_threshold():
    graph = build_reference_graph()
    graph_data = {
        "nodes": [node.model_dump() for node in graph.nodes],
        "edges": [edge.model_dump() for edge in graph.edges],
    }
    result = tissue_threshold_sweep(graph_data, tissues=["Liver"])
    liver_matrix = result.edge_retention_matrix["Liver"]
    # Higher thresholds retain the same number of edges or fewer.
    values = [liver_matrix[t] for t in sorted(liver_matrix.keys())]
    assert all(values[i] >= values[i + 1] for i in range(len(values) - 1))
