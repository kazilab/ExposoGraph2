"""Tests for top-level ExposoGraph public API exports."""

from ExposoGraph import (
    CURATION_SOURCE_MANIFEST,
    REFERENCE_KEGG_PATHWAYS,
    ArchitectureInventoryGroup,
    ArchitectureSummary,
    CytoscapeBundle,
    GraphMode,
    GraphRepository,
    GraphVisibility,
    ViewerLayoutMode,
    apply_viewer_filters,
    build_androgen_module_engine,
    build_androgen_module_graph,
    build_architecture_figure_data,
    build_cytoscape_bundle,
    build_detail_payload,
    build_full_legends_architecture_summary,
    build_full_legends_engine,
    build_full_legends_graph,
    build_reference_architecture_summary,
    build_reference_engine,
    build_reference_graph,
    compute_viewer_positions,
    create_dash_viewer_app,
    filter_knowledge_graph,
    graph_visibility_label,
    launch_dash_viewer,
    load_cytoscape_bundle,
    load_viewer_positions,
    normalize_graph_mode,
    normalize_graph_visibility,
    normalize_viewer_layout_mode,
    paper_architecture_overrides,
    prepare_knowledge_graph,
    render_architecture_figure,
    save_architecture_figure,
    to_plotly_figure,
    to_plotly_html,
    to_plotly_html_string,
    write_cytoscape_bundle,
    write_full_legends_exports,
    write_reference_exports,
    write_viewer_positions,
)


def test_top_level_mode_exports_are_available():
    assert GraphMode.STRICT.value == "strict"
    assert GraphVisibility.VALIDATED_ONLY.value == "validated_only"
    assert normalize_graph_mode("validated") == GraphMode.STRICT
    assert normalize_graph_visibility("provisional") == GraphVisibility.EXPLORATORY_ONLY
    assert ViewerLayoutMode.PRESET.value == "preset"
    assert normalize_viewer_layout_mode("saved") == ViewerLayoutMode.PRESET


def test_top_level_helper_exports_are_available():
    assert ArchitectureInventoryGroup is not None
    assert ArchitectureSummary is not None
    assert isinstance(CURATION_SOURCE_MANIFEST, dict)
    assert isinstance(REFERENCE_KEGG_PATHWAYS, list)
    assert CytoscapeBundle is not None
    assert callable(filter_knowledge_graph)
    assert callable(graph_visibility_label)
    assert callable(prepare_knowledge_graph)
    assert callable(GraphRepository)
    assert callable(build_cytoscape_bundle)
    assert callable(compute_viewer_positions)
    assert callable(write_cytoscape_bundle)
    assert callable(load_cytoscape_bundle)
    assert callable(write_viewer_positions)
    assert callable(load_viewer_positions)
    assert callable(apply_viewer_filters)
    assert callable(build_detail_payload)
    assert callable(build_androgen_module_graph)
    assert callable(build_androgen_module_engine)
    assert callable(build_full_legends_architecture_summary)
    assert callable(build_full_legends_graph)
    assert callable(build_full_legends_engine)
    assert callable(build_reference_architecture_summary)
    assert callable(build_reference_graph)
    assert callable(build_reference_engine)
    assert callable(create_dash_viewer_app)
    assert callable(launch_dash_viewer)
    assert callable(build_architecture_figure_data)
    assert callable(paper_architecture_overrides)
    assert callable(render_architecture_figure)
    assert callable(save_architecture_figure)
    assert callable(to_plotly_figure)
    assert callable(to_plotly_html)
    assert callable(to_plotly_html_string)
    assert callable(write_reference_exports)
    assert callable(write_full_legends_exports)
