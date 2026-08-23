"""Tests for ExposoGraph.exporter."""

import json

from ExposoGraph.config import GraphVisibility
from ExposoGraph.engine import GraphEngine
from ExposoGraph.exporter import (
    export_viewer_bundle,
    parse_graph_artifact,
    parse_graph_data_js,
    parse_graph_data_text,
    parse_graph_html,
    subgraph_to_html_string,
    to_graph_data_js,
    to_interactive_html,
    to_interactive_html_string,
    to_json,
    to_plotly_figure,
    to_plotly_html,
    to_plotly_html_string,
)
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, MatchStatus, Node, NodeType


def _sample_engine() -> GraphEngine:
    engine = GraphEngine()
    kg = KnowledgeGraph(
        nodes=[
            Node(id="BaP", label="Benzo[a]pyrene", type=NodeType.CARCINOGEN),
            Node(id="CYP1A1", label="CYP1A1", type=NodeType.ENZYME),
        ],
        edges=[
            Edge(source="CYP1A1", target="BaP", type=EdgeType.ACTIVATES),
        ],
    )
    engine.load(kg)
    return engine


def _mixed_visibility_engine() -> GraphEngine:
    engine = GraphEngine()
    engine.add_node(
        Node(
            id="CYP1A1",
            label="CYP1A1",
            type=NodeType.ENZYME,
            match_status=MatchStatus.CANONICAL,
        )
    )
    engine.add_node(
        Node(
            id="BaP",
            label="BaP",
            type=NodeType.CARCINOGEN,
            match_status=MatchStatus.ALIAS,
            canonical_id="50-32-8",
            canonical_label="Benzo[a]pyrene",
        )
    )
    engine.add_node(
        Node(
            id="Mystery",
            label="Mystery Chemical",
            type=NodeType.CARCINOGEN,
            match_status=MatchStatus.UNMATCHED,
        )
    )
    engine.add_node(
        Node(
            id="Novel",
            label="Novel Exposure",
            type=NodeType.CARCINOGEN,
            match_status=MatchStatus.CUSTOM,
            custom_type="Exposure",
        )
    )
    engine.add_edge(
        Edge(
            source="CYP1A1",
            target="BaP",
            type=EdgeType.ACTIVATES,
            match_status=MatchStatus.CANONICAL,
        )
    )
    engine.add_edge(
        Edge(
            source="Mystery",
            target="Novel",
            type=EdgeType.CUSTOM,
            match_status=MatchStatus.CUSTOM,
            custom_predicate="CO_OCCURS_WITH",
        )
    )
    return engine


class TestGraphDataJs:
    def test_roundtrip(self, tmp_path):
        engine = _sample_engine()
        out = to_graph_data_js(engine, tmp_path / "graph-data.js")
        assert out.exists()

        content = out.read_text(encoding="utf-8")
        assert "GRAPH_DATA" in content
        assert "Auto-generated" in content

        restored = parse_graph_data_js(out)
        assert len(restored.nodes) == 2
        assert len(restored.edges) == 1

    def test_parse_handles_js_syntax(self, tmp_path):
        js_content = """\
const GRAPH_DATA = {
  nodes: [
    { id: "A", label: "A", type: "Enzyme", detail: "test" },
    { id: "B", label: "B", type: "Metabolite", detail: "test" },
  ],
  edges: [
    { source: "A", target: "B", type: "ACTIVATES" },
  ],
};
"""
        p = tmp_path / "test.js"
        p.write_text(js_content, encoding="utf-8")
        kg = parse_graph_data_js(p)
        assert len(kg.nodes) == 2
        assert len(kg.edges) == 1

    def test_parse_strips_comments(self, tmp_path):
        js_content = """\
// This is a comment
const GRAPH_DATA = {
  nodes: [
    // inline comment
    { id: "X", label: "X", type: "Pathway", detail: "test" },
  ],
  edges: [],
};
"""
        p = tmp_path / "commented.js"
        p.write_text(js_content, encoding="utf-8")
        kg = parse_graph_data_js(p)
        assert len(kg.nodes) == 1

    def test_parse_preserves_urls_inside_strings(self, tmp_path):
        js_content = """\
const GRAPH_DATA = {
  nodes: [
    { id: "X", label: "X", type: "Pathway", detail: "https://example.org/pathway" },
  ],
  edges: [],
};
"""
        p = tmp_path / "url.js"
        p.write_text(js_content, encoding="utf-8")
        kg = parse_graph_data_js(p)
        assert kg.nodes[0].detail == "https://example.org/pathway"


class TestJsonExport:
    def test_export_valid_json(self, tmp_path):
        engine = _sample_engine()
        out = to_json(engine, tmp_path / "kg.json")
        data = json.loads(out.read_text(encoding="utf-8"))
        assert len(data["nodes"]) == 2
        assert len(data["edges"]) == 1

    def test_export_respects_visibility(self, tmp_path):
        engine = _mixed_visibility_engine()

        out = to_json(
            engine,
            tmp_path / "kg_validated.json",
            visibility=GraphVisibility.VALIDATED_ONLY,
        )
        data = json.loads(out.read_text(encoding="utf-8"))

        assert {node["id"] for node in data["nodes"]} == {"CYP1A1", "BaP"}
        assert len(data["edges"]) == 1


class TestNoneCleanup:
    def test_none_values_stripped(self, tmp_path):
        engine = _sample_engine()
        out = to_json(engine, tmp_path / "kg.json")
        data = json.loads(out.read_text(encoding="utf-8"))
        for node in data["nodes"]:
            assert None not in node.values()


class TestViewerBundle:
    def test_export_viewer_bundle_copies_index_template(self, tmp_path):
        engine = _sample_engine()
        template_dir = tmp_path / "template"
        template_dir.mkdir()
        (template_dir / "index.html").write_text("<html>viewer</html>", encoding="utf-8")

        out_dir = export_viewer_bundle(
            engine,
            tmp_path / "bundle",
            template_dir=template_dir,
        )

        assert out_dir.exists()
        assert (out_dir / "index.html").read_text(encoding="utf-8") == "<html>viewer</html>"
        assert (out_dir / "graph-data.js").exists()

    def test_export_viewer_bundle_writes_builtin_template_when_missing(self, tmp_path):
        engine = _sample_engine()

        out_dir = export_viewer_bundle(engine, tmp_path / "bundle")

        assert out_dir.exists()
        assert (out_dir / "index.html").exists()
        assert "ExposoGraph Graph Viewer" in (out_dir / "index.html").read_text(encoding="utf-8")
        assert (out_dir / "graph-data.js").exists()


class TestInteractiveHtml:
    def test_html_export_uses_builtin_template_when_missing(self, tmp_path):
        engine = _sample_engine()

        out = to_interactive_html(engine, tmp_path / "graph.html")

        html = out.read_text(encoding="utf-8")
        assert "GRAPH_DATA" in html
        assert "ExposoGraph Graph Viewer" in html

    def test_html_export_inlines_graph_data(self, tmp_path):
        engine = _sample_engine()
        template = tmp_path / "template.html"
        template.write_text(
            '<html><head><script src="./graph-data.js"></script></head><body></body></html>',
            encoding="utf-8",
        )

        html = to_interactive_html_string(engine, template_path=template)

        assert "GRAPH_DATA" in html
        assert './graph-data.js' not in html

    def test_html_roundtrip(self, tmp_path):
        engine = _sample_engine()
        template = tmp_path / "template.html"
        template.write_text(
            '<html><head><script src="./graph-data.js"></script></head><body></body></html>',
            encoding="utf-8",
        )

        out = to_interactive_html(engine, tmp_path / "graph.html", template_path=template)

        restored = parse_graph_html(out)
        assert len(restored.nodes) == 2
        assert len(restored.edges) == 1

    def test_parse_graph_artifact_supports_html(self, tmp_path):
        engine = _sample_engine()
        template = tmp_path / "template.html"
        template.write_text(
            '<html><head><script src="./graph-data.js"></script></head><body></body></html>',
            encoding="utf-8",
        )
        out = to_interactive_html(engine, tmp_path / "graph.html", template_path=template)

        restored = parse_graph_artifact(out)
        assert len(restored.nodes) == 2
        assert len(restored.edges) == 1

    def test_html_export_respects_visibility(self, tmp_path):
        engine = _mixed_visibility_engine()
        template = tmp_path / "template.html"
        template.write_text(
            '<html><head><script src="./graph-data.js"></script></head><body></body></html>',
            encoding="utf-8",
        )

        html = to_interactive_html_string(
            engine,
            template_path=template,
            visibility=GraphVisibility.EXPLORATORY_ONLY,
        )
        parsed = parse_graph_html(
            to_interactive_html(
                engine,
                tmp_path / "graph.html",
                template_path=template,
                visibility=GraphVisibility.EXPLORATORY_ONLY,
            )
        )

        assert "GRAPH_DATA" in html
        assert {node.id for node in parsed.nodes} == {"Mystery", "Novel"}
        assert len(parsed.edges) == 1


class TestSubgraphToHtmlString:
    """subgraph_to_html_string: the dict-based analog of to_interactive_html_string.

    Used by ui_map_viewer.py to render a filtered {"nodes": ..., "edges": ...}
    dict (from GraphEngine.map_viewer_subgraph/dim_by_tissue_threshold)
    without first assembling a throwaway GraphEngine.
    """

    def test_inlines_arbitrary_dict_as_graph_data(self, tmp_path):
        template = tmp_path / "template.html"
        template.write_text(
            '<html><head><script src="./graph-data.js"></script></head><body></body></html>',
            encoding="utf-8",
        )
        graph = {
            "nodes": [
                {"id": "BaP", "label": "Benzo[a]pyrene", "type": "Carcinogen"},
                {"id": "CYP1A1", "label": "CYP1A1", "type": "Enzyme", "_dimmed": False},
            ],
            "edges": [
                {"source": "BaP", "target": "CYP1A1", "type": "ACTIVATES", "_dimmed": False},
            ],
        }

        html = subgraph_to_html_string(graph, template_path=template)

        assert "GRAPH_DATA" in html
        assert './graph-data.js' not in html
        restored = parse_graph_data_text(html)
        assert {node.id for node in restored.nodes} == {"BaP", "CYP1A1"}
        assert len(restored.edges) == 1

    def test_strips_none_values(self, tmp_path):
        template = tmp_path / "template.html"
        template.write_text(
            '<html><head><script src="./graph-data.js"></script></head><body></body></html>',
            encoding="utf-8",
        )
        graph = {
            "nodes": [{"id": "X", "label": "X", "type": "Enzyme", "group": None}],
            "edges": [],
        }

        html = subgraph_to_html_string(graph, template_path=template)

        start = html.index("const GRAPH_DATA")
        end = html.index(";", start)
        data = json.loads(html[start:end].split("=", 1)[1].strip())
        assert "group" not in data["nodes"][0]

    def test_uses_builtin_template_when_missing(self):
        html = subgraph_to_html_string({"nodes": [], "edges": []})
        assert "GRAPH_DATA" in html


class TestPlotlyExport:
    def test_plotly_figure_contains_graph_content(self):
        engine = _sample_engine()

        fig = to_plotly_figure(engine)
        payload = json.dumps(fig.to_plotly_json())

        assert "Benzo[a]pyrene" in payload
        assert "CYP1A1" in payload
        assert "ACTIVATES" in payload

    def test_plotly_export_respects_visibility(self, tmp_path):
        engine = _mixed_visibility_engine()

        html = to_plotly_html_string(
            engine,
            visibility=GraphVisibility.VALIDATED_ONLY,
        )

        assert "CYP1A1" in html
        assert "Benzo[a]pyrene" in html
        assert "Mystery Chemical" not in html
        assert "Novel Exposure" not in html

        out = to_plotly_html(
            engine,
            tmp_path / "graph_plotly.html",
            visibility=GraphVisibility.EXPLORATORY_ONLY,
        )
        exported = out.read_text(encoding="utf-8")

        assert "Mystery Chemical" in exported
        assert "Novel Exposure" in exported
        assert "Benzo[a]pyrene" not in exported
