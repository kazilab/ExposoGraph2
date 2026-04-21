"""Tests for ExposoGraph.storage."""

import sqlite3

from ExposoGraph.config import GraphVisibility
from ExposoGraph.engine import GraphEngine
from ExposoGraph.exporter import parse_graph_data_text
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, MatchStatus, Node, NodeType
from ExposoGraph.storage import GraphRepository


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
    engine.add_edge(
        Edge(
            source="CYP1A1",
            target="BaP",
            type=EdgeType.ACTIVATES,
            match_status=MatchStatus.CANONICAL,
        )
    )
    return engine


def _template_path(tmp_path):
    template = tmp_path / "template.html"
    template.write_text(
        '<html><head><script src="./graph-data.js"></script></head><body></body></html>',
        encoding="utf-8",
    )
    return template


class TestGraphRepository:
    def test_save_and_load_latest_revision(self, tmp_path):
        repo = GraphRepository(tmp_path / "graphs.sqlite3")
        engine = _sample_engine()

        saved = repo.save_engine(
            graph_key="bap_demo",
            graph_name="BaP Demo",
            engine=engine,
            template_path=_template_path(tmp_path),
            note="initial",
        )

        assert saved.revision_number == 1
        latest = repo.get_latest_revision("bap_demo")
        assert latest is not None
        assert latest.revision_number == 1
        assert "GRAPH_DATA" in latest.html
        assert len(latest.to_knowledge_graph().nodes) == 2

    def test_save_creates_incrementing_revisions(self, tmp_path):
        repo = GraphRepository(tmp_path / "graphs.sqlite3")
        engine = _sample_engine()
        template = _template_path(tmp_path)

        first = repo.save_engine(
            graph_key="bap_demo",
            graph_name="BaP Demo",
            engine=engine,
            template_path=template,
            note="initial",
        )
        engine.add_node(Node(id="EPHX1", label="EPHX1", type=NodeType.ENZYME))
        second = repo.save_engine(
            graph_key="bap_demo",
            graph_name="BaP Demo",
            engine=engine,
            template_path=template,
            note="added EPHX1",
        )

        assert first.revision_number == 1
        assert second.revision_number == 2

        revisions = repo.list_revisions("bap_demo")
        assert [rev.revision_number for rev in revisions] == [2, 1]
        assert revisions[0].node_count == 3

    def test_list_graphs_returns_latest_revision_per_graph(self, tmp_path):
        repo = GraphRepository(tmp_path / "graphs.sqlite3")
        template = _template_path(tmp_path)

        engine = _sample_engine()
        repo.save_engine(
            graph_key="graph_a",
            graph_name="Graph A",
            engine=engine,
            template_path=template,
        )
        engine.add_node(Node(id="XPC", label="XPC", type=NodeType.ENZYME))
        repo.save_engine(
            graph_key="graph_a",
            graph_name="Graph A",
            engine=engine,
            template_path=template,
        )
        repo.save_engine(
            graph_key="graph_b",
            graph_name="Graph B",
            engine=_sample_engine(),
            template_path=template,
        )

        graphs = repo.list_graphs()
        keys = {graph.graph_key for graph in graphs}
        assert keys == {"graph_a", "graph_b"}
        graph_a = next(graph for graph in graphs if graph.graph_key == "graph_a")
        assert graph_a.revision_number == 2

    def test_save_engine_respects_visibility(self, tmp_path):
        repo = GraphRepository(tmp_path / "graphs.sqlite3")
        template = _template_path(tmp_path)

        saved = repo.save_engine(
            graph_key="filtered_demo",
            graph_name="Filtered Demo",
            engine=_mixed_visibility_engine(),
            template_path=template,
            visibility=GraphVisibility.VALIDATED_ONLY,
        )

        latest = repo.get_latest_revision("filtered_demo")
        assert latest is not None
        assert saved.visibility == GraphVisibility.VALIDATED_ONLY
        assert latest.visibility == GraphVisibility.VALIDATED_ONLY
        assert {node.id for node in latest.to_knowledge_graph().nodes} == {"CYP1A1", "BaP"}
        assert len(latest.to_knowledge_graph().edges) == 1
        html_graph = parse_graph_data_text(latest.html)
        assert {node.id for node in html_graph.nodes} == {"CYP1A1", "BaP"}

    def test_context_manager_closes_connection(self, tmp_path):
        with GraphRepository(tmp_path / "graphs.sqlite3") as repo:
            repo.list_graphs()
            assert repo._conn is not None

        assert repo._conn is None

    def test_existing_schema_migrates_visibility_column(self, tmp_path):
        db_path = tmp_path / "legacy.sqlite3"
        conn = sqlite3.connect(db_path)
        conn.executescript(
            """
            CREATE TABLE graphs (
                graph_key TEXT PRIMARY KEY,
                graph_name TEXT NOT NULL,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            );
            CREATE TABLE graph_revisions (
                revision_id INTEGER PRIMARY KEY AUTOINCREMENT,
                graph_key TEXT NOT NULL,
                graph_name TEXT NOT NULL,
                revision_number INTEGER NOT NULL,
                created_at TEXT NOT NULL,
                note TEXT,
                node_count INTEGER NOT NULL,
                edge_count INTEGER NOT NULL,
                graph_json TEXT NOT NULL,
                html TEXT NOT NULL
            );
            """
        )
        conn.close()

        repo = GraphRepository(db_path)
        column_names = {
            row["name"]
            for row in repo.connection.execute("PRAGMA table_info(graph_revisions)").fetchall()
        }

        assert "visibility" in column_names
