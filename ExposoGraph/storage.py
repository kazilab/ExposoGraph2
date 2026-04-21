"""SQLite-backed graph storage with revision history."""

from __future__ import annotations

import json
import sqlite3
from contextlib import suppress
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import TracebackType

from .config import GraphVisibility, normalize_graph_visibility
from .engine import GraphEngine
from .exporter import to_interactive_html_string
from .graph_filters import filter_knowledge_graph
from .models import KnowledgeGraph


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class GraphRevisionSummary:
    revision_id: int
    graph_key: str
    graph_name: str
    revision_number: int
    created_at: str
    node_count: int
    edge_count: int
    note: str | None = None
    visibility: GraphVisibility = GraphVisibility.ALL


@dataclass(frozen=True)
class GraphRevision(GraphRevisionSummary):
    graph_json: str = ""
    html: str = ""

    def to_knowledge_graph(self) -> KnowledgeGraph:
        return KnowledgeGraph(**json.loads(self.graph_json))


class GraphRepository:
    """Persist graphs and their revisions in a local SQLite database."""

    def __init__(self, db_path: str | Path) -> None:
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn: sqlite3.Connection | None = self._create_connection()
        self._initialize()

    def _create_connection(self) -> sqlite3.Connection:
        conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        return conn

    @property
    def connection(self) -> sqlite3.Connection:
        """Return the persistent connection, reconnecting if closed."""
        if self._conn is None:
            self._conn = self._create_connection()
            return self._conn
        try:
            self._conn.execute("SELECT 1")
        except (sqlite3.ProgrammingError, sqlite3.OperationalError):
            self._conn = self._create_connection()
        return self._conn

    def close(self) -> None:
        """Close the persistent connection."""
        conn, self._conn = self._conn, None
        if conn is None:
            return
        with suppress(sqlite3.ProgrammingError):
            conn.close()

    def __enter__(self) -> GraphRepository:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        self.close()

    def __del__(self) -> None:
        with suppress(Exception):
            self.close()

    def _initialize(self) -> None:
        with self.connection as conn:
            conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS graphs (
                    graph_key TEXT PRIMARY KEY,
                    graph_name TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS graph_revisions (
                    revision_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    graph_key TEXT NOT NULL,
                    graph_name TEXT NOT NULL,
                    visibility TEXT NOT NULL DEFAULT 'all',
                    revision_number INTEGER NOT NULL,
                    created_at TEXT NOT NULL,
                    note TEXT,
                    node_count INTEGER NOT NULL,
                    edge_count INTEGER NOT NULL,
                    graph_json TEXT NOT NULL,
                    html TEXT NOT NULL,
                    UNIQUE(graph_key, revision_number),
                    FOREIGN KEY(graph_key) REFERENCES graphs(graph_key) ON DELETE CASCADE
                );

                CREATE INDEX IF NOT EXISTS idx_graph_revisions_graph
                ON graph_revisions(graph_key, revision_number DESC);
                """
            )
            columns = {
                row["name"]
                for row in conn.execute("PRAGMA table_info(graph_revisions)").fetchall()
            }
            if "visibility" not in columns:
                conn.execute(
                    "ALTER TABLE graph_revisions "
                    "ADD COLUMN visibility TEXT NOT NULL DEFAULT 'all'"
                )

    @staticmethod
    def _summary_from_row(row: sqlite3.Row) -> GraphRevisionSummary:
        return GraphRevisionSummary(
            revision_id=row["revision_id"],
            graph_key=row["graph_key"],
            graph_name=row["graph_name"],
            visibility=normalize_graph_visibility(row["visibility"]),
            revision_number=row["revision_number"],
            created_at=row["created_at"],
            node_count=row["node_count"],
            edge_count=row["edge_count"],
            note=row["note"],
        )

    @staticmethod
    def _revision_from_row(row: sqlite3.Row) -> GraphRevision:
        return GraphRevision(
            revision_id=row["revision_id"],
            graph_key=row["graph_key"],
            graph_name=row["graph_name"],
            visibility=normalize_graph_visibility(row["visibility"]),
            revision_number=row["revision_number"],
            created_at=row["created_at"],
            node_count=row["node_count"],
            edge_count=row["edge_count"],
            note=row["note"],
            graph_json=row["graph_json"],
            html=row["html"],
        )

    def save_graph(
        self,
        *,
        graph_key: str,
        graph_name: str,
        kg: KnowledgeGraph,
        html: str,
        visibility: GraphVisibility | str = GraphVisibility.ALL,
        note: str | None = None,
    ) -> GraphRevisionSummary:
        timestamp = _utc_now()
        normalized_visibility = (
            visibility
            if isinstance(visibility, GraphVisibility)
            else normalize_graph_visibility(visibility)
        )
        graph_json = json.dumps(kg.model_dump(mode="json"), indent=2)

        with self.connection as conn:
            existing = conn.execute(
                "SELECT created_at FROM graphs WHERE graph_key = ?",
                (graph_key,),
            ).fetchone()
            created_at = existing["created_at"] if existing else timestamp
            conn.execute(
                """
                INSERT INTO graphs(graph_key, graph_name, created_at, updated_at)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(graph_key) DO UPDATE SET
                    graph_name = excluded.graph_name,
                    updated_at = excluded.updated_at
                """,
                (graph_key, graph_name, created_at, timestamp),
            )

            cursor = conn.execute(
                """
                INSERT INTO graph_revisions(
                    graph_key,
                    graph_name,
                    visibility,
                    revision_number,
                    created_at,
                    note,
                    node_count,
                    edge_count,
                    graph_json,
                    html
                )
                SELECT ?, ?, ?, COALESCE(MAX(revision_number), 0) + 1, ?, ?, ?, ?, ?, ?
                FROM graph_revisions
                WHERE graph_key = ?
                """,
                (
                    graph_key,
                    graph_name,
                    normalized_visibility.value,
                    timestamp,
                    note,
                    len(kg.nodes),
                    len(kg.edges),
                    graph_json,
                    html,
                    graph_key,
                ),
            )

            row = conn.execute(
                "SELECT * FROM graph_revisions WHERE revision_id = ?",
                (cursor.lastrowid,),
            ).fetchone()

        return self._summary_from_row(row)

    def save_engine(
        self,
        *,
        graph_key: str,
        graph_name: str,
        engine: GraphEngine,
        template_path: str | Path | None = None,
        visibility: GraphVisibility | str = GraphVisibility.ALL,
        note: str | None = None,
    ) -> GraphRevisionSummary:
        normalized_visibility = (
            visibility
            if isinstance(visibility, GraphVisibility)
            else normalize_graph_visibility(visibility)
        )
        kg = filter_knowledge_graph(engine.to_knowledge_graph(), normalized_visibility)
        html = to_interactive_html_string(
            engine,
            template_path=template_path,
            visibility=normalized_visibility,
        )
        return self.save_graph(
            graph_key=graph_key,
            graph_name=graph_name,
            kg=kg,
            html=html,
            visibility=normalized_visibility,
            note=note,
        )

    def list_graphs(self) -> list[GraphRevisionSummary]:
        with self.connection as conn:
            rows = conn.execute(
                """
                SELECT r.*
                FROM graph_revisions r
                JOIN (
                    SELECT graph_key, MAX(revision_number) AS max_revision
                    FROM graph_revisions
                    GROUP BY graph_key
                ) latest
                ON latest.graph_key = r.graph_key
                AND latest.max_revision = r.revision_number
                ORDER BY r.created_at DESC, r.graph_name ASC
                """
            ).fetchall()
        return [self._summary_from_row(row) for row in rows]

    def list_revisions(self, graph_key: str) -> list[GraphRevisionSummary]:
        with self.connection as conn:
            rows = conn.execute(
                """
                SELECT *
                FROM graph_revisions
                WHERE graph_key = ?
                ORDER BY revision_number DESC
                """,
                (graph_key,),
            ).fetchall()
        return [self._summary_from_row(row) for row in rows]

    def get_revision(self, revision_id: int) -> GraphRevision | None:
        with self.connection as conn:
            row = conn.execute(
                "SELECT * FROM graph_revisions WHERE revision_id = ?",
                (revision_id,),
            ).fetchone()
        if row is None:
            return None
        return self._revision_from_row(row)

    def get_latest_revision(self, graph_key: str) -> GraphRevision | None:
        with self.connection as conn:
            row = conn.execute(
                """
                SELECT *
                FROM graph_revisions
                WHERE graph_key = ?
                ORDER BY revision_number DESC
                LIMIT 1
                """,
                (graph_key,),
            ).fetchone()
        if row is None:
            return None
        return self._revision_from_row(row)
