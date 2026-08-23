"""Export the knowledge graph to HTML, Plotly HTML, JSON, graph-data.js, and GEXF."""

from __future__ import annotations

import json
import math
import shutil
from collections import defaultdict
from importlib.resources import files as resource_files
from pathlib import Path
from typing import Any

import pyjson5

from .branding import APP_NAME, APP_VERSION
from .config import GraphVisibility, normalize_graph_visibility
from .engine import GraphEngine
from .graph_filters import filtered_engine
from .models import KnowledgeGraph

_NODE_COLORS: dict[str, str] = {
    "Carcinogen": "#e05565",
    "Enzyme": "#4f98a3",
    "Gene": "#3d8b8b",
    "Metabolite": "#e8945a",
    "DNA_Adduct": "#a86fdf",
    "Pathway": "#5591c7",
    "Tissue": "#c2855a",
}

_NODE_SIZES: dict[str, int] = {
    "Carcinogen": 26,
    "Enzyme": 20,
    "Gene": 18,
    "Metabolite": 18,
    "DNA_Adduct": 20,
    "Pathway": 24,
    "Tissue": 17,
}

_NODE_SYMBOLS: dict[str, str] = {
    "Carcinogen": "diamond",
    "Enzyme": "circle",
    "Gene": "circle",
    "Metabolite": "circle",
    "DNA_Adduct": "square",
    "Pathway": "hexagon",
    "Tissue": "triangle-up",
}

_EDGE_COLORS: dict[str, str] = {
    "ACTIVATES": "#e05565",
    "DETOXIFIES": "#6daa45",
    "TRANSPORTS": "#5591c7",
    "FORMS_ADDUCT": "#a86fdf",
    "REPAIRS": "#e8af34",
    "PATHWAY": "#707a8a",
    "EXPRESSED_IN": "#c2855a",
    "INDUCES": "#d4a843",
    "INHIBITS": "#8b4a6b",
    "ENCODES": "#3d8b8b",
    "CUSTOM": "#9ea9bd",
}


def _builtin_viewer_template() -> str:
    """Return a minimal interactive HTML viewer template.

    This fallback keeps HTML export working when no checked-in D3 viewer bundle
    is present in the repo. It intentionally uses only inline CSS/JS so the
    exported file remains self-contained once ``GRAPH_DATA`` is embedded.
    """
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{APP_NAME} Viewer</title>
  <style>
    :root {{
      --bg: #f0f2f5;
      --panel: rgba(255, 255, 255, 0.92);
      --panel-border: rgba(0, 0, 0, 0.12);
      --text: #1a1d23;
      --muted: #5f6673;
      --accent: #01696f;
      --node-stroke: rgba(240, 242, 245, 0.88);
      --edge: rgba(80, 90, 110, 0.35);
      --dim: 0.16;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background:
        radial-gradient(circle at top, rgba(1, 105, 111, 0.08), transparent 42%),
        linear-gradient(180deg, #f0f2f5 0%, #e8eaef 100%);
      color: var(--text);
      min-height: 100vh;
    }}
    .app {{
      max-width: 1600px;
      margin: 0 auto;
      padding: 28px;
    }}
    .header {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      align-items: end;
      margin-bottom: 18px;
    }}
    .title {{
      margin: 0;
      font-size: 1.8rem;
      letter-spacing: 0.02em;
    }}
    .subtitle {{
      margin: 6px 0 0;
      color: var(--muted);
      max-width: 48rem;
      line-height: 1.5;
    }}
    .toolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: center;
      margin-bottom: 18px;
    }}
    .toolbar input {{
      min-width: 260px;
      flex: 1 1 320px;
      border: 1px solid var(--panel-border);
      background: rgba(7, 15, 28, 0.84);
      color: var(--text);
      border-radius: 12px;
      padding: 12px 14px;
      outline: none;
    }}
    .toolbar button {{
      border: 1px solid var(--panel-border);
      background: rgba(11, 25, 44, 0.84);
      color: var(--text);
      border-radius: 12px;
      padding: 12px 14px;
      cursor: pointer;
    }}
    .toolbar button:hover {{
      border-color: rgba(118, 195, 255, 0.48);
    }}
    .summary {{
      color: var(--muted);
      font-size: 0.95rem;
      white-space: nowrap;
    }}
    .layout {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) 320px;
      gap: 18px;
    }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--panel-border);
      border-radius: 18px;
      backdrop-filter: blur(8px);
      box-shadow: 0 22px 60px rgba(0, 0, 0, 0.22);
      overflow: hidden;
    }}
    .graph-wrap {{
      min-height: 760px;
      padding: 10px;
    }}
    #graph {{
      width: 100%;
      height: 760px;
      display: block;
      background:
        radial-gradient(circle at center, rgba(118, 195, 255, 0.08), transparent 46%),
        linear-gradient(180deg, rgba(7, 15, 28, 0.7), rgba(7, 15, 28, 0.92));
      border-radius: 14px;
    }}
    .detail {{
      padding: 18px;
    }}
    .detail h2 {{
      margin: 0 0 8px;
      font-size: 1.15rem;
    }}
    .detail p {{
      margin: 0 0 12px;
      color: var(--muted);
      line-height: 1.5;
    }}
    .detail dl {{
      margin: 0;
      display: grid;
      grid-template-columns: 110px minmax(0, 1fr);
      gap: 8px 10px;
      font-size: 0.92rem;
    }}
    .detail dt {{
      color: var(--muted);
    }}
    .detail dd {{
      margin: 0;
      overflow-wrap: anywhere;
    }}
    .legend {{
      margin-top: 16px;
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }}
    .chip {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      border: 1px solid var(--panel-border);
      border-radius: 999px;
      padding: 6px 10px;
      color: var(--muted);
      font-size: 0.82rem;
    }}
    .chip::before {{
      content: "";
      width: 10px;
      height: 10px;
      border-radius: 50%;
      background: var(--chip-color, var(--accent));
    }}
    .edge {{
      stroke: var(--edge-color, var(--edge));
      stroke-width: 2.2;
      opacity: 0.85;
    }}
    .edge.dimmed,
    .node-group.dimmed,
    .label.dimmed {{
      opacity: var(--dim);
    }}
    .node {{
      stroke: var(--node-stroke);
      stroke-width: 2.5;
      cursor: grab;
    }}
    .node:active {{
      cursor: grabbing;
    }}
    .node.selected {{
      stroke: #ffffff;
      stroke-width: 3;
    }}
    .label {{
      fill: var(--text);
      font-size: 12px;
      font-weight: 600;
      pointer-events: none;
      text-anchor: middle;
    }}
    .footer {{
      margin-top: 14px;
      color: var(--muted);
      font-size: 0.8rem;
    }}
    @media (max-width: 1080px) {{
      .layout {{
        grid-template-columns: 1fr;
      }}
      #graph {{
        height: 560px;
      }}
    }}
  </style>
  <script src="./graph-data.js"></script>
</head>
<body>
  <div class="app">
    <div class="header">
      <div>
        <h1 class="title">{APP_NAME} Graph Viewer</h1>
        <p class="subtitle">Self-contained fallback viewer for exported knowledge graphs. Search, click, and drag nodes to inspect the graph without a checked-in viewer bundle.</p>
      </div>
      <div class="summary" id="summary"></div>
    </div>
    <div class="toolbar">
      <input id="search" type="search" placeholder="Search nodes by id, label, type, detail, or source..." />
      <button id="reset">Reset Focus</button>
    </div>
    <div class="layout">
      <div class="card graph-wrap">
        <svg id="graph" viewBox="0 0 1200 760" role="img" aria-label="Knowledge graph"></svg>
      </div>
      <aside class="card detail">
        <h2 id="detail-title">Graph Overview</h2>
        <p id="detail-text">Click a node to inspect its metadata. Use search to highlight matching nodes and their immediate neighbors.</p>
        <dl id="detail-meta"></dl>
        <div class="legend" id="legend"></div>
      </aside>
    </div>
    <div class="footer">Generated by {APP_NAME} v{APP_VERSION}</div>
  </div>
  <script>
    (() => {{
      const graphData = typeof GRAPH_DATA !== "undefined" ? GRAPH_DATA : {{ nodes: [], edges: [] }};
      const svg = document.getElementById("graph");
      const summary = document.getElementById("summary");
      const search = document.getElementById("search");
      const reset = document.getElementById("reset");
      const detailTitle = document.getElementById("detail-title");
      const detailText = document.getElementById("detail-text");
      const detailMeta = document.getElementById("detail-meta");
      const legend = document.getElementById("legend");
      const width = 1200;
      const height = 760;
      const ns = "http://www.w3.org/2000/svg";
      const colorByType = {{
        Carcinogen: "#e05565",
        Enzyme: "#4f98a3",
        Gene: "#3d8b8b",
        Metabolite: "#e8945a",
        DNA_Adduct: "#a86fdf",
        Pathway: "#5591c7",
        Tissue: "#c2855a",
      }};
      const colorByEdge = {{
        ACTIVATES: "#e05565",
        DETOXIFIES: "#6daa45",
        TRANSPORTS: "#5591c7",
        FORMS_ADDUCT: "#a86fdf",
        REPAIRS: "#e8af34",
        PATHWAY: "#707a8a",
        EXPRESSED_IN: "#c2855a",
        INDUCES: "#d4a843",
        INHIBITS: "#8b4a6b",
        ENCODES: "#3d8b8b",
        CUSTOM: "#9ea9bd",
      }};

      const nodes = (graphData.nodes || []).map((node, index, arr) => {{
        const angle = (Math.PI * 2 * index) / Math.max(arr.length, 1);
        const radius = Math.min(width, height) * 0.34 + (index % 5) * 10;
        return {{
          ...node,
          x: width / 2 + Math.cos(angle) * radius,
          y: height / 2 + Math.sin(angle) * radius,
        }};
      }});
      const edges = (graphData.edges || []).map((edge) => ({{ ...edge }}));
      const nodeById = new Map(nodes.map((node) => [node.id, node]));
      const neighbors = new Map(nodes.map((node) => [node.id, new Set()]));
      edges.forEach((edge) => {{
        neighbors.get(edge.source)?.add(edge.target);
        neighbors.get(edge.target)?.add(edge.source);
      }});

      let searchTerm = "";
      let selectedId = null;
      let draggingId = null;

      function svgPoint(evt) {{
        const point = svg.createSVGPoint();
        point.x = evt.clientX;
        point.y = evt.clientY;
        return point.matrixTransform(svg.getScreenCTM().inverse());
      }}

      function createEl(tag, attrs = {{}}, text = "") {{
        const el = document.createElementNS(ns, tag);
        Object.entries(attrs).forEach(([key, value]) => el.setAttribute(key, String(value)));
        if (text) {{
          el.textContent = text;
        }}
        return el;
      }}

      function clear(el) {{
        while (el.firstChild) {{
          el.removeChild(el.firstChild);
        }}
      }}

      function nodeMatches(node, term) {{
        if (!term) {{
          return true;
        }}
        const haystack = [
          node.id,
          node.label,
          node.type,
          node.detail,
          node.group,
          node.source_db,
          node.match_status,
          node.origin,
        ].filter(Boolean).join(" ").toLowerCase();
        return haystack.includes(term);
      }}

      function highlightedIds() {{
        const term = searchTerm.trim().toLowerCase();
        if (!term) {{
          return null;
        }}
        const matches = new Set(nodes.filter((node) => nodeMatches(node, term)).map((node) => node.id));
        for (const nodeId of Array.from(matches)) {{
          for (const neighbor of neighbors.get(nodeId) || []) {{
            matches.add(neighbor);
          }}
        }}
        return matches;
      }}

      function setDetails(node) {{
        clear(detailMeta);
        if (!node) {{
          detailTitle.textContent = "Graph Overview";
          detailText.textContent = "Click a node to inspect its metadata. Use search to highlight matching nodes and their immediate neighbors.";
          return;
        }}

        detailTitle.textContent = node.label || node.id;
        detailText.textContent = node.detail || "No detail available for this node.";

        const fields = [
          ["ID", node.id],
          ["Type", node.type],
          ["Origin", node.origin],
          ["Match", node.match_status],
          ["Group", node.group],
          ["IARC", node.iarc],
          ["Phase", node.phase],
          ["Role", node.role],
          ["Reactivity", node.reactivity],
          ["Tissue", node.tissue],
          ["Variant", node.variant],
          ["Phenotype", node.phenotype],
          ["Source", node.source_db],
          ["Evidence", node.evidence],
          ["PMID", node.pmid],
          ["Canonical", node.canonical_label || node.canonical_id],
        ].filter(([, value]) => value !== undefined && value !== null && value !== "");

        for (const [label, value] of fields) {{
          const dt = document.createElement("dt");
          dt.textContent = label;
          const dd = document.createElement("dd");
          dd.textContent = String(value);
          detailMeta.appendChild(dt);
          detailMeta.appendChild(dd);
        }}
      }}

      function renderLegend() {{
        clear(legend);
        const usedTypes = Array.from(new Set(nodes.map((node) => node.type))).sort();
        for (const type of usedTypes) {{
          const chip = document.createElement("span");
          chip.className = "chip";
          chip.style.setProperty("--chip-color", colorByType[type] || "#76c3ff");
          chip.textContent = type;
          legend.appendChild(chip);
        }}
      }}

      function render() {{
        clear(svg);

        const defs = createEl("defs");
        const marker = createEl("marker", {{
          id: "arrow",
          markerWidth: 10,
          markerHeight: 10,
          refX: 8,
          refY: 3,
          orient: "auto",
          markerUnits: "strokeWidth",
        }});
        marker.appendChild(createEl("path", {{ d: "M0,0 L0,6 L9,3 z", fill: "#8ea4bb" }}));
        defs.appendChild(marker);
        svg.appendChild(defs);

        const focusIds = highlightedIds();
        summary.textContent = `${{nodes.length}} nodes, ${{edges.length}} edges`;

        edges.forEach((edge) => {{
          const source = nodeById.get(edge.source);
          const target = nodeById.get(edge.target);
          if (!source || !target) {{
            return;
          }}
          const line = createEl("line", {{
            x1: source.x,
            y1: source.y,
            x2: target.x,
            y2: target.y,
            class: "edge",
            "marker-end": "url(#arrow)",
          }});
          line.style.setProperty("--edge-color", colorByEdge[edge.type] || "#8ea4bb");
          if (focusIds && (!focusIds.has(edge.source) || !focusIds.has(edge.target))) {{
            line.classList.add("dimmed");
          }}
          const title = createEl("title", {{}}, `${{edge.type}}: ${{edge.source}} -> ${{edge.target}}`);
          line.appendChild(title);
          svg.appendChild(line);
        }});

        nodes.forEach((node) => {{
          const group = createEl("g", {{ class: "node-group" }});
          if (focusIds && !focusIds.has(node.id)) {{
            group.classList.add("dimmed");
          }}

          const circle = createEl("circle", {{
            cx: node.x,
            cy: node.y,
            r: selectedId === node.id ? 18 : 15,
            class: "node" + (selectedId === node.id ? " selected" : ""),
            fill: colorByType[node.type] || "#76c3ff",
          }});
          circle.addEventListener("mousedown", (evt) => {{
            draggingId = node.id;
            evt.preventDefault();
          }});
          circle.addEventListener("click", () => {{
            selectedId = node.id;
            setDetails(node);
            render();
          }});
          circle.appendChild(createEl("title", {{}}, `${{node.label || node.id}} (${{node.type}})`));

          const label = createEl("text", {{
            x: node.x,
            y: node.y + 30,
            class: "label",
          }}, node.label || node.id);
          if (focusIds && !focusIds.has(node.id)) {{
            label.classList.add("dimmed");
          }}

          group.appendChild(circle);
          group.appendChild(label);
          svg.appendChild(group);
        }});
      }}

      svg.addEventListener("mousemove", (evt) => {{
        if (!draggingId) {{
          return;
        }}
        const node = nodeById.get(draggingId);
        if (!node) {{
          return;
        }}
        const pt = svgPoint(evt);
        node.x = Math.max(24, Math.min(width - 24, pt.x));
        node.y = Math.max(24, Math.min(height - 24, pt.y));
        render();
      }});

      window.addEventListener("mouseup", () => {{
        draggingId = null;
      }});

      search.addEventListener("input", () => {{
        searchTerm = search.value;
        render();
      }});

      reset.addEventListener("click", () => {{
        searchTerm = "";
        selectedId = null;
        search.value = "";
        setDetails(null);
        render();
      }});

      renderLegend();
      setDetails(null);
      render();
    }})();
  </script>
</body>
</html>
"""


def _extract_graph_data_object(raw: str) -> str:
    """Extract the object assigned to ``GRAPH_DATA`` from JS or HTML text."""
    marker = "GRAPH_DATA"
    marker_pos = raw.find(marker)
    if marker_pos == -1:
        start = raw.find("{")
        end = raw.rfind("}")
        if start == -1 or end == -1 or start >= end:
            raise ValueError("Could not locate graph data object")
        return raw[start:end + 1]

    start = raw.find("{", marker_pos)
    if start == -1:
        raise ValueError("Could not locate start of GRAPH_DATA object")

    depth = 0
    in_string = False
    string_char = ""
    escape = False

    for idx in range(start, len(raw)):
        ch = raw[idx]
        if in_string:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == string_char:
                in_string = False
        else:
            if ch in ('"', "'", "`"):
                in_string = True
                string_char = ch
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    return raw[start:idx + 1]

    raise ValueError("Could not find the end of the GRAPH_DATA object")



def parse_graph_data_text(raw: str) -> KnowledgeGraph:
    """Parse JS/HTML text that contains a ``GRAPH_DATA`` object assignment.

    Uses a JSON5 parser to handle unquoted keys, single-line comments,
    and trailing commas that are valid in JavaScript but not in JSON.
    """
    js_obj = _extract_graph_data_object(raw)
    data = pyjson5.loads(js_obj)
    return KnowledgeGraph(**data)


def parse_graph_data_js(path: str | Path) -> KnowledgeGraph:
    """Parse a ``graph-data.js`` file."""
    raw = Path(path).read_text(encoding="utf-8")
    return parse_graph_data_text(raw)


def parse_graph_html(path: str | Path) -> KnowledgeGraph:
    """Parse a standalone HTML graph export with embedded ``GRAPH_DATA``."""
    raw = Path(path).read_text(encoding="utf-8")
    return parse_graph_data_text(raw)


def parse_graph_artifact(path: str | Path) -> KnowledgeGraph:
    """Parse JSON, JS, or HTML graph artifacts."""
    path = Path(path)
    raw = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        return KnowledgeGraph(**json.loads(raw))
    return parse_graph_data_text(raw)


def _clean_for_js(obj: Any) -> Any:
    """Strip None values so the JS side sees clean objects."""
    if isinstance(obj, dict):
        return {k: _clean_for_js(v) for k, v in obj.items() if v is not None}
    if isinstance(obj, list):
        return [_clean_for_js(i) for i in obj]
    return obj


def _export_engine(
    engine: GraphEngine,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> GraphEngine:
    normalized = (
        visibility
        if isinstance(visibility, GraphVisibility)
        else normalize_graph_visibility(visibility)
    )
    if normalized == GraphVisibility.ALL:
        return engine
    return filtered_engine(engine, normalized)


def _require_plotly() -> Any:
    try:
        import plotly.graph_objects as go  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover - exercised only when dependency missing
        raise ImportError(
            "Plotly export requires the `plotly` package. "
            "Install `ExposoGraph[notebook]` or `pip install plotly`."
        ) from exc
    return go


def _plotly_hover_lines(data: dict[str, Any], fields: list[tuple[str, str]]) -> str:
    lines = []
    for label, key in fields:
        value = data.get(key)
        if value not in (None, ""):
            lines.append(f"<b>{label}:</b> {value}")
    return "<br>".join(lines)


def _plotly_node_hover(node: dict[str, Any]) -> str:
    title = node.get("label") or node.get("id") or "Node"
    detail = node.get("detail") or "No detail available."
    meta = _plotly_hover_lines(
        node,
        [
            ("ID", "id"),
            ("Type", "type"),
            ("Origin", "origin"),
            ("Match", "match_status"),
            ("Canonical", "canonical_label"),
            ("Canonical ID", "canonical_id"),
            ("Group", "group"),
            ("IARC", "iarc"),
            ("Phase", "phase"),
            ("Role", "role"),
            ("Reactivity", "reactivity"),
            ("Tissue", "tissue"),
            ("Variant", "variant"),
            ("Phenotype", "phenotype"),
            ("Source", "source_db"),
            ("Evidence", "evidence"),
            ("PMID", "pmid"),
        ],
    )
    return f"<b>{title}</b><br>{detail}" + (f"<br><br>{meta}" if meta else "")


def _plotly_edge_hover(edge: dict[str, Any]) -> str:
    title = f"{edge.get('type', 'EDGE')}: {edge.get('source')} -> {edge.get('target')}"
    meta = _plotly_hover_lines(
        edge,
        [
            ("Label", "label"),
            ("Origin", "origin"),
            ("Match", "match_status"),
            ("Canonical", "canonical_predicate"),
            ("Custom", "custom_predicate"),
            ("Source", "source_db"),
            ("Evidence", "evidence"),
            ("PMID", "pmid"),
        ],
    )
    return title + (f"<br><br>{meta}" if meta else "")


def _plotly_positions(engine: GraphEngine) -> dict[str, tuple[float, float]]:
    import networkx as nx

    graph = engine.G
    if graph.number_of_nodes() == 0:
        return {}
    if graph.number_of_nodes() == 1:
        node_id = next(iter(graph.nodes()))
        return {str(node_id): (0.0, 0.0)}

    layout_graph = nx.Graph()
    layout_graph.add_nodes_from(graph.nodes())
    layout_graph.add_edges_from((u, v) for u, v in graph.edges())
    spring_k = 1.6 / math.sqrt(max(layout_graph.number_of_nodes(), 1))
    pos = nx.spring_layout(layout_graph, seed=42, k=spring_k, iterations=200)
    return {
        str(node_id): (float(coords[0]), float(coords[1]))
        for node_id, coords in pos.items()
    }


def _graph_data_script(
    engine: GraphEngine,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> str:
    data = _clean_for_js(_export_engine(engine, visibility).to_dict())
    return f"const GRAPH_DATA = {json.dumps(data, indent=2)};"


def _graph_data_script_from_dict(graph: dict[str, list[Any]]) -> str:
    """Render an arbitrary ``{"nodes": [...], "edges": [...]}`` dict as a
    ``GRAPH_DATA`` script body -- the dict-based analog of
    :func:`_graph_data_script`, which takes a whole :class:`GraphEngine`
    instead. Used for subgraphs produced by ``GraphEngine`` filtering helpers
    (e.g. :meth:`GraphEngine.map_viewer_subgraph`) that were never assembled
    into a standalone engine.
    """
    data = _clean_for_js(graph)
    return f"const GRAPH_DATA = {json.dumps(data, indent=2)};"


def _inline_graph_data_script(template: str, data_script: str) -> str:
    script_tag = f"<script>\n{data_script.rstrip()}\n</script>"
    external_tag = '<script src="./graph-data.js"></script>'

    if external_tag in template:
        return template.replace(external_tag, script_tag, 1)

    head_close = "</head>"
    if head_close in template:
        return template.replace(head_close, f"{script_tag}\n{head_close}", 1)

    return f"{script_tag}\n{template}"


def _default_template_candidates() -> list[Path]:
    repo_root = Path(__file__).resolve().parent.parent
    return [
        repo_root / "references" / "knowledge-graph" / "index.html",
    ]


def _package_template_text() -> str | None:
    """Load the checked-in D3 HTML template from package resources.

    This makes exports work after installation without relying on a repo-local
    `references/knowledge-graph/` folder.
    """
    try:
        return (resource_files(__package__) / "temp_kg" / "index.html").read_text(encoding="utf-8")
    except Exception:
        return None


def _load_viewer_template(template_path: str | Path | None = None) -> str:
    candidates: list[Path] = []
    if template_path is not None:
        template_path = Path(template_path)
        if template_path.is_dir():
            candidates.append(template_path / "index.html")
        else:
            candidates.append(template_path)
    candidates.extend(_default_template_candidates())

    for candidate in candidates:
        if candidate.exists():
            return candidate.read_text(encoding="utf-8")

    pkg_template = _package_template_text()
    if pkg_template is not None:
        return pkg_template

    return _builtin_viewer_template()


def to_graph_data_js(
    engine: GraphEngine,
    path: str | Path,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> Path:
    """Write a ``graph-data.js`` file consumable by the D3 viewer."""
    path = Path(path)
    header = (
        "// ═══════════════════════════════════════════════════════════════\n"
        f"// Auto-generated by {APP_NAME}\n"
        "// ═══════════════════════════════════════════════════════════════\n\n"
    )

    js = header + _graph_data_script(engine, visibility=visibility) + "\n"
    path.write_text(js, encoding="utf-8")
    return path


def to_plotly_figure(
    engine: GraphEngine,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    title: str | None = None,
) -> Any:
    """Build an interactive Plotly figure for the graph."""
    go = _require_plotly()
    export_engine = _export_engine(engine, visibility)
    graph_data = _clean_for_js(export_engine.to_dict())
    nodes = graph_data.get("nodes", [])
    edges = graph_data.get("edges", [])
    positions = _plotly_positions(export_engine)

    edge_groups: dict[str, dict[str, list[Any]]] = defaultdict(
        lambda: {"x": [], "y": [], "mid_x": [], "mid_y": [], "hover": []}
    )
    for edge in edges:
        source = positions.get(edge.get("source", ""))
        target = positions.get(edge.get("target", ""))
        if source is None or target is None:
            continue
        edge_type = edge.get("type", "EDGE")
        bundle = edge_groups[edge_type]
        bundle["x"].extend([source[0], target[0], None])
        bundle["y"].extend([source[1], target[1], None])
        bundle["mid_x"].append((source[0] + target[0]) / 2)
        bundle["mid_y"].append((source[1] + target[1]) / 2)
        bundle["hover"].append(_plotly_edge_hover(edge))

    node_groups: dict[str, dict[str, list[Any]]] = defaultdict(
        lambda: {"x": [], "y": [], "text": [], "hover": [], "ids": []}
    )
    for node in nodes:
        node_id = node.get("id")
        if not node_id or node_id not in positions:
            continue
        node_type = node.get("type", "Node")
        bundle = node_groups[node_type]
        bundle["x"].append(positions[node_id][0])
        bundle["y"].append(positions[node_id][1])
        bundle["text"].append(node.get("label", node_id))
        bundle["hover"].append(_plotly_node_hover(node))
        bundle["ids"].append(node_id)

    figure = go.Figure()

    for edge_type in sorted(edge_groups):
        color = _EDGE_COLORS.get(edge_type, "#8ea4bb")
        bundle = edge_groups[edge_type]
        figure.add_trace(
            go.Scatter(
                x=bundle["x"],
                y=bundle["y"],
                mode="lines",
                line={"color": color, "width": 2},
                opacity=0.38 if edge_type == "PATHWAY" else 0.62,
                hoverinfo="skip",
                name=f"{edge_type} edge",
                legendgroup=f"edge-{edge_type}",
            )
        )
        if bundle["mid_x"]:
            figure.add_trace(
                go.Scatter(
                    x=bundle["mid_x"],
                    y=bundle["mid_y"],
                    mode="markers",
                    marker={"size": 10, "color": color, "opacity": 0.01},
                    hovertemplate="%{text}<extra></extra>",
                    text=bundle["hover"],
                    showlegend=False,
                    legendgroup=f"edge-{edge_type}",
                )
            )

    for node_type in sorted(node_groups):
        color = _NODE_COLORS.get(node_type, "#76c3ff")
        symbol = _NODE_SYMBOLS.get(node_type, "circle")
        size = _NODE_SIZES.get(node_type, 18)
        bundle = node_groups[node_type]
        figure.add_trace(
            go.Scatter(
                x=bundle["x"],
                y=bundle["y"],
                mode="markers+text",
                text=bundle["text"],
                textposition="bottom center",
                textfont={"size": 11, "color": "#dfe8f4"},
                hovertemplate="%{hovertext}<extra></extra>",
                hovertext=bundle["hover"],
                marker={
                    "size": size,
                    "color": color,
                    "symbol": symbol,
                    "line": {"color": "#f0f2f5", "width": 1.5},
                    "opacity": 0.92,
                },
                name=node_type,
                legendgroup=f"node-{node_type}",
                customdata=bundle["ids"],
            )
        )

    figure.update_layout(
        title=title or f"{APP_NAME} Knowledge Graph",
        paper_bgcolor="#f0f2f5",
        plot_bgcolor="#f0f2f5",
        font={"family": "Inter, system-ui, sans-serif", "color": "#1a1d23"},
        hoverlabel={"bgcolor": "#ffffff", "font": {"color": "#1a1d23"}},
        margin={"l": 24, "r": 24, "t": 68, "b": 24},
        showlegend=True,
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "left",
            "x": 0,
            "bgcolor": "rgba(8, 17, 31, 0.7)",
            "bordercolor": "rgba(142, 176, 205, 0.22)",
            "borderwidth": 1,
            "font": {"size": 11},
        },
        xaxis={"visible": False},
        yaxis={"visible": False, "scaleanchor": "x", "scaleratio": 1},
        dragmode="pan",
        annotations=[
            {
                "text": f"{len(nodes)} nodes · {len(edges)} edges",
                "xref": "paper",
                "yref": "paper",
                "x": 1,
                "y": 1.06,
                "showarrow": False,
                "font": {"size": 12, "color": "#5f6673"},
                "xanchor": "right",
            }
        ],
    )
    figure.update_xaxes(fixedrange=False)
    figure.update_yaxes(fixedrange=False)
    return figure


def to_plotly_html_string(
    engine: GraphEngine,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    title: str | None = None,
    include_plotlyjs: bool | str = True,
) -> str:
    """Render a standalone Plotly HTML document."""
    figure = to_plotly_figure(engine, visibility=visibility, title=title)
    html: str = figure.to_html(
        full_html=True,
        include_plotlyjs=include_plotlyjs,
        config={
            "displaylogo": False,
            "responsive": True,
            "scrollZoom": True,
        },
    )
    return html


def to_plotly_html(
    engine: GraphEngine,
    path: str | Path,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    title: str | None = None,
    include_plotlyjs: bool | str = True,
) -> Path:
    """Write a standalone Plotly HTML file."""
    path = Path(path)
    html = to_plotly_html_string(
        engine,
        visibility=visibility,
        title=title,
        include_plotlyjs=include_plotlyjs,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path


def to_interactive_html_string(
    engine: GraphEngine,
    *,
    template_path: str | Path | None = None,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> str:
    """Render a self-contained interactive HTML document."""
    template = _load_viewer_template(template_path)
    data_script = _graph_data_script(engine, visibility=visibility)
    return _inline_graph_data_script(template, data_script)


def subgraph_to_html_string(
    graph: dict[str, list[Any]],
    *,
    template_path: str | Path | None = None,
) -> str:
    """Render an arbitrary node/edge dict as a self-contained interactive HTML document.

    Dict-based analog of :func:`to_interactive_html_string`, for callers that
    already have a subgraph dict (e.g. from
    :meth:`GraphEngine.map_viewer_subgraph`/:meth:`GraphEngine.dim_by_tissue_threshold`)
    rather than a whole ``GraphEngine`` to select ``visibility=`` on. This is
    the mechanism ``ui_map_viewer.py`` uses to re-render on every filter
    interaction: recompute the filtered dict with the engine's filtering
    helpers, then hand it straight to this function instead of constructing a
    throwaway ``GraphEngine`` just to export it.
    """
    template = _load_viewer_template(template_path)
    data_script = _graph_data_script_from_dict(graph)
    return _inline_graph_data_script(template, data_script)


def bundle_to_html_string(
    template_path: str | Path,
    graph_data_path: str | Path,
) -> str:
    """Render a viewer bundle directory as a self-contained HTML document."""
    template = Path(template_path).read_text(encoding="utf-8")
    data_script = Path(graph_data_path).read_text(encoding="utf-8")
    return _inline_graph_data_script(template, data_script)


def to_interactive_html(
    engine: GraphEngine,
    path: str | Path,
    *,
    template_path: str | Path | None = None,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> Path:
    """Write a standalone interactive HTML file with embedded graph data."""
    path = Path(path)
    html = to_interactive_html_string(
        engine,
        template_path=template_path,
        visibility=visibility,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path


def ensure_viewer_bundle(
    output_dir: str | Path,
    *,
    template_dir: str | Path | None = None,
) -> Path:
    """Ensure *output_dir* is a usable D3 viewer bundle directory.

    When *template_dir* contains an ``index.html`` file, it is copied into the
    output directory if missing so the generated ``graph-data.js`` has a viewer
    to pair with.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if template_dir is not None:
        template_dir = Path(template_dir)
        template_index = template_dir / "index.html"
        target_index = output_dir / "index.html"
        if template_index.exists() and template_index != target_index and not target_index.exists():
            shutil.copy2(template_index, target_index)

    target_index = output_dir / "index.html"
    if not target_index.exists():
        # If a repo-local template isn't present, use the package-embedded
        # one so viewer bundles still work after installation.
        target_index.write_text(_load_viewer_template(template_dir), encoding="utf-8")

    return output_dir


def export_viewer_bundle(
    engine: GraphEngine,
    output_dir: str | Path,
    *,
    template_dir: str | Path | None = None,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> Path:
    """Write a complete viewer bundle directory.

    The resulting folder contains ``graph-data.js`` and, when a template is
    available, ``index.html``.
    """
    output_dir = ensure_viewer_bundle(output_dir, template_dir=template_dir)
    to_graph_data_js(engine, output_dir / "graph-data.js", visibility=visibility)
    return output_dir


def to_json(
    engine: GraphEngine,
    path: str | Path,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> Path:
    """Write a plain JSON export."""
    path = Path(path)
    data = _clean_for_js(_export_engine(engine, visibility).to_dict())
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return path


def to_gexf(
    engine: GraphEngine,
    path: str | Path,
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
) -> Path:
    """Write GEXF (Gephi) format.

    Non-scalar node/edge attributes (lists, dicts) are serialized to
    JSON strings since GEXF only supports primitive attribute types.
    """
    import networkx as nx

    path = Path(path)
    # Deep-copy graph and flatten complex attributes for GEXF compatibility
    G = _export_engine(engine, visibility).G.copy()
    _SCALAR = (str, int, float, bool)
    for _node_id, data in G.nodes(data=True):
        for key, val in list(data.items()):
            if val is not None and not isinstance(val, _SCALAR):
                data[key] = json.dumps(val, default=str)
    for _u, _v, _k, data in G.edges(keys=True, data=True):
        for key, val in list(data.items()):
            if val is not None and not isinstance(val, _SCALAR):
                data[key] = json.dumps(val, default=str)
    nx.write_gexf(G, str(path))
    return path
