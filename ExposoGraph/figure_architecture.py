"""Static architecture-figure helpers for ExposoGraph.

This module builds a publication-style infographic summarizing an
ExposoGraph graph, its curated enzyme panel, and its supported edge
families.  It is intentionally separate from the interactive viewers:
the output is a reproducible Matplotlib figure suitable for PDF/SVG/PNG
export from notebooks and scripts.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from .branding import APP_NAME
from .engine import GraphEngine
from .models import EdgeType, KnowledgeGraph, NodeType
from .reference_data import build_full_panel

_NODE_TYPE_COLORS: dict[str, str] = {
    "Carcinogens": "#ef5345",
    "Enzymes": "#3795de",
    "Metabolites": "#f39a1e",
    "DNA Adducts": "#a25bcc",
    "KEGG Pathways": "#34c96d",
}

_NODE_TYPE_DISPLAY_LABELS: dict[str, str] = {
    "Carcinogens": "CARCINOGEN CONTEXTS",
    "Enzymes": "ENZYMES",
    "Metabolites": "METABOLITES",
    "DNA Adducts": "DNA ADDUCTS",
    "KEGG Pathways": "KEGG PATHWAYS",
}

_EDGE_TYPE_COLORS: dict[str, str] = {
    EdgeType.ACTIVATES.value: "#ef5345",
    EdgeType.DETOXIFIES.value: "#3795de",
    EdgeType.TRANSPORTS.value: "#18a286",
    EdgeType.FORMS_ADDUCT.value: "#a25bcc",
    EdgeType.REPAIRS.value: "#34c96d",
    EdgeType.PATHWAY.value: "#28a25f",
}

_EDGE_TYPE_LABELS: dict[str, tuple[str, str]] = {
    EdgeType.ACTIVATES.value: ("ACTIVATES", "Enzyme -> Metabolite"),
    EdgeType.DETOXIFIES.value: ("DETOXIFIES", "Enzyme -> Metabolite"),
    EdgeType.TRANSPORTS.value: ("TRANSPORTS", "Phase III transporter -> Metabolite"),
    EdgeType.FORMS_ADDUCT.value: ("FORMS_ADDUCT", "Metabolite -> DNA Adduct"),
    EdgeType.REPAIRS.value: ("REPAIRS", "Repair Enzyme -> Adduct"),
    EdgeType.PATHWAY.value: ("PATHWAY", "Entity layer -> \nKEGG pathway"),
}

_EDGE_TYPE_STYLES: dict[str, str] = {
    EdgeType.ACTIVATES.value: "-",
    EdgeType.DETOXIFIES.value: "--",
    EdgeType.TRANSPORTS.value: ":",
    EdgeType.FORMS_ADDUCT.value: "-",
    EdgeType.REPAIRS.value: "-",
    EdgeType.PATHWAY.value: "-.",
}

_DEFAULT_CARCINOGEN_CLASSES: tuple[tuple[str, tuple[str, ...], str, tuple[str, ...]], ...] = (
    ("PAH", ("Benzo[a]pyrene", "DMBA"), "#d64f45", ("PAH",)),
    ("HCA", ("PhIP", "MeIQx"), "#e6861f", ("HCA",)),
    ("Aromatic Amines", ("4-Aminobiphenyl", "Benzidine"), "#dc6a00", ("Aromatic_Amine", "Aromatic Amine")),
    ("Nitrosamines", ("NNK", "NDMA"), "#c24632", ("Nitrosamine",)),
    ("Mycotoxins", ("Aflatoxin B1",), "#9454c6", ("Mycotoxin",)),
    ("Estrogens", ("17beta-Estradiol",), "#dc2b98", ("Estrogen", "Estrogens")),
    ("Androgens", ("Testosterone", "5a-DHT"), "#3281c0", ("Androgen", "Androgens")),
    ("Solvents", ("Benzene", "Vinyl Chloride"), "#1fa287", ("Solvent",)),
    ("Alkylating Agents", ("Ethylene Oxide",), "#7d8a92", ("Alkylating",)),
)


@dataclass(frozen=True)
class ArchitectureListItem:
    title: str
    examples: tuple[str, ...]
    color: str
    count: int | None = None


@dataclass(frozen=True)
class ArchitectureEdgeLegendItem:
    edge_type: str
    label: str
    description: str
    color: str
    linestyle: str
    count: int


@dataclass(frozen=True)
class ArchitectureFigureData:
    title: str
    subtitle: str
    layer_counts: dict[str, int]
    carcinogen_classes: tuple[ArchitectureListItem, ...]
    enzyme_categories: tuple[ArchitectureListItem, ...]
    edge_legend: tuple[ArchitectureEdgeLegendItem, ...]
    summary_lines: tuple[str, ...]


def _require_matplotlib() -> tuple[Any, Any, Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
    except ImportError as exc:  # pragma: no cover - dependency failure path
        raise ImportError(
            "Architecture figure export requires `matplotlib`. "
            "Install `ExposoGraph[notebook]` or `pip install matplotlib`."
        ) from exc
    return plt, FancyBboxPatch, Circle, FancyArrowPatch


def _coerce_engine(graph: GraphEngine | KnowledgeGraph | None) -> GraphEngine:
    if isinstance(graph, GraphEngine):
        return graph
    engine = GraphEngine()
    if graph is None:
        return engine
    for node in graph.nodes:
        engine.add_node(node)
    for edge in graph.edges:
        engine.add_edge(edge)
    return engine


def _example_string(examples: Sequence[str]) -> str:
    return ", ".join(examples)


def _panel_nodes(reference_panel: KnowledgeGraph | None) -> list[Any]:
    panel = reference_panel or build_full_panel()
    return list(panel.nodes)


def _enzyme_category_items(
    reference_panel: KnowledgeGraph | None,
    *,
    count_overrides: Mapping[str, int] | None = None,
) -> tuple[ArchitectureListItem, ...]:
    panel_nodes = _panel_nodes(reference_panel)
    overrides = dict(count_overrides or {})

    phase_i = [node.label for node in panel_nodes if node.phase == "I"]
    phase_ii = [node.label for node in panel_nodes if node.phase == "II"]
    phase_iii = [node.label for node in panel_nodes if node.phase == "III"]
    repair = [node.label for node in panel_nodes if (node.group or "").startswith("DNA Repair")]

    items = (
        ("Phase I -- Activation", phase_i, "#2f80ed"),
        ("Phase II -- Conjugation", phase_ii, "#5d7fa1"),
        ("Phase III -- Transport", phase_iii, "#24a49c"),
        ("DNA Repair", repair, "#8b5bb8"),
    )
    return tuple(
        ArchitectureListItem(
            title=title,
            examples=tuple(labels),
            color=color,
            count=int(overrides.get(title, len(labels))),
        )
        for title, labels, color in items
    )


def _carcinogen_class_items(
    graph: GraphEngine,
    *,
    class_catalog: Sequence[tuple[str, tuple[str, ...], str, tuple[str, ...]]] | None = None,
) -> tuple[ArchitectureListItem, ...]:
    catalog = tuple(class_catalog or _DEFAULT_CARCINOGEN_CLASSES)
    category_counts: Counter[str] = Counter()

    for _, data in graph.G.nodes(data=True):
        if data.get("type") != NodeType.CARCINOGEN.value:
            continue
        group = str(data.get("group") or "").strip()
        if group:
            category_counts[group] += 1

    items: list[ArchitectureListItem] = []
    for title, examples, color, aliases in catalog:
        count = 0
        for alias in aliases:
            count += category_counts.get(alias, 0)
        items.append(
            ArchitectureListItem(
                title=title,
                examples=examples,
                color=color,
                count=count if count > 0 else None,
            )
        )
    return tuple(items)


def _layer_counts(
    graph: GraphEngine,
    *,
    count_overrides: Mapping[str, int] | None = None,
) -> dict[str, int]:
    node_counts = Counter(str(data.get("type", "")) for _, data in graph.G.nodes(data=True))
    counts = {
        "Carcinogens": node_counts.get(NodeType.CARCINOGEN.value, 0),
        "Enzymes": node_counts.get(NodeType.ENZYME.value, 0) + node_counts.get(NodeType.GENE.value, 0),
        "Metabolites": node_counts.get(NodeType.METABOLITE.value, 0),
        "DNA Adducts": node_counts.get(NodeType.DNA_ADDUCT.value, 0),
        "KEGG Pathways": node_counts.get(NodeType.PATHWAY.value, 0),
    }
    for key, value in dict(count_overrides or {}).items():
        counts[str(key)] = int(value)
    return counts


def _edge_legend_items(
    graph: GraphEngine,
    *,
    count_overrides: Mapping[str, int] | None = None,
) -> tuple[ArchitectureEdgeLegendItem, ...]:
    edge_counts = Counter(str(data.get("type", "")) for _, _, data in graph.G.edges(data=True))
    overrides = dict(count_overrides or {})
    ordered = (
        EdgeType.ACTIVATES.value,
        EdgeType.DETOXIFIES.value,
        EdgeType.TRANSPORTS.value,
        EdgeType.FORMS_ADDUCT.value,
        EdgeType.REPAIRS.value,
        EdgeType.PATHWAY.value,
    )
    items: list[ArchitectureEdgeLegendItem] = []
    for edge_type in ordered:
        label, description = _EDGE_TYPE_LABELS[edge_type]
        items.append(
            ArchitectureEdgeLegendItem(
                edge_type=edge_type,
                label=label,
                description=description,
                color=_EDGE_TYPE_COLORS[edge_type],
                linestyle=_EDGE_TYPE_STYLES[edge_type],
                count=int(overrides.get(edge_type, edge_counts.get(edge_type, 0))),
            )
        )
    return tuple(items)


def _summary_lines(
    graph: GraphEngine,
    carcinogen_classes: Sequence[ArchitectureListItem],
    layer_counts: Mapping[str, int],
    edge_legend: Sequence[ArchitectureEdgeLegendItem],
    *,
    summary_overrides: Mapping[str, Any] | None = None,
) -> tuple[str, ...]:
    overrides = dict(summary_overrides or {})
    node_count = int(overrides.get("node_count", graph.node_count))
    edge_count = int(overrides.get("edge_count", graph.edge_count))
    node_type_count = int(
        overrides.get(
            "node_type_count",
            len({str(data.get("type", "")) for _, data in graph.G.nodes(data=True)}),
        )
    )
    edge_type_count = int(
        overrides.get(
            "edge_type_count",
            len({str(data.get("type", "")) for _, _, data in graph.G.edges(data=True)}),
        )
    )
    carcinogen_class_count = int(overrides.get("carcinogen_class_count", len(carcinogen_classes)))
    core_layer_count = int(overrides.get("core_layer_count", len(layer_counts)))
    note = str(
        overrides.get(
            "summary_note",
            "Reproducible Matplotlib export from ExposoGraph",
        )
    )
    edge_total = int(overrides.get("edge_total", sum(item.count for item in edge_legend)))
    return (
        f"{node_count} nodes  ·  {edge_count} edges",
        f"{node_type_count} node types  ·  {edge_type_count} edge types",
        f"{carcinogen_class_count} carcinogen classes  ·  {core_layer_count} architecture layers",
        note if note else f"{edge_total} core-relation edges represented",
    )


def paper_architecture_overrides() -> dict[str, Any]:
    """Return override values for the legacy 133/166 paper architecture snapshot.

    Updated after wiring the procarcinogen-to-intermediate chains for DMBA,
    MeIQx, Benzidine, NDMA, Vinyl Chloride, and 17beta-Estradiol, plus the
    auxiliary Phase I coverage nodes CYP2B6/CYP2C9/CYP2C19/CYP2D6/CYP2F1.
    """
    return {
        "title": "CarcinoGenomic Knowledge Graph -- Platform Architecture",
        "subtitle": "Five-layer data model integrating carcinogen metabolism with pharmacogenomic variation",
        "layer_count_overrides": {
            "Carcinogens": 20,
            "Enzymes": 49,
            "Metabolites": 41,
            "DNA Adducts": 16,
            "KEGG Pathways": 7,
        },
        "enzyme_category_count_overrides": {
            "Phase I -- Activation": 21,
            "Phase II -- Conjugation": 19,
            "Phase III -- Transport": 4,
            "DNA Repair": 5,
        },
        "edge_count_overrides": {
            EdgeType.ACTIVATES.value: 57,
            EdgeType.DETOXIFIES.value: 31,
            EdgeType.TRANSPORTS.value: 9,
            EdgeType.FORMS_ADDUCT.value: 20,
            EdgeType.REPAIRS.value: 19,
            EdgeType.PATHWAY.value: 30,
        },
        "summary_overrides": {
            "node_count": 133,
            "edge_count": 166,
            "node_type_count": 5,
            "edge_type_count": 6,
            "carcinogen_class_count": 10,
            "core_layer_count": 5,
            "edge_total": 166,
            "summary_note": "Interactive D3.js force-directed layout",
        },
    }


def build_architecture_figure_data(
    graph: GraphEngine | KnowledgeGraph | None = None,
    *,
    title: str | None = None,
    subtitle: str = "Five-layer data model integrating carcinogen metabolism with curated reference metadata",
    reference_panel: KnowledgeGraph | None = None,
    layer_count_overrides: Mapping[str, int] | None = None,
    enzyme_category_count_overrides: Mapping[str, int] | None = None,
    edge_count_overrides: Mapping[str, int] | None = None,
    summary_overrides: Mapping[str, Any] | None = None,
    carcinogen_class_catalog: Sequence[tuple[str, tuple[str, ...], str, tuple[str, ...]]] | None = None,
) -> ArchitectureFigureData:
    """Build structured data for the architecture infographic."""
    engine = _coerce_engine(graph)
    carcinogen_classes = _carcinogen_class_items(
        engine,
        class_catalog=carcinogen_class_catalog,
    )
    layer_counts = _layer_counts(engine, count_overrides=layer_count_overrides)
    enzyme_categories = _enzyme_category_items(
        reference_panel,
        count_overrides=enzyme_category_count_overrides,
    )
    edge_legend = _edge_legend_items(engine, count_overrides=edge_count_overrides)
    summary_lines = _summary_lines(
        engine,
        carcinogen_classes,
        layer_counts,
        edge_legend,
        summary_overrides=summary_overrides,
    )
    return ArchitectureFigureData(
        title=title or f"{APP_NAME} -- Platform Architecture",
        subtitle=subtitle,
        layer_counts=layer_counts,
        carcinogen_classes=carcinogen_classes,
        enzyme_categories=enzyme_categories,
        edge_legend=edge_legend,
        summary_lines=summary_lines,
    )


def render_architecture_figure(
    data: ArchitectureFigureData,
    *,
    figsize: tuple[float, float] = (13.84, 10.0),
) -> Any:
    """Render the architecture infographic and return the Matplotlib figure."""
    plt, FancyBboxPatch, Circle, FancyArrowPatch = _require_matplotlib()

    fig, ax = plt.subplots(figsize=figsize, dpi=180)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    fig.patch.set_facecolor("#fbfcfe")
    ax.set_facecolor("#fbfcfe")

    def add_box(
        x: float,
        y: float,
        w: float,
        h: float,
        *,
        face: str,
        edge: str = "#d5deea",
        radius: float = 1.2,
        linewidth: float = 1.2,
        alpha: float = 1.0,
    ) -> Any:
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle=f"round,pad=0.35,rounding_size={radius}",
            linewidth=linewidth,
            edgecolor=edge,
            facecolor=face,
            alpha=alpha,
        )
        ax.add_patch(patch)
        return patch

    def add_text(
        x: float,
        y: float,
        text: str,
        *,
        size: float = 10.0,
        color: str = "#33475b",
        weight: str = "normal",
        ha: str = "left",
        va: str = "center",
        style: str = "normal",
    ) -> None:
        ax.text(
            x,
            y,
            text,
            fontsize=size,
            color=color,
            fontweight=weight,
            ha=ha,
            va=va,
            style=style,
            family="DejaVu Sans",
        )

    def add_arrow(
        start: tuple[float, float],
        end: tuple[float, float],
        *,
        color: str,
        label: str | None = None,
        label_xy: tuple[float, float] | None = None,
        linewidth: float = 1.8,
        linestyle: str = "-",
        mutation_scale: float = 11.0,
        connectionstyle: str = "arc3,rad=0.0",
        label_size: float = 7.0,
        label_ha: str = "center",
    ) -> None:
        patch = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=linewidth,
            linestyle=linestyle,
            color=color,
            connectionstyle=connectionstyle,
        )
        ax.add_patch(patch)
        if label and label_xy:
            add_text(
                label_xy[0],
                label_xy[1],
                label,
                size=label_size,
                color=color,
                style="italic",
                ha=label_ha,
            )

    add_text(50, 96, data.title, size=16, color="#32475b", weight="bold", ha="center")
    add_text(50, 92.5, data.subtitle, size=9.8, color="#627181", ha="center")

    add_box(2, 18, 31, 71, face="#ffffff", edge="#d9e1ec", radius=1.2)
    add_text(17.5, 86.6, f"{len(data.carcinogen_classes)} Carcinogen Classes", size=11.5, weight="bold", ha="center")
    start_y = 81.0
    step = 7.3
    for index, item in enumerate(data.carcinogen_classes):
        y = start_y - index * step
        ax.add_patch(Circle((6, y), 0.45, facecolor=item.color, edgecolor=item.color))
        title = item.title if item.count is None else f"{item.title}  (n={item.count})"
        add_text(8.0, y + 0.6, title, size=9.0, weight="bold")
        add_text(8.0, y - 1.8, _example_string(item.examples), size=7.0, color="#7b8795")

    center_layers = (
        ("Carcinogens", 76.0, 24.0, "#fde1df"),
        ("Enzymes", 61.8, 24.0, "#dceaf9"),
        ("Metabolites", 47.1, 24.0, "#fff0da"),
        ("DNA Adducts", 32.4, 24.0, "#efe3f7"),
        ("KEGG Pathways", 18.2, 24.0, "#ddf7e7"),
    )
    layer_bounds: dict[str, tuple[float, float, float, float]] = {}
    for label, y, height, face in center_layers:
        color = _NODE_TYPE_COLORS[label]
        box_x = 38.0
        box_height = (height / 2.0)-0.5
        box_width = 24.0
        add_box(box_x, y, box_width, box_height, face=face, edge=color, linewidth=1.8, radius=1.1)
        layer_bounds[label] = (box_x, y, box_width, box_height)
        display_label = _NODE_TYPE_DISPLAY_LABELS.get(label, label.upper())
        title_size = 12.3 if label == "Carcinogens" else 13.0
        add_text(50, y + box_height - 3.0, display_label, size=title_size, color=color, weight="bold", ha="center")
        add_text(
            50,
            y + 4.2,
            f"n = {data.layer_counts.get(label, 0)}",
            size=9.2,
            color=color,
            ha="center",
        )
        if label == "Carcinogens":
            add_text(
                50,
                y + 1.8,
                "seed compounds and edge-level context",
                size=6.4,
                color="#9b6c66",
                ha="center",
            )

    _, carcinogen_y, _, carcinogen_h = layer_bounds["Carcinogens"]
    _, enzyme_y, _, enzyme_h = layer_bounds["Enzymes"]
    _, metabolite_y, _, metabolite_h = layer_bounds["Metabolites"]
    _, adduct_y, _, adduct_h = layer_bounds["DNA Adducts"]
    _, pathway_y, _, _ = layer_bounds["KEGG Pathways"]

    # Enzyme-driven reactions converge on the metabolite layer.
    add_arrow(
        (42.1, enzyme_y + enzyme_h - 9.5),
        (44.2, metabolite_y + 5.0),
        color=_EDGE_TYPE_COLORS[EdgeType.ACTIVATES.value],
        label="ACTIVATES",
        label_xy=(42.2, 65.2),
        connectionstyle="arc3,rad=-0.12",
    )
    add_arrow(
        (57.9, enzyme_y + enzyme_h - 9.5),
        (55.8, metabolite_y + 5.0),
        color=_EDGE_TYPE_COLORS[EdgeType.DETOXIFIES.value],
        label="DETOXIFIES",
        label_xy=(57.8, 65.2),
        linestyle=_EDGE_TYPE_STYLES[EdgeType.DETOXIFIES.value],
        connectionstyle="arc3,rad=0.12",
    )
    add_arrow(
        (60.6, enzyme_y + enzyme_h - 4.8),
        (60.6, metabolite_y + 5.3),
        color=_EDGE_TYPE_COLORS[EdgeType.TRANSPORTS.value],
        label="TRANSPORTS",
        label_xy=(60.6, 60.2),
        linestyle=_EDGE_TYPE_STYLES[EdgeType.TRANSPORTS.value],
        mutation_scale=10.0,
        label_size=6.8,
    )

    add_arrow(
        (55.8, metabolite_y + metabolite_h - 9.3),
        (55.8, adduct_y + 5.0),
        color=_EDGE_TYPE_COLORS[EdgeType.FORMS_ADDUCT.value],
        label="FORMS ADDUCT",
        label_xy=(58.8, 45.6),
        mutation_scale=10.5,
        label_size=6.9,
    )
    add_arrow(
        (36.1, enzyme_y + enzyme_h - 4.4),
        (36.1, adduct_y + 4.8),
        color=_EDGE_TYPE_COLORS[EdgeType.REPAIRS.value],
        label="REPAIRS",
        label_xy=(40.4, 50.5),
        mutation_scale=10.5,
        label_size=6.9,
        label_ha="left",
    )

    # Pathway membership spans carcinogen, enzyme, metabolite, and adduct layers.
    pathway_rail_x = 64.4
    pathway_join_y = (
        carcinogen_y + carcinogen_h - 4.0,
        enzyme_y + enzyme_h - 4.0,
        metabolite_y + metabolite_h - 4.0,
        adduct_y + adduct_h - 4.0,
    )
    for join_y in pathway_join_y:
        ax.plot(
            [62.1, pathway_rail_x],
            [join_y, join_y],
            color=_EDGE_TYPE_COLORS[EdgeType.PATHWAY.value],
            linewidth=1.5,
            linestyle=_EDGE_TYPE_STYLES[EdgeType.PATHWAY.value],
        )
    add_arrow(
        (pathway_rail_x, carcinogen_y + carcinogen_h - 2.4),
        (pathway_rail_x, pathway_y + 5.3),
        color=_EDGE_TYPE_COLORS[EdgeType.PATHWAY.value],
        label="PATHWAY",
        label_xy=(60.0, 50.0),
        linestyle=_EDGE_TYPE_STYLES[EdgeType.PATHWAY.value],
        mutation_scale=11.0,
        label_size=6.9,
        label_ha="right",
    )

    add_box(67, 38, 31, 51, face="#ffffff", edge="#d9e1ec", radius=1.2)
    add_text(82.5, 86.6, "Enzyme Categories", size=11.5, weight="bold", ha="center")
    enzyme_box_colors = ("#e8f1fb", "#ecf2f8", "#e3f7f4", "#f1e8f8")
    enzyme_box_y = (78.2, 66.6, 55.0, 43.4)
    for item, box_y, face in zip(data.enzyme_categories, enzyme_box_y, enzyme_box_colors, strict=True):
        add_box(69, box_y, 27, 9.6, face=face, edge="#d9e1ec", radius=0.8, linewidth=1.0)
        count_text = "" if item.count is None else f" (n={item.count})"
        add_text(71, box_y + 6.2, f"{item.title}{count_text}", size=8.8, color=item.color, weight="bold")
        preview = list(item.examples[:5])
        suffix = "..." if len(item.examples) > 5 else ""
        add_text(71, box_y + 2.7, _example_string(preview) + suffix, size=6.9, color="#667788")

    add_box(67, 19, 31, 17.5, face="#f5fbfc", edge="#158c98", radius=1.2, linewidth=1.5)
    add_text(82.5, 33.2, "Graph Summary", size=11.2, color="#157c86", weight="bold", ha="center")
    summary_y = 29.9
    for index, line in enumerate(data.summary_lines):
        add_text(82.5, summary_y - index * 3.1, line, size=7.8, color="#51606f", ha="center")

    add_box(2, 2.3, 96, 14.3, face="#ffffff", edge="#d9e1ec", radius=1.2)
    total_edges = sum(item.count for item in data.edge_legend)
    add_text(50, 15.2, f"{len(data.edge_legend)} Edge Types ({total_edges} Total Edges)", size=10.8, weight="bold", ha="center")

    x_positions = (10, 26, 42, 58, 74, 90)
    for x, legend_item in zip(x_positions, data.edge_legend, strict=True):
        ax.plot(
            [x - 4.2, x - 0.2],
            [10.0, 10.0],
            color=legend_item.color,
            linewidth=2.0,
            linestyle=legend_item.linestyle,
        )
        arrow = FancyArrowPatch(
            (x - 0.4, 10.0),
            (x + 0.6, 10.0),
            arrowstyle="->",
            mutation_scale=8,
            linewidth=1.6,
            color=legend_item.color,
        )
        ax.add_patch(arrow)
        label = (
            legend_item.label
            if legend_item.count <= 0
            else f"{legend_item.label} ({legend_item.count})"
        )
        add_text(x + 1.0, 10.2, label, size=7.0, color=legend_item.color, weight="bold")
        add_text(x + 1.0, 7.6, legend_item.description, size=5.9, color="#7b8795")

    fig.tight_layout()
    return fig


def save_architecture_figure(
    data: ArchitectureFigureData,
    path: str | Path,
    *,
    figsize: tuple[float, float] = (13.84, 10.0),
    dpi: int = 300,
) -> Path:
    """Render and save the architecture figure."""
    plt, _, _, _ = _require_matplotlib()
    figure = render_architecture_figure(data, figsize=figsize)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return path
