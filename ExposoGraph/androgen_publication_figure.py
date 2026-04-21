"""Publication-figure helpers for the seeded androgen module."""

from __future__ import annotations

import math
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from .engine import GraphEngine
from .models import EdgeType, KnowledgeGraph

_NODE_COLORS: dict[str, str] = {
    "Carcinogen": "#e05565",
    "Enzyme": "#4f98a3",
    "Gene": "#3d8b8b",
    "Metabolite": "#e8945a",
    "DNA_Adduct": "#a86fdf",
    "Pathway": "#5591c7",
    "Tissue": "#c2855a",
}

_EDGE_COLORS: dict[str, str] = {
    "ACTIVATES": "#e05565",
    "DETOXIFIES": "#6daa45",
    "FORMS_ADDUCT": "#a86fdf",
    "CUSTOM": "#51606f",
}

_EDGE_STYLES: dict[str, str] = {
    "ACTIVATES": "-",
    "DETOXIFIES": "--",
    "FORMS_ADDUCT": "-",
    "CUSTOM": "-",
}


@dataclass(frozen=True)
class _NodeSpec:
    x: float
    y: float
    label: str
    font_size: float
    radius: float | None = None
    width: float | None = None
    height: float | None = None

    @property
    def center(self) -> tuple[float, float]:
        return (self.x, self.y)

    @property
    def rx(self) -> float:
        if self.width is not None:
            return self.width / 2.0
        return float(self.radius or 0.0)

    @property
    def ry(self) -> float:
        if self.height is not None:
            return self.height / 2.0
        return float(self.radius or 0.0)


_NODE_SPECS: dict[str, _NodeSpec] = {
    "Androstenedione": _NodeSpec(16.0, 66.0, "Androstenedione", 7.1, radius=4.05),
    "Testosterone": _NodeSpec(34.5, 66.0, "Testosterone", 7.3, radius=4.15),
    "DHT": _NodeSpec(52.5, 66.0, "5a-DHT", 7.3, radius=4.05),
    "AR": _NodeSpec(70.3, 66.0, "Androgen receptor\n(AR)", 5.8, width=12.4, height=6.4),
    "AR_signal_program": _NodeSpec(88.7, 66.0, "AR proliferative\nprogram", 6.0, width=11.6, height=6.4),
    "E2": _NodeSpec(34.5, 34.2, "17beta-\nEstradiol", 6.0, width=8.4, height=6.2),
    "HydroxyE2": _NodeSpec(52.5, 34.2, "4-Hydroxy-\nestradiol", 5.8, width=10.8, height=6.2),
    "E2_quinone": _NodeSpec(70.3, 34.2, "Estradiol-3,4-\nquinone", 5.7, width=11.4, height=6.2),
    "E2_Ade": _NodeSpec(88.7, 42.2, "4-OHE2-1-\nN3Ade", 5.7, width=9.6, height=6.1),
    "E2_Gua": _NodeSpec(88.7, 26.0, "4-OHE2-guanine\ndepurinating\nadduct", 4.9, width=13.4, height=6.6),
    "HydroxyTestosterone": _NodeSpec(34.5, 86.2, "6beta-Hydroxy-\ntestosterone", 5.7, width=11.8, height=6.2),
    "Testosterone_gluc": _NodeSpec(52.5, 86.2, "Testosterone-17-\nglucuronide", 5.7, width=12.6, height=6.2),
    "DHT_gluc": _NodeSpec(70.3, 86.2, "DHT-17-\nglucuronide", 5.9, width=10.9, height=6.2),
}

_TISSUE_IDS: tuple[str, ...] = (
    "Prostate",
    "PeripheralHormoneTissues",
    "HormoneResponsiveTissues",
)

_VARIANT_ROWS: tuple[tuple[str, str, str], ...] = (
    ("SRD5A2_V89L", "SRD5A2 V89L", "Reduced 5alpha-reductase"),
    ("SRD5A2_A49T", "SRD5A2 A49T", "Increased 5alpha-reductase"),
    ("CYP19A1_repeat_length", "CYP19A1 repeat-length", "Aromatase activity shift"),
    ("UGT2B17_copy_number_deletion", "UGT2B17 CN deletion", "Reduced glucuronidation"),
)


def _require_matplotlib() -> tuple[Any, Any, Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
    except ImportError as exc:  # pragma: no cover - dependency failure path
        raise ImportError(
            "Androgen publication figure export requires `matplotlib`. "
            "Install `ExposoGraph[notebook]` or `pip install matplotlib`."
        ) from exc
    return plt, Circle, FancyArrowPatch, FancyBboxPatch


def _coerce_graph(graph: GraphEngine | KnowledgeGraph) -> KnowledgeGraph:
    if isinstance(graph, GraphEngine):
        return graph.to_knowledge_graph()
    return graph


def _summary_counts(showcase_summary: Any | None) -> tuple[int | None, int | None]:
    if showcase_summary is None:
        return (None, None)
    node_count = getattr(showcase_summary, "node_count", None)
    edge_count = getattr(showcase_summary, "edge_count", None)
    return (
        int(node_count) if node_count is not None else None,
        int(edge_count) if edge_count is not None else None,
    )


def _wrap_text(text: str, width: int) -> str:
    return textwrap.fill(
        text,
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    )


def _display_gene_label(label: str) -> str:
    return "AR" if label == "Androgen receptor (AR)" else label


def _anchor(spec: _NodeSpec, target: tuple[float, float], *, padding: float = 0.45) -> tuple[float, float]:
    dx = target[0] - spec.x
    dy = target[1] - spec.y
    norm = math.hypot(dx, dy)
    if norm == 0:
        return spec.center

    scale = 1.0 / math.sqrt((dx / max(spec.rx, 0.01)) ** 2 + (dy / max(spec.ry, 0.01)) ** 2)
    unit_x = dx / norm
    unit_y = dy / norm
    return (
        spec.x + dx * scale + padding * unit_x,
        spec.y + dy * scale + padding * unit_y,
    )


def render_androgen_publication_figure(
    graph: GraphEngine | KnowledgeGraph,
    *,
    output_dir: str | Path | None = None,
    stem: str = "androgen_module_publication_figure",
    showcase_summary: Any | None = None,
    figsize: tuple[float, float] = (17.0, 10.0),
    dpi: int = 300,
    formats: Sequence[str] = ("png", "svg", "pdf"),
) -> tuple[Any, dict[str, Path]]:
    """Render the androgen-module publication figure and optionally save it."""
    plt, Circle, FancyArrowPatch, FancyBboxPatch = _require_matplotlib()

    kg = _coerce_graph(graph)
    node_lookup = {node.id: node for node in kg.nodes}
    required_ids = set(_NODE_SPECS) | set(_TISSUE_IDS) | {row[0] for row in _VARIANT_ROWS}
    missing = sorted(required_ids.difference(node_lookup))
    if missing:
        raise ValueError(f"Missing required androgen-module nodes: {missing}")

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    fig.patch.set_facecolor("#f8fafc")
    ax.set_facecolor("#f8fafc")

    def add_box(
        x: float,
        y: float,
        w: float,
        h: float,
        face: str,
        *,
        edge: str = "#d8e2eb",
        radius: float = 2.4,
        lw: float = 1.2,
        alpha: float = 1.0,
    ) -> Any:
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle=f"round,pad=0.5,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            alpha=alpha,
        )
        ax.add_patch(patch)
        return patch

    add_box(4.5, 49.6, 91.0, 25.0, "#eef7fb", edge="#c9dce8", radius=3.0)
    add_box(21.5, 15.8, 74.0, 25.2, "#fff4ef", edge="#ead7c8", radius=3.0)
    add_box(21.5, 77.6, 57.0, 15.6, "#eef8f0", edge="#d4e7d8", radius=3.0)
    add_box(5.8, 15.8, 14.6, 31.0, "#fffdf5", edge="#eadfb4", radius=2.5)
    add_box(79.0, 77.6, 17.0, 16.4, "#f6f1fb", edge="#dfd2ef", radius=2.5)

    ax.text(
        50,
        96.0,
        "Androgen Metabolism Module in ExposoGraph",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
        color="#33475b",
        family="DejaVu Sans",
    )
    ax.text(
        50,
        92.6,
        "Receptor-mediated proliferation, aromatase bridging, detoxification, and variant context rendered from the seeded androgen module",
        ha="center",
        va="center",
        fontsize=9.6,
        color="#617181",
        family="DejaVu Sans",
    )
    ax.text(
        50,
        74.9,
        "Receptor-Mediated Proliferation Axis",
        ha="center",
        va="center",
        fontsize=12.2,
        fontweight="bold",
        color="#35546b",
        family="DejaVu Sans",
    )
    ax.text(
        58.5,
        43.0,
        "Aromatase Bridge to Estrogen-DNA Adduct Formation",
        ha="center",
        va="center",
        fontsize=11.7,
        fontweight="bold",
        color="#8a5031",
        family="DejaVu Sans",
    )
    ax.text(
        50,
        90.8,
        "Clearance and Competing Metabolism",
        ha="center",
        va="center",
        fontsize=11.7,
        fontweight="bold",
        color="#3e6f4f",
        family="DejaVu Sans",
    )
    ax.text(
        13.1,
        43.8,
        "Tissue Context",
        ha="center",
        va="center",
        fontsize=11.4,
        fontweight="bold",
        color="#8b6f18",
        family="DejaVu Sans",
    )
    ax.text(
        87.5,
        89.0,
        "Variant Context",
        ha="center",
        va="center",
        fontsize=11.1,
        fontweight="bold",
        color="#7552a3",
        family="DejaVu Sans",
    )

    for node_id, spec in _NODE_SPECS.items():
        node = node_lookup[node_id]
        node_type = node.type.value
        color = _NODE_COLORS[node_type]
        if spec.width is not None and spec.height is not None:
            rounding = 1.5 if node_type == "DNA_Adduct" else min(spec.height / 2.0, 3.0)
            shape = FancyBboxPatch(
                (spec.x - spec.width / 2.0, spec.y - spec.height / 2.0),
                spec.width,
                spec.height,
                boxstyle=f"round,pad=0.25,rounding_size={rounding}",
                facecolor=color,
                edgecolor=color,
                linewidth=1.4,
                alpha=0.88,
            )
            ax.add_patch(shape)
        else:
            ax.add_patch(
                Circle(
                    spec.center,
                    spec.radius or 3.4,
                    facecolor=color,
                    edgecolor=color,
                    alpha=0.84,
                    linewidth=1.4,
                )
            )
        ax.text(
            spec.x,
            spec.y,
            spec.label,
            ha="center",
            va="center",
            fontsize=spec.font_size,
            color="white",
            fontweight="bold",
            linespacing=0.95,
            family="DejaVu Sans",
        )

    def draw_edge(
        source: str,
        target: str,
        edge_type: str,
        label: str,
        *,
        rad: float = 0.0,
        offset: tuple[float, float] = (0.0, 0.0),
        alpha: float = 0.95,
        linewidth: float = 2.1,
        mutation_scale: float = 14.0,
        label_size: float = 8.0,
    ) -> None:
        source_spec = _NODE_SPECS[source]
        target_spec = _NODE_SPECS[target]
        patch = FancyArrowPatch(
            _anchor(source_spec, target_spec.center),
            _anchor(target_spec, source_spec.center),
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=linewidth,
            linestyle=_EDGE_STYLES[edge_type],
            color=_EDGE_COLORS[edge_type],
            connectionstyle=f"arc3,rad={rad}",
            alpha=alpha,
        )
        ax.add_patch(patch)
        ax.text(
            (source_spec.x + target_spec.x) / 2.0 + offset[0],
            (source_spec.y + target_spec.y) / 2.0 + offset[1],
            label,
            ha="center",
            va="center",
            fontsize=label_size,
            color=_EDGE_COLORS[edge_type],
            fontweight="bold",
            family="DejaVu Sans",
            bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="none", alpha=0.92),
        )

    draw_edge("Androstenedione", "Testosterone", "CUSTOM", "AKR1C3", offset=(0, 4.5), label_size=8.2)
    draw_edge("Testosterone", "DHT", "CUSTOM", "SRD5A1 / SRD5A2", offset=(0, 4.4), label_size=8.0)
    draw_edge("DHT", "AR", "CUSTOM", "High-affinity binding", offset=(0, 4.4), label_size=7.9)
    draw_edge("AR", "AR_signal_program", "CUSTOM", "AR signaling", offset=(0, 4.3), label_size=7.8)
    draw_edge(
        "Testosterone",
        "AR",
        "CUSTOM",
        "Lower-affinity binding",
        rad=-0.32,
        offset=(0, -10.8),
        alpha=0.72,
        label_size=7.7,
    )
    draw_edge("Testosterone", "E2", "CUSTOM", "CYP19A1", offset=(-4.1, 1.5), label_size=7.8)
    draw_edge("E2", "HydroxyE2", "ACTIVATES", "CYP1B1", offset=(0, 4.2), label_size=8.0)
    draw_edge("HydroxyE2", "E2_quinone", "ACTIVATES", "Oxidation", offset=(0, 4.2), label_size=8.0)
    draw_edge(
        "E2_quinone",
        "E2_Ade",
        "FORMS_ADDUCT",
        "Depurinating\nadenine adduct",
        rad=0.14,
        offset=(8.3, 4.6),
        label_size=7.3,
    )
    draw_edge(
        "E2_quinone",
        "E2_Gua",
        "FORMS_ADDUCT",
        "Depurinating\nguanine adduct",
        rad=-0.08,
        offset=(8.6, -1.2),
        label_size=7.3,
    )
    draw_edge("Testosterone", "HydroxyTestosterone", "DETOXIFIES", "CYP3A5", offset=(-2.3, 2.8), label_size=7.6)
    draw_edge("Testosterone", "Testosterone_gluc", "DETOXIFIES", "UGT2B7", offset=(0, 4.4), label_size=7.7)
    draw_edge("DHT", "DHT_gluc", "DETOXIFIES", "UGT2B17", offset=(0, 4.4), label_size=7.7)

    ax.text(
        16.0,
        72.7,
        "CYP17A1-supported\nandrogen precursor branch",
        ha="center",
        va="bottom",
        fontsize=7.6,
        color="#4f98a3",
        family="DejaVu Sans",
    )
    ax.annotate(
        "",
        xy=_anchor(_NODE_SPECS["Androstenedione"], (16.0, 78.4), padding=0.0),
        xytext=(16.0, 78.4),
        arrowprops=dict(arrowstyle="-|>", color="#4f98a3", lw=1.8),
    )

    expression_edges = [edge for edge in kg.edges if edge.type == EdgeType.EXPRESSED_IN]
    tissue_map: dict[str, list[str]] = {}
    for edge in expression_edges:
        tissue_map.setdefault(edge.target, []).append(_display_gene_label(node_lookup[edge.source].label))

    tissue_y = 39.8
    for tissue_id in _TISSUE_IDS:
        tissue = node_lookup[tissue_id]
        genes = ", ".join(sorted(tissue_map.get(tissue_id, [])))
        ax.text(
            8.0,
            tissue_y + 2.4,
            tissue.label,
            ha="left",
            va="center",
            fontsize=8.3,
            fontweight="bold",
            color="#7f6717",
            family="DejaVu Sans",
        )
        ax.text(
            8.0,
            tissue_y - 0.8,
            _wrap_text(genes, 28),
            ha="left",
            va="top",
            fontsize=6.9,
            color="#5f6f7d",
            linespacing=1.15,
            family="DejaVu Sans",
        )
        tissue_y -= 10.1

    variant_y = 86.6
    for variant_id, label, phenotype in _VARIANT_ROWS:
        if variant_id not in node_lookup:
            continue
        ax.text(
            80.7,
            variant_y,
            label,
            ha="left",
            va="center",
            fontsize=7.4,
            fontweight="bold",
            color="#6d5094",
            family="DejaVu Sans",
        )
        ax.text(
            80.7,
            variant_y - 1.85,
            phenotype,
            ha="left",
            va="center",
            fontsize=6.2,
            color="#5f6f7d",
            family="DejaVu Sans",
        )
        variant_y -= 3.55

    legend_y = 8.5
    legend_items = (
        ("ACTIVATES", "Activation"),
        ("DETOXIFIES", "Detoxification"),
        ("FORMS_ADDUCT", "DNA adduct formation"),
        ("CUSTOM", "Curated signaling or conversion"),
    )
    for index, (edge_type, label) in enumerate(legend_items):
        x = 12 + index * 22
        ax.plot([x, x + 5], [legend_y, legend_y], color=_EDGE_COLORS[edge_type], lw=2.2, ls=_EDGE_STYLES[edge_type])
        ax.text(
            x + 6.5,
            legend_y,
            label,
            ha="left",
            va="center",
            fontsize=7.6,
            color="#51606f",
            family="DejaVu Sans",
        )

    module_node_count = len(kg.nodes)
    module_edge_count = len(kg.edges)
    showcase_node_count, showcase_edge_count = _summary_counts(showcase_summary)
    summary = f"Androgen module: {module_node_count} nodes · {module_edge_count} edges"
    if showcase_node_count is not None and showcase_edge_count is not None:
        summary += f"  |  Integrated showcase with module: {showcase_node_count} nodes · {showcase_edge_count} edges"
    ax.text(
        50,
        4.1,
        summary,
        ha="center",
        va="center",
        fontsize=8.2,
        color="#5a6a78",
        family="DejaVu Sans",
    )

    saved_paths: dict[str, Path] = {}
    if output_dir is not None:
        target_dir = Path(output_dir)
        target_dir.mkdir(parents=True, exist_ok=True)
        for suffix in formats:
            path = target_dir / f"{stem}.{suffix}"
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            saved_paths[suffix] = path

    return fig, saved_paths


__all__ = ["render_androgen_publication_figure"]
