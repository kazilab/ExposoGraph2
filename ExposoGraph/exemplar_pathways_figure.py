"""Publication-style helper for Fig. 3 exemplar pathway panels."""

from __future__ import annotations

import math
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from .engine import GraphEngine
from .models import KnowledgeGraph, Node

_NODE_COLORS: dict[str, str] = {
    "Carcinogen": "#e05565",
    "Enzyme": "#4f98a3",
    "Gene": "#3d8b8b",
    "Metabolite": "#e8945a",
    "DNA_Adduct": "#a86fdf",
    "Pathway": "#5591c7",
    "Tissue": "#c2855a",
    "Note": "#f3eee7",
}

_EDGE_COLORS: dict[str, str] = {
    "ACTIVATES": "#e05565",
    "DETOXIFIES": "#6daa45",
    "FORMS_ADDUCT": "#a86fdf",
    "REPAIRS": "#e8af34",
    "CUSTOM": "#51606f",
}

_EDGE_STYLES: dict[str, str] = {
    "ACTIVATES": "-",
    "DETOXIFIES": "--",
    "FORMS_ADDUCT": "-",
    "REPAIRS": ":",
    "CUSTOM": "-",
}


@dataclass(frozen=True)
class _PanelItem:
    key: str
    label: str
    subtitle: str
    kind: str
    x: float
    y: float
    width: float
    height: float
    graph_ids: tuple[str, ...] = ()
    fill: str | None = None
    edge: str | None = None
    text_color: str | None = None
    alpha: float = 0.9

    @property
    def center(self) -> tuple[float, float]:
        return (self.x, self.y)

    @property
    def rx(self) -> float:
        return self.width / 2.0

    @property
    def ry(self) -> float:
        return self.height / 2.0


@dataclass(frozen=True)
class _PanelEdge:
    source: str
    target: str
    relation: str
    label: str = ""
    rad: float = 0.0
    label_offset: tuple[float, float] = (0.0, 0.0)
    alpha: float = 0.94
    linewidth: float = 2.1


@dataclass(frozen=True)
class _PanelSpec:
    letter: str
    title: str
    subtitle: str
    items: tuple[_PanelItem, ...]
    edges: tuple[_PanelEdge, ...]


def _require_matplotlib() -> tuple[Any, Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    except ImportError as exc:  # pragma: no cover - dependency failure path
        raise ImportError(
            "Exemplar pathway figure export requires `matplotlib`. "
            "Install `ExposoGraph[notebook]` or `pip install matplotlib`."
        ) from exc
    return plt, FancyArrowPatch, FancyBboxPatch


def _coerce_graph(graph: GraphEngine | KnowledgeGraph) -> KnowledgeGraph:
    if isinstance(graph, GraphEngine):
        return graph.to_knowledge_graph()
    return graph


def _node_lookup(graph: GraphEngine | KnowledgeGraph) -> dict[str, Node]:
    kg = _coerce_graph(graph)
    return {node.id: node for node in kg.nodes}


def _require_nodes(lookup: dict[str, Node], ids: Sequence[str], *, graph_name: str) -> None:
    missing = [node_id for node_id in ids if node_id not in lookup]
    if missing:
        raise ValueError(f"Missing required nodes in {graph_name}: {missing}")


def _wrap(text: str, width: int) -> str:
    return textwrap.fill(
        text,
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    )


def _anchor(source: _PanelItem, target: _PanelItem, *, padding: float = 0.45) -> tuple[float, float]:
    dx = target.x - source.x
    dy = target.y - source.y
    norm = math.hypot(dx, dy)
    if norm == 0:
        return source.center

    scale = 1.0 / math.sqrt((dx / max(source.rx, 0.01)) ** 2 + (dy / max(source.ry, 0.01)) ** 2)
    unit_x = dx / norm
    unit_y = dy / norm
    return (
        source.x + dx * scale + padding * unit_x,
        source.y + dy * scale + padding * unit_y,
    )


def _style_for_item(item: _PanelItem) -> tuple[str, str, str]:
    if item.kind == "Note":
        face = item.fill or "#f8efe5"
        edge = item.edge or "#e0c6ab"
        text = item.text_color or "#7d5a3d"
        return face, edge, text

    face = item.fill or _NODE_COLORS.get(item.kind, "#d9e4ef")
    edge = item.edge or face
    text = item.text_color or "white"
    return face, edge, text


def _panel_specs() -> tuple[_PanelSpec, ...]:
    return (
        _PanelSpec(
            letter="A",
            title="PAH: Benzo[a]pyrene",
            subtitle="Phase I activation, detoxification, and BPDE adduct formation",
            items=(
                _PanelItem("bap", "Benzo[a]pyrene", "Parent PAH", "Carcinogen", 50, 83, 25, 8, ("BaP",)),
                _PanelItem(
                    "cyp11",
                    "CYP1A1 / CYP1B1",
                    "Phase I activation",
                    "Enzyme",
                    17,
                    74,
                    19,
                    7,
                    ("CYP1A1", "CYP1B1"),
                ),
                _PanelItem("epoxide", "BaP-7,8-epoxide", "Reactive intermediate", "Metabolite", 50, 63, 24, 8, ("BaP_epoxide",)),
                _PanelItem("ephx1", "EPHX1", "Epoxide hydrolase", "Enzyme", 17, 53, 14, 7, ("EPHX1",)),
                _PanelItem("diol", "BaP-7,8-dihydrodiol", "Proximate carcinogen", "Metabolite", 50, 45, 28, 8, ("BaP_diol",)),
                _PanelItem("bpde", "BPDE", "Ultimate carcinogen", "Metabolite", 50, 27, 19, 8, ("BPDE",)),
                _PanelItem("gstm", "GSTM1 / GSTP1", "GSH conjugation", "Enzyme", 84, 33, 18, 7, ("GSTM1", "GSTP1")),
                _PanelItem("gsh", "BPDE-GSH", "Detoxified conjugate", "Metabolite", 84, 18, 18, 8, ("BPDE_GSH",)),
                _PanelItem("adduct", "BPDE-dG Adduct", "N2-guanine lesion", "DNA_Adduct", 50, 10, 22, 8, ("BPDE_dG",)),
                _PanelItem(
                    "risk",
                    "Variant context",
                    _wrap("CYP1A1 high-activity and GSTM1-null interaction tracked as research-use risk context.", 28),
                    "Note",
                    84,
                    52,
                    22,
                    11,
                    fill="#f8efe5",
                    edge="#e3b8a6",
                    text_color="#985140",
                    alpha=0.95,
                ),
            ),
            edges=(
                _PanelEdge("bap", "epoxide", "ACTIVATES", "1st oxidation", (0.0), (0.0, 4.6)),
                _PanelEdge("cyp11", "epoxide", "ACTIVATES"),
                _PanelEdge("epoxide", "diol", "CUSTOM", "hydrolysis", 0.0, (0.0, 4.2)),
                _PanelEdge("ephx1", "diol", "CUSTOM"),
                _PanelEdge("diol", "bpde", "ACTIVATES", "2nd oxidation", 0.0, (0.0, 4.2)),
                _PanelEdge("cyp11", "bpde", "ACTIVATES", rad=-0.18, alpha=0.5, linewidth=1.7),
                _PanelEdge("bpde", "adduct", "FORMS_ADDUCT"),
                _PanelEdge("bpde", "gsh", "DETOXIFIES", "detox", 0.0, (0.0, 3.6)),
                _PanelEdge("gstm", "gsh", "DETOXIFIES"),
            ),
        ),
        _PanelSpec(
            letter="B",
            title="Aromatic Amine: 4-ABP",
            subtitle="CYP1A2 activation with NAT1/NAT2 competition and aromatic-amine adduct formation",
            items=(
                _PanelItem("parent", "4-Aminobiphenyl", "Parent amine", "Carcinogen", 50, 83, 25, 8, ("4ABP",)),
                _PanelItem("cyp1a2", "CYP1A2", "N-oxidation", "Enzyme", 22, 73, 14, 7, ("CYP1A2",)),
                _PanelItem("n-oh", "N-hydroxy-4-ABP", "Proximate carcinogen", "Metabolite", 50, 61, 24, 8, ("NOH_4ABP",)),
                _PanelItem("nat2", "NAT2", "Parent-amine detox", "Enzyme", 20, 50, 12, 7, ("NAT2",)),
                _PanelItem("nat1", "NAT1", "Bladder activation", "Enzyme", 80, 48, 12, 7, ("NAT1",)),
                _PanelItem("gstm1", "GSTM1", "GSH conjugation", "Enzyme", 20, 29, 13, 7, ("GSTM1",)),
                _PanelItem("adduct", "dG-C8-4-ABP", "Aromatic-amine adduct", "DNA_Adduct", 78, 22, 21, 8, ("ABP_dG",)),
                _PanelItem(
                    "risk",
                    "Risk context",
                    _wrap("NAT2 slow-acetylator phenotypes are retained in ExposoGraph activity-score metadata.", 28),
                    "Note",
                    78,
                    9,
                    23,
                    10,
                    fill="#f8efe5",
                    edge="#e3b8a6",
                    text_color="#985140",
                    alpha=0.95,
                ),
            ),
            edges=(
                _PanelEdge("parent", "n-oh", "ACTIVATES", "N-oxidation", 0.0, (0.0, 4.4)),
                _PanelEdge("cyp1a2", "n-oh", "ACTIVATES"),
                _PanelEdge("nat2", "parent", "DETOXIFIES", "N-acetylation", 0.15, (-2.0, 4.2)),
                _PanelEdge("nat1", "n-oh", "ACTIVATES", "O-acetylation", -0.12, (0.0, 4.0)),
                _PanelEdge("gstm1", "n-oh", "DETOXIFIES", "GSH conjugation", 0.18, (0.0, -4.4)),
                _PanelEdge("n-oh", "adduct", "FORMS_ADDUCT", "C8-dG adduct", 0.0, (0.0, 4.3)),
            ),
        ),
        _PanelSpec(
            letter="C",
            title="Androgen Bridge: Testosterone",
            subtitle="Receptor signaling plus aromatase-linked estrogen quinone adduct formation",
            items=(
                _PanelItem("testosterone", "Testosterone", "Primary androgen", "Carcinogen", 50, 83, 22, 8, ("Testosterone",)),
                _PanelItem("srd5", "SRD5A1 / SRD5A2", "5alpha-reduction", "Enzyme", 20, 72, 19, 7, ("SRD5A1", "SRD5A2")),
                _PanelItem("dht", "5a-DHT", "Potent androgen", "Carcinogen", 20, 57, 18, 8, ("DHT",)),
                _PanelItem("ar", "Androgen receptor", "Gene locus", "Gene", 20, 41, 20, 7, ("AR",)),
                _PanelItem("signal", "AR Program", "Proliferative transcription", "Pathway", 20, 24, 21, 8, ("AR_signal_program",)),
                _PanelItem("cyp19", "CYP19A1", "Aromatase bridge", "Enzyme", 56, 72, 16, 7, ("CYP19A1",)),
                _PanelItem("e2", "17beta-Estradiol", "Estrogen bridge", "Carcinogen", 56, 57, 21, 8, ("E2",)),
                _PanelItem("cyp1b1", "CYP1B1", "4-hydroxylation", "Enzyme", 56, 42, 15, 7, ("CYP1B1",)),
                _PanelItem("hydroxye2", "4-Hydroxyestradiol", "Catechol estrogen", "Metabolite", 56, 28, 23, 8, ("HydroxyE2",)),
                _PanelItem("quinone", "Estradiol-3,4-quinone", "Electrophile", "Metabolite", 80, 28, 25, 8, ("E2_quinone",)),
                _PanelItem("adducts", "E2 DNA Adducts", "Depurinating adenine / guanine lesions", "DNA_Adduct", 80, 11, 24, 9, ("E2_Ade", "E2_Gua")),
                _PanelItem(
                    "variants",
                    "Variant context",
                    _wrap("SRD5A2 V89L / A49T, CYP19A1 repeat-length, and UGT2B17 deletion are represented in the androgen module.", 28),
                    "Note",
                    80,
                    57,
                    24,
                    14,
                    fill="#f3eef9",
                    edge="#d8c7eb",
                    text_color="#6d5094",
                    alpha=0.95,
                ),
            ),
            edges=(
                _PanelEdge("testosterone", "dht", "CUSTOM"),
                _PanelEdge("srd5", "dht", "CUSTOM"),
                _PanelEdge("dht", "ar", "CUSTOM", "high-affinity binding", 0.0, (0.0, 4.0)),
                _PanelEdge("ar", "signal", "CUSTOM", "AR signaling", 0.0, (0.0, 4.0)),
                _PanelEdge("testosterone", "e2", "CUSTOM"),
                _PanelEdge("cyp19", "e2", "CUSTOM"),
                _PanelEdge("e2", "hydroxye2", "ACTIVATES", "4-hydroxylation", 0.0, (0.0, 4.2)),
                _PanelEdge("cyp1b1", "hydroxye2", "ACTIVATES"),
                _PanelEdge("hydroxye2", "quinone", "ACTIVATES", "oxidation", 0.0, (0.0, 4.2)),
                _PanelEdge("quinone", "adducts", "FORMS_ADDUCT", "depurinating adducts", 0.0, (0.0, 4.0)),
            ),
        ),
        _PanelSpec(
            letter="D",
            title="Mycotoxin: Aflatoxin B1",
            subtitle="CYP3A4-driven epoxidation with GST detoxification and XPC repair support",
            items=(
                _PanelItem("afb1", "Aflatoxin B1", "Aspergillus mycotoxin", "Carcinogen", 50, 83, 23, 8, ("AFB1",)),
                _PanelItem("cyp", "CYP3A4 / CYP1A2", "8,9-epoxidation", "Enzyme", 22, 73, 20, 7, ("CYP3A4", "CYP1A2")),
                _PanelItem("epoxide", "AFB1-8,9-epoxide", "Highly reactive metabolite", "Metabolite", 50, 60, 25, 8, ("AFB1_epoxide",)),
                _PanelItem("gstm", "GSTM1 / GSTT1", "GSH conjugation", "Enzyme", 83, 58, 19, 7, ("GSTM1", "GSTT1")),
                _PanelItem("gsh", "AFB1-GSH", "Detoxified conjugate", "Metabolite", 83, 42, 18, 8, ("AFB1_GSH",)),
                _PanelItem("ephx1", "EPHX1", "Epoxide hydrolysis", "Enzyme", 18, 44, 14, 7, ("EPHX1",)),
                _PanelItem("adduct", "AFB1-N7-Gua", "Primary adduct", "DNA_Adduct", 50, 31, 20, 8, ("AFB1_Gua",)),
                _PanelItem("xpc", "XPC", "NER recognition", "Enzyme", 50, 15, 12, 7, ("XPC",)),
                _PanelItem(
                    "risk",
                    "Exposure context",
                    _wrap("Detoxification and repair context in the current graph centers on GSTM1 / GSTT1, EPHX1, ABCC2, and XPC.", 28),
                    "Note",
                    83,
                    15,
                    23,
                    12,
                    fill="#eef7f0",
                    edge="#c9e0cf",
                    text_color="#4c7b58",
                    alpha=0.95,
                ),
            ),
            edges=(
                _PanelEdge("afb1", "epoxide", "ACTIVATES", "8,9-epoxidation", 0.0, (0.0, 4.4)),
                _PanelEdge("cyp", "epoxide", "ACTIVATES"),
                _PanelEdge("epoxide", "gsh", "DETOXIFIES", "GSH conjugation", 0.0, (0.0, 4.2)),
                _PanelEdge("gstm", "gsh", "DETOXIFIES"),
                _PanelEdge("ephx1", "epoxide", "DETOXIFIES", "hydrolysis", 0.12, (-1.2, -4.0)),
                _PanelEdge("epoxide", "adduct", "FORMS_ADDUCT", "N7-guanine attack", 0.0, (0.0, 4.2)),
                _PanelEdge("xpc", "adduct", "REPAIRS", "repair support", 0.0, (0.0, 4.0)),
            ),
        ),
    )


def render_exemplar_pathways_figure(
    *,
    showcase_graph: GraphEngine | KnowledgeGraph | None = None,
    androgen_graph: GraphEngine | KnowledgeGraph | None = None,
    output_dir: str | Path | None = None,
    stem: str = "exemplar_pathways_figure",
    figsize: tuple[float, float] = (15.0, 11.8),
    dpi: int = 300,
    formats: Sequence[str] = ("png", "svg", "pdf"),
) -> tuple[Any, dict[str, Path]]:
    """Render a four-panel exemplar pathway figure for manuscript use."""
    if showcase_graph is None or androgen_graph is None:
        from .example_graphs import build_androgen_module_graph, build_full_legends_graph

        showcase_graph = showcase_graph or build_full_legends_graph()
        androgen_graph = androgen_graph or build_androgen_module_graph()

    showcase_lookup = _node_lookup(showcase_graph)
    androgen_lookup = _node_lookup(androgen_graph)

    _require_nodes(
        showcase_lookup,
        (
            "BaP",
            "BaP_epoxide",
            "BaP_diol",
            "BPDE",
            "BPDE_dG",
            "BPDE_GSH",
            "CYP1A1",
            "CYP1B1",
            "EPHX1",
            "GSTM1",
            "GSTP1",
            "4ABP",
            "CYP1A2",
            "NOH_4ABP",
            "NAT1",
            "NAT2",
            "ABP_dG",
            "AFB1",
            "CYP3A4",
            "AFB1_epoxide",
            "AFB1_GSH",
            "AFB1_Gua",
            "GSTT1",
            "XPC",
        ),
        graph_name="showcase graph",
    )
    _require_nodes(
        androgen_lookup,
        (
            "Testosterone",
            "DHT",
            "AR",
            "AR_signal_program",
            "SRD5A1",
            "SRD5A2",
            "CYP19A1",
            "E2",
            "CYP1B1",
            "HydroxyE2",
            "E2_quinone",
            "E2_Ade",
            "E2_Gua",
        ),
        graph_name="androgen graph",
    )

    plt, FancyArrowPatch, FancyBboxPatch = _require_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor("#f7f9fc")
    '''
    fig.suptitle(
        "Figure 3: Exemplar Carcinogen Metabolic Pathways",
        fontsize=18,
        fontweight="bold",
        color="#33475b",
        y=0.975,
    )
    fig.text(
        0.5,
        0.949,
        "Four manuscript-style panels rendered directly from the seeded ExposoGraph showcase and androgen-module graphs.",
        ha="center",
        fontsize=10,
        color="#617181",
    )
    '''

    def add_box(ax: Any, item: _PanelItem) -> None:
        face, edge, text = _style_for_item(item)
        rounding = min(item.height / 2.2, 3.8)
        if item.kind in {"Enzyme", "Gene"}:
            rounding = min(item.height / 1.5, 4.4)
        if item.kind == "Pathway":
            rounding = 2.4
        patch = FancyBboxPatch(
            (item.x - item.width / 2.0, item.y - item.height / 2.0),
            item.width,
            item.height,
            boxstyle=f"round,pad=0.35,rounding_size={rounding}",
            facecolor=face,
            edgecolor=edge,
            linewidth=1.2,
            alpha=item.alpha,
        )
        ax.add_patch(patch)
        if item.kind == "Note":
            note_top = item.y + (item.height / 2.0) - 0.9
            ax.text(
                item.x,
                note_top,
                item.label,
                ha="center",
                va="top",
                fontsize=7.2,
                fontweight="bold",
                color=text,
                linespacing=1.0,
            )
            if item.subtitle:
                ax.text(
                    item.x,
                    note_top - 3.7,
                    item.subtitle,
                    ha="center",
                    va="top",
                    fontsize=5.5,
                    color=text,
                    linespacing=1.05,
                )
            return
        label_y = item.y + (0.9 if item.subtitle else 0.0)
        ax.text(
            item.x,
            label_y,
            item.label,
            ha="center",
            va="center",
            fontsize=8.0 if item.kind != "Note" else 7.5,
            fontweight="bold",
            color=text,
            linespacing=1.0,
        )
        if item.subtitle:
            subtitle_color = "#f4f7fb"
            ax.text(
                item.x,
                item.y - 1.9,
                item.subtitle,
                ha="center",
                va="center",
                fontsize=5.8,
                color=subtitle_color,
                linespacing=1.05,
            )

    def draw_edge(ax: Any, item_by_key: dict[str, _PanelItem], edge: _PanelEdge) -> None:
        source = item_by_key[edge.source]
        target = item_by_key[edge.target]
        color = _EDGE_COLORS[edge.relation]
        patch = FancyArrowPatch(
            _anchor(source, target),
            _anchor(target, source),
            arrowstyle="-|>",
            mutation_scale=12.0,
            linewidth=edge.linewidth,
            linestyle=_EDGE_STYLES[edge.relation],
            color=color,
            connectionstyle=f"arc3,rad={edge.rad}",
            alpha=edge.alpha,
        )
        ax.add_patch(patch)
        if edge.label:
            ax.text(
                (source.x + target.x) / 2.0 + edge.label_offset[0],
                (source.y + target.y) / 2.0 + edge.label_offset[1],
                edge.label,
                ha="center",
                va="center",
                fontsize=6.2,
                color=color,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.92),
            )

    for ax, panel in zip(axes.flat, _panel_specs()):
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.axis("off")
        ax.set_facecolor("#f7f9fc")
        panel_bg = FancyBboxPatch(
            (1.0, 1.5),
            98.0,
            95.0,
            boxstyle="round,pad=0.9,rounding_size=4.5",
            facecolor="white",
            edgecolor="#dde6ef",
            linewidth=1.0,
        )
        ax.add_patch(panel_bg)
        ax.text(3.0, 94.3, panel.letter, fontsize=13, fontweight="bold", color="#33475b", ha="left", va="center")
        ax.text(8.0, 94.6, panel.title, fontsize=11.3, fontweight="bold", color="#e05565", ha="left", va="center")
        ax.text(8.0, 90.3, panel.subtitle, fontsize=7.0, color="#6d7d8d", ha="left", va="center")

        item_by_key = {item.key: item for item in panel.items}
        for item in panel.items:
            add_box(ax, item)
        for edge in panel.edges:
            draw_edge(ax, item_by_key, edge)

    fig.subplots_adjust(left=0.04, right=0.98, top=0.93, bottom=0.05, wspace=0.08, hspace=0.14)

    saved_paths: dict[str, Path] = {}
    if output_dir is not None:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        for fmt in formats:
            path = output_path / f"{stem}.{fmt}"
            fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
            saved_paths[str(fmt)] = path

    return fig, saved_paths
