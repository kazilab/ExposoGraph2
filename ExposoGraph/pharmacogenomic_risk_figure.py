"""Publication-style helpers for pharmacogenomic risk-scoring summaries.

This figure is intentionally conservative. It visualizes the
pharmacogenomic support ExposoGraph already has today:

- representative activity-score annotations on a subset of enzymes
- graph topology that distinguishes activation, detoxification, and repair
- a topology-aware per-gene impact score

It does not claim a full patient-level composite risk engine.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .engine import GraphEngine
from .graph_analysis import variant_impact_score
from .models import KnowledgeGraph, NodeType

_GROUP_ORDER: tuple[str, ...] = (
    "PAH",
    "HCA",
    "Aromatic_Amine",
    "Nitrosamine",
    "Mycotoxin",
    "Estrogen",
    "Androgen",
    "Solvent",
    "Alkylating",
)

_GROUP_LABELS: dict[str, str] = {
    "PAH": "PAH",
    "HCA": "HCA",
    "Aromatic_Amine": "Aromatic amines",
    "Nitrosamine": "Nitrosamines",
    "Mycotoxin": "Mycotoxins",
    "Estrogen": "Estrogens",
    "Androgen": "Androgens",
    "Solvent": "Solvents",
    "Alkylating": "Alkylating agents",
}

_ROLE_COLORS: dict[str, str] = {
    "Activation": "#d65c52",
    "Detoxification": "#5f9f4a",
    "Mixed": "#c89b3f",
    "Repair": "#7f64c1",
    "Other": "#617181",
}

_TEXT_DARK = "#33475b"
_TEXT_MID = "#617181"
_GRID = "#dbe6ef"


@dataclass(frozen=True)
class PharmacogenomicRiskGeneProfile:
    gene_id: str
    label: str
    phase: str | None
    role: str | None
    role_bucket: str
    activity_score: float
    impact_score: float
    downstream_adduct_count: int
    downstream_repair_count: int
    carcinogen_groups: tuple[str, ...]
    class_count: int


@dataclass(frozen=True)
class PharmacogenomicRiskClassProfile:
    carcinogen_group: str
    display_label: str
    activation_score: float
    detoxification_score: float
    mixed_score: float
    repair_score: float
    scored_gene_count: int
    genes: tuple[str, ...]


def _require_matplotlib() -> tuple[Any, Any]:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:  # pragma: no cover - dependency failure path
        raise ImportError(
            "Pharmacogenomic risk figure export requires `matplotlib`. "
            "Install `ExposoGraph[notebook]` or `pip install matplotlib`."
        ) from exc
    return plt, Line2D


def _coerce_engine(graph: GraphEngine | KnowledgeGraph) -> GraphEngine:
    if isinstance(graph, GraphEngine):
        return graph
    engine = GraphEngine()
    for node in graph.nodes:
        engine.add_node(node)
    for edge in graph.edges:
        engine.add_edge(edge)
    return engine


def _role_bucket(node_data: Mapping[str, Any]) -> str:
    role = str(node_data.get("role") or "").strip()
    phase = str(node_data.get("phase") or "").strip()
    if role == "Activation":
        return "Activation"
    if role == "Mixed":
        return "Mixed"
    if role in {"Detoxification", "Transport"} or phase in {"II", "III"}:
        return "Detoxification"
    if role == "Repair" or phase == "Repair":
        return "Repair"
    return "Other"


def _scored_enzyme_ids(engine: GraphEngine) -> list[str]:
    scored: list[str] = []
    for node_id, data in engine.G.nodes(data=True):
        if data.get("type") != NodeType.ENZYME.value:
            continue
        if data.get("activity_score") is None:
            continue
        scored.append(str(node_id))
    return sorted(scored)


def _group_scope_node_ids(engine: GraphEngine, carcinogen_group: str) -> set[str]:
    carcinogen_ids = {
        str(node_id)
        for node_id, data in engine.G.nodes(data=True)
        if data.get("type") == NodeType.CARCINOGEN.value
        and data.get("group") == carcinogen_group
    }
    if not carcinogen_ids:
        return set()

    relevant_node_ids = set(carcinogen_ids)

    for source, target, edge_data in engine.G.edges(data=True):
        if (
            str(source) in carcinogen_ids
            or str(target) in carcinogen_ids
            or edge_data.get("carcinogen") in carcinogen_ids
        ):
            relevant_node_ids.update({str(source), str(target)})

    primary_scope = set(relevant_node_ids)
    for source, target, edge_data in engine.G.edges(data=True):
        if edge_data.get("type") != "FORMS_ADDUCT":
            continue
        if str(source) in primary_scope or str(target) in primary_scope:
            relevant_node_ids.update({str(source), str(target)})

    adduct_scope = set(relevant_node_ids)
    for source, target, edge_data in engine.G.edges(data=True):
        if edge_data.get("type") != "REPAIRS":
            continue
        if str(source) in adduct_scope or str(target) in adduct_scope:
            relevant_node_ids.update({str(source), str(target)})

    return relevant_node_ids


def build_pharmacogenomic_risk_gene_profiles(
    graph: GraphEngine | KnowledgeGraph,
) -> list[PharmacogenomicRiskGeneProfile]:
    """Return scored enzyme profiles for the supplied graph."""
    engine = _coerce_engine(graph)

    group_gene_ids = {
        group: _group_scope_node_ids(engine, group)
        for group in _GROUP_ORDER
    }

    profiles: list[PharmacogenomicRiskGeneProfile] = []
    for gene_id in _scored_enzyme_ids(engine):
        node_data = engine.get_data(gene_id) or {}
        impact = variant_impact_score(engine, gene_id)
        carcinogen_groups = tuple(
            group
            for group in _GROUP_ORDER
            if gene_id in group_gene_ids[group]
        )
        profiles.append(
            PharmacogenomicRiskGeneProfile(
                gene_id=gene_id,
                label=str(node_data.get("label") or gene_id),
                phase=str(node_data.get("phase")) if node_data.get("phase") else None,
                role=str(node_data.get("role")) if node_data.get("role") else None,
                role_bucket=_role_bucket(node_data),
                activity_score=float(node_data["activity_score"]),
                impact_score=float(impact.score if impact is not None else 0.0),
                downstream_adduct_count=int(
                    impact.downstream_adduct_count if impact is not None else 0
                ),
                downstream_repair_count=int(
                    impact.downstream_repair_count if impact is not None else 0
                ),
                carcinogen_groups=carcinogen_groups,
                class_count=len(carcinogen_groups),
            )
        )

    return sorted(
        profiles,
        key=lambda item: (
            -item.impact_score,
            item.activity_score,
            item.gene_id,
        ),
    )


def build_pharmacogenomic_risk_class_profiles(
    graph: GraphEngine | KnowledgeGraph,
) -> list[PharmacogenomicRiskClassProfile]:
    """Summarize scored activation, detoxification, mixed, and repair support."""
    engine = _coerce_engine(graph)
    scored_gene_ids = set(_scored_enzyme_ids(engine))

    profiles: list[PharmacogenomicRiskClassProfile] = []
    for group in _GROUP_ORDER:
        node_ids = _group_scope_node_ids(engine, group)
        totals = {
            "Activation": 0.0,
            "Detoxification": 0.0,
            "Mixed": 0.0,
            "Repair": 0.0,
            "Other": 0.0,
        }
        genes: list[str] = []
        for gene_id in sorted(node_ids.intersection(scored_gene_ids)):
            node_data = engine.get_data(gene_id) or {}
            bucket = _role_bucket(node_data)
            totals[bucket] += float(node_data["activity_score"])
            genes.append(gene_id)

        profiles.append(
            PharmacogenomicRiskClassProfile(
                carcinogen_group=group,
                display_label=_GROUP_LABELS.get(group, group.replace("_", " ")),
                activation_score=round(totals["Activation"], 3),
                detoxification_score=round(totals["Detoxification"], 3),
                mixed_score=round(totals["Mixed"], 3),
                repair_score=round(totals["Repair"], 3),
                scored_gene_count=len(genes),
                genes=tuple(genes),
            )
        )

    return profiles


def pharmacogenomic_risk_gene_rows(
    profiles: Iterable[PharmacogenomicRiskGeneProfile],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for profile in profiles:
        row = asdict(profile)
        row["carcinogen_groups"] = ", ".join(profile.carcinogen_groups)
        rows.append(row)
    return rows


def pharmacogenomic_risk_class_rows(
    profiles: Iterable[PharmacogenomicRiskClassProfile],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for profile in profiles:
        row = asdict(profile)
        row["genes"] = ", ".join(profile.genes)
        rows.append(row)
    return rows


def render_pharmacogenomic_risk_figure(
    graph: GraphEngine | KnowledgeGraph,
    *,
    output_dir: str | Path | None = None,
    stem: str = "pharmacogenomic_risk_scoring_figure",
    figsize: tuple[float, float] = (16.0, 9.2),
    dpi: int = 300,
    formats: Sequence[str] = ("png", "svg", "pdf"),
) -> tuple[Any, dict[str, Path], list[PharmacogenomicRiskGeneProfile], list[PharmacogenomicRiskClassProfile]]:
    """Render a publication-style summary of current PGx risk-scoring support."""
    plt, Line2D = _require_matplotlib()
    engine = _coerce_engine(graph)
    gene_profiles = build_pharmacogenomic_risk_gene_profiles(engine)
    class_profiles = build_pharmacogenomic_risk_class_profiles(engine)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor("#f8fafc")
    grid = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.24), wspace=0.18)
    ax_class = fig.add_subplot(grid[0, 0])
    ax_gene = fig.add_subplot(grid[0, 1])

    for ax in (ax_class, ax_gene):
        ax.set_facecolor("#ffffff")
        for spine in ax.spines.values():
            spine.set_color("#d7e3ec")
            spine.set_linewidth(1.0)

    nonzero_class_profiles = [
        profile for profile in class_profiles if profile.scored_gene_count > 0
    ]
    if not nonzero_class_profiles:
        raise ValueError("The supplied graph has no scored enzyme profiles to plot.")

    max_activation = max(profile.activation_score for profile in nonzero_class_profiles)
    max_detox = max(profile.detoxification_score for profile in nonzero_class_profiles)
    right_text_margin = 2.3
    ax_class.set_xlim(
        -(max_detox + 0.9),
        max_activation + right_text_margin,
    )
    ax_class.axvline(0, color="#9cb2c6", linewidth=1.0)

    y_positions = list(range(len(nonzero_class_profiles)))
    ax_class.barh(
        y_positions,
        [-profile.detoxification_score for profile in nonzero_class_profiles],
        color=_ROLE_COLORS["Detoxification"],
        height=0.56,
        alpha=0.92,
    )
    ax_class.barh(
        y_positions,
        [profile.activation_score for profile in nonzero_class_profiles],
        color=_ROLE_COLORS["Activation"],
        height=0.56,
        alpha=0.92,
    )
    ax_class.set_yticks(y_positions)
    ax_class.set_yticklabels(
        [profile.display_label for profile in nonzero_class_profiles],
        fontsize=10.2,
        color=_TEXT_DARK,
    )
    ax_class.invert_yaxis()
    ax_class.grid(axis="x", color=_GRID, linewidth=0.9, alpha=0.7)
    ax_class.set_axisbelow(True)
    ax_class.tick_params(axis="x", colors=_TEXT_MID, labelsize=9.2)
    ax_class.tick_params(axis="y", length=0)
    ax_class.set_xlabel(
        "Representative activity-score reservoir by carcinogen class",
        fontsize=9.6,
        color=_TEXT_MID,
    )
    ax_class.set_title(
        "Class-level activation and detoxification support",
        fontsize=10,
        fontweight="bold",
        color=_TEXT_DARK,
        loc="left",
        pad=12,
    )
    ax_class.text(
        0.03,
        1.003,
        "Bars summarize currently scored enzymes in each graph-defined class scope.",
        transform=ax_class.transAxes,
        fontsize=7,
        color=_TEXT_MID,
        ha="left",
        va="bottom",
    )

    for idx, profile in enumerate(nonzero_class_profiles):
        if profile.activation_score > 0:
            ax_class.text(
                profile.activation_score + 0.08,
                idx,
                f"{profile.activation_score:.1f}",
                va="center",
                ha="left",
                fontsize=8.8,
                color=_TEXT_DARK,
            )
        if profile.detoxification_score > 0:
            ax_class.text(
                -profile.detoxification_score - 0.08,
                idx,
                f"{profile.detoxification_score:.1f}",
                va="center",
                ha="right",
                fontsize=8.8,
                color=_TEXT_DARK,
            )
        badges: list[str] = []
        if profile.mixed_score > 0:
            badges.append(f"M {profile.mixed_score:.1f}")
        if profile.repair_score > 0:
            badges.append(f"R {profile.repair_score:.1f}")
        if badges:
            ax_class.text(
                max_activation + 0.28,
                idx,
                "  ".join(badges),
                va="center",
                ha="left",
                fontsize=8.8,
                color=_TEXT_MID,
            )

    ax_class.text(
        0.11,
        0.02,
        "Detoxification",
        transform=ax_class.transAxes,
        fontsize=9.0,
        color=_ROLE_COLORS["Detoxification"],
        fontweight="bold",
        ha="center",
    )
    ax_class.text(
        0.73,
        0.02,
        "Activation",
        transform=ax_class.transAxes,
        fontsize=9.0,
        color=_ROLE_COLORS["Activation"],
        fontweight="bold",
        ha="center",
    )

    max_impact = max(profile.impact_score for profile in gene_profiles)
    ax_gene.set_xlim(-0.02, 2.25)
    ax_gene.set_ylim(-0.12, max_impact + 0.75)
    ax_gene.grid(color=_GRID, linewidth=0.9, alpha=0.7)
    ax_gene.set_axisbelow(True)
    ax_gene.tick_params(colors=_TEXT_MID, labelsize=9.2)
    ax_gene.axvline(1.0, color="#9cb2c6", linewidth=1.0, linestyle="--", alpha=0.85)
    ax_gene.axhline(1.0, color="#c1ced9", linewidth=1.0, linestyle=":", alpha=0.75)
    ax_gene.set_xlabel("Representative activity score", fontsize=10.0, color=_TEXT_MID)
    ax_gene.set_ylabel("Topology-aware gene impact score", fontsize=10.0, color=_TEXT_MID)
    ax_gene.set_title(
        "Per-gene impact within the current reference graph",
        fontsize=10,
        fontweight="bold",
        color=_TEXT_DARK,
        loc="left",
        pad=12,
    )
    ax_gene.text(
        0.03,
        1.003,
        "Bubble size tracks how many carcinogen classes each scored gene touches.",
        transform=ax_gene.transAxes,
        fontsize=7,
        color=_TEXT_MID,
        ha="left",
        va="bottom",
    )

    offsets = (
        (8, 8),
        (10, -12),
        (-52, 8),
        (-52, -12),
        (10, 17),
        (-58, 17),
    )
    for idx, gene_profile in enumerate(gene_profiles):
        point_size = 90 + (gene_profile.class_count * 36)
        ax_gene.scatter(
            gene_profile.activity_score,
            gene_profile.impact_score,
            s=point_size,
            color=_ROLE_COLORS.get(gene_profile.role_bucket, _ROLE_COLORS["Other"]),
            alpha=0.92,
            edgecolor="#ffffff",
            linewidth=1.2,
            zorder=3,
        )
        dx, dy = offsets[idx % len(offsets)]
        ax_gene.annotate(
            gene_profile.gene_id,
            (gene_profile.activity_score, gene_profile.impact_score),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8.6,
            color=_TEXT_DARK,
            bbox={
                "boxstyle": "round,pad=0.18",
                "fc": "#ffffff",
                "ec": "#dde7ef",
                "lw": 0.7,
                "alpha": 0.94,
            },
        )

    legend_items = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=_ROLE_COLORS[bucket],
            markeredgecolor="#ffffff",
            markeredgewidth=1.0,
            markersize=8.5,
            label=bucket,
        )
        for bucket in ("Activation", "Detoxification", "Mixed", "Repair")
    ]
    ax_gene.legend(
        handles=legend_items,
        frameon=False,
        fontsize=9.0,
        labelcolor=_TEXT_DARK,
        loc="upper right",
    )

    scored_class_count = sum(
        1 for profile in class_profiles if profile.scored_gene_count > 0
    )
    fig.suptitle(
        "Pharmacogenomic Risk Scoring Foundation in ExposoGraph",
        fontsize=18,
        fontweight="bold",
        color=_TEXT_DARK,
        y=0.98,
    )
    fig.text(
        0.5,
        0.947,
        (
            f"Current reference graph: {engine.G.number_of_nodes()} nodes, "
            f"{engine.G.number_of_edges()} edges, {len(gene_profiles)} scored enzymes "
            f"across {scored_class_count} carcinogen classes."
        ),
        ha="center",
        va="center",
        fontsize=10.0,
        color=_TEXT_MID,
    )
    fig.text(
        0.5,
        0.022,
        (
            "Interpretation: lower activity combined with broader downstream adduct topology yields a higher "
            "gene-level impact score. This figure illustrates current structural support, not a patient-level classifier."
        ),
        ha="center",
        va="center",
        fontsize=9.0,
        color=_TEXT_MID,
    )
    fig.subplots_adjust(top=0.87, bottom=0.1, left=0.07, right=0.98)

    outputs: dict[str, Path] = {}
    if output_dir is not None:
        output_root = Path(output_dir)
        output_root.mkdir(parents=True, exist_ok=True)
        for fmt in formats:
            path = output_root / f"{stem}.{fmt}"
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            outputs[fmt] = path

    return fig, outputs, gene_profiles, class_profiles
