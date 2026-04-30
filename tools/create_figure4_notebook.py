"""Create the Figure 4 glutathione-depletion notebook."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = ROOT / "Figures_Notebook" / "Figure-4-glutathione-depletion-tipping-point.ipynb"


def code_cell(source: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(source.strip() + "\n")


def markdown_cell(source: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(source.strip() + "\n")


def build_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
    }

    nb.cells = [
        markdown_cell(
            """
# Figure 4: Glutathione depletion tipping point

This notebook regenerates the manuscript-style glutathione depletion curve from
the local **ExposoGraph** package.

Package-alignment notes:

- Baseline hepatic GSH, synthesis rate, half-life, and the critical threshold are
  read from `ExposoGraph/data/interaction_parameters.json`.
- The wild-type curve is computed with the same steady-state equation used by
  `ExposoGraph.interaction_engine.gsh_depletion_model()`.
- GCLC and GCLM curves are synthesis-capacity scenarios. They are explicit
  extensions of the same ExposoGraph equation because the public API does not
  currently accept GCLC/GCLM genotypes as direct inputs.
- ExposoGraph uses 20% of baseline GSH as the formal critical threshold and
  flags general GST impairment below 30%.
            """
        ),
        code_cell(
            """
from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Make the notebook robust whether it is launched from repo root or Figures_Notebook/.
for candidate in [Path.cwd(), *Path.cwd().parents]:
    if (candidate / "ExposoGraph").is_dir() and (candidate / "pyproject.toml").exists():
        REPO_ROOT = candidate
        break
else:
    raise RuntimeError("Could not locate the ExposoGraph repository root.")

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

OUTPUT_DIR = REPO_ROOT / "Figures_Notebook" / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PARAMETER_FILE = REPO_ROOT / "ExposoGraph" / "data" / "interaction_parameters.json"
with PARAMETER_FILE.open("r", encoding="utf-8") as handle:
    INTERACTION_PARAMS = json.load(handle)

GSH_PARAMS = INTERACTION_PARAMS["gsh_depletion"]

BASELINE_GSH_MM = float(GSH_PARAMS["baseline_gsh_mM"])
SYNTHESIS_RATE = float(GSH_PARAMS["synthesis_rate_umol_h_g"])
HALF_LIFE_H = float(GSH_PARAMS.get("half_life_h", 2.5))
CRITICAL_THRESHOLD = float(GSH_PARAMS["critical_threshold_fraction"])
IMPAIRMENT_THRESHOLD = 0.30
LIVER_WATER_FRACTION = 0.70

# The manuscript x-axis is relative exposure load. The scale is calibrated so
# that wild-type GSH reaches the 30% impairment threshold at 3 relative units.
REFERENCE_IMPAIRMENT_LOAD = 3.0
NATURAL_TURNOVER = BASELINE_GSH_MM * LIVER_WATER_FRACTION * math.log(2) / HALF_LIFE_H
CONSUMPTION_PER_LOAD_UNIT = (
    (1.0 - IMPAIRMENT_THRESHOLD) * (SYNTHESIS_RATE + NATURAL_TURNOVER)
    / REFERENCE_IMPAIRMENT_LOAD
)

LOAD_GRID = np.linspace(0.0, 5.0, 251)
            """
        ),
        markdown_cell(
            """
## ExposoGraph GSH equation

`gsh_depletion_model()` converts exposure input to a total GSH consumption rate.
For positive net synthesis, ExposoGraph computes steady-state GSH as:

`fraction_normal = 1 - consumption / (synthesis_rate + natural_turnover)`

If consumption equals or exceeds synthesis, the model marks the tipping point as
reached. The helper below keeps that equation unchanged while allowing synthesis
capacity to be scaled for the GCLC/GCLM scenarios.
            """
        ),
        code_cell(
            """
from ExposoGraph.interaction_engine import gsh_depletion_model


@dataclass(frozen=True)
class GSHCurvePoint:
    load: float
    gsh_remaining_percent: float
    consumption_umol_h_g: float
    synthesis_umol_h_g: float
    net_rate_umol_h_g: float
    tipping_point_reached: bool
    time_to_critical_depletion_h: float | None


def compute_gsh_point(load: float, synthesis_scale: float = 1.0) -> GSHCurvePoint:
    baseline_umol_g = BASELINE_GSH_MM * LIVER_WATER_FRACTION
    synthesis = SYNTHESIS_RATE * synthesis_scale
    consumption = load * CONSUMPTION_PER_LOAD_UNIT
    net_rate = synthesis - consumption

    if consumption <= 0:
        fraction_normal = 1.0
    elif net_rate > 0:
        fraction_normal = max(0.0, 1.0 - (consumption / (synthesis + NATURAL_TURNOVER)))
    else:
        fraction_normal = 0.0

    time_to_depletion = None
    if consumption > synthesis:
        gsh_to_lose = baseline_umol_g * (1.0 - CRITICAL_THRESHOLD)
        time_to_depletion = round(gsh_to_lose / (consumption - synthesis), 2)

    return GSHCurvePoint(
        load=round(load, 4),
        gsh_remaining_percent=round(fraction_normal * 100.0, 3),
        consumption_umol_h_g=round(consumption, 4),
        synthesis_umol_h_g=round(synthesis, 4),
        net_rate_umol_h_g=round(net_rate, 4),
        tipping_point_reached=bool(consumption >= synthesis),
        time_to_critical_depletion_h=time_to_depletion,
    )


# Validate that the helper matches gsh_depletion_model() for the wild-type case
# at selected loads when both receive the same direct GSH-consumption rate.
for load in [0.0, 1.0, 2.0, 3.0]:
    local = compute_gsh_point(load, synthesis_scale=1.0)
    package = gsh_depletion_model({"combined_umol_h_g": local.consumption_umol_h_g})
    assert abs(local.gsh_remaining_percent / 100.0 - package.fraction_normal) <= 0.001

print("Validated wild-type helper against ExposoGraph gsh_depletion_model().")
print(f"Baseline GSH: {BASELINE_GSH_MM} mM")
print(f"Synthesis rate: {SYNTHESIS_RATE} umol/h/g liver")
print(f"Half-life: {HALF_LIFE_H} h")
print(f"Critical threshold: {CRITICAL_THRESHOLD:.0%}")
print(f"Impairment threshold used for figure annotation: {IMPAIRMENT_THRESHOLD:.0%}")
print(f"Consumption per x-axis load unit: {CONSUMPTION_PER_LOAD_UNIT:.3f} umol/h/g")
            """
        ),
        markdown_cell(
            """
## Compute curves and export source data
            """
        ),
        code_cell(
            """
SCENARIOS = [
    {
        "scenario": "Wild-type",
        "label": "Wild-type",
        "synthesis_scale": 1.00,
        "color": "#007c89",
        "linestyle": "-",
        "linewidth": 2.2,
    },
    {
        "scenario": "GCLC variant",
        "label": "GCLC variant (70% synthesis)",
        "synthesis_scale": 0.70,
        "color": "#d89c00",
        "linestyle": "--",
        "linewidth": 2.0,
    },
    {
        "scenario": "GCLM variant",
        "label": "GCLM variant (50% synthesis)",
        "synthesis_scale": 0.50,
        "color": "#b65a32",
        "linestyle": "-.",
        "linewidth": 2.0,
    },
]

curve_rows = []
for scenario in SCENARIOS:
    for load in LOAD_GRID:
        point = compute_gsh_point(float(load), float(scenario["synthesis_scale"]))
        curve_rows.append(
            {
                "scenario": scenario["scenario"],
                "label": scenario["label"],
                "synthesis_scale": scenario["synthesis_scale"],
                "relative_exposure_load": point.load,
                "gsh_remaining_percent": point.gsh_remaining_percent,
                "consumption_umol_h_g": point.consumption_umol_h_g,
                "synthesis_umol_h_g": point.synthesis_umol_h_g,
                "net_rate_umol_h_g": point.net_rate_umol_h_g,
                "tipping_point_reached": point.tipping_point_reached,
                "time_to_critical_depletion_h": point.time_to_critical_depletion_h,
            }
        )

curve_csv = OUTPUT_DIR / "figure4_gsh_depletion_tipping_point_curves.csv"
with curve_csv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(curve_rows[0].keys()))
    writer.writeheader()
    writer.writerows(curve_rows)

threshold_rows = []
for scenario in SCENARIOS:
    scale = float(scenario["synthesis_scale"])
    synthesis = SYNTHESIS_RATE * scale
    denominator = synthesis + NATURAL_TURNOVER
    impairment_load = (1.0 - IMPAIRMENT_THRESHOLD) * denominator / CONSUMPTION_PER_LOAD_UNIT
    critical_load = (1.0 - CRITICAL_THRESHOLD) * denominator / CONSUMPTION_PER_LOAD_UNIT
    synthesis_exceeded_load = synthesis / CONSUMPTION_PER_LOAD_UNIT
    threshold_rows.append(
        {
            "scenario": scenario["scenario"],
            "synthesis_scale": scale,
            "load_at_30_percent_gsh": round(impairment_load, 3),
            "load_at_20_percent_gsh": round(critical_load, 3),
            "load_when_consumption_equals_synthesis": round(synthesis_exceeded_load, 3),
        }
    )

threshold_csv = OUTPUT_DIR / "figure4_gsh_depletion_thresholds.csv"
with threshold_csv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(threshold_rows[0].keys()))
    writer.writeheader()
    writer.writerows(threshold_rows)

print(curve_csv)
print(threshold_csv)
for row in threshold_rows:
    print(row)
            """
        ),
        markdown_cell(
            """
## Render Figure 4
            """
        ),
        code_cell(
            """
fig, ax = plt.subplots(figsize=(8.2, 5.35), constrained_layout=False)
fig.subplots_adjust(left=0.12, right=0.96, bottom=0.32, top=0.86)

ax.axhspan(0, IMPAIRMENT_THRESHOLD * 100, color="#f4d7d2", alpha=0.38, zorder=0)
ax.axhline(IMPAIRMENT_THRESHOLD * 100, color="#b44b42", linestyle=":", linewidth=1.1)
ax.axhline(CRITICAL_THRESHOLD * 100, color="#8d1f1f", linestyle=":", linewidth=0.9, alpha=0.75)

for scenario in SCENARIOS:
    rows = [row for row in curve_rows if row["scenario"] == scenario["scenario"]]
    loads = [row["relative_exposure_load"] for row in rows]
    gsh = [row["gsh_remaining_percent"] for row in rows]
    ax.plot(
        loads,
        gsh,
        label=scenario["label"],
        color=scenario["color"],
        linestyle=scenario["linestyle"],
        linewidth=scenario["linewidth"],
    )

threshold_lookup = {row["scenario"]: row for row in threshold_rows}
for scenario in SCENARIOS:
    row = threshold_lookup[scenario["scenario"]]
    impairment_load = row["load_at_30_percent_gsh"]
    if impairment_load <= LOAD_GRID.max():
        ax.plot(
            impairment_load,
            IMPAIRMENT_THRESHOLD * 100,
            marker="o",
            markersize=4.5,
            color=scenario["color"],
            zorder=5,
        )
        ax.annotate(
            scenario["scenario"],
            xy=(impairment_load, IMPAIRMENT_THRESHOLD * 100),
            xytext=(impairment_load + 0.08, IMPAIRMENT_THRESHOLD * 100 + 7),
            fontsize=6.2,
            color=scenario["color"],
            arrowprops={
                "arrowstyle": "-",
                "color": scenario["color"],
                "linewidth": 0.7,
            },
        )

ax.annotate(
    "GST impairment zone (<30%)",
    xy=(3.72, IMPAIRMENT_THRESHOLD * 100),
    xytext=(3.72, IMPAIRMENT_THRESHOLD * 100 + 9),
    fontsize=6.4,
    color="#9b302b",
    ha="center",
    arrowprops={"arrowstyle": "-", "color": "#9b302b", "linewidth": 0.7},
)
ax.text(
    4.18,
    CRITICAL_THRESHOLD * 100 - 3.0,
    "critical threshold\\n(20%)",
    fontsize=6.0,
    color="#8d1f1f",
    ha="center",
    va="top",
)

scenario_annotations = [
    (0.15, "Reference\\nexposure"),
    (1.0, "Single\\ncarcinogen"),
    (2.0, "Moderate\\nco-exposure"),
    (3.0, "GST impairment\\nthreshold"),
    (4.35, "Heavy combined\\noccupational"),
]
for x, label in scenario_annotations:
    ax.text(
        x,
        -0.19,
        label,
        ha="center",
        va="top",
        fontsize=6.0,
        color="#555555",
        transform=ax.get_xaxis_transform(),
        clip_on=False,
    )
    ax.plot(
        [x, x],
        [0.0, -0.035],
        color="#bdbdbd",
        linewidth=0.7,
        transform=ax.get_xaxis_transform(),
        clip_on=False,
    )

ax.set_xlim(0, 5)
ax.set_ylim(0, 105)
ax.set_xticks(np.arange(0, 5.1, 1.0))
ax.set_yticks(np.arange(0, 101, 20))
ax.set_xlabel("Combined GSH-consuming exposure load (relative units)", fontsize=8, labelpad=9)
ax.set_ylabel("Intracellular GSH remaining (%)", fontsize=8)
ax.tick_params(axis="both", labelsize=7)
ax.grid(axis="both", color="#e1e1e1", linewidth=0.55, alpha=0.75)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

ax.legend(loc="upper right", frameon=False, fontsize=7)
fig.suptitle(
    "Glutathione Depletion Under Multi-Carcinogen Exposure",
    fontsize=11,
    weight="bold",
    y=0.965,
)
fig.text(
    0.5,
    0.912,
    "Tipping-point model based on ExposoGraph interaction_engine GSH parameters",
    ha="center",
    va="bottom",
    fontsize=7.5,
)
fig.text(
    0.12,
    0.035,
    "Data source: ExposoGraph/data/interaction_parameters.json and "
    "interaction_engine.gsh_depletion_model(). Baseline hepatic GSH = "
    f"{BASELINE_GSH_MM:g} mM; synthesis rate = {SYNTHESIS_RATE:g} umol/h/g liver; "
    f"half-life = {HALF_LIFE_H:g} h. GCLC/GCLM curves are synthesis-capacity "
    "scenarios using the same equation.",
    fontsize=5.7,
    ha="left",
    va="bottom",
    color="#555555",
)

png_path = OUTPUT_DIR / "figure4_gsh_depletion_tipping_point.png"
pdf_path = OUTPUT_DIR / "figure4_gsh_depletion_tipping_point.pdf"
svg_path = OUTPUT_DIR / "figure4_gsh_depletion_tipping_point.svg"
fig.savefig(png_path, dpi=300, bbox_inches="tight")
fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(svg_path, bbox_inches="tight")
plt.show()

print(png_path)
print(pdf_path)
print(svg_path)
            """
        ),
    ]
    return nb


def main() -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    notebook = build_notebook()
    with NOTEBOOK_PATH.open("w", encoding="utf-8") as handle:
        nbf.write(notebook, handle)
    print(NOTEBOOK_PATH)


if __name__ == "__main__":
    main()
