"""Create the Figure 6 variant-sensitivity notebook."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = (
    ROOT
    / "Figures_Notebook"
    / "Figure-6-variant-contributions-sensitivity-analysis.ipynb"
)


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
# Figure 6: Variant contribution sensitivity analysis

This notebook regenerates a code-aligned Figure 6 from the local
**ExposoGraph** package.

Package-alignment notes:

- Synthetic genotypes are sampled from `ExposoGraph.population_simulation._core`
  using the bundled ancestry-specific genotype tables and haplotype blocks.
- The simulation uses 10,000 synthetic individuals for each of five ancestry
  strata: EUR, EAS, AFR, AMR, and SAS.
- Flux profiles are computed with `ExposoGraph.flux_engine.compute_full_profile()`
  under the same Tier 1/background, liver-context convention used for Figure 5.
- The composite score is the same relative index used in Figure 5: the
  geometric mean of class-wise median-normalized activation/detoxification net
  ratios.
- Variant contribution is estimated by feature neutralization. For each
  modeled variant/phenotype, carriers are reset to the package reference state,
  the ExposoGraph flux profile is recomputed, and the reduction in variance of
  the log2 composite score is recorded. Bars are normalized across modeled
  features with positive variance reduction.

The older draft image included variants such as NQO1*2, UGT1A1*28, XRCC1
R399Q, OGG1 S326C, ABCG2 Q141K, ADH1B*2, and MGMT promoter methylation. Those
features are not currently sampled by the ExposoGraph population simulation
module, so they are exported to an audit CSV but excluded from the code-derived
bar chart.
            """
        ),
        code_cell(
            """
from __future__ import annotations

import csv
import math
import random
import statistics
import sys
import warnings
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

for candidate in [Path.cwd(), *Path.cwd().parents]:
    if (candidate / "ExposoGraph").is_dir() and (candidate / "pyproject.toml").exists():
        REPO_ROOT = candidate
        break
else:
    raise RuntimeError("Could not locate the ExposoGraph repository root.")

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ExposoGraph.flux_engine import compute_full_profile
from ExposoGraph.population_simulation._core import (
    ALLELE_FREQUENCIES,
    HAPLOTYPE_BLOCKS,
    _genes_covered_by_blocks,
    _weighted_pick,
    sample_haplotype_block,
)

OUTPUT_DIR = REPO_ROOT / "Figures_Notebook" / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

N_PER_ANCESTRY = 10_000
SEED = 20260430
TISSUE = "Liver"
EPS = 1e-9

ANCESTRIES = OrderedDict(
    [
        ("EUR", "European"),
        ("EAS", "East Asian"),
        ("AFR", "African"),
        ("AMR", "Admixed American"),
        ("SAS", "South Asian"),
    ]
)

DISPLAY_NAMES = {
    "EUR": "European",
    "EAS": "East Asian",
    "AFR": "African/African American",
    "AMR": "Latino/Admixed",
    "SAS": "South Asian",
}

SAMPLED_GENES = tuple(
    sorted(
        set(ALLELE_FREQUENCIES)
        | {gene for block in HAPLOTYPE_BLOCKS.values() for gene in block["genes"]}
    )
)

# Only these sampled genes currently change compute_full_profile() outputs in
# the population-flux workflow. This key keeps profile caching fast while
# preserving exact model behavior for the active sampled genes.
PROFILE_GENES = (
    "ALDH2",
    "CYP1A1",
    "CYP1A2",
    "CYP2E1",
    "GSTM1",
    "GSTT1",
    "NAT1",
    "NAT2",
)

CATEGORY_COLORS = {
    "Phase I activation": "#007c89",
    "Aldehyde detox": "#d89c00",
    "Phase II acetylation": "#7b5ea7",
    "Phase II GST": "#b65a32",
}

plt.rcParams.update(
    {
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
    }
)

# CYP1A2 *1A/*1F is sampled by the haplotype table but currently defaults to
# normal activity in the package modifier. Suppress repeated warnings so the
# notebook output stays readable.
warnings.filterwarnings("ignore", message="Unrecognized diplotype")

print(f"Repository root: {REPO_ROOT}")
print(f"Sampled genes: {', '.join(SAMPLED_GENES)}")
print(f"Profile-active sampled genes: {', '.join(PROFILE_GENES)}")
            """
        ),
        markdown_cell(
            """
## Simulate genotype profiles

The sampling function mirrors the Figure 5 notebook: haplotype blocks are drawn
first, then independently sampled genes not covered by those blocks are added.
            """
        ),
        code_cell(
            """
def sample_genotypes_for_ancestry(rng: random.Random, ancestry_label: str) -> dict[str, str]:
    genotypes: dict[str, str] = {}
    covered_genes = _genes_covered_by_blocks()

    for block_name in HAPLOTYPE_BLOCKS:
        genotypes.update(sample_haplotype_block(rng, block_name, ancestry_label))

    for gene, per_ancestry in ALLELE_FREQUENCIES.items():
        if gene in covered_genes:
            continue
        genotypes[gene] = _weighted_pick(rng, per_ancestry[ancestry_label])

    return genotypes


profile_cache: dict[tuple[tuple[str, str], ...], dict[str, float]] = {}


def profile_key(genotypes: dict[str, str]) -> tuple[tuple[str, str], ...]:
    return tuple((gene, str(genotypes.get(gene, ""))) for gene in PROFILE_GENES)


def profile_ratios(genotypes: dict[str, str]) -> dict[str, float]:
    key = profile_key(genotypes)
    if key not in profile_cache:
        profile_input = {gene: str(genotypes.get(gene, "")) for gene in PROFILE_GENES}
        profile = compute_full_profile(profile_input, tissue=TISSUE)
        profile_cache[key] = {
            class_name: float(result.net_ratio)
            for class_name, result in profile.per_class_results.items()
        }
    return profile_cache[key]


rng = random.Random(SEED)
rows: list[dict[str, object]] = []

for ancestry_code, ancestry_label in ANCESTRIES.items():
    for index in range(N_PER_ANCESTRY):
        genotypes = sample_genotypes_for_ancestry(rng, ancestry_label)
        ratios = profile_ratios(genotypes)
        rows.append(
            {
                "person_id": f"{ancestry_code}{index:05d}",
                "ancestry_code": ancestry_code,
                "ancestry_label": ancestry_label,
                "display_name": DISPLAY_NAMES[ancestry_code],
                "genotypes": genotypes,
                "ratios": ratios,
            }
        )

print(f"Simulated {len(rows):,} individuals.")
print(f"Unique cached profile states: {len(profile_cache):,}")
print(f"Modeled carcinogen classes: {len(rows[0]['ratios'])}")
            """
        ),
        markdown_cell(
            """
## Compute the relative composite score

Each carcinogen class is median-normalized across the pooled simulated cohort.
The composite score is the geometric mean of the normalized class-specific net
ratios. Variance is computed on the log2 score because the score is
multiplicative.
            """
        ),
        code_cell(
            """
class_names = list(rows[0]["ratios"])

class_medians: dict[str, float] = {}
for class_name in class_names:
    values = [
        float(row["ratios"][class_name])
        for row in rows
        if float(row["ratios"][class_name]) > 0
    ]
    class_medians[class_name] = statistics.median(values) if values else 1.0


def composite_score_from_ratios(ratios: dict[str, float]) -> float:
    components = []
    for class_name in class_names:
        ratio = max(float(ratios[class_name]), EPS)
        median = max(float(class_medians[class_name]), EPS)
        components.append(math.log(ratio / median))
    return math.exp(sum(components) / len(components))


def variance(values: list[float]) -> float:
    if not values:
        return 0.0
    mean = sum(values) / len(values)
    return sum((value - mean) ** 2 for value in values) / len(values)


for row in rows:
    score = composite_score_from_ratios(row["ratios"])
    row["relative_composite_score"] = score
    row["relative_composite_log2"] = math.log2(score)

baseline_log2_scores = [float(row["relative_composite_log2"]) for row in rows]
baseline_variance = variance(baseline_log2_scores)

print(f"Baseline log2 composite variance: {baseline_variance:.6f}")
print(
    "Composite score range: "
    f"{min(row['relative_composite_score'] for row in rows):.3f} to "
    f"{max(row['relative_composite_score'] for row in rows):.3f}"
)
            """
        ),
        markdown_cell(
            """
## Feature-neutralization variance analysis

For each feature, the notebook changes only individuals carrying that feature to
the package reference state, recomputes the ExposoGraph profile, and measures
how much log2 composite-score variance is removed. Because variants can be
correlated through ancestry and haplotype blocks, raw reductions are normalized
to sum to 100% across positive model-supported features for plotting.
            """
        ),
        code_cell(
            """
FEATURES = [
    {
        "key": "GSTM1_null",
        "label": "GSTM1 null",
        "category": "Phase II GST",
        "condition": lambda gt: gt.get("GSTM1") == "null",
        "neutralize": lambda gt: {**gt, "GSTM1": "present"},
        "reference_state": "GSTM1=present",
    },
    {
        "key": "ALDH2_star2",
        "label": "ALDH2*2",
        "category": "Aldehyde detox",
        "condition": lambda gt: "*2" in str(gt.get("ALDH2", "")),
        "neutralize": lambda gt: {**gt, "ALDH2": "*1/*1"},
        "reference_state": "ALDH2=*1/*1",
    },
    {
        "key": "CYP1A1_star2A",
        "label": "CYP1A1*2A",
        "category": "Phase I activation",
        "condition": lambda gt: "*2A" in str(gt.get("CYP1A1", "")),
        "neutralize": lambda gt: {**gt, "CYP1A1": "WT"},
        "reference_state": "CYP1A1=WT",
    },
    {
        "key": "CYP1A2_1F_homozygote",
        "label": "CYP1A2*1F/*1F",
        "category": "Phase I activation",
        "condition": lambda gt: gt.get("CYP1A2") == "*1F/*1F",
        "neutralize": lambda gt: {**gt, "CYP1A2": "*1A/*1A"},
        "reference_state": "CYP1A2=*1A/*1A",
    },
    {
        "key": "NAT2_slow",
        "label": "NAT2 slow",
        "category": "Phase II acetylation",
        "condition": lambda gt: gt.get("NAT2") == "slow",
        "neutralize": lambda gt: {**gt, "NAT2": "rapid"},
        "reference_state": "NAT2=rapid",
    },
    {
        "key": "NAT2_intermediate",
        "label": "NAT2 intermediate",
        "category": "Phase II acetylation",
        "condition": lambda gt: gt.get("NAT2") == "intermediate",
        "neutralize": lambda gt: {**gt, "NAT2": "rapid"},
        "reference_state": "NAT2=rapid",
    },
    {
        "key": "CYP2E1_UM",
        "label": "CYP2E1 UM",
        "category": "Phase I activation",
        "condition": lambda gt: gt.get("CYP2E1") == "UM_c1c1",
        "neutralize": lambda gt: {**gt, "CYP2E1": "NM"},
        "reference_state": "CYP2E1=NM",
    },
    {
        "key": "CYP2E1_IM",
        "label": "CYP2E1 IM",
        "category": "Phase I activation",
        "condition": lambda gt: gt.get("CYP2E1") == "IM",
        "neutralize": lambda gt: {**gt, "CYP2E1": "NM"},
        "reference_state": "CYP2E1=NM",
    },
    {
        "key": "GSTT1_null",
        "label": "GSTT1 null",
        "category": "Phase II GST",
        "condition": lambda gt: gt.get("GSTT1") == "null",
        "neutralize": lambda gt: {**gt, "GSTT1": "present"},
        "reference_state": "GSTT1=present",
    },
]


contribution_rows: list[dict[str, object]] = []

for feature in FEATURES:
    modified_log2_scores: list[float] = []
    carrier_count = 0
    for row in rows:
        genotypes = row["genotypes"]
        if feature["condition"](genotypes):
            carrier_count += 1
            neutralized = feature["neutralize"](genotypes)
            ratios = profile_ratios(neutralized)
            score = composite_score_from_ratios(ratios)
            modified_log2_scores.append(math.log2(score))
        else:
            modified_log2_scores.append(float(row["relative_composite_log2"]))

    neutralized_variance = variance(modified_log2_scores)
    raw_reduction = baseline_variance - neutralized_variance
    positive_reduction = max(0.0, raw_reduction)
    contribution_rows.append(
        {
            "feature_key": feature["key"],
            "label": feature["label"],
            "category": feature["category"],
            "reference_state": feature["reference_state"],
            "carrier_count": carrier_count,
            "carrier_frequency_percent": round(100 * carrier_count / len(rows), 3),
            "baseline_log2_variance": round(baseline_variance, 8),
            "neutralized_log2_variance": round(neutralized_variance, 8),
            "raw_variance_reduction": round(raw_reduction, 8),
            "positive_variance_reduction": positive_reduction,
        }
    )

positive_total = sum(float(row["positive_variance_reduction"]) for row in contribution_rows)
for row in contribution_rows:
    row["normalized_contribution_percent"] = (
        round(100 * float(row["positive_variance_reduction"]) / positive_total, 3)
        if positive_total > 0
        else 0.0
    )

contribution_rows.sort(key=lambda row: float(row["normalized_contribution_percent"]), reverse=True)

for row in contribution_rows:
    print(
        f"{row['label']:<18} "
        f"freq={row['carrier_frequency_percent']:>6.2f}% "
        f"contribution={row['normalized_contribution_percent']:>6.2f}%"
    )
            """
        ),
        markdown_cell(
            """
## Export source tables
            """
        ),
        code_cell(
            """
contribution_csv = OUTPUT_DIR / "figure6_variant_variance_contributions.csv"
with contribution_csv.open("w", newline="", encoding="utf-8") as handle:
    fieldnames = [
        "feature_key",
        "label",
        "category",
        "reference_state",
        "carrier_count",
        "carrier_frequency_percent",
        "baseline_log2_variance",
        "neutralized_log2_variance",
        "raw_variance_reduction",
        "positive_variance_reduction",
        "normalized_contribution_percent",
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(contribution_rows)

population_csv = OUTPUT_DIR / "figure6_variant_sensitivity_population_scores.csv"
with population_csv.open("w", newline="", encoding="utf-8") as handle:
    fieldnames = [
        "person_id",
        "ancestry_code",
        "ancestry_label",
        "display_name",
        *SAMPLED_GENES,
        "relative_composite_score",
        "relative_composite_log2",
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        genotypes = row["genotypes"]
        writer.writerow(
            {
                "person_id": row["person_id"],
                "ancestry_code": row["ancestry_code"],
                "ancestry_label": row["ancestry_label"],
                "display_name": row["display_name"],
                **{gene: genotypes.get(gene, "") for gene in SAMPLED_GENES},
                "relative_composite_score": round(float(row["relative_composite_score"]), 6),
                "relative_composite_log2": round(float(row["relative_composite_log2"]), 6),
            }
        )

audit_rows = [
    {
        "legacy_feature": "NQO1*2",
        "status": "not_sampled_by_population_module",
        "note": (
            "NQO1 is used by the benzene flux model but not currently sampled "
            "in ALLELE_FREQUENCIES or HAPLOTYPE_BLOCKS."
        ),
    },
    {
        "legacy_feature": "GSTP1 Ile105Val",
        "status": "not_sampled_by_population_module",
        "note": "GSTP1 is used by flux models but not currently sampled in the population module.",
    },
    {
        "legacy_feature": "UGT1A1*28",
        "status": "not_active_in_current_population_flux_profile",
        "note": (
            "UGT1A1 is present in reference annotations but is not sampled in "
            "the current population module or active in compute_full_profile()."
        ),
    },
    {
        "legacy_feature": "XRCC1 R399Q",
        "status": "not_sampled_by_population_module",
        "note": (
            "XRCC1 is used by selected proxy repair terms but not sampled in "
            "the population module."
        ),
    },
    {
        "legacy_feature": "OGG1 S326C",
        "status": "not_sampled_by_population_module",
        "note": (
            "OGG1 is used by heavy-metal oxidative repair terms but not sampled "
            "in the population module."
        ),
    },
    {
        "legacy_feature": "ABCG2 Q141K",
        "status": "not_active_in_current_population_flux_profile",
        "note": (
            "ABCG2 is present in reference annotations but not sampled or active "
            "in compute_full_profile()."
        ),
    },
    {
        "legacy_feature": "ADH1B*2",
        "status": "not_sampled_by_population_module",
        "note": (
            "ADH1B is used by the aldehyde/alcohol flux model but not currently "
            "sampled in the population module."
        ),
    },
    {
        "legacy_feature": "MGMT promoter",
        "status": "not_sampled_by_population_module",
        "note": (
            "MGMT repair state is used by selected proxy terms but not sampled "
            "in the population module."
        ),
    },
    {
        "legacy_feature": "CYP3A5*3",
        "status": "sampled_but_not_active_in_current_composite_score",
        "note": (
            "CYP3A5 is sampled but current compute_full_profile() uses CYP3A4 "
            "for aflatoxin activation."
        ),
    },
    {
        "legacy_feature": "CYP2D6*4",
        "status": "sampled_but_not_active_in_current_composite_score",
        "note": (
            "CYP2D6 is sampled and has a modifier, but no current carcinogen "
            "class in compute_full_profile() uses it."
        ),
    },
]

audit_csv = OUTPUT_DIR / "figure6_legacy_variant_audit.csv"
with audit_csv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0].keys()))
    writer.writeheader()
    writer.writerows(audit_rows)

print(contribution_csv)
print(population_csv)
print(audit_csv)
            """
        ),
        markdown_cell(
            """
## Render Figure 6
            """
        ),
        code_cell(
            """
plot_rows = [
    row
    for row in contribution_rows
    if float(row["normalized_contribution_percent"]) > 0
]
plot_rows = sorted(plot_rows, key=lambda row: float(row["normalized_contribution_percent"]))

fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=False)
fig.subplots_adjust(left=0.24, right=0.97, top=0.80, bottom=0.17)

y = np.arange(len(plot_rows))
values = [float(row["normalized_contribution_percent"]) for row in plot_rows]
colors = [CATEGORY_COLORS[str(row["category"])] for row in plot_rows]
labels = [str(row["label"]) for row in plot_rows]

ax.barh(y, values, color=colors, height=0.68)

for ypos, value in zip(y, values):
    ax.text(value + 0.6, ypos, f"{value:.1f}%", va="center", ha="left", fontsize=7)

ax.set_yticks(y)
ax.set_yticklabels(labels)
ax.set_xlabel("Normalized contribution to modeled variant-driven variance (%)")
ax.set_xlim(0, max(values) * 1.18 if values else 1)
ax.grid(axis="x", color="#dddddd", linewidth=0.6)
ax.set_axisbelow(True)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

legend_handles = [
    Patch(facecolor=color, edgecolor="none", label=category)
    for category, color in CATEGORY_COLORS.items()
    if category in {str(row["category"]) for row in plot_rows}
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    frameon=True,
    framealpha=0.95,
    fontsize=6.5,
)

fig.suptitle(
    "Sensitivity Analysis: Variant Contributions to Composite ExposoGraph Score Variance",
    fontsize=9.6,
    weight="bold",
    y=0.965,
)
fig.text(
    0.5,
    0.900,
    "Feature-neutralization analysis across 50,000 ancestry-stratified synthetic genotypes",
    ha="center",
    va="bottom",
    fontsize=7.2,
)
fig.text(
    0.5,
    0.045,
    (
        "Data source: ExposoGraph population_simulation._core genotype tables and "
        "flux_engine.compute_full_profile(); 10,000 individuals per ancestry group. "
        "Bars are normalized across modeled features with positive variance reduction."
    ),
    ha="center",
    va="bottom",
    fontsize=6.2,
    color="#666666",
)

png_path = OUTPUT_DIR / "figure6_variant_contributions_sensitivity_analysis.png"
pdf_path = OUTPUT_DIR / "figure6_variant_contributions_sensitivity_analysis.pdf"
svg_path = OUTPUT_DIR / "figure6_variant_contributions_sensitivity_analysis.svg"
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
    nb = build_notebook()
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(nb, NOTEBOOK_PATH)
    print(NOTEBOOK_PATH)


if __name__ == "__main__":
    main()
