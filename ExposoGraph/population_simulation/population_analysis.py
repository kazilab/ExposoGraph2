#!/usr/bin/env python3
"""
Population-Level Analysis for ExposoGraph Simulations
=======================================================
Aggregates individual patient results into population-level statistics,
ancestry-stratified distributions, and genotype-exposure interaction analyses.

Authors: Kenneth J. Pienta (JHU), Julhash U. Kazi (Lund University)
Date: March 2026
"""

import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, TypeAlias, cast

ResultRow: TypeAlias = dict[str, Any]
ResultRows: TypeAlias = list[ResultRow]
CancerLabel: TypeAlias = dict[str, Any]
CancerLabels: TypeAlias = dict[str, CancerLabel]
JsonDict: TypeAlias = dict[str, Any]
NumericStats: TypeAlias = dict[str, int | float]

# ══════════════════════════════════════════════════════════
#  DATA LOADING
# ══════════════════════════════════════════════════════════

def load_results(filepath: str) -> ResultRows:
    """
    Load simulation results from JSONL file.

    Parameters
    ----------
    filepath : str
        Path to results.jsonl from batch_runner

    Returns
    -------
    list : List of result dicts
    """
    results: ResultRows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line:
                payload = json.loads(line)
                if isinstance(payload, dict):
                    results.append(cast(ResultRow, payload))
    print(f"  Loaded {len(results)} results from {filepath}")
    return results


# ══════════════════════════════════════════════════════════
#  DESCRIPTIVE STATISTICS
# ══════════════════════════════════════════════════════════

def compute_descriptive_stats(results: ResultRows) -> JsonDict:
    """
    Compute population-level descriptive statistics.

    Parameters
    ----------
    results : list
        Participant-level result dicts produced by the batch runner or the
        package simulation helpers.

    Returns
    -------
    dict
        Summary metrics covering cohort size, ancestry and exposure-scenario
        distributions, interaction-factor and GSH statistics, high-risk
        pathway counts, and critical-warning counts.
    """
    n = len(results)
    if n == 0:
        return {"error": "No results to analyze"}

    stats: JsonDict = {"n_total": n}

    # Ancestry distribution
    ancestry_counts = Counter(r.get("ancestry", "unknown") for r in results)
    stats["ancestry_distribution"] = dict(ancestry_counts.most_common())

    # Exposure scenario distribution
    scenario_counts = Counter(r.get("exposure_scenario", "unknown") for r in results)
    stats["exposure_scenario_distribution"] = dict(scenario_counts.most_common())

    # Interaction factor statistics
    if_values = [
        float(r["interaction_factor"])
        for r in results
        if r.get("interaction_factor") is not None
    ]
    if if_values:
        stats["interaction_factor"] = _compute_numeric_stats(if_values)
        stats["interaction_factor"]["n_synergistic"] = sum(
            1 for v in if_values if v > 1.2)
        stats["interaction_factor"]["pct_synergistic"] = (
            stats["interaction_factor"]["n_synergistic"] / len(if_values) * 100
        )

    # GSH statistics
    gsh_values = [
        float(r["gsh_fraction"])
        for r in results
        if r.get("gsh_fraction") is not None
    ]
    if gsh_values:
        stats["gsh_fraction"] = _compute_numeric_stats(gsh_values)
        n_tipping = sum(1 for r in results if r.get("gsh_tipping_point"))
        stats["gsh_fraction"]["n_tipping_point"] = n_tipping
        stats["gsh_fraction"]["pct_tipping_point"] = n_tipping / n * 100

    # High-risk pathway counts
    hr_counts = [int(r.get("n_high_risk_pathways", 0)) for r in results]
    stats["high_risk_pathways"] = _compute_numeric_stats(hr_counts)
    stats["high_risk_pathways"]["n_with_any_high_risk"] = sum(
        1 for v in hr_counts if v > 0)
    stats["high_risk_pathways"]["pct_with_any_high_risk"] = (
        sum(1 for v in hr_counts if v > 0) / n * 100
    )

    # Most common high-risk pathways
    all_hr_pathways = []
    for r in results:
        all_hr_pathways.extend(r.get("high_risk_pathways", []))
    stats["top_high_risk_pathways"] = dict(
        Counter(all_hr_pathways).most_common(10))

    # Critical warnings
    warn_counts = [int(r.get("n_critical_warnings", 0)) for r in results]
    stats["critical_warnings"] = _compute_numeric_stats(warn_counts)
    stats["critical_warnings"]["n_with_warnings"] = sum(
        1 for v in warn_counts if v > 0)

    return stats


def _compute_numeric_stats(values: list[float] | list[int]) -> NumericStats:
    """Compute mean, median, min, max, SD for a list of numbers."""
    if not values:
        return {}
    n = len(values)
    sorted_v = sorted(values)
    mean = sum(values) / n
    median = sorted_v[n // 2] if n % 2 else (sorted_v[n//2-1] + sorted_v[n//2]) / 2
    variance = sum((v - mean) ** 2 for v in values) / n if n > 1 else 0
    sd = math.sqrt(variance)

    return {
        "n": n,
        "mean": round(mean, 4),
        "median": round(median, 4),
        "min": round(min(values), 4),
        "max": round(max(values), 4),
        "sd": round(sd, 4),
        "p25": round(sorted_v[int(n * 0.25)], 4),
        "p75": round(sorted_v[int(n * 0.75)], 4),
        "p95": round(sorted_v[int(n * 0.95)], 4),
    }


# ══════════════════════════════════════════════════════════
#  ANCESTRY-STRATIFIED ANALYSIS
# ══════════════════════════════════════════════════════════

def stratify_by_ancestry(results: ResultRows) -> dict[str, JsonDict]:
    """
    Compute statistics stratified by genetic ancestry.

    Returns
    -------
    dict : {ancestry_label: stats_dict}
    """
    groups: dict[str, ResultRows] = {}
    for r in results:
        ancestry = r.get("ancestry", "unknown")
        if ancestry not in groups:
            groups[ancestry] = []
        groups[ancestry].append(r)

    stratified: dict[str, JsonDict] = {}
    for ancestry, group_results in sorted(groups.items()):
        stratified[ancestry] = compute_descriptive_stats(group_results)

    return stratified


# ══════════════════════════════════════════════════════════
#  GENOTYPE-EXPOSURE INTERACTION ANALYSIS
# ══════════════════════════════════════════════════════════

def analyze_genotype_exposure_interactions(
    results: ResultRows,
    cancer_labels: CancerLabels | None = None
) -> JsonDict:
    """
    Analyze genotype × exposure interactions across the population.

    Computes:
    - Odds ratios for genotype-exposure combinations
    - Stratified risk scores by genotype group
    - Population-attributable risk fractions

    Parameters
    ----------
    results : list
        Simulation results from batch_runner
    cancer_labels : dict, optional
        {person_id: cancer_label_dict} from phenotype_extractor

    Returns
    -------
    dict : interaction analysis results
    """
    analysis = {}

    # ── GSTM1 × Smoking ──────────────────────────────
    gstm1_smoking = _analyze_2x2(
        results,
        gene="gstm1", risk_allele="null",
        exposure_key="exposure_scenario",
        exposure_values=["smoker", "smoker_heavy_drinker",
                         "smoker_industrial_worker",
                         "smoker_industrial_heavy_drinker"],
        outcome_metric="interaction_factor",
        cancer_labels=cancer_labels,
        cancer_type="lung_cancer",
    )
    analysis["GSTM1_null_x_smoking"] = gstm1_smoking

    # ── ALDH2 × Heavy alcohol ────────────────────────
    aldh2_alcohol = _analyze_2x2(
        results,
        gene="aldh2", risk_allele="*1/*2",
        exposure_key="exposure_scenario",
        exposure_values=["heavy_drinker", "smoker_heavy_drinker",
                         "smoker_industrial_heavy_drinker"],
        outcome_metric="interaction_factor",
        cancer_labels=cancer_labels,
        cancer_type="esophageal_cancer",
    )
    analysis["ALDH2_star2_x_alcohol"] = aldh2_alcohol

    # ── NAT2 slow × Smoking ─────────────────────────
    nat2_smoking = _analyze_2x2(
        results,
        gene="nat2", risk_allele="slow",
        exposure_key="exposure_scenario",
        exposure_values=["smoker", "smoker_heavy_drinker",
                         "smoker_industrial_worker"],
        outcome_metric="interaction_factor",
        cancer_labels=cancer_labels,
        cancer_type="bladder_cancer",
    )
    analysis["NAT2_slow_x_smoking"] = nat2_smoking

    # ── GSTT1 × Occupational (TCE) ──────────────────
    gstt1_occupational = _analyze_2x2(
        results,
        gene="gstt1", risk_allele="present",  # NOTE: GSTT1-active = RISK for TCE
        exposure_key="exposure_scenario",
        exposure_values=["industrial_worker", "smoker_industrial_worker"],
        outcome_metric="interaction_factor",
        cancer_labels=cancer_labels,
        cancer_type="kidney_cancer",
    )
    analysis["GSTT1_active_x_TCE_occupational"] = gstt1_occupational

    return analysis


def _analyze_2x2(results: ResultRows, gene: str, risk_allele: str,
                 exposure_key: str, exposure_values: list[str],
                 outcome_metric: str,
                 cancer_labels: CancerLabels | None = None,
                 cancer_type: str | None = None) -> JsonDict:
    """
    Compute a 2×2 genotype × exposure analysis.

    Groups:
      A: gene_risk + exposure_present (high risk)
      B: gene_risk + no_exposure
      C: no_gene_risk + exposure_present
      D: no_gene_risk + no_exposure (reference)
    """
    groups: dict[str, ResultRows] = {"A": [], "B": [], "C": [], "D": []}

    for r in results:
        gene_val = r.get(gene, "")
        has_risk_allele = (str(gene_val).lower() == str(risk_allele).lower())
        has_exposure = r.get(exposure_key, "") in exposure_values

        if has_risk_allele and has_exposure:
            groups["A"].append(r)
        elif has_risk_allele:
            groups["B"].append(r)
        elif has_exposure:
            groups["C"].append(r)
        else:
            groups["D"].append(r)

    analysis: JsonDict = {
        "gene": gene,
        "risk_allele": risk_allele,
        "exposure": exposure_values,
        "group_sizes": {k: len(v) for k, v in groups.items()},
    }

    def _has_cancer_label(row: ResultRow) -> bool:
        if cancer_labels is None or cancer_type is None:
            return False
        pid = row.get("person_id")
        label = cancer_labels.get(str(pid), {}) if pid is not None else {}
        return cancer_type in cast(list[str], label.get("cancer_types", []))

    # Mean outcome metric per group
    for group_name, group_results in groups.items():
        values = [
            float(r[outcome_metric])
            for r in group_results
            if r.get(outcome_metric) is not None
        ]
        if values:
            analysis[f"mean_{outcome_metric}_{group_name}"] = round(
                sum(values) / len(values), 4)

    # Odds ratio if cancer labels provided
    if cancer_labels and cancer_type:
        a_cases = sum(1 for r in groups["A"] if _has_cancer_label(r))
        b_cases = sum(1 for r in groups["B"] if _has_cancer_label(r))
        c_cases = sum(1 for r in groups["C"] if _has_cancer_label(r))
        d_cases = sum(1 for r in groups["D"] if _has_cancer_label(r))

        a_n = len(groups["A"])
        d_n = len(groups["D"])

        analysis["cancer_type"] = cancer_type
        analysis["case_counts"] = {
            "A": a_cases, "B": b_cases, "C": c_cases, "D": d_cases
        }

        # Odds ratio: (A_cases/A_controls) / (D_cases/D_controls)
        a_ctrl = a_n - a_cases
        d_ctrl = d_n - d_cases
        if a_ctrl > 0 and d_cases > 0 and d_ctrl > 0:
            or_val = (a_cases * d_ctrl) / (a_ctrl * d_cases)
            analysis["odds_ratio_AxD"] = round(or_val, 3)

            # 95% CI using Woolf's method
            if a_cases > 0:
                log_or = math.log(or_val)
                se = math.sqrt(1/a_cases + 1/a_ctrl + 1/d_cases + 1/d_ctrl)
                ci_low = math.exp(log_or - 1.96 * se)
                ci_high = math.exp(log_or + 1.96 * se)
                analysis["or_95ci"] = [round(ci_low, 3), round(ci_high, 3)]

        # Population-attributable risk fraction
        if d_n > 0:
            p_exposed_risk = a_n / (a_n + d_n) if (a_n + d_n) > 0 else 0
            or_val = float(analysis.get("odds_ratio_AxD", 1.0))
            if or_val > 0:
                par = p_exposed_risk * (or_val - 1) / (p_exposed_risk * (or_val - 1) + 1)
                analysis["population_attributable_fraction"] = round(par, 4)

    return analysis


# ══════════════════════════════════════════════════════════
#  RISK SCORE DISTRIBUTIONS
# ══════════════════════════════════════════════════════════

def compute_risk_distributions(results: ResultRows) -> JsonDict:
    """
    Compute risk score distributions for the population.

    Returns histograms and percentiles for key metrics.
    """
    distributions: JsonDict = {}

    # Interaction factor distribution
    if_values = [
        float(r["interaction_factor"])
        for r in results
        if r.get("interaction_factor") is not None
    ]
    if if_values:
        distributions["interaction_factor"] = {
            "stats": _compute_numeric_stats(if_values),
            "histogram": _histogram(if_values, bins=20),
            "risk_strata": {
                "synergistic_gt_2": sum(1 for v in if_values if v > 2.0),
                "synergistic_1_2_to_2": sum(1 for v in if_values if 1.2 < v <= 2.0),
                "additive_0_8_to_1_2": sum(1 for v in if_values if 0.8 <= v <= 1.2),
                "protective_lt_0_8": sum(1 for v in if_values if v < 0.8),
            },
        }

    # GSH depletion distribution
    gsh_values = [
        float(r["gsh_fraction"])
        for r in results
        if r.get("gsh_fraction") is not None
    ]
    if gsh_values:
        distributions["gsh_fraction"] = {
            "stats": _compute_numeric_stats(gsh_values),
            "histogram": _histogram(gsh_values, bins=20),
            "risk_strata": {
                "critical_lt_20pct": sum(1 for v in gsh_values if v < 0.20),
                "low_20_50pct": sum(1 for v in gsh_values if 0.20 <= v < 0.50),
                "moderate_50_80pct": sum(1 for v in gsh_values if 0.50 <= v < 0.80),
                "normal_gt_80pct": sum(1 for v in gsh_values if v >= 0.80),
            },
        }

    # Per-carcinogen-class flux ratio distributions
    class_ratios: dict[str, list[float]] = {}
    for r in results:
        flux_classes = r.get("flux_classes", {})
        if not isinstance(flux_classes, dict):
            continue
        for cls_name, cls_data in flux_classes.items():
            if isinstance(cls_data, dict) and cls_data.get("net_ratio") is not None:
                if cls_name not in class_ratios:
                    class_ratios[cls_name] = []
                class_ratios[cls_name].append(float(cls_data["net_ratio"]))

    distributions["per_class_flux_ratios"] = {}
    for cls_name, ratios in sorted(class_ratios.items()):
        distributions["per_class_flux_ratios"][cls_name] = {
            "stats": _compute_numeric_stats(ratios),
            "n_elevated": sum(1 for v in ratios if v > 2.0),
            "n_high": sum(1 for v in ratios if v > 5.0),
        }

    return distributions


def _histogram(values: list[float], bins: int = 20) -> JsonDict:
    """Compute a simple histogram."""
    if not values:
        return {}
    mn, mx = min(values), max(values)
    if mn == mx:
        return {"bins": [mn], "counts": [len(values)]}

    bin_width = (mx - mn) / bins
    bin_edges = [mn + i * bin_width for i in range(bins + 1)]
    counts = [0] * bins
    for v in values:
        idx = min(int((v - mn) / bin_width), bins - 1)
        counts[idx] += 1
    return {
        "bin_edges": [round(e, 4) for e in bin_edges],
        "counts": counts,
        "bin_width": round(bin_width, 4),
    }


# ══════════════════════════════════════════════════════════
#  REPORT GENERATION
# ══════════════════════════════════════════════════════════

def generate_population_report(
    results: ResultRows,
    cancer_labels: CancerLabels | None = None,
    output_dir: str = "./simulation_output"
) -> JsonDict:
    """
    Generate a comprehensive population analysis report.

    Saves to output_dir/population_analysis.json

    Returns
    -------
    dict : Complete analysis results
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    report: JsonDict = {}

    print("\n  Computing descriptive statistics...")
    report["descriptive"] = compute_descriptive_stats(results)

    print("  Computing ancestry-stratified analysis...")
    report["ancestry_stratified"] = stratify_by_ancestry(results)

    print("  Computing risk distributions...")
    report["distributions"] = compute_risk_distributions(results)

    print("  Analyzing genotype-exposure interactions...")
    report["interactions"] = analyze_genotype_exposure_interactions(
        results, cancer_labels)

    # Save
    report_file = out / "population_analysis.json"
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n  Analysis saved to: {report_file}")

    # Print summary
    _print_report_summary(report)

    return report


def _print_report_summary(report: JsonDict) -> None:
    """Print a human-readable summary of the analysis."""
    desc = cast(JsonDict, report.get("descriptive", {}))
    print("\n" + "=" * 70)
    print("  POPULATION ANALYSIS SUMMARY")
    print("=" * 70)
    print(f"\n  Total participants: {desc.get('n_total', 'N/A')}")

    # Ancestry
    anc = cast(dict[str, int], desc.get("ancestry_distribution", {}))
    if anc:
        print("\n  Ancestry distribution:")
        for a, count in anc.items():
            print(f"    {a}: {count} ({count/desc['n_total']*100:.1f}%)")

    # Interaction factors
    if_stats = cast(JsonDict, desc.get("interaction_factor", {}))
    if if_stats:
        print("\n  Interaction factor:")
        print(f"    Mean: {if_stats.get('mean', 'N/A')}")
        print(f"    Median: {if_stats.get('median', 'N/A')}")
        print(f"    Range: [{if_stats.get('min', 'N/A')}, {if_stats.get('max', 'N/A')}]")
        print(f"    Synergistic (>1.2): {if_stats.get('n_synergistic', 0)} "
              f"({if_stats.get('pct_synergistic', 0):.1f}%)")

    # GSH
    gsh_stats = cast(JsonDict, desc.get("gsh_fraction", {}))
    if gsh_stats:
        print("\n  GSH depletion:")
        print(f"    Mean GSH fraction: {gsh_stats.get('mean', 'N/A')}")
        print(f"    At tipping point: {gsh_stats.get('n_tipping_point', 0)} "
              f"({gsh_stats.get('pct_tipping_point', 0):.1f}%)")

    # High-risk pathways
    hr = cast(JsonDict, desc.get("high_risk_pathways", {}))
    if hr:
        print("\n  High-risk pathways:")
        print(f"    Participants with any: {hr.get('n_with_any_high_risk', 0)} "
              f"({hr.get('pct_with_any_high_risk', 0):.1f}%)")

    top_hr = cast(JsonDict, desc.get("top_high_risk_pathways", {}))
    if top_hr:
        print(f"    Top classes: {dict(list(top_hr.items())[:5])}")

    # Genotype-exposure interactions
    interactions = cast(dict[str, JsonDict], report.get("interactions", {}))
    if interactions:
        print("\n  Genotype × Exposure Interactions:")
        for name, inter in interactions.items():
            or_val = inter.get("odds_ratio_AxD")
            if or_val:
                ci = inter.get("or_95ci", [])
                ci_str = f" (95% CI: {ci[0]}-{ci[1]})" if ci else ""
                print(f"    {name}: OR = {or_val}{ci_str}")

    print("\n" + "=" * 70)


# ══════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Population-Level Analysis for ExposoGraph Simulations"
    )
    parser.add_argument("--results", type=str, required=True,
                        help="Path to results.jsonl from batch_runner")
    parser.add_argument("--cancer-labels", type=str,
                        help="Path to cancer_labels.json (optional)")
    parser.add_argument("--output", type=str, default="./simulation_output",
                        help="Output directory")
    args = parser.parse_args()

    results = load_results(args.results)

    cancer_labels: CancerLabels | None = None
    if args.cancer_labels:
        with open(args.cancer_labels) as f:
            loaded = json.load(f)
        if isinstance(loaded, dict):
            cancer_labels = cast(CancerLabels, loaded)
            print(f"  Loaded cancer labels for {len(cancer_labels)} participants")
        else:
            raise ValueError("Cancer labels JSON must be an object keyed by participant id")

    generate_population_report(results, cancer_labels, args.output)


if __name__ == "__main__":
    main()
