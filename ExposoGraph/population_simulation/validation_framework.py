#!/usr/bin/env python3
"""
Validation Framework for ExposoGraph Population Simulations
=============================================================
Validates ExposoGraph risk predictions against observed cancer outcomes
in the All of Us cohort, using epidemiological comparison metrics.

Compares predicted risks against:
1. Published genotype-cancer odds ratios (IARC, meta-analyses)
2. Observed cancer incidence in All of Us EHR data
3. Known exposure-cancer associations

Authors: Kenneth J. Pienta (JHU), Julhash U. Kazi (Lund University)
Date: March 2026
"""

import json
import math
from pathlib import Path
from typing import Any, TypeAlias, cast

ResultRow: TypeAlias = dict[str, Any]
ResultRows: TypeAlias = list[ResultRow]
CancerLabel: TypeAlias = dict[str, Any]
CancerLabels: TypeAlias = dict[str, CancerLabel]
JsonDict: TypeAlias = dict[str, Any]

# ══════════════════════════════════════════════════════════
#  PUBLISHED REFERENCE VALUES
# ══════════════════════════════════════════════════════════

# Known genotype-exposure-cancer odds ratios from epidemiological literature
# Used as gold-standard references for validation
REFERENCE_ORS: dict[str, JsonDict] = {
    "GSTM1_null_smoking_lung": {
        "gene": "GSTM1",
        "risk_allele": "null",
        "exposure": "smoking",
        "cancer": "lung_cancer",
        "published_OR": 1.58,
        "ci_low": 1.21,
        "ci_high": 2.06,
        "source": "Ye et al. 2006, Lung Cancer meta-analysis",
        "pmid": "16893600",
        "n_studies": 41,
    },
    "GSTM1_null_smoking_bladder": {
        "gene": "GSTM1",
        "risk_allele": "null",
        "exposure": "smoking",
        "cancer": "bladder_cancer",
        "published_OR": 1.53,
        "ci_low": 1.11,
        "ci_high": 2.12,
        "source": "Engel et al. 2002, Cancer Epidemiol Biomarkers",
        "pmid": "12163336",
        "n_studies": 15,
    },
    "CYP1A1_2A_smoking_lung": {
        "gene": "CYP1A1",
        "risk_allele": "*2A",
        "exposure": "smoking",
        "cancer": "lung_cancer",
        "published_OR": 2.36,
        "ci_low": 1.16,
        "ci_high": 4.81,
        "source": "Shi et al. 2008, Eur J Cancer meta-analysis",
        "pmid": "18313914",
        "n_studies": 12,
    },
    "GSTM1_CYP1A1_combined_lung": {
        "gene": "GSTM1+CYP1A1",
        "risk_allele": "null+*2A",
        "exposure": "smoking",
        "cancer": "lung_cancer",
        "published_OR": 2.87,
        "ci_low": 1.73,
        "ci_high": 4.76,
        "source": "Vineis et al. 2007, combined meta-analysis",
        "pmid": "17522069",
        "n_studies": 8,
    },
    "ALDH2_star2_alcohol_esophageal": {
        "gene": "ALDH2",
        "risk_allele": "*1/*2",
        "exposure": "heavy_alcohol",
        "cancer": "esophageal_cancer",
        "published_OR": 6.97,
        "ci_low": 4.36,
        "ci_high": 11.12,
        "source": "Yokoyama & Omori 2005, Int J Cancer",
        "pmid": "15386415",
        "n_studies": 10,
    },
    "ALDH2_star2_homozygous_esophageal": {
        "gene": "ALDH2",
        "risk_allele": "*2/*2",
        "exposure": "heavy_alcohol",
        "cancer": "esophageal_cancer",
        "published_OR": 12.5,
        "ci_low": 6.0,
        "ci_high": 26.0,
        "source": "Yokoyama & Omori 2005, Int J Cancer",
        "pmid": "15386415",
        "n_studies": 10,
    },
    "NAT2_slow_smoking_bladder": {
        "gene": "NAT2",
        "risk_allele": "slow",
        "exposure": "smoking",
        "cancer": "bladder_cancer",
        "published_OR": 1.51,
        "ci_low": 1.28,
        "ci_high": 1.78,
        "source": "Garcia-Closas et al. 2005, Lancet",
        "pmid": "16112300",
        "n_studies": 22,
    },
    "GSTT1_active_TCE_kidney": {
        "gene": "GSTT1",
        "risk_allele": "present",
        "exposure": "occupational_TCE",
        "cancer": "kidney_cancer",
        "published_OR": 1.88,
        "ci_low": 1.06,
        "ci_high": 3.33,
        "source": "Karami et al. 2012, Carcinogenesis",
        "pmid": "22610076",
        "n_studies": 5,
        "note": "Inverted genetics: GSTT1-active converts TCE to "
                "DCVC (nephrotoxic); GSTT1-null is protective",
    },
}


# ══════════════════════════════════════════════════════════
#  DISCRIMINATION METRICS (ROC / AUC)
# ══════════════════════════════════════════════════════════

def compute_roc_auc(
    results: ResultRows,
    cancer_labels: CancerLabels,
    score_key: str = "interaction_factor",
    cancer_type: str = "any_cancer"
) -> JsonDict:
    """
    Compute ROC curve and AUC for a risk score vs cancer outcome.

    Parameters
    ----------
    results : list
        Simulation results from batch_runner
    cancer_labels : dict
        {person_id: {cancer_types: [...]}} from phenotype_extractor
    score_key : str
        Which result metric to use as the risk score
    cancer_type : str
        Cancer type to use as outcome

    Returns
    -------
    dict : {auc, roc_points, n_cases, n_controls, thresholds}
    """
    # Build score-outcome pairs
    pairs: list[tuple[float, int]] = []
    for r in results:
        pid = r.get("person_id")
        raw_score = r.get(score_key)
        if raw_score is None or pid is None:
            continue
        score = float(raw_score)

        label = cancer_labels.get(str(pid), {})
        if cancer_type == "any_cancer":
            is_case = bool(label.get("any_cancer", False))
        else:
            is_case = cancer_type in cast(list[str], label.get("cancer_types", []))

        pairs.append((score, 1 if is_case else 0))

    if not pairs:
        return {"error": "No valid score-outcome pairs"}

    n_cases = sum(p[1] for p in pairs)
    n_controls = len(pairs) - n_cases

    if n_cases == 0 or n_controls == 0:
        return {
            "error": "Cannot compute AUC: need both cases and controls",
            "n_cases": n_cases,
            "n_controls": n_controls,
        }

    # Sort by score descending
    pairs.sort(key=lambda x: -x[0])

    # Compute ROC points
    roc_points: list[JsonDict] = []
    tp = fp = 0
    prev_score: float | None = None

    for score, outcome in pairs:
        if score != prev_score and prev_score is not None:
            tpr = tp / n_cases
            fpr = fp / n_controls
            roc_points.append({"fpr": round(fpr, 4), "tpr": round(tpr, 4),
                               "threshold": round(prev_score, 4)})
        if outcome == 1:
            tp += 1
        else:
            fp += 1
        prev_score = score

    # Final point
    final_threshold = 0.0 if prev_score is None else prev_score
    roc_points.append({"fpr": round(fp/n_controls, 4),
                       "tpr": round(tp/n_cases, 4),
                       "threshold": round(final_threshold, 4)})

    # Compute AUC via trapezoidal rule
    auc = 0.0
    for i in range(1, len(roc_points)):
        dx = roc_points[i]["fpr"] - roc_points[i-1]["fpr"]
        dy_avg = (roc_points[i]["tpr"] + roc_points[i-1]["tpr"]) / 2
        auc += dx * dy_avg

    # Standard error (Hanley-McNeil)
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    se = math.sqrt(
        (auc * (1-auc) + (n_cases-1)*(q1-auc**2) + (n_controls-1)*(q2-auc**2))
        / (n_cases * n_controls)
    ) if n_cases > 1 and n_controls > 1 else 0

    return {
        "auc": round(auc, 4),
        "auc_se": round(se, 4),
        "auc_95ci": [round(max(0, auc - 1.96*se), 4),
                     round(min(1, auc + 1.96*se), 4)],
        "n_cases": n_cases,
        "n_controls": n_controls,
        "n_total": len(pairs),
        "score_key": score_key,
        "cancer_type": cancer_type,
        "roc_points": roc_points[::max(1, len(roc_points)//50)],  # subsample
    }


# ══════════════════════════════════════════════════════════
#  CALIBRATION
# ══════════════════════════════════════════════════════════

def compute_calibration(
    results: ResultRows,
    cancer_labels: CancerLabels,
    score_key: str = "interaction_factor",
    cancer_type: str = "any_cancer",
    n_bins: int = 10
) -> JsonDict:
    """
    Compute calibration: predicted risk vs observed frequency.

    Parameters
    ----------
    results, cancer_labels : as above
    score_key : str
    cancer_type : str
    n_bins : int

    Returns
    -------
    dict : {bins: [{predicted_mean, observed_rate, n}], hosmer_lemeshow}
    """
    pairs: list[tuple[float, int]] = []
    for r in results:
        pid = r.get("person_id")
        raw_score = r.get(score_key)
        if raw_score is None or pid is None:
            continue
        score = float(raw_score)

        label = cancer_labels.get(str(pid), {})
        if cancer_type == "any_cancer":
            is_case = bool(label.get("any_cancer", False))
        else:
            is_case = cancer_type in cast(list[str], label.get("cancer_types", []))

        pairs.append((score, 1 if is_case else 0))

    if not pairs:
        return {"error": "No valid pairs"}

    pairs.sort(key=lambda x: x[0])
    bin_size = max(1, len(pairs) // n_bins)

    bins: list[JsonDict] = []
    for i in range(0, len(pairs), bin_size):
        batch = pairs[i:i+bin_size]
        scores = [p[0] for p in batch]
        outcomes = [p[1] for p in batch]
        bins.append({
            "predicted_mean": round(sum(scores) / len(scores), 4),
            "observed_rate": round(sum(outcomes) / len(outcomes), 4),
            "n": len(batch),
            "n_cases": sum(outcomes),
        })

    # Hosmer-Lemeshow-like statistic
    hl = 0.0
    for b in bins:
        predicted_mean = float(b["predicted_mean"])
        count = int(b["n"])
        n_cases = int(b["n_cases"])
        if predicted_mean > 0 and count > 0:
            expected = predicted_mean * count
            if expected > 0:
                hl += (n_cases - expected) ** 2 / expected

    return {
        "bins": bins,
        "hosmer_lemeshow_chi2": round(hl, 4),
        "n_bins": len(bins),
        "score_key": score_key,
        "cancer_type": cancer_type,
    }


# ══════════════════════════════════════════════════════════
#  COMPARISON AGAINST PUBLISHED ORs
# ══════════════════════════════════════════════════════════

def validate_against_published_ors(
    results: ResultRows,
    cancer_labels: CancerLabels
) -> JsonDict:
    """
    Compare observed genotype-exposure-cancer ORs from the simulation
    against published epidemiological references.

    Parameters
    ----------
    results : list
        Simulation results
    cancer_labels : dict
        Cancer outcome labels

    Returns
    -------
    dict : {reference_name: {published_OR, observed_OR, concordant, ...}}
    """
    validations: JsonDict = {}

    for ref_name, ref in REFERENCE_ORS.items():
        gene = str(ref["gene"])
        risk_allele = str(ref["risk_allele"])
        cancer_type = str(ref["cancer"])

        # Handle combined genotypes
        if "+" in gene:
            # Skip combined genotype validation for now
            validations[ref_name] = {
                "published_OR": ref["published_OR"],
                "published_CI": [ref["ci_low"], ref["ci_high"]],
                "source": ref["source"],
                "observed_OR": None,
                "note": "Combined genotype validation requires custom logic",
            }
            continue

        # Map gene name to result field
        gene_field: str | None = gene.lower()
        if gene_field == "gstm1":
            gene_field = "gstm1"
        elif gene_field == "gstt1":
            gene_field = "gstt1"
        elif gene_field == "aldh2":
            gene_field = "aldh2"
        elif gene_field == "nat2":
            gene_field = "nat2"
        elif gene_field == "cyp1a1":
            gene_field = "cyp1a1" if "cyp1a1" in (results[0] if results else {}) else None

        if gene_field is None:
            validations[ref_name] = {
                "published_OR": ref["published_OR"],
                "source": ref["source"],
                "observed_OR": None,
                "note": f"Gene {gene} not available in simulation results",
            }
            continue

        # Build 2x2 table
        a = d = 0  # cases/controls × risk/no-risk
        a_case = d_case = 0

        for r in results:
            pid = r.get("person_id")
            gene_val = str(r.get(gene_field, "")).lower()
            has_risk = gene_val == risk_allele.lower()

            # Check cancer outcome
            label = cancer_labels.get(str(pid), {})
            is_case = cancer_type in cast(list[str], label.get("cancer_types", []))

            if has_risk:
                if is_case:
                    a_case += 1
                a += 1
            else:
                if is_case:
                    d_case += 1
                d += 1

        a_ctrl = a - a_case
        d_ctrl = d - d_case

        observed_or = None
        concordant = None
        if a_ctrl > 0 and d_case > 0 and d_ctrl > 0 and a_case > 0:
            observed_or = round((a_case * d_ctrl) / (a_ctrl * d_case), 3)

            # Check if observed OR falls within published CI
            ci_low = float(ref["ci_low"])
            ci_high = float(ref["ci_high"])
            concordant = ci_low <= observed_or <= ci_high

        validations[ref_name] = {
            "gene": gene,
            "risk_allele": risk_allele,
            "cancer_type": cancer_type,
            "published_OR": ref["published_OR"],
            "published_CI": [ref["ci_low"], ref["ci_high"]],
            "source": ref["source"],
            "pmid": ref.get("pmid"),
            "observed_OR": observed_or,
            "concordant": concordant,
            "n_risk_allele": a,
            "n_cases_risk": a_case,
            "n_reference": d,
            "n_cases_reference": d_case,
        }

    # Summary
    n_tested = sum(1 for v in validations.values()
                   if v.get("observed_OR") is not None)
    n_concordant = sum(1 for v in validations.values()
                       if v.get("concordant") is True)
    validations["_summary"] = {
        "n_references": len(REFERENCE_ORS),
        "n_tested": n_tested,
        "n_concordant": n_concordant,
        "concordance_rate": round(n_concordant / n_tested, 3) if n_tested > 0 else None,
    }

    return validations


# ══════════════════════════════════════════════════════════
#  FULL VALIDATION PIPELINE
# ══════════════════════════════════════════════════════════

def run_full_validation(
    results: ResultRows,
    cancer_labels: CancerLabels,
    output_dir: str = "./simulation_output"
) -> JsonDict:
    """
    Run the complete validation pipeline.

    Parameters
    ----------
    results : list
        From batch_runner
    cancer_labels : dict
        From phenotype_extractor
    output_dir : str

    Returns
    -------
    dict
        Complete validation bundle including ROC/AUC summaries, calibration
        tables, and concordance checks against published epidemiology
        references.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    validation: JsonDict = {}

    print("\n" + "=" * 70)
    print("  ExposoGraph Validation Framework")
    print("=" * 70)

    # 1. ROC/AUC for different score types
    print("\n  Computing ROC/AUC...")
    roc_results: JsonDict = {}
    for score_key in ["interaction_factor", "n_high_risk_pathways",
                      "n_critical_warnings"]:
        for cancer in ["any_cancer", "lung_cancer", "liver_cancer",
                       "bladder_cancer", "esophageal_cancer"]:
            key = f"{score_key}_vs_{cancer}"
            roc = compute_roc_auc(results, cancer_labels, score_key, cancer)
            if "error" not in roc:
                roc_results[key] = roc
                print(f"    {key}: AUC = {roc['auc']} "
                      f"({roc['n_cases']} cases / {roc['n_controls']} controls)")
    validation["roc_auc"] = roc_results

    # 2. Calibration
    print("\n  Computing calibration...")
    cal_results: JsonDict = {}
    for score_key in ["interaction_factor", "n_high_risk_pathways"]:
        cal = compute_calibration(results, cancer_labels, score_key, "any_cancer")
        if "error" not in cal:
            cal_results[score_key] = cal
    validation["calibration"] = cal_results

    # 3. Comparison against published ORs
    print("\n  Validating against published epidemiological ORs...")
    or_validation = validate_against_published_ors(results, cancer_labels)
    validation["published_or_comparison"] = or_validation

    summary = or_validation.get("_summary", {})
    print(f"    Tested: {summary.get('n_tested', 0)}/{summary.get('n_references', 0)} references")
    print(f"    Concordant: {summary.get('n_concordant', 0)}")
    if summary.get("concordance_rate") is not None:
        print(f"    Concordance rate: {summary['concordance_rate']*100:.1f}%")

    for ref_name, ref_result in or_validation.items():
        if ref_name.startswith("_"):
            continue
        if ref_result.get("observed_OR") is not None:
            status = "PASS" if ref_result.get("concordant") else "CHECK"
            print(f"    [{status}] {ref_name}: published OR={ref_result['published_OR']}, "
                  f"observed OR={ref_result['observed_OR']}")

    # Save
    val_file = out / "validation_results.json"
    with open(val_file, "w") as f:
        json.dump(validation, f, indent=2, default=str)
    print(f"\n  Validation results saved to: {val_file}")

    return validation


# ══════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Validation Framework for ExposoGraph Simulations"
    )
    parser.add_argument("--results", type=str, required=True,
                        help="Path to results.jsonl")
    parser.add_argument("--cancer-labels", type=str, required=True,
                        help="Path to cancer_labels.json")
    parser.add_argument("--output", type=str, default="./simulation_output",
                        help="Output directory")
    parser.add_argument("--list-references", action="store_true",
                        help="List all published reference ORs")
    args = parser.parse_args()

    if args.list_references:
        print("\nPublished Reference Odds Ratios:")
        print("=" * 70)
        for name, ref in REFERENCE_ORS.items():
            print(f"\n  {name}")
            print(f"    Gene: {ref['gene']} ({ref['risk_allele']})")
            print(f"    Cancer: {ref['cancer']}")
            print(f"    OR: {ref['published_OR']} "
                  f"(95% CI: {ref['ci_low']}-{ref['ci_high']})")
            print(f"    Source: {ref['source']}")
            if ref.get("note"):
                print(f"    Note: {ref['note']}")
        return

    from .population_analysis import load_results
    results = load_results(args.results)

    with open(args.cancer_labels) as f:
        loaded = json.load(f)
    if not isinstance(loaded, dict):
        raise ValueError("Cancer labels JSON must be an object keyed by participant id")
    cancer_labels = cast(CancerLabels, loaded)

    run_full_validation(results, cancer_labels, args.output)


if __name__ == "__main__":
    main()
