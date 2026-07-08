#!/usr/bin/env python3
"""Run the Module 5 advanced multi-carcinogen interaction workflow locally."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ExposoGraph import compute_interaction_matrix


def main() -> None:
    result = compute_interaction_matrix(
        {"PAH": 3.0, "NNK": 4.0, "benzene": 5.0, "ethanol": 3.0},
        genotypes={"GSTM1": "null", "CYP2E1": "NM"},
        lifestyle={"smoking": True, "alcohol_moderate": True},
        tissue="Liver",
    )
    mechanism_attribution = result.mechanism_attribution or {}
    mechanism_risks = {
        name: {
            "baseline_relative_risk": risk.baseline_relative_risk,
            "adjusted_relative_risk": risk.adjusted_relative_risk,
            "inhibition_status": risk.inhibition_status,
            "review_required": risk.review_required,
            "warnings": risk.warnings,
        }
        for name, risk in result.mechanism_resolved_risks.items()
    }
    summary = {
        "workflow": "Module 5 advanced multi-carcinogen interaction workflow",
        "tissue": result.tissue,
        "interaction_factor": result.interaction_factor,
        "total_independent_risk": result.total_independent_risk,
        "total_interaction_risk": result.total_interaction_risk,
        "gsh_fraction_normal": result.gsh_status.fraction_normal,
        "synergy_pairs": result.synergy_matrix,
        "mechanism_decomposition_basis": mechanism_attribution.get(
            "decomposition_basis"
        ),
        "mechanism_resolved_risks": mechanism_risks,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
