#!/usr/bin/env python3
"""Run the Module 3 simple individual-carcinogen workflow locally."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ExposoGraph import CarcinogenClass, compute_pathway_flux


def main() -> None:
    result = compute_pathway_flux(
        CarcinogenClass.PAH,
        {"CYP1A1": "NM", "GSTM1": "NM", "GSTT1": "NM"},
        tissue="Lung",
    )
    summary = {
        "workflow": "Module 3 simple individual-carcinogen Flux Engine workflow",
        "carcinogen_class": result.carcinogen_class,
        "tissue": result.tissue,
        "net_ratio": result.net_ratio,
        "risk_classification": result.risk_classification,
        "total_activation": result.total_activation,
        "total_detox": result.total_detox,
        "warnings": result.warnings,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
