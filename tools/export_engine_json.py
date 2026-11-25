#!/usr/bin/env python3
"""Write the loaded GraphEngine (reference graph plus baked blocks) as JSON."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ExposoGraph.engine import GraphEngine  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Path of the JSON file to write.")
    args = parser.parse_args()

    engine = GraphEngine()
    engine.load_reference_graph()
    args.output.write_text(engine.to_json(), encoding="utf-8")


if __name__ == "__main__":
    main()
