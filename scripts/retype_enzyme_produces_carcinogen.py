"""One-off curation fix: retype Enzyme--PRODUCES-->Carcinogen edges.

Context (see docs/design/refactor_ui_map_viewer_plan.md follow-up discussion):
an enzyme cannot metabolically *produce* an exogenous carcinogen out of
nothing -- it can only act on one as a substrate, or produce a downstream
*metabolite*. 19 edges in the bundled reference graph were encoded backward:
``Enzyme --PRODUCES--> Carcinogen`` instead of the correct
``Carcinogen --SUBSTRATE_OF--> Enzyme``. This script performs that retype in
place on ``ExposoGraph/map/graph-data.json``, matching the existing
SUBSTRATE_OF edge schema (source/target/type/carcinogen/evidence/origin/
match_status/provenance -- no ``label``/``canonical_predicate`` fields).

Run once: ``python scripts/retype_enzyme_produces_carcinogen.py``. Idempotent
-- if no matching PRODUCES edges remain, it's a no-op.
"""

from __future__ import annotations

import json
from pathlib import Path

GRAPH_DATA_PATH = Path(__file__).resolve().parent.parent / "ExposoGraph" / "map" / "graph-data.json"

# GSTO1 -> AntimonyTrivalent is a genuine bioactivation (Sb(V) -> Sb(III)),
# not a "carcinogen consumed as substrate" case -- flagged distinctly below.
_BIOACTIVATION_NOTE = (
    " [Retyped 2026-08-23 from an Enzyme\u2192Carcinogen PRODUCES edge during "
    "curation review, per explicit instruction to convert all such edges to "
    "SUBSTRATE_OF. CAUTION: this one does not fit the substrate-of pattern "
    "cleanly -- the evidence describes GSTO1 *producing* Sb(III) from Sb(V), "
    "not consuming AntimonyTrivalent as a substrate. The Sb(V) precursor is "
    "not modeled as a separate node in this graph. Recommend follow-up "
    "review rather than treating this retype as validated.]"
)
_STANDARD_NOTE = (
    " [Retyped 2026-08-23 from an Enzyme\u2192Carcinogen PRODUCES edge during "
    "curation review -- an enzyme cannot produce an exogenous carcinogen; "
    "the original edge encoded this enzyme metabolizing this carcinogen as "
    "substrate.]"
)


def main() -> None:
    data = json.loads(GRAPH_DATA_PATH.read_text())
    node_by_id = {n["id"]: n for n in data["nodes"]}

    converted = []
    for edge in data["edges"]:
        if edge.get("type") != "PRODUCES":
            continue
        src_type = node_by_id.get(edge["source"], {}).get("type")
        tgt_type = node_by_id.get(edge["target"], {}).get("type")
        if src_type != "Enzyme" or tgt_type != "Carcinogen":
            continue

        enzyme_id = edge["source"]
        carcinogen_id = edge["target"]
        carcinogen_label = node_by_id[carcinogen_id].get("label", carcinogen_id)

        if edge.get("evidence"):
            base_evidence = edge["evidence"]
        elif edge.get("label"):
            base_evidence = f"{carcinogen_label} is a metabolic substrate of {enzyme_id}."
        else:
            base_evidence = f"{carcinogen_label} is a metabolic substrate of {enzyme_id}."

        note = _BIOACTIVATION_NOTE if (enzyme_id, carcinogen_id) == ("GSTO1", "AntimonyTrivalent") else _STANDARD_NOTE

        edge["source"] = carcinogen_id
        edge["target"] = enzyme_id
        edge["type"] = "SUBSTRATE_OF"
        edge["evidence"] = base_evidence + note
        edge.pop("label", None)
        edge.pop("canonical_predicate", None)
        # carcinogen / origin / match_status / provenance left untouched.

        converted.append(f"{carcinogen_id} --SUBSTRATE_OF--> {enzyme_id}")

    # Match the existing file's exact formatting convention (no ensure_ascii
    # escaping, no added trailing newline) to keep the diff to only the
    # edges that actually changed.
    GRAPH_DATA_PATH.write_text(json.dumps(data, indent=2, ensure_ascii=False))

    print(f"Converted {len(converted)} edges:")
    for line in converted:
        print(" ", line)


if __name__ == "__main__":
    main()
