"""One-off script: add NodeType.SUBSTRATE nodes to map/graph-data.json.

This script is intentionally not part of the package -- it is a scoped,
auditable, single-purpose migration script for the substrate-node-creation
commit. It creates *identity* nodes only (id/label/type/detail). No
Km/Vmax/product/reaction_role/notes data from interaction_parameters.json is
copied in: that data continues to be sourced from interaction_parameters.json
at knowledge-graph instantiation time, per project convention (tissue
expression weights are sourced the same way, from tissue_expression_data.json,
via GraphEngine._apply_tissue_expression).
"""

from __future__ import annotations

import json
from pathlib import Path

GRAPH_DATA_PATH = Path(__file__).resolve().parent.parent / "ExposoGraph" / "map" / "graph-data.json"

# id -> (label, detail). `id` matches the interaction_parameters.json
# substrate key verbatim so a future graph-backed parameter provider can join
# on it without an alias layer. `detail` is general chemical/pharmacological
# identity background (independent of, and not duplicating, the specific
# enzyme-kinetics notes stored in interaction_parameters.json).
NEW_SUBSTRATE_NODES: dict[str, tuple[str, str]] = {
    "1_nitropyrene": (
        "1-Nitropyrene",
        (
            "Nitro-substituted polycyclic aromatic hydrocarbon found in diesel exhaust and "
            "combustion particulate matter."
        ),
    ),
    "3_methylindole": (
        "3-Methylindole",
        (
            "Indole derivative (skatole) formed from tryptophan degradation in the gut; a known "
            "pulmonary toxicant."
        ),
    ),
    "6_aminochrysene": (
        "6-Aminochrysene",
        "Aromatic amine derivative of the polycyclic aromatic hydrocarbon chrysene.",
    ),
    "7_ethoxyresorufin": (
        "7-Ethoxyresorufin",
        (
            "Synthetic resorufin ether widely used as a laboratory probe substrate for CYP1 "
            "O-dealkylase activity (EROD assay)."
        ),
    ),
    "AalphaC": (
        "AαC",
        (
            "Heterocyclic amine (2-amino-9H-pyrido[2,3-b]indole) formed during high-temperature "
            "cooking of protein-rich foods."
        ),
    ),
    "CDNB": (
        "CDNB",
        (
            "1-Chloro-2,4-dinitrobenzene; a synthetic electrophile used as a generic laboratory "
            "substrate for glutathione S-transferase activity assays."
        ),
    ),
    "IQ": (
        "IQ",
        (
            "Heterocyclic amine (2-amino-3-methylimidazo[4,5-f]quinoline) formed in cooked meat "
            "and fish."
        ),
    ),
    "NNAL": (
        "NNAL",
        (
            "4-(Methylnitrosamino)-1-(3-pyridyl)-1-butanol; a major tobacco-specific nitrosamine "
            "metabolite of NNK."
        ),
    ),
    "S_mephenytoin": (
        "S-Mephenytoin",
        "Anticonvulsant drug enantiomer used clinically to phenotype CYP2C19 activity.",
    ),
    "S_warfarin": (
        "S-Warfarin",
        "Anticoagulant drug enantiomer whose clearance is highly CYP2C9-genotype dependent.",
    ),
    "acetaminophen": (
        "Acetaminophen",
        (
            "Widely used analgesic/antipyretic drug; forms the reactive metabolite NAPQI at high "
            "doses."
        ),
    ),
    "benzo_a_anthracene": (
        "Benzo[a]anthracene",
        "Four-ring polycyclic aromatic hydrocarbon and IARC Group 2B combustion product.",
    ),
    "bilirubin": (
        "Bilirubin",
        "Endogenous heme breakdown product cleared primarily by hepatic glucuronidation.",
    ),
    "bufuralol": (
        "Bufuralol",
        (
            "Beta-adrenergic blocking drug used as an alternative clinical probe for CYP2D6 "
            "activity."
        ),
    ),
    "bupropion": (
        "Bupropion",
        (
            "Antidepressant/smoking-cessation drug used clinically as a selective CYP2B6 probe "
            "substrate."
        ),
    ),
    "caffeine": (
        "Caffeine",
        (
            "Widely consumed methylxanthine stimulant; a standard clinical phenotyping probe for "
            "CYP1A2 activity."
        ),
    ),
    "chloroform": (
        "Chloroform",
        (
            "Halogenated solvent and water-disinfection byproduct; oxidatively dehalogenated to "
            "reactive intermediates."
        ),
    ),
    "chrysene": (
        "Chrysene",
        (
            "Four-ring polycyclic aromatic hydrocarbon component of tobacco smoke and combustion "
            "emissions."
        ),
    ),
    "cotinine": (
        "Cotinine",
        (
            "Primary proximate metabolite of nicotine; the standard biomarker of tobacco/nicotine "
            "exposure."
        ),
    ),
    "coumarin": (
        "Coumarin",
        (
            "Naturally occurring benzopyrone compound; the canonical clinical phenotyping probe "
            "for CYP2A6 activity."
        ),
    ),
    "cyclosporine": (
        "Cyclosporine",
        ("Calcineurin-inhibitor immunosuppressant drug with saturable CYP3A4-mediated clearance."),
    ),
    "debrisoquine": (
        "Debrisoquine",
        (
            "Antihypertensive drug whose metabolism defined the original CYP2D6 (debrisoquine "
            "hydroxylase) polymorphism."
        ),
    ),
    "dextromethorphan": (
        "Dextromethorphan",
        ("Over-the-counter antitussive drug used as a standard clinical CYP2D6 phenotyping probe."),
    ),
    "dibenz_ah_anthracene": (
        "Dibenz[a,h]anthracene",
        "Five-ring polycyclic aromatic hydrocarbon and IARC Group 2A carcinogen.",
    ),
    "diclofenac": (
        "Diclofenac",
        (
            "Nonsteroidal anti-inflammatory drug metabolized by CYP2C9 to a reactive "
            "acyl-glucuronide."
        ),
    ),
    "erythromycin": (
        "Erythromycin",
        (
            "Macrolide antibiotic used clinically as a CYP3A4 activity probe (erythromycin breath "
            "test)."
        ),
    ),
    "estradiol_2_OH": (
        "Estradiol (2-hydroxylation substrate)",
        (
            "Represents the 2-hydroxylation branch of estradiol metabolism; a distinct "
            "interaction-parameters substrate key from the parent estrogen node E2."
        ),
    ),
    "estradiol_4_OH": (
        "Estradiol (4-hydroxylation substrate)",
        (
            "Represents the 4-hydroxylation branch of estradiol metabolism, associated with "
            "genotoxic catechol-estrogen formation; a distinct interaction-parameters substrate "
            "key from the parent estrogen node E2."
        ),
    ),
    "ethacrynic_acid": (
        "Ethacrynic acid",
        ("Loop diuretic drug that also acts as a glutathione S-transferase substrate/inhibitor."),
    ),
    "isoniazid": (
        "Isoniazid",
        (
            "First-line antituberculosis drug and the canonical clinical probe for NAT2 "
            "acetylator phenotyping."
        ),
    ),
    "methoxsalen": (
        "Methoxsalen",
        ("Furocoumarin (psoralen) phototherapy agent and potent mechanism-based CYP2A inhibitor."),
    ),
    "methyl_chloride": (
        "Methyl chloride",
        (
            "Small-molecule alkyl halide used industrially as a refrigerant and chemical "
            "intermediate."
        ),
    ),
    "methylene_chloride": (
        "Methylene chloride",
        (
            "Chlorinated solvent (dichloromethane) bioactivated via glutathione conjugation to a "
            "genotoxic intermediate."
        ),
    ),
    "midazolam": (
        "Midazolam",
        "Benzodiazepine sedative widely used as the reference clinical CYP3A4 activity probe.",
    ),
    "naphthalene": (
        "Naphthalene",
        (
            "Bicyclic aromatic hydrocarbon and common combustion/mothball-derived pulmonary "
            "toxicant."
        ),
    ),
    "nicotine": (
        "Nicotine",
        (
            "Primary psychoactive alkaloid of tobacco, metabolized through several minor "
            "CYP-mediated pathways."
        ),
    ),
    "nifedipine": (
        "Nifedipine",
        (
            "Dihydropyridine calcium-channel blocker with well-characterized CYP3A4-mediated "
            "oxidation."
        ),
    ),
    "omeprazole": (
        "Omeprazole",
        "Proton-pump inhibitor drug with clinically significant CYP2C19-dependent metabolism.",
    ),
    "p_aminobenzoic_acid": (
        "p-Aminobenzoic acid",
        (
            "Aromatic amine compound used as the canonical clinical probe substrate for NAT1 "
            "acetylation activity."
        ),
    ),
    "p_nitrophenol": (
        "p-Nitrophenol",
        (
            "Synthetic nitrophenol compound used as the reference laboratory probe for CYP2E1 "
            "hydroxylase activity."
        ),
    ),
    "phenacetin": (
        "Phenacetin",
        (
            "Historic analgesic drug (withdrawn) used as the classic laboratory CYP1A2 "
            "O-deethylation probe."
        ),
    ),
    "resveratrol": (
        "Resveratrol",
        (
            "Dietary stilbene polyphenol found in grapes and red wine with CYP1B1-inhibitory "
            "activity."
        ),
    ),
    "sterigmatocystin": (
        "Sterigmatocystin",
        "Mycotoxin structurally related to aflatoxin B1 and a precursor in its biosynthesis.",
    ),
    "styrene": (
        "Styrene",
        (
            "Aromatic vinyl monomer used in plastics and rubber manufacturing, bioactivated to a "
            "genotoxic epoxide."
        ),
    ),
    "styrene_oxide": (
        "Styrene oxide",
        (
            "Reactive epoxide metabolite of styrene, detoxified primarily via glutathione "
            "conjugation."
        ),
    ),
    "sulfamethazine": (
        "Sulfamethazine",
        (
            "Sulfonamide antibiotic historically used to define the NAT2 acetylator phenotype "
            "ratio test."
        ),
    ),
    "theophylline": (
        "Theophylline",
        "Methylxanthine bronchodilator drug metabolized in part by CYP1A2.",
    ),
    "tolbutamide": (
        "Tolbutamide",
        (
            "First-generation sulfonylurea antidiabetic drug and the canonical clinical CYP2C9 "
            "probe substrate."
        ),
    ),
    "trans_stilbene_oxide": (
        "trans-Stilbene oxide",
        (
            "Synthetic epoxide compound used as a laboratory probe substrate preferentially "
            "conjugated by GSTM1."
        ),
    ),
    "trichloroethylene": (
        "Trichloroethylene",
        ("Chlorinated industrial solvent oxidized primarily via CYP2E1 to reactive intermediates."),
    ),
}


def build_node(node_id: str, label: str, detail: str) -> dict:
    return {
        "id": node_id,
        "label": label,
        "type": "Substrate",
        "detail": detail,
        "group": None,
        "iarc": None,
        "phase": None,
        "role": None,
        "reactivity": None,
        "source_db": None,
        "evidence": None,
        "pmid": None,
        "tissue": None,
        "exposure": None,
        "variant": None,
        "phenotype": None,
        "activity_score": None,
        "tissue_weights": None,
        "exposure_scenarios": None,
        "tier": None,
        "origin": "imported",
        "match_status": "unmatched",
        "canonical_id": None,
        "canonical_label": None,
        "canonical_namespace": None,
        "custom_type": None,
        "provenance": [],
        "curation": None,
    }


def main() -> None:
    with GRAPH_DATA_PATH.open(encoding="utf-8") as f:
        graph_data = json.load(f)

    existing_ids = {n["id"] for n in graph_data["nodes"]}
    added = []
    for node_id, (label, detail) in NEW_SUBSTRATE_NODES.items():
        if node_id in existing_ids:
            raise ValueError(f"Node id already exists, refusing to duplicate: {node_id}")
        graph_data["nodes"].append(build_node(node_id, label, detail))
        added.append(node_id)

    with GRAPH_DATA_PATH.open("w", encoding="utf-8") as f:
        json.dump(graph_data, f, indent=2, ensure_ascii=False)
        f.write("\n")

    print(f"Added {len(added)} Substrate nodes.")
    print(f"Total nodes now: {len(graph_data['nodes'])}")


if __name__ == "__main__":
    main()
