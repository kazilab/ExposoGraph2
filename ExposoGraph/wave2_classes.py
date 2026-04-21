"""Wave 2 carcinogen class expansion: Aldehydes, Dioxins/AhR, Dietary
N-Nitroso compounds, and Chlorinated Solvents.

This module packages the full Wave 2 scientific payload (34 curated source
nodes, 43 edges) and adds one package-native helper node for ``Ethanol`` so the
aldehyde branch resolves cleanly inside ExposoGraph without relying on a
dangling cross-reference.

Classes included:
- Class 11: Aldehydes (Formaldehyde, Acetaldehyde) — IARC Group 1
- Class 12: Dioxins / AhR Ligands & Organochlorines (TCDD, PCB-126/77/118/138/153/169, 2,3,4,7,8-PeCDF, HCB, Lindane, DDT, DDE, PCP, Chlordane, Heptachlor, Toxaphene) — IARC Group 1 / 2A / 2B
- Class 13: Dietary N-Nitroso (NDMA, NDEA) — IARC Group 2A
- Class 14: Chlorinated Solvents (TCE, PCE) — IARC Group 1/2A
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .models import Edge, EdgeType, KnowledgeGraph, Node, NodeType


@dataclass
class CarcinogenClassProfile:
    """Profile for a Wave 2 carcinogen class."""

    class_id: int
    class_name: str
    display_name: str
    carcinogens: list[str]
    key_enzymes: list[str]
    cross_referenced_enzymes: list[str]
    target_organs: list[str]
    pmid_references: list[str] = field(default_factory=list)


# ── Class Profiles ────────────────────────────────────────────────────────

WAVE2_CLASS_PROFILES: dict[str, CarcinogenClassProfile] = {
    "Aldehydes": CarcinogenClassProfile(
        class_id=11,
        class_name="Aldehydes",
        display_name="Aldehydes & Alcohol-Related",
        carcinogens=[
            "Formaldehyde",
            "Acetaldehyde",
            "Urethane",
            "Acrolein",
            "Crotonaldehyde",
            "Furfural",
            "MDA",
            "4_HNE",
        ],
        key_enzymes=["ALDH2", "ADH1B", "ADH1C", "ADH5", "FANCD2"],
        cross_referenced_enzymes=[
            "CYP2E1",
            "GSTP1",
            "GSTM1",
            "GSTT1",
            "SULT1A1",
            "XRCC1",
        ],
        target_organs=[
            "nasopharynx",
            "esophagus",
            "liver",
            "lung",
            "bladder",
            "colon",
        ],
        pmid_references=[
            "PMC7941978",
            "PMC4376259",
            "PMC7758861",
            "PMC5324749",
            "12860588",   # Hoffler et al. 2003 (urethane CYP2E1 activation)
            "19705912",   # Minko et al. 2009 (acrolein-dG adducts)
            "10064842",   # Marnett 1999 (malondialdehyde / M1-dG)
            "1937131",    # Esterbauer et al. 1991 (4-HNE biochemistry)
            "21538843",   # Monien et al. 2011 (HMF/SMF sulfate-ester activation)
            "IARC_96",    # IARC Monograph 96 (urethane)
            "IARC_128",   # IARC Monograph 128 (acrolein)
            "IARC_63",    # IARC Monograph 63 (crotonaldehyde)
        ],
    ),
    "Dioxins_AhR": CarcinogenClassProfile(
        class_id=12,
        class_name="Dioxins_AhR",
        display_name="Dioxins / AhR Ligands & Organochlorines",
        carcinogens=[
            "TCDD",
            "PCB_126",
            "PCB_169",
            "PCB_77",
            "PCB_118",
            "PCB_153",
            "PCB_138",
            "PeCDF_23478",
            "HCB",
            "Lindane",
            "DDT",
            "DDE",
            "PCP",
            "Chlordane",
            "Heptachlor",
            "Toxaphene",
        ],
        key_enzymes=["AHR", "ARNT", "AHRR"],
        cross_referenced_enzymes=["CYP1A1", "CYP1B1"],
        target_organs=["liver", "lung", "soft tissue", "breast", "thyroid", "lymphatic"],
        pmid_references=[
            "PMC3748760",
            "srep34989",
            "IARC_Monograph_107",
            "IARC_Monograph_100F",
            "IARC_Monograph_79",
            "IARC_Monograph_113",
            "IARC_Monograph_117",
            "IARC_Monograph_53",
            "vandenBerg_2006_WHO_TEF",
        ],
    ),
    "Dietary_NNitroso": CarcinogenClassProfile(
        class_id=13,
        class_name="Dietary_NNitroso",
        display_name="Dietary N-Nitroso Compounds",
        carcinogens=["NDMA", "NDEA"],
        key_enzymes=[],
        cross_referenced_enzymes=["CYP2E1", "CYP2A6", "MGMT"],
        target_organs=["liver", "esophagus", "stomach"],
        pmid_references=["40390554"],
    ),
    "Chlorinated_Solvents": CarcinogenClassProfile(
        class_id=14,
        class_name="Chlorinated_Solvents",
        display_name="Chlorinated Solvents",
        carcinogens=["TCE", "PCE"],
        key_enzymes=["CCBL1"],
        cross_referenced_enzymes=["CYP2E1", "GSTT1"],
        target_organs=["kidney", "liver"],
        pmid_references=["20663906", "PMC3867557"],
    ),
    "Alkylating_Agents": CarcinogenClassProfile(
        class_id=15,
        class_name="Alkylating_Agents",
        display_name="Alkylating Carcinogens",
        carcinogens=[
            "Acrylamide",
            "Glycidamide",
            "Cyclophosphamide",
            "Chlorambucil",
            "Sulfur_mustard",
            "Busulfan",
            "MNU",
            "Temozolomide",
        ],
        key_enzymes=[],
        cross_referenced_enzymes=[
            "CYP2E1",
            "CYP2B6",
            "CYP3A4",
            "GSTP1",
            "GSTM1",
            "MGMT",
            "XRCC1",
        ],
        target_organs=["bone marrow", "bladder", "lung", "brain", "liver"],
        pmid_references=[
            "17872912",   # Besaratinia & Pfeifer, Carcinogenesis 2007 (acrylamide/glycidamide)
            "7614537",    # Segerbaeck et al., Carcinogenesis 1995 (glycidamide adducts)
            "10220571",   # Roy et al., Drug Metab Dispos 1999 (cyclophosphamide CYP2B6)
            "7523912",    # Povirk & Shuker, Mutat Res 1994 (nitrogen mustards)
            "8635461",    # Matijasevic et al., Carcinogenesis 1996 (sulfur mustard)
            "12960109",   # Hassan & Ljungman, Mutat Res 1997 (busulfan crosslinks)
            "2185966",    # Kyrtopoulos, Mutat Res 1990 (MNU O6-methylguanine)
            "9327140",    # Newlands et al., Cancer Treat Rev 1997 (temozolomide MTIC)
            "15758009",   # Hegi et al., NEJM 2005 (MGMT/temozolomide)
            "IARC_100A",  # IARC Monograph 100A (cyclophosphamide, chlorambucil, busulfan)
            "IARC_100F",  # IARC Monograph 100F (sulfur mustard)
        ],
    ),
}


# ── Node Data ─────────────────────────────────────────────────────────────

# Cross-referenced node IDs: nodes expected to exist in the core graph.
# The merge function skips these when they already exist.
WAVE2_CROSS_REFERENCES: dict[str, list[str]] = {
    "Aldehydes": [
        "CYP2E1",
        "GSTP1",
        "GSTM1",
        "GSTT1",
        "SULT1A1",
        "XRCC1",
    ],
    "Dioxins_AhR": ["CYP1A1", "CYP1B1"],
    "Dietary_NNitroso": ["CYP2E1", "CYP2A6", "MGMT"],
    "Chlorinated_Solvents": ["CYP2E1", "GSTT1"],
    "Alkylating_Agents": [
        "CYP2E1",
        "CYP2B6",
        "CYP3A4",
        "GSTP1",
        "GSTM1",
        "MGMT",
        "XRCC1",
    ],
}


WAVE2_CARCINOGEN_NODES: list[dict[str, Any]] = [
    # Class 11: Aldehydes
    {
        "id": "Formaldehyde",
        "label": "Formaldehyde",
        "type": "Carcinogen",
        "group": "Aldehyde",
        "iarc": "Group 1",
        "detail": "Building materials, industrial use, vehicle exhaust; also produced endogenously via one-carbon metabolism. IARC Group 1 carcinogen for nasopharyngeal cancer.",
        "exposure": "inhalation (occupational, indoor air), endogenous metabolism",
        "class_name": "Aldehydes",
    },
    {
        "id": "Acetaldehyde",
        "label": "Acetaldehyde",
        "type": "Carcinogen",
        "group": "Aldehyde",
        "iarc": "Group 1",
        "detail": "Primary metabolite of ethanol (ADH1B oxidation); also found in tobacco smoke and fruit ripening. IARC Group 1 carcinogen associated with esophageal squamous cell carcinoma.",
        "exposure": "ethanol metabolism, tobacco smoke, fermented foods",
        "class_name": "Aldehydes",
    },
    {
        "id": "Urethane",
        "label": "Urethane (ethyl carbamate)",
        "type": "Carcinogen",
        "group": "Alcohol",
        "iarc": "Group 2A",
        "detail": (
            "Ethyl carbamate; natural contaminant of fermented foods and "
            "alcoholic beverages, especially stone-fruit brandies, sake, "
            "sherry, and bread crust. CYP2E1 dehydrogenates urethane to "
            "vinyl carbamate, which is then re-epoxidized to vinyl "
            "carbamate epoxide — the ultimate mutagen that forms etheno-A/"
            "etheno-C DNA adducts (same adduct family as vinyl chloride). "
            "References: Hoffler et al., Toxicology 2003 (PMID:12860588); "
            "IARC Monograph Vol 96 (2010)."
        ),
        "exposure": "alcoholic beverages (stone-fruit brandies, sake, sherry), fermented foods, bread crust",
        "class_name": "Aldehydes",
    },
    {
        "id": "Acrolein",
        "label": "Acrolein",
        "type": "Carcinogen",
        "group": "Aldehyde",
        "iarc": "Group 2A",
        "detail": (
            "2-Propenal; smallest alpha,beta-unsaturated aldehyde and "
            "potent Michael acceptor. Major constituents of cigarette "
            "smoke, heated cooking oils, and automobile exhaust; also "
            "generated in vivo as the obligate byproduct of cyclophosphamide "
            "activation (4-hydroxycyclophosphamide fragmentation). Forms "
            "alpha/gamma-OH-1,N2-propano-dG (Acr-dG) exocyclic adducts, "
            "which are mispairing, bulky lesions repaired by NER/BER. "
            "Detoxified by GSTP1 (GS-HPMA) and ALDH2 oxidation to acrylic "
            "acid. References: IARC Monograph Vol 128 (2021); Minko et al., "
            "Chem Res Toxicol 2009 (PMID:19705912)."
        ),
        "exposure": "tobacco smoke, heated cooking oils, wildfire/vehicle exhaust, cyclophosphamide metabolism",
        "class_name": "Aldehydes",
    },
    {
        "id": "Crotonaldehyde",
        "label": "Crotonaldehyde",
        "type": "Carcinogen",
        "group": "Aldehyde",
        "iarc": "Group 2B",
        "detail": (
            "(E)-2-Butenal; alpha,beta-unsaturated aldehyde abundant in "
            "tobacco smoke, alcoholic beverages (congener), and heated "
            "cooking oils. Forms cyclic 1,N2-propano-dG (Cr-dG) exocyclic "
            "adducts via Michael addition/cyclization; adducts are "
            "miscoding and persist in liver, oral mucosa, and lung tissue. "
            "References: Chung et al., Carcinogenesis 1989; IARC Monograph "
            "Vol 63 (1995)."
        ),
        "exposure": "tobacco smoke, alcoholic beverages, heated cooking oils, air pollutant",
        "class_name": "Aldehydes",
    },
    {
        "id": "Furfural",
        "label": "Furfural",
        "type": "Carcinogen",
        "group": "Aldehyde",
        "iarc": "Group 3",
        "detail": (
            "Furan-2-carbaldehyde; Maillard-reaction aldehyde abundant in "
            "roasted coffee, beer, bread, and heated fruit products. "
            "Bioactivated via sequential reduction/hydroxylation to "
            "5-hydroxymethylfurfural (HMF) then sulfated by SULT1A1 to "
            "5-sulfoxymethylfurfural (SMF), a DNA-reactive electrophile "
            "that alkylates dA/dG. Classified IARC Group 3 but included "
            "for mechanistic completeness of dietary aldehyde exposures. "
            "References: Monien et al., Mol Nutr Food Res 2011 "
            "(PMID:21538843); EFSA CONTAM Panel Opinion 2011."
        ),
        "exposure": "roasted coffee, beer, bread, heated fruit/sugar products",
        "class_name": "Aldehydes",
    },
    {
        "id": "MDA",
        "label": "Malondialdehyde (MDA)",
        "type": "Carcinogen",
        "group": "Aldehyde",
        "iarc": "Not classified (endogenous genotoxin)",
        "detail": (
            "Propanedial; ubiquitous lipid-peroxidation end-product and a "
            "direct-acting mutagen. Universal biomarker of oxidative "
            "stress. Reacts with dG to form pyrimido[1,2-a]purin-10(3H)-"
            "one (M1-dG), a miscoding adduct that causes base "
            "substitutions and frameshifts. Elevated in chronic "
            "inflammation, alcoholic liver disease, and hepatocarcinogenesis. "
            "References: Marnett, Mutat Res 1999 (PMID:10064842); "
            "Niedernhofer et al., J Biol Chem 2003."
        ),
        "exposure": "endogenous lipid peroxidation, oxidative stress, chronic inflammation",
        "class_name": "Aldehydes",
    },
    {
        "id": "4_HNE",
        "label": "4-Hydroxynonenal (4-HNE)",
        "type": "Carcinogen",
        "group": "Aldehyde",
        "iarc": "Not classified (endogenous genotoxin)",
        "detail": (
            "(E)-4-Hydroxy-2-nonenal; the most reactive major product of "
            "omega-6 polyunsaturated fatty-acid peroxidation. Bifunctional "
            "Michael acceptor that forms exocyclic 1,N2-propano-dG (HNE-"
            "dG) adducts, preferentially at methylated CpG sites (a "
            "molecular link between oxidative stress and TP53 hotspot "
            "mutations). Detoxified by GSTP1 glutathione conjugation and "
            "ALDH2 oxidation. References: Esterbauer, Schaur & Zollner, "
            "Free Radic Biol Med 1991 (PMID:1937131); Hu et al., PNAS 2002."
        ),
        "exposure": "endogenous omega-6 lipid peroxidation, oxidative stress, alcoholic steatohepatitis",
        "class_name": "Aldehydes",
    },
    # Class 12: Dioxins/AhR
    {
        "id": "TCDD",
        "label": "TCDD",
        "type": "Carcinogen",
        "group": "Dioxin",
        "iarc": "Group 1",
        "detail": "2,3,7,8-Tetrachlorodibenzo-para-dioxin; most potent AhR ligand with picomolar binding affinity. Half-life 7-11 years. Acts as tumor promoter via epigenetic reprogramming, not direct genotoxin.",
        "exposure": "industrial byproduct, incineration, contaminated food (bioaccumulation)",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PCB_126",
        "label": "PCB 126",
        "type": "Carcinogen",
        "group": "PCB",
        "iarc": "Group 1",
        "detail": "3,3',4,4',5-Pentachlorobiphenyl; dioxin-like PCB with WHO-TEF 0.1 relative to TCDD. Persistent organic pollutant that activates AhR pathway.",
        "exposure": "legacy industrial contamination, fish/seafood bioaccumulation",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PCB_169",
        "label": "PCB 169",
        "type": "Carcinogen",
        "group": "PCB",
        "iarc": "Group 1",
        "detail": "3,3',4,4',5,5'-Hexachlorobiphenyl; non-ortho dioxin-like PCB with WHO-TEF 0.03. Potent AhR agonist contributing to total toxic equivalency (TEQ) in serum/food (IARC Monograph 107).",
        "exposure": "legacy transformer/capacitor oils, fatty fish, dairy, breast milk",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PCB_77",
        "label": "PCB 77",
        "type": "Carcinogen",
        "group": "PCB",
        "iarc": "Group 1",
        "detail": "3,3',4,4'-Tetrachlorobiphenyl; coplanar dioxin-like PCB with WHO-TEF 0.0001. Short half-life relative to higher-chlorinated PCBs but remains a marker of recent exposure (IARC Monograph 107).",
        "exposure": "combustion byproduct, contaminated sediments, dietary",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PCB_118",
        "label": "PCB 118",
        "type": "Carcinogen",
        "group": "PCB",
        "iarc": "Group 1",
        "detail": "2,3',4,4',5-Pentachlorobiphenyl; mono-ortho dioxin-like PCB (WHO-TEF 0.00003). Most abundant dioxin-like PCB congener in human serum; used as a biomarker of cumulative PCB exposure (IARC Monograph 107).",
        "exposure": "fish consumption, occupational exposure in electrical industry, background dietary",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PCB_153",
        "label": "PCB 153",
        "type": "Carcinogen",
        "group": "PCB",
        "iarc": "Group 1",
        "detail": "2,2',4,4',5,5'-Hexachlorobiphenyl; non-dioxin-like (NDL) PCB, most abundant PCB in human tissues. Activates CAR/PXR rather than AhR; tumor promoter via oxidative stress and gap-junction inhibition (IARC Monograph 107).",
        "exposure": "fatty fish, meat, dairy; breast milk transfer; highly persistent",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PCB_138",
        "label": "PCB 138",
        "type": "Carcinogen",
        "group": "PCB",
        "iarc": "Group 1",
        "detail": "2,2',3,4,4',5'-Hexachlorobiphenyl; non-dioxin-like PCB co-dominant with PCB-153 in human serum. CAR/PXR activator; commonly reported in NHANES biomonitoring panels (IARC Monograph 107).",
        "exposure": "bioaccumulative dietary exposure, adipose storage, transplacental transfer",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PeCDF_23478",
        "label": "2,3,4,7,8-PeCDF",
        "type": "Carcinogen",
        "group": "Dioxin",
        "iarc": "Group 1",
        "detail": "2,3,4,7,8-Pentachlorodibenzofuran; most potent non-TCDD AhR ligand with WHO-TEF 0.3. Major contributor to Yusho/Yu-cheng poisoning TEQ and dominant PCDF in human adipose (IARC Monograph 100F).",
        "exposure": "PCB combustion, rice-oil contamination (Yusho/Yu-cheng), incineration",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "HCB",
        "label": "Hexachlorobenzene (HCB)",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 2B",
        "detail": "Persistent organochlorine fungicide and industrial byproduct; weak AhR agonist and CYP1A1/CYP1B1 inducer. Linked to hepatic porphyria, thyroid cancer, and elevated breast-cancer risk (IARC Monograph 79).",
        "exposure": "legacy fungicide residues, chlorinated solvent manufacture, dietary fat",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "Lindane",
        "label": "Lindane (γ-HCH)",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 1",
        "detail": "Gamma-hexachlorocyclohexane; organochlorine insecticide classified Group 1 by IARC in 2015 (Monograph 113) based on non-Hodgkin lymphoma evidence. Acts via CAR/PXR-mediated CYP2B/CYP3A induction, GABA-A receptor antagonism, and oxidative stress.",
        "exposure": "legacy agricultural insecticide, head-lice pharmaceutical, persistent in soil and food chain",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "DDT",
        "label": "DDT",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 2A",
        "detail": "Dichlorodiphenyltrichloroethane; classified Group 2A by IARC in 2015 (Monograph 113) based on liver, testicular, and non-Hodgkin lymphoma evidence. Primarily CAR/PXR-mediated tumor promoter; weak AhR activity. Dechlorinated to DDE via CYP2B-mediated metabolism.",
        "exposure": "legacy agricultural/vector-control insecticide, global body-burden biomarker",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "DDE",
        "label": "p,p'-DDE",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 2B",
        "detail": "Principal DDT metabolite and dominant organochlorine in human adipose tissue. Antiandrogenic AR antagonist and CAR/PXR activator; persistent half-life of 10+ years. Biomarker of cumulative DDT exposure (IARC Monograph 113, 2015).",
        "exposure": "bioaccumulation from DDT residues, fatty fish, dairy, adipose storage, transplacental transfer",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "PCP",
        "label": "Pentachlorophenol (PCP)",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 1",
        "detail": "Chlorinated phenol wood preservative and biocide; classified Group 1 by IARC in 2019 (Monograph 117) based on non-Hodgkin lymphoma and multiple myeloma evidence. AhR agonist with redox-cycling tetrachloro-p-benzoquinone metabolite; often co-contaminated with dioxins/furans.",
        "exposure": "sawmill/lumber workers, treated wood products, contaminated soil and groundwater",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "Chlordane",
        "label": "Chlordane (technical)",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 2B",
        "detail": "Technical-grade cyclodiene insecticide mixture; classified Group 2B by IARC in 1991 (Monograph 53). Primarily CAR/PXR activator; liver tumor promoter in rodents; persists in indoor dust from historic termiticide applications.",
        "exposure": "legacy subterranean termiticide, indoor dust residues, dietary bioaccumulation",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "Heptachlor",
        "label": "Heptachlor",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 2B",
        "detail": "Cyclodiene insecticide co-formulated with chlordane; classified Group 2B by IARC in 1991 (Monograph 53). Bioactivated to more toxic heptachlor epoxide by hepatic CYPs; CAR/PXR-mediated tumor promotion.",
        "exposure": "legacy termiticide, soil/sediment residues, dietary exposure via meat and dairy",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "Toxaphene",
        "label": "Toxaphene",
        "type": "Carcinogen",
        "group": "Organochlorine",
        "iarc": "Group 2B",
        "detail": "Complex chlorinated camphene mixture (~670 congeners); classified Group 2B by IARC in 2001 (Monograph 79). Mixed CAR/PXR and weak AhR agonism; GABA-A antagonism; long-range atmospheric transport to Arctic food webs.",
        "exposure": "post-DDT-ban cotton insecticide, persistent in fish/marine mammals, Arctic bioaccumulation",
        "class_name": "Dioxins_AhR",
    },
    # Class 13: Dietary N-Nitroso
    {
        "id": "NDMA",
        "label": "NDMA (N-nitrosodimethylamine)",
        "type": "Carcinogen",
        "group": "Nitrosamine",
        "iarc": "Group 2A",
        "detail": "N-nitrosodimethylamine; potent hepatocarcinogen activated by CYP2E1 (primary) and CYP2A6 (secondary) via alpha-hydroxylation, yielding methylating intermediates that form O6-methyl-dG adducts.",
        "exposure": "processed meats, cured foods, tobacco smoke, contaminated water, drug impurities",
        "class_name": "Dietary_NNitroso",
    },
    {
        "id": "NDEA",
        "label": "NDEA",
        "type": "Carcinogen",
        "group": "Nitrosamine",
        "iarc": "Group 2A",
        "detail": "N-Nitrosodiethylamine; structurally related to NDMA. Found in cured meats, tobacco, and cosmetics. Activated by CYP2E1 alpha-hydroxylation analogous to NDMA pathway.",
        "exposure": "cured meats, tobacco smoke, cosmetics",
        "class_name": "Dietary_NNitroso",
    },
    # Class 14: Chlorinated Solvents
    {
        "id": "TCE",
        "label": "TCE",
        "type": "Carcinogen",
        "group": "Chlorinated_Solvent",
        "iarc": "Group 1",
        "detail": "Trichloroethylene; Group 1 carcinogen for kidney cancer. Unique inverted bioactivation: GSTT1-mediated GSH conjugation is the bioactivation pathway (not detoxification), generating renal-toxic metabolites via CCBL1.",
        "exposure": "dry cleaning, metal degreasing, Superfund sites, contaminated groundwater",
        "class_name": "Chlorinated_Solvents",
    },
    {
        "id": "PCE",
        "label": "PCE",
        "type": "Carcinogen",
        "group": "Chlorinated_Solvent",
        "iarc": "Group 2A",
        "detail": "Tetrachloroethylene (perchloroethylene); current dry cleaning solvent. CYP2E1 oxidation to trichloroacetic acid. IARC Group 2A probable carcinogen.",
        "exposure": "dry cleaning (occupational/residential), contaminated water",
        "class_name": "Chlorinated_Solvents",
    },
    # Class 15: Alkylating Carcinogens
    {
        "id": "Acrylamide",
        "label": "Acrylamide",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 2A",
        "detail": (
            "Water-soluble vinyl monomer generated by Maillard chemistry in "
            "high-temperature cooking of starchy foods (fried/baked potatoes, "
            "coffee roasting) and by polyacrylamide manufacture. Bioactivated "
            "by CYP2E1 epoxidation to glycidamide, the ultimate mutagen. "
            "References: IARC Monograph Vol 60 (1994); Besaratinia & Pfeifer, "
            "Carcinogenesis 2007 (PMID:17872912)."
        ),
        "exposure": "high-temperature cooked starchy foods, tobacco smoke, polyacrylamide industry",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Glycidamide",
        "label": "Glycidamide",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 2A (as ultimate acrylamide metabolite)",
        "detail": (
            "2,3-Epoxypropanamide; direct-acting epoxide and ultimate "
            "mutagenic metabolite of acrylamide. Forms N7-(2-carbamoyl-2-"
            "hydroxyethyl)-dG and N3-GA-dA DNA adducts. Detoxified by GSTP1/"
            "GSTM1 glutathione conjugation. References: Segerbaeck et al., "
            "Carcinogenesis 1995 (PMID:7614537); Besaratinia, Mutat Res 2010."
        ),
        "exposure": "endogenous formation from acrylamide via CYP2E1",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Cyclophosphamide",
        "label": "Cyclophosphamide",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 1",
        "detail": (
            "Oxazaphosphorine prodrug chemotherapeutic. Bioactivated "
            "primarily by CYP2B6 (major) and CYP3A4/2C9 (minor) to "
            "4-hydroxycyclophosphamide, which ring-opens to aldophosphamide "
            "and fragments to phosphoramide mustard (ultimate DNA-crosslinking "
            "agent) plus acrolein. Associated with secondary bladder cancer "
            "and therapy-related AML. References: Roy et al., Drug Metab "
            "Dispos 1999 (PMID:10220571); IARC Monograph Vol 100A (2012)."
        ),
        "exposure": "therapeutic chemotherapy (lymphomas, breast cancer, autoimmune disease)",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Chlorambucil",
        "label": "Chlorambucil",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 1",
        "detail": (
            "Aromatic nitrogen mustard. Forms aziridinium ion via "
            "intramolecular cyclization, yielding N7-guanine monoadducts and "
            "G-G interstrand crosslinks. Detoxified by GSTP1 glutathione "
            "conjugation; GSTP1 Ile105Val (rs1695) modulates treatment "
            "response and therapy-related leukemia risk. References: Povirk "
            "& Shuker, Mutat Res 1994 (PMID:7523912); IARC Monograph Vol "
            "100A (2012)."
        ),
        "exposure": "therapeutic chemotherapy (CLL, Waldenstroem macroglobulinemia, ovarian cancer)",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Sulfur_mustard",
        "label": "Sulfur mustard",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 1",
        "detail": (
            "Bis(2-chloroethyl) sulfide. Chemical-warfare vesicant that "
            "cyclizes to a reactive episulfonium ion, generating N7-mustard-"
            "Gua monoadducts and 5'-d(GNC) interstrand crosslinks. Linked "
            "to lung cancer in exposed workers and veterans. References: "
            "Matijasevic et al., Carcinogenesis 1996 (PMID:8635461); IARC "
            "Monograph Vol 100F (2012)."
        ),
        "exposure": "chemical-weapon exposure, legacy ordnance, historic occupational mustard-gas factories",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Busulfan",
        "label": "Busulfan",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 1",
        "detail": (
            "1,4-Butanediol dimethanesulfonate; bifunctional alkyl sulfonate "
            "used in HSCT conditioning. Hydrolyzes non-enzymatically to a "
            "reactive methanesulfonate intermediate, forming N7-(2,3,4-"
            "trihydroxybutyl)-dG adducts and DNA-DNA interstrand crosslinks. "
            "Detoxified by GSTA1/GSTP1 glutathione conjugation (GSH-Bu). "
            "References: Hassan & Ljungman, Mutat Res 1997 (PMID:12960109); "
            "IARC Monograph Vol 100A (2012)."
        ),
        "exposure": "therapeutic chemotherapy (CML, bone-marrow-transplant conditioning)",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "MNU",
        "label": "MNU (N-methyl-N-nitrosourea)",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 2A",
        "detail": (
            "N-methyl-N-nitrosourea; direct-acting methylating agent "
            "requiring no enzymatic activation. Spontaneously decomposes at "
            "physiological pH to methyldiazohydroxide/methyldiazonium, "
            "alkylating DNA at O6-guanine (G:C to A:T transitions) and "
            "N7-guanine. Canonical model carcinogen for MGMT biology. "
            "References: Kyrtopoulos, Mutat Res 1990 (PMID:2185966); IARC "
            "Monograph Vol 17 (1978)."
        ),
        "exposure": "experimental/research reagent; trace formation in cured meats",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Temozolomide",
        "label": "Temozolomide",
        "type": "Carcinogen",
        "group": "Alkylating",
        "iarc": "Group 2A",
        "detail": (
            "Imidazotetrazine prodrug. Spontaneous pH-dependent hydrolysis "
            "yields MTIC (5-(3-methyltriazen-1-yl)imidazole-4-carboxamide), "
            "which decomposes to methyldiazonium ion. Alkylates O6-guanine; "
            "therapeutic efficacy in glioblastoma tracks MGMT promoter-"
            "methylation status. References: Newlands et al., Cancer Treat "
            "Rev 1997 (PMID:9327140); Hegi et al., NEJM 2005 (PMID:15758009)."
        ),
        "exposure": "therapeutic chemotherapy (glioblastoma, metastatic melanoma)",
        "class_name": "Alkylating_Agents",
    },
]


# Package-native bridge node: the source Wave 2 export expects ``Ethanol`` to
# exist elsewhere. ExposoGraph does not seed that node in the core graph, so
# add it here to keep the acetaldehyde branch traversable.
WAVE2_AUXILIARY_NODES: list[dict[str, Any]] = [
    {
        "id": "Ethanol",
        "label": "Ethanol",
        "type": "Carcinogen",
        "group": "Alcohol",
        "iarc": "Group 1 (alcoholic beverages)",
        "detail": "Upstream alcohol substrate represented to connect ADH-mediated oxidation to intracellular acetaldehyde formation. Included as a package bridge so the aldehyde pathway resolves in standalone and merged graphs.",
        "exposure": "alcoholic beverages",
        "class_name": "Aldehydes",
    }
]


WAVE2_ENZYME_NODES: list[dict[str, Any]] = [
    # Class 11: Aldehydes
    {
        "id": "ALDH2",
        "label": "ALDH2",
        "type": "Enzyme",
        "phase": "II",
        "role": "Detoxification",
        "detail": "Mitochondrial aldehyde dehydrogenase 2; primary acetaldehyde and formaldehyde clearance. rs671 (ALDH2*2): 6-19x acetaldehyde accumulation.",
        "tissue": "liver, esophagus, stomach",
        "variant": "rs671 (Glu504Lys)",
        "class_name": "Aldehydes",
    },
    {
        "id": "ADH1B",
        "label": "ADH1B",
        "type": "Enzyme",
        "phase": "I",
        "role": "Activation",
        "detail": "Alcohol dehydrogenase 1B; oxidizes ethanol to acetaldehyde. ADH1B*2 (rs1229984) fast allele; ESCC OR 2.50 in alcohol consumers.",
        "tissue": "liver, stomach",
        "variant": "rs1229984 (His47Arg)",
        "class_name": "Aldehydes",
    },
    {
        "id": "ADH1C",
        "label": "ADH1C",
        "type": "Enzyme",
        "phase": "I",
        "role": "Activation",
        "detail": "Alcohol dehydrogenase 1C (gamma subunit); contributes to ethanol oxidation. rs698 (Ile349Val) modulates catalytic rate.",
        "tissue": "liver, stomach",
        "variant": "rs698 (Ile349Val)",
        "class_name": "Aldehydes",
    },
    {
        "id": "ADH5",
        "label": "ADH5",
        "type": "Enzyme",
        "phase": "I",
        "role": "Detoxification",
        "detail": "Formaldehyde dehydrogenase (class III ADH); cytosolic formaldehyde clearance via S-hydroxymethylglutathione pathway.",
        "tissue": "liver, ubiquitous",
        "class_name": "Aldehydes",
    },
    {
        "id": "FANCD2",
        "label": "FANCD2",
        "type": "Enzyme",
        "role": "Repair",
        "detail": "Fanconi anemia complementation group D2; repairs acetaldehyde-induced DNA interstrand crosslinks (ICLs).",
        "tissue": "bone marrow, ubiquitous",
        "group": "DNA Repair (FA)",
        "class_name": "Aldehydes",
    },
    # Class 12: Dioxins/AhR
    {
        "id": "AHR",
        "label": "AHR",
        "type": "Enzyme",
        "role": "Transcription Factor",
        "detail": "Aryl hydrocarbon receptor; ligand-activated transcription factor. Binds TCDD/dioxin-like compounds and induces CYP1A1/CYP1B1 via XRE elements.",
        "tissue": "liver, lung, skin, immune cells",
        "variant": "rs2066853 (Arg554Lys)",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "ARNT",
        "label": "ARNT",
        "type": "Enzyme",
        "role": "Transcription Factor",
        "detail": "AhR nuclear translocator (HIF-1beta); obligate co-factor for AhR. Forms AhR-ARNT heterodimer that binds XRE/DRE elements.",
        "tissue": "ubiquitous",
        "variant": "rs2228099 (Ile471Val)",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "AHRR",
        "label": "AHRR",
        "type": "Enzyme",
        "role": "Transcription Factor",
        "detail": "AhR repressor; negative feedback regulator. cg05575921 CpG hypomethylation is a validated blood-based biomarker for tobacco/TCDD exposure.",
        "tissue": "ubiquitous",
        "class_name": "Dioxins_AhR",
    },
    # Class 14: Chlorinated Solvents
    {
        "id": "CCBL1",
        "label": "CCBL1",
        "type": "Enzyme",
        "phase": "II",
        "role": "Bioactivation",
        "detail": "Cysteine conjugate beta-lyase 1; renal proximal tubule enzyme that bioactivates DCVC to reactive DCVT thiol, causing DNA adducts and VHL mutations. rs2293968/rs2280841/rs2259043/rs941960 tagging SNPs.",
        "tissue": "kidney (proximal tubule)",
        "class_name": "Chlorinated_Solvents",
    },
]


WAVE2_METABOLITE_NODES: list[dict[str, Any]] = [
    # Class 11: Aldehydes
    {
        "id": "Acetaldehyde_int",
        "label": "Acetaldehyde (intracellular)",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "Intracellular acetaldehyde; reactive electrophile that forms DNA adducts and protein cross-links.",
        "class_name": "Aldehydes",
    },
    {
        "id": "Formate",
        "label": "Formate",
        "type": "Metabolite",
        "reactivity": "Low",
        "detail": "Terminal product of formaldehyde oxidation; urinary biomarker of formaldehyde exposure.",
        "class_name": "Aldehydes",
    },
    {
        "id": "Vinyl_carbamate",
        "label": "Vinyl carbamate",
        "type": "Metabolite",
        "reactivity": "Intermediate",
        "detail": "CYP2E1-derived dehydrogenation product of urethane; proximal intermediate re-epoxidized to vinyl carbamate epoxide (PMID:12860588).",
        "class_name": "Aldehydes",
    },
    {
        "id": "Vinyl_carbamate_epoxide",
        "label": "Vinyl carbamate epoxide",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "Ultimate mutagenic metabolite of urethane; alkylates dA and dC to yield etheno-DNA adducts equivalent to those formed by vinyl chloride (PMID:12860588).",
        "class_name": "Aldehydes",
    },
    {
        "id": "SMF",
        "label": "5-sulfoxymethylfurfural (SMF)",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "SULT1A1-derived sulfate ester of 5-hydroxymethylfurfural (HMF); DNA-reactive electrophile that alkylates dA/dG. HMF itself arises from hydroxylation of furfural in heated foods (PMID:21538843).",
        "class_name": "Aldehydes",
    },
    # Class 12: Dioxins/AhR
    {
        "id": "AhR_ARNT_complex",
        "label": "AhR-ARNT complex",
        "type": "Metabolite",
        "reactivity": "Low",
        "detail": "Active transcription factor heterodimer; binds XRE/DRE elements to induce CYP1A1/CYP1B1. Not a metabolite per se but modeled as an intermediate complex.",
        "class_name": "Dioxins_AhR",
    },
    # Class 13: Dietary N-Nitroso
    {
        "id": "Hydroxymethylnitrosamine",
        "label": "Hydroxymethylnitrosamine",
        "type": "Metabolite",
        "reactivity": "High",
        "detail": "Alpha-hydroxylation product of NDMA/NDEA; proximate carcinogen that spontaneously decomposes to methyldiazonium + formaldehyde.",
        "class_name": "Dietary_NNitroso",
    },
    {
        "id": "Methyldiazonium",
        "label": "Methyldiazonium ion",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "Ultimate reactive alkylating species produced by nitrosamine alpha-hydroxylation. Methylates guanine O6 to form the primary mutagenic lesion for the NDMA branch.",
        "class_name": "Dietary_NNitroso",
    },
    # Class 14: Chlorinated Solvents
    {
        "id": "DCVG",
        "label": "DCVG",
        "type": "Metabolite",
        "reactivity": "Moderate",
        "detail": "S-(1,2-dichlorovinyl)glutathione; GSTT1 conjugation product of TCE. Transported to kidney for further processing.",
        "class_name": "Chlorinated_Solvents",
    },
    {
        "id": "DCVC",
        "label": "DCVC",
        "type": "Metabolite",
        "reactivity": "High",
        "detail": "S-(1,2-dichlorovinyl)-L-cysteine; renal processing product of DCVG. Substrate for CCBL1 beta-lyase bioactivation.",
        "class_name": "Chlorinated_Solvents",
    },
    {
        "id": "Chloral_hydrate",
        "label": "Chloral hydrate",
        "type": "Metabolite",
        "reactivity": "Moderate",
        "detail": "CYP2E1 oxidative metabolite of TCE; hepatotoxic intermediate oxidized to trichloroacetic acid.",
        "class_name": "Chlorinated_Solvents",
    },
    {
        "id": "TCA",
        "label": "TCA",
        "type": "Metabolite",
        "reactivity": "Low",
        "detail": "Trichloroacetic acid; stable urinary biomarker for TCE/PCE exposure.",
        "class_name": "Chlorinated_Solvents",
    },
    # Class 15: Alkylating Carcinogens
    {
        "id": "4OH_cyclophosphamide",
        "label": "4-Hydroxycyclophosphamide",
        "type": "Metabolite",
        "reactivity": "Intermediate",
        "detail": "Primary CYP2B6/CYP3A4 activation product of cyclophosphamide; equilibrates with aldophosphamide and fragments to phosphoramide mustard plus acrolein (PMID:10220571).",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Phosphoramide_mustard",
        "label": "Phosphoramide mustard",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "Ultimate DNA-crosslinking species from cyclophosphamide activation; bifunctional nitrogen mustard that forms G-G interstrand crosslinks via aziridinium intermediates (PMID:10220571).",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Chlorambucil_aziridinium",
        "label": "Chlorambucil aziridinium",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "Cyclic aziridinium cation formed by intramolecular cyclization of chlorambucil; alkylates N7-guanine and generates interstrand crosslinks (PMID:7523912).",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Mustard_episulfonium",
        "label": "Sulfur-mustard episulfonium ion",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "Cyclic sulfonium intermediate formed by intramolecular cyclization of sulfur mustard; ultimate electrophile for N7-Gua monoadducts and 5'-d(GNC) interstrand crosslinks (PMID:8635461).",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "Busulfan_methanesulfonate",
        "label": "Busulfan methanesulfonate intermediate",
        "type": "Metabolite",
        "reactivity": "Very High",
        "detail": "Reactive bifunctional alkylator generated by hydrolytic displacement of methanesulfonate groups from busulfan; yields N7-THPG adducts and DNA-DNA interstrand crosslinks (PMID:12960109).",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "MTIC",
        "label": "MTIC",
        "type": "Metabolite",
        "reactivity": "High",
        "detail": "5-(3-methyltriazen-1-yl)imidazole-4-carboxamide; monomethyltriazene formed by spontaneous hydrolysis of temozolomide, which in turn decomposes to methyldiazonium (PMID:9327140).",
        "class_name": "Alkylating_Agents",
    },
]


WAVE2_DNA_ADDUCT_NODES: list[dict[str, Any]] = [
    # Class 11: Aldehydes
    {
        "id": "N2_HOMedG",
        "label": "N2-HOMedG",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "N2-hydroxymethyl-dG; formaldehyde-DNA adduct. ~20-fold increase in ALDH2/ADH5 double-knockout mice.",
        "class_name": "Aldehydes",
    },
    {
        "id": "N2_ethylidene_dG",
        "label": "N2-ethylidene-dG",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "Acetaldehyde-DNA adduct; forms interstrand crosslinks (ICLs) requiring Fanconi anemia pathway repair.",
        "class_name": "Aldehydes",
    },
    {
        "id": "Acr_dG",
        "label": "alpha/gamma-OH-PdG (Acr-dG)",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "alpha- and gamma-hydroxy-1,N2-propano-2'-deoxyguanosine; major Michael-addition adduct from acrolein reacting at N1/N2 of guanine. Miscoding in ring-opened form; repaired by NER with XRCC1/BER backup (PMID:19705912).",
        "class_name": "Aldehydes",
    },
    {
        "id": "Cr_dG",
        "label": "Cr-1,N2-PdG (Cr-dG)",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "Cyclic 1,N2-propano-dG adduct from crotonaldehyde Michael addition; miscoding adduct detected in liver, oral mucosa, and lung of smokers (IARC Monograph Vol 63; Chung et al. 1989).",
        "class_name": "Aldehydes",
    },
    {
        "id": "M1_dG",
        "label": "M1-dG (pyrimido-purinone)",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "3-(2'-deoxy-beta-D-erythro-pentofuranosyl)pyrimido[1,2-a]purin-10(3H)-one; principal MDA-DNA adduct. Miscoding lesion causing base substitutions and frameshifts; basal levels in healthy tissues provide a baseline for oxidative-stress-driven mutagenesis (PMID:10064842).",
        "class_name": "Aldehydes",
    },
    {
        "id": "HNE_dG",
        "label": "HNE-1,N2-PdG (HNE-dG)",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "Exocyclic 1,N2-propano-dG adduct from 4-HNE Michael addition/cyclization; formation is strongly enhanced at methylated CpG sites, mechanistically linking lipid peroxidation to TP53 hotspot mutations (Hu et al., PNAS 2002; PMID:1937131).",
        "class_name": "Aldehydes",
    },
    # Class 13: Dietary N-Nitroso
    {
        "id": "O6_methyl_dG",
        "label": "O6-methyl-dG",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "O6-methyldeoxyguanosine; primary mutagenic lesion from nitrosamine alkylation. Causes G:C to A:T transitions. Repaired by MGMT direct reversal.",
        "class_name": "Dietary_NNitroso",
    },
    # Class 14: Chlorinated Solvents
    {
        "id": "Renal_DNA_damage",
        "label": "Renal DNA damage",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "DCVC-derived DCVT thiol forms DNA adducts in renal proximal tubule; associated with VHL tumor suppressor mutations in renal cell carcinoma.",
        "class_name": "Chlorinated_Solvents",
    },
    # Class 15: Alkylating Carcinogens
    {
        "id": "N7_GA_dG",
        "label": "N7-GA-dG",
        "type": "DNA_Adduct",
        "reactivity": "High",
        "detail": "N7-(2-carbamoyl-2-hydroxyethyl)-2'-deoxyguanosine; dominant glycidamide-DNA adduct from acrylamide exposure. Depurinates to abasic site repaired by BER (PMID:17872912).",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "N7_methyl_dG",
        "label": "N7-methyl-dG",
        "type": "DNA_Adduct",
        "reactivity": "Moderate",
        "detail": "N7-methyl-2'-deoxyguanosine; most abundant methylation lesion from MNU and temozolomide via methyldiazonium. Non-mutagenic itself but labile; depurinates to abasic sites repaired by BER (XRCC1 scaffold).",
        "class_name": "Alkylating_Agents",
    },
    {
        "id": "DNA_ICL_mustard",
        "label": "G-G interstrand crosslink (mustard-type)",
        "type": "DNA_Adduct",
        "reactivity": "Very High",
        "detail": "Bifunctional N7-alkyl-dG interstrand crosslink generated by phosphoramide mustard (cyclophosphamide), chlorambucil aziridinium, sulfur-mustard episulfonium, and the busulfan methanesulfonate intermediate. Requires Fanconi-anemia/HR repair; mechanistic basis for therapeutic cytotoxicity (PMID:10220571, PMID:7523912, PMID:8635461, PMID:12960109).",
        "class_name": "Alkylating_Agents",
    },
]


WAVE2_PATHWAY_NODES: list[dict[str, Any]] = [
    {
        "id": "aldehyde_pathway",
        "label": "Aldehyde Carcinogenesis",
        "type": "Pathway",
        "detail": "Composite pathway covering formaldehyde and acetaldehyde metabolism, DNA adduct formation, and Fanconi anemia repair.",
        "class_name": "Aldehydes",
    },
    {
        "id": "dioxin_pathway",
        "label": "AhR & CAR/PXR Transcriptional Induction",
        "type": "Pathway",
        "detail": "Nuclear-receptor-mediated tumor promotion: AhR/ARNT-driven CYP1A1/CYP1B1 induction by dioxin-like ligands, plus CAR/PXR-driven CYP2B6/CYP3A4 induction by non-dioxin-like PCBs and organochlorine pesticides (DDT, DDE, Lindane, Chlordane, Heptachlor, Toxaphene). Epigenetic reprogramming rather than direct DNA damage.",
        "class_name": "Dioxins_AhR",
    },
    {
        "id": "nitrosamine_pathway",
        "label": "N-Nitroso Compound Activation",
        "type": "Pathway",
        "detail": "CYP2E1/CYP2A6-mediated alpha-hydroxylation pathway for dietary nitrosamines; generates alkylating methyldiazonium ions.",
        "class_name": "Dietary_NNitroso",
    },
    {
        "id": "chlorinated_pathway",
        "label": "Chlorinated Solvent Dual Pathway",
        "type": "Pathway",
        "detail": "Dual metabolic pathway: CYP2E1 oxidative (hepatic) and GSTT1 GSH conjugation bioactivation (renal). The GSH pathway is uniquely a bioactivation route.",
        "class_name": "Chlorinated_Solvents",
    },
    {
        "id": "alkylating_pathway",
        "label": "Alkylating Agent Carcinogenesis",
        "type": "Pathway",
        "detail": "Composite pathway covering direct-acting methylators (MNU), prodrug-activated alkylators (cyclophosphamide, temozolomide), chemical-warfare mustards (sulfur mustard), bifunctional sulfonates (busulfan), and CYP2E1-activated acrylamide/glycidamide converging on N7/O6-guanine adducts and DNA interstrand crosslinks.",
        "class_name": "Alkylating_Agents",
    },
]


# ── Edge Data ─────────────────────────────────────────────────────────────

WAVE2_EDGES: list[dict[str, Any]] = [
    # ── Class 11: Aldehydes (13 edges) ────────────────────────────────────
    {"source": "Ethanol", "target": "Acetaldehyde_int", "type": "ACTIVATES",
     "carcinogen": "Acetaldehyde", "evidence": "ADH1B/ADH1C-mediated oxidation of ethanol produces intracellular acetaldehyde; represented explicitly in the source Wave 2 expansion.", "class_name": "Aldehydes"},
    {"source": "ADH1B", "target": "Acetaldehyde_int", "type": "ACTIVATES",
     "carcinogen": "Acetaldehyde", "evidence": "ADH1B oxidizes ethanol to acetaldehyde; *2 (rs1229984) fast allele ESCC OR 2.50 (PMC4122263)", "class_name": "Aldehydes"},
    {"source": "ADH1C", "target": "Acetaldehyde_int", "type": "ACTIVATES",
     "carcinogen": "Acetaldehyde", "evidence": "ADH1C contributes to ethanol oxidation (rs698 Ile349Val)", "class_name": "Aldehydes"},
    {"source": "ALDH2", "target": "Acetaldehyde_int", "type": "DETOXIFIES",
     "carcinogen": "Acetaldehyde", "evidence": "Primary clearance; rs671 ALDH2*2 causes 6-19x acetaldehyde accumulation (PMC4376259)", "class_name": "Aldehydes"},
    {"source": "Acetaldehyde_int", "target": "N2_ethylidene_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "Acetaldehyde", "evidence": "Acetaldehyde forms DNA interstrand crosslinks via N2-ethylidene-dG (PMC5324749)", "class_name": "Aldehydes"},
    {"source": "FANCD2", "target": "N2_ethylidene_dG", "type": "REPAIRS",
     "carcinogen": "Acetaldehyde", "evidence": "Fanconi anemia pathway repairs acetaldehyde-induced ICLs (PMC5324749)", "class_name": "Aldehydes"},
    {"source": "Formaldehyde", "target": "N2_HOMedG", "type": "FORMS_ADDUCT",
     "carcinogen": "Formaldehyde", "evidence": "~20-fold N2-HOMedG increase in ALDH2/ADH5 double-KO mice (PMID:21222454)", "class_name": "Aldehydes"},
    {"source": "ALDH2", "target": "Formaldehyde", "type": "DETOXIFIES",
     "carcinogen": "Formaldehyde", "evidence": "Mitochondrial formaldehyde oxidation to formate (PMC7758861)", "class_name": "Aldehydes"},
    {"source": "ADH5", "target": "Formaldehyde", "type": "DETOXIFIES",
     "carcinogen": "Formaldehyde", "evidence": "Cytosolic formaldehyde dehydrogenase; S-hydroxymethylglutathione pathway (PMC7758861)", "class_name": "Aldehydes"},
    {"source": "ADH5", "target": "Formate", "type": "ACTIVATES",
     "carcinogen": "Formaldehyde", "evidence": "Formaldehyde to S-hydroxymethylglutathione to formate", "class_name": "Aldehydes"},
    {"source": "ALDH2", "target": "Formate", "type": "ACTIVATES",
     "carcinogen": "Formaldehyde", "evidence": "Mitochondrial formaldehyde clearance producing formate", "class_name": "Aldehydes"},
    {"source": "Formaldehyde", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "Formaldehyde", "class_name": "Aldehydes"},
    {"source": "Acetaldehyde", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "Acetaldehyde", "class_name": "Aldehydes"},

    # ── Class 11 expansion: Urethane ──────────────────────────────────────
    {"source": "Urethane", "target": "Vinyl_carbamate", "type": "ACTIVATES",
     "carcinogen": "Urethane", "evidence": "CYP2E1 dehydrogenation of urethane to vinyl carbamate is the first obligatory activation step (PMID:12860588; IARC Monograph 96)", "class_name": "Aldehydes"},
    {"source": "CYP2E1", "target": "Vinyl_carbamate", "type": "ACTIVATES",
     "carcinogen": "Urethane", "evidence": "Primary P450 for urethane dehydrogenation; abolished in CYP2E1-null mice (PMID:12860588)", "class_name": "Aldehydes"},
    {"source": "Vinyl_carbamate", "target": "Vinyl_carbamate_epoxide", "type": "ACTIVATES",
     "carcinogen": "Urethane", "evidence": "Second CYP2E1 epoxidation generates the ultimate electrophile (PMID:12860588)", "class_name": "Aldehydes"},
    {"source": "Vinyl_carbamate_epoxide", "target": "etheno_dA", "type": "FORMS_ADDUCT",
     "carcinogen": "Urethane", "evidence": "Forms 1,N6-etheno-dA adducts mechanistically identical to vinyl-chloride chemistry (PMID:12860588; IARC Monograph 96)", "class_name": "Aldehydes"},
    {"source": "Vinyl_carbamate_epoxide", "target": "etheno_dC", "type": "FORMS_ADDUCT",
     "carcinogen": "Urethane", "evidence": "Forms 3,N4-etheno-dC adducts detected in liver and lung of urethane-exposed rodents (PMID:12860588)", "class_name": "Aldehydes"},
    {"source": "GSTT1", "target": "Vinyl_carbamate_epoxide", "type": "DETOXIFIES",
     "carcinogen": "Urethane", "evidence": "GSH conjugation of vinyl-carbamate epoxide; GSTT1-null allele elevates etheno-adduct burden (IARC Monograph 96)", "class_name": "Aldehydes"},
    {"source": "Urethane", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "Urethane", "class_name": "Aldehydes"},

    # ── Class 11 expansion: Acrolein ──────────────────────────────────────
    {"source": "Acrolein", "target": "Acr_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "Acrolein", "evidence": "Direct Michael addition yields alpha/gamma-hydroxy-1,N2-propano-dG adducts; dominant adduct class in smokers' lung DNA (PMID:19705912; IARC Monograph 128)", "class_name": "Aldehydes"},
    {"source": "ALDH2", "target": "Acrolein", "type": "DETOXIFIES",
     "carcinogen": "Acrolein", "evidence": "ALDH2 oxidation of acrolein to acrylic acid; ALDH2*2 carriers show elevated acrolein-protein adducts (IARC Monograph 128)", "class_name": "Aldehydes"},
    {"source": "GSTP1", "target": "Acrolein", "type": "DETOXIFIES",
     "carcinogen": "Acrolein", "evidence": "GSH conjugation yields 3-hydroxypropylmercapturic acid (HPMA), the principal urinary biomarker (IARC Monograph 128)", "class_name": "Aldehydes"},
    {"source": "XRCC1", "target": "Acr_dG", "type": "REPAIRS",
     "carcinogen": "Acrolein", "evidence": "BER scaffold response to ring-opened Acr-dG adducts with NER backup for bulky ring-closed forms (PMID:19705912)", "class_name": "Aldehydes"},
    {"source": "Cyclophosphamide", "target": "Acrolein", "type": "ACTIVATES",
     "carcinogen": "Cyclophosphamide", "evidence": "Acrolein is the obligate byproduct of aldophosphamide beta-elimination during cyclophosphamide activation; drives hemorrhagic cystitis and contributes to secondary bladder cancer (PMID:10220571)", "class_name": "Aldehydes"},
    {"source": "Acrolein", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "Acrolein", "class_name": "Aldehydes"},

    # ── Class 11 expansion: Crotonaldehyde ────────────────────────────────
    {"source": "Crotonaldehyde", "target": "Cr_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "Crotonaldehyde", "evidence": "Michael addition/cyclization yields 1,N2-propano-dG exocyclic adduct detected in smokers' oral and lung tissue (Chung et al., Carcinogenesis 1989; IARC Monograph 63)", "class_name": "Aldehydes"},
    {"source": "ALDH2", "target": "Crotonaldehyde", "type": "DETOXIFIES",
     "carcinogen": "Crotonaldehyde", "evidence": "ALDH2 oxidation to crotonic acid; ALDH2*2 polymorphism amplifies Cr-dG burden (IARC Monograph 63)", "class_name": "Aldehydes"},
    {"source": "GSTP1", "target": "Crotonaldehyde", "type": "DETOXIFIES",
     "carcinogen": "Crotonaldehyde", "evidence": "GSH conjugation to 3-hydroxy-1-methylpropylmercapturic acid (HMPMA), a validated urinary biomarker for tobacco-smoke crotonaldehyde exposure", "class_name": "Aldehydes"},
    {"source": "XRCC1", "target": "Cr_dG", "type": "REPAIRS",
     "carcinogen": "Crotonaldehyde", "evidence": "NER/BER handling of exocyclic propano-dG adducts (IARC Monograph 63)", "class_name": "Aldehydes"},
    {"source": "Crotonaldehyde", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "Crotonaldehyde", "class_name": "Aldehydes"},

    # ── Class 11 expansion: Furfural ──────────────────────────────────────
    {"source": "Furfural", "target": "SMF", "type": "ACTIVATES",
     "carcinogen": "Furfural", "evidence": "Furfural is reduced/hydroxylated to 5-hydroxymethylfurfural (HMF), then SULT1A1 sulfates HMF to the DNA-reactive 5-sulfoxymethylfurfural (SMF) (PMID:21538843)", "class_name": "Aldehydes"},
    {"source": "SULT1A1", "target": "SMF", "type": "ACTIVATES",
     "carcinogen": "Furfural", "evidence": "SULT1A1 is the principal sulfotransferase generating SMF in human liver and colon (PMID:21538843; EFSA CONTAM Panel 2011)", "class_name": "Aldehydes"},
    {"source": "Furfural", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "Furfural", "class_name": "Aldehydes"},

    # ── Class 11 expansion: Malondialdehyde (MDA) ─────────────────────────
    {"source": "MDA", "target": "M1_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "MDA", "evidence": "Direct reaction of MDA with dG yields pyrimido[1,2-a]purin-10(3H)-one (M1-dG); quantitative biomarker of oxidative-stress mutagenesis (PMID:10064842)", "class_name": "Aldehydes"},
    {"source": "XRCC1", "target": "M1_dG", "type": "REPAIRS",
     "carcinogen": "MDA", "evidence": "NER plus BER handle M1-dG; repair deficiency elevates spontaneous mutation frequency (Niedernhofer et al., JBC 2003)", "class_name": "Aldehydes"},
    {"source": "MDA", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "MDA", "class_name": "Aldehydes"},

    # ── Class 11 expansion: 4-Hydroxynonenal (4-HNE) ──────────────────────
    {"source": "4_HNE", "target": "HNE_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "4_HNE", "evidence": "Michael addition/cyclization yields 1,N2-propano-dG (HNE-dG); selectively enhanced at methylated CpG sites, linking lipid peroxidation to TP53 hotspot spectrum (PMID:1937131; Hu et al., PNAS 2002)", "class_name": "Aldehydes"},
    {"source": "GSTP1", "target": "4_HNE", "type": "DETOXIFIES",
     "carcinogen": "4_HNE", "evidence": "GSH conjugation is the dominant detoxification route; GSTA4 is a preferred isoform but GSTP1 provides substantial activity in human tissue (PMID:1937131)", "class_name": "Aldehydes"},
    {"source": "ALDH2", "target": "4_HNE", "type": "DETOXIFIES",
     "carcinogen": "4_HNE", "evidence": "ALDH2 oxidation to 4-hydroxynon-2-enoic acid; ALDH2*2 carriers accumulate 4-HNE-protein adducts in alcoholic liver disease (PMID:1937131)", "class_name": "Aldehydes"},
    {"source": "4_HNE", "target": "aldehyde_pathway", "type": "PATHWAY",
     "carcinogen": "4_HNE", "class_name": "Aldehydes"},

    # ── Class 12: Dioxins/AhR (10 edges) ──────────────────────────────────
    {"source": "TCDD", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "TCDD", "evidence": "TCDD binds AhR with picomolar Kd; displaces HSP90 (PMC3748760)", "class_name": "Dioxins_AhR"},
    {"source": "PCB_126", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "PCB_126", "evidence": "Dioxin-like PCB; WHO-TEF 0.1 relative to TCDD (PMC3748760)", "class_name": "Dioxins_AhR"},
    {"source": "AHR", "target": "ARNT", "type": "ACTIVATES",
     "carcinogen": "TCDD", "evidence": "Nuclear translocation and dimerization with ARNT", "class_name": "Dioxins_AhR"},
    {"source": "AHR", "target": "AhR_ARNT_complex", "type": "ACTIVATES",
     "carcinogen": "TCDD", "evidence": "AhR-ARNT heterodimerization forming active transcription factor", "class_name": "Dioxins_AhR"},
    {"source": "ARNT", "target": "AhR_ARNT_complex", "type": "ACTIVATES",
     "carcinogen": "TCDD", "evidence": "Obligate co-factor for AhR transcriptional complex", "class_name": "Dioxins_AhR"},
    {"source": "AHR", "target": "CYP1A1", "type": "INDUCES",
     "carcinogen": "TCDD", "evidence": "XRE binding induces CYP1A1; CpG demethylation persists 40+ days (Sci Rep 2016, srep34989)", "class_name": "Dioxins_AhR"},
    {"source": "AHR", "target": "CYP1B1", "type": "INDUCES",
     "carcinogen": "TCDD", "evidence": "Amplifies estrogen/PAH bioactivation via CYP1B1 induction (PMC3748760)", "class_name": "Dioxins_AhR"},
    {"source": "AHRR", "target": "AHR", "type": "INHIBITS",
     "carcinogen": "TCDD", "evidence": "Competes with AhR for ARNT binding; negative feedback loop", "class_name": "Dioxins_AhR"},
    {"source": "TCDD", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "TCDD", "class_name": "Dioxins_AhR"},
    {"source": "PCB_126", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PCB_126", "class_name": "Dioxins_AhR"},
    {"source": "PCB_169", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "PCB_169", "evidence": "Non-ortho dioxin-like PCB; WHO-TEF 0.03 relative to TCDD (van den Berg 2006; IARC Monograph 107, 2016)", "class_name": "Dioxins_AhR"},
    {"source": "PCB_77", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "PCB_77", "evidence": "Coplanar dioxin-like PCB; WHO-TEF 0.0001; CYP1A1 inducer in hepatocytes (IARC Monograph 107, 2016)", "class_name": "Dioxins_AhR"},
    {"source": "PCB_118", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "PCB_118", "evidence": "Mono-ortho dioxin-like PCB; WHO-TEF 0.00003; most abundant dioxin-like congener in human serum (IARC Monograph 107, 2016)", "class_name": "Dioxins_AhR"},
    {"source": "PeCDF_23478", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "PeCDF_23478", "evidence": "Most potent non-TCDD AhR ligand; WHO-TEF 0.3; dominant PCDF in Yusho/Yu-cheng cohorts (IARC Monograph 100F, 2012)", "class_name": "Dioxins_AhR"},
    {"source": "HCB", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "HCB", "evidence": "Weak AhR agonist; induces CYP1A1/CYP1A2/CYP1B1 in rodent liver; associated with hepatic porphyria (IARC Monograph 79, 2001)", "class_name": "Dioxins_AhR"},
    {"source": "PCB_169", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PCB_169", "class_name": "Dioxins_AhR"},
    {"source": "PCB_77", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PCB_77", "class_name": "Dioxins_AhR"},
    {"source": "PCB_118", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PCB_118", "class_name": "Dioxins_AhR"},
    {"source": "PCB_153", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PCB_153", "evidence": "Non-dioxin-like PCB; CAR/PXR-mediated tumor promotion rather than AhR activation (IARC Monograph 107, 2016)", "class_name": "Dioxins_AhR"},
    {"source": "PCB_138", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PCB_138", "evidence": "Non-dioxin-like PCB co-dominant with PCB-153; CAR/PXR activator (IARC Monograph 107, 2016)", "class_name": "Dioxins_AhR"},
    {"source": "PeCDF_23478", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PeCDF_23478", "class_name": "Dioxins_AhR"},
    {"source": "HCB", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "HCB", "class_name": "Dioxins_AhR"},
    {"source": "PCP", "target": "AHR", "type": "ACTIVATES",
     "carcinogen": "PCP", "evidence": "Pentachlorophenol is a moderate AhR agonist; tetrachloro-p-benzoquinone metabolite contributes to redox cycling (IARC Monograph 117, 2019)", "class_name": "Dioxins_AhR"},
    {"source": "Lindane", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "Lindane", "evidence": "Non-AhR organochlorine; CAR/PXR-mediated CYP2B/CYP3A induction and GABA-A antagonism (IARC Monograph 113, 2015)", "class_name": "Dioxins_AhR"},
    {"source": "DDT", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "DDT", "evidence": "Primarily CAR/PXR-mediated tumor promotion; dechlorinated to DDE via hepatic CYP2B (IARC Monograph 113, 2015)", "class_name": "Dioxins_AhR"},
    {"source": "DDE", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "DDE", "evidence": "Antiandrogenic AR antagonist and CAR/PXR activator; dominant DDT body-burden metabolite (IARC Monograph 113, 2015)", "class_name": "Dioxins_AhR"},
    {"source": "PCP", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "PCP", "class_name": "Dioxins_AhR"},
    {"source": "Chlordane", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "Chlordane", "evidence": "CAR-mediated phenobarbital-like liver tumor promotion (IARC Monograph 53, 1991)", "class_name": "Dioxins_AhR"},
    {"source": "Heptachlor", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "Heptachlor", "evidence": "Bioactivated to heptachlor epoxide; CAR/PXR-mediated hepatic CYP induction (IARC Monograph 53, 1991)", "class_name": "Dioxins_AhR"},
    {"source": "Toxaphene", "target": "dioxin_pathway", "type": "PATHWAY",
     "carcinogen": "Toxaphene", "evidence": "Mixed CAR/PXR and weak AhR activity; GABA-A antagonism (IARC Monograph 79, 2001)", "class_name": "Dioxins_AhR"},

    # ── Class 13: Dietary N-Nitroso (9 edges) ─────────────────────────────
    {"source": "NDMA", "target": "Hydroxymethylnitrosamine", "type": "ACTIVATES",
     "carcinogen": "NDMA", "evidence": "CYP2E1 alpha-hydroxylation; confirmed in HepG2-CYP2E1 cells (PMID:40390554)", "class_name": "Dietary_NNitroso"},
    {"source": "CYP2E1", "target": "Hydroxymethylnitrosamine", "type": "ACTIVATES",
     "carcinogen": "NDMA", "evidence": "Primary P450; required for NDMA genotoxicity (PMID:40390554)", "class_name": "Dietary_NNitroso"},
    {"source": "CYP2A6", "target": "NDMA", "type": "ACTIVATES",
     "carcinogen": "NDMA", "evidence": "Secondary activation at higher NDMA concentrations (ATSDR TP-141)", "class_name": "Dietary_NNitroso"},
    {"source": "Hydroxymethylnitrosamine", "target": "Methyldiazonium", "type": "ACTIVATES",
     "carcinogen": "NDMA", "evidence": "Spontaneous decomposition to formaldehyde + methyldiazonium", "class_name": "Dietary_NNitroso"},
    {"source": "Methyldiazonium", "target": "O6_methyl_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "NDMA", "evidence": "Methylates O6-G; causes G:C to A:T transitions (PMID:40390554)", "class_name": "Dietary_NNitroso"},
    {"source": "MGMT", "target": "O6_methyl_dG", "type": "REPAIRS",
     "carcinogen": "NDMA", "evidence": "Direct reversal repair (suicide enzyme mechanism)", "class_name": "Dietary_NNitroso"},
    {"source": "NDEA", "target": "Hydroxymethylnitrosamine", "type": "ACTIVATES",
     "carcinogen": "NDEA", "evidence": "Analogous CYP2E1 alpha-hydroxylation pathway as NDMA", "class_name": "Dietary_NNitroso"},
    {"source": "NDMA", "target": "nitrosamine_pathway", "type": "PATHWAY",
     "carcinogen": "NDMA", "class_name": "Dietary_NNitroso"},
    {"source": "NDEA", "target": "nitrosamine_pathway", "type": "PATHWAY",
     "carcinogen": "NDEA", "class_name": "Dietary_NNitroso"},

    # ── Class 14: Chlorinated Solvents (12 edges) ─────────────────────────
    {"source": "TCE", "target": "Chloral_hydrate", "type": "ACTIVATES",
     "carcinogen": "TCE", "evidence": "CYP2E1 oxidative pathway; hepatic metabolism", "class_name": "Chlorinated_Solvents"},
    {"source": "CYP2E1", "target": "Chloral_hydrate", "type": "ACTIVATES",
     "carcinogen": "TCE", "evidence": "Primary CYP2E1 oxidation of TCE", "class_name": "Chlorinated_Solvents"},
    {"source": "Chloral_hydrate", "target": "TCA", "type": "ACTIVATES",
     "carcinogen": "TCE", "evidence": "Oxidized to trichloroacetic acid; urinary biomarker", "class_name": "Chlorinated_Solvents"},
    {"source": "TCE", "target": "DCVG", "type": "ACTIVATES",
     "carcinogen": "TCE", "evidence": "BIOACTIVATION (not detoxification); GSTT1-active genotype OR 1.88 for RCC (PMID:20663906)", "class_name": "Chlorinated_Solvents"},
    {"source": "GSTT1", "target": "DCVG", "type": "ACTIVATES",
     "carcinogen": "TCE", "evidence": "BIOACTIVATION for TCE (inverted role); GSTT1-active required for renal carcinogenesis. GSTT1-null is PROTECTIVE (OR 0.93) (PMID:20663906)", "class_name": "Chlorinated_Solvents"},
    {"source": "DCVG", "target": "DCVC", "type": "ACTIVATES",
     "carcinogen": "TCE", "evidence": "Renal brush-border processing of DCVG to DCVC", "class_name": "Chlorinated_Solvents"},
    {"source": "CCBL1", "target": "DCVC", "type": "ACTIVATES",
     "carcinogen": "TCE", "evidence": "Beta-lyase generates reactive DCVT thiol from DCVC (PMID:20663906)", "class_name": "Chlorinated_Solvents"},
    {"source": "DCVC", "target": "Renal_DNA_damage", "type": "FORMS_ADDUCT",
     "carcinogen": "TCE", "evidence": "DCVT thiol forms DNA adducts; VHL mutations in RCC (PMID:20663906)", "class_name": "Chlorinated_Solvents"},
    {"source": "PCE", "target": "TCA", "type": "ACTIVATES",
     "carcinogen": "PCE", "evidence": "PCE CYP2E1 oxidation to trichloroacetic acid", "class_name": "Chlorinated_Solvents"},
    {"source": "TCE", "target": "chlorinated_pathway", "type": "PATHWAY",
     "carcinogen": "TCE", "class_name": "Chlorinated_Solvents"},
    {"source": "PCE", "target": "chlorinated_pathway", "type": "PATHWAY",
     "carcinogen": "PCE", "class_name": "Chlorinated_Solvents"},

    # ── Class 15: Alkylating Carcinogens ──────────────────────────────────
    # Acrylamide / Glycidamide branch
    {"source": "Acrylamide", "target": "Glycidamide", "type": "ACTIVATES",
     "carcinogen": "Acrylamide", "evidence": "CYP2E1-dependent epoxidation to glycidamide, the ultimate mutagen; demonstrated in rodent and human hepatic microsomes (Segerbaeck et al. 1995, PMID:7614537)", "class_name": "Alkylating_Agents"},
    {"source": "CYP2E1", "target": "Glycidamide", "type": "ACTIVATES",
     "carcinogen": "Acrylamide", "evidence": "Primary P450 for acrylamide epoxidation (PMID:17872912)", "class_name": "Alkylating_Agents"},
    {"source": "Glycidamide", "target": "N7_GA_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "Glycidamide", "evidence": "Direct alkylation at N7-guanine; principal acrylamide-related adduct in rodent and human tissues (PMID:17872912)", "class_name": "Alkylating_Agents"},
    {"source": "GSTP1", "target": "Glycidamide", "type": "DETOXIFIES",
     "carcinogen": "Glycidamide", "evidence": "GSH conjugation of the epoxide; GSTP1/M1-null genotype elevates adduct burden (Besaratinia, Mutat Res 2010)", "class_name": "Alkylating_Agents"},
    {"source": "GSTM1", "target": "Glycidamide", "type": "DETOXIFIES",
     "carcinogen": "Glycidamide", "evidence": "Secondary GSH conjugation pathway (Besaratinia, Mutat Res 2010)", "class_name": "Alkylating_Agents"},
    {"source": "XRCC1", "target": "N7_GA_dG", "type": "REPAIRS",
     "carcinogen": "Glycidamide", "evidence": "BER scaffold response to depurinated N7-glycidamide adducts (PMID:17872912)", "class_name": "Alkylating_Agents"},
    {"source": "Acrylamide", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "Acrylamide", "class_name": "Alkylating_Agents"},
    {"source": "Glycidamide", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "Glycidamide", "class_name": "Alkylating_Agents"},

    # Cyclophosphamide branch
    {"source": "Cyclophosphamide", "target": "4OH_cyclophosphamide", "type": "ACTIVATES",
     "carcinogen": "Cyclophosphamide", "evidence": "CYP-catalyzed 4-hydroxylation is the obligatory prodrug activation step (IARC Monograph 100A)", "class_name": "Alkylating_Agents"},
    {"source": "CYP2B6", "target": "4OH_cyclophosphamide", "type": "ACTIVATES",
     "carcinogen": "Cyclophosphamide", "evidence": "Major 4-hydroxylation activity; CYP2B6*6 (rs3745274) reduces clearance (PMID:10220571)", "class_name": "Alkylating_Agents"},
    {"source": "CYP3A4", "target": "4OH_cyclophosphamide", "type": "ACTIVATES",
     "carcinogen": "Cyclophosphamide", "evidence": "Secondary 4-hydroxylation with competing N-dechloroethylation producing chloroacetaldehyde (PMID:10220571)", "class_name": "Alkylating_Agents"},
    {"source": "4OH_cyclophosphamide", "target": "Phosphoramide_mustard", "type": "ACTIVATES",
     "carcinogen": "Cyclophosphamide", "evidence": "Spontaneous tautomerization to aldophosphamide then beta-elimination to phosphoramide mustard plus acrolein (PMID:10220571)", "class_name": "Alkylating_Agents"},
    {"source": "Phosphoramide_mustard", "target": "DNA_ICL_mustard", "type": "FORMS_ADDUCT",
     "carcinogen": "Cyclophosphamide", "evidence": "Bifunctional aziridinium forms G-G interstrand crosslinks (PMID:10220571)", "class_name": "Alkylating_Agents"},
    {"source": "GSTP1", "target": "Phosphoramide_mustard", "type": "DETOXIFIES",
     "carcinogen": "Cyclophosphamide", "evidence": "GSH conjugation of aziridinium intermediates; GSTP1 Ile105Val (rs1695) modulates cyclophosphamide sensitivity (Yang et al., Pharmacogenet Genomics 2005)", "class_name": "Alkylating_Agents"},
    {"source": "Cyclophosphamide", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "Cyclophosphamide", "class_name": "Alkylating_Agents"},

    # Chlorambucil branch
    {"source": "Chlorambucil", "target": "Chlorambucil_aziridinium", "type": "ACTIVATES",
     "carcinogen": "Chlorambucil", "evidence": "Intramolecular cyclization displaces chloride to form the cyclic aziridinium ultimate alkylator (PMID:7523912)", "class_name": "Alkylating_Agents"},
    {"source": "Chlorambucil_aziridinium", "target": "DNA_ICL_mustard", "type": "FORMS_ADDUCT",
     "carcinogen": "Chlorambucil", "evidence": "Bifunctional N7-dG interstrand crosslinks drive therapeutic cytotoxicity and therapy-related AML risk (PMID:7523912; IARC 100A)", "class_name": "Alkylating_Agents"},
    {"source": "GSTP1", "target": "Chlorambucil_aziridinium", "type": "DETOXIFIES",
     "carcinogen": "Chlorambucil", "evidence": "GSH conjugation; GSTP1 Ile105Val modulates chlorambucil response in CLL (Meier et al., Clin Cancer Res 2014)", "class_name": "Alkylating_Agents"},
    {"source": "Chlorambucil", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "Chlorambucil", "class_name": "Alkylating_Agents"},

    # Sulfur mustard branch
    {"source": "Sulfur_mustard", "target": "Mustard_episulfonium", "type": "ACTIVATES",
     "carcinogen": "Sulfur_mustard", "evidence": "Intramolecular cyclization to cyclic episulfonium is the rate-limiting step of alkylation (PMID:8635461)", "class_name": "Alkylating_Agents"},
    {"source": "Mustard_episulfonium", "target": "DNA_ICL_mustard", "type": "FORMS_ADDUCT",
     "carcinogen": "Sulfur_mustard", "evidence": "Bifunctional alkylation produces N7-Gua monoadducts and 5'-d(GNC) interstrand crosslinks (PMID:8635461; IARC 100F)", "class_name": "Alkylating_Agents"},
    {"source": "GSTM1", "target": "Mustard_episulfonium", "type": "DETOXIFIES",
     "carcinogen": "Sulfur_mustard", "evidence": "GSH conjugation; GSTM1-null correlates with elevated mustard-adduct burden in exposed veterans (IARC Monograph 100F)", "class_name": "Alkylating_Agents"},
    {"source": "Sulfur_mustard", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "Sulfur_mustard", "class_name": "Alkylating_Agents"},

    # Busulfan branch
    {"source": "Busulfan", "target": "Busulfan_methanesulfonate", "type": "ACTIVATES",
     "carcinogen": "Busulfan", "evidence": "Non-enzymatic hydrolytic displacement of methanesulfonate leaving groups produces the reactive bifunctional alkylator (PMID:12960109)", "class_name": "Alkylating_Agents"},
    {"source": "Busulfan_methanesulfonate", "target": "DNA_ICL_mustard", "type": "FORMS_ADDUCT",
     "carcinogen": "Busulfan", "evidence": "Bifunctional alkylation yields N7-THPG adducts and DNA-DNA interstrand crosslinks (PMID:12960109; IARC 100A)", "class_name": "Alkylating_Agents"},
    {"source": "GSTP1", "target": "Busulfan", "type": "DETOXIFIES",
     "carcinogen": "Busulfan", "evidence": "GSH-Bu conjugation; GSTA1/P1 polymorphisms alter busulfan clearance in pediatric HSCT (Ansari et al., Blood 2010)", "class_name": "Alkylating_Agents"},
    {"source": "Busulfan", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "Busulfan", "class_name": "Alkylating_Agents"},

    # MNU branch (direct methylator)
    {"source": "MNU", "target": "Methyldiazonium", "type": "ACTIVATES",
     "carcinogen": "MNU", "evidence": "Spontaneous pH-dependent decomposition to methyldiazohydroxide and methyldiazonium ion; no enzymatic activation required (PMID:2185966)", "class_name": "Alkylating_Agents"},
    {"source": "Methyldiazonium", "target": "N7_methyl_dG", "type": "FORMS_ADDUCT",
     "carcinogen": "MNU", "evidence": "Dominant guanine methylation product from MNU; repaired by BER after spontaneous depurination (PMID:2185966)", "class_name": "Alkylating_Agents"},
    {"source": "XRCC1", "target": "N7_methyl_dG", "type": "REPAIRS",
     "carcinogen": "MNU", "evidence": "BER scaffold response to abasic sites generated by depurination of N7-methyl-dG", "class_name": "Alkylating_Agents"},
    {"source": "MNU", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "MNU", "class_name": "Alkylating_Agents"},

    # Temozolomide branch
    {"source": "Temozolomide", "target": "MTIC", "type": "ACTIVATES",
     "carcinogen": "Temozolomide", "evidence": "Spontaneous pH-dependent hydrolysis of the imidazotetrazine ring to MTIC; no enzymatic requirement (PMID:9327140)", "class_name": "Alkylating_Agents"},
    {"source": "MTIC", "target": "Methyldiazonium", "type": "ACTIVATES",
     "carcinogen": "Temozolomide", "evidence": "MTIC decomposes to methyldiazonium plus 5-aminoimidazole-4-carboxamide at physiological pH (PMID:9327140)", "class_name": "Alkylating_Agents"},
    {"source": "Temozolomide", "target": "alkylating_pathway", "type": "PATHWAY",
     "carcinogen": "Temozolomide", "class_name": "Alkylating_Agents"},
]


# ── Functions ─────────────────────────────────────────────────────────────


def _all_wave2_node_dicts(
    *,
    class_name: str | None = None,
    include_auxiliary: bool = True,
) -> list[dict[str, Any]]:
    """Return Wave 2 node definition dicts, optionally filtered by class."""
    nodes = (
        WAVE2_CARCINOGEN_NODES
        + WAVE2_ENZYME_NODES
        + WAVE2_METABOLITE_NODES
        + WAVE2_DNA_ADDUCT_NODES
        + WAVE2_PATHWAY_NODES
    )
    if include_auxiliary:
        nodes = nodes + WAVE2_AUXILIARY_NODES
    if class_name is None:
        return nodes
    return [node for node in nodes if node.get("class_name") == class_name]


def _bridge_core_reference_nodes(*, class_name: str | None = None) -> list[Node]:
    """Return minimal core nodes required to make standalone Wave 2 graphs valid."""
    if class_name is None:
        ref_ids = {
            ref_id
            for refs in WAVE2_CROSS_REFERENCES.values()
            for ref_id in refs
        }
    else:
        ref_ids = set(WAVE2_CROSS_REFERENCES.get(class_name, []))
    if not ref_ids:
        return []

    from .expanded_metals import load_core_metallo_reference_graph

    core = load_core_metallo_reference_graph()
    bridge_nodes: list[Node] = []
    for node in core.nodes:
        if node.id in ref_ids:
            bridge_nodes.append(node.model_copy(deep=True))
    return bridge_nodes


def _wave2_node_list(
    *,
    class_name: str | None = None,
    include_auxiliary: bool = True,
    include_bridge_core_refs: bool = False,
) -> list[Node]:
    """Instantiate Wave 2 nodes, optionally filtered by class."""
    nodes: list[Node] = []
    for d in _all_wave2_node_dicts(
        class_name=class_name,
        include_auxiliary=include_auxiliary,
    ):
        ntype = NodeType(d["type"])
        nodes.append(
            Node(
                id=d["id"],
                label=d["label"],
                type=ntype,
                group=d.get("group"),
                iarc=d.get("iarc"),
                phase=d.get("phase"),
                role=d.get("role"),
                detail=d.get("detail", ""),
                tissue=d.get("tissue"),
                exposure=d.get("exposure"),
                reactivity=d.get("reactivity"),
                variant=d.get("variant"),
            )
        )
    if include_bridge_core_refs:
        existing_ids = {node.id for node in nodes}
        for node in _bridge_core_reference_nodes(class_name=class_name):
            if node.id not in existing_ids:
                nodes.append(node)
                existing_ids.add(node.id)
    return nodes


def _wave2_edge_dicts(*, class_name: str | None = None) -> list[dict[str, Any]]:
    """Return Wave 2 edge dicts, optionally filtered by class."""
    if class_name is None:
        return WAVE2_EDGES
    return [edge for edge in WAVE2_EDGES if edge.get("class_name") == class_name]


def _wave2_edges_for_node_ids(
    node_ids: set[str],
    *,
    class_name: str | None = None,
) -> list[Edge]:
    """Instantiate edges whose endpoints and carcinogen refs exist in ``node_ids``."""
    edges: list[Edge] = []
    for d in _wave2_edge_dicts(class_name=class_name):
        if d["source"] not in node_ids or d["target"] not in node_ids:
            continue
        carcinogen = d.get("carcinogen")
        if carcinogen and carcinogen not in node_ids:
            continue
        etype = EdgeType(d["type"])
        edges.append(
            Edge(
                source=d["source"],
                target=d["target"],
                type=etype,
                carcinogen=carcinogen,
                evidence=d.get("evidence"),
            )
        )
    return edges


def merge_wave2_classes_into(
    core_graph: KnowledgeGraph,
    *,
    class_name: str | None = None,
) -> KnowledgeGraph:
    """Merge Wave 2 nodes and edges into an existing graph.

    Uses setdefault to skip nodes that already exist in the core graph
    (e.g., NDMA, Methyldiazonium, O6_methyl_dG, CYP2E1, CYP1A1), preventing duplication while adding
    genuinely new content.
    """
    wave2_nodes = _wave2_node_list(class_name=class_name)
    node_by_id: dict[str, Node] = {n.id: n for n in core_graph.nodes}
    for node in wave2_nodes:
        node_by_id.setdefault(node.id, node)

    combined_nodes = list(node_by_id.values())
    node_ids = set(node_by_id)

    wave2_edges = _wave2_edges_for_node_ids(node_ids, class_name=class_name)
    existing_edges = {(e.source, e.target, e.type.value) for e in core_graph.edges}
    combined_edges = list(core_graph.edges)

    for edge in wave2_edges:
        edge_sig = (edge.source, edge.target, edge.type.value)
        if edge_sig not in existing_edges:
            combined_edges.append(edge)
            existing_edges.add(edge_sig)

    return KnowledgeGraph(nodes=combined_nodes, edges=combined_edges)


def get_class_profile(class_name: str) -> CarcinogenClassProfile | None:
    """Get profile for a specific Wave 2 class.

    Args:
        class_name: Class name (e.g. ``'Aldehydes'``, ``'Dioxins_AhR'``)

    Returns:
        CarcinogenClassProfile if found, None otherwise
    """
    return WAVE2_CLASS_PROFILES.get(class_name)


def get_all_wave2_classes() -> list[CarcinogenClassProfile]:
    """Get all Wave 2 class profiles."""
    return list(WAVE2_CLASS_PROFILES.values())


def get_class_specific_nodes(class_name: str) -> tuple[list[Node], list[Edge]]:
    """Get nodes and edges for a single Wave 2 class.

    Args:
        class_name: Class name (e.g. ``'Aldehydes'``, ``'Chlorinated_Solvents'``)

    Returns:
        Tuple of (nodes, edges) specific to that class
    """
    nodes = _wave2_node_list(
        class_name=class_name,
        include_bridge_core_refs=True,
    )
    node_ids = {n.id for n in nodes}

    class_edges = _wave2_edge_dicts(class_name=class_name)
    edges: list[Edge] = []
    for d in class_edges:
        if d["source"] not in node_ids or d["target"] not in node_ids:
            continue
        carcinogen = d.get("carcinogen")
        if carcinogen and carcinogen not in node_ids:
            continue
        etype = EdgeType(d["type"])
        edges.append(
            Edge(
                source=d["source"],
                target=d["target"],
                type=etype,
                carcinogen=carcinogen,
                evidence=d.get("evidence"),
            )
        )

    return nodes, edges


def build_single_class_graph(
    class_name: str,
    *,
    include_core: bool = False,
) -> KnowledgeGraph:
    """Build graph for a single Wave 2 class.

    Args:
        class_name: Class name (e.g. ``'Aldehydes'``)
        include_core: If ``True``, merge into core reference graph.
    """
    nodes, edges = get_class_specific_nodes(class_name)
    kg = KnowledgeGraph(nodes=nodes, edges=edges)
    if include_core:
        from .expanded_metals import load_core_metallo_reference_graph
        core = load_core_metallo_reference_graph()
        return merge_wave2_classes_into(core, class_name=class_name)
    return kg


def build_wave2_class_graph(*, include_core: bool = True) -> KnowledgeGraph:
    """Build all Wave 2 class pathways.

    Args:
        include_core: If ``True`` (default), merge into the core reference
            graph. If ``False``, return only Wave 2 nodes and edges.
    """
    if include_core:
        from .expanded_metals import load_core_metallo_reference_graph
        core = load_core_metallo_reference_graph()
        return merge_wave2_classes_into(core)

    nodes = _wave2_node_list(include_bridge_core_refs=True)
    node_ids = {n.id for n in nodes}
    edges = _wave2_edges_for_node_ids(node_ids)
    return KnowledgeGraph(nodes=nodes, edges=edges)
