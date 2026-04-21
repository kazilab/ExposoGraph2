"""Reusable seeded example graphs for ExposoGraph."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from .config import GraphMode, GraphVisibility
from .db_clients.iarc import IARCClassifier
from .engine import GraphEngine
from .exporter import (
    parse_graph_data_js,
    to_graph_data_js,
    to_interactive_html,
    to_json,
    to_plotly_html,
)
from .grounding import ground_node
from .models import (
    Edge,
    EdgeType,
    KnowledgeGraph,
    MatchStatus,
    Node,
    NodeType,
    ProvenanceRecord,
)
from .reference_data import ACTIVITY_SCORES, REFERENCE_KEGG_PATHWAYS, build_full_panel

Spec = dict[str, str]


@dataclass(frozen=True)
class ArchitectureInventoryGroup:
    name: str
    members: tuple[str, ...]
    count: int


@dataclass(frozen=True)
class ArchitectureSummary:
    node_count: int
    edge_count: int
    node_type_count: int
    edge_type_count: int
    node_type_counts: dict[str, int]
    edge_type_counts: dict[str, int]
    carcinogens: tuple[str, ...]
    enzymes: tuple[str, ...]
    metabolites: tuple[str, ...]
    dna_adducts: tuple[str, ...]
    pathway_labels: tuple[str, ...]
    carcinogen_classes: tuple[ArchitectureInventoryGroup, ...]
    enzyme_categories: tuple[ArchitectureInventoryGroup, ...]


_SHOWCASE_NAMESPACE = "curated_showcase"
_CURATED_SOURCE = "Curated showcase"
_DEFAULT_BUNDLE_PATH = Path("ExposoGraph/temp_kg/graph-data.js")
_FULL_LEGENDS_REFERENCE_GRAPH = (
    Path(__file__).resolve().parent / "reference_graphs" / "full_legends_graph.js"
)
_IARC = IARCClassifier()

_FULL_LEGENDS_ID_REMAP: dict[str, str] = {
    "8oxodG": "Oxo_dG",
    "AFB1_epox": "AFB1_epoxide",
    "AFB1_N7G": "AFB1_Gua",
    "BENZ": "Benzene",
    "BNZ": "Benzidine",
    "BaP_78diol": "BaP_diol",
    "BaP_78epox": "BaP_epoxide",
    "BenzO": "Benzene_oxide",
    "BQ": "Benzoquinone",
    "CEO": "Chloroethylene_oxide",
    "EO": "EthyleneOxide",
    "TESTO": "Testosterone",
    "VCM": "VinylChloride",
    "NOHPHIP": "NOH_PhIP",
    "NOH4ABP": "NOH_4ABP",
    "NOHBNZ": "NOH_Benzidine",
    "4OHE2": "HydroxyE2",
    "A4": "Androstenedione",
    "6bOHT": "HydroxyTestosterone",
    "T_gluc": "Testosterone_gluc",
    "NDMA_alpha": "NDMA_hydroxyl",
    "NNK_alpha": "NNK_hydroxyl",
    "O6MeG": "O6_methyl_dG",
    "KEGG980": "hsa00980",
    "KEGG982": "hsa00982",
    "KEGG5208": "hsa05208",
    "KEGG5204": "hsa05204",
    "KEGG_ANDRO": "hsa00140",
}

_PHASE_I_HORMONE_GENES = ("CYP17A1", "SRD5A1", "SRD5A2", "CYP19A1", "AKR1C3")
_PHASE_II_HORMONE_GENES = ("COMT", "SULT1E1", "UGT2B17", "UGT2B7")
_CARCINOGEN_CLASS_TITLES: tuple[tuple[str, str], ...] = (
    ("PAH", "PAH"),
    ("HCA", "HCA"),
    ("Aromatic_Amine", "Aromatic Amines"),
    ("Nitrosamine", "Nitrosamines"),
    ("Mycotoxin", "Mycotoxins"),
    ("Estrogen", "Estrogens"),
    ("Androgen", "Androgens"),
    ("Solvent", "Solvents"),
    ("Alkylating", "Alkylating Agents"),
)


def _entity_spec(
    node_id: str,
    label: str,
    detail: str,
    *,
    group: str | None = None,
    reactivity: str | None = None,
) -> Spec:
    spec = {
        "id": node_id,
        "label": label,
        "detail": detail,
    }
    if group is not None:
        spec["group"] = group
    if reactivity is not None:
        spec["reactivity"] = reactivity
    return spec


def _relation_spec(
    source: str,
    target: str,
    edge_type: str,
    carcinogen: str | None,
    label: str,
    *,
    custom_predicate: str | None = None,
) -> Spec:
    spec = {
        "source": source,
        "target": target,
        "type": edge_type,
        "label": label,
    }
    if carcinogen is not None:
        spec["carcinogen"] = carcinogen
    if custom_predicate is not None:
        spec["custom_predicate"] = custom_predicate
    return spec


_CARCINOGEN_SPECS: tuple[Spec, ...] = (
    _entity_spec(
        "BaP",
        "Benzo[a]pyrene",
        "Reference PAH carcinogen used for the canonical activation-to-adduct chain.",
        group="PAH",
    ),
    _entity_spec(
        "DMBA",
        "DMBA",
        (
            "Reference PAH-like polycyclic carcinogen."
        ),
        group="PAH",
    ),
    _entity_spec(
        "PhIP",
        "PhIP",
        "Representative heterocyclic amine produced during high-temperature cooking.",
        group="HCA",
    ),
    _entity_spec(
        "MeIQx",
        "MeIQx",
        (
            "Representative heterocyclic amine activated through "
            "N-hydroxylation and esterification."
        ),
        group="HCA",
    ),
    _entity_spec(
        "4ABP",
        "4-Aminobiphenyl",
        (
            "Representative aromatic amine requiring N-hydroxylation and "
            "acetyl-transfer activation."
        ),
        group="Aromatic_Amine",
    ),
    _entity_spec(
        "Benzidine",
        "Benzidine",
        "Representative aromatic amine included for class-level alignment.",
        group="Aromatic_Amine",
    ),
    _entity_spec(
        "NNK",
        "NNK",
        (
            "Tobacco-specific nitrosamine represented as a "
            "pyridyloxobutylating carcinogen context."
        ),
        group="Nitrosamine",
    ),
    _entity_spec(
        "NDMA",
        "NDMA",
        "Small-molecule nitrosamine represented via methylating chemistry.",
        group="Nitrosamine",
    ),
    _entity_spec(
        "AFB1",
        "Aflatoxin B1",
        (
            "Mycotoxin represented through CYP-mediated epoxidation and "
            "guanine adduct formation."
        ),
        group="Mycotoxin",
    ),
    _entity_spec(
        "E2",
        "17beta-Estradiol",
        "Endogenous estrogen represented in the hormone-metabolism reference layer.",
        group="Estrogen",
    ),
    _entity_spec(
        "Testosterone",
        "Testosterone",
        (
            "Androgen precursor node used to anchor the steroid-hormone "
            "metabolism branch."
        ),
        group="Androgen",
    ),
    _entity_spec(
        "DHT",
        "5a-DHT",
        "Dihydrotestosterone node.",
        group="Androgen",
    ),
    _entity_spec(
        "Benzene",
        "Benzene",
        (
            "Volatile solvent carcinogen represented through CYP2E1 oxidation "
            "and ROS-linked damage."
        ),
        group="Solvent",
    ),
    _entity_spec(
        "VinylChloride",
        "Vinyl Chloride",
        "Solvent carcinogen represented through etheno-adduct chemistry.",
        group="Solvent",
    ),
    _entity_spec(
        "EthyleneOxide",
        "Ethylene Oxide",
        "Alkylating carcinogen represented through hydroxyethyl adduct formation.",
        group="Alkylating",
    ),
)

_METABOLITE_SPECS: tuple[Spec, ...] = (
    _entity_spec(
        "BaP_epoxide",
        "BaP-7,8-epoxide",
        "Initial PAH epoxide intermediate formed during BaP activation.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "BaP_diol",
        "BaP-7,8-diol",
        "Hydrolysis product of the BaP epoxide intermediate.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "BPDE",
        "BPDE",
        "Ultimate diol-epoxide carcinogen formed from BaP.",
        reactivity="High",
    ),
    _entity_spec(
        "BPDE_GSH",
        "BPDE-GSH conjugate",
        "Glutathione-conjugated detoxification product of BPDE.",
        reactivity="Low",
    ),
    _entity_spec(
        "DMBA_epoxide",
        "DMBA-3,4-epoxide",
        "Initial epoxide intermediate formed during DMBA activation.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "DMBA_diol",
        "DMBA-3,4-diol",
        "Hydrolyzed DMBA intermediate prior to diol-epoxide formation.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "DMBA_diol_epoxide",
        "DMBA-diol-epoxide",
        "Ultimate DMBA metabolite used to anchor the PAH class expansion.",
        reactivity="High",
    ),
    _entity_spec(
        "NOH_PhIP",
        "N-hydroxy-PhIP",
        "Phase I activated PhIP intermediate prior to sulfate ester formation.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "PhIP_sulfate",
        "PhIP-N3-sulfate",
        "Reactive sulfate ester representing HCA bioactivation.",
        reactivity="High",
    ),
    _entity_spec(
        "NOH_MeIQx",
        "N-hydroxy-MeIQx",
        "N-hydroxylated MeIQx activation intermediate.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "MeIQx_acetoxy",
        "MeIQx-N-acetoxy ester",
        "Reactive ester representing MeIQx activation.",
        reactivity="High",
    ),
    _entity_spec(
        "NOH_4ABP",
        "N-hydroxy-4-ABP",
        "N-hydroxylated 4-aminobiphenyl activation intermediate.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "ABP_acetoxy",
        "4-ABP-N-acetoxy ester",
        "Reactive aromatic-amine ester leading to guanine adduct formation.",
        reactivity="High",
    ),
    _entity_spec(
        "NOH_Benzidine",
        "N-hydroxy-benzidine",
        "Benzidine activation intermediate represented for class completeness.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "Benzidine_ester",
        "Benzidine-diacetoxy ester",
        "Reactive benzidine ester included as an aromatic-amine endpoint.",
        reactivity="High",
    ),
    _entity_spec(
        "NNK_pob",
        "NNK pyridyloxobutylating intermediate",
        "Representative NNK activation product leading to POB adducts.",
        reactivity="High",
    ),
    _entity_spec(
        "NNK_hydroxyl",
        "NNK alpha-hydroxylated intermediate",
        "Primary NNK activation intermediate for the nitrosamine branch.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "NDMA_hydroxyl",
        "NDMA alpha-hydroxylated intermediate",
        "Primary NDMA activation intermediate.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "Methyldiazonium",
        "Methyldiazonium ion",
        "Reactive methylating species used to represent NDMA alkylation.",
        reactivity="High",
    ),
    _entity_spec(
        "AFB1_epoxide",
        "AFB1-8,9-epoxide",
        "Canonical aflatoxin epoxide intermediate.",
        reactivity="High",
    ),
    _entity_spec(
        "HydroxyE2",
        "4-Hydroxyestradiol",
        (
            "Catechol estrogen intermediate formed during CYP1B1-mediated "
            "estrogen activation."
        ),
        reactivity="Intermediate",
    ),
    _entity_spec(
        "E2_quinone",
        "Estradiol-3,4-quinone",
        "Reactive estrogen quinone that yields depurinating DNA adducts.",
        reactivity="High",
    ),
    _entity_spec(
        "Testosterone_gluc",
        "Testosterone-17-glucuronide",
        "Representative androgen detoxification product.",
        reactivity="Low",
    ),
    _entity_spec(
        "DHT_gluc",
        "DHT-17-glucuronide",
        "Representative DHT conjugate used for androgen detoxification.",
        reactivity="Low",
    ),
    _entity_spec(
        "Benzene_oxide",
        "Benzene oxide",
        "Early benzene oxidation product linked to downstream quinone chemistry.",
        reactivity="Intermediate",
    ),
    _entity_spec(
        "Benzoquinone",
        "Benzoquinone",
        "Representative ROS-linked benzene metabolite.",
        reactivity="High",
    ),
    _entity_spec(
        "Chloroethylene_oxide",
        "Chloroethylene oxide",
        "Reactive vinyl chloride intermediate that yields etheno adducts.",
        reactivity="High",
    ),
    _entity_spec(
        "Hydroxyethylating_species",
        "2-hydroxyethylating species",
        "Reactive ethylene oxide-derived alkylating intermediate.",
        reactivity="High",
    ),
)

_DNA_ADDUCT_SPECS: tuple[Spec, ...] = (
    _entity_spec(
        "BPDE_dG",
        "BPDE-N2-dG",
        "Canonical bulky PAH DNA adduct.",
    ),
    _entity_spec(
        "DMBA_dA",
        "DMBA-dA adduct",
        "Representative DMBA-derived DNA adduct.",
    ),
    _entity_spec(
        "PhIP_dG",
        "PhIP-C8-dG",
        "Representative HCA DNA adduct formed from activated PhIP.",
    ),
    _entity_spec(
        "ABP_dG",
        "dG-C8-4-ABP",
        "Representative aromatic-amine guanine adduct.",
    ),
    _entity_spec(
        "POB_dG",
        "O6-POB-dG",
        "Representative NNK-derived pyridyloxobutyl DNA adduct.",
    ),
    _entity_spec(
        "O6_methyl_dG",
        "O6-methyl-dG",
        "Representative NDMA-derived methyl DNA adduct.",
    ),
    _entity_spec(
        "AFB1_Gua",
        "AFB1-N7-Gua",
        "Canonical aflatoxin-derived guanine adduct.",
    ),
    _entity_spec(
        "E2_Ade",
        "4-OHE2-1-N3Ade",
        "Representative depurinating estrogen-DNA adduct.",
    ),
    _entity_spec(
        "Oxo_dG",
        "8-oxo-dG",
        "Representative oxidative DNA lesion linked to ROS pathway activity.",
    ),
    _entity_spec(
        "EthenoG",
        "N2,3-ethenoguanine",
        "Representative etheno adduct associated with vinyl chloride exposure.",
    ),
    _entity_spec(
        "HEG",
        "N7-(2-hydroxyethyl)guanine",
        "Representative hydroxyethyl DNA adduct associated with ethylene oxide.",
    ),
)

_BASE_CARCINOGEN_SPECS = {spec["id"]: spec for spec in _CARCINOGEN_SPECS}
_BASE_METABOLITE_SPECS = {spec["id"]: spec for spec in _METABOLITE_SPECS}
_BASE_DNA_ADDUCT_SPECS = {spec["id"]: spec for spec in _DNA_ADDUCT_SPECS}

_ANDROGEN_TISSUE_SPECS: tuple[Spec, ...] = (
    _entity_spec(
        "Prostate",
        "Prostate tissue",
        "Primary androgen-responsive tissue represented for the prostate-focused module.",
    ),
    _entity_spec(
        "PeripheralHormoneTissues",
        "Peripheral hormone-metabolism tissues",
        (
            "Peripheral tissues where local androgen conversion contributes "
            "to hormone exposure context."
        ),
    ),
    _entity_spec(
        "HormoneResponsiveTissues",
        "Hormone-responsive tissues",
        (
            "Breast, adipose, and related tissues where aromatase-linked "
            "estrogen formation is relevant."
        ),
    ),
)

_ANDROGEN_ENZYME_SPECS: tuple[Spec, ...] = (
    _entity_spec(
        "CYP3A5",
        "CYP3A5",
        (
            "Testosterone 6beta-hydroxylation enzyme with prostate-relevant "
            "androgen-clearance and androgen-receptor modulation context."
        ),
    ),
)

_ANDROGEN_GENE_SPECS: tuple[Spec, ...] = (
    {
        "id": "AR",
        "label": "Androgen receptor (AR)",
        "detail": (
            "Nuclear androgen receptor represented as the signaling endpoint "
            "for testosterone and DHT binding in prostate-linked proliferation."
        ),
        "tissue": "prostate, hormone-responsive tissues",
        "source_db": _CURATED_SOURCE,
    },
    {
        "id": "SRD5A2_V89L",
        "label": "SRD5A2 V89L",
        "detail": (
            "Variant annotation for the common SRD5A2 V89L polymorphism used to "
            "track reduced DHT-forming activity in tooltip and detail views."
        ),
        "variant": "V89L",
        "phenotype": "Reduced 5alpha-reductase activity",
        "tissue": "prostate",
        "source_db": "PubMed",
    },
    {
        "id": "SRD5A2_A49T",
        "label": "SRD5A2 A49T",
        "detail": (
            "Variant annotation for the SRD5A2 A49T polymorphism used to track "
            "increased DHT-forming activity in the androgen module."
        ),
        "variant": "A49T",
        "phenotype": "Increased 5alpha-reductase activity",
        "tissue": "prostate",
        "source_db": "PubMed",
    },
    {
        "id": "CYP19A1_repeat_length",
        "label": "CYP19A1 repeat-length variant",
        "detail": (
            "Curated placeholder for CYP19A1 repeat-length polymorphisms "
            "associated with altered aromatase expression."
        ),
        "variant": "Repeat-length polymorphism",
        "phenotype": "Expression-dependent aromatase activity shift",
        "tissue": "hormone-responsive tissues",
        "source_db": "PubMed",
    },
    {
        "id": "UGT2B17_copy_number_deletion",
        "label": "UGT2B17 copy number deletion",
        "detail": (
            "Copy-number deletion annotation used to capture loss of UGT2B17-"
            "mediated androgen glucuronidation."
        ),
        "variant": "Copy number deletion",
        "phenotype": "Absent or markedly reduced androgen glucuronidation",
        "tissue": "prostate, liver",
        "source_db": "PubMed",
    },
)

_ANDROGEN_METABOLITE_SPECS: tuple[Spec, ...] = (
    _entity_spec(
        "Androstenedione",
        "Androstenedione",
        (
            "Steroid intermediate represented as the immediate AKR1C3 "
            "substrate for testosterone formation."
        ),
        reactivity="Low",
    ),
    _entity_spec(
        "HydroxyTestosterone",
        "6beta-Hydroxytestosterone",
        "Representative CYP3A5-linked hydroxylated androgen metabolite.",
        reactivity="Low",
    ),
    _BASE_METABOLITE_SPECS["HydroxyE2"],
    _BASE_METABOLITE_SPECS["E2_quinone"],
    _BASE_METABOLITE_SPECS["Testosterone_gluc"],
    _BASE_METABOLITE_SPECS["DHT_gluc"],
)

_ANDROGEN_DNA_ADDUCT_SPECS: tuple[Spec, ...] = (
    _BASE_DNA_ADDUCT_SPECS["E2_Ade"],
    _entity_spec(
        "E2_Gua",
        "4-OHE2-guanine depurinating adduct",
        "Representative depurinating estrogen-DNA adduct at guanine residues.",
    ),
)

_ANDROGEN_PATHWAY_SPECS: tuple[Spec, ...] = (
    {
        "pathway_id": "hsa00140",
        "label": "Steroid hormone biosynthesis",
        "focus": (
            "Hormone-metabolism context for androgen synthesis, DHT "
            "conversion, and aromatase bridging."
        ),
    },
    {
        "pathway_id": "hsa05204",
        "label": "Chemical carcinogenesis - DNA adducts",
        "focus": (
            "Estrogen-quinone DNA-adduct context used to bridge the androgen "
            "and estrogen modules."
        ),
    },
    {
        "pathway_id": "AR_signal_program",
        "label": "AR proliferative transcriptional program",
        "focus": "Curated signaling endpoint for androgen receptor-mediated proliferative drive.",
    },
)

_ANDROGEN_CARCINOGEN_IDS: tuple[str, ...] = ("Testosterone", "DHT", "E2")
_ANDROGEN_PANEL_ENZYME_IDS: tuple[str, ...] = (
    "CYP17A1",
    "SRD5A1",
    "SRD5A2",
    "CYP19A1",
    "AKR1C3",
    "CYP1B1",
    "UGT2B7",
    "UGT2B17",
)

_PATHWAY_MEMBERSHIP: dict[str, tuple[str, ...]] = {
    "hsa00980": (
        "BaP",
        "BPDE",
        "BPDE_GSH",
        "PhIP",
        "PhIP_sulfate",
        "AFB1",
        "AFB1_epoxide",
        "Benzene",
        "CYP1A1",
        "GSTM1",
    ),
    "hsa00140": (
        "E2",
        "Testosterone",
        "DHT",
        "HydroxyE2",
        "E2_quinone",
        "Testosterone_gluc",
        "DHT_gluc",
        "CYP17A1",
        "SRD5A2",
        "CYP19A1",
    ),
    "hsa05204": (
        "BaP",
        "BPDE",
        "BPDE_dG",
        "PhIP_dG",
        "O6_methyl_dG",
        "AFB1_Gua",
        "EthenoG",
        "MGMT",
    ),
    "hsa05208": (
        "Benzene",
        "Benzoquinone",
        "Oxo_dG",
        "E2",
        "COMT",
        "OGG1",
    ),
}

_EDGE_SPECS: tuple[Spec, ...] = (
    # Procarcinogen → reactive-intermediate edges. Without these, the carcinogen
    # nodes themselves appear as isolated singletons in the D3 viewer because
    # the metabolic chain is wired enzyme → intermediate (with `carcinogen=` as
    # metadata only). These edges anchor the parent carcinogen to its chain.
    _relation_spec(
        "DMBA",
        "DMBA_epoxide",
        "ACTIVATES",
        "DMBA",
        "Procarcinogen to 3,4-epoxide",
    ),
    _relation_spec(
        "MeIQx",
        "NOH_MeIQx",
        "ACTIVATES",
        "MeIQx",
        "Procarcinogen to N-hydroxylation",
    ),
    _relation_spec(
        "Benzidine",
        "NOH_Benzidine",
        "ACTIVATES",
        "Benzidine",
        "Procarcinogen to N-hydroxylation",
    ),
    _relation_spec(
        "NDMA",
        "NDMA_hydroxyl",
        "ACTIVATES",
        "NDMA",
        "Procarcinogen to alpha-hydroxylation",
    ),
    _relation_spec(
        "VinylChloride",
        "Chloroethylene_oxide",
        "ACTIVATES",
        "VinylChloride",
        "Procarcinogen to epoxidation",
    ),
    _relation_spec(
        "E2",
        "HydroxyE2",
        "ACTIVATES",
        "E2",
        "Procarcinogen to catechol estrogen",
    ),
    _relation_spec("CYP1A1", "BaP_epoxide", "ACTIVATES", "BaP", "BaP epoxidation"),
    _relation_spec(
        "CYP1A2",
        "BaP_epoxide",
        "ACTIVATES",
        "BaP",
        "Secondary BaP epoxidation",
    ),
    _relation_spec(
        "CYP1B1",
        "BaP_epoxide",
        "ACTIVATES",
        "BaP",
        "Extrahepatic BaP epoxidation",
    ),
    _relation_spec(
        "EPHX1",
        "BaP_diol",
        "DETOXIFIES",
        "BaP",
        "Epoxide hydrolysis",
    ),
    _relation_spec(
        "CYP1A1",
        "BPDE",
        "ACTIVATES",
        "BaP",
        "Second epoxidation to BPDE",
    ),
    _relation_spec(
        "CYP1B1",
        "BPDE",
        "ACTIVATES",
        "BaP",
        "Second epoxidation to BPDE",
    ),
    _relation_spec(
        "GSTM1",
        "BPDE_GSH",
        "DETOXIFIES",
        "BaP",
        "Glutathione conjugation",
    ),
    _relation_spec(
        "GSTP1",
        "BPDE_GSH",
        "DETOXIFIES",
        "BaP",
        "Glutathione conjugation",
    ),
    _relation_spec(
        "ABCB1",
        "BPDE_GSH",
        "TRANSPORTS",
        "BaP",
        "Conjugate efflux",
    ),
    _relation_spec(
        "BPDE",
        "BPDE_dG",
        "FORMS_ADDUCT",
        "BaP",
        "Bulky guanine adduct formation",
    ),
    _relation_spec(
        "XPC",
        "BPDE_dG",
        "REPAIRS",
        "BaP",
        "Damage recognition",
    ),
    _relation_spec(
        "ERCC2",
        "BPDE_dG",
        "REPAIRS",
        "BaP",
        "NER helicase repair",
    ),
    _relation_spec(
        "CYP1A1",
        "DMBA_epoxide",
        "ACTIVATES",
        "DMBA",
        "DMBA epoxidation",
    ),
    _relation_spec(
        "CYP1B1",
        "DMBA_epoxide",
        "ACTIVATES",
        "DMBA",
        "DMBA epoxidation",
    ),
    _relation_spec(
        "EPHX1",
        "DMBA_diol",
        "DETOXIFIES",
        "DMBA",
        "Epoxide hydrolysis",
    ),
    _relation_spec(
        "CYP1A1",
        "DMBA_diol_epoxide",
        "ACTIVATES",
        "DMBA",
        "Diol-epoxide formation",
    ),
    _relation_spec(
        "CYP1B1",
        "DMBA_diol_epoxide",
        "ACTIVATES",
        "DMBA",
        "Diol-epoxide formation",
    ),
    _relation_spec(
        "DMBA_diol_epoxide",
        "DMBA_dA",
        "FORMS_ADDUCT",
        "DMBA",
        "DNA adduct formation",
    ),
    _relation_spec(
        "XPC",
        "DMBA_dA",
        "REPAIRS",
        "DMBA",
        "Bulky adduct recognition",
    ),
    _relation_spec(
        "ERCC2",
        "DMBA_dA",
        "REPAIRS",
        "DMBA",
        "Bulky adduct excision",
    ),
    _relation_spec(
        "CYP1A2",
        "NOH_PhIP",
        "ACTIVATES",
        "PhIP",
        "N-hydroxylation",
    ),
    _relation_spec(
        "SULT1A1",
        "PhIP_sulfate",
        "ACTIVATES",
        "PhIP",
        "Reactive sulfate ester formation",
    ),
    _relation_spec(
        "ABCG2",
        "PhIP_sulfate",
        "TRANSPORTS",
        "PhIP",
        "Efflux of conjugated HCA metabolite",
    ),
    _relation_spec(
        "PhIP_sulfate",
        "PhIP_dG",
        "FORMS_ADDUCT",
        "PhIP",
        "C8-dG adduct formation",
    ),
    _relation_spec(
        "XPC",
        "PhIP_dG",
        "REPAIRS",
        "PhIP",
        "Damage recognition",
    ),
    _relation_spec("ERCC2", "PhIP_dG", "REPAIRS", "PhIP", "NER excision"),
    _relation_spec(
        "CYP1A2",
        "NOH_MeIQx",
        "ACTIVATES",
        "MeIQx",
        "N-hydroxylation",
    ),
    _relation_spec(
        "NAT1",
        "MeIQx_acetoxy",
        "ACTIVATES",
        "MeIQx",
        "Reactive acetoxy ester formation",
    ),
    _relation_spec(
        "MeIQx_acetoxy",
        "PhIP_dG",
        "FORMS_ADDUCT",
        "MeIQx",
        "Representative HCA adduct formation",
    ),
    _relation_spec(
        "CYP1A2",
        "NOH_4ABP",
        "ACTIVATES",
        "4ABP",
        "N-hydroxylation",
    ),
    _relation_spec(
        "NAT2",
        "ABP_acetoxy",
        "ACTIVATES",
        "4ABP",
        "Reactive acetoxy ester formation",
    ),
    _relation_spec(
        "ABP_acetoxy",
        "ABP_dG",
        "FORMS_ADDUCT",
        "4ABP",
        "Guanine adduct formation",
    ),
    _relation_spec(
        "XPC",
        "ABP_dG",
        "REPAIRS",
        "4ABP",
        "Damage recognition",
    ),
    _relation_spec("ERCC2", "ABP_dG", "REPAIRS", "4ABP", "NER excision"),
    _relation_spec(
        "CYP1A2",
        "NOH_Benzidine",
        "ACTIVATES",
        "Benzidine",
        "Benzidine N-hydroxylation",
    ),
    _relation_spec(
        "NAT1",
        "Benzidine_ester",
        "ACTIVATES",
        "Benzidine",
        "Reactive ester formation",
    ),
    _relation_spec(
        "CYP2A6",
        "NNK_hydroxyl",
        "ACTIVATES",
        "NNK",
        "Alpha-hydroxylation",
    ),
    _relation_spec(
        "CYP2A13",
        "NNK_hydroxyl",
        "ACTIVATES",
        "NNK",
        "Pulmonary alpha-hydroxylation",
    ),
    _relation_spec(
        "CYP2D6",
        "NNK_pob",
        "ACTIVATES",
        "NNK",
        "Pyridyloxobutylating intermediate formation",
    ),
    _relation_spec(
        "NNK_pob",
        "POB_dG",
        "FORMS_ADDUCT",
        "NNK",
        "POB adduct formation",
    ),
    _relation_spec(
        "MGMT",
        "POB_dG",
        "REPAIRS",
        "NNK",
        "Direct reversal",
    ),
    _relation_spec(
        "CYP2E1",
        "NDMA_hydroxyl",
        "ACTIVATES",
        "NDMA",
        "Alpha-hydroxylation",
    ),
    _relation_spec(
        "CYP2A6",
        "Methyldiazonium",
        "ACTIVATES",
        "NDMA",
        "Secondary methylating route",
    ),
    _relation_spec(
        "Methyldiazonium",
        "O6_methyl_dG",
        "FORMS_ADDUCT",
        "NDMA",
        "Methyl adduct formation",
    ),
    _relation_spec(
        "MGMT",
        "O6_methyl_dG",
        "REPAIRS",
        "NDMA",
        "Direct reversal",
    ),
    _relation_spec(
        "CYP3A4",
        "AFB1_epoxide",
        "ACTIVATES",
        "AFB1",
        "Aflatoxin epoxidation",
    ),
    _relation_spec(
        "CYP1A2",
        "AFB1_epoxide",
        "ACTIVATES",
        "AFB1",
        "Secondary epoxidation",
    ),
    _relation_spec(
        "AFB1_epoxide",
        "AFB1_Gua",
        "FORMS_ADDUCT",
        "AFB1",
        "AFB1-N7-Gua formation",
    ),
    _relation_spec(
        "XPC",
        "AFB1_Gua",
        "REPAIRS",
        "AFB1",
        "Damage recognition",
    ),
    _relation_spec("ERCC2", "AFB1_Gua", "REPAIRS", "AFB1", "NER excision"),
    _relation_spec(
        "CYP1B1",
        "HydroxyE2",
        "ACTIVATES",
        "E2",
        "Catechol estrogen formation",
    ),
    _relation_spec(
        "CYP1B1",
        "E2_quinone",
        "ACTIVATES",
        "E2",
        "Estrogen quinone formation",
    ),
    _relation_spec(
        "UGT2B7",
        "Testosterone_gluc",
        "DETOXIFIES",
        "Testosterone",
        "Glucuronidation",
    ),
    _relation_spec(
        "UGT2B17",
        "DHT_gluc",
        "DETOXIFIES",
        "DHT",
        "Glucuronidation",
    ),
    _relation_spec(
        "ABCC2",
        "Testosterone_gluc",
        "TRANSPORTS",
        "Testosterone",
        "Conjugate efflux",
    ),
    _relation_spec(
        "ABCC2",
        "DHT_gluc",
        "TRANSPORTS",
        "DHT",
        "Conjugate efflux",
    ),
    _relation_spec(
        "E2_quinone",
        "E2_Ade",
        "FORMS_ADDUCT",
        "E2",
        "Depurinating adduct formation",
    ),
    _relation_spec(
        "XRCC1",
        "E2_Ade",
        "REPAIRS",
        "E2",
        "BER scaffold response",
    ),
    _relation_spec(
        "CYP2E1",
        "Benzene_oxide",
        "ACTIVATES",
        "Benzene",
        "Benzene oxidation",
    ),
    _relation_spec(
        "CYP2E1",
        "Benzoquinone",
        "ACTIVATES",
        "Benzene",
        "Quinone-forming oxidation",
    ),
    _relation_spec(
        "Benzoquinone",
        "Oxo_dG",
        "FORMS_ADDUCT",
        "Benzene",
        "ROS-linked oxidative lesion",
    ),
    _relation_spec(
        "OGG1",
        "Oxo_dG",
        "REPAIRS",
        "Benzene",
        "Base excision repair",
    ),
    _relation_spec(
        "XRCC1",
        "Oxo_dG",
        "REPAIRS",
        "Benzene",
        "BER scaffold response",
    ),
    _relation_spec(
        "CYP2E1",
        "Chloroethylene_oxide",
        "ACTIVATES",
        "VinylChloride",
        "Vinyl chloride oxidation",
    ),
    _relation_spec(
        "Chloroethylene_oxide",
        "EthenoG",
        "FORMS_ADDUCT",
        "VinylChloride",
        "Etheno adduct formation",
    ),
    _relation_spec(
        "ERCC2",
        "EthenoG",
        "REPAIRS",
        "VinylChloride",
        "Damage response",
    ),
    _relation_spec(
        "Hydroxyethylating_species",
        "HEG",
        "FORMS_ADDUCT",
        "EthyleneOxide",
        "Hydroxyethyl adduct formation",
    ),
    _relation_spec(
        "XRCC1",
        "HEG",
        "REPAIRS",
        "EthyleneOxide",
        "BER scaffold response",
    ),
)

_ANDROGEN_EDGE_SPECS: tuple[Spec, ...] = (
    _relation_spec(
        "CYP17A1",
        "Androstenedione",
        "ACTIVATES",
        "Testosterone",
        "Androgen precursor synthesis",
    ),
    _relation_spec(
        "AKR1C3",
        "Testosterone",
        "CUSTOM",
        "Testosterone",
        "Androstenedione reduction to testosterone",
        custom_predicate="REDUCES_TO_TESTOSTERONE",
    ),
    _relation_spec(
        "SRD5A1",
        "DHT",
        "CUSTOM",
        "Testosterone",
        "Peripheral testosterone reduction to DHT",
        custom_predicate="CONVERTS_TO_DHT",
    ),
    _relation_spec(
        "SRD5A2",
        "DHT",
        "CUSTOM",
        "Testosterone",
        "Prostate testosterone reduction to DHT",
        custom_predicate="CONVERTS_TO_DHT",
    ),
    _relation_spec(
        "CYP19A1",
        "E2",
        "CUSTOM",
        "Testosterone",
        "Aromatization to estradiol",
        custom_predicate="AROMATIZES_TO_ESTRADIOL",
    ),
    _relation_spec(
        "CYP1B1",
        "HydroxyE2",
        "ACTIVATES",
        "E2",
        "Catechol estrogen formation",
    ),
    _relation_spec(
        "CYP1B1",
        "E2_quinone",
        "ACTIVATES",
        "E2",
        "Estrogen quinone formation",
    ),
    _relation_spec(
        "E2_quinone",
        "E2_Ade",
        "FORMS_ADDUCT",
        "E2",
        "Depurinating adenine adduct formation",
    ),
    _relation_spec(
        "E2_quinone",
        "E2_Gua",
        "FORMS_ADDUCT",
        "E2",
        "Depurinating guanine adduct formation",
    ),
    _relation_spec(
        "CYP3A5",
        "HydroxyTestosterone",
        "DETOXIFIES",
        "Testosterone",
        "6beta-hydroxylation",
    ),
    _relation_spec(
        "UGT2B7",
        "Testosterone_gluc",
        "DETOXIFIES",
        "Testosterone",
        "Glucuronidation",
    ),
    _relation_spec(
        "UGT2B17",
        "DHT_gluc",
        "DETOXIFIES",
        "DHT",
        "Glucuronidation",
    ),
    _relation_spec(
        "Testosterone",
        "AR",
        "CUSTOM",
        "Testosterone",
        "Lower-affinity androgen receptor binding",
        custom_predicate="BINDS_RECEPTOR",
    ),
    _relation_spec(
        "DHT",
        "AR",
        "CUSTOM",
        "DHT",
        "High-affinity androgen receptor binding",
        custom_predicate="BINDS_RECEPTOR",
    ),
    _relation_spec(
        "AR",
        "AR_signal_program",
        "CUSTOM",
        None,
        "AR-driven proliferative transcription program",
        custom_predicate="ACTIVATES_TRANSCRIPTION",
    ),
    _relation_spec(
        "SRD5A2",
        "Prostate",
        "EXPRESSED_IN",
        None,
        "Prostate-enriched expression context",
    ),
    _relation_spec(
        "AKR1C3",
        "Prostate",
        "EXPRESSED_IN",
        None,
        "Prostate-local androgen activation context",
    ),
    _relation_spec(
        "UGT2B17",
        "Prostate",
        "EXPRESSED_IN",
        None,
        "Prostate androgen-clearance context",
    ),
    _relation_spec(
        "CYP3A5",
        "Prostate",
        "EXPRESSED_IN",
        None,
        "Prostate androgen-clearance context",
    ),
    _relation_spec(
        "SRD5A1",
        "PeripheralHormoneTissues",
        "EXPRESSED_IN",
        None,
        "Peripheral testosterone-conversion context",
    ),
    _relation_spec(
        "CYP19A1",
        "HormoneResponsiveTissues",
        "EXPRESSED_IN",
        None,
        "Hormone-responsive aromatase context",
    ),
    _relation_spec(
        "AR",
        "Prostate",
        "EXPRESSED_IN",
        None,
        "Prostate androgen-receptor signaling context",
    ),
    _relation_spec(
        "SRD5A2_V89L",
        "SRD5A2",
        "ENCODES",
        None,
        "Variant-bearing SRD5A2 locus",
    ),
    _relation_spec(
        "SRD5A2_A49T",
        "SRD5A2",
        "ENCODES",
        None,
        "Variant-bearing SRD5A2 locus",
    ),
    _relation_spec(
        "CYP19A1_repeat_length",
        "CYP19A1",
        "ENCODES",
        None,
        "Variant-bearing CYP19A1 locus",
    ),
    _relation_spec(
        "UGT2B17_copy_number_deletion",
        "UGT2B17",
        "ENCODES",
        None,
        "Variant-bearing UGT2B17 locus",
    ),
    _relation_spec(
        "Testosterone",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "DHT",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "E2",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "CYP17A1",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "SRD5A1",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "SRD5A2",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "CYP19A1",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "AKR1C3",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "CYP3A5",
        "hsa00140",
        "PATHWAY",
        None,
        "Pathway membership: Steroid hormone biosynthesis",
    ),
    _relation_spec(
        "HydroxyE2",
        "hsa05204",
        "PATHWAY",
        None,
        "Pathway membership: Chemical carcinogenesis - DNA adducts",
    ),
    _relation_spec(
        "E2_quinone",
        "hsa05204",
        "PATHWAY",
        None,
        "Pathway membership: Chemical carcinogenesis - DNA adducts",
    ),
    _relation_spec(
        "E2_Ade",
        "hsa05204",
        "PATHWAY",
        None,
        "Pathway membership: Chemical carcinogenesis - DNA adducts",
    ),
    _relation_spec(
        "E2_Gua",
        "hsa05204",
        "PATHWAY",
        None,
        "Pathway membership: Chemical carcinogenesis - DNA adducts",
    ),
    _relation_spec(
        "DHT",
        "AR_signal_program",
        "PATHWAY",
        None,
        "Pathway membership: AR proliferative transcriptional program",
    ),
    _relation_spec(
        "AR",
        "AR_signal_program",
        "PATHWAY",
        None,
        "Pathway membership: AR proliferative transcriptional program",
    ),
)


def _curated_ref(
    record_id: str,
    citation: str,
    evidence: str,
    *,
    source_db: str = _CURATED_SOURCE,
) -> ProvenanceRecord:
    return ProvenanceRecord(
        source_db=source_db,
        record_id=record_id,
        citation=citation,
        evidence=evidence,
    )


def _kegg_ref(pathway_id: str, label: str) -> ProvenanceRecord:
    return ProvenanceRecord(
        source_db="KEGG",
        record_id=pathway_id,
        citation=f"KEGG pathway record for {label}",
        url=f"https://www.kegg.jp/entry/{pathway_id}",
        evidence="Curated KEGG reference pathway tracked by ExposoGraph.",
    )


def _representative_activity_score(gene_id: str) -> float | None:
    rows = ACTIVITY_SCORES.get(gene_id) or []
    numeric_values = [row["value"] for row in rows if isinstance(row.get("value"), (int, float))]
    if not numeric_values:
        return None
    return float(max(numeric_values))


def _annotated_panel_nodes() -> list[Node]:
    nodes: list[Node] = []
    for node in build_full_panel().nodes:
        score = _representative_activity_score(node.id)
        if score is None:
            nodes.append(node)
            continue
        nodes.append(node.model_copy(update={"activity_score": score}))
    return nodes


def _canonical_curated_node(
    *,
    node_id: str,
    label: str,
    node_type: NodeType,
    detail: str,
    provenance: Iterable[ProvenanceRecord],
    group: str | None = None,
    reactivity: str | None = None,
    **extra_fields: Any,
) -> Node:
    return Node(
        id=node_id,
        label=label,
        type=node_type,
        detail=detail,
        group=group,
        reactivity=reactivity,
        match_status=MatchStatus.CANONICAL,
        canonical_id=node_id,
        canonical_label=label,
        canonical_namespace=_SHOWCASE_NAMESPACE,
        provenance=list(provenance),
        **extra_fields,
    )


def _carcinogen_node(spec: Spec) -> Node:
    citation = f"Curated carcinogen entry for {spec['label']}"
    base = Node(
        id=spec["id"],
        label=spec["label"],
        type=NodeType.CARCINOGEN,
        group=spec["group"],
        detail=spec["detail"],
        provenance=[_curated_ref(spec["id"], citation, spec["detail"])],
    )
    grounded = ground_node(base, classifier=_IARC)
    if grounded.match_status == MatchStatus.UNMATCHED:
        return grounded.model_copy(
            update={
                "match_status": MatchStatus.CANONICAL,
                "canonical_id": grounded.id,
                "canonical_label": grounded.label,
                "canonical_namespace": _SHOWCASE_NAMESPACE,
            },
        )
    return grounded


def _metabolite_node(spec: Spec) -> Node:
    citation = f"Curated metabolite entry for {spec['label']}"
    return _canonical_curated_node(
        node_id=spec["id"],
        label=spec["label"],
        node_type=NodeType.METABOLITE,
        detail=spec["detail"],
        reactivity=spec["reactivity"],
        provenance=[_curated_ref(spec["id"], citation, spec["detail"])],
    )


def _adduct_node(spec: Spec) -> Node:
    citation = f"Curated DNA adduct entry for {spec['label']}"
    return _canonical_curated_node(
        node_id=spec["id"],
        label=spec["label"],
        node_type=NodeType.DNA_ADDUCT,
        detail=spec["detail"],
        provenance=[_curated_ref(spec["id"], citation, spec["detail"])],
    )


def _pathway_nodes() -> list[Node]:
    return [
        _canonical_curated_node(
            node_id=pathway["pathway_id"],
            label=pathway["label"],
            node_type=NodeType.PATHWAY,
            detail=pathway["focus"],
            provenance=[_kegg_ref(pathway["pathway_id"], pathway["label"])],
        )
        for pathway in REFERENCE_KEGG_PATHWAYS
    ]


def _gene_node(spec: Spec) -> Node:
    citation = f"Curated gene entry for {spec['label']}"
    source_db = spec.get("source_db", _CURATED_SOURCE)
    return _canonical_curated_node(
        node_id=spec["id"],
        label=spec["label"],
        node_type=NodeType.GENE,
        detail=spec["detail"],
        provenance=[
            _curated_ref(spec["id"], citation, spec["detail"], source_db=source_db),
        ],
        source_db=source_db,
        tissue=spec.get("tissue"),
        variant=spec.get("variant"),
        phenotype=spec.get("phenotype"),
    )


def _tissue_node(spec: Spec) -> Node:
    citation = f"Curated tissue-context entry for {spec['label']}"
    return _canonical_curated_node(
        node_id=spec["id"],
        label=spec["label"],
        node_type=NodeType.TISSUE,
        detail=spec["detail"],
        provenance=[_curated_ref(spec["id"], citation, spec["detail"])],
        tissue=spec["label"],
    )


def _module_enzyme_node(spec: Spec) -> Node:
    citation = f"Curated enzyme entry for {spec['label']}"
    return _canonical_curated_node(
        node_id=spec["id"],
        label=spec["label"],
        node_type=NodeType.ENZYME,
        detail=spec["detail"],
        provenance=[_curated_ref(spec["id"], citation, spec["detail"])],
        phase="I",
        role="Mixed",
        tissue="prostate, liver, intestine",
    )


def _module_pathway_node(spec: Spec) -> Node:
    pathway_id = spec["pathway_id"]
    if pathway_id.startswith("hsa"):
        provenance = [_kegg_ref(pathway_id, spec["label"])]
        source_db = "KEGG"
    else:
        provenance = [
            _curated_ref(pathway_id, f"Curated pathway entry for {spec['label']}", spec["focus"])
        ]
        source_db = _CURATED_SOURCE
    return _canonical_curated_node(
        node_id=pathway_id,
        label=spec["label"],
        node_type=NodeType.PATHWAY,
        detail=spec["focus"],
        provenance=provenance,
        source_db=source_db,
    )


def _relation(spec: Spec) -> Edge:
    record_id = f"{spec['source']}->{spec['target']}"
    edge_type = EdgeType(spec["type"])
    if edge_type == EdgeType.PATHWAY and str(spec["target"]).startswith("hsa"):
        pathway_label = next(
            (
                entry["label"]
                for entry in REFERENCE_KEGG_PATHWAYS
                if entry["pathway_id"] == spec["target"]
            ),
            spec["target"],
        )
        provenance = [_kegg_ref(spec["target"], pathway_label)]
    else:
        provenance = [_curated_ref(record_id, spec["label"], spec["label"])]

    edge_kwargs: dict[str, Any] = {}
    if edge_type == EdgeType.CUSTOM:
        edge_kwargs["custom_predicate"] = spec["custom_predicate"]
        edge_kwargs["match_status"] = MatchStatus.CUSTOM

    return Edge(
        source=spec["source"],
        target=spec["target"],
        type=edge_type,
        label=spec["label"],
        carcinogen=spec.get("carcinogen"),
        provenance=provenance,
        **edge_kwargs,
    )


def _sorted_labels(nodes: Iterable[Node], node_type: NodeType) -> tuple[str, ...]:
    return tuple(sorted(node.label for node in nodes if node.type == node_type))


def _inventory_group(name: str, members: Iterable[str]) -> ArchitectureInventoryGroup:
    ordered_members = tuple(sorted(members))
    return ArchitectureInventoryGroup(
        name=name,
        members=ordered_members,
        count=len(ordered_members),
    )


def _carcinogen_class_groups(nodes: Iterable[Node]) -> tuple[ArchitectureInventoryGroup, ...]:
    labels_by_group: dict[str, list[str]] = {}
    for node in nodes:
        if node.type != NodeType.CARCINOGEN:
            continue
        group = str(node.group or "").strip()
        labels_by_group.setdefault(group, []).append(node.label)

    return tuple(
        _inventory_group(title, labels_by_group.get(group, ()))
        for group, title in _CARCINOGEN_CLASS_TITLES
    )


def _enzyme_category_groups(panel_nodes: Iterable[Node]) -> tuple[ArchitectureInventoryGroup, ...]:
    panel = list(panel_nodes)
    return (
        _inventory_group("Phase I", (node.label for node in panel if node.phase == "I")),
        _inventory_group("Phase II", (node.label for node in panel if node.phase == "II")),
        _inventory_group("Phase III", (node.label for node in panel if node.phase == "III")),
        _inventory_group(
            "DNA Repair",
            (
                node.label
                for node in panel
                if str(node.group or "").startswith("DNA Repair")
            ),
        ),
    )


def _build_legacy_full_legends_graph() -> KnowledgeGraph:
    panel_nodes = _annotated_panel_nodes()
    carcinogen_nodes = [_carcinogen_node(spec) for spec in _CARCINOGEN_SPECS]
    metabolite_nodes = [_metabolite_node(spec) for spec in _METABOLITE_SPECS]
    adduct_nodes = [_adduct_node(spec) for spec in _DNA_ADDUCT_SPECS]
    pathway_nodes = _pathway_nodes()

    pathway_label_by_id = {
        entry["pathway_id"]: entry["label"] for entry in REFERENCE_KEGG_PATHWAYS
    }
    pathway_edges = [
        Edge(
            source=node_id,
            target=pathway_id,
            type=EdgeType.PATHWAY,
            label=f"Pathway membership: {pathway_label_by_id[pathway_id]}",
            source_db="KEGG",
            provenance=[_kegg_ref(pathway_id, pathway_label_by_id[pathway_id])],
        )
        for pathway_id, node_ids in _PATHWAY_MEMBERSHIP.items()
        for node_id in node_ids
    ]

    return KnowledgeGraph(
        nodes=[
            *panel_nodes,
            *carcinogen_nodes,
            *metabolite_nodes,
            *adduct_nodes,
            *pathway_nodes,
        ],
        edges=[*[_relation(spec) for spec in _EDGE_SPECS], *pathway_edges],
    )


def _edge_identity(edge: Edge) -> tuple[str, str, str, str | None, str | None]:
    return (
        edge.source,
        edge.target,
        edge.type.value,
        edge.custom_predicate,
        edge.carcinogen,
    )


def _canonical_showcase_id(node_id: str) -> str:
    return _FULL_LEGENDS_ID_REMAP.get(node_id, node_id)


def _normalize_reference_graph_ids(graph: KnowledgeGraph) -> KnowledgeGraph:
    node_by_id: dict[str, Node] = {}
    for node in graph.nodes:
        node_payload = node.model_dump()
        node_payload["id"] = _canonical_showcase_id(node.id)
        if node_payload.get("canonical_id"):
            node_payload["canonical_id"] = _canonical_showcase_id(
                str(node_payload["canonical_id"])
            )
        normalized_node = Node(**node_payload)
        node_by_id.setdefault(normalized_node.id, normalized_node)

    edge_by_identity: dict[tuple[str, str, str, str | None, str | None], Edge] = {}
    for edge in graph.edges:
        edge_payload = edge.model_dump()
        edge_payload["source"] = _canonical_showcase_id(edge.source)
        edge_payload["target"] = _canonical_showcase_id(edge.target)
        if edge.carcinogen:
            edge_payload["carcinogen"] = _canonical_showcase_id(edge.carcinogen)
        normalized_edge = Edge(**edge_payload)
        edge_by_identity.setdefault(_edge_identity(normalized_edge), normalized_edge)

    return KnowledgeGraph(
        nodes=list(node_by_id.values()),
        edges=list(edge_by_identity.values()),
    )


def _overlay_legacy_showcase_metadata(
    graph: KnowledgeGraph,
    legacy_graph: KnowledgeGraph,
) -> KnowledgeGraph:
    legacy_node_by_id = {node.id: node for node in legacy_graph.nodes}
    legacy_edge_by_identity = {
        _edge_identity(edge): edge for edge in legacy_graph.edges
    }

    return KnowledgeGraph(
        nodes=[legacy_node_by_id.get(node.id, node) for node in graph.nodes],
        edges=[legacy_edge_by_identity.get(_edge_identity(edge), edge) for edge in graph.edges],
    )


def _build_base_full_legends_graph() -> KnowledgeGraph:
    legacy_graph = _build_legacy_full_legends_graph()
    if not _FULL_LEGENDS_REFERENCE_GRAPH.exists():
        return legacy_graph

    reference_graph = parse_graph_data_js(_FULL_LEGENDS_REFERENCE_GRAPH)
    normalized_graph = _normalize_reference_graph_ids(reference_graph)
    return _overlay_legacy_showcase_metadata(normalized_graph, legacy_graph)


def _merge_graphs(*graphs: KnowledgeGraph) -> KnowledgeGraph:
    node_by_id: dict[str, Node] = {}
    for graph in graphs:
        for node in graph.nodes:
            node_by_id.setdefault(node.id, node)

    edge_by_identity: dict[tuple[str, str, str, str | None, str | None], Edge] = {}
    for graph in graphs:
        for edge in graph.edges:
            edge_by_identity.setdefault(_edge_identity(edge), edge)

    return KnowledgeGraph(
        nodes=list(node_by_id.values()),
        edges=list(edge_by_identity.values()),
    )


def _merge_optional_heavy_metal_reference(
    graph: KnowledgeGraph,
    *,
    include_heavy_metals: bool,
) -> KnowledgeGraph:
    if not include_heavy_metals:
        return graph

    from .expanded_metals import load_heavy_metal_reference_graph

    return _merge_graphs(graph, load_heavy_metal_reference_graph())


def build_androgen_module_graph() -> KnowledgeGraph:
    """Return an optional androgen-metabolism module with receptor and variant context."""
    panel_lookup = {node.id: node for node in _annotated_panel_nodes()}
    shared_panel_nodes = [panel_lookup[node_id] for node_id in _ANDROGEN_PANEL_ENZYME_IDS]
    carcinogen_nodes = [
        _carcinogen_node(_BASE_CARCINOGEN_SPECS[node_id])
        for node_id in _ANDROGEN_CARCINOGEN_IDS
    ]
    metabolite_nodes = [_metabolite_node(spec) for spec in _ANDROGEN_METABOLITE_SPECS]
    adduct_nodes = [_adduct_node(spec) for spec in _ANDROGEN_DNA_ADDUCT_SPECS]
    pathway_nodes = [_module_pathway_node(spec) for spec in _ANDROGEN_PATHWAY_SPECS]
    gene_nodes = [_gene_node(spec) for spec in _ANDROGEN_GENE_SPECS]
    tissue_nodes = [_tissue_node(spec) for spec in _ANDROGEN_TISSUE_SPECS]
    extra_enzyme_nodes = [_module_enzyme_node(spec) for spec in _ANDROGEN_ENZYME_SPECS]

    return KnowledgeGraph(
        nodes=[
            *shared_panel_nodes,
            *extra_enzyme_nodes,
            *carcinogen_nodes,
            *metabolite_nodes,
            *adduct_nodes,
            *pathway_nodes,
            *gene_nodes,
            *tissue_nodes,
        ],
        edges=[_relation(spec) for spec in _ANDROGEN_EDGE_SPECS],
    )


def build_androgen_module_engine(
    *,
    mode: GraphMode | str = GraphMode.EXPLORATORY,
) -> GraphEngine:
    """Load :func:`build_androgen_module_graph` into a :class:`GraphEngine`."""
    engine = GraphEngine()
    engine.load(build_androgen_module_graph(), mode=mode)
    return engine


def build_full_legends_graph(
    *,
    include_heavy_metals: bool = False,
    include_androgen_module: bool = False,
) -> KnowledgeGraph:
    """Return a curated graph.

    By default this returns 107-node / 124-edge.
    Set ``include_heavy_metals=True`` to upgrade the showcase with the bundled
    module-01 heavy-metal reference layer (133 nodes / 166 edges).
    """
    graph = _build_base_full_legends_graph()
    graph = _merge_optional_heavy_metal_reference(
        graph,
        include_heavy_metals=include_heavy_metals,
    )
    if include_androgen_module:
        graph = _merge_graphs(graph, build_androgen_module_graph())
    return graph


def build_reference_graph(
    *,
    include_androgen_module: bool = False,
) -> KnowledgeGraph:
    """Return the canonical bundled multi-class reference graph.

    This graph composes the legacy showcase, bundled heavy-metal reference layer,
    and the Wave 2 carcinogen-class expansion into the fullest packaged graph
    currently shipped with ExposoGraph.
    """
    from .wave2_classes import merge_wave2_classes_into

    graph = build_full_legends_graph(include_heavy_metals=True)
    graph = merge_wave2_classes_into(graph)
    graph = _apply_reference_graph_overlay(graph)
    if include_androgen_module:
        graph = _merge_graphs(graph, build_androgen_module_graph())
    return graph


def _apply_reference_graph_overlay(graph: KnowledgeGraph) -> KnowledgeGraph:
    """Restore manuscript-relevant pathway context lost in the bundled overlay.

    The packaged D3 reference graph keeps the curated node inventory but drops a
    substantial subset of the legacy pathway-context edges. For the canonical
    reference graph, retain those pathway edges when both endpoints already
    exist in the current graph and restore the explicit vinyl-chloride
    N2,3-ethenoguanine adduct node from the legacy graph.
    """
    legacy = _build_legacy_full_legends_graph()
    legacy_node_by_id = {node.id: node for node in legacy.nodes}
    current_node_ids = {node.id for node in graph.nodes}
    current_edge_ids = {_edge_identity(edge) for edge in graph.edges}

    overlay_nodes: list[Node] = []
    overlay_edges: list[Edge] = []

    for edge in legacy.edges:
        if edge.type != EdgeType.PATHWAY:
            continue
        if edge.source not in current_node_ids or edge.target not in current_node_ids:
            continue
        identity = _edge_identity(edge)
        if identity in current_edge_ids:
            continue
        overlay_edges.append(edge)
        current_edge_ids.add(identity)

    ethenog = legacy_node_by_id["EthenoG"]
    if ethenog.id not in current_node_ids:
        overlay_nodes.append(ethenog)
        current_node_ids.add(ethenog.id)

    vc_guanine_edge = next(
        edge
        for edge in legacy.edges
        if edge.source == "Chloroethylene_oxide"
        and edge.target == "EthenoG"
        and edge.type == EdgeType.FORMS_ADDUCT
    )
    identity = _edge_identity(vc_guanine_edge)
    if (
        vc_guanine_edge.source in current_node_ids
        and vc_guanine_edge.target in current_node_ids
        and identity not in current_edge_ids
    ):
        overlay_edges.append(vc_guanine_edge)

    if not overlay_nodes and not overlay_edges:
        return graph

    if overlay_edges:
        graph_node_by_id = {node.id: node for node in graph.nodes}
        referenced_ids = {
            node_id
            for edge in overlay_edges
            for node_id in (edge.source, edge.target, edge.carcinogen)
            if node_id is not None
        }
        for node_id in sorted(referenced_ids):
            node = graph_node_by_id.get(node_id)
            if node is not None:
                overlay_nodes.append(node)

    return _merge_graphs(graph, KnowledgeGraph(nodes=overlay_nodes, edges=overlay_edges))


def _build_architecture_summary(graph: KnowledgeGraph) -> ArchitectureSummary:
    """Summarize node/edge inventories for a seeded reference graph."""
    enzyme_nodes = [node for node in graph.nodes if node.type == NodeType.ENZYME]

    node_type_counter = Counter(node.type.value for node in graph.nodes)
    edge_type_counter = Counter(edge.type.value for edge in graph.edges)

    node_type_counts = {
        node_type.value: node_type_counter[node_type.value]
        for node_type in NodeType
        if node_type_counter.get(node_type.value)
    }
    edge_type_counts = {
        edge_type.value: edge_type_counter[edge_type.value]
        for edge_type in EdgeType
        if edge_type_counter.get(edge_type.value)
    }

    return ArchitectureSummary(
        node_count=len(graph.nodes),
        edge_count=len(graph.edges),
        node_type_count=len(node_type_counter),
        edge_type_count=len(edge_type_counter),
        node_type_counts=node_type_counts,
        edge_type_counts=edge_type_counts,
        carcinogens=_sorted_labels(graph.nodes, NodeType.CARCINOGEN),
        enzymes=_sorted_labels(graph.nodes, NodeType.ENZYME),
        metabolites=_sorted_labels(graph.nodes, NodeType.METABOLITE),
        dna_adducts=_sorted_labels(graph.nodes, NodeType.DNA_ADDUCT),
        pathway_labels=_sorted_labels(graph.nodes, NodeType.PATHWAY),
        carcinogen_classes=_carcinogen_class_groups(graph.nodes),
        enzyme_categories=_enzyme_category_groups(enzyme_nodes),
    )


def build_full_legends_architecture_summary(
    *,
    include_heavy_metals: bool = False,
    include_androgen_module: bool = False,
) -> ArchitectureSummary:
    """Return a typed architecture summary derived from the seeded showcase graph."""
    graph = build_full_legends_graph(
        include_heavy_metals=include_heavy_metals,
        include_androgen_module=include_androgen_module,
    )
    return _build_architecture_summary(graph)


def build_reference_architecture_summary(
    *,
    include_androgen_module: bool = False,
) -> ArchitectureSummary:
    """Return a typed architecture summary for the bundled reference graph."""
    graph = build_reference_graph(include_androgen_module=include_androgen_module)
    return _build_architecture_summary(graph)


def build_full_legends_engine(
    *,
    mode: GraphMode | str = GraphMode.EXPLORATORY,
    include_heavy_metals: bool = False,
    include_androgen_module: bool = False,
) -> GraphEngine:
    """Load :func:`build_full_legends_graph` into a :class:`GraphEngine`."""
    engine = GraphEngine()
    engine.load(
        build_full_legends_graph(
            include_heavy_metals=include_heavy_metals,
            include_androgen_module=include_androgen_module,
        ),
        mode=mode,
    )
    return engine


def build_reference_engine(
    *,
    mode: GraphMode | str = GraphMode.EXPLORATORY,
    include_androgen_module: bool = False,
) -> GraphEngine:
    """Load :func:`build_reference_graph` into a :class:`GraphEngine`."""
    engine = GraphEngine()
    engine.load(
        build_reference_graph(include_androgen_module=include_androgen_module),
        mode=mode,
    )
    return engine


def write_full_legends_exports(
    output_dir: str | Path = "exports",
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    include_heavy_metals: bool = False,
    include_androgen_module: bool = False,
    bundle_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Write HTML, JSON, Plotly, and bundled viewer artifacts."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    engine = build_full_legends_engine(
        include_heavy_metals=include_heavy_metals,
        include_androgen_module=include_androgen_module,
    )

    artifact_parts = ["run_kg_full_legends"]
    if include_heavy_metals:
        artifact_parts.append("heavy_metals")
    if include_androgen_module:
        artifact_parts.append("androgen")
    artifact_stem = "_".join(artifact_parts)

    html_path = output_dir / f"{artifact_stem}.html"
    json_path = output_dir / f"{artifact_stem}.json"
    plotly_path = output_dir / f"{artifact_stem}_plotly.html"

    to_interactive_html(engine, html_path, visibility=visibility)
    to_json(engine, json_path, visibility=visibility)
    to_plotly_html(engine, plotly_path, visibility=visibility)

    bundle_target = (
        Path(bundle_dir)
        if bundle_dir is not None
        else (
            output_dir / f"{artifact_stem}.js"
            if include_heavy_metals or include_androgen_module
            else _DEFAULT_BUNDLE_PATH
        )
    )
    if bundle_target.suffix == ".js":
        bundle_target.parent.mkdir(parents=True, exist_ok=True)
        to_graph_data_js(engine, bundle_target, visibility=visibility)

    return {
        "html": html_path,
        "json": json_path,
        "plotly_html": plotly_path,
        "graph_data_js": bundle_target,
    }


def write_reference_exports(
    output_dir: str | Path = "exports",
    *,
    visibility: GraphVisibility | str = GraphVisibility.ALL,
    include_androgen_module: bool = False,
    bundle_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Write HTML, JSON, Plotly, and bundled viewer artifacts for the reference graph."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    engine = build_reference_engine(include_androgen_module=include_androgen_module)

    artifact_parts = ["run_kg_reference"]
    if include_androgen_module:
        artifact_parts.append("androgen")
    artifact_stem = "_".join(artifact_parts)

    html_path = output_dir / f"{artifact_stem}.html"
    json_path = output_dir / f"{artifact_stem}.json"
    plotly_path = output_dir / f"{artifact_stem}_plotly.html"

    to_interactive_html(engine, html_path, visibility=visibility)
    to_json(engine, json_path, visibility=visibility)
    to_plotly_html(engine, plotly_path, visibility=visibility)

    bundle_target = (
        Path(bundle_dir)
        if bundle_dir is not None
        else output_dir / f"{artifact_stem}.js"
    )
    if bundle_target.suffix == ".js":
        bundle_target.parent.mkdir(parents=True, exist_ok=True)
        to_graph_data_js(engine, bundle_target, visibility=visibility)

    return {
        "html": html_path,
        "json": json_path,
        "plotly_html": plotly_path,
        "graph_data_js": bundle_target,
    }
