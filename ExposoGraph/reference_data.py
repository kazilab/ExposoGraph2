"""Curated reference gene panels and activity scores for ExposoGraph.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .models import KnowledgeGraph, Node, NodeType

NCBI_GENE_IDS: dict[str, str] = {
    "ABCB1": "5243",
    "ABCC2": "1244",
    "ABCG2": "9429",
    "AKR1C3": "8644",
    "COMT": "1312",
    "CYP1A1": "1543",
    "CYP1A2": "1544",
    "CYP1B1": "1545",
    "CYP2A6": "1548",
    "CYP2A13": "1553",
    "CYP17A1": "1586",
    "CYP2C19": "1557",
    "CYP2C9": "1559",
    "CYP2D6": "1565",
    "CYP2E1": "1571",
    "CYP19A1": "1588",
    "CYP3A4": "1576",
    "CYP3A5": "1577",
    "EPHX1": "2052",
    "ERCC2": "2068",
    "GSTM1": "2944",
    "GSTP1": "2950",
    "GSTT1": "2952",
    "MLH1": "4292",
    "MGMT": "4255",
    "MSH2": "4436",
    "NAT1": "6530",
    "NAT2": "10",
    "NQO1": "1728",
    "OGG1": "4968",
    "SRD5A1": "6715",
    "SRD5A2": "6716",
    "SULT1A1": "6817",
    "SULT1E1": "6783",
    "UGT1A1": "54658",
    "UGT2B17": "7367",
    "UGT2B7": "7364",
    "XPC": "7508",
    "XRCC1": "7515",
    # Wave 2 class-specific genes
    "ALDH2": "217",
    "ADH1B": "125",
    "ADH1C": "126",
    "ADH5": "128",
    "FANCD2": "2177",
    "AHR": "196",
    "ARNT": "405",
    "AHRR": "57491",
    "CCBL1": "883",
}

# ClinPGx (formerly PharmGKB) accession IDs — verified via PharmGKB API
CLINPGX_ACCESSIONS: dict[str, str] = {
    "ABCB1": "PA267",
    "ABCC2": "PA116",
    "ABCG2": "PA390",
    "AKR1C3": "AKR1C3",
    "COMT": "PA119",
    "CYP1A1": "PA27092",
    "CYP1A2": "PA27093",
    "CYP1B1": "PA27094",
    "CYP2A6": "PA121",
    "CYP2A13": "PA27101",
    "CYP17A1": "CYP17A1",
    "CYP2C9": "PA126",
    "CYP2C19": "PA124",
    "CYP2D6": "PA128",
    "CYP2E1": "PA129",
    "CYP19A1": "CYP19A1",
    "CYP3A4": "PA130",
    "CYP3A5": "CYP3A5",
    "EPHX1": "PA27829",
    "ERCC2": "PA27848",
    "GSTM1": "PA182",
    "GSTP1": "PA29028",
    "GSTT1": "PA183",
    "MLH1": "PA222",
    "MGMT": "PA239",
    "MSH2": "PA283",
    "NAT1": "PA17",
    "NAT2": "PA18",
    "NQO1": "PA31744",
    "OGG1": "PA31912",
    "SRD5A1": "SRD5A1",
    "SRD5A2": "SRD5A2",
    "SULT1A1": "PA343",
    "SULT1E1": "SULT1E1",
    "UGT1A1": "PA420",
    "UGT2B17": "UGT2B17",
    "UGT2B7": "PA361",
    "XPC": "PA37413",
    "XRCC1": "PA369",
    # Wave 2 class-specific genes
    "ALDH2": "PA443",
    "ADH1B": "PA24495",
    "ADH1C": "PA24496",
    "ADH5": "PA24498",
    "FANCD2": "PA28027",
    "AHR": "PA61",
    "ARNT": "PA24717",
    "AHRR": "PA134872392",
    "CCBL1": "PA26127",
}

DNA_REPAIR_CLASSES: dict[str, str] = {
    "XRCC1": "DNA Repair (BER)",
    "OGG1": "DNA Repair (BER)",
    "XPC": "DNA Repair (NER)",
    "ERCC2": "DNA Repair (NER)",
    "MGMT": "DNA Repair (Direct Reversal)",
    "MLH1": "DNA Repair (MMR)",
    "MSH2": "DNA Repair (MMR)",
    "FANCD2": "DNA Repair (FA)",
}


def _ncbi_gene_ref(symbol: str) -> dict[str, str]:
    gene_id = NCBI_GENE_IDS[symbol]
    return {
        "source_db": "NCBI Gene",
        "record_id": gene_id,
        "citation": f"NCBI Gene record for {symbol} (Homo sapiens)",
        "url": f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
        "evidence": "Canonical human gene identifier and nomenclature.",
    }


def _gtex_ref(symbol: str, tissue: str) -> dict[str, str]:
    return {
        "source_db": "GTEx via HPA",
        "record_id": symbol,
        "citation": f"GTEx v8 nTPM expression via Human Protein Atlas v25 for {symbol}",
        "url": "https://www.proteinatlas.org/download/tsv/rna_tissue_detail_gtex.tsv.zip",
        "tissue": tissue,
        "evidence": "Human tissue-expression context from GTEx v8 nTPM values distributed by Human Protein Atlas v25.",
    }


def _clinpgx_ref(symbol: str) -> dict[str, str]:
    accession = CLINPGX_ACCESSIONS.get(symbol, symbol)
    return {
        "source_db": "ClinPGx",
        "record_id": accession,
        "citation": f"ClinPGx gene resource for {symbol}",
        "url": f"https://www.clinpgx.org/gene/{accession}",
        "evidence": "Pharmacogenomics and xenobiotic-response gene resource.",
    }


def _pharmvar_ref(symbol: str) -> dict[str, str]:
    return {
        "source_db": "PharmVar",
        "record_id": symbol,
        "citation": f"PharmVar allele definition resource for {symbol}",
        "url": f"https://www.pharmvar.org/gene/{symbol}",
        "evidence": "Star-allele nomenclature and haplotype definition resource.",
    }


def _cpic_ref(record_id: str, citation: str) -> dict[str, str]:
    return {
        "source_db": "CPIC",
        "record_id": record_id,
        "citation": citation,
        "url": "https://cpicpgx.org/guidelines/",
        "evidence": "Guideline-backed genotype-to-phenotype translation and activity score conventions.",
    }


def _ctd_ref(symbol: str, repair_class: str) -> dict[str, str]:
    return {
        "source_db": "CTD",
        "record_id": symbol,
        "citation": f"Comparative Toxicogenomics Database record for {symbol}",
        "url": f"https://ctdbase.org/detail.go?type=gene&acc={symbol}",
        "evidence": f"Toxicogenomic context supporting the {repair_class} assignment.",
    }


def _pubmed_ref(pmid: str, citation: str, evidence: str) -> dict[str, str]:
    return {
        "source_db": "PubMed",
        "record_id": pmid,
        "pmid": pmid,
        "citation": citation,
        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        "evidence": evidence,
    }


def _gene_provenance(symbol: str, *, tissue: str, repair_class: str | None = None) -> list[dict[str, str]]:
    records = [_ncbi_gene_ref(symbol), _gtex_ref(symbol, tissue)]
    if repair_class:
        records.append(_ctd_ref(symbol, repair_class))
    else:
        records.append(_clinpgx_ref(symbol))
    return records


def _gene_entry(
    symbol: str,
    *,
    detail: str,
    tissue: str,
    role: str,
    phase: str | None = None,
    group: str | None = None,
    label: str | None = None,
) -> dict[str, Any]:
    return {
        "id": symbol,
        "label": label or symbol,
        "detail": detail,
        "tissue": tissue,
        "role": role,
        "phase": phase,
        "group": group,
        "provenance": _gene_provenance(symbol, tissue=tissue, repair_class=group),
    }


def _activity_score_meta(
    evidence_basis: str,
    note: str,
    *references: dict[str, str],
) -> dict[str, object]:
    return {
        "evidence_basis": evidence_basis,
        "note": note,
        "references": list(references),
    }


REFERENCE_KEGG_PATHWAYS: list[dict[str, str]] = [
    {
        "pathway_id": "hsa00980",
        "label": "Xenobiotic metabolism by cytochrome P450",
        "focus": "Core xenobiotic biotransformation pathway for phase I carcinogen metabolism.",
    },
    {
        "pathway_id": "hsa00140",
        "label": "Steroid hormone biosynthesis",
        "focus": "Hormone-metabolism context relevant to estrogenic carcinogen processing.",
    },
    {
        "pathway_id": "hsa05204",
        "label": "Chemical carcinogenesis - DNA adducts",
        "focus": "Reactive-metabolite to DNA-adduct formation and repair context.",
    },
    {
        "pathway_id": "hsa05208",
        "label": "Chemical carcinogenesis - reactive oxygen species",
        "focus": "ROS-linked carcinogenesis and oxidative damage context.",
    },
]


CURATION_SOURCE_MANIFEST: dict[str, dict[str, dict[str, object]]] = {
    "primary_sources": {
        "IARC": {
            "role": "Carcinogen hazard classification",
            "implementation": "Bundled static lookup used to annotate carcinogen nodes and grounding indexes.",
            "notes": (
                "Current implementation stores a curated subset of common carcinogens rather than a "
                "full monograph ETL."
            ),
        },
        "KEGG": {
            "role": "Reference pathway and pathway-membership context",
            "implementation": "KEGG REST seeding plus a curated reference pathway catalog.",
            "reference_pathways": [pathway["pathway_id"] for pathway in REFERENCE_KEGG_PATHWAYS],
        },
        "PharmVar": {
            "role": "Star-allele nomenclature and haplotype definitions for supported pharmacogenes",
            "implementation": "Referenced in activity-score metadata where allele-definition resources are used.",
        },
        "CPIC": {
            "role": "Genotype-to-phenotype translation conventions for guideline-backed activity tables",
            "implementation": (
                "Explicitly referenced for guideline-backed activity metadata such as CYP2D6, CYP2C9, "
                "CYP2C19, and UGT1A1."
            ),
        },
        "CTD": {
            "role": "Chemical-gene interactions and DNA-repair toxicogenomic context",
            "implementation": "Seeder client plus repair-gene provenance records.",
        },
        "GTEx": {
            "role": "Human tissue-expression context for seeded genes",
            "implementation": "Structured gene-level provenance for all seeded panel genes.",
        },
        "PubMed": {
            "role": "Primary literature support for functional activity annotations",
            "implementation": "Per-gene activity metadata references with PMID-backed citations.",
        },
    },
    "supporting_sources": {
        "NCBI Gene": {
            "role": "Stable human gene identifiers and canonical nomenclature",
            "implementation": "Base provenance anchor for all seeded genes.",
        },
        "ClinPGx": {
            "role": "Pharmacogene coverage and implementation-facing gene resources",
            "implementation": (
                "Structured provenance for non-repair seed genes and activity metadata for "
                "supported pharmacogenes."
            ),
        },
    },
}

# ── Tier 1: Core Carcinogen-Metabolizing Enzyme Panel (13 genes) ─────────

TIER1_GENES: list[dict[str, Any]] = [
    _gene_entry(
        "CYP1A1",
        phase="I",
        role="Activation",
        detail="PAH diol-epoxide formation; AhR-inducible; extrahepatic expression",
        tissue="lung, placenta, lymphocytes",
    ),
    _gene_entry(
        "CYP1A2",
        phase="I",
        role="Activation",
        detail="HCA and aromatic amine activation; AhR-inducible; hepatic",
        tissue="liver",
    ),
    _gene_entry(
        "CYP1B1",
        phase="I",
        role="Activation",
        detail="PAH and estrogen activation; 4-OH-estradiol formation",
        tissue="breast, prostate, uterus, lung",
    ),
    _gene_entry(
        "CYP2A6",
        phase="I",
        role="Activation",
        detail="NNK and nitrosamine activation; nicotine metabolism",
        tissue="liver, nasal mucosa",
    ),
    _gene_entry(
        "CYP2E1",
        phase="I",
        role="Activation",
        detail="Small-molecule carcinogen activation; ethanol, benzene, NDMA",
        tissue="liver, lung",
    ),
    _gene_entry(
        "CYP3A4",
        phase="I",
        role="Activation",
        detail="Aflatoxin B1 8,9-epoxidation; broad substrate range",
        tissue="liver, intestine",
    ),
    _gene_entry(
        "GSTM1",
        phase="II",
        role="Detoxification",
        detail="GSH conjugation of PAH diol-epoxides and AFB1-epoxide; null polymorphism common",
        tissue="liver, lung",
    ),
    _gene_entry(
        "GSTT1",
        phase="II",
        role="Detoxification",
        detail="GSH conjugation of small electrophiles; null polymorphism common",
        tissue="liver, kidney, blood cells",
    ),
    _gene_entry(
        "GSTP1",
        phase="II",
        role="Detoxification",
        detail="GSH conjugation of BPDE and other PAH diol-epoxides",
        tissue="lung, brain, placenta",
    ),
    _gene_entry(
        "NAT2",
        phase="II",
        role="Mixed",
        detail="N-acetylation of aromatic amines; rapid vs slow acetylator phenotypes",
        tissue="liver, intestine",
    ),
    _gene_entry(
        "EPHX1",
        phase="II",
        role="Mixed",
        detail="Microsomal epoxide hydrolysis; both activation and detoxification roles",
        tissue="liver, lung",
    ),
    _gene_entry(
        "UGT1A1",
        phase="II",
        role="Detoxification",
        detail="Glucuronidation of PAH metabolites and bilirubin",
        tissue="liver, intestine",
    ),
    _gene_entry(
        "NQO1",
        phase="II",
        role="Detoxification",
        detail="Two-electron quinone reduction; prevents ROS from redox cycling",
        tissue="liver, lung, colon",
    ),
]

# ── Tier 2: Extended Gene Panel (25 genes) ───────────────────────────────

TIER2_GENES: list[dict[str, Any]] = [
    # Phase II
    _gene_entry(
        "SULT1A1",
        phase="II",
        role="Detoxification",
        detail="Sulfation of phenolic compounds, HCA intermediates, PAH metabolites",
        tissue="liver, GI tract, platelets, brain",
    ),
    _gene_entry(
        "NAT1",
        phase="II",
        role="Mixed",
        detail="O-acetylation and N-acetylation of aromatic/heterocyclic amines in peripheral tissues",
        tissue="bladder, colon, breast, lung",
    ),
    _gene_entry(
        "UGT2B7",
        phase="II",
        role="Detoxification",
        detail="Glucuronidation of steroid hormones, carcinogen metabolites",
        tissue="liver, kidney, GI tract, mammary",
    ),
    _gene_entry(
        "UGT2B17",
        phase="II",
        role="Detoxification",
        detail="Glucuronidation of testosterone, DHT, and related androgen metabolites",
        tissue="liver, prostate, intestine",
    ),
    _gene_entry(
        "SULT1E1",
        phase="II",
        role="Detoxification",
        detail="High-affinity estrogen sulfotransferase; sulfation and inactivation of estradiol and catechol estrogens",
        tissue="liver, endometrium, breast, placenta",
    ),
    _gene_entry(
        "COMT",
        phase="II",
        role="Detoxification",
        detail="O-methylation of catechol estrogens and other catechols; limits redox-cycling estrogen metabolites",
        tissue="liver, breast, brain, endometrium",
    ),
    # Phase III transporters
    _gene_entry(
        "ABCB1",
        phase="III",
        role="Transport",
        detail="P-gp; ATP-driven efflux of hydrophobic xenobiotics and carcinogens",
        tissue="liver, intestine, kidney, BBB, placenta",
    ),
    _gene_entry(
        "ABCC2",
        phase="III",
        role="Transport",
        detail="MRP2; export of GSH/glucuronide/sulfate conjugates of carcinogen metabolites",
        tissue="liver, kidney, intestine",
    ),
    _gene_entry(
        "ABCG2",
        phase="III",
        role="Transport",
        detail="BCRP; efflux of PAHs, PhIP, porphyrins",
        tissue="intestine, liver, placenta, mammary, BBB",
    ),
    # DNA repair
    _gene_entry(
        "XRCC1",
        group=DNA_REPAIR_CLASSES["XRCC1"],
        role="Repair",
        detail="BER scaffold protein; coordinates repair of oxidative/alkylation DNA damage",
        tissue="ubiquitous",
    ),
    _gene_entry(
        "XPC",
        group=DNA_REPAIR_CLASSES["XPC"],
        role="Repair",
        detail="GG-NER damage sensor; recognizes bulky DNA adducts from PAHs",
        tissue="ubiquitous",
    ),
    _gene_entry(
        "ERCC2",
        label="ERCC2/XPD",
        group=DNA_REPAIR_CLASSES["ERCC2"],
        role="Repair",
        detail="NER helicase; unwinds DNA at damage sites for bulky adduct excision",
        tissue="ubiquitous",
    ),
    _gene_entry(
        "OGG1",
        group=DNA_REPAIR_CLASSES["OGG1"],
        role="Repair",
        detail="8-oxoguanine DNA glycosylase; excises oxidative DNA damage (8-oxodG)",
        tissue="ubiquitous (nuclear + mitochondrial)",
    ),
    _gene_entry(
        "MGMT",
        group=DNA_REPAIR_CLASSES["MGMT"],
        role="Repair",
        detail="Direct reversal of O6-alkylguanine; single-use suicidal repair enzyme",
        tissue="liver, colon, lung, brain",
    ),
    _gene_entry(
        "MLH1",
        group=DNA_REPAIR_CLASSES["MLH1"],
        role="Repair",
        detail="MutL homolog 1; mismatch repair of replication errors past DNA adducts",
        tissue="ubiquitous",
    ),
    _gene_entry(
        "MSH2",
        group=DNA_REPAIR_CLASSES["MSH2"],
        role="Repair",
        detail="MutS homolog 2; mismatch recognition during post-replicative repair at adduct sites",
        tissue="ubiquitous",
    ),
    # Additional CYPs and hormone-metabolism enzymes
    _gene_entry(
        "CYP2C9",
        phase="I",
        role="Mixed",
        detail="Oxidation of some PAH metabolites; major drug-metabolizing CYP",
        tissue="liver, intestine",
    ),
    _gene_entry(
        "CYP2C19",
        phase="I",
        role="Mixed",
        detail="Minor procarcinogen activation; possible nitrosamine metabolism",
        tissue="liver, intestine",
    ),
    _gene_entry(
        "CYP2D6",
        phase="I",
        role="Mixed",
        detail="Minor NNK activation; important for dual PGx/carcinogen-risk reporting",
        tissue="liver, brain, lung, GI tract",
    ),
    _gene_entry(
        "CYP2A13",
        phase="I",
        role="Activation",
        detail="Primary lung NNK-metabolizing CYP; tobacco-smoke carcinogen activation",
        tissue="lung, nasal mucosa",
    ),
    _gene_entry(
        "CYP17A1",
        phase="I",
        role="Activation",
        detail="Steroid 17alpha-hydroxylase/17,20-lyase; androgen precursor synthesis in hormone-linked carcinogen context",
        tissue="adrenal, gonad, prostate",
    ),
    _gene_entry(
        "SRD5A1",
        phase="I",
        role="Activation",
        detail="5alpha-reductase type 1; converts testosterone to DHT in peripheral hormone-metabolism tissues",
        tissue="skin, liver, prostate, breast",
    ),
    _gene_entry(
        "SRD5A2",
        phase="I",
        role="Activation",
        detail="5alpha-reductase type 2; major DHT-forming enzyme in prostate and urogenital tissues",
        tissue="prostate, seminal vesicle, genital skin, liver",
    ),
    _gene_entry(
        "CYP19A1",
        phase="I",
        role="Activation",
        detail="Aromatase; converts androgen precursors to estrogens in breast and other hormone-responsive tissues",
        tissue="adipose, ovary, breast, placenta, brain",
    ),
    _gene_entry(
        "AKR1C3",
        phase="I",
        role="Mixed",
        detail="17-ketosteroid reductase and quinone reductase; local androgen/estrogen activation with redox-metabolism overlap",
        tissue="prostate, breast, liver, endometrium",
    ),
]


# ── Wave 2: Carcinogen-Class-Specific Enzyme Panel ────────────────────────
# Nine genes required by Wave 2 carcinogen classes (Aldehydes, Dioxins/AhR,
# Dietary N-Nitroso, Chlorinated Solvents).

WAVE2_GENES: list[dict[str, Any]] = [
    _gene_entry(
        "ALDH2",
        phase="II",
        role="Detoxification",
        detail="Mitochondrial aldehyde dehydrogenase 2; primary acetaldehyde and formaldehyde clearance enzyme. rs671 (Glu504Lys, ALDH2*2) causes 6-19x acetaldehyde accumulation; affects ~540 million people globally (East Asian 25-35% allele frequency)",
        tissue="liver, esophagus, stomach",
    ),
    _gene_entry(
        "ADH1B",
        phase="I",
        role="Activation",
        detail="Alcohol dehydrogenase 1B; oxidizes ethanol to acetaldehyde. rs1229984 (His47Arg, ADH1B*2) creates a fast-metabolizer allele (~40x higher Vmax); 70-80% allele frequency in East Asian populations",
        tissue="liver, stomach",
    ),
    _gene_entry(
        "ADH1C",
        phase="I",
        role="Activation",
        detail="Alcohol dehydrogenase 1C (gamma subunit); contributes to ethanol oxidation. rs698 (Ile349Val) modulates catalytic rate",
        tissue="liver, stomach",
    ),
    _gene_entry(
        "ADH5",
        phase="I",
        role="Detoxification",
        detail="Formaldehyde dehydrogenase (class III ADH); cytosolic formaldehyde clearance via S-hydroxymethylglutathione pathway. Dual ADH5+ALDH2 deficiency causes lethal formaldehyde accumulation and bone marrow failure",
        tissue="liver, ubiquitous",
    ),
    _gene_entry(
        "FANCD2",
        phase=None,
        role="Repair",
        group="DNA Repair (FA)",
        detail="Fanconi anemia complementation group D2; repairs acetaldehyde-induced DNA interstrand crosslinks (ICLs). Essential for hematopoietic stem cell protection from endogenous aldehydes",
        tissue="bone marrow, ubiquitous",
    ),
    _gene_entry(
        "AHR",
        phase=None,
        role="Transcription Factor",
        detail="Aryl hydrocarbon receptor; ligand-activated transcription factor that induces CYP1A1/CYP1B1 upon binding TCDD or dioxin-like compounds. rs2066853 (Arg554Lys) modulates transcriptional activity",
        tissue="liver, lung, skin, immune cells",
    ),
    _gene_entry(
        "ARNT",
        phase=None,
        role="Transcription Factor",
        detail="AhR nuclear translocator (HIF-1beta); obligate co-factor for AhR-mediated transcription. Forms AhR-ARNT heterodimer that binds XRE/DRE elements. rs2228099 (Ile471Val)",
        tissue="ubiquitous",
    ),
    _gene_entry(
        "AHRR",
        phase=None,
        role="Transcription Factor",
        detail="AhR repressor; negative feedback regulator that competes with AhR for ARNT binding. cg05575921 CpG hypomethylation is a validated blood-based biomarker for tobacco/TCDD exposure",
        tissue="ubiquitous",
    ),
    _gene_entry(
        "CCBL1",
        phase="II",
        role="Bioactivation",
        detail="Cysteine conjugate beta-lyase 1; renal proximal tubule enzyme that bioactivates DCVC to reactive DCVT thiol, causing DNA adducts and VHL tumor suppressor mutations. rs2293968, rs2280841, rs2259043, rs941960 tagging SNPs associated with elevated RCC risk",
        tissue="kidney (proximal tubule)",
    ),
]


# ── Activity Score Reference (per-allele values, CPIC-format) ────────────
# Key = gene, value = list of {allele, value, phenotype, confidence}

ACTIVITY_SCORES: dict[str, list[dict[str, Any]]] = {
    "CYP2D6": [
        {"allele": "*1, *2, *35", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*9, *17, *29, *41", "value": 0.5, "phenotype": "NM (lower range)", "confidence": "High"},
        {"allele": "*10", "value": 0.25, "phenotype": "NM (lowest); IM", "confidence": "High"},
        {"allele": "*3, *4, *5, *6, *40", "value": 0.0, "phenotype": "Poor Metabolizer", "confidence": "High"},
        {"allele": "*1x2, *2x2 (dupl.)", "value": 2.0, "phenotype": "Ultrarapid Metabolizer", "confidence": "High"},
    ],
    "CYP2C9": [
        {"allele": "*1", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*2", "value": 0.5, "phenotype": "Intermediate Metabolizer", "confidence": "High"},
        {"allele": "*3, *5, *6, *8, *11", "value": 0.0, "phenotype": "IM; PM (AS 0)", "confidence": "High"},
    ],
    "CYP2C19": [
        {"allele": "*1", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*17", "value": 1.5, "phenotype": "Rapid Metabolizer", "confidence": "High"},
        {"allele": "*2, *3, *4, *5", "value": 0.0, "phenotype": "Poor Metabolizer", "confidence": "High"},
        {"allele": "*9", "value": 0.5, "phenotype": "Intermediate Metabolizer", "confidence": "High"},
    ],
    "NAT2": [
        {"allele": "*4, *1 (rapid)", "value": 1.0, "phenotype": "Rapid Acetylator", "confidence": "High"},
        {"allele": "*5, *6, *7, *14 (slow)", "value": 0.5, "phenotype": "Intermediate Acetylator", "confidence": "High"},
        {"allele": "*5/*6, *6/*7, etc.", "value": 0.5, "phenotype": "Slow (Poor) Acetylator", "confidence": "High"},
    ],
    "UGT1A1": [
        {"allele": "*1, *36", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*28, *6, *37", "value": 0.5, "phenotype": "Intermediate Metabolizer", "confidence": "High"},
        {"allele": "*28/*28, *6/*6", "value": 0.5, "phenotype": "Poor Metabolizer", "confidence": "High"},
    ],
    "CYP1A1": [
        {"allele": "*1 (reference)", "value": 1.0, "phenotype": "Normal activity", "confidence": "Moderate"},
        {"allele": "*2A (MspI, 3'-UTR)", "value": 1.2, "phenotype": "Increased inducibility", "confidence": "Moderate"},
        {"allele": "*2B (MspI + Ile462Val)", "value": 1.3, "phenotype": "High activity (UM-like)", "confidence": "Moderate"},
        {"allele": "*4 (Thr461Val)", "value": 1.2, "phenotype": "Increased activity", "confidence": "Moderate"},
    ],
    "CYP1A2": [
        {"allele": "*1A (reference)", "value": 1.0, "phenotype": "Normal (constitutive)", "confidence": "Moderate"},
        {"allele": "*1F (-163C>A) [induced]", "value": 1.5, "phenotype": "Ultrarapid (induced)", "confidence": "Moderate"},
        {"allele": "*1K (-163A + -729T)", "value": 0.5, "phenotype": "Decreased activity", "confidence": "Moderate"},
        {"allele": "*24 (W84X stop)", "value": 0.0, "phenotype": "IM (one null allele)", "confidence": "Moderate"},
    ],
    "CYP1B1": [
        {"allele": "*1 (reference)", "value": 1.0, "phenotype": "Normal activity", "confidence": "Moderate"},
        {"allele": "Val432 (C allele)", "value": 1.3, "phenotype": "Increased 4-OH-E2 formation", "confidence": "Moderate"},
    ],
    "CYP2A6": [
        {"allele": "*1A, *1B (reference)", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*1x2 (gene duplication)", "value": 2.0, "phenotype": "Ultrarapid Metabolizer", "confidence": "High"},
        {"allele": "*9 (-48T>G, TATA box)", "value": 0.5, "phenotype": "Intermediate Metabolizer", "confidence": "High"},
        {"allele": "*2, *4 (null alleles)", "value": 0.0, "phenotype": "Poor Metabolizer", "confidence": "High"},
    ],
    "CYP2E1": [
        {"allele": "c1/c1 (reference; *1A)", "value": 1.0, "phenotype": "Normal activity", "confidence": "Moderate"},
        {"allele": "c2 allele (-1019C>T)", "value": 1.3, "phenotype": "Increased transcription", "confidence": "Moderate"},
        {"allele": "c2/c2 homozygous", "value": 1.3, "phenotype": "High activity", "confidence": "Moderate"},
    ],
    "CYP3A4": [
        {"allele": "*1 (reference)", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*22 (intron 6 C>T)", "value": 0.5, "phenotype": "Intermediate Metabolizer", "confidence": "High"},
        {"allele": "*17, *20 (null-like)", "value": 0.0, "phenotype": "IM (one null allele)", "confidence": "Moderate"},
    ],
    "EPHX1": [
        {"allele": "Tyr113 + His139 (ref)", "value": 1.0, "phenotype": "Normal activity", "confidence": "High"},
        {"allele": "113His (exon 3, slow)", "value": 0.5, "phenotype": "Slow activity (~50%)", "confidence": "High"},
        {"allele": "139Arg (exon 4, fast)", "value": 1.25, "phenotype": "Fast activity (~25% increase)", "confidence": "High"},
    ],
    "GSTM1": [
        {"allele": "Non-null (+/+ or +/null)", "value": 1.0, "phenotype": "Present (functional)", "confidence": "High"},
        {"allele": "Null (homozygous deletion)", "value": 0.0, "phenotype": "Absent (no function)", "confidence": "High"},
    ],
    "GSTT1": [
        {"allele": "Non-null (+/+ or +/null)", "value": 1.0, "phenotype": "Present (functional)", "confidence": "High"},
        {"allele": "Null (homozygous deletion)", "value": 0.0, "phenotype": "Absent (no function)", "confidence": "High"},
    ],
    "NQO1": [
        {"allele": "*1 (Pro187, reference)", "value": 1.0, "phenotype": "Normal activity", "confidence": "High"},
        {"allele": "*2 (Pro187Ser, C609T)", "value": 0.25, "phenotype": "Decreased (~25% hetero)", "confidence": "High"},
        {"allele": "*2/*2 (Ser/Ser)", "value": 0.0, "phenotype": "No function (2-4%)", "confidence": "High"},
    ],
    "XPC": [
        {"allele": "Lys939 (reference)", "value": 1.0, "phenotype": "Normal NER capacity", "confidence": "Moderate"},
        {"allele": "Gln939 (rs2228001)", "value": 0.7, "phenotype": "Reduced NER", "confidence": "Moderate"},
    ],
    "OGG1": [
        {"allele": "Ser326 (reference)", "value": 1.0, "phenotype": "Normal 8-oxoG repair", "confidence": "High"},
        {"allele": "Cys326 (rs1052133)", "value": 0.5, "phenotype": "Reduced repair", "confidence": "High"},
    ],
    "MGMT": [
        {"allele": "Unmethylated promoter", "value": 1.0, "phenotype": "Normal O6-MeG repair", "confidence": "High"},
        {"allele": "Partially methylated", "value": 0.5, "phenotype": "Reduced expression", "confidence": "Moderate"},
        {"allele": "Hypermethylated", "value": 0.0, "phenotype": "Silenced (no repair)", "confidence": "High"},
    ],
    "MLH1": [
        {"allele": "Wild-type", "value": 1.0, "phenotype": "Normal MMR capacity", "confidence": "Moderate"},
        {"allele": "Promoter hypermethylation", "value": 0.0, "phenotype": "Silenced (Lynch-like)", "confidence": "Moderate"},
        {"allele": "Pathogenic variant (heterozygous)", "value": 0.5, "phenotype": "Reduced MMR (Lynch syndrome carrier)", "confidence": "Moderate"},
    ],
    "MSH2": [
        {"allele": "Wild-type", "value": 1.0, "phenotype": "Normal MMR capacity", "confidence": "Moderate"},
        {"allele": "Pathogenic variant (heterozygous)", "value": 0.5, "phenotype": "Reduced MMR (Lynch syndrome carrier)", "confidence": "Moderate"},
        {"allele": "Biallelic loss", "value": 0.0, "phenotype": "Absent MMR (constitutional MMR deficiency)", "confidence": "Moderate"},
    ],
    # Wave 2 class-specific activity scores
    "ALDH2": [
        {"allele": "*1/*1 (Glu504, reference)", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*1/*2 (Glu504Lys, rs671)", "value": 0.1, "phenotype": "Poor Metabolizer (6-19x AcH accumulation)", "confidence": "High"},
        {"allele": "*2/*2 (Lys504Lys)", "value": 0.0, "phenotype": "Null (no measurable activity)", "confidence": "High"},
    ],
    "ADH1B": [
        {"allele": "*1/*1 (Arg48, reference)", "value": 1.0, "phenotype": "Normal Metabolizer", "confidence": "High"},
        {"allele": "*1/*2 (His48, rs1229984)", "value": 20.0, "phenotype": "Rapid Metabolizer (~40x Vmax heterozygote)", "confidence": "High"},
        {"allele": "*2/*2 (His48His)", "value": 40.0, "phenotype": "Ultrarapid Metabolizer (~40x Vmax)", "confidence": "High"},
    ],
}


ACTIVITY_SCORE_METADATA: dict[str, dict[str, object]] = {
    "CYP2D6": _activity_score_meta(
        "Guideline-backed pharmacogene",
        "PharmVar-backed allele definitions are paired with CPIC phenotype-translation and activity-score conventions.",
        _clinpgx_ref("CYP2D6"),
        _pharmvar_ref("CYP2D6"),
        _cpic_ref(
            "CYP2D6-CYP2C19-TCA-2016",
            "Clinical Pharmacogenetics Implementation Consortium Guideline (CPIC) for CYP2D6 and CYP2C19 Genotypes and Dosing of Tricyclic Antidepressants: 2016 Update.",
        ),
        _pubmed_ref(
            "31562822",
            "Clinical Pharmacogenetics Implementation Consortium Guideline (CPIC) for CYP2D6 and CYP2C19 Genotypes and Dosing of Tricyclic Antidepressants: 2016 Update.",
            "CPIC guideline defining CYP2D6 activity score framework and allele-function assignments.",
        ),
    ),
    "CYP2C9": _activity_score_meta(
        "Guideline-backed pharmacogene",
        "Core pharmacogene with PharmVar-backed allele definitions and CPIC activity-score translation.",
        _clinpgx_ref("CYP2C9"),
        _pharmvar_ref("CYP2C9"),
        _cpic_ref(
            "CYP2C9-WARFARIN-2017",
            "Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for Pharmacogenetics-Guided Warfarin Dosing: 2017 Update.",
        ),
        _pubmed_ref(
            "28198005",
            "Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for Pharmacogenetics-Guided Warfarin Dosing: 2017 Update.",
            "CPIC guideline defining CYP2C9 allele-function and activity score assignments.",
        ),
    ),
    "CYP2C19": _activity_score_meta(
        "Guideline-backed pharmacogene",
        "Core pharmacogene with PharmVar-backed allele definitions and CPIC metabolizer-phenotype translation.",
        _clinpgx_ref("CYP2C19"),
        _pharmvar_ref("CYP2C19"),
        _cpic_ref(
            "CYP2C19-PPI",
            "Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for CYP2C19 and Proton Pump Inhibitor Dosing.",
        ),
        _pubmed_ref(
            "29385237",
            "Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for CYP2C19 and Proton Pump Inhibitor Dosing.",
            "CPIC guideline defining CYP2C19 allele-function and metabolizer phenotype assignments.",
        ),
    ),
    "NAT2": _activity_score_meta(
        "Phenotype-backed pharmacogene",
        "Rapid/intermediate/slow acetylator labels are retained for research-use summarization and should not replace diplotype calling.",
        _clinpgx_ref("NAT2"),
        _ncbi_gene_ref("NAT2"),
        _pubmed_ref(
            "30149019",
            "Expression and genotype-dependent catalytic activity of N-acetyltransferase 2 (NAT2) in human peripheral blood mononuclear cells and its modulation by Sirtuin 1.",
            "Supports NAT2 rapid/intermediate/slow acetylator phenotype interpretation.",
        ),
    ),
    "UGT1A1": _activity_score_meta(
        "Guideline-backed pharmacogene",
        "UGT1A1 activity groupings are anchored to implementation-facing gene resources with explicit CPIC phenotype translation.",
        _clinpgx_ref("UGT1A1"),
        _ncbi_gene_ref("UGT1A1"),
        _cpic_ref(
            "UGT1A1-ATAZANAVIR",
            "Clinical Pharmacogenetics Implementation Consortium (CPIC) Guidelines for UGT1A1 and Atazanavir Prescribing.",
        ),
        _pubmed_ref(
            "24296998",
            "Clinical Pharmacogenetics Implementation Consortium (CPIC) Guidelines for UGT1A1 and Atazanavir Prescribing.",
            "CPIC guideline defining UGT1A1 allele-function and activity groupings.",
        ),
    ),
    "CYP2A6": _activity_score_meta(
        "Database-backed pharmacogene",
        "Allele-function groupings are supported by pharmacogene resources and smoking-metabolism literature.",
        _clinpgx_ref("CYP2A6"),
        _pharmvar_ref("CYP2A6"),
        _pubmed_ref(
            "29194389",
            "Variation in CYP2A6 Activity and Personalized Medicine.",
            "Supports CYP2A6 metabolizer-group interpretation in nicotine metabolism.",
        ),
    ),
    "CYP3A4": _activity_score_meta(
        "Database-backed pharmacogene",
        "CYP3A4 activity groupings are simplified from pharmacogene resources and should be treated as research-use summaries.",
        _clinpgx_ref("CYP3A4"),
        _pharmvar_ref("CYP3A4"),
        _pubmed_ref(
            "23665933",
            "CYP3A4 intron 6 C>T polymorphism (CYP3A4*22) is associated with reduced CYP3A4 protein level and function in human liver microsomes.",
            "Supports reduced-function interpretation of CYP3A4*22.",
        ),
    ),
    "CYP1A1": _activity_score_meta(
        "Research-use literature-derived",
        "Current values are heuristic activity summaries for carcinogen-metabolism curation, not CPIC-standardized scores.",
        _ncbi_gene_ref("CYP1A1"),
        _pubmed_ref(
            "15647817",
            "Effect of CYP1A1 gene polymorphisms on estrogen metabolism and bone density.",
            "Supports functional interpretation of common CYP1A1 polymorphisms including Ile462Val-containing haplotypes.",
        ),
    ),
    "CYP1A2": _activity_score_meta(
        "Research-use literature-derived",
        "Current values are heuristic inducibility/activity summaries and should be used only for research curation.",
        _ncbi_gene_ref("CYP1A2"),
        _pubmed_ref(
            "16188490",
            "Influence of the genetic polymorphism in the 5'-noncoding region of the CYP1A2 gene on CYP1A2 phenotype and urinary mutagenicity in smokers.",
            "Supports inducibility interpretation for CYP1A2 promoter polymorphisms.",
        ),
    ),
    "CYP1B1": _activity_score_meta(
        "Research-use literature-derived",
        "Current values summarize literature-reported differences in estradiol hydroxylation rather than a standardized clinical score.",
        _ncbi_gene_ref("CYP1B1"),
        _pubmed_ref(
            "10862525",
            "Polymorphisms in P450 CYP1B1 affect the conversion of estradiol to the potentially carcinogenic metabolite 4-hydroxyestradiol.",
            "Supports functional interpretation of CYP1B1 codon 432 variation.",
        ),
    ),
    "CYP2E1": _activity_score_meta(
        "Research-use literature-derived",
        "CYP2E1 functional polymorphism evidence is mixed; these values are retained as tentative research-use annotations.",
        _ncbi_gene_ref("CYP2E1"),
        _pubmed_ref(
            "10886461",
            "Lack of evidence for a role of cytochrome P450 2E1 genetic polymorphisms in the development of different types of cancer.",
            "Provides a review context showing that CYP2E1 polymorphism effects remain debated and should be interpreted cautiously.",
        ),
    ),
    "EPHX1": _activity_score_meta(
        "Research-use literature-derived",
        "EPHX1 activity values reflect long-used low/intermediate/high genotype conventions rather than a formal clinical standard.",
        _ncbi_gene_ref("EPHX1"),
        _pubmed_ref(
            "21445251",
            "Putative EPHX1 enzyme activity is related with risk of lung and upper aerodigestive tract cancers: a comprehensive meta-analysis.",
            "Supports low-activity and high-activity genotype grouping conventions for Tyr113His and His139Arg.",
        ),
    ),
    "GSTM1": _activity_score_meta(
        "Research-use loss-of-function genotype",
        "The GSTM1 null mapping represents a deletion/absence-of-function convention rather than a CPIC activity score.",
        _ncbi_gene_ref("GSTM1"),
        _pubmed_ref(
            "23506349",
            "Utilization of glutathione S-transferase Mu 1- and Theta 1-null mice as animal models for absorption, distribution, metabolism, excretion and toxicity studies.",
            "Supports absent-function interpretation for GSTM1 null status.",
        ),
    ),
    "GSTT1": _activity_score_meta(
        "Research-use loss-of-function genotype",
        "The GSTT1 null mapping represents a deletion/absence-of-function convention rather than a CPIC activity score.",
        _ncbi_gene_ref("GSTT1"),
        _pubmed_ref(
            "23506349",
            "Utilization of glutathione S-transferase Mu 1- and Theta 1-null mice as animal models for absorption, distribution, metabolism, excretion and toxicity studies.",
            "Supports absent-function interpretation for GSTT1 null status.",
        ),
    ),
    "NQO1": _activity_score_meta(
        "Research-use literature-derived",
        "NQO1 values summarize the functional impact of the common Pro187Ser variant for curation purposes.",
        _ncbi_gene_ref("NQO1"),
        _pubmed_ref(
            "23860519",
            "The NQO1 polymorphism C609T (Pro187Ser) and cancer susceptibility: a comprehensive meta-analysis.",
            "Provides consolidated literature context for the reduced-function NQO1*2 allele.",
        ),
    ),
    "XPC": _activity_score_meta(
        "Research-use literature-derived",
        "XPC repair-capacity values are heuristic and should be treated as exploratory annotations rather than clinical scores.",
        _ncbi_gene_ref("XPC"),
        _ctd_ref("XPC", "DNA Repair (NER)"),
        _pubmed_ref(
            "22592359",
            "XPC Lys939Gln polymorphism contributes to colorectal cancer susceptibility: evidence from a meta-analysis.",
            "Meta-analysis supporting functional interpretation of the XPC Lys939Gln variant on NER capacity.",
        ),
    ),
    "OGG1": _activity_score_meta(
        "Research-use literature-derived",
        "OGG1 values summarize the Ser326Cys literature and should be used only for research interpretation.",
        _ncbi_gene_ref("OGG1"),
        _pubmed_ref(
            "25588927",
            "OGG1 Ser326Cys polymorphism and cancer risk: a meta-analysis of 27 published studies.",
            "Meta-analysis providing consolidated support for functional interpretation of the OGG1 Ser326Cys variant.",
        ),
    ),
    "MGMT": _activity_score_meta(
        "Research-use epigenetic marker",
        "MGMT promoter methylation states are represented as simplified repair-capacity categories for research use.",
        _ncbi_gene_ref("MGMT"),
        _pubmed_ref(
            "20725792",
            "MGMT promoter methylation in malignant gliomas.",
            "Supports the reduced-expression and silencing interpretation of MGMT promoter methylation states.",
        ),
    ),
    "MLH1": _activity_score_meta(
        "Research-use mismatch repair marker",
        "MLH1 activity categories reflect promoter methylation and pathogenic variant status relevant to Lynch syndrome and post-replicative repair at adduct sites.",
        _ncbi_gene_ref("MLH1"),
        _pubmed_ref(
            "25559809",
            "Revised guidelines for the clinical management of Lynch syndrome.",
            "Supports MLH1 germline variant and promoter methylation interpretation for mismatch repair deficiency.",
        ),
    ),
    "MSH2": _activity_score_meta(
        "Research-use mismatch repair marker",
        "MSH2 activity categories reflect pathogenic variant status relevant to Lynch syndrome and mismatch recognition at replication errors past DNA adducts.",
        _ncbi_gene_ref("MSH2"),
        _pubmed_ref(
            "25559809",
            "Revised guidelines for the clinical management of Lynch syndrome.",
            "Supports MSH2 germline variant interpretation for mismatch repair deficiency.",
        ),
    ),
    # Wave 2 class-specific activity score metadata
    "ALDH2": _activity_score_meta(
        "Guideline-backed pharmacogene",
        "ALDH2*2 (rs671) is the most clinically significant variant; near-complete loss of enzymatic activity causes 6-19x acetaldehyde accumulation after ethanol ingestion.",
        _ncbi_gene_ref("ALDH2"),
        _clinpgx_ref("ALDH2"),
        _pubmed_ref(
            "PMC7941978",
            "ALDH2 rs671 polymorphism and esophageal cancer risk: a meta-analysis.",
            "Meta-analysis supporting ALDH2*2 as a major risk factor for esophageal squamous cell carcinoma in alcohol consumers.",
        ),
    ),
    "ADH1B": _activity_score_meta(
        "Research-use literature-derived",
        "ADH1B*2 (rs1229984, His47Arg) encodes a ~40x faster enzyme; high prevalence in East Asian populations (70-80%).",
        _ncbi_gene_ref("ADH1B"),
        _clinpgx_ref("ADH1B"),
        _pubmed_ref(
            "PMC4122263",
            "ADH1B Arg47His polymorphism and cancer risk.",
            "Supports functional interpretation of ADH1B*2 rapid-metabolizer allele in ethanol-related cancer risk.",
        ),
    ),
}


def build_tier1_panel() -> KnowledgeGraph:
    """Return a KnowledgeGraph containing all 13 Tier 1 genes as Enzyme nodes."""
    nodes = [
        Node(id=g["id"], label=g["label"], type=NodeType.ENZYME, tier=1, **{
            k: v for k, v in g.items() if k not in ("id", "label")
        })
        for g in TIER1_GENES
    ]
    return KnowledgeGraph(nodes=nodes, edges=[])


def build_tier2_panel() -> KnowledgeGraph:
    """Return a KnowledgeGraph containing all 25 Tier 2 genes as Enzyme nodes."""
    nodes = [
        Node(id=g["id"], label=g["label"], type=NodeType.ENZYME, tier=2, **{
            k: v for k, v in g.items() if k not in ("id", "label")
        })
        for g in TIER2_GENES
    ]
    return KnowledgeGraph(nodes=nodes, edges=[])


def build_wave2_panel() -> KnowledgeGraph:
    """Return a KnowledgeGraph containing all 9 Wave 2 class-specific genes as Enzyme nodes."""
    nodes = [
        Node(id=g["id"], label=g["label"], type=NodeType.ENZYME, tier=3, **{
            k: v for k, v in g.items() if k not in ("id", "label")
        })
        for g in WAVE2_GENES
    ]
    return KnowledgeGraph(nodes=nodes, edges=[])


def build_full_panel(*, include_wave2: bool = False) -> KnowledgeGraph:
    """Return a KnowledgeGraph with all Tier 1 + Tier 2 genes (38 total).

    Args:
        include_wave2: If ``True``, also include the 9 Wave 2 class-specific
            genes (tier 3). Defaults to ``False`` for backward compatibility.
    """
    t1 = build_tier1_panel()
    t2 = build_tier2_panel()
    nodes = t1.nodes + t2.nodes
    if include_wave2:
        nodes = nodes + build_wave2_panel().nodes
    return KnowledgeGraph(nodes=nodes, edges=[])


def _reference_graph_data_path() -> Path:
    """Return the canonical JSON source of the bundled reference graph.

    ``graph-data.js`` remains checked in alongside this file purely as the
    asset consumed by the bundled Streamlit/D3 viewer (see ``ui_map_viewer``
    and ``exporter.to_graph_data_js``); all Python-side reference graph
    loading now reads the JSON duplicate instead.
    """
    return Path(__file__).resolve().parent / "map" / "graph-data.json"


def build_reference_graph() -> KnowledgeGraph:
    """Return the bundled full reference graph as a package-native model."""
    from .exporter import parse_graph_artifact

    return parse_graph_artifact(_reference_graph_data_path())


def build_reference_engine() -> "GraphEngine":
    """Return a GraphEngine loaded with the bundled full reference graph."""
    from .engine import GraphEngine

    engine = GraphEngine()
    engine.load(build_reference_graph())
    return engine


def get_activity_scores(gene: str) -> list[dict[str, Any]] | None:
    """Look up activity score entries for a gene. Returns None if not found."""
    return ACTIVITY_SCORES.get(gene)


def get_activity_score_metadata(gene: str) -> dict[str, object] | None:
    """Return evidence metadata for a gene's activity score table."""
    return ACTIVITY_SCORE_METADATA.get(gene)


def get_activity_score_references(gene: str) -> list[dict[str, Any]] | None:
    """Return supporting references for a gene's activity score table."""
    metadata = get_activity_score_metadata(gene)
    if metadata is None:
        return None
    return metadata.get("references")  # type: ignore[return-value]
