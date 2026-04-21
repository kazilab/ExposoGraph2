Gene Panels & Activity Scores
=============================

ExposoGraph ships with curated reference gene panels and CPIC-standardized
activity scores from the CarcinoGenomic Platform (Tables 3–6).

Each seed gene now includes structured database provenance, and the package
exposes a formal source manifest through
:data:`~ExposoGraph.reference_data.CURATION_SOURCE_MANIFEST`.
The manuscript-aligned primary curation sources are IARC, KEGG, PharmVar,
CPIC, CTD, GTEx, and PubMed. ExposoGraph also uses NCBI Gene and ClinPGx as
supporting implementation sources for stable identifiers and pharmacogene
coverage.

The curated KEGG pathway catalog exposed through
:data:`~ExposoGraph.reference_data.REFERENCE_KEGG_PATHWAYS` currently tracks
``hsa00980`` (xenobiotic metabolism by cytochrome P450), ``hsa00140``
(steroid hormone biosynthesis), ``hsa05204`` (chemical carcinogenesis - DNA
adducts), and ``hsa05208`` (chemical carcinogenesis - reactive oxygen species).

Terminology:

- ``Tier`` refers to panel priority or inclusion level.
- ``Phase`` is reserved for Phase I, II, and III metabolism/transport labels.
- DNA repair genes are classified by repair class such as BER, NER, or Direct Reversal.

Tier 1: Core Panel (13 genes)
------------------------------

The core carcinogen-metabolizing enzyme panel:

.. list-table::
   :header-rows: 1
   :widths: 12 8 15 65

   * - Gene
     - Classification
     - Role
     - Detail
   * - CYP1A1
     - Phase I
     - Activation
     - PAH diol-epoxide formation; AhR-inducible; extrahepatic expression
   * - CYP1A2
     - Phase I
     - Activation
     - HCA and aromatic amine activation; AhR-inducible; hepatic
   * - CYP1B1
     - Phase I
     - Activation
     - PAH and estrogen activation; 4-OH-estradiol formation
   * - CYP2A6
     - Phase I
     - Activation
     - NNK and nitrosamine activation; nicotine metabolism
   * - CYP2E1
     - Phase I
     - Activation
     - Small-molecule carcinogen activation; ethanol, benzene, NDMA
   * - CYP3A4
     - Phase I
     - Activation
     - Aflatoxin B1 8,9-epoxidation; broad substrate range
   * - GSTM1
     - Phase II
     - Detoxification
     - GSH conjugation of PAH diol-epoxides; null polymorphism common
   * - GSTT1
     - Phase II
     - Detoxification
     - GSH conjugation of small electrophiles; null polymorphism common
   * - GSTP1
     - Phase II
     - Detoxification
     - GSH conjugation of BPDE and PAH diol-epoxides
   * - NAT2
     - Phase II
     - Mixed
     - N-acetylation of aromatic amines; rapid vs slow acetylator
   * - EPHX1
     - Phase II
     - Mixed
     - Microsomal epoxide hydrolysis; activation and detoxification
   * - UGT1A1
     - Phase II
     - Detoxification
     - Glucuronidation of PAH metabolites and bilirubin
   * - NQO1
     - Phase II
     - Detoxification
     - Two-electron quinone reduction; prevents ROS from redox cycling

Tier 2: Extended Panel (23 genes)
----------------------------------

Additional Phase I, II, III, and DNA repair genes:

.. list-table::
   :header-rows: 1
   :widths: 12 8 15 65

   * - Gene
     - Classification
     - Role
     - Detail
   * - SULT1A1
     - Phase II
     - Detoxification
     - Sulfation of phenolic compounds, HCA intermediates
   * - NAT1
     - Phase II
     - Mixed
     - O-acetylation of aromatic/heterocyclic amines in peripheral tissues
   * - UGT2B7
     - Phase II
     - Detoxification
     - Glucuronidation of steroid hormones, carcinogen metabolites
   * - UGT2B17
     - Phase II
     - Detoxification
     - Glucuronidation of testosterone, DHT, and related androgen metabolites
   * - SULT1E1
     - Phase II
     - Detoxification
     - High-affinity estrogen sulfotransferase; estradiol and catechol-estrogen inactivation
   * - COMT
     - Phase II
     - Detoxification
     - O-methylation of catechol estrogens; limits redox-cycling estrogen metabolites
   * - ABCB1
     - Phase III
     - Transport
     - P-gp; ATP-driven efflux of hydrophobic xenobiotics
   * - ABCC2
     - Phase III
     - Transport
     - MRP2; export of GSH/glucuronide/sulfate conjugates
   * - ABCG2
     - Phase III
     - Transport
     - BCRP; efflux of PAHs, PhIP, porphyrins
   * - XRCC1
     - DNA Repair (BER)
     - Repair
     - BER scaffold protein; oxidative/alkylation DNA damage
   * - XPC
     - DNA Repair (NER)
     - Repair
     - GG-NER damage sensor; bulky DNA adducts from PAHs
   * - ERCC2/XPD
     - DNA Repair (NER)
     - Repair
     - NER helicase; unwinds DNA at damage sites
   * - OGG1
     - DNA Repair (BER)
     - Repair
     - 8-oxoguanine DNA glycosylase; oxidative DNA damage
   * - MGMT
     - DNA Repair (Direct Reversal)
     - Repair
     - Direct reversal of O6-alkylguanine; suicidal repair enzyme
   * - CYP2C9
     - Phase I
     - Mixed
     - Oxidation of PAH metabolites; major drug-metabolizing CYP
   * - CYP2C19
     - Phase I
     - Mixed
     - Minor procarcinogen activation; nitrosamine metabolism
   * - CYP2D6
     - Phase I
     - Mixed
     - Minor NNK activation; dual PGx/carcinogen-risk reporting
   * - CYP2A13
     - Phase I
     - Activation
     - Primary lung NNK-metabolizing CYP; tobacco-smoke activation
   * - CYP17A1
     - Phase I
     - Activation
     - Steroid 17alpha-hydroxylase/17,20-lyase; androgen precursor synthesis
   * - SRD5A1
     - Phase I
     - Activation
     - 5alpha-reductase type 1; peripheral testosterone-to-DHT conversion
   * - SRD5A2
     - Phase I
     - Activation
     - 5alpha-reductase type 2; major DHT-forming enzyme in prostate tissues
   * - CYP19A1
     - Phase I
     - Activation
     - Aromatase; converts androgen precursors to estrogens
   * - AKR1C3
     - Phase I
     - Mixed
     - Local androgen/estrogen activation with quinone-reductase overlap

Activity Scores
---------------

ExposoGraph currently ships activity-score tables for 18 genes. These tables
mix two evidence classes:

- guideline-backed pharmacogene resources, primarily surfaced through ClinPGx/PharmVar
- research-use literature-derived mappings for carcinogen metabolism and DNA repair genes

They should therefore be treated as referenced curation aids, not as a full
clinical pharmacogenomics engine. Each gene now has supporting evidence
metadata available through
:data:`~ExposoGraph.reference_data.ACTIVITY_SCORE_METADATA`.

Each per-allele entry has:

- **allele** — Star allele or variant name
- **value** — Numeric activity score (0.0 = no function, 1.0 = normal, 2.0 = ultrarapid)
- **phenotype** — Functional interpretation
- **confidence** — High or Moderate

Example: CYP2D6
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 10 30 15

   * - Allele
     - Score
     - Phenotype
     - Confidence
   * - \*1, \*2, \*35
     - 1.0
     - Normal Metabolizer
     - High
   * - \*9, \*17, \*29, \*41
     - 0.5
     - NM (lower range)
     - High
   * - \*10
     - 0.25
     - NM (lowest); IM
     - High
   * - \*3, \*4, \*5, \*6, \*40
     - 0.0
     - Poor Metabolizer
     - High
   * - \*1x2, \*2x2 (dupl.)
     - 2.0
     - Ultrarapid Metabolizer
     - High

The full set covers 18 genes. See
:data:`~ExposoGraph.reference_data.ACTIVITY_SCORES` for the score table and
:data:`~ExposoGraph.reference_data.ACTIVITY_SCORE_METADATA` for the references.

Usage
^^^^^

.. code-block:: python

   from ExposoGraph import (
       build_full_panel,
       get_activity_scores,
       get_activity_score_references,
   )

   # Load all 36 genes into a KnowledgeGraph
   kg = build_full_panel()

   # Look up scores for a specific gene
   scores = get_activity_scores("CYP1A1")
   for entry in scores:
       print(f"{entry['allele']}: {entry['value']} — {entry['phenotype']}")

   # Inspect the supporting references
   for ref in get_activity_score_references("CYP1A1") or []:
       print(ref["source_db"], ref.get("pmid") or ref.get("record_id"), ref["url"])
