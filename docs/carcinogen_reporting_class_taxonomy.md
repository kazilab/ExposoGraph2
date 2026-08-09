# Carcinogen Reporting Class Taxonomy

**Branch:** `feature_extend_knowledge_graph`
**Source data:** `ExposoGraph/map/graph-data.js` (commit `d6c373f`)
**Date:** 2026-08-09

---

## 1. Overview

This document presents a two-dimensional taxonomy of the 21 carcinogen reporting classes in the ExposoGraph knowledge graph. The two dimensions are:

1. **IARC Agent Formalism** — the real-world nature and exposure source of the agent, following the IARC Monographs Preamble agent-type categories ([IARC Preamble, 2019](https://monographs.iarc.who.int/wp-content/uploads/2019/01/Preamble-2019.pdf)).
2. **Mechanistic Processing** — what the body does to the agent before it can exert carcinogenic damage, derived from the enzyme-metabolite-adduct chains encoded in the graph edges.

A third section tabulates all 66 individual carcinogens against the IARC [Key Characteristics of Carcinogens (KCC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4898322/) framework (Smith et al., *Env Health Perspect*, 2016), demonstrating that KCC is orthogonal to both the reporting classes and the mechanistic processing tiers.

The core thesis: **KCC describes cellular effects (what the agent does to the cell); ExposoGraph's reporting classes describe exposure identity; ExposoGraph's mechanistic tiers describe bodily processing chains. The three frameworks are intentionally orthogonal rather than one-to-one.**

---

## 2. Enumeration of 21 Reporting Classes

| # | Group | Members | Count | IARC Group(s) |
|---|-------|---------|-------|----------------|
| 1 | Alcohol | Ethanol, Urethane | 2 | 1, 2A |
| 2 | Aldehyde | Formaldehyde, Acrolein, Crotonaldehyde, Furfural, MDA | 5 | 1, 2A, 2B, 3, 3 |
| 3 | Alkylating | EthyleneOxide, Acrylamide, Cyclophosphamide, Chlorambucil, Sulfur_mustard, Busulfan, MNU, Temozolomide | 8 | 1, 2A, 1, 1, 1, 1, 2A, — |
| 4 | Androgen | AndrogenicAnabolicSteroids | 1 | 2A |
| 5 | Aromatic_Amine | 4ABP, Benzidine | 2 | 1, 1 |
| 6 | Benzene | Benzene | 1 | 1 |
| 7 | Chlorinated_Solvent | TCE, PCE | 2 | 1, 2A |
| 8 | Dioxin | TCDD, PeCDF_23478 | 2 | 1, 1 |
| 9 | Estrogen | EstrogenProgestogenTherapy | 1 | 1 |
| 10 | HCA | PhIP, MeIQx | 2 | 2B, 2B |
| 11 | Heavy_Metal | ArsenicInorganic, Cd, CrVI, NickelCompounds, NickelMetallic, Be, LeadInorganicCompounds, LeadMetallic, LeadOrganicCompounds, MethylmercuryCompounds, MercuryInorganicCompounds, CobaltMetal, CobaltOxide, AntimonyTrivalent, AntimonyPentavalent | 15 | 1, 1, 1, 1, 2B, 1, 2A, 2B, 3, 2B, 3, 2A, 2B, 2A, 3 |
| 12 | Mycotoxin | AFB1 | 1 | 1 |
| 13 | Nitrosamine | NNK, NDMA, NDEA | 3 | 1, 2A, 2A |
| 14 | Organochlorine | HCB, Lindane, DDT, DDE, PCP, Chlordane, Heptachlor, Toxaphene | 8 | 2B, 1, 2A, 2B, 1, 2B, 2B, 2B |
| 15 | PAH | BaP, DMBA | 2 | 1, — |
| 16 | PCB | PCB_126, PCB_169, PCB_77, PCB_118, PCB_153, PCB_138 | 6 | 1×4, —, — |
| 17 | PFAS | PFOA, PFOS | 2 | 1, 2B |
| 18 | Plant_Alkaloid | AristolochicAcid | 1 | 1 |
| 19 | UV_Radiation | UVRadiation | 1 | 1 |
| 20 | Ionizing_Radiation | IonizingRadiation, Radon | 2 | 1, 1 |
| 21 | Vinyl Chloride | VinylChloride | 1 | 1 |

**Notes on demoted nodes (now Metabolites, not Carcinogens):**
- **Acetaldehyde**: now a Metabolite (demoted from Carcinogen in the Aldehyde class). Intracellular metabolite of ethanol. Not separately IARC-evaluated.
- **Glycidamide**: now a Metabolite (demoted from Carcinogen in the Alkylating class). Intracellular metabolite of acrylamide. Not separately IARC-evaluated.
- **4-HNE**: now a Metabolite (demoted from Carcinogen in the Aldehyde class). Endogenous lipid peroxidation product. Not IARC-evaluated.
- **MDA**: IARC Group 3 (Vol. 71, 1999). Endogenous genotoxin but IARC-evaluated.
- **E2 and 4-OHE2**: now Metabolites (demoted from Carcinogen). Endogenous estrogen metabolites. Not separately IARC-evaluated.
- **DHT**: now a Metabolite (demoted from Carcinogen). Endogenous androgen metabolite. Not separately IARC-evaluated.
- **Temozolomide**: Not IARC-evaluated. Chemotherapy drug.
- **DMBA**: experimental carcinogen, not a human IARC agent.
- **DDE**: evaluated within the DDT monograph (IARC Vol. 53, 1991). DDT upgraded to Group 2A in Vol. 113 (2018) but DDE not separately re-evaluated.

---

## 3. Dimension 1: IARC Agent Formalism

IARC's Preamble classifies agents by their real-world nature and exposure source — what kind of thing is it, and where does it come from? The monograph volumes themselves are organized by agent type (Vol. 100A = Pharmaceuticals, Vol. 100C = Metals/particles/fibers, Vol. 100D = Radiation, Vol. 100E = Lifestyle factors, Vol. 100F = Chemical agents and related occupations).

| IARC Agent Type | ExposoGraph Groups | Rationale |
|---|---|---|
| Industrial/Environmental Chemicals | Benzene, Chlorinated_Solvent, Dioxin, Organochlorine, PCB, PFAS, Vinyl Chloride | Synthetic compounds from industrial manufacture, environmental contamination, or persistent organic pollutants |
| Metals & Elemental Substances | Heavy_Metal | Inorganic elements; IARC Vol. 100C |
| Combustion/Pyrolysis Products | PAH, HCA | Formed during incomplete combustion (PAH) or high-temperature cooking of food (HCA) |
| Natural Toxins & Food Contaminants | Mycotoxin, Plant_Alkaloid, Nitrosamine | Naturally occurring (fungal metabolites, plant constituents) or formed in preserved foods/tobacco |
| Pharmaceuticals & Hormones | Androgen, Estrogen, Alkylating (partial) | Steroid hormones and oncology drugs (cyclophosphamide, busulfan, temozolomide, chlorambucil) |
| Reactive Chemical Classes | Aldehyde, Aromatic_Amine | Functional chemistry groupings that span multiple exposure sources |
| Lifestyle/Habit | Alcohol | IARC Vol. 100E classifies alcoholic beverages as a personal habit |
| Physical Agents (Non-ionizing) | UV_Radiation | Energy transfer via photons; IARC Vol. 100D |
| Physical Agents (Ionizing) | Ionizing_Radiation | Particle/radiation energy transfer; IARC Vol. 100D |

Several groups do not map cleanly to a single IARC category. **Alkylating** spans pharmaceuticals (cyclophosphamide, busulfan) and military/industrial chemicals (sulfur mustard, ethylene oxide). **Nitrosamine** spans tobacco-specific (NNK) and dietary/industrial (NDMA, NDEA). This is inherent to IARC's source-based formalism — it classifies by where the agent comes from, not by how it behaves in the body.

---

## 4. Dimension 2: Mechanistic Processing by the Body

This is the dimension that distinguishes ExposoGraph from IARC. The question here is: what must the body do to the agent before it can exert carcinogenic damage? This is orthogonal to IARC's agent type and to the KCC formalism.

The KCC describes what the agent **does to the cell** (10 characteristics). ExposoGraph's mechanistic processing describes what the **body does to the agent** — the processing chain from exposure to damage. Two agents can share the same KCC (e.g., both "is genotoxic") but require completely different bodily processing.

### Tier 1: Direct-acting (no metabolic activation required)

The agent is already the ultimate carcinogen. DNA damage or cellular disruption occurs without enzymatic bioactivation.

| Group | Mechanism | Graph evidence |
|---|---|---|
| Alkylating | Direct electrophilic attack on DNA bases | Most members FORMS_ADDUCT directly; no metabolite intermediary (except cyclophosphamide → CYP3A4) |
| Aldehyde | Direct DNA crosslinking and adduct formation | Formaldehyde FORMS_ADDUCT to N2_ethylidene_dG; no metabolite node in chain |
| Heavy_Metal | Direct oxidative stress, enzyme inhibition, metal-DNA binding | CrVI/As/Cd FORMS_ADDUCT to Oxo_dG; no CYP activation step |
| UV_Radiation | Photodimer formation (CPD, 6-4PP) | UVRadiation FORMS_ADDUCT to CPD/64PP — no enzyme intermediary |
| Ionizing_Radiation | DSB/SSB via radiolysis | IonizingRadiation FORMS_ADDUCT to DSB/SSB; Radon alpha particles — no enzyme intermediary |

### Tier 2: Metabolic bioactivation required (procarcinogen → ultimate carcinogen)

The parent compound is relatively inert. The body enzymatically transforms it (primarily via CYP450, sometimes ADH) into a reactive electrophile that binds DNA.

| Group | Activating enzyme | Bioactivation chain in graph |
|---|---|---|
| PAH | CYP1A1 | BaP → BaP_epoxide → BaP_diol → BPDE → BPDE_dG |
| Aromatic_Amine | CYP1A2 → NAT | 4ABP → NOH_4ABP → ABP_dG |
| HCA | CYP1A2 → NAT | PhIP → NOH_PhIP → PhIP_dG |
| Nitrosamine | CYP2A6 / CYP2A13 | NNK → NNK_hydroxyl → POB_dG |
| Mycotoxin | CYP3A4 | AFB1 → AFB1_epoxide → AFB1_Gua |
| Plant_Alkaloid | CYP1A1 / CYP1A2 | AristolochicAcid → AL_nitrenium → dA_AL_I |
| Benzene | CYP2E1 | Benzene → Benzene_oxide → HQ → Benzoquinone → BQ_dG |
| Vinyl Chloride | CYP2E1 | VinylChloride → CEO (implied) → DNA adduct |
| Chlorinated_Solvent | CYP2E1 | TCE/PCE → reactive metabolites → Oxo_dG |
| Alcohol | ADH (not CYP) | Ethanol → Acetaldehyde → N2_ethylidene_dG |

### Tier 3: Receptor-mediated (non-genotoxic)

The agent acts through binding to cellular receptors — nuclear receptors, AHR, hormone receptors. No DNA-reactive metabolite is formed. The mechanism is altered gene expression, cell proliferation, or immunosuppression.

| Group | Primary receptor | Edge type in graph |
|---|---|---|
| Dioxin | AHR | TCDD ACTIVATES AHR; INDUCES CYP1A1 |
| PFAS | PPARα, CAR, PXR | PFOA/PFOS ACTIVATES PPARA |
| Estrogen | ESR1 (ERα) | EstrogenProgestogenTherapy ACTIVATES ESR1 |
| Androgen | AR | AndrogenicAnabolicSteroids ACTIVATES AR |

### Tier 4: Mixed mechanism (dual processing pathways)

These groups have evidence for more than one processing tier. The graph captures both pathways as separate edges.

| Group | Tier 2 pathway | Tier 3 pathway | Notes |
|---|---|---|---|
| PCB | CYP-mediated quinone metabolites → DNA adducts (non-dioxin-like) | AHR activation (dioxin-like: PCB_126, PCB_169, PCB_77) | Congener-dependent mechanism |
| Organochlorine | Oxidative stress → Oxo_dG (several members) | AHR/ER mediated (DDT, lindane) | PCP is direct-acting (mitochondrial uncoupler) |
| Alkylating (partial) | Cyclophosphamide requires CYP3A4 activation | — | Most members are Tier 1 direct-acting |

---

## 5. Cross-Tabulation: IARC Agent Type × Mechanistic Tier

| | Tier 1: Direct | Tier 2: Bioactivation | Tier 3: Receptor-mediated | Tier 4: Mixed |
|---|---|---|---|---|
| Industrial/Environmental Chemicals | — | Benzene, Chlorinated_Solvent, Vinyl Chloride | Dioxin, PFAS | PCB, Organochlorine |
| Metals & Elemental | Heavy_Metal | — | — | — |
| Combustion/Pyrolysis Products | — | PAH, HCA | — | — |
| Natural Toxins & Food Contaminants | — | Mycotoxin, Plant_Alkaloid, Nitrosamine | — | — |
| Pharmaceuticals & Hormones | — | Alkylating (partial) | Androgen, Estrogen | — |
| Reactive Chemical Classes | Aldehyde, Alkylating | Aromatic_Amine | — | — |
| Lifestyle/Habit | — | Alcohol | — | — |
| Physical Agents (Non-ionizing) | UV_Radiation | — | — | — |
| Physical Agents (Ionizing) | Ionizing_Radiation | — | — | — |

---

## 6. How This Distinguishes ExposoGraph from IARC/KCC

### Three key distinctions

**1. Orthogonal axes.** IARC's agent formalism asks "what kind of thing is it?" The KCC asks "what does it do to the cell?" ExposoGraph's mechanistic processing asks "what does the body do to it first?" These are independent questions. PAH and Aromatic_Amine are in different IARC categories (combustion product vs. reactive chemical class) but share the same mechanistic processing (CYP450 bioactivation). Conversely, Aldehyde and Alkylating are both direct-acting but come from entirely different exposure sources.

**2. Processing chain visibility.** The KCC is a flat list of characteristics — an agent either "is genotoxic" or isn't. ExposoGraph captures the chain: procarcinogen → activating enzyme → reactive metabolite → DNA adduct. This chain identifies intervention points (e.g., CYP2E1 polymorphisms affect benzene, VC, and TCE risk simultaneously; NAT2 slow acetylators affect aromatic amine and HCA risk).

**3. Shared enzyme infrastructure.** The cross-tabulation reveals that multiple IARC agent categories share the same CYP enzymes: CYP1A1 (PAH, plant alkaloid), CYP1A2 (aromatic amine, HCA), CYP2E1 (benzene, vinyl chloride, chlorinated solvent), CYP3A4 (mycotoxin, alkylating-pharma). IARC's source-based grouping obscures these shared metabolic dependencies; ExposoGraph's mechanistic grouping makes them visible because the enzyme nodes are shared across groups in the graph.

---

## 7. KCC Tabulation: Individual Carcinogens × 10 Key Characteristics

### KCC definitions

The 10 Key Characteristics of Carcinogens, as defined by [Smith et al., 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4898322/) and adopted in the [IARC Preamble (2019)](https://monographs.iarc.who.int/wp-content/uploads/2019/01/Preamble-2019.pdf):

| KC | Characteristic |
|----|----------------|
| KC1 | Is electrophilic or can be metabolically activated |
| KC2 | Is genotoxic |
| KC3 | Alters DNA repair or causes genomic instability |
| KC4 | Induces epigenetic alterations |
| KC5 | Induces oxidative stress |
| KC6 | Induces chronic inflammation |
| KC7 | Is immunosuppressive |
| KC8 | Modulates receptor-mediated effects |
| KC9 | Causes immortalization |
| KC10 | Alters cell proliferation, cell death, or nutrient supply |

### Rating scheme

| Symbol | Meaning |
|--------|---------|
| ● | Strong/sufficient evidence from IARC, council synthesis, or graph edges |
| ◐ | Moderate/limited evidence; class-supported but not cleanly established for this node |
| — | Not assessed or no clear support in current literature |

### Basis column

| Label | Meaning |
|-------|---------|
| graph edge | Directly supported by ACTIVATES/FORMS_ADDUCT/INDUCES edges in graph-data.js |
| IARC KC | Supported by IARC Monograph KC evaluation or Preamble |
| council | Supported by ExposoGraph Model Council synthesis (Radiation or PFAS) |
| literature | Supported by published toxicology literature, class-level evidence |

---

### 7.1 Alcohol

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ethanol in alcoholic beverages | ● | ◐ | — | — | ● | ● | ◐ | ◐ | — | ● | graph edge, literature |
| Urethane | ● | ● | — | — | ◐ | — | — | — | — | ◐ | literature |

*Within-group heterogeneity:* Ethanol is metabolically activated (ADH → acetaldehyde, Tier 2) while urethane is a direct alkylating agent (Tier 1). Same reporting class, different mechanistic tiers. Ethanol spans 6 KCs; urethane spans 3.

### 7.2 Aldehyde

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Formaldehyde | ● | ● | ◐ | — | ◐ | ● | — | — | — | ● | graph edge, IARC KC |
| Acrolein | ● | ● | — | — | ● | ● | — | — | — | — | literature |
| Crotonaldehyde | ● | ● | — | — | ◐ | — | — | — | — | — | literature |
| Furfural | ◐ | ◐ | — | — | ● | — | — | — | — | — | literature |
| MDA | ● | ● | — | — | ● | — | — | — | — | — | graph edge, IARC Vol. 71 |

*Within-group heterogeneity:* All 5 members are direct-acting electrophiles (Tier 1), but KC profiles vary. Formaldehyde is a chronic inflammation agent (KC6); MDA is an endogenous genotoxin now IARC-evaluated as Group 3 (IARC Vol. 71, 1999).

### 7.3 Alkylating

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EthyleneOxide | ● | ● | — | — | — | — | — | — | — | — | graph edge |
| Acrylamide | ● | ● | — | — | ◐ | ● | ◐ | — | — | ◐ | literature |
| Cyclophosphamide | ● | ● | — | — | ● | — | ● | — | — | — | graph edge, literature |
| Chlorambucil | ● | ● | — | — | — | — | — | — | — | — | literature |
| Sulfur_mustard | ● | ● | — | — | ◐ | ● | — | — | — | — | literature |
| Busulfan | ● | ● | — | — | — | — | — | — | — | ◐ | literature |
| MNU | ● | ● | — | — | — | — | — | — | — | — | literature |
| Temozolomide | ● | ● | — | — | — | — | — | — | — | — | literature |

*Within-group heterogeneity:* Most of the 8 members are Tier 1 direct-acting, but cyclophosphamide is a Tier 2 prodrug requiring CYP3A4 activation and is also immunosuppressive (KC7). Acrylamide requires metabolic activation to glycidamide (now a Metabolite node). Sulfur mustard causes chronic lung inflammation (KC6). Temozolomide is a chemotherapy drug not separately IARC-evaluated. The "alkylating" class name suggests a single mechanism, but KC profiles vary substantially.

### 7.4 Androgen

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AndrogenicAnabolicSteroids | — | — | — | — | ◐ | — | — | ● | — | ● | graph edge |

*Exposure-anchor pattern — AndrogenicAnabolicSteroids (Group 2A) is the IARC-listed exogenous exposure. Endogenous androgens (testosterone, DHT) are now Metabolite nodes without a reporting class.*

### 7.5 Aromatic_Amine

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4ABP | ● | ● | — | — | ◐ | ● | — | — | — | — | graph edge |
| Benzidine | ● | ● | — | — | — | — | — | — | — | — | literature |

*Within-group heterogeneity:* Both require CYP1A2 → NAT bioactivation (Tier 2), but 4ABP has stronger evidence for chronic inflammation (KC6, bladder carcinogenesis). Narrow KC profiles for a chemical class.

### 7.6 Benzene

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Benzene | ● | ● | ◐ | — | ● | ● | ● | — | — | ◐ | graph edge, literature |

*Single-member class.* Benzene has one of the broadest KC profiles (7 of 10 KCs): metabolically activated (KC1), genotoxic (KC2), oxidative stress (KC5), chronic inflammation (KC6), immunosuppression (KC7, bone marrow toxicity). This breadth is unusual — most carcinogens span 2–4 KCs.

### 7.7 Chlorinated_Solvent

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| TCE | ● | ◐ | — | — | ● | ● | — | ◐ | — | ◐ | graph edge, literature |
| PCE | ● | ◐ | — | — | ● | — | — | — | — | — | literature |

*Within-group heterogeneity:* Both require CYP2E1 activation (Tier 2), but TCE has broader KC coverage including receptor-mediated effects (KC8, PPARα in kidney) and chronic inflammation (KC6).

### 7.8 Dioxin

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| TCDD | — | — | — | ◐ | ● | — | ● | ● | — | ● | graph edge, IARC KC |
| PeCDF_23478 | — | — | — | — | ◐ | — | ◐ | ● | — | ◐ | graph edge |

*Within-group heterogeneity:* Both are purely receptor-mediated (Tier 3, AHR). TCDD has stronger evidence for immunosuppression (KC7), epigenetic alterations (KC4), and cell proliferation effects (KC10). No genotoxicity (KC2) — dioxins are the archetypal non-genotoxic carcinogens.

### 7.9 Estrogen

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EstrogenProgestogenTherapy | — | — | — | — | ◐ | ● | — | ● | — | ● | graph edge, literature |

*Exposure-anchor pattern — EstrogenProgestogenTherapy (Group 1) is the IARC-listed exogenous exposure. Endogenous estrogens (E2, 4-OHE2) are now Metabolite nodes without a reporting class.*

### 7.10 HCA

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PhIP | ● | ● | — | — | ◐ | — | — | — | — | — | graph edge |
| MeIQx | ● | ● | — | — | — | — | — | — | — | — | literature |

*Within-group heterogeneity:* Both require CYP1A2 → NAT bioactivation (Tier 2). PhIP has some oxidative stress evidence (KC5). Narrow, homogeneous KC profiles for a food-borne mutagen class.

### 7.11 Heavy_Metal

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ArsenicInorganic | — | ◐ | ● | ● | ● | ● | — | — | — | ● | graph edge, IARC KC |
| Cd | — | ◐ | ● | ● | ● | — | — | ◐ | — | ● | graph edge, literature |
| CrVI | ◐ | ● | ● | — | ● | — | — | — | — | — | graph edge, literature |
| NickelCompounds | — | ◐ | ● | ● | ● | — | — | — | — | — | literature |
| NickelMetallic | — | — | — | — | ◐ | — | — | — | — | — | literature |
| Be | — | — | — | — | ◐ | ● | ◐ | — | — | — | literature |
| LeadInorganicCompounds | — | ◐ | ◐ | ◐ | ● | — | — | — | — | — | literature |
| LeadMetallic | — | — | — | — | ◐ | — | — | — | — | — | literature |
| LeadOrganicCompounds | — | — | — | — | ◐ | — | — | — | — | — | literature |
| MethylmercuryCompounds | — | — | — | — | ● | ● | ◐ | — | — | — | literature |
| MercuryInorganicCompounds | — | — | — | — | ● | ● | ◐ | — | — | — | literature |
| CobaltMetal | — | ◐ | — | — | ● | — | — | ◐ | — | ◐ | graph edge, literature |
| CobaltOxide | — | ◐ | — | — | ● | — | — | — | — | — | literature |
| AntimonyTrivalent | — | — | — | — | ● | ● | — | — | — | — | literature |
| AntimonyPentavalent | — | — | — | — | ◐ | — | — | — | — | — | literature |

*Within-group heterogeneity:* Most heterogeneous class in the graph. 15 members spanning IARC Groups 1–3. ArsenicInorganic has 6 KCs including epigenetic alterations (KC4) and DNA repair inhibition (KC3). Beryllium has 3 KCs, dominated by chronic inflammation (KC6, berylliosis). CrVI is the only metal with strong KC1 (electrophilic after intracellular Cr(VI)→Cr(III) reduction). CobaltMetal has receptor-mediated effects (KC8, HIF-1α stabilization). The expanded species nodes reveal that within-metal heterogeneity is substantial: nickel compounds (Group 1) have broader KC profiles than metallic nickel (Group 2B); lead inorganic compounds (Group 2A) show DNA repair and epigenetic effects absent from metallic lead (Group 2B) and organic lead (Group 3).

### 7.12 Mycotoxin

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AFB1 | ● | ● | — | — | ● | — | — | — | — | — | graph edge |

*Single-member class.* Requires CYP3A4 bioactivation (Tier 2). Classic genotoxic hepatocarcinogen with narrow KC profile (KC1, KC2, KC5).

### 7.13 Nitrosamine

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| NNK | ● | ● | — | — | ◐ | ● | — | ◐ | — | — | graph edge, literature |
| NDMA | ● | ● | — | — | — | — | — | — | — | — | literature |
| NDEA | ● | ● | — | — | — | — | — | — | — | — | literature |

*Within-group heterogeneity:* All require CYP bioactivation (Tier 2). NNK is the broadest — it also activates β-adrenergic receptors in lung tissue (KC8) and causes chronic inflammation (KC6). NDMA and NDEA are narrower, primarily genotoxic.

### 7.14 Organochlorine

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HCB | — | — | — | — | ● | ● | ◐ | ◐ | — | ◐ | literature |
| Lindane | — | — | — | — | ● | ● | — | ◐ | — | — | literature |
| DDT | — | — | — | ◐ | ● | — | — | ● | — | ● | graph edge, literature |
| DDE | — | — | — | — | ◐ | — | — | ● | — | ◐ | graph edge |
| PCP | — | ◐ | — | — | ● | ● | — | — | — | — | graph edge, literature |
| Chlordane | — | — | — | — | ● | ● | ◐ | ◐ | — | — | literature |
| Heptachlor | — | — | — | — | ● | — | — | ◐ | — | — | literature |
| Toxaphene | — | — | — | — | ● | ● | — | ◐ | — | — | literature |

*Within-group heterogeneity:* **Second most heterogeneous class.** DDT is receptor-mediated (KC8, ER; KC10, endocrine disruption) while PCP is direct-acting (KC2, KC5, mitochondrial uncoupling). HCB and lindane cause chronic inflammation (KC6). Most members have weak receptor-mediated effects (KC8) alongside oxidative stress (KC5). The class spans Tiers 1, 3, and 4.

### 7.15 PAH

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| BaP | ● | ● | — | — | ● | ● | — | ● | — | ◐ | graph edge, literature |
| DMBA \* | ● | ● | — | — | ● | — | — | ● | — | — | literature |

*\* Experimental carcinogen — not a human IARC agent.*

*Within-group heterogeneity:* Both require CYP1A1 bioactivation (Tier 2) and activate AHR (Tier 3). BaP has broader evidence for chronic inflammation (KC6) and cell proliferation (KC10). The PAH class is actually a Tier 4 mixed class (metabolic activation + receptor-mediated), not purely Tier 2.

### 7.16 PCB

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PCB_126 | — | — | — | — | ● | — | ● | ● | — | ● | graph edge, IARC KC |
| PCB_169 | — | — | — | — | ◐ | — | ◐ | ● | — | ◐ | graph edge |
| PCB_77 | — | — | — | — | ◐ | — | — | ● | — | — | graph edge |
| PCB_118 | ◐ | ◐ | — | — | ● | — | — | ◐ | — | — | graph edge |
| PCB_153 | ◐ | ◐ | — | — | ● | — | — | — | — | — | graph edge |
| PCB_138 | ◐ | ◐ | — | — | ● | — | — | — | — | — | graph edge |

*Within-group heterogeneity:* **Cleanest illustration of within-class KC divergence.** Dioxin-like PCBs (126, 169, 77) are purely receptor-mediated (Tier 3, AHR). Non-dioxin-like PCBs (118, 153, 138) require metabolic activation to quinone metabolites (Tier 2, KC1 + KC2). Same reporting class, different mechanistic tiers — exactly the distinction KCC flattens.

### 7.17 PFAS

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PFOA | — | — | — | ● | ● | — | ● | ● | — | ● | council, graph edge |
| PFOS | — | — | — | ● | ● | — | ● | ● | — | ● | council, graph edge |

*Within-group heterogeneity:* Per the [ExposoGraph PFAS Council Synthesis](council-reports/Council_PFAS_Claim_Synthesis.md), IARC Vol. 135 cited five KCs for PFOA (KC4, KC5, KC7, KC8, KC10), with the primary in-human basis being KC4 (epigenetic) and KC7 (immunosuppression) — not KC2 (genotoxicity) or KC8 (receptor-mediated effects). Group 1 was reached **in the absence of strong evidence for KC2 and KC8 in exposed humans**. This is a case where the KCC evidence hierarchy matters: the same KC8 is marked ● (strong in experimental systems) but was explicitly insufficient for the human classification.

### 7.18 Plant_Alkaloid

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| AristolochicAcid | ● | ● | — | — | ◐ | — | — | — | — | — | graph edge |

*Single-member class.* Requires CYP1A1/2 bioactivation (Tier 2). Characteristic A→T transversion mutational signature. Narrow KC profile (KC1, KC2, KC5).

### 7.19 UV_Radiation

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| UVRadiation | — | ● | ◐ | — | ● | ● | — | — | — | ● | council, graph edge |

*Within-group heterogeneity:* Per the [ExposoGraph Radiation Council Synthesis](council-reports/Council_Radiation_Claim_Synthesis.md), all three are physical agents (Tier 1) that directly damage DNA (KC2) and generate oxidative stress (KC5, via radiolysis/photodimerization). **KC1 does not apply** — radiation is not electrophilic and is not metabolically activated; it is a physical energy transfer. Ionizing radiation has the broadest profile including immunosuppression (KC7, radiation-induced immune suppression) and genomic instability (KC3). UV causes chronic inflammation (KC6, skin inflammation) and cell proliferation (KC10, UV-induced hyperplasia).

### 7.20 Ionizing_Radiation

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| IonizingRadiation | — | ● | ◐ | — | ● | ● | ◐ | — | — | ● | council, graph edge |
| Radon | — | ● | — | — | ● | ● | — | — | — | — | council, graph edge |

*Within-group heterogeneity:* Per the [ExposoGraph Radiation Council Synthesis](council-reports/Council_Radiation_Claim_Synthesis.md), all three are physical agents (Tier 1) that directly damage DNA (KC2) and generate oxidative stress (KC5, via radiolysis/photodimerization). **KC1 does not apply** — radiation is not electrophilic and is not metabolically activated; it is a physical energy transfer. Ionizing radiation has the broadest profile including immunosuppression (KC7, radiation-induced immune suppression) and genomic instability (KC3). UV causes chronic inflammation (KC6, skin inflammation) and cell proliferation (KC10, UV-induced hyperplasia).

### 7.21 Vinyl Chloride

| Carcinogen | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | Basis |
|---|---|---|---|---|---|---|---|---|---|---|---|
| VinylChloride | ● | ● | — | — | ● | — | — | — | — | ◐ | graph edge, literature |

*Single-member class.* Requires CYP2E1 bioactivation (Tier 2). Characteristic angiosarcoma of the liver. Narrow KC profile (KC1, KC2, KC5, KC10).

---

## 8. KCC Rollup: Reporting Class × Key Characteristics

This rollup shows, for each reporting class, which KCs are represented by at least one member (● or ◐), and the degree of within-class heterogeneity.

| # | Group | KC1 | KC2 | KC3 | KC4 | KC5 | KC6 | KC7 | KC8 | KC9 | KC10 | KCs spanned | Heterogeneity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Alcohol | ● | ◐ | — | — | ● | ● | ◐ | ◐ | — | ● | 6/10 | High — ethanol spans 6 KCs; urethane only 3 |
| 2 | Aldehyde | ● | ● | ◐ | — | ● | ● | — | — | — | ● | 6/10 | Moderate — 5 members; formaldehyde broadest, furfural narrowest |
| 3 | Alkylating | ● | ● | — | — | ◐ | ◐ | ● | — | — | ◐ | 5/10 | High — 8 members; cyclophosphamide unique KC7; most others KC1+KC2 only |
| 4 | Androgen | — | — | — | — | ◐ | — | — | ● | — | ● | 2/10 | Low — single member (AndrogenicAnabolicSteroids); narrowest profile, purely receptor-mediated |
| 5 | Aromatic_Amine | ● | ● | — | — | ◐ | ● | — | — | — | — | 3/10 | Low — homogeneous Tier 2 class |
| 6 | Benzene | ● | ● | ◐ | — | ● | ● | ● | — | — | ◐ | 7/10 | N/A — single member |
| 7 | Chlorinated_Solvent | ● | ◐ | — | — | ● | ● | — | ◐ | — | ◐ | 5/10 | Moderate — TCE broader than PCE |
| 8 | Dioxin | — | — | — | ◐ | ● | — | ● | ● | — | ● | 5/10 | Moderate — TCDD broader than PeCDF |
| 9 | Estrogen | — | — | — | — | ◐ | ● | — | ● | — | ● | 4/10 | N/A — single member (EstrogenProgestogenTherapy); exogenous exposure anchor |
| 10 | HCA | ● | ● | — | — | ◐ | — | — | — | — | — | 3/10 | Low — homogeneous Tier 2 class |
| 11 | Heavy_Metal | ◐ | ● | ● | ● | ● | ● | ◐ | ◐ | — | ● | 8/10 | **Highest** — 15 members; no two metals share the same KC profile |
| 12 | Mycotoxin | ● | ● | — | — | ● | — | — | — | — | — | 3/10 | N/A — single member |
| 13 | Nitrosamine | ● | ● | — | — | ◐ | ● | — | ◐ | — | — | 4/10 | Moderate — NNK broader (KC6, KC8) |
| 14 | Organochlorine | ◐ | ◐ | — | ◐ | ● | ● | ◐ | ● | — | ◐ | 7/10 | **Second highest** — DDT (KC8+KC10) vs. PCP (KC2+KC5) |
| 15 | PAH | ● | ● | — | — | ● | ● | — | ● | — | ◐ | 5/10 | Moderate — actually Tier 4 mixed, not pure Tier 2 |
| 16 | PCB | ◐ | ◐ | — | — | ● | — | ● | ● | — | ● | 5/10 | **Cleanest divergence** — dioxin-like (Tier 3) vs. non-dioxin (Tier 2) |
| 17 | PFAS | — | — | — | ● | ● | — | ● | ● | — | ● | 5/10 | Low — homogeneous receptor-mediated class; Group 1 despite no KC2 |
| 18 | Plant_Alkaloid | ● | ● | — | — | ◐ | — | — | — | — | — | 3/10 | N/A — single member |
| 19 | UV_Radiation | — | ● | ◐ | — | ● | ● | — | — | — | ● | 5/10 | N/A — single member |
| 20 | Ionizing_Radiation | — | ● | ◐ | — | ● | ● | ◐ | — | — | ● | 5/10 | Moderate — ionizing broader than radon |
| 21 | Vinyl Chloride | ● | ● | — | — | ● | — | — | — | — | ◐ | 4/10 | N/A — single member |

### KC ubiquity across classes

| KC | Characteristic | Classes with ● or ◐ | Ubiquity |
|----|----------------|---------------------|----------|
| KC1 | Electrophilic/metabolically activated | 14/21 | High — but absent from all receptor-mediated and physical classes |
| KC2 | Genotoxic | 15/21 | **Highest** — but absent from PFAS, dioxin, androgen, radiation (physical) |
| KC3 | Alters DNA repair/genomic instability | 5/21 | Moderate — primarily metals |
| KC4 | Epigenetic alterations | 4/21 | Low — arsenic, cadmium, nickel, PFAS |
| KC5 | Oxidative stress | 18/21 | **Near-universal** — nearly non-discriminating |
| KC6 | Chronic inflammation | 12/21 | High — but mechanism varies (tissue injury vs. immune dysregulation) |
| KC7 | Immunosuppressive | 7/21 | Moderate — selective (dioxin, PFAS, benzene, cyclophosphamide, metals) |
| KC8 | Receptor-mediated effects | 9/21 | Moderate — defines Tier 3 classes but also appears in mixed classes |
| KC9 | Immortalization | 0/21 | **Absent** — rarely useful for chemical/physical agents |
| KC10 | Cell proliferation/death/nutrient supply | 11/21 | High — overlaps heavily with KC8 (receptor-mediated agents) |

---

## 9. Why KCC Does Not Map Cleanly to the 21 Reporting Classes

### Problem 1: Near-universal KCs are non-discriminating

KC2 (genotoxic) appears in 15 of 21 classes and KC5 (oxidative stress) in 18 of 21. These KCs cannot distinguish between a direct-acting aldehyde and a metabolically activated PAH — both are "genotoxic" and both "induce oxidative stress," but they require entirely different bodily processing chains and activate different CYP enzymes.

### Problem 2: Within-class KC heterogeneity exceeds between-class differences

Heavy_Metal spans 8 of 10 KCs across its 15 members. ArsenicInorganic alone covers 6 KCs (KC2, KC3, KC4, KC5, KC6, KC10) while beryllium covers only 3 (KC5, KC6, KC7). Two metals in the same reporting class have less KC overlap than arsenic does with benzene (both have KC2, KC5, KC6, KC7). The reporting class captures shared chemical/processing properties; KCC captures downstream cellular effects that vary by metal.

### Problem 3: Same KC, different mechanistic tier

PCB_126 (dioxin-like) and PCB_153 (non-dioxin) are in the same reporting class. PCB_126 has KC8 (AHR activation) without KC1/KC2. PCB_153 has KC1/KC2 (metabolic activation to quinone adducts) without strong KC8. They share KC5 (oxidative stress). KCC records "oxidative stress" for both but cannot distinguish the AHR-mediated pathway from the CYP450-activation pathway. ExposoGraph's mechanistic tier (Tier 3 vs. Tier 2) captures this distinction.

### Problem 4: KCC evidence hierarchy is lost in a flat matrix

The PFAS council synthesis revealed that PFOA's Group 1 classification was based on KC4 (epigenetic) and KC7 (immunosuppression) in exposed humans, despite KC8 (receptor-mediated) being marked ● based on experimental evidence. The KCC matrix records both as ● but cannot represent that one was sufficient for human hazard identification and the other was explicitly insufficient. ExposoGraph's council synthesis and provenance fields capture this evidence hierarchy.

### Problem 5: Processing chain invisible

KCC says BaP "is genotoxic" (KC2) and "is metabolically activated" (KC1). It does not specify that the activation requires CYP1A1, proceeds through three sequential metabolites (BaP epoxide → diol → BPDE), and produces a specific DNA adduct (BPDE-dG). ExposoGraph's graph edges encode this full chain, making it possible to identify shared enzymatic dependencies (e.g., CYP2E1 processes benzene, vinyl chloride, and TCE — a connection invisible in both IARC's agent type and KCC).

---

## 10. Conclusion

The three frameworks serve complementary purposes:

| Framework | Question answered | Granularity |
|-----------|-------------------|-------------|
| IARC Agent Formalism | What kind of thing is it? Where does it come from? | Exposure source and agent identity |
| KCC (Key Characteristics) | What does it do to the cell? | 10 binary-ish cellular effects |
| ExposoGraph Mechanistic Tier | What does the body do to it first? | Enzyme → metabolite → adduct chain |

ExposoGraph's 21 reporting classes are organized by mechanistic processing similarity (shared CYP enzymes, shared adduct types, shared receptor targets). This grouping captures information that neither IARC's source-based classification nor KCC's effect-based characteristics can represent: the bodily processing chain from exposure to DNA damage, including shared enzymatic dependencies that create correlated risk profiles across seemingly unrelated agents.

---

## References

1. IARC Preamble (2019). [Preamble to the IARC Monographs](https://monographs.iarc.who.int/wp-content/uploads/2019/01/Preamble-2019.pdf). International Agency for Research on Cancer.
2. Smith MT, et al. (2016). [Key Characteristics of Carcinogens](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4898322/). *Environ Health Perspect* 124(6):713–722.
3. IARC Monographs Vol. 100D (2012). [Radiation](https://publications.iarc.who.int/Book-And-Report-Series/Iarc-Monographs-On-The-Identification-Of-Carcinogenic-Hazards-To-Humans/Radiation-2012).
4. IARC Monographs Vol. 135 (2025). PFOA and PFOS. International Agency for Research on Cancer.
5. ExposoGraph Model Council. Radiation Claim Synthesis. `council-reports/Council_Radiation_Claim_Synthesis.md`.
6. ExposoGraph Model Council. PFAS Claim Synthesis. `council-reports/Council_PFAS_Claim_Synthesis.md`.
7. IARC Monographs Vol. 53 (1991). [Occupational Exposures in Insecticide Application, and Some Pesticides](https://publications.iarc.who.int/Book-And-Report-Series/Iarc-Monographs-On-The-Identification-Of-Carcinogenic-Hazards-To-Humans/Occupational-Exposures-In-Insecticide-Application-And-Some-Pesticides-1991). International Agency for Research on Cancer. (DDT and associated compounds, including DDE.)
8. IARC Monographs Vol. 58 (1993). [Beryllium, Cadmium, Mercury, and Exposures in the Glass Manufacturing Industry](https://publications.iarc.who.int/Book-And-Report-Series/Iarc-Monographs-On-The-Identification-Of-Carcinogenic-Hazards-To-Humans/Beryllium-Cadmium-Mercury-And-Exposures-In-The-Glass-Manufacturing-Industry-1993). International Agency for Research on Cancer.
9. IARC Monographs Vol. 71 (1999). [Re-evaluation of Some Organic Chemicals, Hydrazine and Hydrogen Peroxide](https://publications.iarc.who.int/Book-And-Report-Series/Iarc-Monographs-On-The-Identification-Of-Carcinogenic-Hazards-To-Humans/Re-evaluation-Of-Some-Organic-Chemicals-Hydrazine-And-Hydrogen-Peroxide-Part-1-Part-2-Part-3--1999). International Agency for Research on Cancer. (Includes 4,4′-methylenedianiline, Group 3.)
10. IARC Monographs Vol. 87 (2006). [Inorganic and Organic Lead Compounds](https://publications.iarc.who.int/Book-And-Report-Series/Iarc-Monographs-On-The-Identification-Of-Carcinogenic-Hazards-To-Humans/Inorganic-And-Organic-Lead-Compounds-2006). International Agency for Research on Cancer.
11. IARC Monographs Vol. 100C (2012). [Arsenic, Metals, Fibres, and Dusts](https://publications.iarc.who.int/Book-And-Report-Series/Iarc-Monographs-On-The-Identification-Of-Carcinogenic-Hazards-To-Humans/Arsenic-Metals-Fibres-And-Dusts-2012). International Agency for Research on Cancer. (Arsenic, cadmium, chromium(VI), nickel, beryllium.)
12. IARC Monographs Vol. 131 (2023). [Cobalt, Antimony Compounds, and Weapons-grade Tungsten Alloy](https://monographs.iarc.who.int/iarc-monographs-volume-131/). International Agency for Research on Cancer.
