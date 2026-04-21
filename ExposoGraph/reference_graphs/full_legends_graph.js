// ═══════════════════════════════════════════════════════════════════════
// CarcinoGenomic Platform — Layer 2: Carcinogen-Gene Knowledge Graph Data
// Node types: Carcinogen, Enzyme, Metabolite, DNA_Adduct, Pathway
// Edge types: ACTIVATES, DETOXIFIES, TRANSPORTS, FORMS_ADDUCT, REPAIRS
// ═══════════════════════════════════════════════════════════════════════

const GRAPH_DATA = {
  nodes: [
    // ── CARCINOGENS (IARC Group 1) ───────────────────────────────────
    { id: "BaP", label: "Benzo[a]pyrene", type: "Carcinogen", group: "PAH", iarc: "Group 1", detail: "Prototypical PAH; tobacco smoke, grilled food, coal tar" },
    { id: "DMBA", label: "DMBA", type: "Carcinogen", group: "PAH", iarc: "Group 2A", detail: "7,12-Dimethylbenz[a]anthracene; experimental PAH" },
    { id: "PhIP", label: "PhIP", type: "Carcinogen", group: "HCA", iarc: "Group 2B", detail: "2-Amino-1-methyl-6-phenylimidazo[4,5-b]pyridine; cooked meat" },
    { id: "MeIQx", label: "MeIQx", type: "Carcinogen", group: "HCA", iarc: "Group 2B", detail: "2-Amino-3,8-dimethylimidazo[4,5-f]quinoxaline; cooked meat" },
    { id: "4ABP", label: "4-Aminobiphenyl", type: "Carcinogen", group: "Aromatic_Amine", iarc: "Group 1", detail: "Aromatic amine; tobacco smoke, dye industry" },
    { id: "BNZ", label: "Benzidine", type: "Carcinogen", group: "Aromatic_Amine", iarc: "Group 1", detail: "Aromatic amine; dye manufacturing" },
    { id: "NNK", label: "NNK", type: "Carcinogen", group: "Nitrosamine", iarc: "Group 1", detail: "4-(Methylnitrosamino)-1-(3-pyridyl)-1-butanone; tobacco" },
    { id: "NDMA", label: "NDMA", type: "Carcinogen", group: "Nitrosamine", iarc: "Group 2A", detail: "N-Nitrosodimethylamine; processed meat, contaminated water" },
    { id: "AFB1", label: "Aflatoxin B1", type: "Carcinogen", group: "Mycotoxin", iarc: "Group 1", detail: "Aspergillus flavus mycotoxin; contaminated grains/nuts" },
    { id: "E2", label: "17β-Estradiol", type: "Carcinogen", group: "Estrogen", iarc: "Group 1", detail: "Endogenous estrogen; breast/endometrial carcinogenesis" },
    { id: "TESTO", label: "Testosterone", type: "Carcinogen", group: "Androgen", iarc: "Group 2A", detail: "Primary androgen; prostate carcinogenesis via AR signaling and aromatization to estradiol" },
    { id: "DHT", label: "5α-DHT", type: "Carcinogen", group: "Androgen", iarc: "Group 2A", detail: "5α-Dihydrotestosterone; most potent androgen; non-aromatizable; AR-mediated proliferation" },
    { id: "BENZ", label: "Benzene", type: "Carcinogen", group: "Solvent", iarc: "Group 1", detail: "Industrial solvent; petroleum products, tobacco smoke" },
    { id: "VCM", label: "Vinyl Chloride", type: "Carcinogen", group: "Solvent", iarc: "Group 1", detail: "PVC production; liver angiosarcoma" },
    { id: "EO", label: "Ethylene Oxide", type: "Carcinogen", group: "Alkylating", iarc: "Group 1", detail: "Sterilization agent; direct alkylator" },

    // ── PHASE I ENZYMES (Activation) ──────────────────────────────────
    { id: "CYP1A1", label: "CYP1A1", type: "Enzyme", phase: "I", role: "Activation", detail: "PAH → diol-epoxides; extrahepatic expression; AhR-inducible" },
    { id: "CYP1B1", label: "CYP1B1", type: "Enzyme", phase: "I", role: "Activation", detail: "PAH activation; estradiol 4-hydroxylation; mammary/ovary" },
    { id: "CYP1A2", label: "CYP1A2", type: "Enzyme", phase: "I", role: "Activation", detail: "HCA/aromatic amine N-oxidation; hepatic; AhR-inducible" },
    { id: "CYP2A6", label: "CYP2A6", type: "Enzyme", phase: "I", role: "Activation", detail: "NNK α-hydroxylation (hepatic); nicotine metabolism" },
    { id: "CYP2A13", label: "CYP2A13", type: "Enzyme", phase: "I", role: "Activation", detail: "NNK α-hydroxylation (pulmonary); respiratory-specific" },
    { id: "CYP2B6", label: "CYP2B6", type: "Enzyme", phase: "I", role: "Activation", detail: "Auxiliary xenobiotic CYP represented for expanded Phase I competition coverage, including NNK-linked activation context" },
    { id: "CYP2C9", label: "CYP2C9", type: "Enzyme", phase: "I", role: "Activation", detail: "Broad hepatic xenobiotic CYP represented for expanded Phase I competition coverage and secondary aromatic-amine/NNK context" },
    { id: "CYP2C19", label: "CYP2C19", type: "Enzyme", phase: "I", role: "Activation", detail: "Hepatic CYP represented for expanded Phase I competition coverage with NNK-linked activation context" },
    { id: "CYP2D6", label: "CYP2D6", type: "Enzyme", phase: "I", role: "Activation", detail: "Phase I CYP represented for expanded competition coverage, including NNK pyridyloxobutylation context" },
    { id: "CYP2E1", label: "CYP2E1", type: "Enzyme", phase: "I", role: "Activation", detail: "Benzene/NDMA/vinyl chloride activation; ethanol-inducible" },
    { id: "CYP2F1", label: "CYP2F1", type: "Enzyme", phase: "I", role: "Activation", detail: "Pulmonary xenobiotic CYP represented for lung-specific benzene bioactivation and inhalation-focused competition context" },
    { id: "CYP3A4", label: "CYP3A4", type: "Enzyme", phase: "I", role: "Activation", detail: "AFB1 8,9-epoxidation; most abundant hepatic CYP" },
    { id: "EPHX1", label: "EPHX1", type: "Enzyme", phase: "I", role: "Mixed", detail: "Epoxide hydrolase; converts epoxides to diols (can activate or detoxify)" },
    { id: "CYP17A1", label: "CYP17A1", type: "Enzyme", phase: "I", role: "Activation", detail: "Steroid 17α-hydroxylase/17,20-lyase; androgen biosynthesis from C21 precursors" },
    { id: "SRD5A2", label: "SRD5A2", type: "Enzyme", phase: "I", role: "Activation", detail: "Steroid 5α-reductase type 2; converts testosterone to DHT in prostate" },
    { id: "SRD5A1", label: "SRD5A1", type: "Enzyme", phase: "I", role: "Activation", detail: "Steroid 5α-reductase type 1; peripheral testosterone-to-DHT conversion" },
    { id: "CYP19A1", label: "CYP19A1", type: "Enzyme", phase: "I", role: "Activation", detail: "Aromatase; converts testosterone to 17β-estradiol; genotoxic bridging" },
    { id: "CYP3A5", label: "CYP3A5", type: "Enzyme", phase: "I", role: "Mixed", detail: "6β-hydroxylation of testosterone; AR nuclear translocation in prostate" },
    { id: "AKR1C3", label: "AKR1C3", type: "Enzyme", phase: "I", role: "Activation", detail: "Aldo-keto reductase 1C3; converts androstenedione to testosterone; CRPC driver" },

    // ── PHASE II ENZYMES (Detoxification) ─────────────────────────────
    { id: "GSTM1", label: "GSTM1", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Glutathione conjugation of PAH diol-epoxides, AFB1-epoxide" },
    { id: "GSTT1", label: "GSTT1", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Glutathione conjugation; small reactive epoxides" },
    { id: "GSTP1", label: "GSTP1", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Glutathione conjugation; PAH metabolites; lung expression" },
    { id: "NAT2", label: "NAT2", type: "Enzyme", phase: "II", role: "Mixed", detail: "N-acetylation (detoxifies amines) / O-acetylation (activates N-OH-amines)" },
    { id: "NAT1", label: "NAT1", type: "Enzyme", phase: "II", role: "Activation", detail: "O-acetylation in bladder epithelium; activates N-OH aromatic amines" },
    { id: "SULT1A1", label: "SULT1A1", type: "Enzyme", phase: "II", role: "Mixed", detail: "Sulfotransferase; can activate or detoxify depending on substrate" },
    { id: "UGT1A1", label: "UGT1A1", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Glucuronidation of HCA metabolites; bilirubin metabolism" },
    { id: "UGT2B7", label: "UGT2B7", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Glucuronidation of carcinogen metabolites" },
    { id: "NQO1", label: "NQO1", type: "Enzyme", phase: "II", role: "Detoxification", detail: "NAD(P)H:quinone reductase; benzene quinone detoxification" },
    { id: "COMT", label: "COMT", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Catechol-O-methyltransferase; methylates catechol estrogens" },
    { id: "UGT2B17", label: "UGT2B17", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Primary glucuronidation of testosterone/DHT; high deletion polymorphism" },
    { id: "UGT2B15", label: "UGT2B15", type: "Enzyme", phase: "II", role: "Detoxification", detail: "Glucuronidation of DHT/androsterone; compensates for UGT2B17 deletion" },
    { id: "HSD3B2", label: "HSD3B2", type: "Enzyme", phase: "I", role: "Activation", detail: "3β-Hydroxysteroid dehydrogenase type 2; DHEA→androstenedione conversion" },
    { id: "AKR1C2", label: "AKR1C2", type: "Enzyme", phase: "II", role: "Detoxification", detail: "3α-Hydroxysteroid dehydrogenase; DHT→3α-androstanediol inactivation" },

    // ── PHASE III TRANSPORTERS ────────────────────────────────────────
    { id: "ABCB1", label: "ABCB1 (P-gp)", type: "Enzyme", phase: "III", role: "Transport", detail: "P-glycoprotein; efflux of xenobiotic conjugates" },
    { id: "ABCC2", label: "ABCC2 (MRP2)", type: "Enzyme", phase: "III", role: "Transport", detail: "MRP2; biliary/renal efflux of glutathione conjugates" },
    { id: "ABCG2", label: "ABCG2 (BCRP)", type: "Enzyme", phase: "III", role: "Transport", detail: "BCRP; efflux of sulfate/glucuronide conjugates" },

    // ── DNA REPAIR ENZYMES ────────────────────────────────────────────
    { id: "XRCC1", label: "XRCC1", type: "Enzyme", phase: "Repair", role: "Repair", detail: "Base excision repair scaffold; oxidative/alkylation damage" },
    { id: "XPC", label: "XPC", type: "Enzyme", phase: "Repair", role: "Repair", detail: "Nucleotide excision repair; bulky adduct recognition" },
    { id: "ERCC2", label: "ERCC2/XPD", type: "Enzyme", phase: "Repair", role: "Repair", detail: "NER helicase; opens DNA at bulky adduct damage sites" },
    { id: "OGG1", label: "OGG1", type: "Enzyme", phase: "Repair", role: "Repair", detail: "8-oxoguanine glycosylase; repairs oxidative DNA damage" },
    { id: "MGMT", label: "MGMT", type: "Enzyme", phase: "Repair", role: "Repair", detail: "O6-methylguanine-DNA methyltransferase; removes alkyl adducts" },

    // ── METABOLITES (Reactive Intermediates) ──────────────────────────
    { id: "BPDE", label: "BPDE", type: "Metabolite", reactivity: "High", detail: "Benzo[a]pyrene-7,8-diol-9,10-epoxide; ultimate PAH carcinogen" },
    { id: "BaP_78diol", label: "BaP-7,8-diol", type: "Metabolite", reactivity: "Intermediate", detail: "Trans-7,8-dihydrodiol; proximate PAH carcinogen" },
    { id: "BaP_78epox", label: "BaP-7,8-epoxide", type: "Metabolite", reactivity: "High", detail: "Initial epoxide; substrate for EPHX1" },
    { id: "NOHPHIP", label: "N-OH-PhIP", type: "Metabolite", reactivity: "High", detail: "N-hydroxy-PhIP; reactive intermediate for HCA" },
    { id: "NOH4ABP", label: "N-OH-4-ABP", type: "Metabolite", reactivity: "High", detail: "N-hydroxy-4-aminobiphenyl; bladder carcinogen intermediate" },
    { id: "NOHBNZ", label: "N-OH-Benzidine", type: "Metabolite", reactivity: "High", detail: "N-hydroxy-benzidine; reactive intermediate" },
    { id: "NNK_alpha", label: "α-OH-NNK", type: "Metabolite", reactivity: "High", detail: "α-Hydroxylated NNK; methyl/pyridyloxobutyl diazohydroxide" },
    { id: "AFB1_epox", label: "AFB1-8,9-epoxide", type: "Metabolite", reactivity: "High", detail: "Exo-AFB1-8,9-epoxide; ultimate aflatoxin carcinogen" },
    { id: "4OHE2", label: "4-OH-Estradiol", type: "Metabolite", reactivity: "High", detail: "Catechol estrogen; forms depurinating DNA adducts" },
    { id: "2OHE2", label: "2-OH-Estradiol", type: "Metabolite", reactivity: "Low", detail: "Catechol estrogen; less genotoxic pathway" },
    { id: "E2_quinone", label: "E2-3,4-Quinone", type: "Metabolite", reactivity: "High", detail: "Estrogen quinone; forms depurinating DNA adducts" },
    { id: "BenzO", label: "Benzene oxide", type: "Metabolite", reactivity: "High", detail: "Initial benzene metabolite; rearranges to phenol" },
    { id: "HQ", label: "Hydroquinone", type: "Metabolite", reactivity: "High", detail: "Benzene metabolite; oxidized to benzoquinone" },
    { id: "BQ", label: "p-Benzoquinone", type: "Metabolite", reactivity: "High", detail: "Highly reactive quinone; protein/DNA adducts" },
    { id: "CEO", label: "Chloroethylene oxide", type: "Metabolite", reactivity: "High", detail: "Vinyl chloride epoxide; etheno-DNA adducts" },
    { id: "PhIP_NAc", label: "PhIP-N-acetyl", type: "Metabolite", reactivity: "Low", detail: "N-acetylated PhIP; detoxified conjugate" },
    { id: "BPDE_GSH", label: "BPDE-glutathione", type: "Metabolite", reactivity: "Low", detail: "Glutathione conjugate of BPDE; water-soluble" },
    { id: "AFB1_GSH", label: "AFB1-GSH", type: "Metabolite", reactivity: "Low", detail: "Glutathione conjugate of AFB1 epoxide; detoxified" },
    { id: "PhIP_gluc", label: "PhIP-glucuronide", type: "Metabolite", reactivity: "Low", detail: "Glucuronide conjugate; renal excretion" },
    { id: "E2_methyl", label: "Methoxy-E2", type: "Metabolite", reactivity: "Low", detail: "Methylated catechol estrogen; non-genotoxic" },
    { id: "NDMA_alpha", label: "α-OH-NDMA", type: "Metabolite", reactivity: "High", detail: "Methyldiazohydroxide; forms O6-methylguanine" },
    // ── DMBA METABOLISM CHAIN (added) ────────────────────────────────
    { id: "DMBA_epoxide", label: "DMBA-3,4-epoxide", type: "Metabolite", reactivity: "High", detail: "Initial DMBA epoxide; substrate for EPHX1 diol formation" },
    { id: "DMBA_diol", label: "DMBA-3,4-diol", type: "Metabolite", reactivity: "Intermediate", detail: "Trans-3,4-dihydrodiol DMBA; proximate PAH carcinogen" },
    { id: "DMBA_diol_epoxide", label: "DMBA-diol-epoxide", type: "Metabolite", reactivity: "High", detail: "Ultimate DMBA metabolite; bulky DNA adduct precursor" },
    // ── MeIQx METABOLISM CHAIN (added) ───────────────────────────────
    { id: "NOH_MeIQx", label: "N-OH-MeIQx", type: "Metabolite", reactivity: "High", detail: "N-hydroxylated MeIQx activation intermediate" },
    { id: "MeIQx_acetoxy", label: "MeIQx-N-acetoxy ester", type: "Metabolite", reactivity: "High", detail: "Reactive acetoxy ester representing MeIQx activation" },
    // ── ANDROGEN METABOLITES ─────────────────────────────────────────
    { id: "DHEA", label: "DHEA", type: "Metabolite", reactivity: "Low", detail: "Dehydroepiandrosterone; adrenal androgen precursor" },
    { id: "A4", label: "Androstenedione", type: "Metabolite", reactivity: "Low", detail: "Δ4-Androstenedione; immediate testosterone precursor" },
    { id: "6bOHT", label: "6β-OH-Testosterone", type: "Metabolite", reactivity: "Low", detail: "CYP3A4/5-catalyzed hydroxylated testosterone; inactive metabolite" },
    { id: "T_gluc", label: "Testosterone-glucuronide", type: "Metabolite", reactivity: "Low", detail: "UGT2B17-catalyzed glucuronide conjugate; renal excretion" },
    { id: "DHT_gluc", label: "DHT-glucuronide", type: "Metabolite", reactivity: "Low", detail: "UGT2B17/2B15-catalyzed glucuronide; DHT inactivation" },
    { id: "3aAdiol", label: "3α-Androstanediol", type: "Metabolite", reactivity: "Low", detail: "AKR1C2-catalyzed DHT metabolite; weak ERβ agonist" },
    { id: "E2_from_T", label: "E2 (from T)", type: "Metabolite", reactivity: "High", detail: "Estradiol produced by aromatization of testosterone; genotoxic bridge" },

    // ── DNA ADDUCTS ───────────────────────────────────────────────────
    { id: "BPDE_dG", label: "BPDE-N2-dG", type: "DNA_Adduct", detail: "Major PAH-DNA adduct; G→T transversions in p53 codon 249" },
    { id: "PhIP_dG", label: "PhIP-C8-dG", type: "DNA_Adduct", detail: "HCA-DNA adduct at C8 of guanine; frameshift mutations" },
    { id: "ABP_dG", label: "4-ABP-C8-dG", type: "DNA_Adduct", detail: "Aromatic amine-DNA adduct; bladder mutagenesis" },
    { id: "O6MeG", label: "O6-MeG", type: "DNA_Adduct", detail: "O6-methylguanine; G:C→A:T transitions; alkylating damage" },
    { id: "POB_dG", label: "POB-DNA", type: "DNA_Adduct", detail: "Pyridyloxobutyl-DNA adducts from NNK metabolism" },
    { id: "AFB1_N7G", label: "AFB1-N7-dG", type: "DNA_Adduct", detail: "Major aflatoxin adduct; p53 codon 249 AGG→AGT hotspot" },
    { id: "etheno_dA", label: "εdA (etheno-dA)", type: "DNA_Adduct", detail: "1,N6-etheno-deoxyadenosine; vinyl chloride marker" },
    { id: "etheno_dC", label: "εdC (etheno-dC)", type: "DNA_Adduct", detail: "3,N4-etheno-deoxycytidine; vinyl chloride marker" },
    { id: "E2_depurin", label: "E2-depurinating", type: "DNA_Adduct", detail: "4-OHE2 depurinating adducts; apurinic sites → mutations" },
    { id: "8oxodG", label: "8-oxo-dG", type: "DNA_Adduct", detail: "Oxidative DNA damage marker; G→T transversions" },
    { id: "BQ_dG", label: "BQ-DNA", type: "DNA_Adduct", detail: "Benzoquinone-DNA adducts; hematotoxic lesions" },
    // ── DMBA DNA ADDUCT (added) ──────────────────────────────────────
    { id: "DMBA_dA", label: "DMBA-dA", type: "DNA_Adduct", detail: "Representative DMBA-derived deoxyadenosine adduct" },

    // ── PATHWAYS (KEGG References) ────────────────────────────────────
    { id: "KEGG980", label: "Xenobiotic CYP450\n(KEGG:00980)", type: "Pathway", detail: "Metabolism of xenobiotics by cytochrome P450" },
    { id: "KEGG5208", label: "Chemical Carcino.\nROS (KEGG:05208)", type: "Pathway", detail: "Chemical carcinogenesis — reactive oxygen species" },
    { id: "KEGG5204", label: "Chemical Carcino.\nDNA (KEGG:05204)", type: "Pathway", detail: "Chemical carcinogenesis — DNA adducts" },
    { id: "KEGG982", label: "Drug Metabolism\nCYP450 (KEGG:00982)", type: "Pathway", detail: "Drug metabolism — cytochrome P450" },
    { id: "KEGG980G", label: "Glutathione\nMetabolism (KEGG:00480)", type: "Pathway", detail: "Glutathione metabolism pathway" },
    { id: "KEGG_ANDRO", label: "Steroid Hormone\nBiosynthesis (KEGG:00140)", type: "Pathway", detail: "Androgen/estrogen biosynthesis and metabolism pathway" },
  ],

  edges: [
    // ═══ PAH METABOLISM CHAIN ═══════════════════════════════════════
    // BaP → CYP1A1 → BaP-7,8-epoxide → EPHX1 → BaP-7,8-diol → CYP1A1 → BPDE → BPDE-N2-dG
    { source: "CYP1A1", target: "BaP_78epox", type: "ACTIVATES", label: "epoxidation", carcinogen: "BaP" },
    { source: "CYP1B1", target: "BaP_78epox", type: "ACTIVATES", label: "epoxidation", carcinogen: "BaP" },
    { source: "EPHX1", target: "BaP_78diol", type: "ACTIVATES", label: "hydrolysis", carcinogen: "BaP" },
    { source: "CYP1A1", target: "BPDE", type: "ACTIVATES", label: "second epoxidation", carcinogen: "BaP" },
    { source: "CYP1B1", target: "BPDE", type: "ACTIVATES", label: "second epoxidation", carcinogen: "BaP" },
    { source: "BPDE", target: "BPDE_dG", type: "FORMS_ADDUCT", label: "N2-dG adduct" },
    // PAH detoxification
    { source: "GSTM1", target: "BPDE_GSH", type: "DETOXIFIES", label: "GSH conjugation", carcinogen: "BaP" },
    { source: "GSTP1", target: "BPDE_GSH", type: "DETOXIFIES", label: "GSH conjugation", carcinogen: "BaP" },
    { source: "ABCB1", target: "BPDE_GSH", type: "TRANSPORTS", label: "efflux" },
    { source: "ABCC2", target: "BPDE_GSH", type: "TRANSPORTS", label: "biliary efflux" },
    // PAH DNA repair
    { source: "XPC", target: "BPDE_dG", type: "REPAIRS", label: "NER recognition" },
    { source: "ERCC2", target: "BPDE_dG", type: "REPAIRS", label: "NER helicase" },
    // PAH oxidative damage
    { source: "BaP", target: "8oxodG", type: "FORMS_ADDUCT", label: "ROS generation" },
    { source: "OGG1", target: "8oxodG", type: "REPAIRS", label: "8-oxoG excision" },

    // ═══ HCA METABOLISM CHAIN ═══════════════════════════════════════
    // PhIP → CYP1A2 (N-oxidation) → N-OH-PhIP → NAT2/SULT → reactive ester → DNA
    { source: "CYP1A2", target: "NOHPHIP", type: "ACTIVATES", label: "N-oxidation", carcinogen: "PhIP" },
    { source: "NAT2", target: "NOHPHIP", type: "ACTIVATES", label: "O-acetylation", carcinogen: "PhIP" },
    { source: "SULT1A1", target: "NOHPHIP", type: "ACTIVATES", label: "O-sulfonation", carcinogen: "PhIP" },
    { source: "NOHPHIP", target: "PhIP_dG", type: "FORMS_ADDUCT", label: "C8-dG adduct" },
    // HCA detoxification
    { source: "UGT1A1", target: "PhIP_gluc", type: "DETOXIFIES", label: "glucuronidation", carcinogen: "PhIP" },
    { source: "NAT2", target: "PhIP_NAc", type: "DETOXIFIES", label: "N-acetylation", carcinogen: "PhIP" },
    { source: "ABCG2", target: "PhIP_gluc", type: "TRANSPORTS", label: "efflux" },
    // HCA repair
    { source: "XRCC1", target: "PhIP_dG", type: "REPAIRS", label: "BER" },

    // ═══ AROMATIC AMINE METABOLISM ═══════════════════════════════════
    // 4-ABP → CYP1A2 (N-oxidation) → N-OH-4-ABP → NAT1 (O-acetylation in bladder) → DNA
    { source: "CYP1A2", target: "NOH4ABP", type: "ACTIVATES", label: "N-oxidation", carcinogen: "4ABP" },
    { source: "NAT1", target: "NOH4ABP", type: "ACTIVATES", label: "O-acetylation (bladder)", carcinogen: "4ABP" },
    { source: "NOH4ABP", target: "ABP_dG", type: "FORMS_ADDUCT", label: "C8-dG adduct" },
    // Benzidine → CYP1A2 → N-OH-Benzidine → bladder adducts
    { source: "CYP1A2", target: "NOHBNZ", type: "ACTIVATES", label: "N-oxidation", carcinogen: "BNZ" },
    // Aromatic amine detoxification
    { source: "NAT2", target: "4ABP", type: "DETOXIFIES", label: "N-acetylation (deactivation)", carcinogen: "4ABP" },
    { source: "GSTM1", target: "NOH4ABP", type: "DETOXIFIES", label: "GSH conjugation", carcinogen: "4ABP" },
    { source: "XRCC1", target: "ABP_dG", type: "REPAIRS", label: "BER" },

    // ═══ NITROSAMINE METABOLISM ══════════════════════════════════════
    // NNK → CYP2A6/2A13 (α-hydroxylation) → methyldiazohydroxide → O6-MeG
    { source: "CYP2A6", target: "NNK_alpha", type: "ACTIVATES", label: "α-hydroxylation", carcinogen: "NNK" },
    { source: "CYP2A13", target: "NNK_alpha", type: "ACTIVATES", label: "α-hydroxylation (lung)", carcinogen: "NNK" },
    { source: "NNK_alpha", target: "O6MeG", type: "FORMS_ADDUCT", label: "methylation" },
    { source: "NNK_alpha", target: "POB_dG", type: "FORMS_ADDUCT", label: "POB adduction" },
    // NDMA → CYP2E1 → methyldiazohydroxide → O6-MeG
    { source: "CYP2E1", target: "NDMA_alpha", type: "ACTIVATES", label: "α-hydroxylation", carcinogen: "NDMA" },
    { source: "NDMA_alpha", target: "O6MeG", type: "FORMS_ADDUCT", label: "methylation" },
    // Nitrosamine detoxification/repair
    { source: "UGT2B7", target: "NNK_alpha", type: "DETOXIFIES", label: "glucuronidation", carcinogen: "NNK" },
    { source: "MGMT", target: "O6MeG", type: "REPAIRS", label: "methyl transfer" },
    { source: "ABCB1", target: "NNK_alpha", type: "TRANSPORTS", label: "efflux" },

    // ═══ AFLATOXIN B1 METABOLISM ═════════════════════════════════════
    // AFB1 → CYP3A4 → AFB1-8,9-epoxide → DNA adduct
    { source: "CYP3A4", target: "AFB1_epox", type: "ACTIVATES", label: "8,9-epoxidation", carcinogen: "AFB1" },
    { source: "CYP1A2", target: "AFB1_epox", type: "ACTIVATES", label: "3α-hydroxylation", carcinogen: "AFB1" },
    { source: "AFB1_epox", target: "AFB1_N7G", type: "FORMS_ADDUCT", label: "N7-dG adduct" },
    // AFB1 detoxification
    { source: "GSTM1", target: "AFB1_GSH", type: "DETOXIFIES", label: "GSH conjugation", carcinogen: "AFB1" },
    { source: "GSTT1", target: "AFB1_GSH", type: "DETOXIFIES", label: "GSH conjugation", carcinogen: "AFB1" },
    { source: "EPHX1", target: "AFB1_epox", type: "DETOXIFIES", label: "hydrolysis", carcinogen: "AFB1" },
    { source: "ABCC2", target: "AFB1_GSH", type: "TRANSPORTS", label: "biliary efflux" },
    { source: "XPC", target: "AFB1_N7G", type: "REPAIRS", label: "NER" },

    // ═══ ESTROGEN METABOLISM ═════════════════════════════════════════
    // E2 → CYP1B1 (4-hydroxylation) → 4-OH-E2 → quinone → depurinating adducts
    { source: "CYP1B1", target: "4OHE2", type: "ACTIVATES", label: "4-hydroxylation", carcinogen: "E2" },
    { source: "CYP1A1", target: "2OHE2", type: "ACTIVATES", label: "2-hydroxylation", carcinogen: "E2" },
    { source: "4OHE2", target: "E2_quinone", type: "ACTIVATES", label: "auto-oxidation" },
    { source: "E2_quinone", target: "E2_depurin", type: "FORMS_ADDUCT", label: "depurinating adduct" },
    // Estrogen detoxification
    { source: "COMT", target: "E2_methyl", type: "DETOXIFIES", label: "O-methylation", carcinogen: "E2" },
    { source: "NQO1", target: "E2_quinone", type: "DETOXIFIES", label: "quinone reduction", carcinogen: "E2" },
    { source: "SULT1A1", target: "4OHE2", type: "DETOXIFIES", label: "sulfation", carcinogen: "E2" },
    { source: "UGT1A1", target: "4OHE2", type: "DETOXIFIES", label: "glucuronidation", carcinogen: "E2" },

    // ═══ ANDROGEN METABOLISM ═════════════════════════════════════════════
    // Biosynthesis: CYP17A1 → DHEA → HSD3B2 → A4 → AKR1C3 → Testosterone
    { source: "CYP17A1", target: "DHEA", type: "ACTIVATES", label: "17,20-lyase", carcinogen: "TESTO" },
    { source: "HSD3B2", target: "A4", type: "ACTIVATES", label: "3β-oxidation", carcinogen: "TESTO" },
    { source: "AKR1C3", target: "TESTO", type: "ACTIVATES", label: "17β-reduction", carcinogen: "TESTO" },
    // Testosterone → DHT (5α-reduction)
    { source: "SRD5A2", target: "DHT", type: "ACTIVATES", label: "5α-reduction (prostate)", carcinogen: "TESTO" },
    { source: "SRD5A1", target: "DHT", type: "ACTIVATES", label: "5α-reduction (peripheral)", carcinogen: "TESTO" },
    // Testosterone → Estradiol (aromatization — genotoxic bridge)
    { source: "CYP19A1", target: "E2_from_T", type: "ACTIVATES", label: "aromatization", carcinogen: "TESTO" },
    // E2 from T links to existing estrogen genotoxic pathway
    { source: "E2_from_T", target: "4OHE2", type: "ACTIVATES", label: "CYP1B1 → catechol" },
    { source: "E2_from_T", target: "E2_depurin", type: "FORMS_ADDUCT", label: "via quinone pathway" },
    // Testosterone hydroxylation (inactivation)
    { source: "CYP3A4", target: "6bOHT", type: "DETOXIFIES", label: "6β-hydroxylation", carcinogen: "TESTO" },
    { source: "CYP3A5", target: "6bOHT", type: "DETOXIFIES", label: "6β-hydroxylation", carcinogen: "TESTO" },
    // Testosterone glucuronidation (Phase II detox)
    { source: "UGT2B17", target: "T_gluc", type: "DETOXIFIES", label: "glucuronidation", carcinogen: "TESTO" },
    // DHT inactivation
    { source: "AKR1C2", target: "3aAdiol", type: "DETOXIFIES", label: "3α-reduction", carcinogen: "DHT" },
    { source: "UGT2B17", target: "DHT_gluc", type: "DETOXIFIES", label: "glucuronidation", carcinogen: "DHT" },
    { source: "UGT2B15", target: "DHT_gluc", type: "DETOXIFIES", label: "glucuronidation", carcinogen: "DHT" },
    // Transport of androgen conjugates
    { source: "ABCG2", target: "T_gluc", type: "TRANSPORTS", label: "efflux" },
    { source: "ABCB1", target: "DHT_gluc", type: "TRANSPORTS", label: "efflux" },
    // Androgen-driven oxidative DNA damage
    { source: "TESTO", target: "8oxodG", type: "FORMS_ADDUCT", label: "ROS via AR signaling" },

    // ═══ BENZENE METABOLISM ══════════════════════════════════════════
    // Benzene → CYP2E1 → benzene oxide → phenol → hydroquinone → benzoquinone → DNA
    { source: "CYP2E1", target: "BenzO", type: "ACTIVATES", label: "epoxidation", carcinogen: "BENZ" },
    { source: "BenzO", target: "HQ", type: "ACTIVATES", label: "rearrangement/hydroxylation" },
    { source: "HQ", target: "BQ", type: "ACTIVATES", label: "auto-oxidation" },
    { source: "BQ", target: "BQ_dG", type: "FORMS_ADDUCT", label: "quinone adducts" },
    // Benzene detoxification
    { source: "NQO1", target: "BQ", type: "DETOXIFIES", label: "quinone reduction", carcinogen: "BENZ" },
    { source: "GSTT1", target: "BenzO", type: "DETOXIFIES", label: "GSH conjugation", carcinogen: "BENZ" },
    { source: "XRCC1", target: "BQ_dG", type: "REPAIRS", label: "BER" },

    // ═══ VINYL CHLORIDE METABOLISM ═══════════════════════════════════
    // VC → CYP2E1 → chloroethylene oxide → etheno-DNA adducts
    { source: "CYP2E1", target: "CEO", type: "ACTIVATES", label: "epoxidation", carcinogen: "VCM" },
    { source: "CEO", target: "etheno_dA", type: "FORMS_ADDUCT", label: "εdA formation" },
    { source: "CEO", target: "etheno_dC", type: "FORMS_ADDUCT", label: "εdC formation" },
    { source: "GSTT1", target: "CEO", type: "DETOXIFIES", label: "GSH conjugation", carcinogen: "VCM" },
    { source: "XRCC1", target: "etheno_dA", type: "REPAIRS", label: "BER" },

    // ═══ PROCARCINOGEN → REACTIVE INTERMEDIATE (added) ════════════════
    // Wire each procarcinogen node to its first reactive intermediate so the
    // carcinogen is visually connected to its metabolism chain in the D3 viewer.
    { source: "BNZ", target: "NOHBNZ", type: "ACTIVATES", label: "procarcinogen → N-hydroxylation", carcinogen: "BNZ" },
    { source: "NDMA", target: "NDMA_alpha", type: "ACTIVATES", label: "procarcinogen → α-hydroxylation", carcinogen: "NDMA" },
    { source: "VCM", target: "CEO", type: "ACTIVATES", label: "procarcinogen → epoxidation", carcinogen: "VCM" },
    { source: "E2", target: "4OHE2", type: "ACTIVATES", label: "procarcinogen → catechol estrogen", carcinogen: "E2" },

    // ═══ DMBA METABOLISM CHAIN (added) ════════════════════════════════
    // DMBA → CYP1A1/1B1 → DMBA-3,4-epoxide → EPHX1 → DMBA-3,4-diol
    //      → CYP1A1/1B1 → DMBA-diol-epoxide → DNA adduct
    { source: "DMBA", target: "DMBA_epoxide", type: "ACTIVATES", label: "procarcinogen → 3,4-epoxidation", carcinogen: "DMBA" },
    { source: "CYP1A1", target: "DMBA_epoxide", type: "ACTIVATES", label: "3,4-epoxidation", carcinogen: "DMBA" },
    { source: "CYP1B1", target: "DMBA_epoxide", type: "ACTIVATES", label: "3,4-epoxidation", carcinogen: "DMBA" },
    { source: "EPHX1", target: "DMBA_diol", type: "DETOXIFIES", label: "epoxide hydrolysis", carcinogen: "DMBA" },
    { source: "CYP1A1", target: "DMBA_diol_epoxide", type: "ACTIVATES", label: "second epoxidation", carcinogen: "DMBA" },
    { source: "CYP1B1", target: "DMBA_diol_epoxide", type: "ACTIVATES", label: "second epoxidation", carcinogen: "DMBA" },
    { source: "DMBA_diol_epoxide", target: "DMBA_dA", type: "FORMS_ADDUCT", label: "bulky dA adduct" },
    { source: "XPC", target: "DMBA_dA", type: "REPAIRS", label: "NER recognition" },
    { source: "ERCC2", target: "DMBA_dA", type: "REPAIRS", label: "NER excision" },

    // ═══ MeIQx METABOLISM CHAIN (added) ═══════════════════════════════
    // MeIQx → CYP1A2 → N-OH-MeIQx → NAT1 → reactive acetoxy ester → C8-dG adduct
    { source: "MeIQx", target: "NOH_MeIQx", type: "ACTIVATES", label: "procarcinogen → N-hydroxylation", carcinogen: "MeIQx" },
    { source: "CYP1A2", target: "NOH_MeIQx", type: "ACTIVATES", label: "N-hydroxylation", carcinogen: "MeIQx" },
    { source: "NAT1", target: "MeIQx_acetoxy", type: "ACTIVATES", label: "O-acetylation", carcinogen: "MeIQx" },
    { source: "MeIQx_acetoxy", target: "PhIP_dG", type: "FORMS_ADDUCT", label: "HCA C8-dG adduct" },

    // ═══ CARCINOGEN → PATHWAY LINKS ═════════════════════════════════
    { source: "BaP", target: "KEGG5204", type: "PATHWAY", label: "PAH → DNA adducts" },
    { source: "BaP", target: "KEGG980", type: "PATHWAY", label: "CYP450 metabolism" },
    { source: "BENZ", target: "KEGG5208", type: "PATHWAY", label: "ROS carcinogenesis" },
    { source: "BENZ", target: "KEGG980", type: "PATHWAY", label: "CYP450 metabolism" },
    { source: "NNK", target: "KEGG5204", type: "PATHWAY", label: "DNA adduct pathway" },
    { source: "AFB1", target: "KEGG5204", type: "PATHWAY", label: "DNA adduct pathway" },
    { source: "PhIP", target: "KEGG980", type: "PATHWAY", label: "CYP450 metabolism" },
    { source: "GSTM1", target: "KEGG980G", type: "PATHWAY", label: "glutathione pathway" },
    { source: "GSTT1", target: "KEGG980G", type: "PATHWAY", label: "glutathione pathway" },
    { source: "CYP1A1", target: "KEGG980", type: "PATHWAY", label: "member" },
    { source: "CYP1A2", target: "KEGG982", type: "PATHWAY", label: "member" },
    { source: "CYP2B6", target: "KEGG982", type: "PATHWAY", label: "member" },
    { source: "CYP2C9", target: "KEGG982", type: "PATHWAY", label: "member" },
    { source: "CYP2C19", target: "KEGG982", type: "PATHWAY", label: "member" },
    { source: "CYP2D6", target: "KEGG982", type: "PATHWAY", label: "member" },
    { source: "CYP2E1", target: "KEGG980", type: "PATHWAY", label: "member" },
    { source: "CYP2F1", target: "KEGG980", type: "PATHWAY", label: "member" },
    { source: "CYP3A4", target: "KEGG982", type: "PATHWAY", label: "member" },
    // Androgen pathway links
    { source: "TESTO", target: "KEGG_ANDRO", type: "PATHWAY", label: "androgen signaling" },
    { source: "DHT", target: "KEGG_ANDRO", type: "PATHWAY", label: "AR activation" },
    { source: "CYP17A1", target: "KEGG_ANDRO", type: "PATHWAY", label: "steroidogenesis" },
    { source: "SRD5A2", target: "KEGG_ANDRO", type: "PATHWAY", label: "5α-reduction" },
    { source: "CYP19A1", target: "KEGG_ANDRO", type: "PATHWAY", label: "aromatization" },
    { source: "TESTO", target: "KEGG5208", type: "PATHWAY", label: "ROS carcinogenesis" },
  ]
};
