"""Population-scale synthetic simulation utilities.

Provides synthetic cohort generation, participant-level pipeline execution,
population summaries, genotype-exposure interaction analysis, and simple
validation helpers built on top of the package's flux, exposure, and
interaction engines.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from typing import Any, cast

from .phenotype_extractor import CANCER_PHENOTYPES as _PHENOTYPE_EXTRACTOR_CANCER_PHENOTYPES
from .validation_framework import REFERENCE_ORS as _VALIDATION_FRAMEWORK_REFERENCE_ORS

# ── Public constants ───────────────────────────────────────────────────────

ALLELE_FREQUENCIES: dict[str, dict[str, dict[str, float]]] = {
    "CYP2D6": {
        "European": {"*1/*1": 0.40, "*1/*4": 0.20, "*2/*4": 0.10, "*4/*4": 0.05, "*1/*2": 0.15, "*1/*41": 0.10},
        "African": {"*1/*1": 0.30, "*1/*17": 0.25, "*17/*17": 0.10, "*1/*4": 0.10, "*2/*17": 0.15, "*1/*2": 0.10},
        "East Asian": {"*1/*1": 0.20, "*1/*10": 0.35, "*10/*10": 0.15, "*1/*2": 0.15, "*1/*4": 0.05, "*2/*10": 0.10},
        "Admixed American": {"*1/*1": 0.34, "*1/*4": 0.18, "*2/*4": 0.10, "*4/*4": 0.04, "*1/*2": 0.18, "*1/*41": 0.16},
        "South Asian": {"*1/*1": 0.33, "*1/*4": 0.16, "*2/*4": 0.08, "*4/*4": 0.04, "*1/*2": 0.22, "*1/*41": 0.17},
        "Middle Eastern": {"*1/*1": 0.36, "*1/*4": 0.16, "*2/*4": 0.08, "*4/*4": 0.05, "*1/*2": 0.19, "*1/*41": 0.16},
        "Other": {"*1/*1": 0.35, "*1/*4": 0.17, "*2/*4": 0.09, "*4/*4": 0.04, "*1/*2": 0.18, "*1/*41": 0.17},
    },
    "CYP3A5": {
        "European": {"*3/*3": 0.85, "*1/*3": 0.12, "*1/*1": 0.03},
        "African": {"*1/*1": 0.35, "*1/*3": 0.40, "*3/*3": 0.25},
        "East Asian": {"*3/*3": 0.70, "*1/*3": 0.25, "*1/*1": 0.05},
        "Admixed American": {"*3/*3": 0.72, "*1/*3": 0.22, "*1/*1": 0.06},
        "South Asian": {"*3/*3": 0.68, "*1/*3": 0.24, "*1/*1": 0.08},
        "Middle Eastern": {"*3/*3": 0.75, "*1/*3": 0.20, "*1/*1": 0.05},
        "Other": {"*3/*3": 0.73, "*1/*3": 0.21, "*1/*1": 0.06},
    },
    "CYP2C19": {
        "European": {"*1/*1": 0.60, "*1/*2": 0.20, "*1/*17": 0.12, "*2/*2": 0.03, "*17/*17": 0.05},
        "African": {"*1/*1": 0.55, "*1/*2": 0.15, "*1/*17": 0.20, "*2/*2": 0.02, "*17/*17": 0.08},
        "East Asian": {"*1/*1": 0.30, "*1/*2": 0.30, "*1/*3": 0.15, "*2/*2": 0.10, "*2/*3": 0.10, "*3/*3": 0.05},
        "Admixed American": {"*1/*1": 0.50, "*1/*2": 0.20, "*1/*17": 0.14, "*2/*2": 0.04, "*17/*17": 0.12},
        "South Asian": {"*1/*1": 0.48, "*1/*2": 0.22, "*1/*17": 0.14, "*2/*2": 0.06, "*17/*17": 0.10},
        "Middle Eastern": {"*1/*1": 0.54, "*1/*2": 0.18, "*1/*17": 0.16, "*2/*2": 0.04, "*17/*17": 0.08},
        "Other": {"*1/*1": 0.51, "*1/*2": 0.19, "*1/*17": 0.16, "*2/*2": 0.04, "*17/*17": 0.10},
    },
    "CYP2C9": {
        "European": {"*1/*1": 0.70, "*1/*2": 0.14, "*1/*3": 0.10, "*2/*3": 0.03, "*3/*3": 0.03},
        "African": {"*1/*1": 0.72, "*1/*2": 0.05, "*1/*3": 0.03, "*1/*8": 0.12, "*8/*8": 0.08},
        "East Asian": {"*1/*1": 0.85, "*1/*3": 0.10, "*3/*3": 0.03, "*1/*13": 0.01, "*13/*13": 0.01},
        "Admixed American": {"*1/*1": 0.74, "*1/*2": 0.10, "*1/*3": 0.08, "*2/*3": 0.04, "*3/*3": 0.04},
        "South Asian": {"*1/*1": 0.68, "*1/*2": 0.10, "*1/*3": 0.15, "*2/*3": 0.04, "*3/*3": 0.03},
        "Middle Eastern": {"*1/*1": 0.71, "*1/*2": 0.12, "*1/*3": 0.11, "*2/*3": 0.03, "*3/*3": 0.03},
        "Other": {"*1/*1": 0.72, "*1/*2": 0.11, "*1/*3": 0.10, "*2/*3": 0.04, "*3/*3": 0.03},
    },
    "CYP2B6": {
        "European": {"*1/*1": 0.55, "*1/*6": 0.30, "*6/*6": 0.08, "*1/*4": 0.07},
        "African": {"*1/*1": 0.30, "*1/*6": 0.35, "*6/*6": 0.20, "*1/*18": 0.15},
        "East Asian": {"*1/*1": 0.50, "*1/*6": 0.35, "*6/*6": 0.10, "*1/*4": 0.05},
        "Admixed American": {"*1/*1": 0.48, "*1/*6": 0.33, "*6/*6": 0.12, "*1/*4": 0.07},
        "South Asian": {"*1/*1": 0.46, "*1/*6": 0.36, "*6/*6": 0.11, "*1/*4": 0.07},
        "Middle Eastern": {"*1/*1": 0.50, "*1/*6": 0.32, "*6/*6": 0.10, "*1/*4": 0.08},
        "Other": {"*1/*1": 0.49, "*1/*6": 0.33, "*6/*6": 0.11, "*1/*4": 0.07},
    },
    "GSTM1": {
        "European": {"null": 0.50, "present": 0.50},
        "African": {"null": 0.27, "present": 0.73},
        "East Asian": {"null": 0.55, "present": 0.45},
        "Admixed American": {"null": 0.40, "present": 0.60},
        "South Asian": {"null": 0.33, "present": 0.67},
        "Middle Eastern": {"null": 0.42, "present": 0.58},
        "Other": {"null": 0.41, "present": 0.59},
    },
    "GSTT1": {
        "European": {"null": 0.20, "present": 0.80},
        "African": {"null": 0.24, "present": 0.76},
        "East Asian": {"null": 0.47, "present": 0.53},
        "Admixed American": {"null": 0.23, "present": 0.77},
        "South Asian": {"null": 0.22, "present": 0.78},
        "Middle Eastern": {"null": 0.21, "present": 0.79},
        "Other": {"null": 0.24, "present": 0.76},
    },
    "ALDH2": {
        "European": {"*1/*1": 0.99, "*1/*2": 0.01, "*2/*2": 0.00},
        "African": {"*1/*1": 0.99, "*1/*2": 0.01, "*2/*2": 0.00},
        "East Asian": {"*1/*1": 0.55, "*1/*2": 0.35, "*2/*2": 0.10},
        "Admixed American": {"*1/*1": 0.95, "*1/*2": 0.04, "*2/*2": 0.01},
        "South Asian": {"*1/*1": 0.98, "*1/*2": 0.02, "*2/*2": 0.00},
        "Middle Eastern": {"*1/*1": 0.99, "*1/*2": 0.01, "*2/*2": 0.00},
        "Other": {"*1/*1": 0.96, "*1/*2": 0.03, "*2/*2": 0.01},
    },
    "NAT2": {
        "European": {"rapid": 0.40, "intermediate": 0.35, "slow": 0.25},
        "African": {"rapid": 0.55, "intermediate": 0.30, "slow": 0.15},
        "East Asian": {"rapid": 0.65, "intermediate": 0.25, "slow": 0.10},
        "Admixed American": {"rapid": 0.46, "intermediate": 0.33, "slow": 0.21},
        "South Asian": {"rapid": 0.52, "intermediate": 0.28, "slow": 0.20},
        "Middle Eastern": {"rapid": 0.44, "intermediate": 0.33, "slow": 0.23},
        "Other": {"rapid": 0.48, "intermediate": 0.31, "slow": 0.21},
    },
    "CYP1A1": {
        "European": {"WT": 0.82, "*1/*2A": 0.13, "*2A/*2A": 0.05},
        "African": {"WT": 0.92, "*1/*2A": 0.07, "*2A/*2A": 0.01},
        "East Asian": {"WT": 0.70, "*1/*2A": 0.20, "*2A/*2A": 0.10},
        "Admixed American": {"WT": 0.80, "*1/*2A": 0.14, "*2A/*2A": 0.06},
        "South Asian": {"WT": 0.76, "*1/*2A": 0.16, "*2A/*2A": 0.08},
        "Middle Eastern": {"WT": 0.81, "*1/*2A": 0.14, "*2A/*2A": 0.05},
        "Other": {"WT": 0.80, "*1/*2A": 0.14, "*2A/*2A": 0.06},
    },
    "CYP2E1": {
        "European": {"NM": 0.80, "UM_c1c1": 0.12, "IM": 0.08},
        "African": {"NM": 0.78, "UM_c1c1": 0.10, "IM": 0.12},
        "East Asian": {"NM": 0.68, "UM_c1c1": 0.22, "IM": 0.10},
        "Admixed American": {"NM": 0.78, "UM_c1c1": 0.14, "IM": 0.08},
        "South Asian": {"NM": 0.76, "UM_c1c1": 0.16, "IM": 0.08},
        "Middle Eastern": {"NM": 0.79, "UM_c1c1": 0.13, "IM": 0.08},
        "Other": {"NM": 0.78, "UM_c1c1": 0.14, "IM": 0.08},
    },
    # CYP1B1 L432V (*3) is a common gain-of-function variant that elevates
    # catalysis of PAH and estrogen-2/4-hydroxylation. Frequencies follow
    # 1000G phase 3 / Bailey 1998. Mapped onto the standard PM/IM/NM/RM scale:
    # NM = Leu/Leu reference, RM = Leu/Val (one *3 copy), UM = Val/Val (*3/*3).
    "CYP1B1": {
        "European": {"NM": 0.30, "RM": 0.50, "UM": 0.20},
        "African": {"NM": 0.10, "RM": 0.45, "UM": 0.45},
        "East Asian": {"NM": 0.55, "RM": 0.38, "UM": 0.07},
        "Admixed American": {"NM": 0.28, "RM": 0.48, "UM": 0.24},
        "South Asian": {"NM": 0.32, "RM": 0.49, "UM": 0.19},
        "Middle Eastern": {"NM": 0.31, "RM": 0.49, "UM": 0.20},
        "Other": {"NM": 0.31, "RM": 0.48, "UM": 0.21},
    },
    # EPHX1 Y113H/H139R confer slow / fast epoxide hydrolase phenotypes
    # (Hassett et al. 1994; Smith & Harrison 1997). PM = slow (Y113H/Y113H),
    # IM = intermediate, NM = wild type, RM = fast (H139R carriers).
    "EPHX1": {
        "European": {"PM": 0.06, "IM": 0.30, "NM": 0.50, "RM": 0.14},
        "African": {"PM": 0.04, "IM": 0.26, "NM": 0.55, "RM": 0.15},
        "East Asian": {"PM": 0.10, "IM": 0.36, "NM": 0.42, "RM": 0.12},
        "Admixed American": {"PM": 0.07, "IM": 0.31, "NM": 0.48, "RM": 0.14},
        "South Asian": {"PM": 0.07, "IM": 0.30, "NM": 0.49, "RM": 0.14},
        "Middle Eastern": {"PM": 0.06, "IM": 0.30, "NM": 0.50, "RM": 0.14},
        "Other": {"PM": 0.07, "IM": 0.30, "NM": 0.49, "RM": 0.14},
    },
    # NQO1*2 (Pro187Ser, rs1800566) homozygotes have ~5% activity (PM);
    # heterozygotes ~50% (IM). Frequencies from 1000G / Eguchi-Ishimae 2005.
    "NQO1": {
        "European": {"NM": 0.62, "IM": 0.32, "PM": 0.06},
        "African": {"NM": 0.66, "IM": 0.30, "PM": 0.04},
        "East Asian": {"NM": 0.30, "IM": 0.50, "PM": 0.20},
        "Admixed American": {"NM": 0.50, "IM": 0.40, "PM": 0.10},
        "South Asian": {"NM": 0.55, "IM": 0.36, "PM": 0.09},
        "Middle Eastern": {"NM": 0.60, "IM": 0.33, "PM": 0.07},
        "Other": {"NM": 0.57, "IM": 0.36, "PM": 0.07},
    },
    # GSTP1 Ile105Val (*B/*C alleles): Val carriers have reduced thermal
    # stability and 2-3x lower activity for several PAH-diol epoxide
    # substrates (Watson et al. 1998; Hu et al. 1997). NM = Ile/Ile,
    # IM = Ile/Val, PM = Val/Val.
    "GSTP1": {
        "European": {"NM": 0.43, "IM": 0.45, "PM": 0.12},
        "African": {"NM": 0.36, "IM": 0.48, "PM": 0.16},
        "East Asian": {"NM": 0.62, "IM": 0.33, "PM": 0.05},
        "Admixed American": {"NM": 0.45, "IM": 0.43, "PM": 0.12},
        "South Asian": {"NM": 0.48, "IM": 0.42, "PM": 0.10},
        "Middle Eastern": {"NM": 0.44, "IM": 0.44, "PM": 0.12},
        "Other": {"NM": 0.46, "IM": 0.43, "PM": 0.11},
    },
}

ANCESTRY_DISTRIBUTION: dict[str, float] = {
    "European": 0.45,
    "African": 0.22,
    "Admixed American": 0.18,
    "East Asian": 0.05,
    "South Asian": 0.04,
    "Middle Eastern": 0.03,
    "Other": 0.03,
}


HAPLOTYPE_BLOCKS: dict[str, dict[str, Any]] = {
    "CYP1_cluster": {
        "genes": ("CYP1A1", "CYP1A2"),
        "chromosome": "15q24.1",
        "reference": (
            "Landi et al. 2005; 1000 Genomes phase 3. CYP1A1 and CYP1A2 are "
            "adjacent (~25 kb) and share AhR-responsive regulatory elements; "
            "moderate-to-high LD reported across ancestries."
        ),
        "haplotypes": {
            "European": [
                {"p": 0.70, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1A"}},
                {"p": 0.12, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1F"}},
                {"p": 0.10, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1A/*1F"}},
                {"p": 0.05, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1F/*1F"}},
                {"p": 0.03, "genotypes": {"CYP1A1": "*2A/*2A", "CYP1A2": "*1F/*1F"}},
            ],
            "African": [
                {"p": 0.85, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1A"}},
                {"p": 0.07, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1F"}},
                {"p": 0.05, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1A/*1F"}},
                {"p": 0.02, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1F/*1F"}},
                {"p": 0.01, "genotypes": {"CYP1A1": "*2A/*2A", "CYP1A2": "*1F/*1F"}},
            ],
            "East Asian": [
                {"p": 0.55, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1A"}},
                {"p": 0.10, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1F"}},
                {"p": 0.15, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1A/*1F"}},
                {"p": 0.12, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1F/*1F"}},
                {"p": 0.08, "genotypes": {"CYP1A1": "*2A/*2A", "CYP1A2": "*1F/*1F"}},
            ],
            "Admixed American": [
                {"p": 0.67, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1A"}},
                {"p": 0.12, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1F"}},
                {"p": 0.11, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1A/*1F"}},
                {"p": 0.06, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1F/*1F"}},
                {"p": 0.04, "genotypes": {"CYP1A1": "*2A/*2A", "CYP1A2": "*1F/*1F"}},
            ],
            "South Asian": [
                {"p": 0.62, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1A"}},
                {"p": 0.13, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1F"}},
                {"p": 0.12, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1A/*1F"}},
                {"p": 0.08, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1F/*1F"}},
                {"p": 0.05, "genotypes": {"CYP1A1": "*2A/*2A", "CYP1A2": "*1F/*1F"}},
            ],
            "Middle Eastern": [
                {"p": 0.69, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1A"}},
                {"p": 0.12, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1F"}},
                {"p": 0.11, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1A/*1F"}},
                {"p": 0.05, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1F/*1F"}},
                {"p": 0.03, "genotypes": {"CYP1A1": "*2A/*2A", "CYP1A2": "*1F/*1F"}},
            ],
            "Other": [
                {"p": 0.68, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1A"}},
                {"p": 0.12, "genotypes": {"CYP1A1": "WT", "CYP1A2": "*1A/*1F"}},
                {"p": 0.11, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1A/*1F"}},
                {"p": 0.05, "genotypes": {"CYP1A1": "*1/*2A", "CYP1A2": "*1F/*1F"}},
                {"p": 0.04, "genotypes": {"CYP1A1": "*2A/*2A", "CYP1A2": "*1F/*1F"}},
            ],
        },
    },

    "GSTM_cluster": {
        "genes": ("GSTM1", "GSTM3"),
        "chromosome": "1p13.3",
        "reference": (
            "Inskip et al. 1995; Mitrunen et al. 2001. GSTM1 and GSTM3 are "
            "tightly linked (~140 kb) with coordinated null-deletion haplotypes "
            "in Europeans; frequency of linked null differs by ancestry."
        ),
        "haplotypes": {
            "European": [
                {"p": 0.40, "genotypes": {"GSTM1": "present", "GSTM3": "AA"}},
                {"p": 0.10, "genotypes": {"GSTM1": "present", "GSTM3": "AB"}},
                {"p": 0.35, "genotypes": {"GSTM1": "null", "GSTM3": "AA"}},
                {"p": 0.10, "genotypes": {"GSTM1": "null", "GSTM3": "AB"}},
                {"p": 0.05, "genotypes": {"GSTM1": "null", "GSTM3": "BB"}},
            ],
            "African": [
                {"p": 0.60, "genotypes": {"GSTM1": "present", "GSTM3": "AA"}},
                {"p": 0.13, "genotypes": {"GSTM1": "present", "GSTM3": "AB"}},
                {"p": 0.20, "genotypes": {"GSTM1": "null", "GSTM3": "AA"}},
                {"p": 0.05, "genotypes": {"GSTM1": "null", "GSTM3": "AB"}},
                {"p": 0.02, "genotypes": {"GSTM1": "null", "GSTM3": "BB"}},
            ],
            "East Asian": [
                {"p": 0.36, "genotypes": {"GSTM1": "present", "GSTM3": "AA"}},
                {"p": 0.09, "genotypes": {"GSTM1": "present", "GSTM3": "AB"}},
                {"p": 0.40, "genotypes": {"GSTM1": "null", "GSTM3": "AA"}},
                {"p": 0.10, "genotypes": {"GSTM1": "null", "GSTM3": "AB"}},
                {"p": 0.05, "genotypes": {"GSTM1": "null", "GSTM3": "BB"}},
            ],
            "Admixed American": [
                {"p": 0.50, "genotypes": {"GSTM1": "present", "GSTM3": "AA"}},
                {"p": 0.10, "genotypes": {"GSTM1": "present", "GSTM3": "AB"}},
                {"p": 0.28, "genotypes": {"GSTM1": "null", "GSTM3": "AA"}},
                {"p": 0.08, "genotypes": {"GSTM1": "null", "GSTM3": "AB"}},
                {"p": 0.04, "genotypes": {"GSTM1": "null", "GSTM3": "BB"}},
            ],
            "South Asian": [
                {"p": 0.54, "genotypes": {"GSTM1": "present", "GSTM3": "AA"}},
                {"p": 0.13, "genotypes": {"GSTM1": "present", "GSTM3": "AB"}},
                {"p": 0.22, "genotypes": {"GSTM1": "null", "GSTM3": "AA"}},
                {"p": 0.08, "genotypes": {"GSTM1": "null", "GSTM3": "AB"}},
                {"p": 0.03, "genotypes": {"GSTM1": "null", "GSTM3": "BB"}},
            ],
            "Middle Eastern": [
                {"p": 0.48, "genotypes": {"GSTM1": "present", "GSTM3": "AA"}},
                {"p": 0.10, "genotypes": {"GSTM1": "present", "GSTM3": "AB"}},
                {"p": 0.29, "genotypes": {"GSTM1": "null", "GSTM3": "AA"}},
                {"p": 0.09, "genotypes": {"GSTM1": "null", "GSTM3": "AB"}},
                {"p": 0.04, "genotypes": {"GSTM1": "null", "GSTM3": "BB"}},
            ],
            "Other": [
                {"p": 0.49, "genotypes": {"GSTM1": "present", "GSTM3": "AA"}},
                {"p": 0.10, "genotypes": {"GSTM1": "present", "GSTM3": "AB"}},
                {"p": 0.28, "genotypes": {"GSTM1": "null", "GSTM3": "AA"}},
                {"p": 0.09, "genotypes": {"GSTM1": "null", "GSTM3": "AB"}},
                {"p": 0.04, "genotypes": {"GSTM1": "null", "GSTM3": "BB"}},
            ],
        },
    },

    "NAT_cluster": {
        "genes": ("NAT1", "NAT2"),
        "chromosome": "8p22",
        "reference": (
            "Hein 2002; Sabbagh et al. 2008. NAT1 and NAT2 are ~170 kb apart "
            "and in moderate LD; co-inheritance affects aromatic-amine "
            "acetylation of bladder carcinogens."
        ),
        "haplotypes": {
            "European": [
                {"p": 0.22, "genotypes": {"NAT1": "*10/*10", "NAT2": "rapid"}},
                {"p": 0.20, "genotypes": {"NAT1": "*4/*10", "NAT2": "intermediate"}},
                {"p": 0.18, "genotypes": {"NAT1": "*4/*4", "NAT2": "intermediate"}},
                {"p": 0.25, "genotypes": {"NAT1": "*4/*4", "NAT2": "slow"}},
                {"p": 0.15, "genotypes": {"NAT1": "*14/*14", "NAT2": "slow"}},
            ],
            "African": [
                {"p": 0.35, "genotypes": {"NAT1": "*10/*10", "NAT2": "rapid"}},
                {"p": 0.22, "genotypes": {"NAT1": "*4/*10", "NAT2": "intermediate"}},
                {"p": 0.18, "genotypes": {"NAT1": "*4/*4", "NAT2": "intermediate"}},
                {"p": 0.15, "genotypes": {"NAT1": "*4/*4", "NAT2": "slow"}},
                {"p": 0.10, "genotypes": {"NAT1": "*14/*14", "NAT2": "slow"}},
            ],
            "East Asian": [
                {"p": 0.45, "genotypes": {"NAT1": "*10/*10", "NAT2": "rapid"}},
                {"p": 0.22, "genotypes": {"NAT1": "*4/*10", "NAT2": "intermediate"}},
                {"p": 0.15, "genotypes": {"NAT1": "*4/*4", "NAT2": "intermediate"}},
                {"p": 0.12, "genotypes": {"NAT1": "*4/*4", "NAT2": "slow"}},
                {"p": 0.06, "genotypes": {"NAT1": "*14/*14", "NAT2": "slow"}},
            ],
            "Admixed American": [
                {"p": 0.27, "genotypes": {"NAT1": "*10/*10", "NAT2": "rapid"}},
                {"p": 0.21, "genotypes": {"NAT1": "*4/*10", "NAT2": "intermediate"}},
                {"p": 0.17, "genotypes": {"NAT1": "*4/*4", "NAT2": "intermediate"}},
                {"p": 0.22, "genotypes": {"NAT1": "*4/*4", "NAT2": "slow"}},
                {"p": 0.13, "genotypes": {"NAT1": "*14/*14", "NAT2": "slow"}},
            ],
            "South Asian": [
                {"p": 0.33, "genotypes": {"NAT1": "*10/*10", "NAT2": "rapid"}},
                {"p": 0.22, "genotypes": {"NAT1": "*4/*10", "NAT2": "intermediate"}},
                {"p": 0.15, "genotypes": {"NAT1": "*4/*4", "NAT2": "intermediate"}},
                {"p": 0.20, "genotypes": {"NAT1": "*4/*4", "NAT2": "slow"}},
                {"p": 0.10, "genotypes": {"NAT1": "*14/*14", "NAT2": "slow"}},
            ],
            "Middle Eastern": [
                {"p": 0.25, "genotypes": {"NAT1": "*10/*10", "NAT2": "rapid"}},
                {"p": 0.21, "genotypes": {"NAT1": "*4/*10", "NAT2": "intermediate"}},
                {"p": 0.18, "genotypes": {"NAT1": "*4/*4", "NAT2": "intermediate"}},
                {"p": 0.23, "genotypes": {"NAT1": "*4/*4", "NAT2": "slow"}},
                {"p": 0.13, "genotypes": {"NAT1": "*14/*14", "NAT2": "slow"}},
            ],
            "Other": [
                {"p": 0.29, "genotypes": {"NAT1": "*10/*10", "NAT2": "rapid"}},
                {"p": 0.21, "genotypes": {"NAT1": "*4/*10", "NAT2": "intermediate"}},
                {"p": 0.17, "genotypes": {"NAT1": "*4/*4", "NAT2": "intermediate"}},
                {"p": 0.21, "genotypes": {"NAT1": "*4/*4", "NAT2": "slow"}},
                {"p": 0.12, "genotypes": {"NAT1": "*14/*14", "NAT2": "slow"}},
            ],
        },
    },

    "CYP2C_cluster": {
        "genes": ("CYP2C8", "CYP2C9", "CYP2C19"),
        "chromosome": "10q24",
        "reference": (
            "Goldstein 2001; Schaeffeler et al. 2002; 1000 Genomes phase 3. "
            "The 10q24 CYP2C cluster spans ~400 kb with haplotype-level LD "
            "across CYP2C8, CYP2C9, and CYP2C19."
        ),
        "haplotypes": {
            "European": [
                {"p": 0.40, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*1"}},
                {"p": 0.12, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*17"}},
                {"p": 0.12, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*2", "CYP2C19": "*1/*1"}},
                {"p": 0.15, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.11, "genotypes": {"CYP2C8": "*3/*3", "CYP2C9": "*2/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*1/*4", "CYP2C9": "*3/*3", "CYP2C19": "*2/*2"}},
            ],
            "African": [
                {"p": 0.50, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*1"}},
                {"p": 0.15, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*17"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*1/*2", "CYP2C9": "*1/*8", "CYP2C19": "*1/*1"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*2/*2", "CYP2C9": "*8/*8", "CYP2C19": "*1/*2"}},
                {"p": 0.08, "genotypes": {"CYP2C8": "*1/*4", "CYP2C9": "*1/*2", "CYP2C19": "*1/*17"}},
                {"p": 0.07, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*3", "CYP2C19": "*2/*17"}},
            ],
            "East Asian": [
                {"p": 0.55, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*1"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*3", "CYP2C19": "*1/*3"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*3/*3", "CYP2C19": "*2/*2"}},
                {"p": 0.08, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*13", "CYP2C19": "*2/*3"}},
                {"p": 0.07, "genotypes": {"CYP2C8": "*3/*3", "CYP2C9": "*13/*13", "CYP2C19": "*3/*3"}},
            ],
            "Admixed American": [
                {"p": 0.44, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*1"}},
                {"p": 0.13, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*17"}},
                {"p": 0.12, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*2", "CYP2C19": "*1/*1"}},
                {"p": 0.13, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*3/*3", "CYP2C9": "*2/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.08, "genotypes": {"CYP2C8": "*1/*4", "CYP2C9": "*3/*3", "CYP2C19": "*2/*2"}},
            ],
            "South Asian": [
                {"p": 0.40, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*1"}},
                {"p": 0.12, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*17"}},
                {"p": 0.14, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*2", "CYP2C19": "*1/*2"}},
                {"p": 0.14, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.12, "genotypes": {"CYP2C8": "*3/*3", "CYP2C9": "*2/*3", "CYP2C19": "*2/*2"}},
                {"p": 0.08, "genotypes": {"CYP2C8": "*1/*4", "CYP2C9": "*3/*3", "CYP2C19": "*17/*17"}},
            ],
            "Middle Eastern": [
                {"p": 0.42, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*1"}},
                {"p": 0.13, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*17"}},
                {"p": 0.12, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*2", "CYP2C19": "*1/*1"}},
                {"p": 0.14, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.11, "genotypes": {"CYP2C8": "*3/*3", "CYP2C9": "*2/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.08, "genotypes": {"CYP2C8": "*1/*4", "CYP2C9": "*3/*3", "CYP2C19": "*2/*2"}},
            ],
            "Other": [
                {"p": 0.43, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*1"}},
                {"p": 0.13, "genotypes": {"CYP2C8": "*1/*1", "CYP2C9": "*1/*1", "CYP2C19": "*1/*17"}},
                {"p": 0.12, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*2", "CYP2C19": "*1/*1"}},
                {"p": 0.13, "genotypes": {"CYP2C8": "*1/*3", "CYP2C9": "*1/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.10, "genotypes": {"CYP2C8": "*3/*3", "CYP2C9": "*2/*3", "CYP2C19": "*1/*2"}},
                {"p": 0.09, "genotypes": {"CYP2C8": "*1/*4", "CYP2C9": "*3/*3", "CYP2C19": "*2/*2"}},
            ],
        },
    },
}

LIFESTYLE_PROBABILITIES: dict[str, float] = {
    "smoking": 0.14,
    "heavy_alcohol": 0.06,
    "moderate_alcohol": 0.25,
    "occupational_exposure": 0.08,
}

EXPOSURE_SCENARIO_RULES: dict[str, str] = {
    "smoking+heavy_alcohol+occupational_exposure": "smoker_industrial_heavy_drinker",
    "smoking+heavy_alcohol": "smoker_heavy_drinker",
    "smoking+occupational_exposure": "smoker_industrial_worker",
    "smoking+moderate_alcohol": "smoker_moderate_drinker",
    "smoking": "smoker",
    "heavy_alcohol": "heavy_drinker",
    "moderate_alcohol": "moderate_drinker",
    "occupational_exposure": "industrial_worker",
    "baseline": "general_population",
}

CANCER_PHENOTYPES: dict[str, dict[str, Any]] = {
    "lung_cancer": {
        "icd10_codes": ["C34", "C34.0", "C34.1", "C34.9", "C33"],
        "tissue": "Lung",
        "carcinogen_classes": ["PAH", "Benzene", "Dioxin", "Nitrosamine", "HeavyMetal_Chromium", "HeavyMetal_Nickel", "HeavyMetal_Arsenic"],
        "base_rate": 0.04,
    },
    "liver_cancer": {
        "icd10_codes": ["C22", "C22.0", "C22.1", "C22.9"],
        "tissue": "Liver",
        "carcinogen_classes": ["Aflatoxin", "Aldehyde", "ChlorinatedSolvent", "HeavyMetal_Arsenic"],
        "base_rate": 0.008,
    },
    "bladder_cancer": {
        "icd10_codes": ["C67", "C67.0", "C67.9"],
        "tissue": "Bladder",
        "carcinogen_classes": ["PAH", "HCA", "HeavyMetal_Arsenic", "HeavyMetal_Cadmium"],
        "base_rate": 0.015,
    },
    "kidney_cancer": {
        "icd10_codes": ["C64", "C64.1", "C64.9", "C65"],
        "tissue": "Kidney",
        "carcinogen_classes": ["ChlorinatedSolvent_TCE", "HeavyMetal_Cadmium", "HeavyMetal_Arsenic"],
        "base_rate": 0.012,
    },
    "colorectal_cancer": {
        "icd10_codes": ["C18", "C19", "C20", "C18.9"],
        "tissue": "Colon",
        "carcinogen_classes": ["HCA", "PAH", "NDMA", "Nitrosamine"],
        "base_rate": 0.035,
    },
    "esophageal_cancer": {
        "icd10_codes": ["C15", "C15.4", "C15.5", "C15.9"],
        "tissue": "Esophagus",
        "carcinogen_classes": ["Aldehyde_Acetaldehyde", "Nitrosamine", "NDMA", "PAH"],
        "base_rate": 0.005,
    },
    "breast_cancer": {
        "icd10_codes": ["C50", "C50.1", "C50.9"],
        "tissue": "Breast",
        "carcinogen_classes": ["PAH", "Dioxin", "HCA"],
        "base_rate": 0.030,
    },
    "prostate_cancer": {
        "icd10_codes": ["C61"],
        "tissue": "Prostate",
        "carcinogen_classes": ["PAH", "HCA", "HeavyMetal_Cadmium", "HeavyMetal_Arsenic"],
        "base_rate": 0.028,
    },
    "head_neck_cancer": {
        "icd10_codes": ["C01", "C02", "C09", "C10", "C32"],
        "tissue": "Esophagus",
        "carcinogen_classes": ["PAH", "Aldehyde_Acetaldehyde", "Nitrosamine", "HeavyMetal_Nickel"],
        "base_rate": 0.010,
    },
    "leukemia": {
        "icd10_codes": ["C91", "C92", "C93", "C95"],
        "tissue": "Liver",
        "carcinogen_classes": ["Benzene", "Dioxin"],
        "base_rate": 0.009,
    },
}

GENOTYPE_CANCER_ORS: dict[str, dict[str, Any]] = {
    "GSTM1_null_smoking_lung": {
        "gene": "GSTM1",
        "risk_allele": "null",
        "exposure": "smoking",
        "cancer": "lung_cancer",
        "published_OR": 1.58,
        "ci": (1.21, 2.06),
        "source": "Ye et al. 2006",
    },
    "GSTM1_null_smoking_bladder": {
        "gene": "GSTM1",
        "risk_allele": "null",
        "exposure": "smoking",
        "cancer": "bladder_cancer",
        "published_OR": 1.53,
        "ci": (1.11, 2.12),
        "source": "Engel et al. 2002",
    },
    "CYP1A1_2A_smoking_lung": {
        "gene": "CYP1A1",
        "risk_allele": "*2A",
        "exposure": "smoking",
        "cancer": "lung_cancer",
        "published_OR": 2.36,
        "ci": (1.16, 4.81),
        "source": "Shi et al. 2008",
    },
    "GSTM1_CYP1A1_combined_lung": {
        "gene": "GSTM1+CYP1A1",
        "risk_allele": "null+*2A",
        "exposure": "smoking",
        "cancer": "lung_cancer",
        "published_OR": 2.87,
        "ci": (1.73, 4.76),
        "source": "Vineis et al. 2007",
    },
    "ALDH2_star2_alcohol_esophageal": {
        "gene": "ALDH2",
        "risk_allele": "*1/*2",
        "exposure": "heavy_alcohol",
        "cancer": "esophageal_cancer",
        "published_OR": 6.97,
        "ci": (4.36, 11.12),
        "source": "Yokoyama & Omori 2005",
    },
    "ALDH2_star2_homozygous_esophageal": {
        "gene": "ALDH2",
        "risk_allele": "*2/*2",
        "exposure": "heavy_alcohol",
        "cancer": "esophageal_cancer",
        "published_OR": 12.5,
        "ci": (6.0, 26.0),
        "source": "Yokoyama & Omori 2005",
    },
    "NAT2_slow_smoking_bladder": {
        "gene": "NAT2",
        "risk_allele": "slow",
        "exposure": "smoking",
        "cancer": "bladder_cancer",
        "published_OR": 1.51,
        "ci": (1.28, 1.78),
        "source": "Garcia-Closas et al. 2005",
    },
    "GSTT1_active_TCE_kidney": {
        "gene": "GSTT1",
        "risk_allele": "present",
        "exposure": "occupational_TCE",
        "cancer": "kidney_cancer",
        "published_OR": 1.88,
        "ci": (1.06, 3.33),
        "source": "Karami et al. 2012",
    },
}

_CANCER_BASE_RATES: dict[str, float] = {
    "lung_cancer": 0.04,
    "liver_cancer": 0.008,
    "bladder_cancer": 0.015,
    "kidney_cancer": 0.012,
    "colorectal_cancer": 0.035,
    "esophageal_cancer": 0.005,
    "breast_cancer": 0.030,
    "prostate_cancer": 0.028,
    "head_neck_cancer": 0.010,
    "leukemia": 0.009,
}


def _build_canonical_cancer_phenotypes() -> dict[str, dict[str, Any]]:
    """Use phenotype_extractor as the canonical catalog, with base rates layered in."""
    merged: dict[str, dict[str, Any]] = {}
    for name, definition in _PHENOTYPE_EXTRACTOR_CANCER_PHENOTYPES.items():
        record = dict(definition)
        base_rate = _CANCER_BASE_RATES.get(name)
        if base_rate is not None:
            record["base_rate"] = base_rate
        merged[name] = record
    return merged


def _build_canonical_reference_ors() -> dict[str, dict[str, Any]]:
    """Normalize validation-framework references into the typed simulation shape."""
    normalized: dict[str, dict[str, Any]] = {}
    for name, reference in _VALIDATION_FRAMEWORK_REFERENCE_ORS.items():
        normalized[name] = {
            "gene": reference["gene"],
            "risk_allele": reference["risk_allele"],
            "exposure": reference["exposure"],
            "cancer": reference["cancer"],
            "published_OR": float(reference["published_OR"]),
            "ci": (float(reference["ci_low"]), float(reference["ci_high"])),
            "source": str(reference["source"]),
        }
    return normalized


# Use the richer phenotype/reference catalogs as the package canonical source.
CANCER_PHENOTYPES = _build_canonical_cancer_phenotypes()
GENOTYPE_CANCER_ORS = _build_canonical_reference_ors()


# ── Dataclasses ────────────────────────────────────────────────────────────


@dataclass
class SyntheticParticipant:
    """A synthetic participant with genotype, ancestry, and lifestyle."""

    person_id: str
    genotypes: dict[str, str]
    ancestry: str
    ancestry_label: str
    lifestyle: dict[str, bool]
    exposure_scenario: str
    tissue: str = "Liver"


@dataclass
class ParticipantRiskSummary:
    """Compact risk summary for a single participant."""

    person_id: str
    ancestry: str
    exposure_scenario: str
    interaction_factor: float | None
    total_independent_risk: float | None
    total_interaction_risk: float | None
    gsh_fraction: float | None
    gsh_tipping_point: bool
    flux_classes: dict[str, dict[str, Any]]
    critical_warnings: list[str]
    key_genotypes: dict[str, str]
    error: str | None = None


@dataclass
class PopulationStats:
    """Descriptive statistics for a population simulation."""

    n_total: int
    ancestry_distribution: dict[str, int]
    exposure_scenario_distribution: dict[str, int]
    interaction_factor_stats: dict[str, float]
    gsh_fraction_stats: dict[str, float]
    synergistic_count: int
    synergistic_fraction: float
    tipping_point_count: int
    tipping_point_fraction: float
    high_risk_pathway_counts: dict[str, int]


@dataclass
class GxEAnalysisResult:
    """2x2 genotype x exposure interaction analysis."""

    interaction_name: str
    gene: str
    risk_allele: str
    exposure: str
    cancer_type: str
    group_a_count: int
    group_b_count: int
    group_c_count: int
    group_d_count: int
    observed_or: float | None
    or_95ci: tuple[float, float] | None
    published_or: float
    published_ci: tuple[float, float]
    concordant: bool | None
    source: str


@dataclass
class ValidationResult:
    """ROC/AUC and calibration validation results."""

    roc_auc: float
    roc_auc_se: float
    roc_auc_95ci: tuple[float, float]
    n_cases: int
    n_controls: int
    calibration_bins: list[dict[str, float]]
    hosmer_lemeshow_chi2: float | None


@dataclass
class PopulationSimulationResult:
    """Complete population simulation output."""

    n_participants: int
    n_completed: int
    n_errors: int
    participant_summaries: list[ParticipantRiskSummary]
    population_stats: PopulationStats
    gxe_analyses: list[GxEAnalysisResult]
    validation: dict[str, ValidationResult] | None
    or_comparison: list[GxEAnalysisResult] | None
    elapsed_seconds: float


# ── Private helpers ────────────────────────────────────────────────────────


_SUMMARY_GENES = ("GSTM1", "GSTT1", "ALDH2", "NAT2", "CYP1A1", "CYP2E1")


def _weighted_pick(rng: random.Random, distribution: dict[str, float]) -> str:
    """Pick one label from a normalized categorical distribution."""
    threshold = rng.random()
    cumulative = 0.0
    last_key = next(iter(distribution))
    for key, probability in distribution.items():
        cumulative += probability
        if threshold <= cumulative:
            return key
        last_key = key
    return last_key


def sample_haplotype_block(
    rng: random.Random,
    block_name: str,
    ancestry: str,
) -> dict[str, str]:
    """Sample a curated haplotype for a gene cluster from ``HAPLOTYPE_BLOCKS``."""
    block = HAPLOTYPE_BLOCKS[block_name]
    per_ancestry = block["haplotypes"]
    candidates = per_ancestry.get(ancestry) or per_ancestry.get("European") or []
    if not candidates:
        return {}
    total = sum(float(entry.get("p", 0.0)) for entry in candidates) or 1.0
    threshold = rng.random() * total
    cumulative = 0.0
    last_entry = candidates[-1]
    for entry in candidates:
        cumulative += float(entry.get("p", 0.0))
        if threshold <= cumulative:
            return dict(entry["genotypes"])
        last_entry = entry
    return dict(last_entry["genotypes"])


def _genes_covered_by_blocks() -> set[str]:
    """Return the set of gene names covered by any haplotype block."""
    covered: set[str] = set()
    for block in HAPLOTYPE_BLOCKS.values():
        for gene in block["genes"]:
            covered.add(gene)
    return covered


def _percentile(sorted_values: list[float], q: float) -> float:
    """Compute a simple percentile for a sorted sequence."""
    if not sorted_values:
        return 0.0
    index = min(len(sorted_values) - 1, max(0, int(round((len(sorted_values) - 1) * q))))
    return sorted_values[index]


def _compute_numeric_stats(values: list[float]) -> dict[str, float]:
    """Compute n, mean, median, sd, min, max, p25, p75, p95."""
    if not values:
        return {}
    sorted_values = sorted(values)
    n = len(sorted_values)
    mean = sum(sorted_values) / n
    mid = n // 2
    if n % 2:
        median = sorted_values[mid]
    else:
        median = (sorted_values[mid - 1] + sorted_values[mid]) / 2
    variance = sum((value - mean) ** 2 for value in sorted_values) / n if n > 1 else 0.0
    return {
        "n": float(n),
        "mean": round(mean, 4),
        "median": round(median, 4),
        "sd": round(math.sqrt(variance), 4),
        "min": round(sorted_values[0], 4),
        "max": round(sorted_values[-1], 4),
        "p25": round(_percentile(sorted_values, 0.25), 4),
        "p75": round(_percentile(sorted_values, 0.75), 4),
        "p95": round(_percentile(sorted_values, 0.95), 4),
    }


def _derive_exposure_scenario(lifestyle: dict[str, bool]) -> str:
    """Map lifestyle factors to a synthetic exposure scenario label."""
    smoking = lifestyle.get("smoking", False)
    heavy_alcohol = lifestyle.get("alcohol_heavy", False)
    moderate_alcohol = lifestyle.get("alcohol_moderate", False)
    occupational = lifestyle.get("occupational_exposure", False)

    if smoking and heavy_alcohol and occupational:
        return "smoker_industrial_heavy_drinker"
    if smoking and heavy_alcohol:
        return "smoker_heavy_drinker"
    if smoking and occupational:
        return "smoker_industrial_worker"
    if smoking and moderate_alcohol:
        return "smoker_moderate_drinker"
    if smoking:
        return "smoker"
    if heavy_alcohol:
        return "heavy_drinker"
    if occupational:
        return "industrial_worker"
    if moderate_alcohol:
        return "moderate_drinker"
    return "general_population"


def _build_exposure_answers(
    exposure_scenario: str,
    lifestyle: dict[str, bool],
) -> dict[str, Any]:
    """Build minimal exposure-questionnaire answers for the exposure engine."""
    smoking = lifestyle.get("smoking", False)
    heavy_alcohol = lifestyle.get("alcohol_heavy", False)
    moderate_alcohol = lifestyle.get("alcohol_moderate", False)
    occupational = lifestyle.get("occupational_exposure", False)

    answers: dict[str, Any] = {
        "smoking_status": "current" if smoking else "never",
        "cigarettes_per_day": 25 if smoking else 0,
        "secondhand_smoke": False,
        "alcohol_status": "heavy" if heavy_alcohol else "moderate" if moderate_alcohol else "nondrinker",
        "drinks_per_week": 20 if heavy_alcohol else 8 if moderate_alcohol else 0,
        "grilled_meat_frequency": "weekly",
        "processed_meat_frequency": "monthly",
        "red_meat_servings": 4,
        "seafood_frequency": 1,
        "aflatoxin_exposure": False,
        "occupation_category": "manufacturing_industrial" if occupational else "office_indoor",
        "occupational_chemicals": [],
        "urban_rural": "suburban",
        "proximity_industrial": occupational,
        "water_source": "municipal_treated",
        "arsenic_endemic_area": False,
        "sun_exposure": "moderate",
    }

    if occupational:
        answers["occupational_chemicals"] = ["benzene_petroleum", "formaldehyde", "chromium_vi", "TCE_solvents"]
    if "heavy_drinker" in exposure_scenario:
        answers["processed_meat_frequency"] = "daily"
    if "smoker" in exposure_scenario:
        answers["secondhand_smoke"] = True
        answers["grilled_meat_frequency"] = "daily"

    return answers


def _build_interaction_exposure_profile(
    exposure_scenario: str,
    lifestyle: dict[str, bool],
) -> dict[str, float]:
    """Build the compact carcinogen multiplier profile used by interaction_engine."""
    profile: dict[str, float] = {}

    if lifestyle.get("smoking", False) or "smoker" in exposure_scenario:
        profile.update(
            {
                "PAH": 3.0,
                "NNK": 4.0,
                "HCA": 1.5,
                "benzene": 6.0,
                "formaldehyde": 2.0,
                "cadmium": 2.0,
                "acrolein": 5.0,
            }
        )

    if lifestyle.get("alcohol_moderate", False) or "moderate_drinker" in exposure_scenario:
        profile["acetaldehyde"] = max(profile.get("acetaldehyde", 1.0), 1.5)
        profile["NDMA"] = max(profile.get("NDMA", 1.0), 1.2)
        profile["ethanol"] = max(profile.get("ethanol", 1.0), 3.0)

    if lifestyle.get("alcohol_heavy", False) or "heavy_drinker" in exposure_scenario:
        profile["acetaldehyde"] = max(profile.get("acetaldehyde", 1.0), 3.0)
        profile["NDMA"] = max(profile.get("NDMA", 1.0), 2.0)
        profile["benzene"] = max(profile.get("benzene", 1.0), 1.5)
        profile["ethanol"] = max(profile.get("ethanol", 1.0), 8.0)
        profile["acrolein"] = max(profile.get("acrolein", 1.0), 1.5)

    if lifestyle.get("occupational_exposure", False) or "industrial" in exposure_scenario:
        profile["chromium_VI"] = max(profile.get("chromium_VI", 1.0), 10.0)
        profile["benzene"] = max(profile.get("benzene", 1.0), 10.0)
        profile["formaldehyde"] = max(profile.get("formaldehyde", 1.0), 5.0)
        profile["PAH"] = max(profile.get("PAH", 1.0), 2.0)
        profile["cadmium"] = max(profile.get("cadmium", 1.0), 3.0)
        profile["vinyl_chloride"] = max(profile.get("vinyl_chloride", 1.0), 2.0)

    return profile or {"PAH": 1.0}


def _extract_flux_classes(profile_result: Any) -> dict[str, dict[str, Any]]:
    """Convert the flux profile dataclasses into a compact JSON-like summary."""
    if profile_result is None:
        return {}
    flux_classes: dict[str, dict[str, Any]] = {}
    for class_name, result in profile_result.per_class_results.items():
        flux_classes[class_name] = {
            "net_ratio": result.net_ratio,
            "susceptibility_score_log2": result.susceptibility_score_log2,
            "reactive_intermediate_uM": result.steady_state_concentrations_uM.get(
                "reactive_intermediate_uM"
            ),
            "time_to_steady_state_days": result.steady_state_model.get(
                "time_to_steady_state_days"
            ),
            "risk_classification": result.risk_classification.value,
            "model_kind": getattr(result, "model_kind", None),
            "parameter_source": getattr(result, "parameter_source", None),
        }
    return flux_classes


def _score_value(summary: ParticipantRiskSummary, score_key: str) -> float | None:
    """Extract a numeric validation score from a participant summary."""
    if score_key == "interaction_factor":
        return summary.interaction_factor
    if score_key == "total_interaction_risk":
        return summary.total_interaction_risk
    if score_key == "total_independent_risk":
        return summary.total_independent_risk
    if score_key == "high_risk_pathway_count":
        return float(
            sum(
                1
                for data in summary.flux_classes.values()
                if data.get("risk_classification") in {"ELEVATED", "HIGH"}
            )
        )
    return None


def _scenario_has_exposure(exposure_scenario: str, exposure_name: str) -> bool:
    """Determine whether a scenario implies the reference exposure."""
    scenario = exposure_scenario.lower()
    if exposure_name == "smoking":
        return "smoker" in scenario
    if exposure_name == "heavy_alcohol":
        return "heavy_drinker" in scenario
    if exposure_name == "moderate_alcohol":
        return "moderate_drinker" in scenario
    if exposure_name == "occupational_TCE":
        return "industrial" in scenario
    return exposure_name.lower() in scenario


def _matches_risk_allele(genotype_value: str, risk_allele: str) -> bool:
    """Return whether a genotype string matches a risk-allele definition."""
    genotype = genotype_value.lower()
    risk = risk_allele.lower()
    if genotype == risk:
        return True
    if "+" in risk:
        return False
    if risk.startswith("*") and risk in genotype:
        return True
    return False


def _summary_has_risk(summary: ParticipantRiskSummary, gene: str, risk_allele: str) -> bool:
    """Return whether a summary carries the requested single or combined risk."""
    if "+" in gene:
        genes = gene.split("+")
        alleles = risk_allele.split("+")
        return all(
            _matches_risk_allele(summary.key_genotypes.get(sub_gene, ""), sub_risk)
            for sub_gene, sub_risk in zip(genes, alleles)
        )
    return _matches_risk_allele(summary.key_genotypes.get(gene, ""), risk_allele)


def _label_has_cancer(
    cancer_labels: dict[str, dict[str, Any]],
    person_id: str,
    cancer_type: str,
) -> bool:
    """Return whether the label dict marks a given participant as a case."""
    label = cancer_labels.get(person_id, {})
    if cancer_type == "any_cancer":
        return bool(label.get("any_cancer", False))
    return cancer_type in label.get("cancer_types", [])


def _analyze_2x2(
    summaries: list[ParticipantRiskSummary],
    cancer_labels: dict[str, dict[str, Any]] | None,
    *,
    interaction_name: str,
    gene: str,
    risk_allele: str,
    exposure: str,
    cancer_type: str,
    published_or: float,
    published_ci: tuple[float, float],
    source: str,
) -> GxEAnalysisResult:
    """Compute a 2x2 GxE analysis with OR and 95% CI (Woolf's method)."""
    group_a: list[ParticipantRiskSummary] = []
    group_b: list[ParticipantRiskSummary] = []
    group_c: list[ParticipantRiskSummary] = []
    group_d: list[ParticipantRiskSummary] = []

    for summary in summaries:
        has_risk = _summary_has_risk(summary, gene, risk_allele)
        has_exposure = _scenario_has_exposure(summary.exposure_scenario, exposure)
        if has_risk and has_exposure:
            group_a.append(summary)
        elif has_risk:
            group_b.append(summary)
        elif has_exposure:
            group_c.append(summary)
        else:
            group_d.append(summary)

    observed_or: float | None = None
    interval: tuple[float, float] | None = None
    concordant: bool | None = None

    if cancer_labels is not None and group_a and group_d:
        a_cases = sum(1 for summary in group_a if _label_has_cancer(cancer_labels, summary.person_id, cancer_type))
        d_cases = sum(1 for summary in group_d if _label_has_cancer(cancer_labels, summary.person_id, cancer_type))
        a_controls = len(group_a) - a_cases
        d_controls = len(group_d) - d_cases

        a = a_cases + 0.5
        b = a_controls + 0.5
        c = d_cases + 0.5
        d = d_controls + 0.5
        odds_ratio = (a * d) / (b * c)
        log_or = math.log(odds_ratio)
        standard_error = math.sqrt((1 / a) + (1 / b) + (1 / c) + (1 / d))
        interval = (
            round(math.exp(log_or - 1.96 * standard_error), 3),
            round(math.exp(log_or + 1.96 * standard_error), 3),
        )
        observed_or = round(odds_ratio, 3)
        concordant = published_ci[0] <= observed_or <= published_ci[1]

    return GxEAnalysisResult(
        interaction_name=interaction_name,
        gene=gene,
        risk_allele=risk_allele,
        exposure=exposure,
        cancer_type=cancer_type,
        group_a_count=len(group_a),
        group_b_count=len(group_b),
        group_c_count=len(group_c),
        group_d_count=len(group_d),
        observed_or=observed_or,
        or_95ci=interval,
        published_or=published_or,
        published_ci=published_ci,
        concordant=concordant,
        source=source,
    )


# ── Public API ────────────────────────────────────────────────────────────


def generate_synthetic_cohort(
    n: int = 1000,
    *,
    seed: int = 42,
    use_haplotypes: bool = False,
) -> list[SyntheticParticipant]:
    """Generate synthetic All of Us-like cohort using published allele frequencies.

    When ``use_haplotypes=True`` the four curated LD/haplotype blocks in
    ``HAPLOTYPE_BLOCKS`` (CYP1_cluster, GSTM_cluster, NAT_cluster,
    CYP2C_cluster) are sampled jointly — preserving within-block co-inheritance
    and introducing additional linked genes (CYP1A2, GSTM3, NAT1, CYP2C8) that
    are not present in the independent ``ALLELE_FREQUENCIES`` table. Genes
    outside those blocks continue to be sampled marginally.
    """
    rng = random.Random(seed)
    participants: list[SyntheticParticipant] = []
    haplotype_block_names = list(HAPLOTYPE_BLOCKS.keys())
    covered_genes = _genes_covered_by_blocks() if use_haplotypes else set()

    for index in range(n):
        ancestry = _weighted_pick(rng, ANCESTRY_DISTRIBUTION)
        if use_haplotypes:
            genotypes: dict[str, str] = {}
            for block_name in haplotype_block_names:
                genotypes.update(sample_haplotype_block(rng, block_name, ancestry))
            for gene, per_ancestry in ALLELE_FREQUENCIES.items():
                if gene in covered_genes:
                    continue
                genotypes[gene] = _weighted_pick(rng, per_ancestry[ancestry])
        else:
            genotypes = {
                gene: _weighted_pick(rng, per_ancestry[ancestry])
                for gene, per_ancestry in ALLELE_FREQUENCIES.items()
            }
        smoking = rng.random() < LIFESTYLE_PROBABILITIES["smoking"]
        heavy_alcohol = rng.random() < LIFESTYLE_PROBABILITIES["heavy_alcohol"]
        moderate_alcohol = (not heavy_alcohol) and rng.random() < LIFESTYLE_PROBABILITIES["moderate_alcohol"]
        occupational = rng.random() < LIFESTYLE_PROBABILITIES["occupational_exposure"]
        lifestyle = {
            "smoking": smoking,
            "alcohol_heavy": heavy_alcohol,
            "alcohol_moderate": moderate_alcohol,
            "occupational_exposure": occupational,
        }
        participants.append(
            SyntheticParticipant(
                person_id=f"SYN{100000 + index}",
                genotypes=genotypes,
                ancestry=ancestry,
                ancestry_label=ancestry,
                lifestyle=lifestyle,
                exposure_scenario=_derive_exposure_scenario(lifestyle),
            )
        )

    return participants


def generate_synthetic_cancer_labels(
    participants: list[SyntheticParticipant],
    *,
    base_rate: float = 0.10,
    genotype_effects: bool = True,
    seed: int = 42,
) -> dict[str, dict[str, Any]]:
    """Generate synthetic cancer labels using published genotype-cancer ORs."""
    rng = random.Random(seed)
    labels: dict[str, dict[str, Any]] = {}
    scale = base_rate / 0.10 if base_rate > 0 else 0.0

    refs_by_cancer: dict[str, list[dict[str, Any]]] = {}
    for ref in GENOTYPE_CANCER_ORS.values():
        refs_by_cancer.setdefault(ref["cancer"], []).append(ref)

    for participant in participants:
        cancers: list[str] = []
        for cancer_type, definition in CANCER_PHENOTYPES.items():
            base_rate_value = definition.get("base_rate")
            if cancer_type == "any_cancer" or base_rate_value is None:
                continue
            rate = float(base_rate_value) * scale
            if genotype_effects:
                for ref in refs_by_cancer.get(cancer_type, []):
                    has_risk = _matches_risk_allele(
                        participant.genotypes.get(ref["gene"].split("+")[0], ""),
                        ref["risk_allele"].split("+")[0],
                    ) if "+" not in ref["gene"] else all(
                        _matches_risk_allele(participant.genotypes.get(gene_name, ""), risk)
                        for gene_name, risk in zip(ref["gene"].split("+"), ref["risk_allele"].split("+"))
                    )
                    has_exposure = _scenario_has_exposure(participant.exposure_scenario, ref["exposure"])
                    if has_risk and has_exposure:
                        rate *= ref["published_OR"]
                    elif has_risk:
                        rate *= ref["published_OR"] ** 0.5
                    elif has_exposure:
                        rate *= ref["published_OR"] ** 0.7
            if rng.random() < min(rate, 0.95):
                cancers.append(cancer_type)

        labels[participant.person_id] = {
            "any_cancer": bool(cancers),
            "cancer_types": cancers,
            "n_cancer_diagnoses": len(cancers),
        }

    return labels


def run_participant_risk(
    participant: SyntheticParticipant,
) -> ParticipantRiskSummary:
    """Run full risk assessment for a single participant."""
    error_messages: list[str] = []
    flux_result: Any | None = None
    exposure_result: Any | None = None
    interaction_result: Any | None = None
    warnings: list[str] = []

    try:
        from ..flux_engine import compute_full_profile

        flux_result = compute_full_profile(participant.genotypes, participant.tissue)
    except Exception as exc:  # pragma: no cover - defensive path
        error_messages.append(f"flux:{type(exc).__name__}:{exc}")

    try:
        from ..exposure_engine import apply_exposure_profile

        exposure_answers = _build_exposure_answers(participant.exposure_scenario, participant.lifestyle)
        exposure_result = apply_exposure_profile(
            patient_genotypes=participant.genotypes,
            exposure_answers=exposure_answers,
            tissue=participant.tissue,
        )
    except Exception as exc:  # pragma: no cover - defensive path
        error_messages.append(f"exposure:{type(exc).__name__}:{exc}")

    try:
        from ..interaction_engine import compute_interaction_matrix, identify_critical_interactions

        interaction_result = compute_interaction_matrix(
            cast(
                dict[str, float | dict[str, Any]],
                _build_interaction_exposure_profile(
                    participant.exposure_scenario,
                    participant.lifestyle,
                ),
            ),
            genotypes=participant.genotypes,
            tissue=participant.tissue,
            lifestyle=participant.lifestyle,
        )
        warnings = [f"[{item.severity}] {item.interaction}" for item in identify_critical_interactions(participant.genotypes)]
    except Exception as exc:  # pragma: no cover - defensive path
        error_messages.append(f"interaction:{type(exc).__name__}:{exc}")

    flux_classes = _extract_flux_classes(flux_result)
    if exposure_result is not None:
        for risk in exposure_result.risk_scores:
            flux_classes.setdefault(
                risk.carcinogen_class,
                {
                    "net_ratio": None,
                    "risk_classification": risk.risk_category.value,
                    "model_kind": "exposure_only",
                    "parameter_source": "",
                },
            )

    return ParticipantRiskSummary(
        person_id=participant.person_id,
        ancestry=participant.ancestry_label,
        exposure_scenario=participant.exposure_scenario,
        interaction_factor=interaction_result.interaction_factor if interaction_result is not None else None,
        total_independent_risk=interaction_result.total_independent_risk if interaction_result is not None else None,
        total_interaction_risk=interaction_result.total_interaction_risk if interaction_result is not None else None,
        gsh_fraction=interaction_result.gsh_status.fraction_normal if interaction_result is not None else None,
        gsh_tipping_point=interaction_result.gsh_status.tipping_point_reached if interaction_result is not None else False,
        flux_classes=flux_classes,
        critical_warnings=warnings,
        key_genotypes={gene: participant.genotypes.get(gene, "") for gene in _SUMMARY_GENES},
        error="; ".join(error_messages) if error_messages else None,
    )


def compute_population_stats(
    summaries: list[ParticipantRiskSummary],
) -> PopulationStats:
    """Compute population-level descriptive statistics."""
    ancestry_distribution: dict[str, int] = {}
    scenario_distribution: dict[str, int] = {}
    interaction_values: list[float] = []
    gsh_values: list[float] = []
    high_risk_pathway_counts: dict[str, int] = {}

    for summary in summaries:
        ancestry_distribution[summary.ancestry] = ancestry_distribution.get(summary.ancestry, 0) + 1
        scenario_distribution[summary.exposure_scenario] = scenario_distribution.get(summary.exposure_scenario, 0) + 1
        if summary.interaction_factor is not None:
            interaction_values.append(summary.interaction_factor)
        if summary.gsh_fraction is not None:
            gsh_values.append(summary.gsh_fraction)
        for class_name, data in summary.flux_classes.items():
            if data.get("risk_classification") in {"ELEVATED", "HIGH"}:
                high_risk_pathway_counts[class_name] = high_risk_pathway_counts.get(class_name, 0) + 1

    synergistic_count = sum(1 for value in interaction_values if value > 1.2)
    tipping_point_count = sum(1 for summary in summaries if summary.gsh_tipping_point)
    top_high_risk = dict(
        sorted(high_risk_pathway_counts.items(), key=lambda item: (-item[1], item[0]))[:10]
    )

    return PopulationStats(
        n_total=len(summaries),
        ancestry_distribution=dict(sorted(ancestry_distribution.items(), key=lambda item: (-item[1], item[0]))),
        exposure_scenario_distribution=dict(sorted(scenario_distribution.items(), key=lambda item: (-item[1], item[0]))),
        interaction_factor_stats=_compute_numeric_stats(interaction_values),
        gsh_fraction_stats=_compute_numeric_stats(gsh_values),
        synergistic_count=synergistic_count,
        synergistic_fraction=round(synergistic_count / len(interaction_values), 4) if interaction_values else 0.0,
        tipping_point_count=tipping_point_count,
        tipping_point_fraction=round(tipping_point_count / len(summaries), 4) if summaries else 0.0,
        high_risk_pathway_counts=top_high_risk,
    )


def analyze_gxe_interactions(
    summaries: list[ParticipantRiskSummary],
    cancer_labels: dict[str, dict[str, Any]] | None = None,
) -> list[GxEAnalysisResult]:
    """Analyze genotype x exposure interactions with 2x2 tables and ORs."""
    results: list[GxEAnalysisResult] = []
    for name, ref in GENOTYPE_CANCER_ORS.items():
        results.append(
            _analyze_2x2(
                summaries,
                cancer_labels,
                interaction_name=name,
                gene=ref["gene"],
                risk_allele=ref["risk_allele"],
                exposure=ref["exposure"],
                cancer_type=ref["cancer"],
                published_or=float(ref["published_OR"]),
                published_ci=tuple(ref["ci"]),
                source=str(ref["source"]),
            )
        )
    return results


def validate_against_published_ors(
    summaries: list[ParticipantRiskSummary],
    cancer_labels: dict[str, dict[str, Any]],
) -> list[GxEAnalysisResult]:
    """Compare observed ORs from simulation against published references."""
    return analyze_gxe_interactions(summaries, cancer_labels)


def compute_roc_auc(
    summaries: list[ParticipantRiskSummary],
    cancer_labels: dict[str, dict[str, Any]],
    *,
    score_key: str = "interaction_factor",
    cancer_type: str = "any_cancer",
) -> ValidationResult:
    """Compute ROC/AUC for a risk score vs cancer outcome."""
    pairs: list[tuple[float, int]] = []
    for summary in summaries:
        score = _score_value(summary, score_key)
        if score is None:
            continue
        outcome = 1 if _label_has_cancer(cancer_labels, summary.person_id, cancer_type) else 0
        pairs.append((score, outcome))

    n_cases = sum(label for _, label in pairs)
    n_controls = len(pairs) - n_cases
    if not pairs or n_cases == 0 or n_controls == 0:
        return ValidationResult(
            roc_auc=0.5,
            roc_auc_se=0.0,
            roc_auc_95ci=(0.5, 0.5),
            n_cases=n_cases,
            n_controls=n_controls,
            calibration_bins=[],
            hosmer_lemeshow_chi2=None,
        )

    ranked = sorted(pairs, key=lambda item: item[0], reverse=True)
    roc_points: list[tuple[float, float]] = [(0.0, 0.0)]
    true_positive = 0
    false_positive = 0
    previous_score: float | None = None
    for score, label in ranked:
        if previous_score is not None and score != previous_score:
            roc_points.append((false_positive / n_controls, true_positive / n_cases))
        if label == 1:
            true_positive += 1
        else:
            false_positive += 1
        previous_score = score
    roc_points.append((false_positive / n_controls, true_positive / n_cases))

    auc = 0.0
    for index in range(1, len(roc_points)):
        x0, y0 = roc_points[index - 1]
        x1, y1 = roc_points[index]
        auc += (x1 - x0) * (y0 + y1) / 2

    q1 = auc / (2 - auc) if auc != 2 else 0.0
    q2 = (2 * auc * auc) / (1 + auc) if auc != -1 else 0.0
    se = math.sqrt(
        max(
            0.0,
            (auc * (1 - auc) + (n_cases - 1) * (q1 - auc * auc) + (n_controls - 1) * (q2 - auc * auc))
            / (n_cases * n_controls),
        )
    )

    raw_scores = [score for score, _ in pairs]
    max_score = max(raw_scores)
    min_score = min(raw_scores)
    normalised: list[tuple[float, int]] = []
    for score, label in sorted(pairs, key=lambda item: item[0]):
        if max_score == min_score:
            normalised_score = 0.5
        else:
            normalised_score = (score - min_score) / (max_score - min_score)
        normalised.append((normalised_score, label))

    bin_size = max(1, len(normalised) // min(10, max(1, len(normalised))))
    calibration_bins: list[dict[str, float]] = []
    for index in range(0, len(normalised), bin_size):
        bin_batch = normalised[index : index + bin_size]
        predicted_mean = sum(score for score, _ in bin_batch) / len(bin_batch)
        observed_rate = sum(label for _, label in bin_batch) / len(bin_batch)
        calibration_bins.append(
            {
                "predicted_mean": round(predicted_mean, 4),
                "observed_rate": round(observed_rate, 4),
                "n": float(len(bin_batch)),
            }
        )

    hosmer_lemeshow = 0.0
    for calibration_bin in calibration_bins:
        expected = calibration_bin["predicted_mean"] * calibration_bin["n"]
        observed = calibration_bin["observed_rate"] * calibration_bin["n"]
        if expected > 0:
            hosmer_lemeshow += (observed - expected) ** 2 / expected

    return ValidationResult(
        roc_auc=round(auc, 4),
        roc_auc_se=round(se, 4),
        roc_auc_95ci=(round(max(0.0, auc - 1.96 * se), 4), round(min(1.0, auc + 1.96 * se), 4)),
        n_cases=n_cases,
        n_controls=n_controls,
        calibration_bins=calibration_bins,
        hosmer_lemeshow_chi2=round(hosmer_lemeshow, 4),
    )


def get_cancer_phenotype_definitions() -> dict[str, dict[str, Any]]:
    """Return all cancer phenotype definitions."""
    return {name: dict(definition) for name, definition in CANCER_PHENOTYPES.items()}


def run_population_simulation(
    participants: list[SyntheticParticipant] | None = None,
    *,
    n_synthetic: int = 1000,
    seed: int = 42,
    cancer_labels: dict[str, dict[str, Any]] | None = None,
    validate: bool = True,
) -> PopulationSimulationResult:
    """Run population-scale simulation."""
    start = time.monotonic()
    if participants is None:
        participants = generate_synthetic_cohort(n=n_synthetic, seed=seed)

    summaries = [run_participant_risk(participant) for participant in participants]
    n_errors = sum(1 for summary in summaries if summary.error)
    population_stats = compute_population_stats(summaries)

    labels = cancer_labels
    if validate and labels is None:
        labels = generate_synthetic_cancer_labels(participants, seed=seed)

    gxe_results = analyze_gxe_interactions(summaries, labels if validate else None)
    validation: dict[str, ValidationResult] | None = None
    or_comparison: list[GxEAnalysisResult] | None = None

    if validate and labels is not None:
        validation = {
            "interaction_factor": compute_roc_auc(summaries, labels, score_key="interaction_factor", cancer_type="any_cancer"),
            "total_interaction_risk": compute_roc_auc(summaries, labels, score_key="total_interaction_risk", cancer_type="any_cancer"),
            "high_risk_pathway_count": compute_roc_auc(summaries, labels, score_key="high_risk_pathway_count", cancer_type="any_cancer"),
        }
        or_comparison = validate_against_published_ors(summaries, labels)

    return PopulationSimulationResult(
        n_participants=len(participants),
        n_completed=len(summaries) - n_errors,
        n_errors=n_errors,
        participant_summaries=summaries,
        population_stats=population_stats,
        gxe_analyses=gxe_results,
        validation=validation,
        or_comparison=or_comparison,
        elapsed_seconds=round(time.monotonic() - start, 4),
    )
