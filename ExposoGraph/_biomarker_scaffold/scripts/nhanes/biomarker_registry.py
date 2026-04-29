"""NHANES variable-to-biomarker mapping registry.

Keep this file limited to measurement metadata. Model parameters stay in
biomarker_mapping.json or the registry resolver layer.
"""

BIOMARKER_REGISTRY = {
    "PAH": {
        "URXP10": {
            "biomarker": "urinary_1_hydroxypyrene",
            "parent_compound": "pyrene",
            "mw_g_mol": 218.25,
            "units_expected": "ng_per_L",
        },
    },
    "COTININE": {
        "LBXCOT": {
            "biomarker": "serum_cotinine",
            "parent_compound": "nicotine",
            "units_expected": "ng_per_mL",
        },
    },
    "VOC_BLOOD": {
        "LBXVBZ": {
            "biomarker": "blood_benzene",
            "parent_compound": "benzene",
            "mw_g_mol": 78.11,
            "units_expected": "ng_per_mL",
        },
    },
    "METALS_URINE": {
        "URXUAS": {
            "biomarker": "urinary_total_arsenic",
            "parent_compound": "arsenic",
            "units_expected": "ug_per_L",
        },
    },
}
