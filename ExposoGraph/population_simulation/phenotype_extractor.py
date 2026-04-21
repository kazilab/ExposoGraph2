#!/usr/bin/env python3
"""
Cancer Phenotype Extractor for All of Us EHR Data
====================================================
Extracts cancer diagnoses from All of Us OMOP CDM condition_occurrence
table using validated ICD-10-CM and SNOMED concept mappings.

Provides outcome labels for validating ExposoGraph risk predictions
against observed cancer incidence.

Authors: Kenneth J. Pienta (JHU), Julhash U. Kazi (Lund University)
Date: March 2026
"""

from typing import Any, TypeAlias, cast

JsonDict: TypeAlias = dict[str, Any]
ConditionRecord: TypeAlias = dict[str, Any]
ConditionRecords: TypeAlias = list[ConditionRecord]
CancerLabel: TypeAlias = dict[str, Any]
CancerLabels: TypeAlias = dict[str, CancerLabel]
ParticipantRecord: TypeAlias = dict[str, Any]
ParticipantRecords: TypeAlias = list[ParticipantRecord]


# ══════════════════════════════════════════════════════════
#  CANCER PHENOTYPE DEFINITIONS
# ══════════════════════════════════════════════════════════

# SNOMED ancestor concept for all malignant neoplasms
SNOMED_MALIGNANT_NEOPLASM = 443392  # 363346000 "Malignant neoplastic disease"

# ICD-10-CM cancer codes organized by organ system
# These map to ExposoGraph tissue targets
CANCER_PHENOTYPES: dict[str, JsonDict] = {
    "lung_cancer": {
        "icd10_codes": ["C34", "C34.0", "C34.1", "C34.2", "C34.3",
                        "C34.8", "C34.9", "C33"],
        "snomed_concept_ids": [
            253717, 254637, 4110588, 4110589,  # Malignant neoplasm of lung
            258375,  # Malignant neoplasm of trachea
        ],
        "tissue": "Lung",
        "carcinogen_classes": ["PAH", "Benzene", "Dioxin", "Nitrosamine",
                               "HeavyMetal_Chromium", "HeavyMetal_Nickel",
                               "HeavyMetal_Arsenic"],
        "description": "Primary lung cancer (bronchus and lung)",
    },
    "liver_cancer": {
        "icd10_codes": ["C22", "C22.0", "C22.1", "C22.2", "C22.3",
                        "C22.4", "C22.7", "C22.8", "C22.9"],
        "snomed_concept_ids": [
            4245614, 197029, 4091513,  # HCC, cholangiocarcinoma
        ],
        "tissue": "Liver",
        "carcinogen_classes": ["Aflatoxin", "Aldehyde", "ChlorinatedSolvent",
                               "HeavyMetal_Arsenic"],
        "description": "Primary liver cancer (HCC and cholangiocarcinoma)",
    },
    "bladder_cancer": {
        "icd10_codes": ["C67", "C67.0", "C67.1", "C67.2", "C67.3",
                        "C67.4", "C67.5", "C67.6", "C67.7", "C67.8",
                        "C67.9"],
        "snomed_concept_ids": [
            197508, 4181483,  # Malignant neoplasm of bladder
        ],
        "tissue": "Bladder",
        "carcinogen_classes": ["PAH", "HCA", "HeavyMetal_Arsenic",
                               "HeavyMetal_Cadmium"],
        "description": "Primary bladder cancer",
    },
    "kidney_cancer": {
        "icd10_codes": ["C64", "C64.1", "C64.2", "C64.9", "C65"],
        "snomed_concept_ids": [
            197320, 4091344,  # Renal cell carcinoma
        ],
        "tissue": "Kidney",
        "carcinogen_classes": ["ChlorinatedSolvent_TCE",
                               "HeavyMetal_Cadmium", "HeavyMetal_Arsenic"],
        "description": "Primary kidney cancer (RCC)",
    },
    "colorectal_cancer": {
        "icd10_codes": ["C18", "C19", "C20", "C18.0", "C18.1", "C18.2",
                        "C18.3", "C18.4", "C18.5", "C18.6", "C18.7",
                        "C18.8", "C18.9"],
        "snomed_concept_ids": [
            432837, 438981, 4181351,  # Colorectal cancer
        ],
        "tissue": "Colon",
        "carcinogen_classes": ["HCA", "PAH", "NDMA", "Nitrosamine"],
        "description": "Primary colorectal cancer",
    },
    "esophageal_cancer": {
        "icd10_codes": ["C15", "C15.3", "C15.4", "C15.5", "C15.8",
                        "C15.9"],
        "snomed_concept_ids": [
            24576, 4109054,  # Malignant neoplasm of esophagus
        ],
        "tissue": "Esophagus",
        "carcinogen_classes": ["Aldehyde_Acetaldehyde", "Nitrosamine",
                               "NDMA", "PAH"],
        "description": "Primary esophageal cancer (SCC and adenocarcinoma)",
    },
    "breast_cancer": {
        "icd10_codes": ["C50", "C50.0", "C50.1", "C50.2", "C50.3",
                        "C50.4", "C50.5", "C50.6", "C50.8", "C50.9"],
        "snomed_concept_ids": [
            137809, 4112853, 4048384,  # Breast cancer
        ],
        "tissue": "Breast",
        "carcinogen_classes": ["PAH", "Dioxin", "HCA"],
        "description": "Primary breast cancer",
    },
    "prostate_cancer": {
        "icd10_codes": ["C61"],
        "snomed_concept_ids": [
            200962, 4163261,  # Malignant neoplasm of prostate
        ],
        "tissue": "Prostate",
        "carcinogen_classes": ["PAH", "HCA", "HeavyMetal_Cadmium",
                               "HeavyMetal_Arsenic"],
        "description": "Primary prostate cancer",
    },
    "head_neck_cancer": {
        "icd10_codes": ["C01", "C02", "C03", "C04", "C05", "C06",
                        "C09", "C10", "C11", "C12", "C13", "C14",
                        "C30", "C31", "C32"],
        "snomed_concept_ids": [
            4119138, 4113900,  # Head and neck SCC
        ],
        "tissue": "Esophagus",  # closest available ExposoGraph tissue
        "carcinogen_classes": ["PAH", "Aldehyde_Acetaldehyde",
                               "Nitrosamine", "HeavyMetal_Nickel"],
        "description": "Head and neck squamous cell carcinoma",
    },
    "leukemia": {
        "icd10_codes": ["C91", "C92", "C93", "C94", "C95"],
        "snomed_concept_ids": [
            440383, 434298,  # Leukemia
        ],
        "tissue": "Liver",  # systemic; use liver as metabolic target
        "carcinogen_classes": ["Benzene", "Dioxin"],
        "description": "Leukemia (all types)",
    },
    "any_cancer": {
        "icd10_prefix": "C",
        "snomed_ancestor_id": SNOMED_MALIGNANT_NEOPLASM,
        "tissue": None,
        "carcinogen_classes": [],
        "description": "Any malignant neoplasm (umbrella phenotype)",
    },
}

# Phenotypes to EXCLUDE (benign, uncertain, in-situ)
EXCLUSION_ICD10_PREFIXES = ["D0", "D1", "D2", "D3", "D4"]


# ══════════════════════════════════════════════════════════
#  BIGQUERY SQL FOR PHENOTYPE EXTRACTION
# ══════════════════════════════════════════════════════════

PHENOTYPE_SQL: dict[str, str] = {
    "cancer_conditions": """
    -- Extract cancer diagnoses from condition_occurrence
    -- Uses SNOMED standard concepts mapped from ICD-10-CM
    -- Run in All of Us Researcher Workbench (Controlled Tier)
    SELECT
        co.person_id,
        co.condition_concept_id,
        c.concept_name AS condition_name,
        c.vocabulary_id,
        c.concept_code,
        co.condition_start_date,
        co.condition_source_value,
        cs.concept_name AS source_concept_name,
        co.condition_status_concept_id,
        cst.concept_name AS condition_status
    FROM `{CDR}.condition_occurrence` co
    JOIN `{CDR}.concept` c
        ON co.condition_concept_id = c.concept_id
    LEFT JOIN `{CDR}.concept` cs
        ON co.condition_source_concept_id = cs.concept_id
    LEFT JOIN `{CDR}.concept` cst
        ON co.condition_status_concept_id = cst.concept_id
    WHERE co.condition_concept_id IN (
        -- Get all descendants of "Malignant neoplastic disease"
        SELECT descendant_concept_id
        FROM `{CDR}.concept_ancestor`
        WHERE ancestor_concept_id = 443392
    )
    -- Only participants with WGS data
    AND co.person_id IN (
        SELECT DISTINCT person_id
        FROM `{CDR}.cb_search_person`
        WHERE has_whole_genome_variant = 1
    )
    ORDER BY co.person_id, co.condition_start_date
    """,

    "cancer_by_icd10": """
    -- Alternative: Extract using ICD-10 source codes directly
    SELECT
        co.person_id,
        co.condition_source_value AS icd10_code,
        cs.concept_name AS source_diagnosis,
        co.condition_start_date,
        co.condition_concept_id,
        c.concept_name AS standard_name
    FROM `{CDR}.condition_occurrence` co
    JOIN `{CDR}.concept` cs
        ON co.condition_source_concept_id = cs.concept_id
    JOIN `{CDR}.concept` c
        ON co.condition_concept_id = c.concept_id
    WHERE co.condition_source_value LIKE 'C%'
    AND co.condition_source_value NOT LIKE 'C44%'  -- exclude non-melanoma skin
    AND co.person_id IN (
        SELECT DISTINCT person_id
        FROM `{CDR}.cb_search_person`
        WHERE has_whole_genome_variant = 1
    )
    ORDER BY co.person_id, co.condition_start_date
    """,

    "cancer_counts_by_type": """
    -- Summary: cancer case counts by type
    SELECT
        c.concept_name,
        COUNT(DISTINCT co.person_id) AS n_patients,
        MIN(co.condition_start_date) AS earliest,
        MAX(co.condition_start_date) AS latest
    FROM `{CDR}.condition_occurrence` co
    JOIN `{CDR}.concept` c
        ON co.condition_concept_id = c.concept_id
    WHERE co.condition_concept_id IN (
        SELECT descendant_concept_id
        FROM `{CDR}.concept_ancestor`
        WHERE ancestor_concept_id = 443392
    )
    GROUP BY c.concept_name
    HAVING COUNT(DISTINCT co.person_id) >= 20
    ORDER BY n_patients DESC
    """,
}


# ══════════════════════════════════════════════════════════
#  PHENOTYPE CLASSIFICATION FUNCTIONS
# ══════════════════════════════════════════════════════════

def classify_cancer_from_icd10(icd10_code: str) -> list[str]:
    """
    Classify an ICD-10-CM code into ExposoGraph cancer phenotypes.

    Parameters
    ----------
    icd10_code : str
        ICD-10-CM code (e.g., "C34.1", "C22.0")

    Returns
    -------
    list : matching phenotype names (can match multiple)
    """
    code = icd10_code.strip().upper()
    matches: list[str] = []

    for pheno_name, pheno_def in CANCER_PHENOTYPES.items():
        if pheno_name == "any_cancer":
            if code.startswith("C") and not any(
                code.startswith(ex) for ex in EXCLUSION_ICD10_PREFIXES
            ):
                matches.append(pheno_name)
            continue

        icd_codes = cast(list[str], pheno_def.get("icd10_codes", []))
        for ref_code in icd_codes:
            if code.startswith(ref_code):
                matches.append(pheno_name)
                break

    return matches


def classify_cancer_from_snomed(concept_id: int) -> list[str]:
    """
    Classify a SNOMED concept_id into ExposoGraph cancer phenotypes.

    Parameters
    ----------
    concept_id : int
        OMOP standard concept_id

    Returns
    -------
    list : matching phenotype names
    """
    matches: list[str] = []
    for pheno_name, pheno_def in CANCER_PHENOTYPES.items():
        if pheno_name == "any_cancer":
            continue
        snomed_ids = cast(list[int], pheno_def.get("snomed_concept_ids", []))
        if concept_id in snomed_ids:
            matches.append(pheno_name)
    return matches


def build_cancer_labels(condition_records: ConditionRecords) -> CancerLabels:
    """
    Build binary cancer outcome labels from a list of condition records.

    Parameters
    ----------
    condition_records : list of dict
        Each dict has: person_id, condition_source_value (ICD-10),
        condition_concept_id (SNOMED), condition_start_date

    Returns
    -------
    dict
        Mapping of ``person_id`` to binary cancer labels, matched cancer
        types, per-diagnosis details, earliest cancer date, and diagnosis
        counts.
    """
    labels: CancerLabels = {}

    for record in condition_records:
        pid = str(record.get("person_id", ""))
        if not pid:
            continue

        icd10 = str(record.get("condition_source_value", "")).strip()
        snomed_id = record.get("condition_concept_id")
        date = str(record.get("condition_start_date", ""))

        # Classify
        phenotypes = set()
        if icd10:
            phenotypes.update(classify_cancer_from_icd10(icd10))
        if snomed_id:
            phenotypes.update(classify_cancer_from_snomed(int(snomed_id)))

        if not phenotypes:
            continue

        if pid not in labels:
            labels[pid] = {
                "any_cancer": False,
                "cancer_types": set(),
                "cancer_details": [],
                "earliest_cancer_date": None,
                "n_cancer_diagnoses": 0,
            }

        label_record = labels[pid]
        label_record["any_cancer"] = True
        cast(set[str], label_record["cancer_types"]).update(phenotypes)
        cast(list[JsonDict], label_record["cancer_details"]).append({
            "type": list(phenotypes),
            "icd10": icd10,
            "date": date,
        })
        label_record["n_cancer_diagnoses"] = int(label_record["n_cancer_diagnoses"]) + 1

        earliest = cast(str | None, label_record["earliest_cancer_date"])
        if date and (earliest is None or date < earliest):
            label_record["earliest_cancer_date"] = date

    # Convert sets to lists for JSON serialization
    for pid in labels:
        labels[pid]["cancer_types"] = sorted(cast(set[str], labels[pid]["cancer_types"]))

    return labels


def generate_synthetic_cancer_labels(
    participants: ParticipantRecords,
    base_rate: float = 0.10,
    genotype_effects: bool = True,
    seed: int = 42
) -> CancerLabels:
    """
    Generate synthetic cancer outcome labels for testing.

    Uses published genotype-cancer ORs to create realistic case/control ratios:
    - GSTM1-null + smoking → lung cancer OR ~1.5-2.0
    - ALDH2*2 + alcohol → esophageal cancer OR ~6-12
    - NAT2 slow + occupational → bladder cancer OR ~1.5-2.0

    Parameters
    ----------
    participants : list
        From allofus_adapter.generate_synthetic_cohort()
    base_rate : float
        Background cancer incidence rate
    genotype_effects : bool
        If True, modulate rates by genotype/exposure
    seed : int
        Random seed

    Returns
    -------
    dict : {person_id: cancer_label_dict}
    """
    from ._core import SyntheticParticipant
    from ._core import generate_synthetic_cancer_labels as _generate_core_synthetic_cancer_labels
    from .allofus_adapter import _derive_exposure_scenario

    normalized_participants: list[SyntheticParticipant] = []
    for participant in participants:
        if isinstance(participant, SyntheticParticipant):
            normalized_participants.append(participant)
            continue
        lifestyle = dict(participant.get("lifestyle", {}))
        occupational = bool(
            lifestyle.get("occupational_exposure", lifestyle.get("occupational_risk", False))
        )
        ancestry_value = str(
            participant.get("ancestry_label", participant.get("ancestry", "Unknown"))
        )
        normalized_participants.append(
            SyntheticParticipant(
                person_id=str(participant.get("person_id")),
                genotypes=dict(participant.get("genotypes", {})),
                ancestry=ancestry_value,
                ancestry_label=ancestry_value,
                lifestyle={
                    "smoking": bool(lifestyle.get("smoking", False)),
                    "alcohol_heavy": bool(lifestyle.get("alcohol_heavy", False)),
                    "alcohol_moderate": bool(lifestyle.get("alcohol_moderate", False)),
                    "occupational_exposure": occupational,
                },
                exposure_scenario=str(
                    participant.get(
                        "exposure_scenario",
                        _derive_exposure_scenario(
                            {
                                "smoking": bool(lifestyle.get("smoking", False)),
                                "alcohol_heavy": bool(lifestyle.get("alcohol_heavy", False)),
                                "alcohol_moderate": bool(lifestyle.get("alcohol_moderate", False)),
                                "occupational_risk": occupational,
                                "occupational_exposure": occupational,
                            }
                        ),
                    )
                ),
                tissue=str(participant.get("tissue", "Liver")),
            )
        )

    labels: CancerLabels = _generate_core_synthetic_cancer_labels(
        normalized_participants,
        base_rate=base_rate,
        genotype_effects=genotype_effects,
        seed=seed,
    )

    for label in labels.values():
        cancers = cast(list[str], label.get("cancer_types", []))
        label["cancer_details"] = [
            {"type": [cancer_type], "icd10": "SYN", "date": "2025-01-01"}
            for cancer_type in cancers
        ]
        label["earliest_cancer_date"] = "2025-01-01" if cancers else None

    n_cases = sum(1 for value in labels.values() if value["any_cancer"])
    print(f"  Generated cancer labels: {n_cases}/{len(labels)} cases "
          f"({n_cases/len(labels)*100:.1f}%)")

    type_counts: dict[str, int] = {}
    for value in labels.values():
        for cancer_type in cast(list[str], value["cancer_types"]):
            type_counts[cancer_type] = type_counts.get(cancer_type, 0) + 1
    for cancer_type, count in sorted(type_counts.items(), key=lambda item: -item[1]):
        print(f"    {cancer_type}: {count}")

    return labels


def get_phenotype_sql(query_name: str, cdr_dataset: str | None = None) -> str:
    """Get a phenotype extraction SQL template."""
    sql = PHENOTYPE_SQL.get(query_name)
    if sql is None:
        raise ValueError(f"Unknown query: {query_name}. "
                         f"Available: {list(PHENOTYPE_SQL.keys())}")
    if cdr_dataset:
        sql = sql.replace("{CDR}", cdr_dataset)
    return sql


def get_cancer_phenotype_definitions() -> dict[str, JsonDict]:
    """Return all cancer phenotype definitions."""
    return CANCER_PHENOTYPES


# ══════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Cancer Phenotype Extractor for All of Us"
    )
    parser.add_argument("--list-phenotypes", action="store_true",
                        help="List all cancer phenotype definitions")
    parser.add_argument("--classify-icd10", type=str,
                        help="Classify an ICD-10 code")
    parser.add_argument("--show-sql", type=str,
                        help="Print SQL template (cancer_conditions, "
                             "cancer_by_icd10, cancer_counts_by_type)")
    parser.add_argument("--generate-synthetic", type=int, metavar="N",
                        help="Generate synthetic labels for N participants")
    args = parser.parse_args()

    if args.list_phenotypes:
        print("\nCancer Phenotype Definitions:")
        print("=" * 60)
        for name, defn in CANCER_PHENOTYPES.items():
            print(f"\n  {name}")
            print(f"    Tissue: {defn.get('tissue', 'N/A')}")
            print(f"    ICD-10: {defn.get('icd10_codes', defn.get('icd10_prefix', 'N/A'))}")
            print(f"    Carcinogens: {defn.get('carcinogen_classes', [])}")
            print(f"    {defn.get('description', '')}")
        return

    if args.classify_icd10:
        matches = classify_cancer_from_icd10(args.classify_icd10)
        print(f"ICD-10 '{args.classify_icd10}' → {matches}")
        return

    if args.show_sql:
        print(get_phenotype_sql(args.show_sql, "YOUR_CDR_DATASET"))
        return

    if args.generate_synthetic:
        from .allofus_adapter import generate_synthetic_cohort
        participants = generate_synthetic_cohort(args.generate_synthetic)
        labels = generate_synthetic_cancer_labels(participants)
        print(
            f"Generated {len(participants)} synthetic participants with "
            f"{len(labels)} cancer label records."
        )
        return

    parser.print_help()


if __name__ == "__main__":
    main()
