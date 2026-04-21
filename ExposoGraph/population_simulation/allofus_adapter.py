#!/usr/bin/env python3
"""
All of Us → ExposoGraph Adapter
=================================
Maps All of Us Researcher Workbench CDR tables into ExposoGraph-compatible
genotype, lifestyle, and exposure input formats.

Data sources consumed:
    1. PGx auxiliary dataset (per-gene TSVs) → ExposoGraph genotype dict
    2. Survey / observation tables (OMOP) → lifestyle + exposure scenario
    3. Ancestry predictions (ancestry_preds.tsv) → ancestry label

Designed to run INSIDE the All of Us Researcher Workbench Jupyter environment,
where CDR tables are accessible via BigQuery SQL or pandas.

Authors: Kenneth J. Pienta (JHU), Julhash U. Kazi (Lund University)
Date: March 2026
"""

import csv
import json
from pathlib import Path
from typing import Any, Mapping, TypeAlias, cast

JsonDict: TypeAlias = dict[str, Any]
PGXGeneRecord: TypeAlias = dict[str, str]
PGXParticipantRecord: TypeAlias = dict[str, PGXGeneRecord]
PGXDataset: TypeAlias = dict[str, PGXParticipantRecord]
AncestryRecord: TypeAlias = dict[str, Any]
AncestryDataset: TypeAlias = dict[str, AncestryRecord]
LifestyleRecord: TypeAlias = dict[str, Any]
WGSRecord: TypeAlias = dict[str, Any]
WGSDataset: TypeAlias = dict[str, WGSRecord]
ParticipantRecord: TypeAlias = dict[str, Any]
ParticipantRecords: TypeAlias = list[ParticipantRecord]

# ══════════════════════════════════════════════════════════
#  CONSTANTS
# ══════════════════════════════════════════════════════════

# All of Us PGx genes and their mapping to ExposoGraph gene names
ALLOFUS_PGX_GENES: dict[str, str] = {
    "CYP2D6":       "CYP2D6",
    "CYP3A5":       "CYP3A5",
    "CYP2B6":       "CYP2B6",
    "CYP2C9":       "CYP2C9",
    "CYP2C19":      "CYP2C19",
    "UGT1A1":       "UGT1A1",
    "DPYD":         "DPYD",
    "TPMT":         "TPMT",
    "SLCO1B1":      "SLCO1B1",
    "G6PD":         "G6PD",
    "ABCG2":        "ABCG2",
    "CYP4F2":       "CYP4F2",
    "VKORC1":       "VKORC1",
    "NUDT15":       "NUDT15",
    "CYP2C_CLUSTER": "CYP2C_CLUSTER",
    "CACNA1S":      "CACNA1S",
    "CFTR":         "CFTR",
    "RYR1":         "RYR1",
}

# Genes directly used by ExposoGraph flux/interaction engines
EXPOSOGRAPH_RELEVANT_GENES = {
    "CYP2D6", "CYP3A5", "CYP2B6", "CYP2C9", "CYP2C19",
    "UGT1A1", "CYP4F2",
}

# Genes NOT in All of Us PGx but needed by ExposoGraph — must infer or default
EXPOSOGRAPH_EXTRA_GENES: dict[str, str] = {
    "CYP1A1": "NM",       # Not in AoU PGx panel; default to NM
    "CYP1A2": "NM",       # Not in AoU PGx panel
    "CYP2E1": "NM",       # Not in AoU PGx panel
    "CYP3A4": "NM",       # Not in AoU PGx panel
    "GSTM1":  "present",  # Not in AoU PGx panel; can be derived from WGS
    "GSTP1":  "NM",       # Not in AoU PGx panel
    "GSTT1":  "present",  # Not in AoU PGx panel; can be derived from WGS
    "NAT2":   "rapid",    # Not in AoU PGx panel
    "ALDH2":  "*1/*1",    # Not in AoU PGx panel
}

# Phenotype → ExposoGraph activity mapping
PHENOTYPE_ACTIVITY_MAP: dict[str, str] = {
    # Metabolizer phenotypes
    "ultrarapid metabolizer":   "UM",
    "rapid metabolizer":        "RM",
    "normal metabolizer":       "NM",
    "intermediate metabolizer": "IM",
    "poor metabolizer":         "PM",
    "indeterminate":            "NM",  # default to normal if unknown
    # Function phenotypes
    "normal function":          "NM",
    "decreased function":       "IM",
    "poor function":            "PM",
    "possible decreased function": "IM",
    "possible poor function":   "PM",
    # CYP2C cluster
    "variant present":          "variant",
    "variant absent":           "NM",
    # G6PD
    "normal":                   "NM",
    "variable":                 "IM",
    "deficient":                "PM",
}

# All of Us survey concept IDs for lifestyle factors (OMOP observation table)
# These are the LOINC/PPI concept_ids used in All of Us surveys
SURVEY_CONCEPTS: dict[str, JsonDict] = {
    "smoking_status": {
        # "Tobacco smoking status" LOINC 72166-2 → OMOP concept_id 40766929
        "concept_ids": [40766929, 1585860, 1585864, 1585870],
        "positive_values": [
            "Current every day smoker",
            "Current some day smoker",
            "Yes",
        ],
        "negative_values": [
            "Never smoker",
            "Former smoker",
            "No",
        ],
    },
    "alcohol_use": {
        # "How often do you have a drink containing alcohol" LOINC 68518-0
        "concept_ids": [40771103, 40771104],
        "heavy_values": [
            "4 or more times a week",
            "Daily or almost daily",
        ],
        "moderate_values": [
            "2-3 times a week",
            "2-4 times a month",
        ],
    },
    "occupational_exposure": {
        # Occupation-related survey fields
        "concept_ids": [40771090, 1585952],
        "high_risk_values": [
            "Construction",
            "Manufacturing",
            "Mining",
            "Agriculture",
        ],
    },
}

# Ancestry label mapping
ANCESTRY_MAP: dict[str, str] = {
    "afr": "African",
    "amr": "Admixed American",
    "eas": "East Asian",
    "eur": "European",
    "mid": "Middle Eastern",
    "sas": "South Asian",
    "oth": "Other/Admixed",
}


# ══════════════════════════════════════════════════════════
#  PGx GENOTYPE ADAPTER
# ══════════════════════════════════════════════════════════

def parse_pgx_tsv(filepath: str) -> dict[str, PGXGeneRecord]:
    """
    Parse a single All of Us PGx TSV file into a dict of
    {person_id: {gene, genotype, phenotype, activity_score, ...}}.

    Parameters
    ----------
    filepath : str
        Path to a per-gene PGx TSV file from the AoU Controlled CDR.

    Returns
    -------
    dict : {person_id: row_dict}
    """
    results: dict[str, PGXGeneRecord] = {}
    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = row.get("person_id", "").strip()
            if pid:
                results[pid] = {
                    "gene": row.get("gene", "").strip(),
                    "genotype": row.get("genotype", "").strip(),
                    "phenotype": row.get("phenotype", "").strip(),
                    "activity_score": row.get("activity_score", ""),
                    "copy_number": row.get("copy_number", ""),
                    "cyrius_filter": row.get("cyrius_filter", ""),
                }
    return results


def load_all_pgx_files(pgx_directory: str) -> PGXDataset:
    """
    Load ALL per-gene PGx TSV files from the AoU PGx auxiliary directory.

    Parameters
    ----------
    pgx_directory : str
        Path to directory containing per-gene TSV files.
        In AoU Workbench: access via Controlled CDR Directory paths.

    Returns
    -------
    dict : {person_id: {gene_name: {genotype, phenotype, ...}}}
    """
    all_data: PGXDataset = {}
    pgx_dir = Path(pgx_directory)

    for tsv_file in sorted(pgx_dir.glob("*.tsv")):
        gene_data = parse_pgx_tsv(str(tsv_file))
        for pid, data in gene_data.items():
            if pid not in all_data:
                all_data[pid] = {}
            gene = data["gene"]
            all_data[pid][gene] = data

    print(f"  Loaded PGx data for {len(all_data)} participants across "
          f"{len(list(pgx_dir.glob('*.tsv')))} gene files")
    return all_data


def pgx_to_exposograph_genotype(pgx_record: Mapping[str, Mapping[str, str]]) -> dict[str, str]:
    """
    Convert a single participant's All of Us PGx record into
    ExposoGraph genotype format.

    Parameters
    ----------
    pgx_record : dict
        {gene_name: {genotype, phenotype, activity_score, ...}}
        from load_all_pgx_files()[person_id]

    Returns
    -------
    dict
        ExposoGraph-compatible genotype mapping, for example
        ``{"CYP2D6": "*2/*4", "CYP3A5": "*3/*3", "GSTM1": "present"}``.
    """
    genotypes: dict[str, str] = {}

    for gene, data in pgx_record.items():
        # Map AoU gene name to ExposoGraph name
        eg_gene = ALLOFUS_PGX_GENES.get(gene, gene)

        diplotype = data.get("genotype", "")
        phenotype = data.get("phenotype", "").strip().lower()

        # For genes with star allele diplotypes (CYP family)
        if diplotype.startswith("*") and "/" in diplotype:
            genotypes[eg_gene] = diplotype
        # For genes with phenotype-based representation
        elif phenotype in PHENOTYPE_ACTIVITY_MAP:
            activity = PHENOTYPE_ACTIVITY_MAP[phenotype]
            genotypes[eg_gene] = activity
        # Fallback: store raw diplotype
        elif diplotype:
            genotypes[eg_gene] = diplotype
        else:
            genotypes[eg_gene] = "NM"

    # Add ExposoGraph-specific genes not in AoU panel with defaults
    for gene, default in EXPOSOGRAPH_EXTRA_GENES.items():
        if gene not in genotypes:
            genotypes[gene] = default

    return genotypes


def derive_gst_from_wgs(person_id: str, wgs_data: WGSDataset | None = None) -> dict[str, str]:
    """
    Derive GSTM1/GSTT1 deletion status from WGS copy number data.

    In All of Us, GSTM1 and GSTT1 are NOT in the PGx panel but can be
    derived from WGS structural variant calls or copy number estimation.

    Parameters
    ----------
    person_id : str
        All of Us research_id
    wgs_data : dict, optional
        Pre-loaded SV/CNV data. If None, returns defaults.

    Returns
    -------
    dict : {"GSTM1": "null"|"present", "GSTT1": "null"|"present"}
    """
    if wgs_data and person_id in wgs_data:
        sv = wgs_data[person_id]
        return {
            "GSTM1": "null" if sv.get("GSTM1_cn", 2) == 0 else "present",
            "GSTT1": "null" if sv.get("GSTT1_cn", 2) == 0 else "present",
        }

    # Default: cannot determine without SV data
    return {"GSTM1": "present", "GSTT1": "present"}


# ══════════════════════════════════════════════════════════
#  LIFESTYLE / EXPOSURE ADAPTER
# ══════════════════════════════════════════════════════════

def extract_lifestyle_from_surveys(
    person_id: str,
    observation_df: Any = None,
    survey_df: Mapping[str, str] | None = None,
) -> LifestyleRecord:
    """
    Extract lifestyle factors from All of Us survey/observation data.

    In the AoU Workbench, data comes from the OMOP observation table
    queried via BigQuery. The function accepts either a pandas DataFrame of
    observation records for one participant or a pre-extracted survey-answer
    dict.

    Parameters
    ----------
    person_id : str
    observation_df : pandas.DataFrame, optional
        Filtered observation table for this person_id. Expected columns are
        ``observation_concept_id``, ``value_as_concept_id``,
        ``value_source_value``, and ``value_as_string``.
    survey_df : dict, optional
        Pre-extracted survey responses such as
        ``{"smoking": "Current every day smoker", "alcohol": "4 or more times a week"}``.

    Returns
    -------
    dict
        Lifestyle mapping with smoking, alcohol, occupational exposure,
        and derived ``exposure_scenario`` fields in the package risk format.
    """
    del person_id

    lifestyle: LifestyleRecord = {
        "smoking": False,
        "alcohol_heavy": False,
        "alcohol_moderate": False,
        "occupational_risk": False,
        "occupational_exposure": False,
    }

    if survey_df and isinstance(survey_df, Mapping):
        # Direct survey answers
        smoking_val = survey_df.get("smoking", "").lower()
        lifestyle["smoking"] = any(
            v.lower() in smoking_val
            for v in ["current", "every day", "some day", "yes"]
        )

        alcohol_val = survey_df.get("alcohol", "").lower()
        lifestyle["alcohol_heavy"] = any(
            v.lower() in alcohol_val
            for v in ["4 or more", "daily", "almost daily"]
        )
        lifestyle["alcohol_moderate"] = any(
            v.lower() in alcohol_val
            for v in ["2-3 times", "2-4 times"]
        )

        occupation_val = survey_df.get("occupation", "").lower()
        lifestyle["occupational_risk"] = any(
            v.lower() in occupation_val
            for v in ["construction", "manufacturing", "mining", "agriculture"]
        )
        lifestyle["occupational_exposure"] = lifestyle["occupational_risk"]

    elif observation_df is not None:
        # Pandas DataFrame from OMOP observation table
        try:
            for _, row in observation_df.iterrows():
                cid = int(row.get("observation_concept_id", 0))
                val = str(row.get("value_source_value",
                          row.get("value_as_string", ""))).lower()

                # Smoking
                if cid in cast(list[int], SURVEY_CONCEPTS["smoking_status"]["concept_ids"]):
                    smoking_values = cast(
                        list[str], SURVEY_CONCEPTS["smoking_status"]["positive_values"]
                    )
                    lifestyle["smoking"] = any(v.lower() in val for v in smoking_values)

                # Alcohol
                if cid in cast(list[int], SURVEY_CONCEPTS["alcohol_use"]["concept_ids"]):
                    lifestyle["alcohol_heavy"] = any(
                        v.lower() in val
                        for v in cast(list[str], SURVEY_CONCEPTS["alcohol_use"]["heavy_values"])
                    )
                    lifestyle["alcohol_moderate"] = any(
                        v.lower() in val
                        for v in cast(list[str], SURVEY_CONCEPTS["alcohol_use"]["moderate_values"])
                    )

                # Occupation
                if cid in cast(list[int], SURVEY_CONCEPTS["occupational_exposure"]["concept_ids"]):
                    lifestyle["occupational_risk"] = any(
                        v.lower() in val
                        for v in cast(
                            list[str],
                            SURVEY_CONCEPTS["occupational_exposure"]["high_risk_values"],
                        )
                    )
                    lifestyle["occupational_exposure"] = lifestyle["occupational_risk"]
        except ImportError:
            pass  # pandas not available outside Workbench

    lifestyle["occupational_exposure"] = bool(
        lifestyle.get("occupational_exposure", lifestyle.get("occupational_risk", False))
    )
    # Derive exposure scenario from lifestyle
    lifestyle["exposure_scenario"] = _derive_exposure_scenario(lifestyle)
    return lifestyle


def _derive_exposure_scenario(lifestyle: Mapping[str, Any]) -> str:
    """
    Map lifestyle factors to the closest ExposoGraph exposure scenario.

    Returns one of the pre-defined scenario keys from the exposure engine.
    """
    smoking = lifestyle.get("smoking", False)
    heavy_alcohol = lifestyle.get("alcohol_heavy", False)
    occupational = bool(
        lifestyle.get("occupational_exposure", lifestyle.get("occupational_risk", False))
    )

    if smoking and heavy_alcohol and occupational:
        return "smoker_industrial_heavy_drinker"
    elif smoking and heavy_alcohol:
        return "smoker_heavy_drinker"
    elif smoking and occupational:
        return "smoker_industrial_worker"
    elif smoking:
        return "smoker"
    elif heavy_alcohol:
        return "heavy_drinker"
    elif occupational:
        return "industrial_worker"
    elif lifestyle.get("alcohol_moderate", False):
        return "moderate_drinker"
    else:
        return "general_population"


# ══════════════════════════════════════════════════════════
#  ANCESTRY ADAPTER
# ══════════════════════════════════════════════════════════

def load_ancestry_predictions(filepath: str) -> AncestryDataset:
    """
    Load All of Us genetic ancestry predictions.

    Parameters
    ----------
    filepath : str
        Path to ancestry_preds.tsv from AoU Controlled CDR.

    Returns
    -------
    dict : {person_id: {"ancestry": str, "ancestry_label": str,
                        "probabilities": list}}
    """
    results: AncestryDataset = {}
    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = row.get("research_id", "").strip()
            pred = row.get("ancestry_pred", "").strip().lower()
            pred_other = row.get("ancestry_pred_other", "").strip().lower()

            # Parse probability array
            prob_str = row.get("probabilities", "[]")
            try:
                probabilities = json.loads(prob_str.replace("'", '"'))
            except (json.JSONDecodeError, ValueError):
                probabilities = []

            # Use the "with other" prediction if it differs
            ancestry = pred_other if pred_other else pred
            results[pid] = {
                "ancestry": ancestry,
                "ancestry_label": ANCESTRY_MAP.get(ancestry, ancestry),
                "probabilities": probabilities,
            }

    print(f"  Loaded ancestry for {len(results)} participants")
    return results


# ══════════════════════════════════════════════════════════
#  UNIFIED PARTICIPANT BUILDER
# ══════════════════════════════════════════════════════════

def build_participant_record(
    person_id: str,
    pgx_data: PGXParticipantRecord,
    ancestry_data: AncestryDataset | None = None,
    lifestyle_data: LifestyleRecord | None = None,
    wgs_sv_data: WGSDataset | None = None,
) -> ParticipantRecord:
    """
    Assemble a complete ExposoGraph-ready participant record.

    Parameters
    ----------
    person_id : str
    pgx_data : dict
        {gene: {genotype, phenotype, ...}} from load_all_pgx_files
    ancestry_data : dict, optional
        Ancestry predictions dict from load_ancestry_predictions
    lifestyle_data : dict, optional
        Lifestyle dict from extract_lifestyle_from_surveys
    wgs_sv_data : dict, optional
        SV/CNV data for GST derivation

    Returns
    -------
    dict
        Complete participant record containing participant id, genotypes,
        lifestyle, exposure scenario, ancestry labels, and a default target
        tissue.
    """
    # 1. Genotypes
    genotypes = pgx_to_exposograph_genotype(pgx_data)

    # Override GSTM1/GSTT1 if SV data available
    gst = derive_gst_from_wgs(person_id, wgs_sv_data)
    genotypes.update(gst)

    # 2. Lifestyle
    if lifestyle_data is None:
        lifestyle_data = {
            "smoking": False,
            "alcohol_heavy": False,
            "alcohol_moderate": False,
            "occupational_risk": False,
            "occupational_exposure": False,
            "exposure_scenario": "general_population",
        }

    # 3. Ancestry
    ancestry_info: AncestryRecord = {}
    if ancestry_data and person_id in ancestry_data:
        ancestry_info = ancestry_data[person_id]

    return {
        "person_id": person_id,
        "genotypes": genotypes,
        "lifestyle": {
            "smoking": lifestyle_data.get("smoking", False),
            "alcohol_heavy": lifestyle_data.get("alcohol_heavy", False),
            "alcohol_moderate": lifestyle_data.get("alcohol_moderate", False),
            "occupational_risk": lifestyle_data.get("occupational_risk", False),
            "occupational_exposure": lifestyle_data.get(
                "occupational_exposure",
                lifestyle_data.get("occupational_risk", False),
            ),
        },
        "exposure_scenario": lifestyle_data.get(
            "exposure_scenario", "general_population"
        ),
        "ancestry": ancestry_info.get("ancestry", "unknown"),
        "ancestry_label": ancestry_info.get("ancestry_label", "Unknown"),
        "tissue": "Liver",  # default; override per analysis
    }


# ══════════════════════════════════════════════════════════
#  BIGQUERY SQL TEMPLATES
# ══════════════════════════════════════════════════════════

# These SQL templates run inside the All of Us Researcher Workbench
# to extract the data needed for population simulation.

BIGQUERY_SQL: dict[str, str] = {
    "smoking_status": """
    -- Extract smoking status for all participants
    -- Run in All of Us Researcher Workbench (Controlled Tier)
    SELECT
        o.person_id,
        c.concept_name AS smoking_status,
        o.observation_date
    FROM `{CDR}.observation` o
    JOIN `{CDR}.concept` c
        ON o.value_as_concept_id = c.concept_id
    WHERE o.observation_concept_id IN (
        40766929,  -- Tobacco smoking status LOINC 72166-2
        1585860,   -- Smoking status PPI
        1585864,   -- Cigarettes per day PPI
        1585870    -- Smoking frequency PPI
    )
    AND o.observation_date = (
        SELECT MAX(o2.observation_date)
        FROM `{CDR}.observation` o2
        WHERE o2.person_id = o.person_id
        AND o2.observation_concept_id = o.observation_concept_id
    )
    """,

    "alcohol_use": """
    -- Extract alcohol consumption for all participants
    SELECT
        o.person_id,
        c.concept_name AS alcohol_frequency,
        o.value_as_number AS drinks_per_occasion,
        o.observation_date
    FROM `{CDR}.observation` o
    JOIN `{CDR}.concept` c
        ON o.value_as_concept_id = c.concept_id
    WHERE o.observation_concept_id IN (
        40771103,  -- How often do you have a drink LOINC 68518-0
        40771104   -- How many standard drinks LOINC 68519-8
    )
    AND o.observation_date = (
        SELECT MAX(o2.observation_date)
        FROM `{CDR}.observation` o2
        WHERE o2.person_id = o.person_id
        AND o2.observation_concept_id = o.observation_concept_id
    )
    """,

    "occupation": """
    -- Extract occupational status for all participants
    SELECT
        o.person_id,
        c.concept_name AS occupation_status,
        o.value_source_value AS occupation_detail
    FROM `{CDR}.observation` o
    JOIN `{CDR}.concept` c
        ON o.value_as_concept_id = c.concept_id
    WHERE o.observation_concept_id IN (
        40771090,  -- Current occupational status LOINC 68505-7
        1585952    -- Occupation detail PPI
    )
    """,

    "demographics": """
    -- Extract demographics for all participants with WGS
    SELECT
        p.person_id,
        p.year_of_birth,
        p.gender_concept_id,
        gc.concept_name AS gender,
        p.race_concept_id,
        rc.concept_name AS race,
        p.ethnicity_concept_id,
        ec.concept_name AS ethnicity
    FROM `{CDR}.person` p
    LEFT JOIN `{CDR}.concept` gc ON p.gender_concept_id = gc.concept_id
    LEFT JOIN `{CDR}.concept` rc ON p.race_concept_id = rc.concept_id
    LEFT JOIN `{CDR}.concept` ec ON p.ethnicity_concept_id = ec.concept_id
    WHERE p.person_id IN (
        -- Only participants with WGS data
        SELECT DISTINCT person_id
        FROM `{CDR}.cb_search_person`
        WHERE has_whole_genome_variant = 1
    )
    """,
}


def get_bigquery_sql(query_name: str, cdr_dataset: str | None = None) -> str:
    """
    Get a BigQuery SQL template for All of Us data extraction.

    Parameters
    ----------
    query_name : str
        One of: smoking_status, alcohol_use, occupation, demographics
    cdr_dataset : str, optional
        CDR dataset name. In the Workbench, use:
        os.environ.get("WORKSPACE_CDR")

    Returns
    -------
    str : SQL query string
    """
    sql = BIGQUERY_SQL.get(query_name)
    if sql is None:
        raise ValueError(f"Unknown query: {query_name}. "
                         f"Available: {list(BIGQUERY_SQL.keys())}")

    if cdr_dataset:
        sql = sql.replace("{CDR}", cdr_dataset)
    return sql


# ══════════════════════════════════════════════════════════
#  SYNTHETIC DATA GENERATOR (for testing outside Workbench)
# ══════════════════════════════════════════════════════════

def generate_synthetic_cohort(
    n: int = 1000,
    seed: int = 42,
    output_dir: str | None = None
) -> ParticipantRecords:
    """
    Generate a synthetic All of Us-like cohort for testing the pipeline
    outside the Researcher Workbench.

    Uses published allele frequencies from PharmVar and CPIC to create
    realistic genotype distributions across ancestries.

    Parameters
    ----------
    n : int
        Number of synthetic participants
    seed : int
        Random seed for reproducibility
    output_dir : str, optional
        If provided, saves synthetic PGx TSVs and ancestry file

    Returns
    -------
    list : List of participant dicts ready for batch_runner
    """
    from ._core import generate_synthetic_cohort as _generate_core_synthetic_cohort

    label_to_code = {label: code for code, label in ANCESTRY_MAP.items()}
    typed_participants = _generate_core_synthetic_cohort(n=n, seed=seed)
    participants: ParticipantRecords = []
    for participant in typed_participants:
        ancestry_code = label_to_code.get(
            participant.ancestry_label, participant.ancestry_label.lower()
        )
        lifestyle = {
            "smoking": bool(participant.lifestyle.get("smoking", False)),
            "alcohol_heavy": bool(participant.lifestyle.get("alcohol_heavy", False)),
            "alcohol_moderate": bool(participant.lifestyle.get("alcohol_moderate", False)),
            "occupational_risk": bool(
                participant.lifestyle.get(
                    "occupational_risk",
                    participant.lifestyle.get("occupational_exposure", False),
                )
            ),
            "occupational_exposure": bool(
                participant.lifestyle.get("occupational_exposure", False)
            ),
        }
        participants.append(
            {
                "person_id": participant.person_id,
                "genotypes": dict(participant.genotypes),
                "lifestyle": lifestyle,
                "exposure_scenario": participant.exposure_scenario,
                "ancestry": ancestry_code,
                "ancestry_label": participant.ancestry_label,
                "tissue": participant.tissue,
            }
        )

    print(f"  Generated {n} synthetic participants")
    print(
        "  Ancestry distribution: "
        + ", ".join(
            f"{code}: {sum(1 for participant in participants if participant['ancestry'] == code)}"
            for code in sorted(ANCESTRY_MAP.keys())
        )
    )

    if output_dir:
        _save_synthetic_files(participants, output_dir)

    return participants


def _save_synthetic_files(participants: ParticipantRecords, output_dir: str) -> None:
    """Save synthetic participant data to files mimicking AoU format."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Save participant records
    with open(out / "synthetic_cohort.json", "w") as f:
        json.dump(participants, f, indent=2)

    # Save summary stats
    ancestry_counts: dict[str, int] = {}
    summary: JsonDict = {
        "n_participants": len(participants),
        "ancestry_counts": ancestry_counts,
        "smoking_count": sum(1 for p in participants if p["lifestyle"]["smoking"]),
        "heavy_alcohol_count": sum(1 for p in participants
                                    if p["lifestyle"]["alcohol_heavy"]),
    }
    for p in participants:
        ancestry = str(p["ancestry"])
        ancestry_counts[ancestry] = ancestry_counts.get(ancestry, 0) + 1

    with open(out / "synthetic_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Saved synthetic data to {output_dir}")


# ══════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="All of Us → ExposoGraph Adapter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m ExposoGraph.population_simulation.allofus_adapter \\
      --generate-synthetic 1000 --output ./synthetic_data/
  python -m ExposoGraph.population_simulation.allofus_adapter \\
      --pgx-dir /path/to/pgx/tsvs --ancestry /path/to/ancestry_preds.tsv
  python -m ExposoGraph.population_simulation.allofus_adapter \\
      --show-sql smoking_status
        """
    )
    parser.add_argument("--generate-synthetic", type=int, metavar="N",
                        help="Generate N synthetic participants for testing")
    parser.add_argument("--output", type=str, default="./synthetic_data",
                        help="Output directory for synthetic data")
    parser.add_argument("--pgx-dir", type=str,
                        help="Path to All of Us PGx TSV directory")
    parser.add_argument("--ancestry", type=str,
                        help="Path to ancestry_preds.tsv")
    parser.add_argument("--show-sql", type=str, metavar="QUERY_NAME",
                        help="Print BigQuery SQL template")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for synthetic data")
    args = parser.parse_args()

    if args.show_sql:
        print(get_bigquery_sql(args.show_sql, "YOUR_CDR_DATASET"))
        return

    if args.generate_synthetic:
        participants = generate_synthetic_cohort(
            n=args.generate_synthetic,
            seed=args.seed,
            output_dir=args.output
        )
        print("\n  Example participant:")
        print(json.dumps(participants[0], indent=2))
        return

    if args.pgx_dir:
        pgx_data = load_all_pgx_files(args.pgx_dir)
        ancestry = {}
        if args.ancestry:
            ancestry = load_ancestry_predictions(args.ancestry)

        # Convert first 5 as examples
        for pid in list(pgx_data.keys())[:5]:
            record = build_participant_record(
                person_id=pid,
                pgx_data=pgx_data[pid],
                ancestry_data=ancestry,
            )
            print(json.dumps(record, indent=2))
        return

    parser.print_help()


if __name__ == "__main__":
    main()
