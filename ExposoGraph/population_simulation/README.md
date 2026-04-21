# Module 07: Population-Scale Simulation

## Overview
Runs ExposoGraph risk predictions across the entire All of Us Research Program cohort (~414,000 participants with WGS) and validates predictions against observed cancer outcomes in linked EHR data.

## Architecture

```
ExposoGraph/population_simulation/
├── __init__.py              # Unified package surface and exports
├── _core.py                 # Canonical typed synthetic cohort API
├── allofus_adapter.py       # AoU/source-compatible wrappers and SQL templates
├── phenotype_extractor.py   # Canonical phenotype catalog + source-style wrappers
├── batch_runner.py          # Parallelized simulation with checkpointing
├── population_analysis.py   # Aggregate statistics, ancestry stratification
├── validation_framework.py  # Canonical reference ORs, ROC/AUC, calibration
└── README.md
```

## Canonical Layout

This package no longer treats the old extension folder and the typed package
API as separate implementations.

- `_core.py` is the canonical source for typed synthetic participants, typed
  summaries, and the top-level helpers exported from `ExposoGraph`.
- `phenotype_extractor.py` is the canonical source for the richer cancer
  phenotype catalog, including SNOMED mappings and the `any_cancer` umbrella.
- `validation_framework.py` is the canonical source for the published
  genotype-exposure-cancer reference set.
- `allofus_adapter.generate_synthetic_cohort()` and
  `phenotype_extractor.generate_synthetic_cancer_labels()` are source-style
  wrappers over that canonical package logic, so the dict-based workflow stays
  aligned with the typed API.

## Quick Start (Testing Outside Workbench)

```bash
# Generate 100 synthetic participants and run full pipeline
python -m ExposoGraph.population_simulation.batch_runner --synthetic 100 --workers 1 --tissue Liver

# Analyze results
python -m ExposoGraph.population_simulation.population_analysis \
    --results ./simulation_output/results.jsonl --output ./simulation_output/

# Run validation against synthetic cancer labels
python -m ExposoGraph.population_simulation.validation_framework \
    --results ./simulation_output/results.jsonl \
    --cancer-labels ./simulation_output/cancer_labels.json
```

`results.jsonl` stores one compact per-participant summary per line. Each
`flux_classes` entry now preserves `net_ratio`, risk label, `model_kind`, and
`parameter_source`, and the batch summary also separates
`measured_high_risk_pathways` from `proxy_high_risk_pathways` so downstream
analysis can distinguish literature-backed kinetics from proxy-backed classes.

## Running in All of Us Researcher Workbench

### Step 1: Set Up Workspace
1. Log into [Researcher Workbench](https://workbench.researchallofus.org)
2. Create a new workspace (Controlled Tier access required for genomic data)
3. Upload or install the `ExposoGraph/` package in your workspace environment
4. Launch a Jupyter Notebook (Python 3, recommend 16 CPUs / 104 GB RAM)

### Step 2: Extract Data
```python
import os

# The CDR dataset name is available as an environment variable
CDR = os.environ.get("WORKSPACE_CDR")

# Extract PGx data
from ExposoGraph.population_simulation.allofus_adapter import get_bigquery_sql
import pandas as pd

# Smoking status
smoking_sql = get_bigquery_sql("smoking_status", CDR)
smoking_df = pd.read_gbq(smoking_sql, dialect="standard")

# Alcohol use  
alcohol_sql = get_bigquery_sql("alcohol_use", CDR)
alcohol_df = pd.read_gbq(alcohol_sql, dialect="standard")

# Cancer diagnoses
from ExposoGraph.population_simulation.phenotype_extractor import get_phenotype_sql
cancer_sql = get_phenotype_sql("cancer_conditions", CDR)
cancer_df = pd.read_gbq(cancer_sql, dialect="standard")
```

### Step 3: Build Participant Records
```python
from ExposoGraph.population_simulation.allofus_adapter import (
    load_all_pgx_files,
    load_ancestry_predictions,
    build_participant_record,
    extract_lifestyle_from_surveys,
)

# Load PGx data (path from Controlled CDR Directory)
pgx_data = load_all_pgx_files("/path/to/pgx/tsvs/")
ancestry = load_ancestry_predictions("/path/to/ancestry_preds.tsv")

# Build participant records
participants = []
for person_id, pgx in pgx_data.items():
    # Get lifestyle from survey data
    person_surveys = smoking_df[smoking_df.person_id == int(person_id)]
    lifestyle = extract_lifestyle_from_surveys(person_id, person_surveys)
    
    record = build_participant_record(
        person_id=person_id,
        pgx_data=pgx,
        ancestry_data=ancestry,
        lifestyle_data=lifestyle,
    )
    participants.append(record)
```

### Step 4: Run Simulation
```python
from ExposoGraph.population_simulation.batch_runner import PopulationSimulation

config = {
    "max_workers": 16,       # Match your VM CPU count
    "batch_size": 1000,
    "output_dir": "./simulation_output",
    "checkpoint_interval": 500,
}

sim = PopulationSimulation(participants, config)
sim.run(tissue="Liver")      # Run for Liver
sim.run(tissue="Lung")       # Run for Lung (resumes, skips completed)
sim.save_run_report()
```

### Step 5: Analyze and Validate
```python
from ExposoGraph.population_simulation.population_analysis import generate_population_report, load_results
from ExposoGraph.population_simulation.phenotype_extractor import build_cancer_labels
from ExposoGraph.population_simulation.validation_framework import run_full_validation

# Load results
results = load_results("./simulation_output/results.jsonl")

# Build cancer outcome labels from EHR
cancer_labels = build_cancer_labels(cancer_df.to_dict('records'))

# Population analysis
report = generate_population_report(results, cancer_labels)

# Validation against published ORs
validation = run_full_validation(results, cancer_labels)
```

## All of Us Data Requirements

| Data Source | CDR Table / File | Tier | Used By |
|-------------|-----------------|------|---------|
| PGx star alleles | pgx/ auxiliary TSVs (18 genes) | Controlled | allofus_adapter |
| Genetic ancestry | ancestry_preds.tsv | Controlled | allofus_adapter |
| Smoking status | observation (concept 40766929) | Registered | allofus_adapter |
| Alcohol use | observation (concept 40771103) | Registered | allofus_adapter |
| Occupation | observation (concept 40771090) | Registered | allofus_adapter |
| Cancer diagnoses | condition_occurrence (SNOMED 443392 descendants) | Controlled | phenotype_extractor |
| Demographics | person table | Registered | allofus_adapter |
| SV/CNV (GSTM1/GSTT1) | structural_variants | Controlled | allofus_adapter (optional) |

## ExposoGraph Genes Covered

| Gene | All of Us Source | ExposoGraph Engine |
|------|-----------------|-------------------|
| CYP2D6 | PGx auxiliary (Cyrius) | Flux model |
| CYP3A5 | PGx auxiliary (Stargazer) | Flux model |
| CYP2B6 | PGx auxiliary (Stargazer) | Flux model |
| CYP2C9 | PGx auxiliary (Stargazer) | Flux model |
| CYP2C19 | PGx auxiliary (Stargazer) | Flux model |
| UGT1A1 | PGx auxiliary (Stargazer) | Flux model |
| GSTM1 | WGS SV calls (derived) | Flux, Interaction |
| GSTT1 | WGS SV calls (derived) | Flux, Interaction |
| CYP1A1 | Not in panel (default NM) | Flux, Interaction |
| CYP1A2 | Not in panel (default NM) | Interaction |
| CYP2E1 | Not in panel (default NM) | Interaction |
| ALDH2 | Not in panel (default *1/*1) | Interaction |
| NAT2 | Not in panel (default rapid) | Flux |

**Note**: CYP1A1, CYP1A2, CYP2E1, ALDH2, and NAT2 are not in the All of Us PGx panel.
For these genes, the adapter defaults to wild-type. However, their star alleles CAN be called
from the raw WGS VCF using Stargazer or PharmCAT. A future enhancement would extract these
directly from the srWGS variant data.

## Validation References

The validation framework compares observed ORs against 8 published genotype-exposure-cancer
associations from IARC and meta-analyses:

| Association | Published OR | Source |
|-------------|-------------|--------|
| GSTM1-null × smoking → lung | 1.58 (1.21-2.06) | Ye et al. 2006 |
| CYP1A1*2A × smoking → lung | 2.36 (1.16-4.81) | Shi et al. 2008 |
| GSTM1+CYP1A1 combined → lung | 2.87 (1.73-4.76) | Vineis et al. 2007 |
| ALDH2*1/*2 × alcohol → esophageal | 6.97 (4.36-11.12) | Yokoyama & Omori 2005 |
| ALDH2*2/*2 × alcohol → esophageal | 12.5 (6.0-26.0) | Yokoyama & Omori 2005 |
| NAT2 slow × smoking → bladder | 1.51 (1.28-1.78) | Garcia-Closas et al. 2005 |
| GSTM1-null × smoking → bladder | 1.53 (1.11-2.12) | Engel et al. 2002 |
| GSTT1-active × TCE → kidney | 1.88 (1.06-3.33) | Karami et al. 2012 |

## Compute Budget Estimates

| Cohort Size | Workers | Estimated Time | Estimated GCP Cost |
|-------------|---------|---------------|-------------------|
| 1,000 (test) | 4 | ~5 min | $0.10 |
| 10,000 | 8 | ~30 min | $1-2 |
| 100,000 | 16 | ~4 hours | $10-20 |
| 414,000 (full) | 32 | ~12-16 hours | $40-80 |

Based on 16 CPUs / 104 GB RAM VM at ~$2.50/hr on GCP.
