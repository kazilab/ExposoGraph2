"""NHANES transformation pipeline for long-format biomarker tables."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterable

import pandas as pd

from .biomarker_registry import BIOMARKER_REGISTRY
from .catalog import class_available, get_file_url
from .class_registry import NHANES_CLASS_REGISTRY
from .downloader import read_xpt
from .smoking import current_smoking_flag, recent_smoking_flag
from .weights import add_weight_metadata, combine_2yr_weights


def load_base_covariates(cycle: str, cache_dir: str | Path = "data/nhanes/raw") -> pd.DataFrame:
    demo = read_xpt(get_file_url(cycle, "DEMO"), cache_dir=cache_dir)
    keep = [col for col in ["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "RIDRETH3", "WTMEC2YR"] if col in demo.columns]
    out = demo[keep].copy()
    try:
        smq = read_xpt(get_file_url(cycle, "SMQ"), cache_dir=cache_dir)
        smq_keep = [col for col in ["SEQN", "SMQ020", "SMQ040", "SMD641"] if col in smq.columns]
        out = out.merge(smq[smq_keep], on="SEQN", how="left")
    except Exception as exc:
        warnings.warn(f"Smoking questionnaire unavailable for {cycle}: {exc}", stacklevel=2)
    try:
        creat = read_xpt(get_file_url(cycle, "UCREAT"), cache_dir=cache_dir)
        creat_cols = [col for col in ["SEQN", "URXUCR", "URDUCRLC"] if col in creat.columns]
        out = out.merge(creat[creat_cols], on="SEQN", how="left")
    except Exception as exc:
        warnings.warn(f"Urine creatinine unavailable for {cycle}: {exc}", stacklevel=2)
    return out


def load_exposure_class(cycle: str, class_name: str, cache_dir: str | Path = "data/nhanes/raw") -> pd.DataFrame:
    if class_name not in NHANES_CLASS_REGISTRY:
        raise KeyError(f"Unsupported NHANES class: {class_name}")
    if not class_available(cycle, class_name):
        raise KeyError(f"Class {class_name!r} is not available for cycle {cycle}")
    info = NHANES_CLASS_REGISTRY[class_name]
    return read_xpt(get_file_url(cycle, info["default_file_key"]), cache_dir=cache_dir)


def _creatinine_adjust(value: pd.Series, creatinine_mg_dl: pd.Series, mw_g_mol: float | None) -> pd.Series:
    if mw_g_mol is None:
        return pd.Series(pd.NA, index=value.index, dtype="Float64")
    # For ng/L: (ng/L) / MW = nmol/L. Creatinine mg/dL * 0.0884 = mmol/L.
    creatinine_mmol_l = pd.to_numeric(creatinine_mg_dl, errors="coerce") * 0.0884
    nmol_l = pd.to_numeric(value, errors="coerce") / mw_g_mol
    return (nmol_l / creatinine_mmol_l).astype("Float64")


def class_to_long(df: pd.DataFrame, class_name: str) -> pd.DataFrame:
    if class_name not in NHANES_CLASS_REGISTRY:
        raise KeyError(f"Unsupported NHANES class: {class_name}")
    class_info = NHANES_CLASS_REGISTRY[class_name]
    registry = BIOMARKER_REGISTRY.get(class_name, {})
    rows: list[pd.DataFrame] = []
    for variable, meta in registry.items():
        if variable not in df.columns:
            continue
        piece = pd.DataFrame({
            "SEQN": df["SEQN"],
            "source": "NHANES",
            "nhanes_class": class_name,
            "chemical_class": class_info["chemical_class"],
            "matrix": class_info["matrix"],
            "nhanes_variable": variable,
            "biomarker": meta["biomarker"],
            "parent_compound": meta.get("parent_compound"),
            "value_raw": df[variable],
            "units_expected": meta.get("units_expected"),
            "creatinine_mg_dl": df["URXUCR"] if "URXUCR" in df.columns else pd.NA,
            "age_years": df["RIDAGEYR"] if "RIDAGEYR" in df.columns else pd.NA,
            "sex_code": df["RIAGENDR"] if "RIAGENDR" in df.columns else pd.NA,
            "race_ethnicity_code": df["RIDRETH3"] if "RIDRETH3" in df.columns else df["RIDRETH1"] if "RIDRETH1" in df.columns else pd.NA,
        })
        if class_info["requires_creatinine"]:
            if meta.get("mw_g_mol") is None:
                warnings.warn(f"Missing molecular weight for {meta['biomarker']}; creatinine adjustment skipped", stacklevel=2)
                piece["value_creatinine_adjusted"] = pd.NA
                piece["creatinine_adjusted_units"] = pd.NA
            else:
                piece["value_creatinine_adjusted"] = _creatinine_adjust(piece["value_raw"], piece["creatinine_mg_dl"], meta.get("mw_g_mol"))
                piece["creatinine_adjusted_units"] = "nmol_per_mmol_creatinine"
        else:
            piece["value_creatinine_adjusted"] = pd.NA
            piece["creatinine_adjusted_units"] = pd.NA
        piece["current_smoking"] = current_smoking_flag(df)
        piece["recent_smoking"] = recent_smoking_flag(df)
        piece = add_weight_metadata(piece.join(df[[c for c in df.columns if c.startswith("WT") and c not in piece.columns]], how="left"), class_name)
        rows.append(piece)
    if not rows:
        return pd.DataFrame(columns=[
            "SEQN", "source", "cycle", "nhanes_class", "chemical_class", "matrix", "nhanes_variable",
            "biomarker", "parent_compound", "value_raw", "units_expected", "creatinine_mg_dl",
            "value_creatinine_adjusted", "creatinine_adjusted_units", "current_smoking", "recent_smoking",
            "age_years", "sex_code", "race_ethnicity_code", "nhanes_weight_variable", "nhanes_weight",
            "combined_nhanes_weight",
        ])
    return pd.concat(rows, ignore_index=True)


def build_class_database(cycles: Iterable[str], classes: Iterable[str], cache_dir: str | Path = "data/nhanes/raw", output_csv: str | Path | None = None) -> pd.DataFrame:
    cycles = list(cycles)
    outputs: list[pd.DataFrame] = []
    for cycle in cycles:
        covariates = load_base_covariates(cycle, cache_dir=cache_dir)
        for class_name in classes:
            if class_name not in NHANES_CLASS_REGISTRY:
                warnings.warn(f"Skipping unsupported NHANES class {class_name!r}", stacklevel=2)
                continue
            if not class_available(cycle, class_name):
                warnings.warn(f"Skipping unavailable class/cycle combination: {class_name} @ {cycle}", stacklevel=2)
                continue
            exposure = load_exposure_class(cycle, class_name, cache_dir=cache_dir)
            merged = exposure.merge(covariates, on="SEQN", how="left", suffixes=("", "_cov"))
            long = class_to_long(merged, class_name)
            long["cycle"] = cycle
            outputs.append(long)
    out = pd.concat(outputs, ignore_index=True) if outputs else pd.DataFrame()
    out = combine_2yr_weights(out, n_cycles=len(cycles)) if not out.empty else out
    if output_csv:
        path = Path(output_csv)
        path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(path, index=False)
    return out
