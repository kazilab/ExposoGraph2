"""NHANES survey weight selection helpers."""

from __future__ import annotations

import warnings

import pandas as pd

CLASS_WEIGHT_CANDIDATES = {
    "PAH": ["WTSPAH2YR", "WTSA2YR", "WTMEC2YR"],
    "COTININE": ["WTSS2YR", "WTMEC2YR"],
    "METALS_BLOOD": ["WTSB2YR", "WTMEC2YR"],
    "METALS_URINE": ["WTSA2YR", "WTMEC2YR"],
    "PFAS": ["WTSA2YR", "WTMEC2YR"],
    "PHTHALATES": ["WTSA2YR", "WTMEC2YR"],
    "PHENOLS": ["WTSA2YR", "WTMEC2YR"],
    "VOC_METABOLITES": ["WTSA2YR", "WTMEC2YR"],
    "VOC_BLOOD": ["WTSB2YR", "WTMEC2YR"],
}
DEFAULT_WEIGHT_CANDIDATES = ["WTSA2YR", "WTSB2YR", "WTSPAH2YR", "WTSS2YR", "WTMEC2YR"]


def select_weight_variable(df: pd.DataFrame, class_name: str | None = None) -> str | None:
    candidates = CLASS_WEIGHT_CANDIDATES.get(class_name or "", DEFAULT_WEIGHT_CANDIDATES)
    for candidate in candidates:
        if candidate in df.columns:
            if candidate == "WTMEC2YR" and class_name:
                warnings.warn(f"Using MEC weight for {class_name}; check whether a class-specific subsample weight is required", stacklevel=2)
            return candidate
    return None


def add_weight_metadata(df: pd.DataFrame, class_name: str | None = None) -> pd.DataFrame:
    out = df.copy()
    weight = select_weight_variable(out, class_name)
    out["nhanes_weight_variable"] = weight
    out["nhanes_weight"] = out[weight] if weight else pd.NA
    return out


def combine_2yr_weights(df: pd.DataFrame, n_cycles: int) -> pd.DataFrame:
    if n_cycles <= 0:
        raise ValueError("n_cycles must be positive")
    out = df.copy()
    if "nhanes_weight" not in out.columns:
        out = add_weight_metadata(out)
    out["combined_nhanes_weight"] = pd.to_numeric(out["nhanes_weight"], errors="coerce") / n_cycles
    return out
