"""Smoking status helpers for NHANES covariates."""

from __future__ import annotations

import pandas as pd


def _nullable_false(index) -> pd.Series:
    return pd.Series(pd.NA, index=index, dtype="boolean")


def current_smoking_flag(df: pd.DataFrame) -> pd.Series:
    """Return nullable boolean for current smoking.

    SMQ020: smoked at least 100 cigarettes in life (1=yes, 2=no)
    SMQ040: now smokes cigarettes (1=every day, 2=some days, 3=not at all)
    """
    if "SMQ020" not in df or "SMQ040" not in df:
        return _nullable_false(df.index)
    ever = df["SMQ020"].eq(1)
    current = df["SMQ040"].isin([1, 2])
    result = ever & current
    missing = df[["SMQ020", "SMQ040"]].isna().any(axis=1)
    return result.mask(missing, pd.NA).astype("boolean")


def recent_smoking_flag(df: pd.DataFrame) -> pd.Series:
    """Return nullable boolean for cigarette smoking in the past 5 days."""
    if "SMD641" in df:
        # NHANES coding varies by cycle; treat positive numeric counts as recent use.
        result = pd.to_numeric(df["SMD641"], errors="coerce").gt(0)
        return result.mask(df["SMD641"].isna(), pd.NA).astype("boolean")
    return current_smoking_flag(df)
