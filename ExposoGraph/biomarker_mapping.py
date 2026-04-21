"""Supplementary Table S1 biomarker -> [S]/Km mapping for ExposoGraph 2.0.

Provides a thin helper around ``data/biomarker_mapping.json`` so downstream
modules (interaction_engine, exposure_engine, validation notebooks) can
derive tissue-level [S]/Km ratios from reported urinary / blood / serum
biomarker concentrations without duplicating the curated reference ranges
and partition coefficients.

Public entry points:

* :func:`get_biomarker_catalog` - raw catalogue (metadata + entries).
* :func:`get_biomarker_entries` - list of :class:`BiomarkerEntry` dataclasses.
* :func:`get_entries_for_lifestyle_factor` - filter by Tier 2 lifestyle factor.
* :func:`get_entry_by_biomarker` - lookup by biomarker identifier.
* :func:`compute_s_over_km` - convert a measured biomarker value (in the
  biomarker's reference units) into a tissue-level [S]/Km ratio using the
  partition coefficient and Km stored for that biomarker.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

_DATA_FILE = Path(__file__).with_name("data") / "biomarker_mapping.json"
_CATALOG_CACHE: dict[str, Any] | None = None


@dataclass(frozen=True)
class BiomarkerEntry:
    """A single Supplementary Table S1 biomarker -> [S]/Km row."""

    lifestyle_factor: str
    biomarker: str
    matrix: str
    reference_range: tuple[float, float]
    reference_units: str
    partition_coefficient: float
    target_tissue: str
    target_enzyme: str
    carcinogen_class: str
    Km_uM: float
    S_over_Km_central: float
    S_over_Km_range: tuple[float, float]
    tier2_multiplier: float
    references: tuple[str, ...] = field(default_factory=tuple)


def _load_catalog() -> dict[str, Any]:
    global _CATALOG_CACHE
    if _CATALOG_CACHE is None:
        with _DATA_FILE.open("r", encoding="utf-8") as handle:
            _CATALOG_CACHE = json.load(handle)
    return _CATALOG_CACHE


def _to_entry(row: dict[str, Any]) -> BiomarkerEntry:
    ref_range = row.get("reference_range") or [0.0, 0.0]
    sok_range = row.get("S_over_Km_range") or [0.0, 0.0]
    return BiomarkerEntry(
        lifestyle_factor=str(row["lifestyle_factor"]),
        biomarker=str(row["biomarker"]),
        matrix=str(row["matrix"]),
        reference_range=(float(ref_range[0]), float(ref_range[1])),
        reference_units=str(row["reference_units"]),
        partition_coefficient=float(row["partition_coefficient"]),
        target_tissue=str(row["target_tissue"]),
        target_enzyme=str(row["target_enzyme"]),
        carcinogen_class=str(row.get("carcinogen_class", "")),
        Km_uM=float(row["Km_uM"]),
        S_over_Km_central=float(row["S_over_Km_central"]),
        S_over_Km_range=(float(sok_range[0]), float(sok_range[1])),
        tier2_multiplier=float(row["tier2_multiplier"]),
        references=tuple(row.get("references", []) or []),
    )


def get_biomarker_catalog() -> dict[str, Any]:
    """Return a deep-copy-safe view of the Supplementary Table S1 JSON."""
    catalog: dict[str, Any] = json.loads(json.dumps(_load_catalog()))
    return catalog


def get_biomarker_entries() -> list[BiomarkerEntry]:
    """Return every biomarker row as a :class:`BiomarkerEntry`."""
    catalog = _load_catalog()
    return [_to_entry(row) for row in catalog.get("entries", [])]


def get_entries_for_lifestyle_factor(
    lifestyle_factor: str,
) -> list[BiomarkerEntry]:
    """Return the biomarker rows backing a Tier 2 lifestyle multiplier."""
    target = lifestyle_factor.strip().lower()
    return [
        entry
        for entry in get_biomarker_entries()
        if entry.lifestyle_factor.lower() == target
    ]


def get_entry_by_biomarker(biomarker: str) -> BiomarkerEntry | None:
    """Return the first entry matching ``biomarker`` (case-insensitive)."""
    target = biomarker.strip().lower()
    for entry in get_biomarker_entries():
        if entry.biomarker.lower() == target:
            return entry
    return None


def get_lifestyle_factors() -> list[str]:
    """Return the sorted set of Tier 2 lifestyle factors in the catalogue."""
    return sorted({entry.lifestyle_factor for entry in get_biomarker_entries()})


def compute_s_over_km(
    biomarker: str,
    measured_value: float,
    *,
    entry: BiomarkerEntry | None = None,
) -> float:
    """Convert a measured biomarker concentration into a tissue [S]/Km ratio.

    The conversion scales the curated central [S]/Km value by the ratio of
    the measured concentration to the midpoint of the biomarker's reference
    range, then clamps the result to the published [S]/Km_range so a single
    outlier measurement cannot produce a physiologically implausible ratio::

        midpoint = mean(reference_range)
        scale    = measured_value / midpoint
        raw      = S_over_Km_central * scale
        result   = clip(raw, S_over_Km_range.min, S_over_Km_range.max)

    Parameters
    ----------
    biomarker:
        Biomarker identifier (ignored when ``entry`` is provided).
    measured_value:
        Measurement in the biomarker's ``reference_units``.
    entry:
        Optional pre-resolved :class:`BiomarkerEntry`; useful when the
        caller has already looked the row up.

    Returns
    -------
    Tissue-level [S]/Km ratio clipped to the curated plausible range.

    Raises
    ------
    KeyError
        If ``biomarker`` cannot be resolved and ``entry`` is not supplied.
    ValueError
        If ``measured_value`` is negative or the reference midpoint is zero.
    """
    if measured_value < 0:
        raise ValueError(f"measured_value must be non-negative, got {measured_value}")
    row = entry or get_entry_by_biomarker(biomarker)
    if row is None:
        raise KeyError(f"Unknown biomarker: {biomarker!r}")

    low, high = row.reference_range
    midpoint = (low + high) / 2.0
    if midpoint <= 0:
        raise ValueError(
            f"Biomarker {row.biomarker!r} has non-positive reference midpoint "
            f"({low}, {high}); cannot scale [S]/Km."
        )

    scale = measured_value / midpoint
    raw = row.S_over_Km_central * scale
    sok_low, sok_high = row.S_over_Km_range
    return max(sok_low, min(sok_high, raw))


def biomarker_entry_to_dict(entry: BiomarkerEntry) -> dict[str, Any]:
    """Convert a :class:`BiomarkerEntry` to a JSON-friendly dict."""
    payload = asdict(entry)
    payload["reference_range"] = list(entry.reference_range)
    payload["S_over_Km_range"] = list(entry.S_over_Km_range)
    payload["references"] = list(entry.references)
    return payload


__all__ = [
    "BiomarkerEntry",
    "biomarker_entry_to_dict",
    "compute_s_over_km",
    "get_biomarker_catalog",
    "get_biomarker_entries",
    "get_entries_for_lifestyle_factor",
    "get_entry_by_biomarker",
    "get_lifestyle_factors",
]
