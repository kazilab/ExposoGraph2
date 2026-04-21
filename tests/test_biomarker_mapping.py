"""Smoke tests for the Supplementary Table S1 biomarker -> [S]/Km mapping."""

import pytest

from ExposoGraph import (
    BiomarkerEntry,
    biomarker_entry_to_dict,
    compute_s_over_km,
    get_biomarker_catalog,
    get_biomarker_entries,
    get_entries_for_lifestyle_factor,
    get_entry_by_biomarker,
    get_lifestyle_factors,
)


def test_catalog_metadata_and_entry_list():
    catalog = get_biomarker_catalog()
    assert "_metadata" in catalog
    assert "entries" in catalog
    assert isinstance(catalog["entries"], list)
    assert catalog["entries"], "catalog should contain at least one entry"


def test_entries_are_typed_biomarker_entries():
    entries = get_biomarker_entries()
    assert all(isinstance(e, BiomarkerEntry) for e in entries)
    assert len(entries) >= 10


def test_lifestyle_factor_filter_returns_smoking_biomarkers():
    smoking = get_entries_for_lifestyle_factor("current_smoking")
    names = {e.biomarker for e in smoking}
    # Canonical smoking biomarkers in the manuscript
    assert {"urinary_1_hydroxypyrene", "serum_cotinine"}.issubset(names)
    assert all(e.lifestyle_factor == "current_smoking" for e in smoking)


def test_lifestyle_factor_filter_is_case_insensitive():
    upper = get_entries_for_lifestyle_factor("CURRENT_SMOKING")
    lower = get_entries_for_lifestyle_factor("current_smoking")
    assert [e.biomarker for e in upper] == [e.biomarker for e in lower]


def test_lookup_by_biomarker_hits_and_misses():
    hit = get_entry_by_biomarker("urinary_1_hydroxypyrene")
    assert hit is not None
    assert hit.target_tissue == "Lung"
    assert hit.target_enzyme == "CYP1A1"
    miss = get_entry_by_biomarker("not_a_biomarker")
    assert miss is None


def test_lifestyle_factors_sorted_and_unique():
    factors = get_lifestyle_factors()
    assert factors == sorted(set(factors))
    assert "current_smoking" in factors


def test_compute_s_over_km_at_midpoint_returns_central_value():
    entry = get_entry_by_biomarker("urinary_1_hydroxypyrene")
    low, high = entry.reference_range
    midpoint = (low + high) / 2.0
    result = compute_s_over_km(entry.biomarker, midpoint)
    assert result == pytest.approx(entry.S_over_Km_central, rel=1e-6)


def test_compute_s_over_km_clamps_to_published_range():
    entry = get_entry_by_biomarker("urinary_1_hydroxypyrene")
    # Huge value clamps to S_over_Km_range max
    high_clip = compute_s_over_km(entry.biomarker, 1_000.0)
    assert high_clip == pytest.approx(entry.S_over_Km_range[1])
    # Near-zero value clamps to S_over_Km_range min
    low_clip = compute_s_over_km(entry.biomarker, 1e-9)
    assert low_clip == pytest.approx(entry.S_over_Km_range[0])


def test_compute_s_over_km_accepts_explicit_entry():
    entry = get_entry_by_biomarker("blood_benzene")
    mid = sum(entry.reference_range) / 2.0
    result = compute_s_over_km("ignored", mid, entry=entry)
    assert result == pytest.approx(entry.S_over_Km_central, rel=1e-6)


def test_compute_s_over_km_rejects_unknown_biomarker():
    with pytest.raises(KeyError):
        compute_s_over_km("nonexistent_biomarker", 1.0)


def test_compute_s_over_km_rejects_negative_measurement():
    with pytest.raises(ValueError):
        compute_s_over_km("urinary_1_hydroxypyrene", -0.01)


def test_biomarker_entry_to_dict_is_json_friendly():
    entry = get_entry_by_biomarker("serum_cotinine")
    payload = biomarker_entry_to_dict(entry)
    assert payload["biomarker"] == "serum_cotinine"
    assert isinstance(payload["reference_range"], list)
    assert len(payload["reference_range"]) == 2
    assert isinstance(payload["S_over_Km_range"], list)
    assert isinstance(payload["references"], list)


def test_every_entry_declares_required_manuscript_fields():
    required_positive = (
        "partition_coefficient",
        "Km_uM",
        "S_over_Km_central",
        "tier2_multiplier",
    )
    for entry in get_biomarker_entries():
        for field in required_positive:
            assert getattr(entry, field) > 0, f"{entry.biomarker}.{field} must be > 0"
        assert entry.reference_range[0] <= entry.reference_range[1]
        assert entry.S_over_Km_range[0] <= entry.S_over_Km_range[1]
