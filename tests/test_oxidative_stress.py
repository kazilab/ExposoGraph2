"""Smoke tests for oxidative stress marker integration."""

from ExposoGraph import (
    OxidativeStressMarker,
    get_markers_by_metal,
)


def test_lead_has_known_oxidative_stress_markers():
    markers = get_markers_by_metal("Lead")

    assert len(markers) > 0
    assert all(isinstance(m, OxidativeStressMarker) for m in markers)


def test_unknown_metal_returns_empty_marker_list():
    assert get_markers_by_metal("Unobtanium") == []
