"""Smoke tests for the expanded heavy-metal catalogue."""

from ExposoGraph import (
    MetalProfile,
    get_all_expanded_metals,
    get_metal_profile,
)


def test_lead_profile_is_registered_and_typed_correctly():
    lead = get_metal_profile("Lead")

    assert isinstance(lead, MetalProfile)


def test_expanded_metals_catalogue_is_non_empty():
    metals = get_all_expanded_metals()
    assert len(metals) > 0
