"""Smoke tests for cross-species extrapolation."""

from ExposoGraph import allometric_scale


def test_rat_to_human_allometric_scaling_increases_dose_with_bw_ratio():
    scaled = allometric_scale(1.0, from_species_bw=0.25, to_species_bw=70.0)

    # With the default exponent 0.75, a rat→human scale on body-weight ratio
    # 280× should produce a finite multiplier well above 1.
    assert scaled > 1.0
    assert scaled < 500.0


def test_allometric_scale_is_identity_for_equal_body_weights():
    assert allometric_scale(1.0, from_species_bw=70.0, to_species_bw=70.0) == 1.0
