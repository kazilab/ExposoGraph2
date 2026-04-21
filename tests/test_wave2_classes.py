"""Smoke tests for Wave 2 carcinogen class profiles."""

from ExposoGraph import (
    WAVE2_CLASS_PROFILES,
    CarcinogenClassProfile,
    get_all_wave2_classes,
)


def test_wave2_profiles_export_five_classes():
    assert len(WAVE2_CLASS_PROFILES) == 5
    for profile in WAVE2_CLASS_PROFILES.values():
        assert isinstance(profile, CarcinogenClassProfile)


def test_get_all_wave2_classes_matches_profile_keys():
    returned = get_all_wave2_classes()
    assert len(returned) == len(WAVE2_CLASS_PROFILES)
    returned_names = {profile.class_name for profile in returned}
    assert returned_names == set(WAVE2_CLASS_PROFILES.keys())
