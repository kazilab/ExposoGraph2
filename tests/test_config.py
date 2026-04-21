"""Tests for ExposoGraph.config."""

from ExposoGraph.config import (
    AppMode,
    GraphMode,
    GraphVisibility,
    get_app_mode,
    normalize_app_mode,
    normalize_graph_mode,
    normalize_graph_visibility,
    persistence_enabled,
)


class TestNormalizeAppMode:
    def test_default_is_stateless(self):
        assert normalize_app_mode(None) == AppMode.STATELESS

    def test_public_aliases_map_to_stateless(self):
        for value in ["stateless", "public", "web", "unknown"]:
            assert normalize_app_mode(value) == AppMode.STATELESS

    def test_local_aliases_map_to_local(self):
        for value in ["local", "curation", "persistent"]:
            assert normalize_app_mode(value) == AppMode.LOCAL


class TestGetAppMode:
    def test_reads_from_environment_mapping(self):
        assert get_app_mode({"ExposoGraph_MODE": "local"}) == AppMode.LOCAL
        assert get_app_mode({"ExposoGraph_MODE": "public"}) == AppMode.STATELESS


class TestNormalizeGraphMode:
    def test_default_is_exploratory(self):
        assert normalize_graph_mode(None) == GraphMode.EXPLORATORY

    def test_exploratory_aliases_map_to_exploratory(self):
        for value in ["exploratory", "draft", "flexible", "unknown"]:
            assert normalize_graph_mode(value) == GraphMode.EXPLORATORY

    def test_strict_aliases_map_to_strict(self):
        for value in ["strict", "validated", "canonical"]:
            assert normalize_graph_mode(value) == GraphMode.STRICT


class TestNormalizeGraphVisibility:
    def test_default_is_all(self):
        assert normalize_graph_visibility(None) == GraphVisibility.ALL

    def test_all_aliases_map_to_all(self):
        for value in ["all", "full", "unknown"]:
            assert normalize_graph_visibility(value) == GraphVisibility.ALL

    def test_validated_aliases_map_to_validated_only(self):
        for value in ["validated", "validated_only", "canonical", "strict"]:
            assert normalize_graph_visibility(value) == GraphVisibility.VALIDATED_ONLY

    def test_exploratory_aliases_map_to_exploratory_only(self):
        for value in ["exploratory", "exploratory_only", "provisional"]:
            assert normalize_graph_visibility(value) == GraphVisibility.EXPLORATORY_ONLY


class TestPersistenceEnabled:
    def test_enabled_only_for_local_mode(self):
        assert persistence_enabled(AppMode.LOCAL) is True
        assert persistence_enabled(AppMode.STATELESS) is False
