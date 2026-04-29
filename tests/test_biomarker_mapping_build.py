"""Tests for rebuilding and comparing the biomarker mapping snapshot."""

from pathlib import Path

from ExposoGraph._biomarker_scaffold.scripts.registries.build_mapping import (
    build_biomarker_mapping_document,
    compare_biomarker_mapping_documents,
    load_biomarker_mapping_manifest,
)
from ExposoGraph._biomarker_scaffold.scripts.registries.loader import load_json_mapping


def test_biomarker_mapping_source_rebuild_matches_old_snapshot():
    root = Path(__file__).resolve().parents[1]
    source_path = (
        root
        / "ExposoGraph"
        / "_biomarker_scaffold"
        / "data"
        / "registries"
        / "biomarkers_master.yaml"
    )
    old_path = root / "ExposoGraph" / "data" / "biomarker_mapping_old.json"

    manifest_metadata, manifest_entries = load_biomarker_mapping_manifest(source_path)
    rebuilt = build_biomarker_mapping_document(source_path)
    old_document = load_json_mapping(old_path)
    report = compare_biomarker_mapping_documents(old_document, rebuilt)
    source_counts = {
        Path(item["path"]).name: item["entry_count"]
        for item in manifest_metadata["source_documents"]
    }

    assert manifest_metadata["source_registry_kind"] == "manifest"
    assert manifest_metadata["source_registry"] == (
        "ExposoGraph/_biomarker_scaffold/data/registries/biomarkers_master.yaml"
    )
    assert manifest_metadata["source_documents"][0]["path"] == (
        "ExposoGraph/_biomarker_scaffold/data/registries/source_nhanes.yaml"
    )
    assert len(manifest_entries) == 15
    assert source_counts["source_nhanes.yaml"] == 6
    assert source_counts["source_literature.yaml"] == 9
    assert source_counts["source_brenda.yaml"] == 0
    assert source_counts["source_pubchem.yaml"] == 0
    assert source_counts["source_comptox.yaml"] == 0
    assert rebuilt["entries"], "rebuilt mapping should contain biomarker entries"
    assert report["mapped_biomarkers_unchanged"] is True
    assert report["added_count"] == 0
    assert report["removed_count"] == 0
    assert report["changed_count"] == 0
