"""Regression tests for data-driven proxy flux configuration."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ExposoGraph import FluxTissueWeightSource, flux_engine

DATA_DIR = Path(__file__).resolve().parents[1] / "ExposoGraph" / "data"


def _load_json(name: str) -> dict:
    with open(DATA_DIR / name, "r") as fh:
        return json.load(fh)


def _resolve_ref(doc: dict, ref: str) -> dict:
    node = doc
    for part in ref.split("."):
        node = node[part]
    return node


def test_proxy_flux_classes_match_dispatch_table():
    params = _load_json("proxy_flux_parameters.json")["classes"]
    assert set(params).issubset(flux_engine._DISPATCH)


def test_proxy_flux_terms_have_resolvable_provenance():
    params = _load_json("proxy_flux_parameters.json")["classes"]
    provenance = _load_json("proxy_flux_provenance.json")

    for class_name, class_cfg in params.items():
        for section in ("activation_terms", "detox_terms", "repair_terms"):
            for term_name, term_cfg in class_cfg.get(section, {}).items():
                ref = term_cfg.get("provenance_ref")
                assert ref, f"{class_name}.{section}.{term_name} is missing provenance_ref"
                entry = _resolve_ref(provenance, ref)
                assert entry["confidence"] == term_cfg["confidence"]
                assert entry["sources"], f"{class_name}.{section}.{term_name} has no sources"

        signal_cfg = class_cfg.get("signal")
        if signal_cfg is not None:
            ref = signal_cfg.get("provenance_ref")
            assert ref, f"{class_name}.signal is missing provenance_ref"
            entry = _resolve_ref(provenance, ref)
            assert entry["confidence"] == signal_cfg["confidence"]
            assert entry["sources"], f"{class_name}.signal has no sources"


def test_proxy_mm_term_preserves_vmax_without_capacity_label():
    flux, gm, tw = flux_engine._compute_proxy_mm_term(
        {"gene": "CYP2E1", "vmax": 10.0, "km": 5.0},
        {"CYP2E1": "NM"},
        tissue="Liver",
        S=1.0,
        tissue_weight_source=FluxTissueWeightSource.CURATED,
    )

    expected = flux_engine.michaelis_menten(1.0, 10.0 * gm * tw, 5.0)
    assert flux == pytest.approx(expected)
