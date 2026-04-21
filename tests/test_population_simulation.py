"""Smoke tests for the synthetic population simulation."""

import random

import pytest

from ExposoGraph import (
    HAPLOTYPE_BLOCKS,
    SyntheticParticipant,
    compute_full_profile,
    generate_synthetic_cohort,
    patient_risk_query,
    sample_haplotype_block,
)
from ExposoGraph.population_simulation._core import _extract_flux_classes
from ExposoGraph.population_simulation.batch_runner import _extract_summary_metrics


def test_generate_synthetic_cohort_produces_requested_size_with_fixed_seed():
    cohort = generate_synthetic_cohort(n=100, seed=42)

    assert len(cohort) == 100
    assert all(isinstance(p, SyntheticParticipant) for p in cohort)
    # Deterministic under fixed seed.
    cohort_repeat = generate_synthetic_cohort(n=100, seed=42)
    assert [p.person_id for p in cohort] == [p.person_id for p in cohort_repeat]


def test_haplotype_block_probabilities_sum_to_one_per_ancestry():
    for block_name, block in HAPLOTYPE_BLOCKS.items():
        for ancestry, entries in block["haplotypes"].items():
            total = sum(entry["p"] for entry in entries)
            assert total == pytest.approx(1.0, abs=0.01), (
                f"{block_name}/{ancestry} sums to {total}"
            )


def test_sample_haplotype_block_returns_all_genes_in_block():
    rng = random.Random(7)
    for block_name, block in HAPLOTYPE_BLOCKS.items():
        sample = sample_haplotype_block(rng, block_name, "European")
        assert set(sample.keys()) == set(block["genes"])


def test_cohort_with_haplotypes_includes_linked_genes():
    cohort = generate_synthetic_cohort(n=20, seed=5, use_haplotypes=True)
    assert cohort
    sample = cohort[0].genotypes
    # Linked genes not present in default ALLELE_FREQUENCIES should be set
    for linked in ("CYP1A2", "GSTM3", "NAT1", "CYP2C8"):
        assert linked in sample


def test_cohort_default_is_backward_compatible():
    cohort = generate_synthetic_cohort(n=20, seed=5)
    sample = cohort[0].genotypes
    # Without haplotype sampling, those linked genes should not appear
    for linked in ("CYP1A2", "GSTM3", "NAT1", "CYP2C8"):
        assert linked not in sample


def test_extract_flux_classes_preserves_model_metadata():
    flux_profile = compute_full_profile(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Liver",
    )

    flux_classes = _extract_flux_classes(flux_profile)

    assert flux_classes["PAH"]["model_kind"] == "measured_kinetics"
    assert flux_classes["PAH"]["parameter_source"] == "kinetic_parameters.json"
    assert flux_classes["Dioxin"]["model_kind"] == "receptor_mediated_proxy"
    assert flux_classes["Dioxin"]["parameter_source"] == "proxy_flux_parameters.json"


def test_batch_runner_summary_keeps_measured_and_proxy_flux_support_split():
    participant = {
        "person_id": "P1",
        "ancestry": "EUR",
        "exposure_scenario": "general_population",
        "genotypes": {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
    }
    profile = patient_risk_query(participant["genotypes"], tissue="Liver")

    summary = _extract_summary_metrics(profile, participant)

    assert summary["flux_classes"]["PAH"]["model_kind"] == "measured_kinetics"
    assert summary["flux_classes"]["PAH"]["parameter_source"] == "kinetic_parameters.json"
    assert summary["flux_classes"]["Dioxin"]["model_kind"] == "receptor_mediated_proxy"
    assert summary["flux_classes"]["Dioxin"]["parameter_source"] == "proxy_flux_parameters.json"
    assert "PAH" in summary["measured_high_risk_pathways"]
    assert "Dioxin" in summary["proxy_high_risk_pathways"]
