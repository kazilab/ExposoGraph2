"""Smoke tests for population allele-frequency lookups."""

from ExposoGraph import (
    AlleleFrequency,
    Population,
    get_allele_frequencies,
)


def test_allele_frequencies_span_multiple_ancestries_for_cyp1a1():
    freqs = get_allele_frequencies("CYP1A1")

    assert len(freqs) > 0
    assert all(isinstance(f, AlleleFrequency) for f in freqs)
    populations = {f.population for f in freqs}
    assert Population.EUR in populations
    assert any(p != Population.EUR for p in populations)


def test_allele_frequencies_return_empty_for_unknown_gene():
    assert get_allele_frequencies("NOT_A_REAL_GENE") == []
