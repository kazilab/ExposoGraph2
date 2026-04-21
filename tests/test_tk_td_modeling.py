"""Smoke tests for toxicokinetic / toxicodynamic modeling."""

import math

from ExposoGraph import (
    ExposureRoute,
    ExposureScenario,
    Species,
    TKParameters,
    TKResults,
    one_compartment_model,
)


def _lead_params() -> TKParameters:
    return TKParameters(
        metal="Lead",
        species=Species.HUMAN,
        oral_bioavailability=0.15,
        inhalation_absorption=0.3,
        dermal_absorption=0.01,
        volume_of_distribution=0.3,
        tissue_partition_coeffs={"Liver": 0.5, "Bone": 20.0},
        elimination_rate_constant=0.01,
        clearance=0.05,
        half_life_days=30.0,
        exposure_specific_params={},
    )


def test_one_compartment_pk_curve_rises_from_zero_and_stays_finite():
    scenario = ExposureScenario(
        route=ExposureRoute.ORAL,
        dose_per_day=0.1,
        duration_days=10.0,
        body_weight_kg=70.0,
    )

    results = one_compartment_model(_lead_params(), scenario)

    assert isinstance(results, TKResults)
    assert len(results.time_days) > 0
    assert all(math.isfinite(c) for c in results.blood_concentration)
    assert results.blood_concentration[0] == 0.0
    assert results.blood_concentration[-1] > 0.0
