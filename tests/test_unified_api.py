"""Smoke tests for the unified high-level risk API."""

from dataclasses import asdict

from ExposoGraph import (
    FluxClassEvidence,
    PatientRiskProfile,
    build_reference_graph,
    enrich_knowledge_graph,
    patient_risk_query,
)


def test_patient_risk_query_returns_profile_for_tier1_genotype():
    profile = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Liver",
    )

    assert isinstance(profile, PatientRiskProfile)
    assert profile.tissue == "Liver"
    assert profile.genotypes.get("CYP1A1") == "NM"


def test_patient_risk_query_exposes_flux_model_evidence_and_summary_labels():
    profile = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Liver",
    )

    assert profile.flux_class_evidence
    assert all(isinstance(item, FluxClassEvidence) for item in profile.flux_class_evidence)

    evidence_by_class = {
        item.carcinogen_class: item
        for item in profile.flux_class_evidence
    }
    assert evidence_by_class["PAH"].model_kind == "measured_kinetics"
    assert evidence_by_class["PAH"].parameter_source == "kinetic_parameters.json"
    assert evidence_by_class["Dioxin"].model_kind == "receptor_mediated_proxy"
    assert evidence_by_class["Dioxin"].parameter_source == "proxy_flux_parameters.json"

    payload = asdict(profile)
    assert payload["flux_class_evidence"]
    assert "Measured-kinetics flux support:" in profile.summary
    assert "Proxy-backed flux support:" in profile.summary
    assert "Dioxin (receptor-mediated proxy)" in profile.summary


def test_patient_risk_query_uses_biomarker_substrate_inputs_for_flux_and_exposure():
    profile = patient_risk_query(
        {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
        tissue="Lung",
        biomarker_measurements={"urinary_1_hydroxypyrene": 0.175},
        include_interactions=False,
        include_tissue_report=False,
    )

    assert profile.biomarker_dose_estimates
    pah_flux = profile.flux_profile.per_class_results["PAH"]
    assert pah_flux.substrate_concentration_uM == profile.biomarker_dose_estimates[0]["tissue_conc_uM"]

    pah_risk = next(
        risk for risk in profile.exposure_profile.risk_scores if risk.carcinogen_class == "PAH"
    )
    assert pah_risk.exposure_tier == 3
    assert pah_risk.biomarker_dose_estimate["biomarker"] == "urinary_1_hydroxypyrene"


def _find_edge(graph, source: str, target: str):
    for edge in graph.edges:
        if edge.source == source and edge.target == target:
            return edge
    raise AssertionError(f"Could not find edge {source!r} -> {target!r}")


def test_enrich_knowledge_graph_attaches_wave2_kinetics_to_reference_edges():
    enriched = enrich_knowledge_graph(
        build_reference_graph(),
        attach_tissue_weights=False,
        attach_exposure_data=False,
    )

    tce_edge = _find_edge(enriched, "CYP2E1", "Chloral_hydrate")
    assert tce_edge.kinetics is not None
    assert tce_edge.kinetics["ChlorinatedSolvent"]["pathway_type"] == "oxidation"
    assert tce_edge.kinetics["ChlorinatedSolvent"]["matched_enzyme_key"] == "CYP2E1_TCE"

    tcdd_edge = _find_edge(enriched, "TCDD", "AHR")
    assert tcdd_edge.kinetics is not None
    assert tcdd_edge.kinetics["Dioxin"]["pathway_type"] == "AhR_signaling"
    assert tcdd_edge.kinetics["Dioxin"]["matched_enzyme_key"] == "AhR_binding"


def test_dioxin_binding_kinetics_do_not_leak_to_downstream_signaling_edges():
    enriched = enrich_knowledge_graph(
        build_reference_graph(),
        attach_tissue_weights=False,
        attach_exposure_data=False,
    )

    downstream_edge = _find_edge(enriched, "AHR", "ARNT")
    assert downstream_edge.kinetics is None or "Dioxin" not in downstream_edge.kinetics
