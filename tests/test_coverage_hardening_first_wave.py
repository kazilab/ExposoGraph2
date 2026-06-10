import pytest

from ExposoGraph._biomarker_scaffold.scripts.registries.mapping_document import (
    apply_update_list,
    build_mapping_document,
    validate_mapping_document,
)
from ExposoGraph.effective_burden import (
    compute_effective_carcinogenic_burden,
    couple_gsh_consumption_to_activation_burden,
)
from ExposoGraph.endpoint_toxic_flux import interpret_competitive_endpoint_flux
from ExposoGraph.exporter import parse_graph_artifact, parse_graph_data_text
from ExposoGraph.gsh_redox_capacity import (
    compute_gsh_redox_capacity,
    compute_quasi_steady_gsh_fraction,
)
from ExposoGraph.interaction_schema import (
    AssumptionWarning,
    EvidenceGrade,
    EvidenceRecord,
    ReactionRole,
    RiskEndpoint,
    SMEReviewStatus,
    enum_from_value,
)
from ExposoGraph.mechanism_attribution import (
    MechanismStateValue,
    compute_mechanism_attribution,
    generate_mechanism_states,
    validate_eight_state_values,
)
from ExposoGraph.model_transparency import (
    AssumptionCategory,
    ReviewSeverity,
    build_transparency_report,
    collect_assumption_warnings,
)
from ExposoGraph.models import Edge, EdgeType, KnowledgeGraph, Node, NodeType
from ExposoGraph.reaction_role_semantics import (
    interpret_competitive_flux_ratio,
    make_unknown_reaction_role_annotation,
)
from ExposoGraph.storage import GraphRepository
from ExposoGraph.unified_api import PatientRiskProfile, summarize_risk_profile


def warning_codes(record):
    return {warning.code for warning in record.warnings}


def test_interaction_schema_enum_fallback_and_nested_serialization():
    assert enum_from_value(ReactionRole, "bioactivation", ReactionRole.UNKNOWN) is ReactionRole.BIOACTIVATION
    assert enum_from_value(ReactionRole, ReactionRole.DETOXIFICATION, ReactionRole.UNKNOWN) is ReactionRole.DETOXIFICATION
    assert enum_from_value(ReactionRole, "not-a-role", ReactionRole.UNKNOWN) is ReactionRole.UNKNOWN
    assert enum_from_value(ReactionRole, None, ReactionRole.UNKNOWN) is ReactionRole.UNKNOWN

    record = EvidenceRecord(
        source="h3-local",
        grade=EvidenceGrade.CURATED,
        metadata={
            "role": ReactionRole.DETOXIFICATION,
            "values": (RiskEndpoint.ROS, {"grade": EvidenceGrade.PLACEHOLDER}),
        },
    )

    assert record.to_dict() == {
        "source": "h3-local",
        "grade": "curated",
        "confidence": None,
        "provenance_ref": None,
        "notes": None,
        "metadata": {
            "role": "detoxification",
            "values": ["ROS", {"grade": "placeholder"}],
        },
    }


def test_reaction_role_unknown_fallback_stays_neutral_and_warned():
    annotation = make_unknown_reaction_role_annotation(
        "CYPX",
        "substrate-y",
        tissue="liver",
        endpoint=RiskEndpoint.ROS,
    )

    assert annotation.active is False
    assert annotation.reaction_role is ReactionRole.UNKNOWN
    assert annotation.sme_notes[0].status is SMEReviewStatus.UNKNOWN
    assert annotation.sme_notes[0].endpoint_context == "ROS"
    assert warning_codes(annotation) == {"reaction_role_unmatched"}

    interpretation = interpret_competitive_flux_ratio(0.25, annotation)

    assert interpretation.burden_multiplier == 1.0
    assert interpretation.flux.flux_decreased is True
    assert warning_codes(interpretation) >= {"reaction_role_unmatched", "reaction_role_inactive"}

    with pytest.raises(ValueError, match="greater than zero"):
        interpret_competitive_flux_ratio(0.0, annotation)


def test_endpoint_toxic_flux_neutral_metadata_and_error_guard():
    annotation = make_unknown_reaction_role_annotation(
        "CYPX",
        "substrate-y",
        tissue="lung",
        endpoint=RiskEndpoint.DNA_ADDUCT,
    )

    result = interpret_competitive_endpoint_flux(
        0.2,
        annotation,
        enzyme="CYPX",
        substrate="substrate-y",
        tissue="lung",
        endpoint=RiskEndpoint.DNA_ADDUCT,
        metadata={"caller": "h3-first-wave"},
    )

    assert result.activation_burden_ratio == 1.0
    assert result.detox_failure_ratio == 1.0
    assert result.endpoint_toxic_flux_ratio == 1.0
    assert result.metadata["interpretation_channel"] == "neutral"
    assert result.metadata["caller_metadata"] == {"caller": "h3-first-wave"}
    assert warning_codes(result) >= {"reaction_role_unmatched", "endpoint_role_inactive"}

    with pytest.raises(ValueError, match="greater than zero"):
        interpret_competitive_endpoint_flux(-1.0, annotation)


def test_gsh_redox_capacity_bounds_unknown_tissue_and_invalid_inputs():
    assert compute_quasi_steady_gsh_fraction(0.0, 0.0) == 0.0
    assert compute_quasi_steady_gsh_fraction(1.0, 0.0) == 1.0
    assert compute_quasi_steady_gsh_fraction(1.0, 3.0) == 0.25

    result = compute_gsh_redox_capacity(
        tissue="not-a-tissue",
        consumption_load=-2.0,
        synthesis_capacity="bad-number",
        turnover_capacity=-1.0,
        baseline_capacity=0.0,
        metadata={"caller": "h3"},
    )

    codes = warning_codes(result)
    assert result.tissue == "default"
    assert result.gsh_fraction == 1.0
    assert result.redox_capacity_ratio == 0.0
    assert result.detox_penalty_multiplier > 1.0
    assert result.clamped is True
    assert result.metadata["caller_metadata"] == {"caller": "h3"}
    assert codes >= {
        "unknown_tissue_default_preset",
        "synthesis_capacity_invalid",
        "consumption_load_negative_clamped",
        "turnover_capacity_negative_clamped",
        "baseline_capacity_zero",
    }


def test_effective_burden_neutral_defaults_and_gsh_coupling_branches():
    result = compute_effective_carcinogenic_burden(
        activation_burden_ratio="bad-number",
        detox_failure_ratio=-2.0,
        susceptibility_modifier=float("inf"),
        gsh_relevant=True,
    )

    codes = warning_codes(result)
    assert result.activation_burden_ratio == 1.0
    assert result.detox_failure_ratio == 0.0
    assert result.susceptibility_modifier == 1.0
    assert result.gsh_detox_penalty_ratio == 1.0
    assert result.effective_carcinogenic_burden_ratio == 0.0
    assert codes >= {
        "activation_burden_ratio_invalid_neutral_default",
        "detox_failure_ratio_negative_clamped",
        "susceptibility_modifier_invalid_neutral_default",
        "gsh_detox_penalty_missing_neutral",
    }

    neutral_coupling = couple_gsh_consumption_to_activation_burden(
        gsh_relevant=False,
        base_gsh_consumption_load=-5.0,
    )
    assert neutral_coupling.scaling_source == "not_gsh_relevant_neutral"
    assert neutral_coupling.gsh_redox_capacity_result is None
    assert neutral_coupling.detox_penalty_multiplier == 1.0
    assert warning_codes(neutral_coupling) == {"base_gsh_consumption_load_negative_clamped"}

    coupled = couple_gsh_consumption_to_activation_burden(
        gsh_relevant=True,
        base_gsh_consumption_load=2.0,
        d_factor=2.0,
        k_factor=3.0,
        tissue="not-a-tissue",
    )
    assert coupled.scaling_source == "d_times_k_approximation"
    assert coupled.upstream_activation_burden_ratio == 6.0
    assert coupled.gsh_consumption_load_scaled == 12.0
    assert coupled.gsh_redox_capacity_result is not None
    assert "unknown_tissue_default_preset" in warning_codes(coupled)


def test_mechanism_attribution_complete_states_and_validation_guards():
    states = generate_mechanism_states()
    state_values = {state: float(index) for index, state in enumerate(states)}

    assert len(states) == 8
    assert validate_eight_state_values(state_values) == state_values

    result = compute_mechanism_attribution(state_values)
    assert result.baseline_value == 0.0
    assert result.full_value == 7.0
    assert result.total_effect == 7.0
    assert result.residuals_are_zero_within_tolerance is True
    assert result.metadata["not_biological_causality"] is True
    assert [item.state_key for item in result.state_values] == [state.key for state in states]

    missing_one = dict(state_values)
    missing_one.pop(states[-1])
    with pytest.raises(ValueError, match="Missing Phase 9 mechanism state values"):
        validate_eight_state_values(missing_one)

    duplicate_values = [
        MechanismStateValue(state=states[0], value=1.0),
        MechanismStateValue(state=states[0], value=2.0),
    ]
    with pytest.raises(ValueError, match="Duplicate Phase 9 mechanism state value"):
        validate_eight_state_values(duplicate_values)

    with pytest.raises(ValueError, match="tolerance"):
        compute_mechanism_attribution(state_values, tolerance=-0.1)


def test_model_transparency_collects_records_and_keeps_no_accepted_caveats():
    warning = AssumptionWarning(
        code="product_carcinogenic_guardrail",
        message="product_carcinogenic evidence is metadata only",
        review_status=SMEReviewStatus.PENDING_TEAM_AGREEMENT,
    )

    records = collect_assumption_warnings({"warnings": [warning]})
    assert len(records) == 1
    assert records[0].category is AssumptionCategory.PRODUCT_CARCINOGENIC_GUARDRAIL
    assert records[0].severity is ReviewSeverity.REVIEW_REQUIRED
    assert records[0].review_status is SMEReviewStatus.PENDING_TEAM_AGREEMENT

    report = build_transparency_report(
        {"warnings": [warning]},
        validation_summary={"targeted_tests": "h3-first-wave"},
    )
    assert report.accepted_non_blocking_caveats == []
    assert report.unresolved_blockers == []
    assert report.metadata["public_adjusted_risk_output"] is False
    assert report.model_card_summary.validation_summary["targeted_tests"] == "h3-first-wave"
    assert report.warning_counts_by_category["product_carcinogenic_guardrail"] == 1


def test_unified_api_summary_reports_public_warning_and_error_contracts():
    profile = PatientRiskProfile(
        tissue="lung",
        genotypes={"CYP1A1": "wildtype"},
        lifestyle={"smoking": True, "alcohol_heavy": False, "alcohol_moderate": False, "occupational_exposure": True},
        exposure_scenario="smoker_industrial_worker",
        exposure_answers={"smoking_status": "current"},
        tissue_report="local deterministic tissue report",
        tissue_weight_count=2,
        top_tissue_genes=[("CYP1A1", 0.876), ("NQO1", 0.321)],
        flux_profile=None,
        exposure_profile=None,
        interactions=None,
        biomarker_dose_estimates=[{"biomarker": "1-OHP", "tissue_conc_uM": 0.01234}],
        pipeline_warnings=["manual review retained"],
        pipeline_errors={"flux": "not available"},
    )

    summary = summarize_risk_profile(profile)

    assert "Patient risk assessment for lung" in summary
    assert "Exposure scenario: smoker_industrial_worker" in summary
    assert "Top tissue-weighted genes: CYP1A1 0.88, NQO1 0.32" in summary
    assert "Biomarker-derived substrate inputs: 1-OHP -> 0.0123 uM" in summary
    assert "Pipeline warnings: manual review retained" in summary
    assert "Pipeline errors: flux=not available" in summary


def test_mapping_document_updates_validation_contracts():
    base_entry = {
        "biomarker": "1-OHP",
        "matrix": "urine",
        "reference_range": "0-1",
        "reference_units": "ug/L",
        "source_status": "measured",
        "lifestyle_factor": "smoking",
    }
    document = build_mapping_document([base_entry], metadata={"owner": "h3"})

    assert validate_mapping_document(document) == []
    assert document["_metadata"]["owner"] == "h3"
    assert document["entries"][0]["entry_id"] == "1-OHP::smoking"
    assert document["entries"][0]["trace"]["created_index"] == 0

    updated = apply_update_list(
        document,
        [
            {"op": "upsert", "entry_id": "1-OHP::smoking", "entry": {"reference_units": "ng/mL"}},
            {"op": "upsert", "entry": {**base_entry, "biomarker": "t,t-MA", "lifestyle_factor": "default"}},
            {"op": "remove", "entry_id": "missing-entry"},
        ],
    )
    assert len(updated["entries"]) == 2
    assert updated["entries"][0]["reference_units"] == "ng/mL"
    assert len(updated["_update_list"]) == 3

    duplicate_document = build_mapping_document([base_entry, dict(base_entry)])
    duplicate_issues = validate_mapping_document(duplicate_document)
    assert duplicate_issues == ["Duplicate biomarker/lifestyle combination not allowed: 1-OHP::smoking"]



def test_exporter_and_storage_round_trip_local_tmp_artifacts(tmp_path):
    js_text = """
    window.GRAPH_DATA = {
      nodes: [
        {id: 'c1', label: 'Benzo[a]pyrene', type: 'Carcinogen'},
        {id: 'g1', label: 'CYP1A1', type: 'Gene'},
      ],
      edges: [
        {source: 'c1', target: 'g1', type: 'ACTIVATES', label: 'activates'},
      ],
    };
    """
    parsed = parse_graph_data_text(js_text)
    assert [node.id for node in parsed.nodes] == ["c1", "g1"]
    assert parsed.edges[0].type is EdgeType.ACTIVATES

    json_path = tmp_path / "graph.json"
    json_path.write_text(parsed.model_dump_json(), encoding="utf-8")
    assert parse_graph_artifact(json_path).nodes[0].label == "Benzo[a]pyrene"

    html_path = tmp_path / "graph.html"
    html_path.write_text(f"<script>{js_text}</script>", encoding="utf-8")
    assert parse_graph_artifact(html_path).edges[0].target == "g1"

    kg = KnowledgeGraph(
        nodes=[
            Node(id="c1", label="Benzo[a]pyrene", type=NodeType.CARCINOGEN),
            Node(id="g1", label="CYP1A1", type=NodeType.GENE),
        ],
        edges=[Edge(source="c1", target="g1", type=EdgeType.ACTIVATES)],
    )

    repository = GraphRepository(tmp_path / "graphs.db")
    first = repository.save_graph(
        graph_key="h3",
        graph_name="H3 Local Graph",
        kg=kg,
        html="<html>first</html>",
        note="first revision",
    )
    second = repository.save_graph(
        graph_key="h3",
        graph_name="H3 Local Graph",
        kg=kg,
        html="<html>second</html>",
    )

    assert first.revision_number == 1
    assert second.revision_number == 2
    assert repository.list_graphs()[0].revision_number == 2
    assert [item.revision_number for item in repository.list_revisions("h3")] == [2, 1]
    latest = repository.get_latest_revision("h3")
    assert latest is not None
    assert latest.to_knowledge_graph().nodes[1].label == "CYP1A1"
    assert repository.get_revision(first.revision_id).note == "first revision"
    repository.close()
    assert repository.connection.execute("SELECT COUNT(*) FROM graphs").fetchone()[0] == 1




