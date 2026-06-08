from ExposoGraph.interaction_schema import (
    AssumptionWarning,
    EvidenceGrade,
    EvidenceRecord,
    ReleaseTarget,
    SMEReviewStatus,
)
import ExposoGraph.model_transparency as model_transparency
from ExposoGraph.model_transparency import (
    AssumptionCategory,
    ReviewSeverity,
    TransparencyRecord,
    build_deferral_records,
    build_model_card_summary,
    build_sme_review_queue,
    build_transparency_report,
    collect_assumption_warnings,
    summarize_model_boundaries,
    transparency_report_from_phase_outputs,
)


def _warning(code, message=None, *, status=SMEReviewStatus.UNKNOWN, severity="warning"):
    return AssumptionWarning(
        code=code,
        message=message or code,
        severity=severity,
        review_status=status,
    )


def _messages(items):
    return [item.message for item in items]


def test_transparency_records_serialize_cleanly():
    record = TransparencyRecord(
        category=AssumptionCategory.REACTION_ROLE,
        message="Unknown reaction role remains unresolved.",
        severity=ReviewSeverity.REVIEW_REQUIRED,
        source_phase="phase_4_reaction_role_semantics",
        review_status=SMEReviewStatus.UNKNOWN,
        release_target=ReleaseTarget.V2_0,
        code="reaction_role_unmatched",
        evidence=EvidenceRecord(source="local test", grade=EvidenceGrade.PLACEHOLDER),
    )

    payload = record.to_dict()

    assert payload["category"] == "reaction_role"
    assert payload["severity"] == "review_required"
    assert payload["review_status"] == "unknown"
    assert payload["release_target"] == "v2_0"
    assert payload["evidence"]["grade"] == "placeholder"


def test_warnings_are_grouped_by_category_and_severity():
    records = collect_assumption_warnings(
        [
            _warning("reaction_role_unmatched", "No role."),
            _warning("km_used_as_ki_proxy", "Km proxy.", status=SMEReviewStatus.CURATED),
            _warning("ic50_conversion_unavailable", "IC50 unavailable."),
        ]
    )
    by_category = {record.category for record in records}
    by_severity = {record.severity for record in records}

    assert AssumptionCategory.REACTION_ROLE in by_category
    assert AssumptionCategory.KINETIC_PARAMETER in by_category
    assert AssumptionCategory.IC50_CONVERSION in by_category
    assert ReviewSeverity.REVIEW_REQUIRED in by_severity
    assert ReviewSeverity.WARNING in by_severity


def test_sme_review_queue_includes_unknown_reaction_role_warning():
    queue = build_sme_review_queue(
        [_warning("reaction_role_unmatched", "No explicit reaction-role record matched.")]
    )

    assert any(item.code == "reaction_role_unmatched" for item in queue)
    assert any(item.category is AssumptionCategory.REACTION_ROLE for item in queue)


def test_sme_review_queue_includes_tce_cyp2e1_pending_candidate_item():
    queue = build_sme_review_queue()

    tce = [item for item in queue if item.code == "spyros_tce_cyp2e1_candidate_pending"]

    assert len(tce) == 1
    assert tce[0].review_status is SMEReviewStatus.PENDING_TEAM_AGREEMENT
    assert tce[0].metadata["candidate_pending"] is True
    assert "TCE x CYP2E1" in tce[0].message


def test_deferred_nat2_aldh_gst_ephx1_items_are_visible():
    deferrals = build_deferral_records()
    subjects = {item.subject for item in deferrals}

    assert "NAT2 broader detox-arm semantics" in subjects
    assert "ALDH broader detox-arm semantics" in subjects
    assert "GST broader detox-arm semantics" in subjects
    assert "EPHX1 broader detox-arm semantics" in subjects
    assert all(item.release_target in {ReleaseTarget.V3_0, ReleaseTarget.FUTURE} for item in deferrals)


def test_km_proxy_warning_appears_as_kinetic_parameter_review_item():
    queue = build_sme_review_queue(
        [_warning("km_used_as_ki_proxy", "Km_uM was used as a low-confidence Ki proxy.")]
    )
    km_items = [item for item in queue if item.code == "km_used_as_ki_proxy"]

    assert len(km_items) == 1
    assert km_items[0].category is AssumptionCategory.KINETIC_PARAMETER
    assert km_items[0].metadata["low_confidence_proxy"] is True


def test_ic50_guard_warning_appears_and_no_conversion_is_claimed():
    records = collect_assumption_warnings(
        [_warning("ic50_conversion_unavailable", "IC50-to-Ki conversion was requested but unavailable.")]
    )
    ic50 = [record for record in records if record.category is AssumptionCategory.IC50_CONVERSION]

    assert len(ic50) == 1
    assert ic50[0].metadata["conversion_claimed"] is False


def test_affinity_rdkit_tanimoto_unavailable_status_appears_without_conversion_claim():
    report = build_transparency_report(
        [_warning("affinity_fallback_unavailable", "RDKit ECFP4 Tanimoto affinity fallback unavailable.")]
    )
    records = [record for record in report.records if record.category is AssumptionCategory.AFFINITY_FALLBACK]
    deferrals = [item for item in report.deferrals if item.subject == "RDKit/ECFP4/Tanimoto affinity fallback"]

    assert records[0].metadata["conversion_claimed"] is False
    assert records[0].metadata["rdkit_ecfp4_tanimoto_active"] is False
    assert deferrals[0].metadata["active_in_2_0"] is False


def test_gsh_tissue_presets_appear_as_local_2_0_preset_caveats():
    queue = build_sme_review_queue()
    preset_items = [item for item in queue if item.code == "local_2_0_gsh_preset"]

    assert len(preset_items) == 1
    assert preset_items[0].category is AssumptionCategory.GSH_TISSUE_PRESET
    assert preset_items[0].metadata["local_2_0_preset"] is True


def test_model_boundary_states_gsh_is_not_pbpk_ode_gsh_gssg_nrf2_validated():
    boundaries = summarize_model_boundaries()
    text = " ".join(boundary.statement for boundary in boundaries)

    assert "not PBPK/ODE/GSH-GSSG/Nrf2 validated" in text


def test_model_boundary_states_shapley_is_model_output_attribution_not_causality():
    boundaries = summarize_model_boundaries()
    text = " ".join(boundary.statement for boundary in boundaries)

    assert "model-output attribution only" in text
    assert "not proof of biological causality" in text


def test_model_card_includes_context_limitations_assumptions_review_queue_and_deferrals():
    card = build_model_card_summary(
        [_warning("reaction_role_unmatched", "Unknown reaction role.")],
        validation_summary={"targeted_tests": "not run in unit fixture"},
    )

    assert "research-oriented comparison" in card.intended_context_of_use
    assert card.known_limitations
    assert card.major_assumptions
    assert card.sme_review_items
    assert card.deferred_items
    assert card.validation_summary["targeted_tests"] == "not run in unit fixture"


def test_product_carcinogenic_guardrail_appears_in_model_boundary_statements():
    boundaries = summarize_model_boundaries()
    guardrails = [item for item in boundaries if item.category is AssumptionCategory.PRODUCT_CARCINOGENIC_GUARDRAIL]

    assert len(guardrails) == 1
    assert "not the final sign rule" in guardrails[0].statement


def test_no_adjusted_risk_public_output_is_produced():
    report = build_transparency_report()
    payload = report.to_dict()

    assert "adjusted_risk" not in payload
    assert payload["metadata"]["public_adjusted_risk_output"] is False
    assert report.model_card_summary.metadata["public_adjusted_risk_output"] is False


def test_no_interaction_engine_integration_occurs():
    assert "interaction_engine" not in model_transparency.__dict__
    assert "compute_interaction_matrix" not in model_transparency.__dict__
    assert "competitive_inhibition_flux" not in model_transparency.__dict__
    assert "get_ki" not in model_transparency.__dict__


def test_no_phase11_release_packaging_or_github_behavior_occurs():
    report = transparency_report_from_phase_outputs()

    assert report.metadata["release_packaging_behavior"] is False
    assert report.metadata["github_behavior"] is False
    assert report.metadata["phase11_behavior"] is False
    assert "github" not in " ".join(dir(model_transparency)).lower()


def test_deterministic_output_ordering():
    warnings = [
        _warning("ic50_conversion_unavailable", "IC50 unavailable."),
        _warning("reaction_role_unmatched", "Unknown role."),
        _warning("affinity_fallback_unavailable", "Affinity unavailable."),
    ]

    first = build_transparency_report(warnings).to_dict()
    second = build_transparency_report(list(reversed(warnings))).to_dict()

    assert first["records"] == second["records"]
    assert first["sme_review_queue"] == second["sme_review_queue"]
    assert first["deferrals"] == second["deferrals"]
