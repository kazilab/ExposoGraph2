import json
from types import SimpleNamespace

import pytest

from ExposoGraph._biomarker_scaffold.scripts.registries.build_mapping import (
    build_biomarker_mapping_document,
    compare_biomarker_mapping_documents,
    load_biomarker_mapping_manifest,
)
from ExposoGraph._biomarker_scaffold.scripts.registries.check_mapping import (
    validate_biomarker_mapping_document,
)
from ExposoGraph._biomarker_scaffold.scripts.registries.loader import (
    load_json_mapping,
    load_registry_document,
    write_json_mapping,
)
from ExposoGraph.exporter import (
    bundle_to_html_string,
    parse_graph_artifact,
    parse_graph_data_text,
)
from ExposoGraph.interaction_schema import (
    EvidenceGrade,
    ReactionRole,
    ReleaseTarget,
    RiskDirectionIfFluxDecreases,
    SMEReviewStatus,
)
from ExposoGraph.model_transparency import (
    AssumptionCategory,
    ReviewSeverity,
    build_model_card_summary,
    build_sme_review_queue,
    build_transparency_report,
    collect_assumption_warnings,
)
from ExposoGraph.parameter_provider import (
    JSONInteractionParameterProvider,
    KGInteractionParameterProvider,
)
from ExposoGraph.unified_api import PatientRiskProfile, compare_patient_profiles


def _warning_codes(record):
    return {warning.code for warning in record.warnings or []}


def test_model_transparency_dedupes_sorts_counts_and_preserves_no_accepted_caveats():
    warnings = [
        {
            "code": "z_gsh_preset",
            "message": "GSH tissue preset needs review",
            "review_status": SMEReviewStatus.UNKNOWN.value,
            "release_target": ReleaseTarget.V2_0.value,
        },
        {
            "code": "a_ic50_conversion",
            "message": "IC50 conversion error requires review",
            "severity": "error",
            "review_status": SMEReviewStatus.PENDING_TEAM_AGREEMENT.value,
            "release_target": ReleaseTarget.FUTURE.value,
        },
        {
            "code": "z_gsh_preset",
            "message": "GSH tissue preset needs review",
            "review_status": SMEReviewStatus.UNKNOWN.value,
            "release_target": ReleaseTarget.V2_0.value,
        },
    ]
    phase_output = {"warnings": warnings}

    records = collect_assumption_warnings(phase_output)

    assert len(records) == 2
    assert [(record.category, record.code) for record in records] == [
        (AssumptionCategory.GSH_TISSUE_PRESET, "z_gsh_preset"),
        (AssumptionCategory.IC50_CONVERSION, "a_ic50_conversion"),
    ]
    assert records[0].severity is ReviewSeverity.REVIEW_REQUIRED
    assert records[1].severity is ReviewSeverity.BLOCKER

    queue = build_sme_review_queue(
        phase_output,
        include_phase4_registry=False,
        include_model_boundary_caveats=False,
    )
    assert [item.code for item in queue] == ["z_gsh_preset", "a_ic50_conversion"]
    assert {item.release_target for item in queue} == {ReleaseTarget.V2_0, ReleaseTarget.FUTURE}

    card = build_model_card_summary(phase_output, validation_summary={"h4a": "local"})
    assert card.validation_summary == {"h4a": "local"}
    assert card.warning_counts_by_category["gsh_tissue_preset"] == 1
    assert card.warning_counts_by_category["ic50_conversion"] == 1
    assert card.release_target_summary["future"] >= 1
    assert card.metadata["public_adjusted_risk_output"] is False

    report = build_transparency_report(phase_output, validation_summary={"suite": "h4a"})
    assert report.accepted_non_blocking_caveats == []
    assert report.unresolved_blockers == []
    assert report.warning_counts_by_severity["blocker"] == 1
    assert report.model_card_summary.validation_summary["suite"] == "h4a"


def _profile(
    *,
    tissue="lung",
    risk=None,
    factor=None,
    flux_classes=None,
    exposure_high=None,
    exposure_elevated=None,
):
    interactions = None
    if risk is not None or factor is not None:
        interactions = SimpleNamespace(total_interaction_risk=risk, interaction_factor=factor)
    flux_profile = SimpleNamespace(elevated_or_high_risk_classes=flux_classes or [])
    exposure_profile = SimpleNamespace(
        high_risk_classes=exposure_high or [],
        elevated_risk_classes=exposure_elevated or [],
    )
    return PatientRiskProfile(
        tissue=tissue,
        genotypes={},
        lifestyle={},
        exposure_scenario="local_fixture",
        exposure_answers={},
        tissue_report=None,
        tissue_weight_count=None,
        top_tissue_genes=[],
        flux_profile=flux_profile,
        exposure_profile=exposure_profile,
        interactions=interactions,
    )


def test_unified_api_compare_patient_profiles_uses_public_objects_only():
    left = _profile(
        tissue="lung",
        risk=2.5,
        factor=1.2,
        flux_classes=["PAH", "Benzene"],
        exposure_high=["AflatoxinB1"],
        exposure_elevated=["PAH"],
    )
    right = _profile(
        tissue="liver",
        risk=1.0,
        factor=1.5,
        flux_classes=["Benzene"],
        exposure_elevated=["HeavyMetals"],
    )

    comparison = compare_patient_profiles(left, right, left_label="Left", right_label="Right")

    assert comparison.tissue == "lung vs liver"
    assert comparison.total_interaction_risk_delta == 1.5
    assert comparison.interaction_factor_delta == -0.3
    assert comparison.shared_high_risk_classes == ["Benzene"]
    assert comparison.left_only_high_risk_classes == ["AflatoxinB1", "PAH"]
    assert comparison.right_only_high_risk_classes == ["HeavyMetals"]
    assert comparison.more_concerning_profile == "Left"
    assert "Overall more concerning profile: Left" in comparison.summary

    class_only = compare_patient_profiles(
        _profile(flux_classes=["PAH", "Benzene"]),
        _profile(flux_classes=["PAH"]),
        left_label="Class-rich",
        right_label="Class-sparse",
    )
    assert class_only.more_concerning_profile == "Class-rich"
    assert class_only.total_interaction_risk_delta is None


def test_exporter_local_bundle_and_error_paths_are_deterministic(tmp_path):
    graph_data = tmp_path / "graph-data.js"
    graph_data.write_text("const GRAPH_DATA = {nodes: [], edges: []};\n", encoding="utf-8")

    template_with_external = tmp_path / "index.html"
    template_with_external.write_text(
        '<html><head><script src="./graph-data.js"></script></head><body>Omega Ω</body></html>',
        encoding="utf-8",
    )
    html = bundle_to_html_string(template_with_external, graph_data)
    assert '<script src="./graph-data.js"></script>' not in html
    assert "const GRAPH_DATA = {nodes: [], edges: []};" in html
    assert "Omega Ω" in html

    template_with_head = tmp_path / "head.html"
    template_with_head.write_text("<html><head></head><body>fallback</body></html>", encoding="utf-8")
    head_html = bundle_to_html_string(template_with_head, graph_data)
    assert head_html.index("const GRAPH_DATA") < head_html.index("</head>")

    with pytest.raises(ValueError, match="Could not locate graph data object"):
        parse_graph_data_text("no graph payload here")
    with pytest.raises(ValueError, match="Could not locate start of GRAPH_DATA object"):
        parse_graph_data_text("const GRAPH_DATA = null;")
    with pytest.raises(ValueError, match="Could not find the end of the GRAPH_DATA object"):
        parse_graph_data_text("const GRAPH_DATA = {nodes: [")
    with pytest.raises(FileNotFoundError):
        parse_graph_artifact(tmp_path / "missing.json")


def test_registry_loader_build_compare_and_validation_contracts(tmp_path):
    mapping_path = tmp_path / "nested" / "mapping.json"
    write_json_mapping(
        mapping_path,
        {
            "_metadata": {"schema_version": "fixture"},
            "_update_list": [],
            "entries": [{"biomarker": "β-HCH", "lifestyle_factor": "default"}],
        },
    )
    assert load_registry_document(mapping_path)["entries"][0]["biomarker"] == "β-HCH"
    assert mapping_path.read_text(encoding="utf-8").endswith("\n")
    assert "β-HCH" in mapping_path.read_text(encoding="utf-8")

    list_path = tmp_path / "list.json"
    list_path.write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="Expected JSON object"):
        load_json_mapping(list_path)
    with pytest.raises(ValueError, match="Unsupported registry file type"):
        load_registry_document(tmp_path / "mapping.txt")

    registry_path = tmp_path / "registry.json"
    write_json_mapping(
        registry_path,
        {
            "_metadata": {"source_family": "urinary", "registry_phase": "fixture"},
            "entries": [
                {
                    "biomarker": "1-OHP",
                    "lifestyle_factor": "smoking",
                    "matrix": "urine",
                    "reference_range": "0-1",
                    "reference_units": "ug/L",
                    "source_status": "measured",
                }
            ],
        },
    )
    manifest_path = tmp_path / "manifest.json"
    write_json_mapping(
        manifest_path,
        {
            "_metadata": {"schema_version": "manifest-fixture"},
            "source_documents": [
                {
                    "path": str(registry_path),
                    "source_family": "urinary",
                    "registry_phase": "fixture",
                    "registry_tier": "tier2",
                    "source_note": "local fixture",
                }
            ],
        },
    )

    metadata, entries = load_biomarker_mapping_manifest(manifest_path)
    assert metadata["source_registry_entry_count"] == 1
    assert entries[0]["biomarker"] == "1-OHP"

    built = build_biomarker_mapping_document(manifest_path)
    built_entry = built["entries"][0]
    assert built["_metadata"]["forward_update_compatible"] is True
    assert built_entry["entry_id"] == "1-OHP::smoking"
    assert built_entry["trace"]["source_registry"] == str(registry_path)
    assert built_entry["trace"]["registry_tier"] == "tier2"

    old_document = {
        "entries": [
            {"biomarker": "A", "lifestyle_factor": "smoking", "reference_units": "u", "trace": {"created_index": 0}},
            {"biomarker": "B", "lifestyle_factor": "default", "reference_units": "u", "trace": {"created_index": 1}},
        ]
    }
    new_document = {
        "entries": [
            {"biomarker": "B", "lifestyle_factor": "default", "reference_units": "ng", "trace": {"created_index": 99}},
            {"biomarker": "C", "lifestyle_factor": "default", "reference_units": "u", "trace": {"created_index": 2}},
        ]
    }
    comparison = compare_biomarker_mapping_documents(old_document, new_document)
    assert comparison == {
        "mapped_biomarkers_unchanged": False,
        "old_count": 2,
        "new_count": 2,
        "added_count": 1,
        "removed_count": 1,
        "changed_count": 1,
        "added": ["C::default"],
        "removed": ["A::smoking"],
        "changed": ["B::default"],
    }

    trace_only_change = compare_biomarker_mapping_documents(
        {"entries": [{"biomarker": "A", "lifestyle_factor": "default", "trace": {"created_index": 0}}]},
        {"entries": [{"biomarker": "A", "lifestyle_factor": "default", "trace": {"created_index": 100}}]},
    )
    assert trace_only_change["mapped_biomarkers_unchanged"] is True

    fixed, errors = validate_biomarker_mapping_document(
        {
            "_metadata": {},
            "_update_list": [],
            "entries": [{"biomarker": "1-OHP", "lifestyle_factor": "smoking", "trace": {"source_registry": "fixture"}}],
        },
        fix=True,
    )
    assert errors == []
    assert fixed["_metadata"]["forward_update_compatible"] is True
    assert fixed["entries"][0]["entry_id"] == "1-OHP::smoking"
    assert fixed["entries"][0]["trace"]["created_index"] == 0

    _, duplicate_errors = validate_biomarker_mapping_document(
        {
            "_metadata": {},
            "_update_list": [],
            "entries": [
                {"biomarker": "D", "lifestyle_factor": "default", "entry_id": "D::default", "trace": {"source_registry": "fixture"}},
                {"biomarker": "D", "lifestyle_factor": "default", "entry_id": "D::default", "trace": {"source_registry": "fixture"}},
            ],
        }
    )
    assert "duplicate biomarker entry_id: D::default" in duplicate_errors


def test_json_parameter_provider_preserves_local_provenance_and_scaffold_guardrails(tmp_path):
    data_dir = tmp_path / "parameters"
    data_dir.mkdir()
    (data_dir / "interaction_parameters.json").write_text(
        json.dumps(
            {
                "competitive_inhibition": {
                    "CYP1A1": {
                        "substrates": {
                            "BaP": {
                                "Km_uM": "2.5",
                                "Ki_uM": "not-numeric",
                                "Vmax_relative": "1.25",
                                "hill_coefficient": "bad",
                                "relative_priority": "2",
                                "assumed_ki": True,
                                "product": "BPDE",
                                "product_carcinogenic": True,
                                "reaction_role": "unknown",
                                "risk_direction_if_flux_decreases": "unknown",
                                "provenance_ref": "fixture-provenance",
                                "notes": "local fixture",
                            },
                            "skip_me": "not an object",
                        }
                    },
                    "not_an_enzyme": "skip",
                },
                "gsh_depletion": {
                    "consumers": {
                        "quinone": {
                            "enzyme": "GSTP1",
                            "substrate_class": "quinone",
                            "gsh_per_umol_substrate": "2.5",
                            "tissue": "liver",
                        },
                        "lung_case": {
                            "enzyme": "GSTM1",
                            "substrate_class": "epoxide",
                            "gsh_per_umol_substrate": "bad",
                            "tissue": "lung",
                        },
                    }
                },
                "enzyme_induction": {
                    "smoke": {
                        "CYP1A1": {
                            "fold_induction": "3.2",
                            "range_min": "1.1",
                            "range_max": "bad",
                            "mechanism": "AhR",
                            "tissue_specificity": "liver, lung",
                        },
                        "_notes": {"ignored": True},
                        "bad_rule": "skip",
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    (data_dir / "parameter_provenance.json").write_text(
        json.dumps(
            {
                "pairs": {
                    "CYP1A1": {
                        "bap": {
                            "ki_status": "curated",
                            "ki_reference": "local Ki review",
                            "notes": "case-insensitive provenance key",
                        }
                    }
                }
            }
        ),
        encoding="utf-8",
    )

    provider = JSONInteractionParameterProvider(data_dir)
    interactions = provider.get_competitive_interactions("cyp1a1")

    assert len(interactions) == 1
    interaction = interactions[0]
    assert interaction.enzyme == "CYP1A1"
    assert interaction.substrate == "BaP"
    assert interaction.reaction_role is ReactionRole.UNKNOWN
    assert interaction.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.UNKNOWN
    assert _warning_codes(interaction) == {"reaction_role_unknown", "risk_direction_unknown"}
    assert interaction.evidence.source == "local Ki review"
    assert interaction.evidence.grade is EvidenceGrade.CURATED
    assert interaction.evidence.provenance_ref == "parameter_provenance.json#pairs/CYP1A1/BaP"
    assert interaction.kinetic_parameters.km_uM == 2.5
    assert interaction.kinetic_parameters.ki_uM is None
    assert interaction.kinetic_parameters.hill_coefficient is None
    assert interaction.kinetic_parameters.relative_priority == 2
    assert interaction.kinetic_parameters.product_hazard == {"product": "BPDE", "product_carcinogenic": True}
    assert interaction.kinetic_parameters.uncertainty.confidence == "curated"

    assert provider.get_competitive_interactions("missing") == []
    assert provider.get_parameter_evidence("CYP1A1", "missing") is None
    assert provider.get_reactions_for_enzyme("CYP1A1")[0].product == "BPDE"
    assert provider.get_reactions_for_carcinogen("bap")[0].substrate == "BaP"

    liver_consumers = provider.get_gsh_consumers("liver")
    assert len(liver_consumers) == 1
    assert liver_consumers[0].gsh_per_umol_substrate == 2.5
    assert provider.get_gsh_consumers("kidney") == []

    lung_rules = provider.get_induction_rules("lung")
    assert len(lung_rules) == 1
    assert lung_rules[0].fold_induction == 3.2
    assert lung_rules[0].range_max is None
    assert provider.get_induction_rules("kidney") == []

    with pytest.raises(NotImplementedError, match="Phase 3 scaffold"):
        KGInteractionParameterProvider().get_competitive_interactions()
