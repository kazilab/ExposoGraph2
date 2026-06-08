from ExposoGraph.interaction_schema import (
    AssumptionWarning,
    CompetitiveInteraction,
    EvidenceGrade,
    EvidenceRecord,
    KineticParameterSet,
    MetabolicReaction,
    ReactionRole,
    ReleaseTarget,
    RiskDirectionIfFluxDecreases,
    RiskEndpoint,
    SMEReviewNote,
    SMEReviewStatus,
    TissueContext,
)


def test_phase3_enums_include_mechanism_placeholders_without_behavior():
    assert ReactionRole.BIOACTIVATION.value == "bioactivation"
    assert ReactionRole.DETOXIFICATION.value == "detoxification"
    assert ReactionRole.UNKNOWN.value == "unknown"
    assert RiskDirectionIfFluxDecreases.INCREASE.value == "increase"
    assert RiskDirectionIfFluxDecreases.UNKNOWN.value == "unknown"
    assert RiskEndpoint.DNA_ADDUCT.value == "DNA_adduct"
    assert SMEReviewStatus.DEFERRED_3_0.value == "deferred_3_0"


def test_competitive_interaction_defaults_leave_mechanism_unresolved():
    interaction = CompetitiveInteraction(
        enzyme="CYP2E1",
        substrate="benzene",
        kinetic_parameters=KineticParameterSet(
            km_uM=10.0,
            product="phenol",
            product_hazard={"product": "phenol", "product_carcinogenic": False},
        ),
    )

    assert interaction.reaction_role is ReactionRole.UNKNOWN
    assert interaction.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.UNKNOWN
    assert interaction.kinetic_parameters.product_hazard["product_carcinogenic"] is False


def test_sme_notes_can_be_carried_without_driving_sign_logic():
    interaction = CompetitiveInteraction(
        enzyme="CYP1A2",
        substrate="candidate substrate",
        sme_notes=[
            SMEReviewNote(
                status=SMEReviewStatus.CANDIDATE,
                release_target=ReleaseTarget.V3_0,
                notes="Spyros draft input captured as SME note only.",
            )
        ],
    )

    assert interaction.sme_notes[0].status is SMEReviewStatus.CANDIDATE
    assert interaction.sme_notes[0].release_target is ReleaseTarget.V3_0
    assert interaction.reaction_role is ReactionRole.UNKNOWN
    assert interaction.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.UNKNOWN


def test_schema_serialization_uses_enum_values_and_preserves_context():
    reaction = MetabolicReaction(
        enzyme="CYP1A1",
        substrate="benzo[a]pyrene",
        product="reactive metabolite",
        reaction_role=ReactionRole.UNKNOWN,
        risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.UNKNOWN,
        risk_endpoints=[RiskEndpoint.DNA_ADDUCT],
        tissue_context=TissueContext(
            tissue="liver",
            species="human",
            endpoint=RiskEndpoint.DNA_ADDUCT,
        ),
        evidence=[
            EvidenceRecord(
                source="local metadata",
                grade=EvidenceGrade.INFERRED_FROM_LOCAL_METADATA,
                confidence="medium",
            )
        ],
        warnings=[
            AssumptionWarning(
                code="reaction_role_unknown",
                message="No curated role was supplied.",
                field="reaction_role",
            )
        ],
        sme_notes=[SMEReviewNote(status=SMEReviewStatus.PENDING_TEAM_AGREEMENT)],
    )

    payload = reaction.to_dict()

    assert payload["reaction_role"] == "unknown"
    assert payload["risk_direction_if_flux_decreases"] == "unknown"
    assert payload["risk_endpoints"] == ["DNA_adduct"]
    assert payload["tissue_context"]["endpoint"] == "DNA_adduct"
    assert payload["evidence"][0]["grade"] == "inferred_from_local_metadata"
    assert payload["warnings"][0]["code"] == "reaction_role_unknown"
    assert payload["sme_notes"][0]["status"] == "pending_team_agreement"


def test_metabolic_reaction_serializes_single_evidence_record():
    reaction = MetabolicReaction(
        enzyme="CYP2E1",
        substrate="benzene",
        evidence=EvidenceRecord(
            source="parameter_provenance.json",
            grade=EvidenceGrade.INFERRED_FROM_LOCAL_METADATA,
            provenance_ref="parameter_provenance.json#pairs/CYP2E1/benzene",
        ),
    )

    payload = reaction.to_dict()

    assert payload["evidence"]["source"] == "parameter_provenance.json"
    assert payload["evidence"]["grade"] == "inferred_from_local_metadata"
    assert payload["reaction_role"] == "unknown"
    assert payload["risk_direction_if_flux_decreases"] == "unknown"


def test_metabolic_reaction_serializes_multiple_evidence_records():
    reaction = MetabolicReaction(
        enzyme="CYP2E1",
        substrate="benzene",
        evidence=[
            EvidenceRecord(source="interaction_parameters.json", grade=EvidenceGrade.PLACEHOLDER),
            EvidenceRecord(source="parameter_provenance.json", grade=EvidenceGrade.INFERRED_FROM_LOCAL_METADATA),
        ],
    )

    payload = reaction.to_dict()

    assert [item["source"] for item in payload["evidence"]] == [
        "interaction_parameters.json",
        "parameter_provenance.json",
    ]
    assert payload["evidence"][0]["grade"] == "placeholder"
    assert payload["evidence"][1]["grade"] == "inferred_from_local_metadata"
    assert payload["reaction_role"] == "unknown"
    assert payload["risk_direction_if_flux_decreases"] == "unknown"
