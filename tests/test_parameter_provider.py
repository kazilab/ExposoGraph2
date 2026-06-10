import pytest

from ExposoGraph.interaction_schema import (
    CompetitiveInteraction,
    ReactionRole,
    ReleaseTarget,
    RiskDirectionIfFluxDecreases,
    SMEReviewNote,
    SMEReviewStatus,
)
from ExposoGraph.parameter_provider import (
    JSONInteractionParameterProvider,
    KGInteractionParameterProvider,
)


def test_json_provider_loads_competitive_interactions_from_local_files():
    provider = JSONInteractionParameterProvider()

    all_interactions = provider.get_competitive_interactions()
    cyp2e1_interactions = provider.get_competitive_interactions("CYP2E1")

    assert all_interactions
    assert cyp2e1_interactions
    assert all(isinstance(item, CompetitiveInteraction) for item in all_interactions)


def test_json_provider_preserves_local_kinetics_and_marks_missing_mechanism_unknown():
    provider = JSONInteractionParameterProvider()

    benzene = next(
        item
        for item in provider.get_competitive_interactions("CYP2E1")
        if item.substrate == "benzene"
    )

    assert benzene.enzyme == "CYP2E1"
    assert benzene.reaction_role is ReactionRole.UNKNOWN
    assert benzene.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.UNKNOWN
    assert benzene.kinetic_parameters.km_uM is not None
    assert benzene.kinetic_parameters.product == "benzene_oxide"
    assert "product_carcinogenic" in benzene.kinetic_parameters.product_hazard
    assert any(warning.code == "reaction_role_unknown" for warning in benzene.warnings)
    assert any(
        warning.code == "risk_direction_unknown"
        for warning in benzene.warnings
    )


def test_json_provider_returns_provenance_as_evidence_record():
    provider = JSONInteractionParameterProvider()

    evidence = provider.get_parameter_evidence("CYP2E1", "benzene")

    assert evidence is not None
    assert evidence.provenance_ref == "parameter_provenance.json#pairs/CYP2E1/benzene"
    assert evidence.grade.value in {"curated", "inferred_from_local_metadata"}
    assert evidence.source is not None


def test_reaction_lookup_by_carcinogen_uses_existing_json_without_tissue_policy():
    provider = JSONInteractionParameterProvider()

    reactions = provider.get_reactions_for_carcinogen("benzene", tissue="liver")

    assert reactions
    assert all(item.substrate == "benzene" for item in reactions)
    assert all(item.reaction_role is ReactionRole.UNKNOWN for item in reactions)


def test_gsh_and_induction_provider_surfaces_are_read_only_scaffolds():
    provider = JSONInteractionParameterProvider()

    consumers = provider.get_gsh_consumers()
    rules = provider.get_induction_rules()

    assert consumers
    assert rules
    assert all(consumer.name for consumer in consumers)
    assert all(rule.enzyme for rule in rules)


def test_kg_provider_is_explicitly_deferred_scaffold():
    provider = KGInteractionParameterProvider()

    with pytest.raises(NotImplementedError, match="Phase 3 scaffold"):
        provider.get_competitive_interactions()


def test_spyros_draft_fields_are_representable_as_review_notes_only():
    interaction = CompetitiveInteraction(
        enzyme="CYP3A4",
        substrate="future SME substrate",
        sme_notes=[
            SMEReviewNote(
                status=SMEReviewStatus.CANDIDATE,
                release_target=ReleaseTarget.V3_0,
                notes="Candidate mechanism note retained for later agreement.",
            )
        ],
    )

    assert interaction.sme_notes[0].status is SMEReviewStatus.CANDIDATE
    assert interaction.sme_notes[0].release_target is ReleaseTarget.V3_0
    assert interaction.reaction_role is ReactionRole.UNKNOWN
    assert interaction.risk_direction_if_flux_decreases is RiskDirectionIfFluxDecreases.UNKNOWN


