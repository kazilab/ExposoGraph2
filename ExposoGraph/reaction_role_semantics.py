"""Reaction-role semantics for competitive inhibition.

Phase 4 keeps toxicological interpretation separate from the existing
competitive-flux calculation. The records in this module are explicit SME
annotations; they do not change interaction_engine public outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .interaction_schema import (
    AssumptionWarning,
    EvidenceGrade,
    EvidenceRecord,
    JsonDict,
    ReactionRole,
    ReleaseTarget,
    RiskDirectionIfFluxDecreases,
    RiskEndpoint,
    SMEReviewNote,
    SMEReviewStatus,
    SerializableRecord,
    TissueContext,
)


REACTION_ROLE_SME_SOURCE = "SME amendment captured for Phase 4 reaction-role semantics"


@dataclass
class ReactionRoleAnnotation(SerializableRecord):
    """Typed role annotation for one enzyme/substrate interpretation case."""

    enzyme: str
    substrate: str
    reaction_role: ReactionRole = ReactionRole.UNKNOWN
    risk_direction_if_flux_decreases: RiskDirectionIfFluxDecreases = (
        RiskDirectionIfFluxDecreases.UNKNOWN
    )
    evidence: EvidenceRecord | None = None
    warnings: list[AssumptionWarning] | None = None
    sme_notes: list[SMEReviewNote] | None = None
    tissue_context: TissueContext | None = None
    endpoint_context: str | None = None
    release_target: ReleaseTarget = ReleaseTarget.UNKNOWN
    review_status: SMEReviewStatus = SMEReviewStatus.UNKNOWN
    active: bool = True
    record_id: str | None = None
    metadata: JsonDict | None = None

    @property
    def is_context_specific(self) -> bool:
        return self.tissue_context is not None or self.endpoint_context is not None


@dataclass
class CompetitiveFluxInterpretation(SerializableRecord):
    """Flux-ratio facts passed into the standalone interpretation helper."""

    flux_ratio: float
    flux_decreased: bool
    annotation: ReactionRoleAnnotation


@dataclass
class CompetitiveRiskInterpretation(SerializableRecord):
    """Toxic-burden interpretation derived from an explicit role annotation."""

    flux: CompetitiveFluxInterpretation
    burden_multiplier: float
    risk_direction_if_flux_decreases: RiskDirectionIfFluxDecreases
    warnings: list[AssumptionWarning]
    sme_notes: list[SMEReviewNote]
    metadata: JsonDict | None = None


class ReactionRoleRegistry:
    """Deterministic registry for explicit reaction-role annotations."""

    def __init__(self, records: Iterable[ReactionRoleAnnotation]):
        self._records = tuple(records)

    @property
    def records(self) -> tuple[ReactionRoleAnnotation, ...]:
        return self._records

    def lookup(
        self,
        enzyme: str,
        substrate: str,
        *,
        tissue: str | None = None,
        endpoint: str | RiskEndpoint | None = None,
        release_target: ReleaseTarget | str | None = None,
        review_status: SMEReviewStatus | str | None = None,
    ) -> ReactionRoleAnnotation:
        """Return the best explicit annotation or an unknown neutral fallback."""

        release_filter = _coerce_release_target(release_target)
        review_filter = _coerce_review_status(review_status)
        matches: list[tuple[int, int, ReactionRoleAnnotation]] = []

        for index, record in enumerate(self._records):
            if not _matches_pair(enzyme, record.enzyme, record.metadata, "enzyme_aliases"):
                continue
            if not _matches_pair(substrate, record.substrate, record.metadata, "substrate_aliases"):
                continue
            if release_filter is not None and record.release_target is not release_filter:
                continue
            if review_filter is not None and record.review_status is not review_filter:
                continue
            if not _matches_context(record, tissue, endpoint):
                continue
            matches.append((_match_score(record, tissue, endpoint), -index, record))

        if not matches:
            return make_unknown_reaction_role_annotation(enzyme, substrate, tissue=tissue, endpoint=endpoint)

        return max(matches, key=lambda item: (item[0], item[1]))[2]


def interpret_competitive_flux_ratio(
    flux_ratio: float,
    annotation: ReactionRoleAnnotation,
) -> CompetitiveRiskInterpretation:
    """Interpret a competitive flux ratio under an explicit role annotation."""

    if flux_ratio <= 0:
        raise ValueError("flux_ratio must be greater than zero")

    warnings = list(annotation.warnings or [])
    sme_notes = list(annotation.sme_notes or [])
    flux = CompetitiveFluxInterpretation(
        flux_ratio=flux_ratio,
        flux_decreased=flux_ratio < 1.0,
        annotation=annotation,
    )
    direction = annotation.risk_direction_if_flux_decreases

    if not annotation.active:
        warnings.append(_warning("reaction_role_inactive", "Role annotation is not active curated behavior."))
        return _risk(flux, 1.0, direction, warnings, sme_notes)
    if annotation.reaction_role is ReactionRole.UNKNOWN:
        warnings.append(_warning("reaction_role_unknown", "Reaction role is unknown; using neutral burden."))
        return _risk(flux, 1.0, direction, warnings, sme_notes)
    if direction is RiskDirectionIfFluxDecreases.UNKNOWN:
        warnings.append(_warning("risk_direction_unknown", "Risk direction is unknown; using neutral burden."))
        return _risk(flux, 1.0, direction, warnings, sme_notes)
    if annotation.reaction_role is ReactionRole.PROBE_ONLY:
        warnings.append(_warning("probe_only_neutral", "Probe-only reactions are not interpreted as risk."))
        return _risk(flux, 1.0, direction, warnings, sme_notes)
    if direction is RiskDirectionIfFluxDecreases.NEUTRAL:
        warnings.append(_warning("risk_direction_neutral", "Neutral risk direction leaves burden unchanged."))
        return _risk(flux, 1.0, direction, warnings, sme_notes)
    if direction is RiskDirectionIfFluxDecreases.MIXED:
        warnings.append(_warning("risk_direction_mixed", "Mixed risk direction requires context-specific resolution."))
        return _risk(flux, 1.0, direction, warnings, sme_notes)
    if direction is RiskDirectionIfFluxDecreases.INCREASE:
        return _risk(flux, 1.0 / flux_ratio, direction, warnings, sme_notes)

    return _risk(flux, flux_ratio, direction, warnings, sme_notes)


def get_default_reaction_role_registry() -> ReactionRoleRegistry:
    """Build the default explicit Phase 4 registry."""

    return ReactionRoleRegistry(get_reaction_role_sme_records())


def get_reaction_role_sme_records() -> list[ReactionRoleAnnotation]:
    """Return SME-derived reaction-role records encoded as explicit annotations."""

    records = [
        _record(
            record_id="reaction_role_benzene_cyp2e1_v2_direction_adjustment",
            enzyme="CYP2E1",
            substrate="benzene",
            role=ReactionRole.DETOXIFICATION,
            direction=RiskDirectionIfFluxDecreases.INCREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            notes="Phase 4 encodes reaction-role-based direction interpretation relative to prior K-factor behavior.",
            warning_code="reaction_role_direction_adjustment",
            warning_message="Benzene CYP2E1 flux decrease is interpreted as higher burden in the scoped 2.0 SME record.",
            metadata={"prior_k_factor_direction_adjustment": True, "product_carcinogenic": False},
        ),
        _record(
            record_id="reaction_role_benzene_cyp2f1_bone_marrow_direction_adjustment",
            enzyme="CYP2F1",
            substrate="benzene",
            role=ReactionRole.DETOXIFICATION,
            direction=RiskDirectionIfFluxDecreases.INCREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            tissue="bone marrow",
            endpoint=RiskEndpoint.DNA_ADDUCT,
            endpoint_context="hematopoietic benzene toxicity",
            notes="Direction adjustment applies only in hematopoietic or bone-marrow benzene context.",
            warning_code="context_specific_direction_adjustment",
            warning_message="CYP2F1 benzene interpretation is context specific, not global.",
            metadata={"context_required": "hematopoietic_or_bone_marrow"},
        ),
        _record(
            record_id="reaction_role_benzene_cyp2f1_outside_context_unknown",
            enzyme="CYP2F1",
            substrate="benzene",
            role=ReactionRole.UNKNOWN,
            direction=RiskDirectionIfFluxDecreases.UNKNOWN,
            status=SMEReviewStatus.PENDING_TEAM_AGREEMENT,
            release=ReleaseTarget.V2_0,
            notes="Outside hematopoietic or bone-marrow benzene context the role remains unresolved.",
            warning_code="context_missing",
            warning_message="No global CYP2F1 benzene direction adjustment is inferred outside the approved context.",
        ),
        _record(
            record_id="reaction_role_benzene_cyp2a13_bone_marrow_direction_adjustment",
            enzyme="CYP2A13",
            substrate="benzene",
            role=ReactionRole.DETOXIFICATION,
            direction=RiskDirectionIfFluxDecreases.INCREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            tissue="bone marrow",
            endpoint=RiskEndpoint.DNA_ADDUCT,
            endpoint_context="hematopoietic benzene toxicity",
            notes="Direction adjustment applies only in hematopoietic or bone-marrow benzene context.",
            warning_code="context_specific_direction_adjustment",
            warning_message="CYP2A13 benzene interpretation is context specific, not global.",
            metadata={"context_required": "hematopoietic_or_bone_marrow"},
        ),
        _record(
            record_id="reaction_role_benzene_cyp2a13_outside_context_unknown",
            enzyme="CYP2A13",
            substrate="benzene",
            role=ReactionRole.UNKNOWN,
            direction=RiskDirectionIfFluxDecreases.UNKNOWN,
            status=SMEReviewStatus.PENDING_TEAM_AGREEMENT,
            release=ReleaseTarget.V2_0,
            notes="Outside hematopoietic or bone-marrow benzene context the role remains unresolved.",
            warning_code="context_missing",
            warning_message="No global CYP2A13 benzene direction adjustment is inferred outside the approved context.",
        ),
        _record(
            record_id="reaction_role_benzene_cyp1a1_default_direction",
            enzyme="CYP1A1",
            substrate="benzene",
            role=ReactionRole.BIOACTIVATION,
            direction=RiskDirectionIfFluxDecreases.DECREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            notes="Explicit default-direction benzene CYP1A1 record.",
            metadata={"default_direction_preserved": True},
        ),
        _record(
            record_id="reaction_role_ndma_cyp2e1_default_direction",
            enzyme="CYP2E1",
            substrate="NDMA",
            role=ReactionRole.BIOACTIVATION,
            direction=RiskDirectionIfFluxDecreases.DECREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            notes="Explicit default-direction NDMA CYP2E1 record.",
            metadata={"default_direction_preserved": True},
        ),
        _record(
            record_id="reaction_role_vinyl_chloride_cyp2e1_default_direction",
            enzyme="CYP2E1",
            substrate="vinyl_chloride",
            role=ReactionRole.BIOACTIVATION,
            direction=RiskDirectionIfFluxDecreases.DECREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            notes="Explicit default-direction vinyl chloride CYP2E1 record.",
            metadata={"default_direction_preserved": True, "substrate_aliases": ["vinyl chloride"]},
        ),
        _record(
            record_id="reaction_role_hca_cyp1a1_default_direction",
            enzyme="CYP1A1",
            substrate="HCA",
            role=ReactionRole.BIOACTIVATION,
            direction=RiskDirectionIfFluxDecreases.DECREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            notes="Explicit default-direction heterocyclic amine CYP1A1 class record.",
            metadata={"default_direction_preserved": True, "substrate_aliases": ["heterocyclic amine"]},
        ),
        _record(
            record_id="reaction_role_hca_cyp1a2_default_direction",
            enzyme="CYP1A2",
            substrate="HCA",
            role=ReactionRole.BIOACTIVATION,
            direction=RiskDirectionIfFluxDecreases.DECREASE,
            status=SMEReviewStatus.CURATED,
            release=ReleaseTarget.V2_0,
            notes="Explicit default-direction heterocyclic amine CYP1A2 class record.",
            metadata={"default_direction_preserved": True, "substrate_aliases": ["heterocyclic amine"]},
        ),
        _record(
            record_id="reaction_role_tce_cyp2e1_candidate_pending",
            enzyme="CYP2E1",
            substrate="trichloroethylene",
            role=ReactionRole.UNKNOWN,
            direction=RiskDirectionIfFluxDecreases.UNKNOWN,
            status=SMEReviewStatus.PENDING_TEAM_AGREEMENT,
            release=ReleaseTarget.V2_0,
            active=False,
            notes="Candidate 2.0 TCE CYP2E1 addition pending team agreement; not active behavior.",
            warning_code="pending_team_agreement",
            warning_message="TCE CYP2E1 role is represented as pending, not curated behavior.",
            metadata={"substrate_aliases": ["TCE"], "candidate_addition": True},
        ),
    ]

    records.extend(_deferred_detox_arm_records())
    return records


def make_unknown_reaction_role_annotation(
    enzyme: str,
    substrate: str,
    *,
    tissue: str | None = None,
    endpoint: str | RiskEndpoint | None = None,
) -> ReactionRoleAnnotation:
    """Create an explicit unknown fallback instead of inferring behavior."""

    endpoint_text = endpoint.value if isinstance(endpoint, RiskEndpoint) else endpoint
    return ReactionRoleAnnotation(
        enzyme=enzyme,
        substrate=substrate,
        reaction_role=ReactionRole.UNKNOWN,
        risk_direction_if_flux_decreases=RiskDirectionIfFluxDecreases.UNKNOWN,
        warnings=[
            _warning(
                "reaction_role_unmatched",
                "No explicit reaction-role record matched; using neutral unknown behavior.",
                field="reaction_role",
            )
        ],
        sme_notes=[
            SMEReviewNote(
                status=SMEReviewStatus.UNKNOWN,
                release_target=ReleaseTarget.UNKNOWN,
                tissue_context=tissue,
                endpoint_context=endpoint_text,
                notes="SME review required before assigning toxicological direction.",
            )
        ],
        review_status=SMEReviewStatus.UNKNOWN,
        release_target=ReleaseTarget.UNKNOWN,
        active=False,
        record_id="unknown_reaction_role",
    )


def _record(
    *,
    record_id: str,
    enzyme: str,
    substrate: str,
    role: ReactionRole,
    direction: RiskDirectionIfFluxDecreases,
    status: SMEReviewStatus,
    release: ReleaseTarget,
    notes: str,
    active: bool = True,
    tissue: str | None = None,
    endpoint: RiskEndpoint = RiskEndpoint.UNKNOWN,
    endpoint_context: str | None = None,
    warning_code: str | None = None,
    warning_message: str | None = None,
    metadata: JsonDict | None = None,
) -> ReactionRoleAnnotation:
    warnings = None
    if warning_code and warning_message:
        warnings = [_warning(warning_code, warning_message, review_status=status)]

    return ReactionRoleAnnotation(
        enzyme=enzyme,
        substrate=substrate,
        reaction_role=role,
        risk_direction_if_flux_decreases=direction,
        tissue_context=(
            TissueContext(tissue=tissue, species="human", endpoint=endpoint)
            if tissue or endpoint is not RiskEndpoint.UNKNOWN
            else None
        ),
        endpoint_context=endpoint_context,
        evidence=EvidenceRecord(
            source=REACTION_ROLE_SME_SOURCE,
            grade=EvidenceGrade.CURATED if status is SMEReviewStatus.CURATED else EvidenceGrade.PLACEHOLDER,
            confidence="SME_scoped" if status is SMEReviewStatus.CURATED else "pending",
            notes=notes,
            metadata=metadata,
        ),
        warnings=warnings,
        sme_notes=[
            SMEReviewNote(
                status=status,
                source=REACTION_ROLE_SME_SOURCE,
                release_target=release,
                tissue_context=tissue,
                endpoint_context=endpoint_context,
                team_agreement_status=status.value,
                notes=notes,
                assumption_warnings=warnings,
            )
        ],
        release_target=release,
        review_status=status,
        active=active,
        record_id=record_id,
        metadata=metadata,
    )


def _deferred_detox_arm_records() -> list[ReactionRoleAnnotation]:
    deferred = []
    for enzyme in ("NAT2", "ALDH", "GST", "EPHX1"):
        deferred.append(
            _record(
                record_id=f"reaction_role_{_normalize_key(enzyme)}_detox_arm_deferred_v3",
                enzyme=enzyme,
                substrate="broader_detox_arm_cases",
                role=ReactionRole.UNKNOWN,
                direction=RiskDirectionIfFluxDecreases.UNKNOWN,
                status=SMEReviewStatus.DEFERRED_3_0,
                release=ReleaseTarget.V3_0,
                active=False,
                notes=f"{enzyme} broader detox-arm cases are deferred to 3.0 and inactive in Phase 4.",
                warning_code="deferred_3_0",
                warning_message=f"{enzyme} detox-arm semantics require separate 3.0 approval.",
                metadata={"substrate_aliases": ["broader detox-arm cases", "detox arm"]},
            )
        )
    return deferred


def _matches_pair(query: str, canonical: str, metadata: JsonDict | None, alias_field: str) -> bool:
    query_key = _normalize_key(query)
    if query_key == _normalize_key(canonical):
        return True
    aliases = (metadata or {}).get(alias_field, [])
    return any(query_key == _normalize_key(alias) for alias in aliases)


def _matches_context(
    record: ReactionRoleAnnotation,
    tissue: str | None,
    endpoint: str | RiskEndpoint | None,
) -> bool:
    if not record.is_context_specific:
        return True
    query_parts = [_normalize_optional(tissue), _normalize_optional(endpoint)]
    if not any(query_parts):
        return False

    record_parts = [
        _normalize_optional(record.tissue_context.tissue if record.tissue_context else None),
        _normalize_optional(record.endpoint_context),
        _normalize_optional(record.tissue_context.endpoint if record.tissue_context else None),
    ]
    return any(_context_overlaps(query, record_value) for query in query_parts for record_value in record_parts)


def _match_score(
    record: ReactionRoleAnnotation,
    tissue: str | None,
    endpoint: str | RiskEndpoint | None,
) -> int:
    score = 0
    if record.active:
        score += 8
    if record.review_status is SMEReviewStatus.CURATED:
        score += 4
    if record.is_context_specific and (tissue is not None or endpoint is not None):
        score += 16
    return score


def _coerce_release_target(value: ReleaseTarget | str | None) -> ReleaseTarget | None:
    if value is None or isinstance(value, ReleaseTarget):
        return value
    try:
        return ReleaseTarget(str(value))
    except ValueError:
        return None


def _coerce_review_status(value: SMEReviewStatus | str | None) -> SMEReviewStatus | None:
    if value is None or isinstance(value, SMEReviewStatus):
        return value
    try:
        return SMEReviewStatus(str(value))
    except ValueError:
        return None


def _risk(
    flux: CompetitiveFluxInterpretation,
    burden_multiplier: float,
    direction: RiskDirectionIfFluxDecreases,
    warnings: list[AssumptionWarning],
    sme_notes: list[SMEReviewNote],
) -> CompetitiveRiskInterpretation:
    return CompetitiveRiskInterpretation(
        flux=flux,
        burden_multiplier=burden_multiplier,
        risk_direction_if_flux_decreases=direction,
        warnings=warnings,
        sme_notes=sme_notes,
    )


def _warning(
    code: str,
    message: str,
    *,
    field: str | None = "risk_direction_if_flux_decreases",
    review_status: SMEReviewStatus = SMEReviewStatus.UNKNOWN,
) -> AssumptionWarning:
    return AssumptionWarning(
        code=code,
        message=message,
        field=field,
        review_status=review_status,
    )


def _context_overlaps(query: str | None, record_value: str | None) -> bool:
    if not query or not record_value:
        return False
    return query in record_value or record_value in query


def _normalize_optional(value: str | RiskEndpoint | None) -> str | None:
    if value is None:
        return None
    if isinstance(value, RiskEndpoint):
        value = value.value
    return _normalize_key(value)


def _normalize_key(value: object) -> str:
    return "".join(character.lower() for character in str(value) if character.isalnum())
