"""Standalone transparency, provenance, SME-review, and model-card layer.

Phase 10 makes assumptions, caveats, review status, release boundaries, and
model limits visible. It consumes already-created records and warnings; it does
not compute mechanism outputs or integrate with public engine behavior.
"""

from __future__ import annotations

from dataclasses import dataclass, field, is_dataclass
from typing import Any, Iterable, Mapping

from .interaction_schema import (
    AssumptionWarning,
    EvidenceGrade,
    EvidenceRecord,
    JsonDict,
    ReleaseTarget,
    SMEReviewNote,
    SMEReviewStatus,
    SerializableRecord,
    ValueEnum,
)
from .reaction_role_semantics import get_reaction_role_sme_records


class AssumptionCategory(ValueEnum):
    """Phase 10 transparency categories."""

    REACTION_ROLE = "reaction_role"
    KINETIC_PARAMETER = "kinetic_parameter"
    IC50_CONVERSION = "ic50_conversion"
    AFFINITY_FALLBACK = "affinity_fallback"
    GSH_TISSUE_PRESET = "gsh_tissue_preset"
    GSH_MODEL_BOUNDARY = "gsh_model_boundary"
    EFFECTIVE_BURDEN_BOUNDARY = "effective_burden_boundary"
    SHAPLEY_ATTRIBUTION_BOUNDARY = "shapley_attribution_boundary"
    RELEASE_DEFERRAL = "release_deferral"
    UNKNOWN_BIOLOGY = "unknown_biology"
    PRODUCT_CARCINOGENIC_GUARDRAIL = "product_carcinogenic_guardrail"
    PUBLIC_API_BOUNDARY = "public_api_boundary"
    GENERAL = "general"


class ReviewSeverity(ValueEnum):
    """Review severity for transparency and SME queue records."""

    INFO = "info"
    WARNING = "warning"
    REVIEW_REQUIRED = "review_required"
    BLOCKER = "blocker"


@dataclass
class TransparencyRecord(SerializableRecord):
    """Serializable assumption, warning, provenance, or caveat record."""

    category: AssumptionCategory
    message: str
    severity: ReviewSeverity = ReviewSeverity.WARNING
    source_phase: str = "phase_10"
    review_status: SMEReviewStatus = SMEReviewStatus.UNKNOWN
    release_target: ReleaseTarget = ReleaseTarget.UNKNOWN
    code: str | None = None
    evidence: EvidenceRecord | list[EvidenceRecord] | None = None
    metadata: JsonDict | None = None


@dataclass
class SMEReviewQueueItem(SerializableRecord):
    """Serializable item requiring SME or team review."""

    category: AssumptionCategory
    message: str
    severity: ReviewSeverity
    source_phase: str
    review_status: SMEReviewStatus
    release_target: ReleaseTarget
    reason: str
    code: str | None = None
    evidence: EvidenceRecord | list[EvidenceRecord] | None = None
    metadata: JsonDict | None = None


@dataclass
class DeferralRecord(SerializableRecord):
    """Deferred item with an explicit release target."""

    subject: str
    category: AssumptionCategory
    message: str
    source_phase: str
    release_target: ReleaseTarget
    review_status: SMEReviewStatus = SMEReviewStatus.DEFERRED_3_0
    evidence: EvidenceRecord | None = None
    metadata: JsonDict | None = None


@dataclass
class ModelBoundaryStatement(SerializableRecord):
    """Deterministic model-boundary statement for a phase."""

    phase: str
    category: AssumptionCategory
    statement: str
    context: str
    severity: ReviewSeverity = ReviewSeverity.INFO
    release_target: ReleaseTarget = ReleaseTarget.V2_0
    evidence: EvidenceRecord | None = None
    metadata: JsonDict | None = None


@dataclass
class ModelCardSummary(SerializableRecord):
    """Compact model-card-like summary for local review packages."""

    intended_context_of_use: str
    out_of_scope_use: list[str]
    major_assumptions: list[str]
    known_limitations: list[str]
    sme_review_items: list[SMEReviewQueueItem]
    deferred_items: list[DeferralRecord]
    provenance_evidence_summary: JsonDict
    validation_summary: JsonDict
    warning_counts_by_category: JsonDict
    release_target_summary: JsonDict
    model_boundaries: list[ModelBoundaryStatement]
    metadata: JsonDict | None = None


@dataclass
class TransparencyReport(SerializableRecord):
    """Complete Phase 10 transparency report payload."""

    records: list[TransparencyRecord]
    sme_review_queue: list[SMEReviewQueueItem]
    deferrals: list[DeferralRecord]
    model_boundaries: list[ModelBoundaryStatement]
    model_card_summary: ModelCardSummary
    warning_counts_by_category: JsonDict
    warning_counts_by_severity: JsonDict
    unresolved_unknown_count: int
    accepted_non_blocking_caveats: list[str] = field(default_factory=list)
    unresolved_blockers: list[str] = field(default_factory=list)
    metadata: JsonDict | None = None


REVIEW_QUEUE_STATUSES = {
    SMEReviewStatus.UNKNOWN,
    SMEReviewStatus.PENDING_TEAM_AGREEMENT,
    SMEReviewStatus.CANDIDATE,
    SMEReviewStatus.DEFERRED_3_0,
}

INTENDED_CONTEXT_OF_USE = (
    "Mechanism-resolved, PBK/PBPK-compatible, semi-mechanistic, "
    "relative-risk interaction index for research-oriented comparison, "
    "prioritization, and hypothesis generation."
)

OUT_OF_SCOPE_USE = [
    "clinical risk prediction",
    "validated absolute risk estimation",
    "fully validated PBPK/ODE simulation",
    "proof of biological causality from attribution outputs",
    "resolution of unknown biology without curated evidence or SME agreement",
]

MAJOR_ASSUMPTIONS = [
    "Unknown biology remains unknown, neutral, warned, and queued for review.",
    "Exact curated Ki is preferred; Km proxy use is low-confidence and visible.",
    "IC50 conversion and affinity/RDKit/ECFP4/Tanimoto fallback are guarded and not claimed when unavailable.",
    "GSH tissue presets are local 2.0 relative-capacity placeholders unless curated later.",
    "Effective burden remains an internal semi-mechanistic relative-burden ratio, not public adjusted risk.",
    "Shapley and mechanism interactions are model-output attribution only, not causal proof.",
]

KNOWN_LIMITATIONS = [
    "No clinical prediction claim is made.",
    "No fully validated PBPK, ODE, GSH-GSSG, or Nrf2 claim is made.",
    "No causal proof claim is made for Shapley attribution.",
    "Unknown reaction roles, endpoint signs, and biological mechanisms are not resolved by this layer.",
    "Release 2.0 boundaries keep deferred NAT2, ALDH, GST, EPHX1, and optional affinity fallback work inactive.",
]

MODULE5_MECHANISM_MODEL_VERSION = "module5_mechanism_resolved_v2"
MODULE5_RISK_CALCULATION_BASIS = "mechanism_resolved_adjusted_relative_risk"
MODULE5_SYNERGY_DECOMPOSITION_BASIS = "eight_state_shapley"
MODULE5_INHIBITION_PARAMETER_POLICY = (
    "curated_Ki_preferred; Km_proxy_low_confidence_warned; IC50_conversion_guarded"
)
MODULE5_REACTION_ROLE_POLICY = (
    "explicit_role_and_risk_direction_required; unknown_pending_or_deferred_biology_neutral_warned"
)
MODULE5_DIAGNOSTIC_OUTPUT_POLICY = (
    "mechanism_resolved_adjusted_risk_authoritative; biological_output_marks_selected_or_diagnostic_effects"
)


def build_module5_model_card_summary(
    *,
    mechanism_model_version: str = MODULE5_MECHANISM_MODEL_VERSION,
    risk_calculation_basis: str = MODULE5_RISK_CALCULATION_BASIS,
    gsh_model_version: str | None = None,
    inhibition_parameter_policy: str = MODULE5_INHIBITION_PARAMETER_POLICY,
    reaction_role_policy: str = MODULE5_REACTION_ROLE_POLICY,
    synergy_decomposition_basis: str = MODULE5_SYNERGY_DECOMPOSITION_BASIS,
    diagnostic_output_policy: str = MODULE5_DIAGNOSTIC_OUTPUT_POLICY,
    review_required_count: int = 0,
    warning_count: int = 0,
    unresolved_or_deferred_count: int = 0,
    evidence_summary_available: bool = True,
    detailed_records_location: Mapping[str, Any] | None = None,
    ki_resolver_statuses: Iterable[str] | None = None,
    reaction_role_review_statuses: Iterable[str] | None = None,
    model_boundary_warnings_present: bool = False,
) -> JsonDict:
    """Build the compact public Module 5 model-card summary."""

    return {
        "mechanism_model_version": mechanism_model_version,
        "risk_calculation_basis": risk_calculation_basis,
        "gsh_model_version": gsh_model_version,
        "inhibition_parameter_policy": inhibition_parameter_policy,
        "reaction_role_policy": reaction_role_policy,
        "synergy_decomposition_basis": synergy_decomposition_basis,
        "diagnostic_output_policy": diagnostic_output_policy,
        "review_required_count": int(review_required_count),
        "warning_count": int(warning_count),
        "unresolved_or_deferred_count": int(unresolved_or_deferred_count),
        "evidence_summary_available": bool(evidence_summary_available),
        "detailed_records_location": dict(detailed_records_location or {}),
        "ki_resolver_statuses": sorted({str(item) for item in (ki_resolver_statuses or []) if item}),
        "reaction_role_review_statuses": sorted(
            {str(item) for item in (reaction_role_review_statuses or []) if item}
        ),
        "model_boundary_warnings_present": bool(model_boundary_warnings_present),
    }


def collect_assumption_warnings(*phase_outputs: Any) -> list[TransparencyRecord]:
    """Collect warnings, SME notes, and transparency records from supplied outputs."""

    records: list[TransparencyRecord] = []
    for item in _iter_items(phase_outputs):
        _collect_from_item(item, records, seen=set())
    return _sort_records(_dedupe_records(records))


def build_sme_review_queue(
    *phase_outputs: Any,
    include_phase4_registry: bool = True,
    include_model_boundary_caveats: bool = True,
) -> list[SMEReviewQueueItem]:
    """Build a deterministic SME review queue from warnings and fixed phase caveats."""

    records = collect_assumption_warnings(*phase_outputs)
    items = [_queue_item_from_record(record) for record in records if _needs_review(record)]

    if include_phase4_registry:
        items.extend(_phase4_registry_queue_items())
    if include_model_boundary_caveats:
        items.extend(_default_boundary_queue_items())

    return _sort_queue_items(_dedupe_queue_items(items))


def summarize_model_boundaries() -> list[ModelBoundaryStatement]:
    """Return deterministic boundary statements for accepted Module 5 phases."""

    statements = [
        ModelBoundaryStatement(
            phase="phase_4_reaction_role_semantics",
            category=AssumptionCategory.PRODUCT_CARCINOGENIC_GUARDRAIL,
            statement="product_carcinogenic is evidence and metadata only, not the final sign rule.",
            context="Reaction-role semantics require explicit role and risk direction evidence.",
            metadata={"public_risk_behavior_changed": False},
        ),
        ModelBoundaryStatement(
            phase="phase_4_reaction_role_semantics",
            category=AssumptionCategory.UNKNOWN_BIOLOGY,
            statement="Unknown roles remain neutral and warned until curated evidence or SME agreement resolves them.",
            context="Unknown biology stays unknown.",
        ),
        ModelBoundaryStatement(
            phase="phase_5_kinetic_resolution",
            category=AssumptionCategory.KINETIC_PARAMETER,
            statement="Curated Ki is preferred; Km proxy is low-confidence and must be warned.",
            context="Ki, Km, IC50, docking scores, and affinity scores are not interchangeable.",
        ),
        ModelBoundaryStatement(
            phase="phase_5_kinetic_resolution",
            category=AssumptionCategory.IC50_CONVERSION,
            statement="IC50 conversion is guarded and is not claimed when assay context is insufficient.",
            context="No hidden IC50-to-Ki conversion is performed by the transparency layer.",
        ),
        ModelBoundaryStatement(
            phase="phase_5_kinetic_resolution",
            category=AssumptionCategory.AFFINITY_FALLBACK,
            statement="Affinity/RDKit/ECFP4/Tanimoto conversion is not implemented or claimed active.",
            context="Future affinity fallback remains optional, bounded, and inactive until validated and approved.",
            release_target=ReleaseTarget.FUTURE,
        ),
        ModelBoundaryStatement(
            phase="phase_5_kinetic_resolution",
            category=AssumptionCategory.KINETIC_PARAMETER,
            statement="Module 3 Km remains static for 2.0.",
            context="No dynamic Module 3 Km update is introduced.",
        ),
        ModelBoundaryStatement(
            phase="phase_6_endpoint_toxic_flux",
            category=AssumptionCategory.PUBLIC_API_BOUNDARY,
            statement="Endpoint toxic flux accepts precomputed flux ratios and produces no public adjusted-risk output.",
            context="Flux calculation remains separate from toxicological interpretation.",
            metadata={"adjusted_risk_output": False},
        ),
        ModelBoundaryStatement(
            phase="phase_7_gsh_redox_capacity",
            category=AssumptionCategory.GSH_MODEL_BOUNDARY,
            statement="GSH redox capacity is quasi-steady relative capacity, not PBPK/ODE/GSH-GSSG/Nrf2 validated.",
            context="The GSH model must not overclaim biological validation.",
            metadata={"validated_pbpk_ode_gsh_gssg_nrf2": False},
        ),
        ModelBoundaryStatement(
            phase="phase_7_gsh_redox_capacity",
            category=AssumptionCategory.GSH_TISSUE_PRESET,
            statement="GSH tissue presets are local 2.0 presets unless curated later.",
            context="Tissue constants are not invented or externally validated by this layer.",
        ),
        ModelBoundaryStatement(
            phase="phase_8_effective_burden",
            category=AssumptionCategory.EFFECTIVE_BURDEN_BOUNDARY,
            statement="Effective burden is a standalone semi-mechanistic relative-burden ratio, not public adjusted risk.",
            context="No clinical risk claim or public adjusted-risk output change is made.",
            metadata={"clinical_risk_claim": False, "adjusted_risk_output": False},
        ),
        ModelBoundaryStatement(
            phase="phase_9_shapley_interactions",
            category=AssumptionCategory.SHAPLEY_ATTRIBUTION_BOUNDARY,
            statement="Shapley results are model-output attribution only, not proof of biological causality.",
            context="Attribution explains configured model outputs under supplied state values.",
            metadata={"biological_causality_claim": False},
        ),
        ModelBoundaryStatement(
            phase="phase_9_shapley_interactions",
            category=AssumptionCategory.SHAPLEY_ATTRIBUTION_BOUNDARY,
            statement="Complete eight-state inputs have no unexplained residual beyond tolerance.",
            context="Residual fields are numerical verification, not extra biology.",
        ),
    ]
    return sorted(statements, key=lambda item: (item.phase, item.category.value, item.statement))


def build_model_card_summary(
    *phase_outputs: Any,
    validation_summary: Mapping[str, Any] | None = None,
) -> ModelCardSummary:
    """Build a compact deterministic model-card-like summary."""

    records = collect_assumption_warnings(*phase_outputs)
    queue = build_sme_review_queue(*phase_outputs)
    deferrals = build_deferral_records()
    boundaries = summarize_model_boundaries()
    warning_counts = _count_by_category(records)

    return ModelCardSummary(
        intended_context_of_use=INTENDED_CONTEXT_OF_USE,
        out_of_scope_use=list(OUT_OF_SCOPE_USE),
        major_assumptions=list(MAJOR_ASSUMPTIONS),
        known_limitations=list(KNOWN_LIMITATIONS),
        sme_review_items=queue,
        deferred_items=deferrals,
        provenance_evidence_summary=_summarize_evidence(records, boundaries),
        validation_summary=dict(validation_summary or {}),
        warning_counts_by_category=warning_counts,
        release_target_summary=_release_target_summary(queue, deferrals),
        model_boundaries=boundaries,
        metadata={
            "model_card_type": "phase_10_transparency_summary",
            "clinical_risk_prediction_claim": False,
            "fully_validated_pbpk_ode_claim": False,
            "shapley_causality_proof_claim": False,
            "unknown_biology_resolved_claim": False,
            "public_adjusted_risk_output": False,
            "github_or_release_packaging_behavior": False,
        },
    )


def build_transparency_report(
    *phase_outputs: Any,
    validation_summary: Mapping[str, Any] | None = None,
) -> TransparencyReport:
    """Build the full Phase 10 report payload from supplied phase outputs."""

    records = collect_assumption_warnings(*phase_outputs)
    queue = build_sme_review_queue(*phase_outputs)
    deferrals = build_deferral_records()
    boundaries = summarize_model_boundaries()
    card = build_model_card_summary(*phase_outputs, validation_summary=validation_summary)

    return TransparencyReport(
        records=records,
        sme_review_queue=queue,
        deferrals=deferrals,
        model_boundaries=boundaries,
        model_card_summary=card,
        warning_counts_by_category=_count_by_category(records),
        warning_counts_by_severity=_count_by_severity(records),
        unresolved_unknown_count=sum(1 for item in queue if item.review_status is SMEReviewStatus.UNKNOWN),
        accepted_non_blocking_caveats=[],
        unresolved_blockers=[],
        metadata={
            "phase": "phase_10_transparency_reporting",
            "transparency_model_card_layer_only": True,
            "scientific_calculations_changed": False,
            "engine_integration": False,
            "interaction_engine_integration": False,
            "unified_api_integration": False,
            "public_adjusted_risk_output": False,
            "release_packaging_behavior": False,
            "github_behavior": False,
            "phase11_behavior": False,
        },
    )


def transparency_report_from_phase_outputs(
    phase_outputs: Mapping[str, Any] | Iterable[Any] | None = None,
    *,
    validation_summary: Mapping[str, Any] | None = None,
    **named_phase_outputs: Any,
) -> TransparencyReport:
    """Build a transparency report from caller-supplied phase outputs."""

    supplied: list[Any] = []
    if isinstance(phase_outputs, Mapping):
        supplied.extend(phase_outputs.values())
    elif phase_outputs is not None:
        supplied.append(phase_outputs)
    supplied.extend(named_phase_outputs.values())
    return build_transparency_report(*supplied, validation_summary=validation_summary)


def build_deferral_records(include_phase4_registry: bool = True) -> list[DeferralRecord]:
    """Return visible 3.0/future release deferrals."""

    deferrals: list[DeferralRecord] = []
    if include_phase4_registry:
        for record in get_reaction_role_sme_records():
            if record.review_status is SMEReviewStatus.DEFERRED_3_0:
                subject = _deferred_subject(record.enzyme)
                deferrals.append(
                    DeferralRecord(
                        subject=subject,
                        category=AssumptionCategory.RELEASE_DEFERRAL,
                        message=_first_sme_note(record) or f"{record.enzyme} broader detox-arm semantics are deferred to 3.0.",
                        source_phase="phase_4_reaction_role_semantics",
                        release_target=record.release_target,
                        review_status=record.review_status,
                        evidence=record.evidence,
                        metadata={
                            "record_id": record.record_id,
                            "enzyme": record.enzyme,
                            "active": record.active,
                            "reaction_role": record.reaction_role.value,
                        },
                    )
                )

    deferrals.append(
        DeferralRecord(
            subject="RDKit/ECFP4/Tanimoto affinity fallback",
            category=AssumptionCategory.AFFINITY_FALLBACK,
            message="Future RDKit/ECFP4/Tanimoto affinity fallback remains optional, bounded, and inactive unless validated and approved.",
            source_phase="phase_5_kinetic_resolution",
            release_target=ReleaseTarget.FUTURE,
            review_status=SMEReviewStatus.DEFERRED_3_0,
            evidence=EvidenceRecord(
                source="Phase 5 kinetic resolver boundary",
                grade=EvidenceGrade.PLACEHOLDER,
                confidence="future_optional",
                notes="Affinity fallback unavailable in 2.0.",
            ),
            metadata={"active_in_2_0": False, "conversion_claimed": False},
        )
    )
    return _sort_deferrals(_dedupe_deferrals(deferrals))


def _iter_items(values: Iterable[Any]) -> Iterable[Any]:
    for value in values:
        if value is None:
            continue
        if isinstance(value, (list, tuple, set, frozenset)):
            for nested in value:
                yield nested
        else:
            yield value


def _collect_from_item(item: Any, records: list[TransparencyRecord], seen: set[int]) -> None:
    if item is None:
        return
    item_id = id(item)
    if item_id in seen:
        return
    seen.add(item_id)

    if isinstance(item, TransparencyRecord):
        records.append(item)
        return
    if _looks_like_warning(item):
        records.append(_record_from_warning(item))
        return
    if isinstance(item, SMEReviewNote):
        records.extend(_records_from_sme_note(item))
        return
    if isinstance(item, Mapping):
        if _mapping_looks_like_warning(item):
            records.append(_record_from_warning_mapping(item))
        for key in ("warnings", "assumption_warnings", "sme_notes", "records"):
            if key in item:
                for nested in _iter_items((item[key],)):
                    _collect_from_item(nested, records, seen)
        return

    for attr in ("warnings", "assumption_warnings", "sme_notes", "records"):
        nested = getattr(item, attr, None)
        if nested is not None:
            for value in _iter_items((nested,)):
                _collect_from_item(value, records, seen)

    metadata = getattr(item, "metadata", None)
    if isinstance(metadata, Mapping):
        _metadata_boundary_records(metadata, records)

    if is_dataclass(item):
        note = getattr(item, "review_status", None)
        if note in REVIEW_QUEUE_STATUSES:
            message = getattr(item, "record_id", None) or getattr(item, "notes", None) or item.__class__.__name__
            category = _category_from_code_message(str(message), str(message))
            records.append(
                TransparencyRecord(
                    category=category,
                    message=str(message),
                    severity=ReviewSeverity.REVIEW_REQUIRED,
                    source_phase=_source_phase_from_category(category),
                    review_status=note,
                    release_target=getattr(item, "release_target", ReleaseTarget.UNKNOWN),
                    code=str(message),
                    evidence=getattr(item, "evidence", None),
                    metadata={"source_type": item.__class__.__name__},
                )
            )


def _records_from_sme_note(note: SMEReviewNote) -> list[TransparencyRecord]:
    records: list[TransparencyRecord] = []
    for warning in note.assumption_warnings or []:
        records.append(
            _record_from_warning(
                warning,
                review_status=note.status,
                release_target=note.release_target,
                metadata={"sme_note_source": note.source, "team_agreement_status": note.team_agreement_status},
            )
        )
    if note.status in REVIEW_QUEUE_STATUSES and note.notes:
        category = _category_from_code_message(note.team_agreement_status, note.notes)
        records.append(
            TransparencyRecord(
                category=category,
                message=note.notes,
                severity=ReviewSeverity.REVIEW_REQUIRED,
                source_phase=_source_phase_from_category(category),
                review_status=note.status,
                release_target=note.release_target,
                code=note.team_agreement_status,
                evidence=None,
                metadata={"sme_note_source": note.source, "tissue_context": note.tissue_context, "endpoint_context": note.endpoint_context},
            )
        )
    return records


def _record_from_warning(
    warning: Any,
    *,
    review_status: SMEReviewStatus | None = None,
    release_target: ReleaseTarget | None = None,
    metadata: JsonDict | None = None,
) -> TransparencyRecord:
    code = str(getattr(warning, "code", None) or getattr(warning, "source_field", None) or "warning")
    message = str(getattr(warning, "message", None) or code)
    category = _category_from_code_message(code, message)
    status = review_status or _coerce_review_status(getattr(warning, "review_status", SMEReviewStatus.UNKNOWN))
    meta: JsonDict = {"field": getattr(warning, "field", None), "source_field": getattr(warning, "source_field", None)}
    meta.update(metadata or {})
    meta.update(_status_metadata_for_category(category, code, message))
    return TransparencyRecord(
        category=category,
        message=message,
        severity=_severity_from_warning(warning, status),
        source_phase=_source_phase_from_category(category),
        review_status=status,
        release_target=release_target or ReleaseTarget.UNKNOWN,
        code=code,
        evidence=getattr(warning, "evidence", None),
        metadata=meta,
    )


def _record_from_warning_mapping(mapping: Mapping[str, Any]) -> TransparencyRecord:
    code = str(mapping.get("code") or mapping.get("source_field") or "warning")
    message = str(mapping.get("message") or code)
    category = _category_from_code_message(code, message)
    status = _coerce_review_status(mapping.get("review_status", SMEReviewStatus.UNKNOWN))
    return TransparencyRecord(
        category=category,
        message=message,
        severity=_severity_from_text(mapping.get("severity"), status),
        source_phase=_source_phase_from_category(category),
        review_status=status,
        release_target=_coerce_release_target(mapping.get("release_target", ReleaseTarget.UNKNOWN)),
        code=code,
        evidence=mapping.get("evidence"),
        metadata=_status_metadata_for_category(category, code, message),
    )


def _metadata_boundary_records(metadata: Mapping[str, Any], records: list[TransparencyRecord]) -> None:
    if metadata.get("public_risk_output") == "not_produced_or_modified" or metadata.get("adjusted_risk_output") is False:
        records.append(
            TransparencyRecord(
                category=AssumptionCategory.PUBLIC_API_BOUNDARY,
                message="Public adjusted-risk output was not produced or modified.",
                severity=ReviewSeverity.INFO,
                source_phase="phase_10",
                review_status=SMEReviewStatus.CURATED,
                release_target=ReleaseTarget.V2_0,
                code="public_adjusted_risk_output_not_modified",
                metadata={"adjusted_risk_output": False},
            )
        )
    if metadata.get("validated_pbpk_ode_gsh_gssg_nrf2") is False or metadata.get("validated_pbpk_ode_model") is False:
        records.append(
            TransparencyRecord(
                category=AssumptionCategory.GSH_MODEL_BOUNDARY,
                message="GSH or burden model metadata states no validated PBPK/ODE claim.",
                severity=ReviewSeverity.INFO,
                source_phase="phase_7_or_8_model_boundary",
                review_status=SMEReviewStatus.UNKNOWN,
                release_target=ReleaseTarget.V2_0,
                code="no_validated_pbpk_ode_claim",
                metadata={"validated_pbpk_ode_claim": False},
            )
        )
    if metadata.get("not_biological_causality") is True:
        records.append(
            TransparencyRecord(
                category=AssumptionCategory.SHAPLEY_ATTRIBUTION_BOUNDARY,
                message="Shapley attribution metadata states model-output attribution, not biological causality.",
                severity=ReviewSeverity.INFO,
                source_phase="phase_9_shapley_interactions",
                review_status=SMEReviewStatus.UNKNOWN,
                release_target=ReleaseTarget.V2_0,
                code="shapley_not_biological_causality",
                metadata={"biological_causality_claim": False},
            )
        )


def _phase4_registry_queue_items() -> list[SMEReviewQueueItem]:
    items: list[SMEReviewQueueItem] = []
    for record in get_reaction_role_sme_records():
        if record.review_status not in REVIEW_QUEUE_STATUSES:
            continue
        message = _first_sme_note(record) or f"{record.enzyme} {record.substrate} requires review."
        category = AssumptionCategory.RELEASE_DEFERRAL if record.review_status is SMEReviewStatus.DEFERRED_3_0 else AssumptionCategory.REACTION_ROLE
        if record.record_id and "tce" in record.record_id:
            message = "TCE x CYP2E1 candidate/pending status is visible and inactive pending team agreement."
        items.append(
            SMEReviewQueueItem(
                category=category,
                message=message,
                severity=ReviewSeverity.REVIEW_REQUIRED,
                source_phase="phase_4_reaction_role_semantics",
                review_status=record.review_status,
                release_target=record.release_target,
                reason="phase4_registry_review_status",
                code=record.record_id,
                evidence=record.evidence,
                metadata={
                    "enzyme": record.enzyme,
                    "substrate": record.substrate,
                    "active": record.active,
                    "record_id": record.record_id,
                    "candidate_pending": record.review_status in {SMEReviewStatus.CANDIDATE, SMEReviewStatus.PENDING_TEAM_AGREEMENT},
                },
            )
        )
    return items


def _default_boundary_queue_items() -> list[SMEReviewQueueItem]:
    return [
        SMEReviewQueueItem(
            category=AssumptionCategory.GSH_TISSUE_PRESET,
            message="GSH tissue presets are local 2.0 placeholders unless curated later.",
            severity=ReviewSeverity.REVIEW_REQUIRED,
            source_phase="phase_7_gsh_redox_capacity",
            review_status=SMEReviewStatus.UNKNOWN,
            release_target=ReleaseTarget.V2_0,
            reason="local_2_0_preset_caveat",
            code="local_2_0_gsh_preset",
            metadata={"local_2_0_preset": True},
        ),
        SMEReviewQueueItem(
            category=AssumptionCategory.GSH_MODEL_BOUNDARY,
            message="GSH redox-capacity model is not PBPK/ODE/GSH-GSSG/Nrf2 validated.",
            severity=ReviewSeverity.REVIEW_REQUIRED,
            source_phase="phase_7_gsh_redox_capacity",
            review_status=SMEReviewStatus.UNKNOWN,
            release_target=ReleaseTarget.V2_0,
            reason="model_boundary_caveat",
            code="gsh_not_pbpk_ode_validated",
            metadata={"validated_pbpk_ode_gsh_gssg_nrf2": False},
        ),
        SMEReviewQueueItem(
            category=AssumptionCategory.SHAPLEY_ATTRIBUTION_BOUNDARY,
            message="Shapley output is model-output attribution only and not biological causality proof.",
            severity=ReviewSeverity.REVIEW_REQUIRED,
            source_phase="phase_9_shapley_interactions",
            review_status=SMEReviewStatus.UNKNOWN,
            release_target=ReleaseTarget.V2_0,
            reason="model_boundary_caveat",
            code="shapley_not_causality",
            metadata={"biological_causality_claim": False},
        ),
    ]


def _queue_item_from_record(record: TransparencyRecord) -> SMEReviewQueueItem:
    return SMEReviewQueueItem(
        category=record.category,
        message=record.message,
        severity=record.severity if record.severity is not ReviewSeverity.INFO else ReviewSeverity.REVIEW_REQUIRED,
        source_phase=record.source_phase,
        review_status=record.review_status,
        release_target=record.release_target,
        reason="warning_or_unknown_status",
        code=record.code,
        evidence=record.evidence,
        metadata=record.metadata,
    )


def _needs_review(record: TransparencyRecord) -> bool:
    return record.review_status in REVIEW_QUEUE_STATUSES or record.category in {
        AssumptionCategory.RELEASE_DEFERRAL,
        AssumptionCategory.UNKNOWN_BIOLOGY,
        AssumptionCategory.GSH_MODEL_BOUNDARY,
        AssumptionCategory.SHAPLEY_ATTRIBUTION_BOUNDARY,
    }


def _looks_like_warning(item: Any) -> bool:
    return hasattr(item, "code") and hasattr(item, "message")


def _mapping_looks_like_warning(item: Mapping[str, Any]) -> bool:
    return "message" in item and ("code" in item or "source_field" in item)


def _category_from_code_message(code: object, message: object = "") -> AssumptionCategory:
    text = f"{code or ''} {message or ''}".lower()
    if "product_carcinogenic" in text:
        return AssumptionCategory.PRODUCT_CARCINOGENIC_GUARDRAIL
    if "ic50" in text:
        return AssumptionCategory.IC50_CONVERSION
    if any(token in text for token in ("affinity", "rdkit", "ecfp4", "tanimoto", "docking")):
        return AssumptionCategory.AFFINITY_FALLBACK
    if "km" in text or "ki" in text:
        return AssumptionCategory.KINETIC_PARAMETER
    if "deferred" in text or "3.0" in text or "v3_0" in text:
        return AssumptionCategory.RELEASE_DEFERRAL
    if "reaction_role" in text or "endpoint_role" in text or "risk_direction" in text or "tce" in text or "cyp2e1" in text:
        return AssumptionCategory.REACTION_ROLE
    if "gsh" in text and any(token in text for token in ("pbpk", "ode", "gssg", "nrf2", "boundary", "validated")):
        return AssumptionCategory.GSH_MODEL_BOUNDARY
    if "gsh" in text or "tissue" in text or "preset" in text:
        return AssumptionCategory.GSH_TISSUE_PRESET
    if "effective" in text or "burden" in text:
        return AssumptionCategory.EFFECTIVE_BURDEN_BOUNDARY
    if "shapley" in text or "attribution" in text or "causality" in text:
        return AssumptionCategory.SHAPLEY_ATTRIBUTION_BOUNDARY
    if "adjusted_risk" in text or "public" in text or "api" in text:
        return AssumptionCategory.PUBLIC_API_BOUNDARY
    if "unknown" in text or "unresolved" in text:
        return AssumptionCategory.UNKNOWN_BIOLOGY
    return AssumptionCategory.GENERAL


def _source_phase_from_category(category: AssumptionCategory) -> str:
    mapping = {
        AssumptionCategory.REACTION_ROLE: "phase_4_reaction_role_semantics",
        AssumptionCategory.KINETIC_PARAMETER: "phase_5_kinetic_resolution",
        AssumptionCategory.IC50_CONVERSION: "phase_5_kinetic_resolution",
        AssumptionCategory.AFFINITY_FALLBACK: "phase_5_kinetic_resolution",
        AssumptionCategory.GSH_TISSUE_PRESET: "phase_7_gsh_redox_capacity",
        AssumptionCategory.GSH_MODEL_BOUNDARY: "phase_7_gsh_redox_capacity",
        AssumptionCategory.EFFECTIVE_BURDEN_BOUNDARY: "phase_8_effective_burden",
        AssumptionCategory.SHAPLEY_ATTRIBUTION_BOUNDARY: "phase_9_shapley_interactions",
        AssumptionCategory.RELEASE_DEFERRAL: "phase_4_reaction_role_semantics",
        AssumptionCategory.PRODUCT_CARCINOGENIC_GUARDRAIL: "phase_4_reaction_role_semantics",
        AssumptionCategory.PUBLIC_API_BOUNDARY: "phase_10_public_api_boundary",
    }
    return mapping.get(category, "phase_10_transparency")


def _severity_from_warning(warning: Any, review_status: SMEReviewStatus) -> ReviewSeverity:
    return _severity_from_text(getattr(warning, "severity", None), review_status)


def _severity_from_text(value: object, review_status: SMEReviewStatus) -> ReviewSeverity:
    text = str(value or "warning").lower()
    if "block" in text or "error" in text:
        return ReviewSeverity.BLOCKER
    if review_status in REVIEW_QUEUE_STATUSES:
        return ReviewSeverity.REVIEW_REQUIRED
    if "info" in text:
        return ReviewSeverity.INFO
    return ReviewSeverity.WARNING


def _coerce_review_status(value: object) -> SMEReviewStatus:
    if isinstance(value, SMEReviewStatus):
        return value
    try:
        return SMEReviewStatus(str(value))
    except ValueError:
        return SMEReviewStatus.UNKNOWN


def _coerce_release_target(value: object) -> ReleaseTarget:
    if isinstance(value, ReleaseTarget):
        return value
    try:
        return ReleaseTarget(str(value))
    except ValueError:
        return ReleaseTarget.UNKNOWN


def _status_metadata_for_category(category: AssumptionCategory, code: str, message: str) -> JsonDict:
    text = f"{code} {message}".lower()
    metadata: JsonDict = {}
    if category is AssumptionCategory.IC50_CONVERSION:
        metadata["conversion_claimed"] = False
    if category is AssumptionCategory.AFFINITY_FALLBACK:
        metadata["conversion_claimed"] = False
        metadata["rdkit_ecfp4_tanimoto_active"] = False
    if category is AssumptionCategory.GSH_TISSUE_PRESET:
        metadata["local_2_0_preset"] = True
    if category is AssumptionCategory.PUBLIC_API_BOUNDARY:
        metadata["adjusted_risk_output"] = False
    if "km" in text and "proxy" in text:
        metadata["low_confidence_proxy"] = True
    return metadata


def _sort_records(records: list[TransparencyRecord]) -> list[TransparencyRecord]:
    return sorted(records, key=lambda item: (item.category.value, item.severity.value, item.source_phase, item.code or "", item.message))


def _dedupe_records(records: list[TransparencyRecord]) -> list[TransparencyRecord]:
    seen: set[tuple[str, str, str, str | None]] = set()
    unique: list[TransparencyRecord] = []
    for record in records:
        key = (record.category.value, record.source_phase, record.message, record.code)
        if key not in seen:
            seen.add(key)
            unique.append(record)
    return unique


def _sort_queue_items(items: list[SMEReviewQueueItem]) -> list[SMEReviewQueueItem]:
    return sorted(items, key=lambda item: (item.category.value, item.release_target.value, item.source_phase, item.code or "", item.message))


def _dedupe_queue_items(items: list[SMEReviewQueueItem]) -> list[SMEReviewQueueItem]:
    seen: set[tuple[str, str, str | None]] = set()
    unique: list[SMEReviewQueueItem] = []
    for item in items:
        key = (item.category.value, item.message, item.code)
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique


def _sort_deferrals(items: list[DeferralRecord]) -> list[DeferralRecord]:
    return sorted(items, key=lambda item: (item.release_target.value, item.subject, item.message))


def _dedupe_deferrals(items: list[DeferralRecord]) -> list[DeferralRecord]:
    seen: set[tuple[str, str]] = set()
    unique: list[DeferralRecord] = []
    for item in items:
        key = (item.subject, item.release_target.value)
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique


def _count_by_category(records: list[TransparencyRecord]) -> JsonDict:
    counts = {category.value: 0 for category in AssumptionCategory}
    for record in records:
        counts[record.category.value] += 1
    return {key: counts[key] for key in sorted(counts)}


def _count_by_severity(records: list[TransparencyRecord]) -> JsonDict:
    counts = {severity.value: 0 for severity in ReviewSeverity}
    for record in records:
        counts[record.severity.value] += 1
    return {key: counts[key] for key in sorted(counts)}


def _release_target_summary(queue: list[SMEReviewQueueItem], deferrals: list[DeferralRecord]) -> JsonDict:
    counts: dict[str, int] = {}
    for item in list(queue) + list(deferrals):
        key = item.release_target.value
        counts[key] = counts.get(key, 0) + 1
    return {key: counts[key] for key in sorted(counts)}


def _summarize_evidence(records: list[TransparencyRecord], boundaries: list[ModelBoundaryStatement]) -> JsonDict:
    grades: dict[str, int] = {}
    sources: set[str] = set()
    for evidence in _iter_evidence(records, boundaries):
        grades[evidence.grade.value] = grades.get(evidence.grade.value, 0) + 1
        if evidence.source:
            sources.add(evidence.source)
    return {
        "evidence_grade_counts": {key: grades[key] for key in sorted(grades)},
        "sources": sorted(sources),
        "provenance_is_visible": True,
    }


def _iter_evidence(records: list[TransparencyRecord], boundaries: list[ModelBoundaryStatement]) -> Iterable[EvidenceRecord]:
    for record in records:
        evidence = record.evidence
        if isinstance(evidence, EvidenceRecord):
            yield evidence
        elif isinstance(evidence, list):
            for item in evidence:
                if isinstance(item, EvidenceRecord):
                    yield item
    for boundary in boundaries:
        if boundary.evidence is not None:
            yield boundary.evidence


def _deferred_subject(enzyme: str) -> str:
    return f"{enzyme} broader detox-arm semantics"


def _first_sme_note(record: Any) -> str | None:
    notes = getattr(record, "sme_notes", None) or []
    if not notes:
        return None
    return notes[0].notes
