"""Pydantic data models for the Carcinogen-Gene Knowledge Graph."""

from __future__ import annotations

import hashlib
import re
from enum import Enum
from typing import Any, Optional, Dict

from pydantic import (
    AliasChoices,
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)

# ── Enums ────────────────────────────────────────────────────────────────


class NodeType(str, Enum):
    CARCINOGEN = "Carcinogen"
    ENZYME = "Enzyme"
    GENE = "Gene"
    METABOLITE = "Metabolite"
    DNA_ADDUCT = "DNA_Adduct"
    PATHWAY = "Pathway"
    TISSUE = "Tissue"


class EdgeType(str, Enum):
    ACTIVATES = "ACTIVATES"
    DETOXIFIES = "DETOXIFIES"
    TRANSPORTS = "TRANSPORTS"
    FORMS_ADDUCT = "FORMS_ADDUCT"
    REPAIRS = "REPAIRS"
    PATHWAY = "PATHWAY"
    EXPRESSED_IN = "EXPRESSED_IN"
    INDUCES = "INDUCES"
    INHIBITS = "INHIBITS"
    ENCODES = "ENCODES"
    CUSTOM = "CUSTOM"
    # new edge types
    ACCUMULATES = "ACCUMULATES"
    GENERATES = "GENERATES"
    DEPLETES = "DEPLETES"
    ANTAGONIZES = "ANTAGONIZES"
    REGULATES = "REGULATES"
    STABILIZES = "STABILIZES"
    METABOLIZES = "METABOLIZES"
    CO_EXPOSED = "CO_EXPOSED"


class CurationStatus(str, Enum):
    DRAFT = "Draft"
    IN_REVIEW = "In Review"
    REVIEWED = "Reviewed"
    APPROVED = "Approved"
    REJECTED = "Rejected"


class CurationConfidence(str, Enum):
    LOW = "Low"
    MEDIUM = "Medium"
    HIGH = "High"


class RecordOrigin(str, Enum):
    IMPORTED = "imported"
    SEEDED = "seeded"
    USER = "user"
    LLM = "llm"


class MatchStatus(str, Enum):
    UNKNOWN = "unknown"
    CANONICAL = "canonical"
    ALIAS = "alias"
    UNMATCHED = "unmatched"
    CUSTOM = "custom"


class ProvenanceRecord(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    source_db: Optional[str] = None
    record_id: Optional[str] = Field(
        default=None,
        validation_alias=AliasChoices("record_id", "accession"),
    )
    evidence: Optional[str] = None
    pmid: Optional[str] = None
    tissue: Optional[str] = None
    citation: Optional[str] = None
    url: Optional[str] = None


class CurationRecord(BaseModel):
    status: CurationStatus = CurationStatus.DRAFT
    confidence: Optional[CurationConfidence] = None
    curator: Optional[str] = None
    reviewed_by: Optional[str] = None
    reviewed_at: Optional[str] = None
    notes: Optional[str] = None

    @field_validator("reviewed_at")
    @classmethod
    def _validate_reviewed_at(cls, v: str | None) -> str | None:
        if v is None:
            return v
        import datetime

        for fmt in ("%Y-%m-%d", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%dT%H:%M:%S%z"):
            try:
                datetime.datetime.strptime(v, fmt)
                return v
            except ValueError:
                continue
        raise ValueError(
            f"reviewed_at must be a date string (YYYY-MM-DD or ISO 8601), got {v!r}"
        )


def _join_unique(values: list[str]) -> str | None:
    unique = dict.fromkeys(v for value in values if (v := value.strip()))
    return "; ".join(unique) if unique else None


def _first_nonempty(values: list[str]) -> str | None:
    for value in values:
        cleaned = value.strip()
        if cleaned:
            return cleaned
    return None


def _normalize_provenance_fields(
    owner: BaseModel, *, summary_only_fields: tuple[str, ...]
) -> None:
    provenance = list(getattr(owner, "provenance", []))
    if not provenance:
        legacy = ProvenanceRecord(
            source_db=getattr(owner, "source_db", None),
            evidence=getattr(owner, "evidence", None),
            pmid=getattr(owner, "pmid", None),
            tissue=getattr(owner, "tissue", None),
        )
        if any(legacy.model_dump(exclude_none=True).values()):
            provenance = [legacy]
            setattr(owner, "provenance", provenance)

    if not provenance:
        return

    for field in ("source_db", "pmid", "tissue"):
        if getattr(owner, field, None) is None:
            joined = _join_unique(
                [getattr(item, field) for item in provenance if getattr(item, field)]
            )
            setattr(owner, field, joined)

    for field in summary_only_fields:
        if getattr(owner, field, None) is None:
            first_value = _first_nonempty(
                [getattr(item, field) for item in provenance if getattr(item, field)]
            )
            setattr(owner, field, first_value)


# ── Node ─────────────────────────────────────────────────────────────────


class Node(BaseModel):
    id: str
    label: str
    type: NodeType
    detail: str = ""
    group: Optional[str] = None
    iarc: Optional[str] = None
    phase: Optional[str] = None
    role: Optional[str] = None
    reactivity: Optional[str] = None
    source_db: Optional[str] = None
    evidence: Optional[str] = None
    pmid: Optional[str] = None
    tissue: Optional[str] = None
    exposure: Optional[str] = None
    variant: Optional[str] = None
    phenotype: Optional[str] = None
    activity_score: Optional[float] = None
    tissue_weights: Optional[dict[str, float]] = None
    exposure_scenarios: Optional[dict[str, Any]] = None
    tier: Optional[int] = None
    origin: RecordOrigin = RecordOrigin.IMPORTED
    match_status: MatchStatus = MatchStatus.UNKNOWN
    canonical_id: Optional[str] = None
    canonical_label: Optional[str] = None
    canonical_namespace: Optional[str] = None
    custom_type: Optional[str] = None
    provenance: list[ProvenanceRecord] = Field(default_factory=list)
    curation: Optional[CurationRecord] = None

    @classmethod
    def generate_id(cls, label: str) -> str:
        """Derive a safe node ID from a label.

        Simple labels (e.g. ``"CYP1A1"``) pass through unchanged.  Labels
        containing special characters that are stripped during sanitisation
        get a short hash suffix to avoid collisions (e.g.
        ``"Benzo[a]pyrene"`` → ``"Benzo_a_pyrene_a4f2c1"``).
        """
        sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_.")
        if not sanitized:
            raise ValueError("Cannot generate an ID from an empty label")
        # Append hash only when characters were actually stripped
        if sanitized != label:
            short_hash = hashlib.sha256(label.encode()).hexdigest()[:6]
            sanitized = f"{sanitized}_{short_hash}"
        return sanitized

    @model_validator(mode="after")
    def _normalize(self) -> "Node":
        if not self.id:
            self.id = self.generate_id(self.label)
        _normalize_provenance_fields(self, summary_only_fields=("evidence",))
        if self.match_status == MatchStatus.CANONICAL:
            self.canonical_id = self.canonical_id or self.id
            self.canonical_label = self.canonical_label or self.label
        elif self.match_status == MatchStatus.ALIAS:
            if not self.canonical_id:
                raise ValueError("Alias-matched nodes must define canonical_id")
            if not self.canonical_label:
                raise ValueError("Alias-matched nodes must define canonical_label")
        elif self.match_status == MatchStatus.CUSTOM and not self.custom_type:
            raise ValueError("Custom nodes must define custom_type")
        return self


# ── Edge ─────────────────────────────────────────────────────────────────


class Edge(BaseModel):
    source: str
    target: str
    type: EdgeType
    label: Optional[str] = None
    carcinogen: Optional[str] = None
    source_db: Optional[str] = None
    evidence: Optional[str] = None
    pmid: Optional[str] = None
    tissue: Optional[str] = None
    kinetics: Optional[dict[str, Any]] = None
    origin: RecordOrigin = RecordOrigin.IMPORTED
    match_status: MatchStatus = MatchStatus.UNKNOWN
    canonical_predicate: Optional[str] = None
    canonical_namespace: Optional[str] = None
    custom_predicate: Optional[str] = None
    provenance: list[ProvenanceRecord] = Field(default_factory=list)
    curation: Optional[CurationRecord] = None

    @model_validator(mode="after")
    def _normalize(self) -> "Edge":
        _normalize_provenance_fields(self, summary_only_fields=("evidence",))
        if self.type == EdgeType.CUSTOM and not self.custom_predicate:
            raise ValueError("Edges with type CUSTOM must define custom_predicate")
        if self.type == EdgeType.CUSTOM and self.match_status in (
            MatchStatus.CANONICAL,
            MatchStatus.ALIAS,
        ):
            raise ValueError(
                "Edges with type CUSTOM cannot be canonical or alias-matched"
            )
        if self.match_status in (MatchStatus.CANONICAL, MatchStatus.ALIAS):
            self.canonical_predicate = self.canonical_predicate or self.type.value
        elif self.match_status == MatchStatus.CUSTOM and not self.custom_predicate:
            raise ValueError("Custom edges must define custom_predicate")
        return self


# ── Top-level graph container ────────────────────────────────────────────


class KnowledgeGraph(BaseModel):
    nodes: list[Node] = Field(default_factory=list)
    edges: list[Edge] = Field(default_factory=list)

    @model_validator(mode="after")
    def _check_edge_references(self) -> "KnowledgeGraph":
        node_ids = {n.id for n in self.nodes}
        bad: list[str] = []
        for edge in self.edges:
            if edge.source not in node_ids:
                bad.append(f"Edge references missing source node: {edge.source!r}")
            if edge.target not in node_ids:
                bad.append(f"Edge references missing target node: {edge.target!r}")
            if edge.carcinogen and edge.carcinogen not in node_ids:
                bad.append(
                    f"Edge references missing carcinogen node: {edge.carcinogen!r}"
                )
        if bad:
            raise ValueError(
                f"Referential integrity errors ({len(bad)}):\n  " + "\n  ".join(bad)
            )
        return self


# --- Data classes for interactions, modifiers, and parameters added after base load
class GenotypeModifiers(BaseModel):
    activity_multiplier: float
    frequency: float
    alleles: Optional[str] = None
    note: Optional[str] = None
    model_config = {"extra": "ignore"}


class ModifierMetrics(BaseModel):
    activity_multiplier: float
    frequency: float
    alleles: Optional[str] = None
    note: Optional[str] = None
    model_config = {"extra": "ignore"}


class GenotypeModifiersContainer(BaseModel):
    """Maps genotype keys (e.g., 'UM', 'active') to their metric data."""

    model_config = {"extra": "allow"}

    @property
    def modifiers(self) -> Dict[str, ModifierMetrics]:
        """Helper to extract just the validated data dictionaries, ignoring description."""
        print(self.model_dump(by_alias=True))
        return {
            k: ModifierMetrics(**v)
            for k, v in self.model_dump(by_alias=True).items()
            if k != "_description"
        }


class InteractionEdgeDataBase(BaseModel):
    # key of outer dictionary is enzyme, key of inner dictionary is substrate
    Km_uM: float
    Vmax_relative: float
    relative_priority: int
    product: str
    product_carcinogenic: bool
    assumed_ki: bool = False
    parameter_provenance_ref: Optional[str] = None
    notes: Optional[str] = None
    Ki_uM: Optional[float] = None


# edge data between enzyme and substrate. Eventually connect back to carcinogen
class CompetitiveInhibitionSubstrate(InteractionEdgeDataBase):
    pass


class Phase2ConjugationSubstrate(InteractionEdgeDataBase):
    pass


class EnzymeInduction(BaseModel):
    fold_induction: float = (Field(ge=0),)
    range_min: float = Field(ge=0)
    range_max: float = Field(ge=0)
    mechanism: str
    tissue_specificity: str
    notes: str
