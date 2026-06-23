from typing import List, Optional, Literal, Union, Annotated, Dict, Any
from pydantic import (
    BaseModel,
    Field,
    field_validator,
    TypeAdapter,
    Discriminator,
    model_validator,
)

# --- SHARED TYPES ---


"""
need to indicate which data classes are affiliated with which 
node or edge
"""

# --- Non-visible Nodes ---


class Lifesyle(BaseModel):
    type: Literal["Lifestyle"]
    visible: bool = False
    pass


class Patient(BaseModel):
    type: Literal["Patient"]
    visible: bool = False


class Tissue(BaseModel):
    type: Literal["Tissue"]
    visible: bool = False


class Substrate(BaseModel):
    type: Literal["Substrate"]
    visible: bool = False


class CarcinogenClass(BaseModel):
    type: str = None
    class_id: int
    class_label: str
    index_carcinogen: str
    iarc_classification: str
    primary_targets: List[str]
    key_enzymes_activation: List[str]
    key_enzymes_detox: List[str]
    visible: bool = False
    epa_cancer_slope_factor: Optional[CancerSlopeFactor] = None
    epa_inhalation_unit_risk: Optional[InhalationUnitRisk] = None
    regulatory_limits: Dict[str, str] = Field(default_factory=dict)
    population_prevalence_notes: Optional[str] = None

    @model_validator(mode="before")
    @classmethod
    def populate_type_field(cls, data: Any) -> Any:
        if isinstance(data, dict):
            # If 'type' is not explicitly provided, default it to 'class_label'
            if "type" not in data or data["type"] is None:
                data["type"] = data.get("class_label")
        return data

    model_config = {"validate_assignment": True}


class ExposureScenario(BaseModel):
    label: str
    estimated_tissue_conc_uM: float = Field(ge=0)
    daily_intake_ug_kg: float = Field(ge=0)
    multiplier_vs_baseline: float = Field(ge=0)
    source: str
    type: str = None
    biomarker: Optional[str] = None
    air_concentration_ng_m3: Optional[float] = Field(default=None, ge=0)
    dietary_bap_ug_day: Optional[float] = Field(default=None, ge=0)
    note: Optional[str] = None

    @model_validator(mode="before")
    @classmethod
    def populate_type_field(cls, data: Any) -> Any:
        if isinstance(data, dict):
            # If 'type' is not explicitly provided, default it to 'label'
            if "type" not in data or data["type"] is None:
                data["type"] = data.get("label")
        return data


class CancerSlopeFactor(BaseModel):
    value: float
    unit: str
    adaf_adjusted_value: Optional[float] = None
    source: str


class InhalationUnitRisk(BaseModel):
    value: float
    unit: str
    source: str


class CarcinogenClass(BaseModel):
    class_id: int
    class_label: str
    index_carcinogen: str
    iarc_classification: str = Field(pattern=r"^Group (1|2A|2B|3|4)$")
    primary_targets: List[str]
    key_enzymes_activation: List[str]
    key_enzymes_detox: List[str]
    exposure_scenarios: Dict[str, ExposureScenario]
    epa_cancer_slope_factor: Optional[CancerSlopeFactor] = None
    epa_inhalation_unit_risk: Optional[InhalationUnitRisk] = None
    regulatory_limits: Dict[str, str] = Field(default_factory=dict)
    population_prevalence_notes: Optional[str] = None


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


class Phase2Conjugation(InteractionEdgeDataBase):
    pass


# per enzyme within competitive_inhitibition dict in interaction_parameters.json
class CompetitiveInduction(BaseModel):
    fold: float
    range_min: float
    range_max: float
    mechanism: str


class GshDepletionBiologyModelData(BaseModel):
    Ki_feedback_mM: float
    n_feedback: float
    Km_GST_mM: float
    _param_notes: str


# class GshDepletionBaseData(BaseModel):
#     baseline_gsh_mM: float
#     baseline_gsh_range_mM: List[float]
#     critical_threshold_fraction: float
#     critical_threshold_mM: float
#     synthesis_rate_umol_h_g: float
#     half_life_h: float
#     liver_volume_g: float
#     notes: str
#     biology_model: GshDepletionBiologyModelData


# class GshDepletionModel(BaseModel):
#     gsh_base_data = GshDepletionBaseData
#     pass  # come back to, interaction parameters are not clear here


class EnzymeInductionFromLifestyle(BaseModel):
    fold_induction: float
    range_min: float
    range_max: float
    mechanism: str
    tissue_specificity: Optional[str]
    notes: Optional[str]


class TissueWeights(BaseModel):
    """Enforces strict, explicit tissue key verification."""

    Liver: float = 0.0
    Lung: float = 0.0
    Prostate: float = 0.0
    Bladder: float = 0.0
    Colon: float = 0.0
    Breast: float = 0.0
    Kidney: float = 0.0
    Esophagus: float = 0.0
    visible: bool = False
    model_config = {"extra": "allow"}


# --- Base Node Data ---


class NodeProvenance(BaseModel):
    source_db: Optional[str] = None
    record_id: Optional[str] = None
    evidence: Optional[str] = None
    citation: Optional[str] = None


class BaseNode(BaseModel):
    """Fields that every single node in your graph must have."""

    id: str
    label: str
    detail: str
    origin: Optional[str] = None
    match_status: Optional[str] = None
    provenance: List[NodeProvenance] = Field(default_factory=list)
    # begin optional fields
    source_db: Optional[str] = None
    canonical_id: Optional[str] = None
    canonical_label: Optional[str] = None
    canonical_namespace: Optional[str] = None
    evidence: Optional[str] = None
    visible: Optional[bool] = True


# --- TYPE-SPECIFIC NODE MODELS ---


class CarcinogenNode(BaseNode):
    type: Literal["Carcinogen"]
    group: str
    exposure: Optional[str] = None
    iarc: Optional[str] = None
    model_config = {"validate_assignment": True}


class EnzymeNode(BaseNode):
    type: Literal["Enzyme"]
    role: str
    phase: Optional[str] = None
    activity_score: Optional[float] = None
    tier: Optional[int] = None
    variant: Optional[str] = None
    group: Optional[str] = None
    tissue: Optional[str] = None
    tissue_weights: TissueWeights = Field(default_factory=TissueWeights)
    genotype_modifiers: Dict[str, ModifierMetrics] = Field(default_factory=dict)
    model_config = {"validate_assignment": True}


class MetaboliteNode(BaseNode):
    type: Literal["Metabolite"]
    reactivity: str
    model_config = {"validate_assignment": True}


class DnaAdductNode(BaseNode):
    type: Literal["DNA_Adduct"]
    reactivity: Optional[str] = None
    model_config = {"validate_assignment": True}


class PathwayNode(BaseNode):
    type: Literal["Pathway"]
    model_config = {"validate_assignment": True}


class GeneNode(BaseNode):
    type: Literal["Gene"]
    model_config = {"validate_assignment": True}


# --- Edge Data --- #


class EdgeProvenance(BaseModel):
    source_db: Optional[str] = None
    record_id: Optional[str] = None
    evidence: Optional[str] = None
    citation: Optional[str] = None
    url: Optional[str] = None  # Optional since some entries have URLs


# --- BASE EDGE MODEL ---
class BaseEdge(BaseModel):
    """Core parameters required on every single relationship link."""

    visible: Optional[bool] = True
    source: str = Field(..., description="The ID of the source node")
    target: str = Field(..., description="The ID of the target node")
    match_status: Optional[str] = None
    origin: Optional[str] = None
    provenance: List[EdgeProvenance] = Field(default_factory=list)
    label: Optional[str] = None
    source_db: Optional[str] = None
    evidence: Optional[str] = None
    canonical_predicate: Optional[str] = None
    canonical_namespace: Optional[str] = None
    model_config = {"validate_assignment": True}


# -- Specific Edge Models ---


class VisibleEdge(BaseEdge):
    type: Literal[
        "PATHWAY",
        "METABOLIZES",
        "ACTIVATES",
        "FORMS_ADDUCT",
        "DETOXIFIES",
        "TRANSPORTS",
        "REPAIRS",
        "INDUCES",
        "INHIBITS",
        "DEPLETES",  # NEW
        "GENERATES",
        "ANTAGONIZES",
        "REGULATES",
        "STABILIZES",
        "CO_EXPOSED",
        "ACCUMULATES",
    ]
    model_config = {"validate_assignment": True}


AnyNode = Union[
    CarcinogenNode, EnzymeNode, MetaboliteNode, DnaAdductNode, PathwayNode, GeneNode
]
AnnotatedNode = Annotated[AnyNode, Discriminator("type")]
AnyEdge = Union[VisibleEdge]
AnnotatedEdge = Annotated[AnyEdge, Discriminator("type")]
node_adapter = TypeAdapter(List[AnnotatedNode])
edge_adapter = TypeAdapter(List[AnnotatedEdge])
