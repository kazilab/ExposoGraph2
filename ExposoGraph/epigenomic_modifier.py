"""Epigenomic (promoter-methylation) modifier layer.

Implements the epigenomic modifier module anticipated in ExposoGraph 2.0
(Manuscript v6, Discussion / Future development):

    "An epigenomic modifier layer (CYP1A1/GSTP1/MGMT promoter methylation)
    will be developed as a parallel somatic modifier analogous to the CHIP
    module."

The layer models aberrant promoter methylation of carcinogen-relevant genes as
a post-hoc modifier applied on top of germline flux and interaction results,
analogous to :mod:`ExposoGraph.chip_modifier`.

Inputs may be either beta values (typical EPIC/450K bisulfite output, 0-1) or
binary hyper/hypo-methylated flags. Hyper- vs hypo-methylation has opposite
functional consequences depending on whether the gene is a bioactivator
(CYP1A1/CYP1B1) or a detoxifier / repair enzyme (GSTP1, MGMT, OGG1).

Mechanistic grounding (Manuscript v6, CHIP / Discussion sections):
* GSTP1 promoter hypermethylation is a well-documented event in prostate,
  breast, and hepatic tumours; it reduces Phase II conjugation capacity.
* MGMT promoter hypermethylation silences O6-methylguanine repair, amplifying
  the mutagenic burden of nitrosamines and alkylating agents.
* CYP1A1 hypermethylation typically reduces AhR-induced bioactivation and is
  therefore protective for PAH and estrogen activation.
* AHRR hypomethylation is the canonical smoking biomarker (illumina cg05575921)
  and is treated here as an exposure-intensity correlate rather than a risk
  multiplier.

All multipliers are conservative model estimates consistent with published
epidemiology and are expected to be refined as methylation-stratified enzyme
activity data become available.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable

# ── Constants ──────────────────────────────────────────────────────────────

# Default beta-value thresholds. Beta > HYPER_BETA_THRESHOLD is treated as
# hypermethylated, beta < HYPO_BETA_THRESHOLD as hypomethylated. Values between
# are "normal" and do not trigger a modifier.
HYPER_BETA_THRESHOLD: float = 0.50
HYPO_BETA_THRESHOLD: float = 0.20

# Default aggregate modifier when *any* detrimental methylation event is
# present without an applicable gene-specific rule (fallback, conservative).
AGGREGATE_METHYLATION_MODIFIER: float = 1.2


# Per-gene methylation effect tables. Each entry specifies:
#   * direction: "hyper" or "hypo" methylation that triggers the effect
#   * role: "bioactivator" (silencing protective) or "detox_repair"
#     (silencing harmful)
#   * carcinogen_class_multipliers: multiplicative risk modifiers applied to
#     the matching carcinogen classes when the trigger direction is observed
#
# A gene with direction="hyper" role="detox_repair" means: when the promoter
# is hypermethylated, detox/repair capacity drops, so risk *increases* (values
# > 1). For bioactivators, hypermethylation *decreases* risk (values < 1).
EPIGENETIC_GENE_EFFECTS: dict[str, dict[str, Any]] = {
    "GSTP1": {
        "direction": "hyper",
        "role": "detox_repair",
        "description": (
            "GSTP1 promoter hypermethylation silences pi-class glutathione "
            "S-transferase, a canonical Phase II conjugation enzyme for PAH "
            "diol-epoxides and oxidative-stress electrophiles."
        ),
        "carcinogen_class_multipliers": {
            "PAH": 1.5,
            "Estrogen": 1.3,
            "Aflatoxin": 1.3,
            "HeavyMetal": 1.2,
        },
        "mechanism": "phase2_conjugation_loss",
    },
    "MGMT": {
        "direction": "hyper",
        "role": "detox_repair",
        "description": (
            "MGMT hypermethylation silences O6-methylguanine-DNA "
            "methyltransferase, the dedicated repair enzyme for alkylated "
            "guanine adducts produced by nitrosamines and alkylating agents."
        ),
        "carcinogen_class_multipliers": {
            "Nitrosamine": 1.7,
            "NDMA": 1.7,
            "AlkylatingAgent": 2.0,
            "PAH": 1.1,
        },
        "mechanism": "repair_impairment",
    },
    "OGG1": {
        "direction": "hyper",
        "role": "detox_repair",
        "description": (
            "OGG1 hypermethylation reduces base-excision repair of 8-oxoG "
            "lesions produced by oxidative stress from PAH, heavy metals, "
            "and benzene metabolites."
        ),
        "carcinogen_class_multipliers": {
            "PAH": 1.3,
            "Benzene": 1.4,
            "HeavyMetal": 1.3,
            "Aflatoxin": 1.2,
        },
        "mechanism": "oxidative_repair_loss",
    },
    "CYP1A1": {
        "direction": "hyper",
        "role": "bioactivator",
        "description": (
            "CYP1A1 promoter hypermethylation reduces AhR-inducible "
            "bioactivation of PAH, estrogens, and dioxins (protective)."
        ),
        "carcinogen_class_multipliers": {
            "PAH": 0.8,
            "Estrogen": 0.85,
            "Dioxin": 0.8,
        },
        "mechanism": "bioactivation_loss_protective",
    },
    "CYP1B1": {
        "direction": "hyper",
        "role": "bioactivator",
        "description": (
            "CYP1B1 promoter hypermethylation reduces 4-hydroxylation of "
            "17-beta-estradiol, lowering catechol-estrogen bioactivation in "
            "mammary and prostate tissue."
        ),
        "carcinogen_class_multipliers": {
            "Estrogen": 0.8,
            "PAH": 0.9,
        },
        "mechanism": "bioactivation_loss_protective",
    },
    "CDKN2A": {
        "direction": "hyper",
        "role": "tumor_suppressor",
        "description": (
            "CDKN2A (p16INK4a) promoter hypermethylation silences a core "
            "cell-cycle tumor suppressor; amplifies mutagenic consequences "
            "of any DNA-reactive carcinogen."
        ),
        "carcinogen_class_multipliers": {
            "PAH": 1.4,
            "Nitrosamine": 1.4,
            "Aflatoxin": 1.4,
            "AlkylatingAgent": 1.4,
            "HCA": 1.3,
            "HeavyMetal": 1.3,
            "Benzene": 1.3,
        },
        "mechanism": "tumor_suppressor_silencing",
    },
    "AHRR": {
        "direction": "hypo",
        "role": "exposure_biomarker",
        "description": (
            "AHRR hypomethylation (cg05575921) is the canonical DNA "
            "methylation biomarker of tobacco smoke exposure. Treated as an "
            "exposure-intensity correlate, not a direct risk multiplier."
        ),
        "carcinogen_class_multipliers": {},
        "mechanism": "exposure_biomarker",
    },
}


# ── Dataclasses ────────────────────────────────────────────────────────────


@dataclass
class MethylationMark:
    """A single gene-level methylation observation."""

    gene: str
    beta_value: float | None = None
    status: str | None = None  # "hyper", "hypo", or "normal"
    probe_id: str | None = None

    def resolved_status(
        self,
        *,
        hyper_threshold: float = HYPER_BETA_THRESHOLD,
        hypo_threshold: float = HYPO_BETA_THRESHOLD,
    ) -> str:
        """Resolve the effective hyper/hypo/normal call for this mark."""
        if self.status:
            return self.status.lower()
        if self.beta_value is None:
            return "unknown"
        if self.beta_value >= hyper_threshold:
            return "hyper"
        if self.beta_value <= hypo_threshold:
            return "hypo"
        return "normal"


@dataclass
class MethylationStatus:
    """Collected methylation marks for a single individual."""

    present: bool
    marks: list[MethylationMark] = field(default_factory=list)

    def detrimental_genes(
        self,
        *,
        hyper_threshold: float = HYPER_BETA_THRESHOLD,
        hypo_threshold: float = HYPO_BETA_THRESHOLD,
    ) -> list[str]:
        """Return genes whose methylation matches a curated trigger direction."""
        seen: list[str] = []
        for mark in self.marks:
            gene = mark.gene.upper()
            effect = EPIGENETIC_GENE_EFFECTS.get(gene)
            if effect is None:
                continue
            observed = mark.resolved_status(
                hyper_threshold=hyper_threshold,
                hypo_threshold=hypo_threshold,
            )
            if observed == effect["direction"] and gene not in seen:
                seen.append(gene)
        return seen


@dataclass
class EpigeneticEffect:
    """Result of applying the epigenomic modifier layer."""

    methylation_positive: bool
    aggregate_modifier: float
    per_class_multipliers: dict[str, float]
    affected_genes: list[str]
    affected_pathways: dict[str, list[str]]
    mechanism_summary: str
    notes: list[str] = field(default_factory=list)


# ── Parsing helpers ────────────────────────────────────────────────────────


def parse_methylation_calls(
    calls: Iterable[dict[str, Any]] | None,
) -> MethylationStatus:
    """Build a :class:`MethylationStatus` from a list of methylation records.

    Each call may provide ``gene`` plus either a numeric ``beta_value`` / ``beta``
    or a string ``status`` ("hyper"/"hypo"/"normal"). Calls for genes not in
    :data:`EPIGENETIC_GENE_EFFECTS` are retained for transparency but do not
    drive effect computation.
    """
    if not calls:
        return MethylationStatus(present=False, marks=[])

    marks: list[MethylationMark] = []
    for raw in calls:
        gene = str(raw.get("gene", "")).upper()
        if not gene:
            continue
        beta_raw = raw.get("beta_value", raw.get("beta"))
        try:
            beta = float(beta_raw) if beta_raw is not None else None
        except (TypeError, ValueError):
            beta = None
        status_val = raw.get("status")
        marks.append(
            MethylationMark(
                gene=gene,
                beta_value=round(beta, 4) if beta is not None else None,
                status=str(status_val).lower() if status_val else None,
                probe_id=str(raw["probe_id"]) if raw.get("probe_id") else None,
            )
        )

    status = MethylationStatus(present=False, marks=marks)
    status.present = bool(status.detrimental_genes())
    return status


def derive_methylation_status(
    marks: Iterable[dict[str, Any]] | None = None,
    *,
    is_present: bool | None = None,
) -> MethylationStatus:
    """Build MethylationStatus from calls or a binary flag (mirrors CHIP API)."""
    status = parse_methylation_calls(marks)
    if is_present is not None and not status.marks:
        status = MethylationStatus(present=bool(is_present), marks=[])
    return status


# ── Core modifier computation ──────────────────────────────────────────────


def compute_epigenetic_modifier(
    status: MethylationStatus | None,
    *,
    aggregate_modifier: float = AGGREGATE_METHYLATION_MODIFIER,
    hyper_threshold: float = HYPER_BETA_THRESHOLD,
    hypo_threshold: float = HYPO_BETA_THRESHOLD,
) -> EpigeneticEffect:
    """Translate a methylation status into per-carcinogen-class multipliers.

    Multipliers combine multiplicatively across genes, mirroring the CHIP
    modifier. A bioactivator silencing produces a multiplier < 1 (protective);
    a detox/repair silencing produces a multiplier > 1 (harmful); the
    aggregate of both directions is a single net per-class multiplier.
    """
    if status is None or not status.present:
        return EpigeneticEffect(
            methylation_positive=False,
            aggregate_modifier=1.0,
            per_class_multipliers={},
            affected_genes=[],
            affected_pathways={},
            mechanism_summary="No detrimental methylation calls; identity modifier applied.",
        )

    per_class: dict[str, float] = {}
    affected_pathways: dict[str, list[str]] = {}
    mechanisms: list[str] = []
    notes: list[str] = []
    detrimental_count = 0

    for gene in status.detrimental_genes(
        hyper_threshold=hyper_threshold,
        hypo_threshold=hypo_threshold,
    ):
        effect = EPIGENETIC_GENE_EFFECTS.get(gene)
        if effect is None:
            notes.append(f"Gene {gene} flagged but not in curated effect table.")
            continue
        multipliers = effect.get("carcinogen_class_multipliers", {})
        if not multipliers:
            notes.append(
                f"{gene}: biomarker-only entry ({effect.get('mechanism', 'n/a')}); "
                "no multiplier applied."
            )
            continue
        detrimental_count += 1
        for class_name, multiplier in multipliers.items():
            per_class[class_name] = round(
                per_class.get(class_name, 1.0) * float(multiplier), 4
            )
            pathway_list = affected_pathways.setdefault(class_name, [])
            if gene not in pathway_list:
                pathway_list.append(gene)
        mechanisms.append(f"{gene} ({effect['mechanism']})")

    if detrimental_count == 0:
        return EpigeneticEffect(
            methylation_positive=False,
            aggregate_modifier=1.0,
            per_class_multipliers={},
            affected_genes=[],
            affected_pathways={},
            mechanism_summary=(
                "Methylation marks detected but none map to curated "
                "risk-modifying entries."
            ),
            notes=notes,
        )

    has_harmful = any(m > 1.0 for m in per_class.values())
    aggregate = round(float(aggregate_modifier), 4) if has_harmful else 1.0

    summary_bits: list[str] = []
    if mechanisms:
        summary_bits.append(f"Active marks: {', '.join(mechanisms)}.")
    summary_bits.append(
        f"Aggregate modifier: {aggregate:.2f}x applied to germline risk."
    )
    if per_class:
        top = max(per_class.items(), key=lambda item: item[1])
        bottom = min(per_class.items(), key=lambda item: item[1])
        summary_bits.append(
            f"Largest amplification: {top[0]} ({top[1]:.2f}x); "
            f"largest protection: {bottom[0]} ({bottom[1]:.2f}x)."
        )

    return EpigeneticEffect(
        methylation_positive=True,
        aggregate_modifier=aggregate,
        per_class_multipliers=per_class,
        affected_genes=list(status.detrimental_genes(
            hyper_threshold=hyper_threshold,
            hypo_threshold=hypo_threshold,
        )),
        affected_pathways=affected_pathways,
        mechanism_summary=" ".join(summary_bits),
        notes=notes,
    )


# ── Applying the modifier to existing risk outputs ─────────────────────────


def apply_methylation_to_risk_score(
    risk_score: float,
    effect: EpigeneticEffect,
    carcinogen_class: str | None = None,
) -> float:
    """Return a methylation-adjusted risk score for one carcinogen class."""
    if not effect.methylation_positive:
        return round(float(risk_score), 4)

    adjusted = float(risk_score) * float(effect.aggregate_modifier)
    if carcinogen_class is not None:
        class_multiplier = effect.per_class_multipliers.get(carcinogen_class, 1.0)
        adjusted *= float(class_multiplier)
    return round(adjusted, 4)


def apply_methylation_to_risk_map(
    risks: dict[str, float],
    effect: EpigeneticEffect,
) -> dict[str, float]:
    """Apply :func:`apply_methylation_to_risk_score` across a {class: score} map."""
    return {
        class_name: apply_methylation_to_risk_score(score, effect, class_name)
        for class_name, score in risks.items()
    }


# ── JSON-friendly serialization ────────────────────────────────────────────


def epigenetic_effect_to_dict(effect: EpigeneticEffect) -> dict[str, Any]:
    """Convert :class:`EpigeneticEffect` to a JSON-serializable dictionary."""
    payload = asdict(effect)
    payload["affected_pathways"] = {
        cls: list(genes) for cls, genes in effect.affected_pathways.items()
    }
    return payload


def methylation_status_to_dict(status: MethylationStatus) -> dict[str, Any]:
    """Convert :class:`MethylationStatus` to a JSON-serializable dictionary."""
    return {
        "present": status.present,
        "marks": [asdict(mark) for mark in status.marks],
        "detrimental_genes": status.detrimental_genes(),
    }


def get_epigenetic_gene_effects() -> dict[str, dict[str, Any]]:
    """Return a defensive copy of the curated methylation effect table."""
    return deepcopy(EPIGENETIC_GENE_EFFECTS)


__all__ = [
    "AGGREGATE_METHYLATION_MODIFIER",
    "EPIGENETIC_GENE_EFFECTS",
    "EpigeneticEffect",
    "HYPER_BETA_THRESHOLD",
    "HYPO_BETA_THRESHOLD",
    "MethylationMark",
    "MethylationStatus",
    "apply_methylation_to_risk_map",
    "apply_methylation_to_risk_score",
    "compute_epigenetic_modifier",
    "derive_methylation_status",
    "epigenetic_effect_to_dict",
    "get_epigenetic_gene_effects",
    "methylation_status_to_dict",
    "parse_methylation_calls",
]
