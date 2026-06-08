"""Provider layer for typed interaction parameter data.

The JSON provider reads only bundled local JSON files and preserves local values
as typed records. It does not infer reaction-role semantics or risk direction
from product hazard metadata.
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from .interaction_schema import (
    AssumptionWarning,
    CompetitiveInteraction,
    EvidenceGrade,
    EvidenceRecord,
    GSHConsumer,
    InductionRule,
    KineticParameterSet,
    MetabolicReaction,
    ParameterUncertainty,
    ReactionRole,
    RiskDirectionIfFluxDecreases,
    TissueContext,
    enum_from_value,
)


class InteractionParameterProvider(ABC):
    """Abstract access layer for mechanism-resolved interaction parameters."""

    @abstractmethod
    def get_competitive_interactions(self, enzyme: str | None = None) -> list[CompetitiveInteraction]:
        raise NotImplementedError

    @abstractmethod
    def get_reactions_for_enzyme(self, enzyme: str) -> list[MetabolicReaction]:
        raise NotImplementedError

    @abstractmethod
    def get_reactions_for_carcinogen(self, carcinogen: str, tissue: str | None = None) -> list[MetabolicReaction]:
        raise NotImplementedError

    @abstractmethod
    def get_gsh_consumers(self, tissue: str | None = None) -> list[GSHConsumer]:
        raise NotImplementedError

    @abstractmethod
    def get_induction_rules(self, tissue: str | None = None) -> list[InductionRule]:
        raise NotImplementedError

    @abstractmethod
    def get_parameter_evidence(self, enzyme: str, substrate: str) -> EvidenceRecord | None:
        raise NotImplementedError


class JSONInteractionParameterProvider(InteractionParameterProvider):
    """Read typed records from bundled local interaction JSON files."""

    def __init__(self, data_dir: str | Path | None = None) -> None:
        self.data_dir = Path(data_dir) if data_dir is not None else Path(__file__).resolve().parent / "data"
        self.interaction_path = self.data_dir / "interaction_parameters.json"
        self.provenance_path = self.data_dir / "parameter_provenance.json"
        self._interaction_data = self._read_json(self.interaction_path)
        self._provenance_data = self._read_json(self.provenance_path)

    @staticmethod
    def _read_json(path: Path) -> dict[str, Any]:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    @property
    def _competitive_data(self) -> dict[str, Any]:
        return self._interaction_data.get("competitive_inhibition", {})

    @property
    def _provenance_pairs(self) -> dict[str, Any]:
        return self._provenance_data.get("pairs", {})

    def get_competitive_interactions(self, enzyme: str | None = None) -> list[CompetitiveInteraction]:
        interactions: list[CompetitiveInteraction] = []
        for enzyme_name, enzyme_data in self._iter_enzyme_data(enzyme):
            for substrate, substrate_data in self._iter_substrate_data(enzyme_data):
                interactions.append(self._build_competitive_interaction(enzyme_name, substrate, substrate_data))
        return interactions

    def get_reactions_for_enzyme(self, enzyme: str) -> list[MetabolicReaction]:
        return [self._interaction_to_reaction(item) for item in self.get_competitive_interactions(enzyme)]

    def get_reactions_for_carcinogen(self, carcinogen: str, tissue: str | None = None) -> list[MetabolicReaction]:
        del tissue  # Phase 3 carries tissue context but does not infer tissue semantics.
        target = carcinogen.lower()
        return [
            self._interaction_to_reaction(item)
            for item in self.get_competitive_interactions()
            if item.substrate.lower() == target
        ]

    def get_gsh_consumers(self, tissue: str | None = None) -> list[GSHConsumer]:
        consumers = self._interaction_data.get("gsh_depletion", {}).get("consumers", {})
        parsed: list[GSHConsumer] = []
        for name, data in consumers.items():
            tissue_context = TissueContext(tissue=data.get("tissue") if isinstance(data, dict) else None)
            if tissue and tissue_context.tissue and tissue_context.tissue.lower() != tissue.lower():
                continue
            parsed.append(
                GSHConsumer(
                    name=name,
                    enzyme=data.get("enzyme"),
                    substrate_class=data.get("substrate_class"),
                    gsh_per_umol_substrate=self._float_or_none(data.get("gsh_per_umol_substrate")),
                    tissue_context=tissue_context,
                    metadata=dict(data),
                )
            )
        return parsed

    def get_induction_rules(self, tissue: str | None = None) -> list[InductionRule]:
        induction_data = self._interaction_data.get("enzyme_induction", {})
        rules: list[InductionRule] = []
        for exposure_context, exposure_data in induction_data.items():
            if not isinstance(exposure_data, dict):
                continue
            for enzyme, rule_data in exposure_data.items():
                if enzyme.startswith("_") or not isinstance(rule_data, dict):
                    continue
                tissue_specificity = rule_data.get("tissue_specificity")
                if tissue and tissue_specificity and tissue.lower() not in str(tissue_specificity).lower():
                    continue
                rules.append(
                    InductionRule(
                        exposure_context=exposure_context,
                        enzyme=enzyme,
                        fold_induction=self._float_or_none(rule_data.get("fold_induction")),
                        range_min=self._float_or_none(rule_data.get("range_min")),
                        range_max=self._float_or_none(rule_data.get("range_max")),
                        mechanism=rule_data.get("mechanism"),
                        tissue_context=TissueContext(tissue=tissue_specificity),
                        metadata=dict(rule_data),
                    )
                )
        return rules

    def get_parameter_evidence(self, enzyme: str, substrate: str) -> EvidenceRecord | None:
        provenance = self._get_provenance(enzyme, substrate)
        if provenance is None:
            return None
        return self._evidence_from_provenance(enzyme, substrate, provenance)

    def _iter_enzyme_data(self, enzyme: str | None = None) -> list[tuple[str, dict[str, Any]]]:
        if enzyme is None:
            return [(name, data) for name, data in self._competitive_data.items() if isinstance(data, dict)]
        for name, data in self._competitive_data.items():
            if name.lower() == enzyme.lower() and isinstance(data, dict):
                return [(name, data)]
        return []

    @staticmethod
    def _iter_substrate_data(enzyme_data: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
        substrates = enzyme_data.get("substrates", {})
        return [(name, data) for name, data in substrates.items() if isinstance(data, dict)]

    def _build_competitive_interaction(
        self, enzyme: str, substrate: str, substrate_data: dict[str, Any]
    ) -> CompetitiveInteraction:
        evidence = self.get_parameter_evidence(enzyme, substrate)
        reaction_role = enum_from_value(
            ReactionRole,
            substrate_data.get("reaction_role"),
            ReactionRole.UNKNOWN,
        )
        risk_direction = enum_from_value(
            RiskDirectionIfFluxDecreases,
            substrate_data.get("risk_direction_if_flux_decreases"),
            RiskDirectionIfFluxDecreases.UNKNOWN,
        )
        warnings = self._unknown_warnings(reaction_role, risk_direction)
        return CompetitiveInteraction(
            enzyme=enzyme,
            substrate=substrate,
            kinetic_parameters=self._kinetic_parameters(substrate_data, evidence),
            reaction_role=reaction_role,
            risk_direction_if_flux_decreases=risk_direction,
            evidence=evidence,
            warnings=warnings,
            metadata={"source": "interaction_parameters.json", "raw": dict(substrate_data)},
        )

    def _interaction_to_reaction(self, item: CompetitiveInteraction) -> MetabolicReaction:
        product = item.kinetic_parameters.product if item.kinetic_parameters else None
        return MetabolicReaction(
            enzyme=item.enzyme,
            substrate=item.substrate,
            product=product,
            kinetic_parameters=item.kinetic_parameters,
            reaction_role=item.reaction_role,
            risk_direction_if_flux_decreases=item.risk_direction_if_flux_decreases,
            evidence=item.evidence,
            warnings=item.warnings,
            sme_notes=item.sme_notes,
            metadata=item.metadata,
        )

    def _kinetic_parameters(self, data: dict[str, Any], evidence: EvidenceRecord | None) -> KineticParameterSet:
        product = data.get("product")
        product_carcinogenic = data.get("product_carcinogenic")
        product_hazard = {
            "product": product,
            "product_carcinogenic": product_carcinogenic,
        }
        uncertainty = None
        if evidence and evidence.confidence:
            uncertainty = ParameterUncertainty(confidence=evidence.confidence, notes=evidence.notes)
        return KineticParameterSet(
            km_uM=self._float_or_none(data.get("Km_uM")),
            ki_uM=self._float_or_none(data.get("Ki_uM")),
            vmax_relative=self._float_or_none(data.get("Vmax_relative")),
            hill_coefficient=self._float_or_none(data.get("hill_coefficient")),
            relative_priority=self._int_or_none(data.get("relative_priority")),
            assumed_ki=data.get("assumed_ki"),
            product=product,
            product_hazard=product_hazard,
            evidence=evidence,
            uncertainty=uncertainty,
            metadata={
                "notes": data.get("notes"),
                "provenance_ref": data.get("provenance_ref"),
                "local_fields": dict(data),
            },
        )

    @staticmethod
    def _unknown_warnings(
        reaction_role: ReactionRole,
        risk_direction: RiskDirectionIfFluxDecreases,
    ) -> list[AssumptionWarning]:
        warnings: list[AssumptionWarning] = []
        if reaction_role is ReactionRole.UNKNOWN:
            warnings.append(
                AssumptionWarning(
                    code="reaction_role_unknown",
                    message="Reaction role is not explicitly curated in local data.",
                    field="reaction_role",
                )
            )
        if risk_direction is RiskDirectionIfFluxDecreases.UNKNOWN:
            warnings.append(
                AssumptionWarning(
                    code="risk_direction_unknown",
                    message="Risk direction if flux decreases is not explicitly curated in local data.",
                    field="risk_direction_if_flux_decreases",
                )
            )
        return warnings

    def _get_provenance(self, enzyme: str, substrate: str) -> dict[str, Any] | None:
        enzyme_pairs = self._provenance_pairs.get(enzyme, {})
        if substrate in enzyme_pairs:
            return enzyme_pairs[substrate]
        target = substrate.lower()
        for name, data in enzyme_pairs.items():
            if name.lower() == target:
                return data
        return None

    @staticmethod
    def _evidence_from_provenance(enzyme: str, substrate: str, provenance: dict[str, Any]) -> EvidenceRecord:
        ki_status = provenance.get("ki_status")
        grade = EvidenceGrade.CURATED if ki_status == "curated" else EvidenceGrade.INFERRED_FROM_LOCAL_METADATA
        source = provenance.get("ki_reference") or provenance.get("km_source") or provenance.get("vmax_source")
        confidence = provenance.get("ki_status") or provenance.get("km_confidence") or provenance.get("vmax_confidence")
        return EvidenceRecord(
            source=source,
            grade=grade,
            confidence=confidence,
            provenance_ref=f"parameter_provenance.json#pairs/{enzyme}/{substrate}",
            notes=provenance.get("notes"),
            metadata=dict(provenance),
        )

    @staticmethod
    def _float_or_none(value: Any) -> float | None:
        if value is None:
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _int_or_none(value: Any) -> int | None:
        if value is None:
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            return None


class KGInteractionParameterProvider(InteractionParameterProvider):
    """Scaffold for future KG-backed parameter population.

    Phase 3 must not traverse graph internals. These methods fail explicitly so
    callers cannot mistake the scaffold for an implemented provider.
    """

    _MESSAGE = "KGInteractionParameterProvider is a Phase 3 scaffold; KG traversal is reserved for later phases."

    def _not_implemented(self) -> None:
        raise NotImplementedError(self._MESSAGE)

    def get_competitive_interactions(self, enzyme: str | None = None) -> list[CompetitiveInteraction]:
        del enzyme
        self._not_implemented()

    def get_reactions_for_enzyme(self, enzyme: str) -> list[MetabolicReaction]:
        del enzyme
        self._not_implemented()

    def get_reactions_for_carcinogen(self, carcinogen: str, tissue: str | None = None) -> list[MetabolicReaction]:
        del carcinogen, tissue
        self._not_implemented()

    def get_gsh_consumers(self, tissue: str | None = None) -> list[GSHConsumer]:
        del tissue
        self._not_implemented()

    def get_induction_rules(self, tissue: str | None = None) -> list[InductionRule]:
        del tissue
        self._not_implemented()

    def get_parameter_evidence(self, enzyme: str, substrate: str) -> EvidenceRecord | None:
        del enzyme, substrate
        self._not_implemented()
