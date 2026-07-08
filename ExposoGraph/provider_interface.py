"""Typed local facade for packaged v2 data and public workflow entry points."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Mapping

from .models import KnowledgeGraph


def _json_safe(value: Any) -> Any:
    """Return a JSON-compatible copy of common project result objects."""
    if value is None or isinstance(value, str | int | bool):
        return value
    if isinstance(value, float):
        return value
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return _json_safe(asdict(value))
    if hasattr(value, "model_dump"):
        try:
            return _json_safe(value.model_dump(mode="json"))
        except TypeError:
            return _json_safe(value.model_dump())
    if hasattr(value, "to_dict"):
        return _json_safe(value.to_dict())
    if isinstance(value, Mapping):
        return {str(_json_safe(key)): _json_safe(item) for key, item in value.items()}
    if isinstance(value, list | tuple | set | frozenset):
        return [_json_safe(item) for item in value]
    return value


class LocalV2DataProvider:
    """Stable local accessors over bundled graph, JSON, and workflow data."""

    @property
    def package_root(self) -> Path:
        return Path(__file__).resolve().parent

    @property
    def data_dir(self) -> Path:
        return self.package_root / "data"

    @property
    def map_dir(self) -> Path:
        return self.package_root / "map"

    @property
    def reference_graph_path(self) -> Path:
        return self.map_dir / "graph-data.js"

    def runtime_paths(self) -> dict[str, Path]:
        """Return packaged local runtime paths used by v2 workflows."""
        return {
            "kinetic_parameters": self.data_dir / "kinetic_parameters.json",
            "interaction_parameters": self.data_dir / "interaction_parameters.json",
            "exposure_database": self.data_dir / "exposure_database.json",
            "proxy_flux_parameters": self.data_dir / "proxy_flux_parameters.json",
            "proxy_flux_provenance": self.data_dir / "proxy_flux_provenance.json",
            "parameter_provenance": self.data_dir / "parameter_provenance.json",
            "tissue_expression_data": self.data_dir / "tissue_expression_data.json",
            "biomarker_mapping": self.data_dir / "biomarker_mapping.json",
            "reference_graph": self.reference_graph_path,
        }

    def build_reference_graph(self) -> KnowledgeGraph:
        from .reference_data import build_reference_graph

        return build_reference_graph()

    def build_reference_engine(self) -> Any:
        from .reference_data import build_reference_engine

        return build_reference_engine()

    def get_activity_scores(self, gene: str) -> list[dict[str, Any]] | None:
        from .reference_data import get_activity_scores

        return get_activity_scores(gene)

    def get_activity_score_metadata(self, gene: str) -> dict[str, object] | None:
        from .reference_data import get_activity_score_metadata

        return get_activity_score_metadata(gene)

    def compute_module3_pathway_flux(self, *args: Any, **kwargs: Any) -> Any:
        from .flux_engine import compute_pathway_flux

        return compute_pathway_flux(*args, **kwargs)

    def compute_module5_interaction_matrix(self, *args: Any, **kwargs: Any) -> Any:
        from .interaction_engine import compute_interaction_matrix

        return compute_interaction_matrix(*args, **kwargs)

    def json_safe(self, value: Any) -> Any:
        return _json_safe(value)


def get_default_v2_provider() -> LocalV2DataProvider:
    """Return the default packaged local provider facade."""
    return LocalV2DataProvider()

