from kg_interface import NodeSpec, EdgeSpec
from typing import Any, Dict, List


class NodeMetadataMatch(NodeSpec):
    def __init__(self, key: str, value: Any):
        self.key = key
        self.value = value

    def is_satisfied_by(self, id: str, attributes: Dict[str, Any]) -> bool:
        return attributes.get(self.key) == self.value


class EdgeMetadataMatch(EdgeSpec):
    def __init__(self, key: str, value: Any):
        self.key = key
        self.value = value

    def is_satisfied_by(
        self, source: str, target: str, attributes: Dict[str, Any]
    ) -> bool:
        return attributes.get(self.key) == self.value


class IdNameSpec(NodeSpec):
    """Filters nodes based on naming conventions or prefixes (e.g., 'Enzyme_')."""

    def __init__(self, prefix: str):
        self.prefix = prefix

    def is_satisfied_by(self, id: str, attributes: Dict[str, Any]) -> bool:
        return node_id.startswith(self.prefix)


class AndSpec(NodeSpec):
    """Combines multiple specifications. All must evaluate to True."""

    def __init__(self, specs: List[NodeSpec]):
        self.specs = specs

    def is_satisfied_by(self, node_id: str, attributes: Dict[str, Any]) -> bool:
        return all(spec.is_satisfied_by(node_id, attributes) for spec in self.specs)
