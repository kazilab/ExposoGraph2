from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

# --- Specifications (The "What") ---


class NodeSpec(ABC):
    @abstractmethod
    def is_satisfied_by(self, node_id: str, attributes: Dict[str, Any]) -> bool:
        pass


class EdgeSpec(ABC):
    @abstractmethod
    def is_satisfied_by(
        self, source: str, target: str, attributes: Dict[str, Any]
    ) -> bool:
        pass


@dataclass
class SubgraphQuery:
    node_spec: Optional[NodeSpec] = None
    edge_spec: Optional[EdgeSpec] = None
    include_neighbors: bool = False


# --- Unified Graph Interface ---


class KnowledgeGraph(ABC):
    @abstractmethod
    def get_subgraph(self, query: SubgraphQuery) -> "KnowledgeGraph":
        pass

    @abstractmethod
    def get_node_metadata(self, node_id: str) -> Dict[str, Any]:
        """Retrieves all metadata attributes for a specific node ID."""
        pass

    @abstractmethod
    def get_edge_metadata(self, source_id: str, target_id: str) -> Dict[str, Any]:
        pass

    @abstractmethod
    def get_all_nodes(self) -> List[Tuple[str, Dict[str, Any]]]:
        pass
