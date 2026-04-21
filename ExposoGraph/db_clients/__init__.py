"""Public database clients for enriching the knowledge graph.

Submodules:
- :mod:`.kegg` — KEGG REST API client for pathway and enzyme lookups.
- :mod:`.ctd` — CTD chemical-gene interaction parser.
- :mod:`.iarc` — Bundled IARC carcinogen classification data.
"""

from .ctd import CTDClient
from .iarc import IARCClassifier, IARCGroup
from .kegg import KEGGClient

__all__ = [
    "CTDClient",
    "IARCClassifier",
    "IARCGroup",
    "KEGGClient",
]
