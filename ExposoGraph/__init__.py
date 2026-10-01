"""ExposoGraph.

Build, curate, and export carcinogen metabolism knowledge graphs with
literature extraction, manual curation, and a quantitative risk stack
(Michaelis-Menten flux, tissue-specific subgraphs, multi-carcinogen
interaction, exposure integration, and population-scale simulation) for
14 IARC carcinogen classes including heavy metals, alcohol/acetaldehyde,
dioxins/PCBs, dietary nitrosamines, and chlorinated solvents.

"""

from ._version import __version__

__all__ = ["__version__"]
