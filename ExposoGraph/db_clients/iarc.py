"""Bundled IARC carcinogen classification data.

Provides a static lookup of IARC monograph classifications for common
carcinogens relevant to the carcinogen metabolism knowledge graph.
No external API calls are required — data is embedded as a Python dict.
"""

from __future__ import annotations

from enum import Enum
from typing import Optional


class IARCGroup(str, Enum):
    """IARC carcinogen classification groups."""

    GROUP_1 = "Group 1"
    GROUP_2A = "Group 2A"
    GROUP_2B = "Group 2B"
    GROUP_3 = "Group 3"


# ── Static classification data ────────────────────────────────────────────
# Source: IARC Monographs on the Identification of Carcinogenic Hazards to Humans
# Last updated from IARC Monographs Volumes 1-135

_IARC_CLASSIFICATIONS: dict[str, dict[str, str]] = {
    # Group 1 — Carcinogenic to humans
    "Benzo[a]pyrene": {"group": "Group 1", "cas": "50-32-8", "category": "PAH", "monograph": "Vol. 92 (2010); Vol. 100F (2012)"},
    "BaP": {"group": "Group 1", "cas": "50-32-8", "category": "PAH", "monograph": "Vol. 92 (2010); Vol. 100F (2012)"},
    "Aflatoxin B1": {"group": "Group 1", "cas": "1162-65-8", "category": "Mycotoxin", "monograph": "Vol. 82 (2002); Vol. 100F (2012)"},
    "AFB1": {"group": "Group 1", "cas": "1162-65-8", "category": "Mycotoxin", "monograph": "Vol. 82 (2002); Vol. 100F (2012)"},
    "Benzene": {"group": "Group 1", "cas": "71-43-2", "category": "Solvent", "monograph": "Vol. 120 (2018)"},
    "Formaldehyde": {"group": "Group 1", "cas": "50-00-0", "category": "Aldehyde", "monograph": "Vol. 100F (2012)"},
    "Vinyl chloride": {"group": "Group 1", "cas": "75-01-4", "category": "Halogenated", "monograph": "Vol. 100F (2012)"},
    "2-Naphthylamine": {"group": "Group 1", "cas": "91-59-8", "category": "Aromatic_Amine", "monograph": "Vol. 100F (2012)"},
    "4-Aminobiphenyl": {"group": "Group 1", "cas": "92-67-1", "category": "Aromatic_Amine", "monograph": "Vol. 100F (2012)"},
    "Benzidine": {"group": "Group 1", "cas": "92-87-5", "category": "Aromatic_Amine", "monograph": "Vol. 100F (2012)"},
    "Cyclophosphamide": {"group": "Group 1", "cas": "50-18-0", "category": "Alkylating", "monograph": "Vol. 100A (2012)"},
    "Ethylene oxide": {"group": "Group 1", "cas": "75-21-8", "category": "Alkylating", "monograph": "Vol. 100F (2012)"},
    "Tobacco smoke": {"group": "Group 1", "cas": "", "category": "Mixture", "monograph": "Vol. 100E (2012)"},
    "Aristolochic acid": {"group": "Group 1", "cas": "313-67-7", "category": "Plant_Toxin", "monograph": "Vol. 100A (2012)"},
    "Asbestos": {"group": "Group 1", "cas": "", "category": "Mineral", "monograph": "Vol. 100C (2012)"},
    "Acetaldehyde": {"group": "Group 1", "cas": "75-07-0", "category": "Aldehyde", "monograph": "Vol. 100E (2012)"},
    "2,3,7,8-Tetrachlorodibenzo-para-dioxin": {"group": "Group 1", "cas": "1746-01-6", "category": "Dioxin", "monograph": "Vol. 100F (2012)"},
    "TCDD": {"group": "Group 1", "cas": "1746-01-6", "category": "Dioxin", "monograph": "Vol. 100F (2012)"},
    "PCB 126": {"group": "Group 1", "cas": "57465-28-8", "category": "PCB", "monograph": "Vol. 107 (2016)"},
    "PCB_126": {"group": "Group 1", "cas": "57465-28-8", "category": "PCB", "monograph": "Vol. 107 (2016)"},
    "Trichloroethylene": {"group": "Group 1", "cas": "79-01-6", "category": "Chlorinated_Solvent", "monograph": "Vol. 106 (2014)"},
    "TCE": {"group": "Group 1", "cas": "79-01-6", "category": "Chlorinated_Solvent", "monograph": "Vol. 106 (2014)"},
    # Group 2A — Probably carcinogenic
    "Acrylamide": {"group": "Group 2A", "cas": "79-06-1", "category": "Alkylating", "monograph": "Vol. 60 (1994)"},
    "Benz[a]anthracene": {"group": "Group 2A", "cas": "56-55-3", "category": "PAH", "monograph": "Vol. 92 (2010)"},
    "Dibenz[a,h]anthracene": {"group": "Group 2A", "cas": "53-70-3", "category": "PAH", "monograph": "Vol. 92 (2010)"},
    "N-Nitrosodiethylamine": {"group": "Group 2A", "cas": "55-18-5", "category": "Nitrosamine", "monograph": "Suppl. 7 (1987)"},
    "NDEA": {"group": "Group 2A", "cas": "55-18-5", "category": "Nitrosamine", "monograph": "Suppl. 7 (1987)"},
    "N-Nitrosodimethylamine": {"group": "Group 2A", "cas": "62-75-9", "category": "Nitrosamine", "monograph": "Vol. 17 (1978); Suppl. 7 (1987)"},
    "NDMA": {"group": "Group 2A", "cas": "62-75-9", "category": "Nitrosamine", "monograph": "Vol. 17 (1978); Suppl. 7 (1987)"},
    "Styrene": {"group": "Group 2A", "cas": "100-42-5", "category": "Solvent", "monograph": "Vol. 121 (2019)"},
    "Tetrachloroethylene": {"group": "Group 2A", "cas": "127-18-4", "category": "Chlorinated_Solvent", "monograph": "Vol. 106 (2014)"},
    "PCE": {"group": "Group 2A", "cas": "127-18-4", "category": "Chlorinated_Solvent", "monograph": "Vol. 106 (2014)"},
    # Group 2B — Possibly carcinogenic
    "Chrysene": {"group": "Group 2B", "cas": "218-01-9", "category": "PAH", "monograph": "Vol. 92 (2010)"},
    "Benzo[b]fluoranthene": {"group": "Group 2B", "cas": "205-99-2", "category": "PAH", "monograph": "Vol. 92 (2010)"},
    "Naphthalene": {"group": "Group 2B", "cas": "91-20-3", "category": "PAH", "monograph": "Vol. 82 (2002)"},
    "Styrene-7,8-oxide": {"group": "Group 2B", "cas": "96-09-3", "category": "Epoxide", "monograph": "Vol. 121 (2019)"},
    "4-(Methylnitrosamino)-1-(3-pyridyl)-1-butanone": {"group": "Group 1", "cas": "64091-91-4", "category": "Nitrosamine", "monograph": "Vol. 89 (2007); Vol. 100E (2012)"},
    "NNK": {"group": "Group 1", "cas": "64091-91-4", "category": "Nitrosamine", "monograph": "Vol. 89 (2007); Vol. 100E (2012)"},
    "PhIP": {"group": "Group 2B", "cas": "105650-23-5", "category": "HCA", "monograph": "Vol. 56 (1993)"},
    "MeIQx": {"group": "Group 2B", "cas": "77500-04-0", "category": "HCA", "monograph": "Vol. 56 (1993)"},
    # Group 3 — Not classifiable
    "Pyrene": {"group": "Group 3", "cas": "129-00-0", "category": "PAH", "monograph": "Vol. 92 (2010)"},
    "Fluoranthene": {"group": "Group 3", "cas": "206-44-0", "category": "PAH", "monograph": "Vol. 92 (2010)"},
    "Anthracene": {"group": "Group 3", "cas": "120-12-7", "category": "PAH", "monograph": "Vol. 92 (2010)"},
}


class IARCClassifier:
    """Look up IARC classifications from the bundled static dataset.

    Example
    -------
    >>> clf = IARCClassifier()
    >>> clf.classify("Benzo[a]pyrene")
    IARCGroup.GROUP_1
    """

    def __init__(self, extra: Optional[dict[str, dict[str, str]]] = None) -> None:
        self._data = dict(_IARC_CLASSIFICATIONS)
        if extra:
            self._data.update(extra)
        self._normalized_index = {
            self._normalize_name(name): name
            for name in self._data
        }

    @staticmethod
    def _normalize_name(chemical_name: str) -> str:
        return "".join(ch.lower() for ch in chemical_name if ch.isalnum())

    def _resolve_name(self, chemical_name: str) -> Optional[str]:
        if chemical_name in self._data:
            return chemical_name
        return self._normalized_index.get(self._normalize_name(chemical_name))

    def classify(self, chemical_name: str) -> Optional[IARCGroup]:
        """Return the IARC group for a chemical, or ``None`` if not found."""
        resolved_name = self._resolve_name(chemical_name)
        if resolved_name is None:
            return None
        entry = self._data.get(resolved_name)
        if entry is None:
            return None
        group_str = entry["group"]
        for member in IARCGroup:
            if member.value == group_str:
                return member
        return None

    def get_entry(self, chemical_name: str) -> Optional[dict[str, str]]:
        """Return the full classification entry (group, CAS, category)."""
        resolved_name = self._resolve_name(chemical_name)
        if resolved_name is None:
            return None
        return self._data.get(resolved_name)

    def list_by_group(self, group: IARCGroup) -> list[str]:
        """Return all chemical names in a given IARC group."""
        return [
            name
            for name, entry in self._data.items()
            if entry["group"] == group.value
        ]

    def list_by_category(self, category: str) -> list[str]:
        """Return all chemical names in a given category (e.g. ``'PAH'``)."""
        return [
            name
            for name, entry in self._data.items()
            if entry.get("category") == category
        ]

    @property
    def all_chemicals(self) -> list[str]:
        """Return all known chemical names."""
        return list(self._data.keys())
