"""CTD (Comparative Toxicogenomics Database) chemical-gene interaction client.

Queries the CTD public data via their batch query API to retrieve
chemical-gene interactions relevant to carcinogen metabolism.
"""

from __future__ import annotations

import csv
import io
import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

_BASE_URL = "https://ctdbase.org/tools/batchQuery.go"


@dataclass
class ChemicalGeneInteraction:
    """A single chemical-gene interaction from CTD."""

    chemical_name: str
    chemical_id: str
    gene_symbol: str
    gene_id: str
    organism: str = ""
    interaction: str = ""
    pubmed_ids: list[str] = field(default_factory=list)


class CTDClient:
    """Client for querying CTD chemical-gene interactions.

    Parameters
    ----------
    base_url:
        Override the CTD batch query URL (useful for testing).
    timeout:
        HTTP request timeout in seconds.
    """

    def __init__(
        self,
        base_url: str = _BASE_URL,
        timeout: int = 60,
    ) -> None:
        try:
            import requests as _requests  # noqa: F401
        except ModuleNotFoundError as exc:  # pragma: no cover
            raise RuntimeError(
                "The 'requests' package is required for CTD lookups. "
                "Install with: pip install ExposoGraph[db]"
            ) from exc
        self.base_url = base_url
        self.timeout = timeout

    def get_chemical_gene_interactions(
        self,
        chemical_name: str,
        *,
        organism: str = "Homo sapiens",
    ) -> list[ChemicalGeneInteraction]:
        """Fetch chemical-gene interactions for a given chemical.

        Parameters
        ----------
        chemical_name:
            Chemical name to query (e.g. ``"Benzo(a)pyrene"``).
        organism:
            Organism filter. Defaults to ``"Homo sapiens"``.
        """
        import requests

        params = {
            "inputType": "chem",
            "inputTerms": chemical_name,
            "report": "cgixns",
            "format": "tsv",
        }
        resp = requests.get(self.base_url, params=params, timeout=self.timeout)
        resp.raise_for_status()

        return self._parse_interactions(resp.text, organism=organism)

    def _parse_interactions(
        self,
        tsv_text: str,
        *,
        organism: str = "Homo sapiens",
    ) -> list[ChemicalGeneInteraction]:
        """Parse CTD TSV response into interaction objects."""
        interactions: list[ChemicalGeneInteraction] = []

        # Skip comment lines starting with #
        lines = [line for line in tsv_text.splitlines() if not line.startswith("#")]
        if not lines:
            return interactions

        reader = csv.reader(io.StringIO("\n".join(lines)), delimiter="\t")
        for row in reader:
            if len(row) < 6:
                continue

            row_organism = row[4].strip() if len(row) > 4 else ""
            if organism and row_organism != organism:
                continue

            pmids = row[7].split("|") if len(row) > 7 and row[7] else []

            interactions.append(
                ChemicalGeneInteraction(
                    chemical_name=row[0].strip(),
                    chemical_id=row[1].strip(),
                    gene_symbol=row[2].strip(),
                    gene_id=row[3].strip(),
                    organism=row_organism,
                    interaction=row[5].strip() if len(row) > 5 else "",
                    pubmed_ids=[p.strip() for p in pmids],
                )
            )

        return interactions

    def get_gene_interactions(
        self,
        gene_symbol: str,
        *,
        organism: str = "Homo sapiens",
    ) -> list[ChemicalGeneInteraction]:
        """Fetch chemical-gene interactions for a given gene.

        Parameters
        ----------
        gene_symbol:
            Gene symbol to query (e.g. ``"CYP1A1"``).
        organism:
            Organism filter. Defaults to ``"Homo sapiens"``.
        """
        import requests

        params = {
            "inputType": "gene",
            "inputTerms": gene_symbol,
            "report": "cgixns",
            "format": "tsv",
        }
        resp = requests.get(self.base_url, params=params, timeout=self.timeout)
        resp.raise_for_status()

        return self._parse_interactions(resp.text, organism=organism)
