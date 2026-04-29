"""PubChem identity enrichment connector with cache and dry-run support."""

from __future__ import annotations

import json
from pathlib import Path
from urllib.parse import quote
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from .source_base import EvidenceSource


class PubChemConnector(EvidenceSource):
    def __init__(self, cache_dir: str | Path = "data/cache/pubchem", dry_run: bool = False):
        super().__init__(source_name="PubChem", version=None)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.dry_run = dry_run

    def search(self, query: str) -> list[dict]:
        if self.dry_run:
            return []
        encoded = quote(query)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/JSON"
        data = self._get_json(url, self.cache_dir / f"search_{encoded}.json")
        return [{"cid": cid, "query": query} for cid in data.get("IdentifierList", {}).get("CID", [])]

    def fetch(self, identifier: str) -> dict:
        if self.dry_run:
            return {"pubchem_cid": identifier, "source_status": "dry_run"}
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{quote(identifier)}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey/JSON"
        raw = self._get_json(url, self.cache_dir / f"cid_{identifier}.json")
        return self.normalize(raw)

    def normalize(self, raw: dict) -> dict:
        props = raw.get("PropertyTable", {}).get("Properties", [{}])[0]
        return {
            "pubchem_cid": str(props.get("CID")) if props.get("CID") is not None else None,
            "formula": props.get("MolecularFormula"),
            "mw_g_mol": props.get("MolecularWeight"),
            "smiles": props.get("CanonicalSMILES"),
            "inchikey": props.get("InChIKey"),
            "source_status": "candidate_identity",
        }

    def _get_json(self, url: str, cache_path: Path) -> dict:
        if cache_path.exists():
            return json.loads(cache_path.read_text(encoding="utf-8"))
        try:
            with urlopen(url, timeout=60) as response:
                data = json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            raise RuntimeError(f"PubChem request failed with HTTP {exc.code}: {url}") from exc
        except URLError as exc:
            raise RuntimeError(f"PubChem request failed: {url}: {exc}") from exc
        cache_path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
        return data
