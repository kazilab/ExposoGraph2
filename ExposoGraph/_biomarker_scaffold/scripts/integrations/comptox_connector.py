"""EPA CompTox connector skeleton with cache-first behavior."""

from __future__ import annotations

import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen

from .source_base import EvidenceSource


class CompToxConnector(EvidenceSource):
    def __init__(self, cache_dir: str | Path = "data/cache/comptox", dry_run: bool = False):
        super().__init__(source_name="EPA CompTox", version=None)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.dry_run = dry_run

    def search(self, query: str) -> list[dict]:
        if self.dry_run:
            return []
        cached = self.cache_dir / f"search_{query.lower().replace(' ', '_')}.json"
        if cached.exists():
            data = json.loads(cached.read_text(encoding="utf-8"))
            return data if isinstance(data, list) else data.get("results", [])
        # Optional live API lookup when cache is missing.
        url = f"https://comptox.epa.gov/dashboard-api/ccdapp1/search/chemical/equal/{quote(query)}"
        try:
            data = self._get_json(url)
        except RuntimeError:
            return []
        payload = data if isinstance(data, list) else data.get("results", [])
        cached.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        return payload

    def fetch(self, identifier: str) -> dict:
        if self.dry_run:
            return {"comptox_dtxsid": identifier, "source_status": "dry_run"}
        cached = self.cache_dir / f"{identifier}.json"
        if cached.exists():
            return self.normalize(json.loads(cached.read_text(encoding="utf-8")))
        url = f"https://comptox.epa.gov/dashboard-api/ccdapp2/chemical-detail/search/by-dsstoxsid?id={quote(identifier)}"
        try:
            data = self._get_json(url)
        except RuntimeError:
            return {"comptox_dtxsid": identifier, "source_status": "missing_cached_record"}
        cached.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
        return self.normalize(data)

    def normalize(self, raw: dict) -> dict:
        return {
            "comptox_dtxsid": raw.get("dtxsid") or raw.get("DTXSID") or raw.get("comptox_dtxsid"),
            "comptox_dtxcid": raw.get("dtxcid") or raw.get("DTXCID"),
            "cas": raw.get("casrn") or raw.get("cas"),
            "preferred_name": raw.get("preferredName") or raw.get("preferred_name"),
            "source_status": "candidate_identity",
            "raw": raw,
        }

    def _get_json(self, url: str) -> dict:
        try:
            with urlopen(url, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            raise RuntimeError(f"CompTox request failed with HTTP {exc.code}: {url}") from exc
        except URLError as exc:
            raise RuntimeError(f"CompTox request failed: {url}: {exc}") from exc
