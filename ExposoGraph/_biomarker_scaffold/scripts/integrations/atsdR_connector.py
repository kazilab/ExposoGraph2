"""Curated toxicological evidence loader connector."""

from __future__ import annotations

import csv
from pathlib import Path

import yaml

from .source_base import EvidenceSource


class CuratedToxicologyConnector(EvidenceSource):
    def __init__(self, source_name: str, curated_path: str | Path | None = None, version: str | None = None):
        super().__init__(source_name=source_name, version=version)
        self.curated_path = Path(curated_path) if curated_path else None
        self._records = self._load_records()

    def _load_records(self) -> list[dict]:
        if not self.curated_path or not self.curated_path.exists():
            return []
        suffix = self.curated_path.suffix.lower()
        if suffix in {".yaml", ".yml"}:
            data = yaml.safe_load(self.curated_path.read_text(encoding="utf-8")) or []
            if isinstance(data, dict):
                data = data.get("records", data.get("entries", []))
            return data if isinstance(data, list) else []
        with self.curated_path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def search(self, query: str) -> list[dict]:
        q = query.lower()
        return [self.normalize(row) for row in self._records if q in " ".join(str(v).lower() for v in row.values())]

    def fetch(self, identifier: str) -> dict:
        for row in self._records:
            if identifier in {row.get("chemical"), row.get("identifier"), row.get("cas")}:
                return self.normalize(row)
        return {"identifier": identifier, "source_status": "missing_curated_record"}

    def normalize(self, raw: dict) -> dict:
        return {
            "carcinogen_class": raw.get("carcinogen_class"),
            "target_organs": raw.get("target_organs"),
            "cancer_sites": raw.get("cancer_sites"),
            "reference_dose": raw.get("reference_dose"),
            "slope_factor": raw.get("slope_factor"),
            "classification": raw.get("classification"),
            "source_reference": raw.get("source_reference") or raw.get("reference"),
            "source_status": raw.get("source_status") or "candidate_toxicology_evidence",
            "raw": raw,
        }
