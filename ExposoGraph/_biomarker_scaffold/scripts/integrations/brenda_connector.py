"""BRENDA curated export connector for enzyme-substrate kinetic evidence."""

from __future__ import annotations

import csv
from pathlib import Path

from .source_base import EvidenceSource


class BrendaConnector(EvidenceSource):
    def __init__(self, curated_path: str | Path | None = None):
        super().__init__(source_name="BRENDA", version=None)
        self.curated_path = Path(curated_path) if curated_path else None
        self._records = self._load_records()

    def _load_records(self) -> list[dict]:
        if not self.curated_path or not self.curated_path.exists():
            return []
        with self.curated_path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def search(self, query: str) -> list[dict]:
        q = query.lower()
        return [self.normalize(row) for row in self._records if q in " ".join(str(v).lower() for v in row.values())]

    def fetch(self, identifier: str) -> dict:
        for row in self._records:
            if identifier in {row.get("enzyme"), row.get("EC"), row.get("substrate"), row.get("id")}:
                return self.normalize(row)
        return {"identifier": identifier, "source_status": "missing_curated_record"}

    def normalize(self, raw: dict) -> dict:
        km = raw.get("Km_uM") or raw.get("km_um") or raw.get("Km")
        confidence = raw.get("confidence", 0.5)
        return {
            "enzyme": raw.get("enzyme"),
            "substrate": raw.get("substrate"),
            "Km_uM": float(km) if km not in (None, "") else None,
            "organism": raw.get("organism"),
            "tissue": raw.get("tissue") or raw.get("source"),
            "reference": raw.get("reference"),
            "source_status": "candidate_kinetic_evidence",
            "confidence": float(confidence) if confidence not in (None, "") else 0.5,
        }
