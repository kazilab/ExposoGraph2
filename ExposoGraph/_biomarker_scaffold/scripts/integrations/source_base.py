"""Base classes for optional external evidence connectors."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class EvidenceSource:
    source_name: str
    version: str | None = None

    def search(self, query: str) -> list[dict]:
        raise NotImplementedError

    def fetch(self, identifier: str) -> dict:
        raise NotImplementedError

    def normalize(self, raw: dict) -> dict:
        return dict(raw)
