"""Runtime configuration helpers for app deployment modes."""

from __future__ import annotations

import os
from enum import Enum
from typing import Mapping


class AppMode(str, Enum):
    STATELESS = "stateless"
    LOCAL = "local"


class LLMProvider(str, Enum):
    OPENAI = "openai"
    OLLAMA = "ollama"


class GraphMode(str, Enum):
    EXPLORATORY = "exploratory"
    STRICT = "strict"


class GraphVisibility(str, Enum):
    ALL = "all"
    VALIDATED_ONLY = "validated_only"
    EXPLORATORY_ONLY = "exploratory_only"


def normalize_app_mode(value: str | None) -> AppMode:
    """Normalize environment/config values into a supported app mode."""
    if value is None:
        return AppMode.STATELESS

    cleaned = value.strip().lower()
    aliases = {
        "stateless": AppMode.STATELESS,
        "public": AppMode.STATELESS,
        "web": AppMode.STATELESS,
        "local": AppMode.LOCAL,
        "curation": AppMode.LOCAL,
        "persistent": AppMode.LOCAL,
    }
    return aliases.get(cleaned, AppMode.STATELESS)


def normalize_graph_mode(value: str | None) -> GraphMode:
    """Normalize graph-validation mode values."""
    if value is None:
        return GraphMode.EXPLORATORY

    cleaned = value.strip().lower()
    aliases = {
        "exploratory": GraphMode.EXPLORATORY,
        "draft": GraphMode.EXPLORATORY,
        "flexible": GraphMode.EXPLORATORY,
        "strict": GraphMode.STRICT,
        "validated": GraphMode.STRICT,
        "canonical": GraphMode.STRICT,
    }
    return aliases.get(cleaned, GraphMode.EXPLORATORY)


def normalize_graph_visibility(value: str | None) -> GraphVisibility:
    """Normalize graph-view visibility mode values."""
    if value is None:
        return GraphVisibility.ALL

    cleaned = value.strip().lower()
    aliases = {
        "all": GraphVisibility.ALL,
        "full": GraphVisibility.ALL,
        "validated": GraphVisibility.VALIDATED_ONLY,
        "validated_only": GraphVisibility.VALIDATED_ONLY,
        "canonical": GraphVisibility.VALIDATED_ONLY,
        "strict": GraphVisibility.VALIDATED_ONLY,
        "exploratory": GraphVisibility.EXPLORATORY_ONLY,
        "exploratory_only": GraphVisibility.EXPLORATORY_ONLY,
        "provisional": GraphVisibility.EXPLORATORY_ONLY,
    }
    return aliases.get(cleaned, GraphVisibility.ALL)


def get_app_mode(env: Mapping[str, str] | None = None) -> AppMode:
    """Resolve the current app mode from the environment."""
    source = env if env is not None else os.environ
    return normalize_app_mode(source.get("ExposoGraph_MODE"))


def persistence_enabled(mode: AppMode | str) -> bool:
    """Whether server-side persistence features are allowed."""
    normalized = mode if isinstance(mode, AppMode) else normalize_app_mode(mode)
    return normalized == AppMode.LOCAL
