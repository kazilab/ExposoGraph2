"""Pluggable LLM backends for knowledge graph extraction.

Each backend accepts a text prompt and system prompt, returning raw JSON
that is validated into a :class:`KnowledgeGraph` by the caller.
"""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Protocol

from pydantic import BaseModel, ConfigDict, Field

logger = logging.getLogger(__name__)

# ── Usage tracking ────────────────────────────────────────────────────────


@dataclass
class UsageRecord:
    """Token usage and cost metadata from a single LLM call."""

    provider: str
    model: str
    prompt_tokens: int = 0
    completion_tokens: int = 0
    total_tokens: int = 0
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    duration_ms: int = 0


class _ExtractionNode(BaseModel):
    """Lightweight schema for raw LLM extraction output."""

    model_config = ConfigDict(extra="ignore")

    id: str
    label: str
    type: str
    detail: str = ""
    group: str | None = None
    iarc: str | None = None
    phase: str | None = None
    role: str | None = None
    reactivity: str | None = None
    source_db: str | None = None
    evidence: str | None = None
    pmid: str | None = None
    tissue: str | None = None
    variant: str | None = None
    phenotype: str | None = None
    activity_score: float | None = None
    tier: int | None = None


class _ExtractionEdge(BaseModel):
    model_config = ConfigDict(extra="ignore")

    source: str
    target: str
    type: str
    label: str | None = None
    carcinogen: str | None = None
    source_db: str | None = None
    evidence: str | None = None
    pmid: str | None = None
    tissue: str | None = None


class _ExtractionKnowledgeGraph(BaseModel):
    model_config = ConfigDict(extra="ignore")

    nodes: list[_ExtractionNode] = Field(default_factory=list)
    edges: list[_ExtractionEdge] = Field(default_factory=list)


# ── Backend protocol ──────────────────────────────────────────────────────


class LLMBackend(Protocol):
    """Protocol that all LLM backends must satisfy."""

    def extract_json(
        self,
        text: str,
        system_prompt: str,
        model: str,
    ) -> tuple[dict[str, Any], UsageRecord]:
        """Return parsed JSON dict and usage metadata."""
        ...  # pragma: no cover


# ── Retry helper ──────────────────────────────────────────────────────────


def _retry(
    fn: Any,
    *,
    max_retries: int = 3,
    base_delay: float = 1.0,
    retryable: tuple[type[Exception], ...] = (),
) -> Any:
    """Call *fn* with exponential backoff on retryable exceptions."""
    last_exc: Exception | None = None
    for attempt in range(max_retries + 1):
        try:
            return fn()
        except retryable as exc:
            last_exc = exc
            if attempt < max_retries:
                delay = base_delay * (2 ** attempt)
                logger.warning(
                    "Attempt %d/%d failed (%s), retrying in %.1fs…",
                    attempt + 1, max_retries + 1, exc, delay,
                )
                time.sleep(delay)
    raise last_exc  # type: ignore[misc]


# ── OpenAI backend ────────────────────────────────────────────────────────


class OpenAIBackend:
    """Backend using the OpenAI API (or any OpenAI-compatible endpoint)."""

    def __init__(
        self,
        api_key: str | None = None,
        base_url: str | None = None,
        max_retries: int = 3,
    ) -> None:
        try:
            from openai import OpenAI
        except ModuleNotFoundError as exc:  # pragma: no cover
            raise RuntimeError(
                "The 'openai' package is required for the OpenAI backend. "
                "Install with: pip install ExposoGraph[llm]"
            ) from exc

        self.client = OpenAI(
            api_key=api_key or os.environ.get("OPENAI_API_KEY"),
            base_url=base_url or os.environ.get("OPENAI_BASE_URL"),
        )
        self.max_retries = max_retries

    def extract_json(
        self,
        text: str,
        system_prompt: str,
        model: str,
    ) -> tuple[dict[str, Any], UsageRecord]:
        import openai as _openai

        messages: list[dict[str, str]] = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": text},
        ]

        retryable = (
            _openai.APITimeoutError,
            _openai.APIConnectionError,
            _openai.RateLimitError,
            _openai.InternalServerError,
        )

        start = time.monotonic()

        try:
            def _parse() -> Any:
                return self.client.beta.chat.completions.parse(
                    model=model,
                    messages=messages,  # type: ignore[arg-type]
                    response_format=_ExtractionKnowledgeGraph,
                )

            response = _retry(_parse, max_retries=self.max_retries, retryable=retryable)
            duration = int((time.monotonic() - start) * 1000)
            kg = response.choices[0].message.parsed
            if kg is not None:
                usage = response.usage
                return kg.model_dump(mode="json"), UsageRecord(
                    provider="openai",
                    model=model,
                    prompt_tokens=usage.prompt_tokens if usage else 0,
                    completion_tokens=usage.completion_tokens if usage else 0,
                    total_tokens=usage.total_tokens if usage else 0,
                    duration_ms=duration,
                )
        except (
            _openai.BadRequestError,
            _openai.APIResponseValidationError,
            _openai.LengthFinishReasonError,
            _openai.ContentFilterFinishReasonError,
            AttributeError,
            TypeError,
            ValueError,
        ) as exc:
            logger.warning("Structured output failed, falling back to JSON mode: %s", exc)

        # Fallback: JSON mode
        def _create() -> Any:
            return self.client.chat.completions.create(  # type: ignore[call-overload]
                model=model,
                messages=messages,
                response_format={"type": "json_object"},
            )

        fallback = _retry(_create, max_retries=self.max_retries, retryable=retryable)
        duration = int((time.monotonic() - start) * 1000)
        content = fallback.choices[0].message.content or "{}"
        usage = fallback.usage
        return json.loads(content), UsageRecord(
            provider="openai",
            model=model,
            prompt_tokens=usage.prompt_tokens if usage else 0,
            completion_tokens=usage.completion_tokens if usage else 0,
            total_tokens=usage.total_tokens if usage else 0,
            duration_ms=duration,
        )


# ── Ollama backend ────────────────────────────────────────────────────────


class OllamaBackend:
    """Backend using a local Ollama instance.

    Communicates via the Ollama HTTP API (``/api/chat``).  No extra
    dependencies beyond ``requests`` are required.
    """

    DEFAULT_BASE_URL = "http://localhost:11434"

    def __init__(
        self,
        base_url: str | None = None,
        max_retries: int = 3,
    ) -> None:
        try:
            import requests as _requests  # noqa: F401
        except ModuleNotFoundError as exc:  # pragma: no cover
            raise RuntimeError(
                "The 'requests' package is required for the Ollama backend. "
                "Install with: pip install requests"
            ) from exc
        self.base_url = (
            base_url
            or os.environ.get("OLLAMA_BASE_URL")
            or self.DEFAULT_BASE_URL
        )
        self.max_retries = max_retries

    def extract_json(
        self,
        text: str,
        system_prompt: str,
        model: str,
    ) -> tuple[dict[str, Any], UsageRecord]:
        import requests

        url = f"{self.base_url.rstrip('/')}/api/chat"
        payload = {
            "model": model,
            "messages": [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": text},
            ],
            "format": "json",
            "stream": False,
        }

        retryable = (
            requests.ConnectionError,
            requests.Timeout,
        )

        start = time.monotonic()

        def _call() -> requests.Response:
            resp = requests.post(url, json=payload, timeout=120)
            resp.raise_for_status()
            return resp

        response = _retry(_call, max_retries=self.max_retries, retryable=retryable)
        duration = int((time.monotonic() - start) * 1000)
        body = response.json()

        content = body.get("message", {}).get("content", "{}")
        data = json.loads(content)

        prompt_tokens = body.get("prompt_eval_count", 0)
        completion_tokens = body.get("eval_count", 0)

        return data, UsageRecord(
            provider="ollama",
            model=model,
            prompt_tokens=prompt_tokens,
            completion_tokens=completion_tokens,
            total_tokens=prompt_tokens + completion_tokens,
            duration_ms=duration,
        )
