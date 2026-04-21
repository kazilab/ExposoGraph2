"""Tests for ExposoGraph.llm_backend (mocked backends and retry logic)."""

import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from ExposoGraph.llm_backend import (
    OllamaBackend,
    OpenAIBackend,
    UsageRecord,
    _retry,
)
from ExposoGraph.models import KnowledgeGraph

# ── Sample data ───────────────────────────────────────────────────────────

SAMPLE_KG_DICT = {
    "nodes": [
        {"id": "BaP", "label": "Benzo[a]pyrene", "type": "Carcinogen", "group": "PAH"},
        {"id": "CYP1A1", "label": "CYP1A1", "type": "Enzyme", "phase": "I", "role": "Activation"},
    ],
    "edges": [
        {"source": "CYP1A1", "target": "BaP", "type": "ACTIVATES"},
    ],
}


# ── UsageRecord ───────────────────────────────────────────────────────────


class TestUsageRecord:
    def test_defaults(self):
        rec = UsageRecord(provider="openai", model="gpt-4o")
        assert rec.prompt_tokens == 0
        assert rec.completion_tokens == 0
        assert rec.total_tokens == 0
        assert rec.duration_ms == 0
        assert rec.timestamp  # non-empty ISO string

    def test_custom_values(self):
        rec = UsageRecord(
            provider="ollama",
            model="llama3.1",
            prompt_tokens=100,
            completion_tokens=50,
            total_tokens=150,
            duration_ms=1234,
        )
        assert rec.total_tokens == 150
        assert rec.duration_ms == 1234


# ── Retry helper ──────────────────────────────────────────────────────────


class TestRetry:
    def test_success_first_try(self):
        fn = MagicMock(return_value="ok")
        assert _retry(fn, retryable=(ValueError,)) == "ok"
        fn.assert_called_once()

    def test_retries_on_retryable(self):
        fn = MagicMock(side_effect=[ValueError("err"), ValueError("err"), "ok"])
        result = _retry(fn, max_retries=3, base_delay=0.01, retryable=(ValueError,))
        assert result == "ok"
        assert fn.call_count == 3

    def test_raises_after_max_retries(self):
        fn = MagicMock(side_effect=ValueError("persistent"))
        with pytest.raises(ValueError, match="persistent"):
            _retry(fn, max_retries=2, base_delay=0.01, retryable=(ValueError,))
        assert fn.call_count == 3  # initial + 2 retries

    def test_non_retryable_raises_immediately(self):
        fn = MagicMock(side_effect=TypeError("not retryable"))
        with pytest.raises(TypeError, match="not retryable"):
            _retry(fn, max_retries=3, base_delay=0.01, retryable=(ValueError,))
        fn.assert_called_once()


# ── OpenAI backend ────────────────────────────────────────────────────────


class TestOpenAIBackend:
    @patch("openai.OpenAI")
    def test_structured_output_path(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client

        kg = KnowledgeGraph(**SAMPLE_KG_DICT)
        usage_ns = SimpleNamespace(prompt_tokens=10, completion_tokens=5, total_tokens=15)
        parsed_msg = SimpleNamespace(parsed=kg)
        choice = SimpleNamespace(message=parsed_msg)
        response = SimpleNamespace(choices=[choice], usage=usage_ns)
        mock_client.beta.chat.completions.parse.return_value = response

        backend = OpenAIBackend(api_key="test-key")
        data, usage = backend.extract_json("text", "system prompt", "gpt-4o")

        assert data["nodes"][0]["id"] == "BaP"
        assert usage.provider == "openai"
        assert usage.prompt_tokens == 10
        assert usage.total_tokens == 15
        call_kwargs = mock_client.beta.chat.completions.parse.call_args.kwargs
        response_format = call_kwargs["response_format"]
        assert response_format is not KnowledgeGraph

    @patch("openai.OpenAI")
    def test_fallback_to_json_mode(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client

        # Structured output fails
        mock_client.beta.chat.completions.parse.side_effect = AttributeError("no .parse")

        # JSON fallback succeeds
        usage_ns = SimpleNamespace(prompt_tokens=20, completion_tokens=10, total_tokens=30)
        json_msg = SimpleNamespace(content=json.dumps(SAMPLE_KG_DICT))
        choice = SimpleNamespace(message=json_msg)
        fallback = SimpleNamespace(choices=[choice], usage=usage_ns)
        mock_client.chat.completions.create.return_value = fallback

        backend = OpenAIBackend(api_key="test-key")
        data, usage = backend.extract_json("text", "system prompt", "gpt-4o")

        assert len(data["nodes"]) == 2
        assert usage.prompt_tokens == 20
        mock_client.chat.completions.create.assert_called_once()

    @patch("openai.OpenAI")
    def test_fallback_none_content(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client

        mock_client.beta.chat.completions.parse.side_effect = AttributeError("no .parse")
        json_msg = SimpleNamespace(content=None)
        choice = SimpleNamespace(message=json_msg)
        fallback = SimpleNamespace(choices=[choice], usage=None)
        mock_client.chat.completions.create.return_value = fallback

        backend = OpenAIBackend(api_key="test-key")
        data, usage = backend.extract_json("text", "system prompt", "gpt-4o")

        assert data == {}
        assert usage.prompt_tokens == 0

    @patch("openai.OpenAI")
    def test_tracks_duration(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client

        kg = KnowledgeGraph(**SAMPLE_KG_DICT)
        usage_ns = SimpleNamespace(prompt_tokens=0, completion_tokens=0, total_tokens=0)
        parsed_msg = SimpleNamespace(parsed=kg)
        choice = SimpleNamespace(message=parsed_msg)
        response = SimpleNamespace(choices=[choice], usage=usage_ns)
        mock_client.beta.chat.completions.parse.return_value = response

        backend = OpenAIBackend(api_key="test-key")
        _data, usage = backend.extract_json("text", "system prompt", "gpt-4o")

        assert usage.duration_ms >= 0


# ── Ollama backend ────────────────────────────────────────────────────────


class TestOllamaBackend:
    @patch("requests.post")
    def test_successful_extraction(self, mock_post):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "message": {"content": json.dumps(SAMPLE_KG_DICT)},
            "prompt_eval_count": 50,
            "eval_count": 25,
        }
        mock_resp.raise_for_status = MagicMock()
        mock_post.return_value = mock_resp

        backend = OllamaBackend(base_url="http://localhost:11434")
        data, usage = backend.extract_json("text", "system prompt", "llama3.1")

        assert data["nodes"][0]["id"] == "BaP"
        assert usage.provider == "ollama"
        assert usage.model == "llama3.1"
        assert usage.prompt_tokens == 50
        assert usage.completion_tokens == 25
        assert usage.total_tokens == 75

        # Verify request payload
        call_kwargs = mock_post.call_args
        payload = call_kwargs.kwargs["json"]
        assert payload["model"] == "llama3.1"
        assert payload["format"] == "json"
        assert payload["stream"] is False

    @patch("requests.post")
    def test_uses_env_base_url(self, mock_post):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "message": {"content": "{}"},
            "prompt_eval_count": 0,
            "eval_count": 0,
        }
        mock_resp.raise_for_status = MagicMock()
        mock_post.return_value = mock_resp

        with patch.dict("os.environ", {"OLLAMA_BASE_URL": "http://custom:1234"}):
            backend = OllamaBackend()
            backend.extract_json("text", "system prompt", "mistral")

        url = mock_post.call_args.args[0]
        assert url == "http://custom:1234/api/chat"

    @patch("requests.post")
    def test_default_base_url(self, mock_post):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "message": {"content": "{}"},
            "prompt_eval_count": 0,
            "eval_count": 0,
        }
        mock_resp.raise_for_status = MagicMock()
        mock_post.return_value = mock_resp

        with patch.dict("os.environ", {}, clear=True):
            backend = OllamaBackend()
            backend.extract_json("text", "system prompt", "llama3.1")

        url = mock_post.call_args.args[0]
        assert url == "http://localhost:11434/api/chat"

    @patch("requests.post")
    def test_empty_response_content(self, mock_post):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "message": {"content": "{}"},
        }
        mock_resp.raise_for_status = MagicMock()
        mock_post.return_value = mock_resp

        backend = OllamaBackend()
        data, usage = backend.extract_json("text", "system prompt", "llama3.1")

        assert data == {}
        assert usage.prompt_tokens == 0


# ── Integration: extract_graph_with_usage ─────────────────────────────────


class TestExtractGraphWithUsage:
    @patch("openai.OpenAI")
    def test_returns_kg_and_usage(self, mock_openai_cls):
        from ExposoGraph.llm_extractor import extract_graph_with_usage

        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client

        kg = KnowledgeGraph(**SAMPLE_KG_DICT)
        usage_ns = SimpleNamespace(prompt_tokens=10, completion_tokens=5, total_tokens=15)
        parsed_msg = SimpleNamespace(parsed=kg)
        choice = SimpleNamespace(message=parsed_msg)
        response = SimpleNamespace(choices=[choice], usage=usage_ns)
        mock_client.beta.chat.completions.parse.return_value = response

        result_kg, usage = extract_graph_with_usage("text", api_key="test-key")

        assert isinstance(result_kg, KnowledgeGraph)
        assert len(result_kg.nodes) == 2
        assert isinstance(usage, UsageRecord)
        assert usage.total_tokens == 15

    def test_with_custom_backend(self):
        from ExposoGraph.llm_extractor import extract_graph_with_usage

        mock_backend = MagicMock()
        mock_backend.extract_json.return_value = (
            SAMPLE_KG_DICT,
            UsageRecord(provider="test", model="test-model", total_tokens=42),
        )

        result_kg, usage = extract_graph_with_usage(
            "text", backend=mock_backend, model="test-model",
        )

        assert len(result_kg.nodes) == 2
        assert usage.provider == "test"
        assert usage.total_tokens == 42
        mock_backend.extract_json.assert_called_once()
