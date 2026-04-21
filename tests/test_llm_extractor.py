"""Tests for ExposoGraph.llm_extractor (mocked OpenAI client)."""

import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from ExposoGraph.config import GraphMode
from ExposoGraph.llm_extractor import SYSTEM_PROMPT, extract_graph
from ExposoGraph.models import KnowledgeGraph, MatchStatus, RecordOrigin

# ── Fixtures ─────────────────────────────────────────────────────────────

SAMPLE_KG_DICT = {
    "nodes": [
        {
            "id": "BaP",
            "label": "Benzo[a]pyrene",
            "type": "Carcinogen",
            "group": "PAH",
            "iarc": "Group 1",
        },
        {"id": "CYP1A1", "label": "CYP1A1", "type": "Enzyme", "phase": "I", "role": "Activation"},
        {"id": "BPDE", "label": "BPDE", "type": "Metabolite", "reactivity": "High"},
    ],
    "edges": [
        {"source": "CYP1A1", "target": "BaP", "type": "ACTIVATES", "label": "epoxidation"},
        {"source": "CYP1A1", "target": "BPDE", "type": "ACTIVATES", "label": "second epoxidation"},
    ],
}

_USAGE = SimpleNamespace(prompt_tokens=10, completion_tokens=5, total_tokens=15)


class _StubBackend:
    def __init__(self, payload: dict):
        self.payload = payload

    def extract_json(self, text: str, system_prompt: str, model: str):
        return self.payload, _USAGE


def _make_structured_response(kg, usage=_USAGE):
    """Build a mock structured-output response."""
    message = SimpleNamespace(parsed=kg)
    choice = SimpleNamespace(message=message)
    return SimpleNamespace(choices=[choice], usage=usage)


def _make_json_response(data: dict, usage=_USAGE):
    """Build a mock JSON-mode response."""
    message = SimpleNamespace(content=json.dumps(data))
    choice = SimpleNamespace(message=message)
    return SimpleNamespace(choices=[choice], usage=usage)


# ── Tests ────────────────────────────────────────────────────────────────


class TestStructuredOutputPath:
    """Tests for the primary structured-output (.parse) code path."""

    @patch("openai.OpenAI")
    def test_returns_parsed_kg(self, mock_openai_cls):
        kg = KnowledgeGraph(**SAMPLE_KG_DICT)
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.return_value = _make_structured_response(kg)

        result = extract_graph("some text", api_key="test-key")

        assert isinstance(result, KnowledgeGraph)
        assert len(result.nodes) == 3
        assert len(result.edges) == 2
        assert result.nodes[0].id == "BaP"

    @patch("openai.OpenAI")
    def test_passes_model_and_api_key(self, mock_openai_cls):
        kg = KnowledgeGraph(**SAMPLE_KG_DICT)
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.return_value = _make_structured_response(kg)

        extract_graph("text", model="gpt-4o-mini", api_key="sk-test")

        mock_openai_cls.assert_called_once_with(api_key="sk-test", base_url=None)
        call_kwargs = mock_client.beta.chat.completions.parse.call_args
        assert call_kwargs.kwargs["model"] == "gpt-4o-mini"

    @patch("openai.OpenAI")
    def test_system_prompt_sent(self, mock_openai_cls):
        kg = KnowledgeGraph(**SAMPLE_KG_DICT)
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.return_value = _make_structured_response(kg)

        extract_graph("describe BaP", api_key="sk-test")

        call_kwargs = mock_client.beta.chat.completions.parse.call_args
        messages = call_kwargs.kwargs["messages"]
        assert messages[0]["role"] == "system"
        assert messages[0]["content"] == SYSTEM_PROMPT
        assert messages[1]["role"] == "user"
        assert messages[1]["content"] == "describe BaP"


class TestJsonFallbackPath:
    """Tests for the JSON-mode fallback when structured output fails."""

    @patch("openai.OpenAI")
    def test_falls_back_on_attribute_error(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.side_effect = AttributeError("no .parse")
        mock_client.chat.completions.create.return_value = _make_json_response(SAMPLE_KG_DICT)

        result = extract_graph("text", api_key="sk-test")

        assert isinstance(result, KnowledgeGraph)
        assert len(result.nodes) == 3
        mock_client.chat.completions.create.assert_called_once()

    @patch("openai.OpenAI")
    def test_falls_back_on_value_error(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.side_effect = ValueError("parse failed")
        mock_client.chat.completions.create.return_value = _make_json_response(SAMPLE_KG_DICT)

        result = extract_graph("text", api_key="sk-test")

        assert isinstance(result, KnowledgeGraph)
        assert len(result.nodes) == 3

    @patch("openai.OpenAI")
    def test_falls_back_when_parsed_is_none(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.return_value = _make_structured_response(None)
        mock_client.chat.completions.create.return_value = _make_json_response(SAMPLE_KG_DICT)

        result = extract_graph("text", api_key="sk-test")

        assert isinstance(result, KnowledgeGraph)
        assert len(result.nodes) == 3

    @patch("openai.OpenAI")
    def test_fallback_empty_content_returns_empty_graph(self, mock_openai_cls):
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.side_effect = AttributeError("no .parse")
        empty_response = SimpleNamespace(
            choices=[SimpleNamespace(message=SimpleNamespace(content=None))],
            usage=None,
        )
        mock_client.chat.completions.create.return_value = empty_response

        result = extract_graph("text", api_key="sk-test")

        assert isinstance(result, KnowledgeGraph)
        assert len(result.nodes) == 0
        assert len(result.edges) == 0


class TestErrorPropagation:
    """Verify that non-fallback errors propagate instead of being swallowed."""

    @patch("openai.OpenAI")
    def test_auth_error_propagates(self, mock_openai_cls):
        import openai

        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.side_effect = openai.AuthenticationError(
            message="invalid key",
            response=MagicMock(status_code=401),
            body=None,
        )

        with pytest.raises(openai.AuthenticationError):
            extract_graph("text", api_key="bad-key")

    @patch("openai.OpenAI")
    def test_rate_limit_error_propagates(self, mock_openai_cls):
        import openai

        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.side_effect = openai.RateLimitError(
            message="rate limited",
            response=MagicMock(status_code=429),
            body=None,
        )

        with pytest.raises(openai.RateLimitError):
            extract_graph("text", api_key="sk-test")

    @patch("openai.OpenAI")
    def test_timeout_error_propagates(self, mock_openai_cls):
        import openai

        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.side_effect = openai.APITimeoutError(
            request=MagicMock(),
        )

        with pytest.raises(openai.APITimeoutError):
            extract_graph("text", api_key="sk-test")


class TestModelValidation:
    """Verify that the KnowledgeGraph model validator catches bad LLM output."""

    @patch("openai.OpenAI")
    def test_bad_node_references_in_fallback_raises(self, mock_openai_cls):
        bad_data = {
            "nodes": [{"id": "A", "label": "A", "type": "Enzyme"}],
            "edges": [{"source": "A", "target": "MISSING", "type": "ACTIVATES"}],
        }
        mock_client = MagicMock()
        mock_openai_cls.return_value = mock_client
        mock_client.beta.chat.completions.parse.side_effect = AttributeError("no .parse")
        mock_client.chat.completions.create.return_value = _make_json_response(bad_data)

        with pytest.raises(Exception, match="MISSING"):
            extract_graph("text", api_key="sk-test")


class TestModeAwareExtraction:
    def test_exploratory_mode_marks_llm_origin_and_grounds_nodes(self):
        payload = {
            "nodes": [
                {"id": "n1", "label": "CYP1A1", "type": "Enzyme"},
                {"id": "n2", "label": "BaP", "type": "Carcinogen"},
            ],
            "edges": [
                {"source": "n1", "target": "n2", "type": "ACTIVATES"},
            ],
        }

        result = extract_graph(
            "text",
            backend=_StubBackend(payload),
            mode=GraphMode.EXPLORATORY,
        )

        assert all(node.origin == RecordOrigin.LLM for node in result.nodes)
        assert all(edge.origin == RecordOrigin.LLM for edge in result.edges)
        assert result.nodes[0].match_status == MatchStatus.CANONICAL
        assert result.nodes[1].match_status == MatchStatus.ALIAS
        assert result.edges[0].match_status == MatchStatus.CANONICAL

    def test_strict_mode_filters_unmatched_llm_content(self):
        payload = {
            "nodes": [
                {"id": "n1", "label": "CYP1A1", "type": "Enzyme"},
                {"id": "n2", "label": "Unknown Chemical", "type": "Carcinogen"},
            ],
            "edges": [
                {"source": "n1", "target": "n2", "type": "ACTIVATES"},
            ],
        }

        result = extract_graph(
            "text",
            backend=_StubBackend(payload),
            mode=GraphMode.STRICT,
        )

        assert [node.id for node in result.nodes] == ["n1"]
        assert result.edges == []
