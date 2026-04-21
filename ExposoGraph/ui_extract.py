"""Tab 1: LLM-Assisted Extraction."""

from __future__ import annotations

import streamlit as st

from ._app_shared import get_pending_extraction, get_secret, load_example_text
from .config import GraphMode, LLMProvider
from .engine import GraphEngine
from .llm_backend import OllamaBackend, OpenAIBackend
from .llm_extractor import EXAMPLE_INPUT, extract_graph_with_usage

_OPENAI_MODELS = ["gpt-4o", "gpt-4o-mini", "gpt-4.1", "gpt-4.1-mini", "gpt-4.1-nano"]
_OLLAMA_MODELS = ["llama3.1", "llama3.2", "mistral", "gemma2", "deepseek-r1", "phi4"]


def render(engine: GraphEngine) -> None:
    """Render the LLM Extract tab."""
    st.markdown("#### Extract Knowledge Graph from Text")
    st.caption(
        "Describe a carcinogen metabolism pathway in plain English. "
        "The LLM will extract nodes and edges automatically."
    )
    pending_warnings = st.session_state.pop("extract_merge_warnings", [])
    if pending_warnings:
        st.warning("\n".join(pending_warnings))

    col_mode1, col_mode2 = st.columns(2)
    with col_mode1:
        provider = st.selectbox(
            "LLM Provider",
            [p.value for p in LLMProvider],
            index=0,
            help="Choose OpenAI (cloud) or Ollama (local).",
        )
    with col_mode2:
        graph_mode = st.selectbox(
            "Graph Mode",
            [mode.value for mode in GraphMode],
            index=0,
            help="Exploratory keeps provisional entities. Strict keeps only canonically grounded content.",
        )

    if provider == LLMProvider.OPENAI:
        col_cfg1, col_cfg2 = st.columns(2)
        with col_cfg1:
            api_key = st.text_input(
                "OpenAI API Key",
                type="password",
                help="Leave blank to use the OPENAI_API_KEY env var.",
            )
        with col_cfg2:
            model = st.selectbox("Model", _OPENAI_MODELS, index=0)
    else:
        col_cfg1, col_cfg2 = st.columns(2)
        with col_cfg1:
            ollama_url = st.text_input(
                "Ollama URL",
                value="http://localhost:11434",
                help="Base URL of your running Ollama instance.",
            )
        with col_cfg2:
            model = st.selectbox("Model", _OLLAMA_MODELS, index=0)
        api_key = ""

    user_text = st.text_area(
        "Pathway description",
        key="extract_text",
        height=220,
        placeholder=EXAMPLE_INPUT,
    )

    col1, col2 = st.columns([1, 3])
    with col1:
        run_extract = st.button("Extract", type="primary", use_container_width=True)
    with col2:
        st.button(
            "Load example text",
            use_container_width=True,
            on_click=load_example_text,
        )

    if run_extract and user_text.strip():
        # Build the appropriate backend
        backend: OpenAIBackend | OllamaBackend
        if provider == LLMProvider.OPENAI:
            resolved_key = api_key or get_secret("OPENAI_API_KEY")
            if not resolved_key:
                st.error("Provide an OpenAI API key or set OPENAI_API_KEY before extracting.")
                st.stop()
            backend = OpenAIBackend(api_key=resolved_key or None)
        else:
            backend = OllamaBackend(base_url=ollama_url)

        with st.spinner("Calling LLM…"):
            try:
                kg, usage = extract_graph_with_usage(
                    user_text, model=model, backend=backend, mode=graph_mode,
                )
                st.success(
                    f"Extracted **{len(kg.nodes)}** nodes and "
                    f"**{len(kg.edges)}** edges "
                    f"({usage.total_tokens} tokens, {usage.duration_ms}ms)"
                )
                st.session_state.pending_extraction = kg.model_dump(mode="json")
                st.session_state.pending_extraction_mode = graph_mode
                st.rerun()

            except Exception as exc:
                st.error(f"Extraction failed: {exc}")

    pending_kg = get_pending_extraction()
    if pending_kg is not None:
        pending_mode = st.session_state.get("pending_extraction_mode", GraphMode.EXPLORATORY.value)
        with st.expander("Preview extracted data", expanded=True):
            st.caption(f"Prepared in `{pending_mode}` mode")
            st.json(pending_kg.model_dump(mode="json"))

        col_merge, col_discard = st.columns(2)
        with col_merge:
            if st.button("Merge into graph", type="primary", use_container_width=True):
                try:
                    warnings = engine.merge(pending_kg, mode=pending_mode)
                    if warnings:
                        st.session_state.extract_merge_warnings = warnings
                    st.session_state.pop("pending_extraction", None)
                    st.session_state.pop("pending_extraction_mode", None)
                    st.rerun()
                except Exception as exc:
                    st.error(f"Merge failed: {exc}")
        with col_discard:
            if st.button("Discard extraction", use_container_width=True):
                st.session_state.pop("pending_extraction", None)
                st.session_state.pop("pending_extraction_mode", None)
                st.rerun()
