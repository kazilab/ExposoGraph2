``ExposoGraph.llm_extractor``
=============================

LLM-powered entity and relation extraction using pluggable backends.

``extract_graph(...)`` and ``extract_graph_with_usage(...)`` accept a
``mode`` argument:

- ``exploratory`` keeps provisional nodes and edges after grounding
- ``strict`` returns only canonically grounded content

.. autofunction:: ExposoGraph.llm_extractor.extract_graph

.. data:: ExposoGraph.llm_extractor.SYSTEM_PROMPT

   The system prompt sent to the LLM, containing the full schema specification
   and extraction guidelines.

.. data:: ExposoGraph.llm_extractor.EXAMPLE_INPUT

   Example input text describing the Benzo[a]pyrene metabolism pathway,
   useful for testing and demonstration.
