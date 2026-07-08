Changelog
=========

Unreleased
----------

- Restored the documented ``ExposoGraph._biomarker_scaffold`` package,
  including split YAML registry sources, mapping rebuild/compare helpers, and
  the scaffold checker used by the Makefile and console scripts
- Promoted all local quality checks to a single required gate: tests, docs,
  biomarker mapping validation, coverage, Ruff, and strict mypy now pass via
  ``make ci``
- Added a root ``Makefile`` and an in-tree GitHub Actions workflow
  (``.github/workflows/ci.yml``) so tests, docs, lint, typing, and coverage
  checks are codified in the repository rather than only described in prose
- Added a new ``production-readiness`` documentation page describing the staged
  CI model, current gate status, and the remaining blockers before a
  production-ready release
- Reduced the Ruff backlog from mixed import/unused/lambda/style failures down
  to line-length-only findings by fixing the non-``E501`` issues directly and
  scoping figure/UI modules whose authored geometry literals are clearer left
  unwrapped
- Updated deployment and README guidance to reflect the current local workflow:
  all configured quality gates now run together through ``make ci``
- Clarified the standalone D3 viewer documentation so ``map/index.html`` is
  described as the bundled reference ``graph-data.js`` payload rather than the
  full quantitative interaction engine
- Documented the new interaction provenance helpers:
  ``get_parameter_provenance()``, ``get_interaction_source_catalog()``, and
  ``get_interaction_expansion_backlog()``
- Added public API reference coverage for the interaction engine, including
  its source catalog, expansion backlog, and synergy-analysis helpers

0.0.4 (2026-04-15)
------------------

Current release.

- Preserved parallel edges in the NetworkX ``MultiDiGraph`` engine so distinct
  evidence edges with the same source, predicate, and target are no longer overwritten
- Corrected KEGG fixed-width record parsing for multi-line ``GENE`` and
  ``PATHWAY`` sections used by the seeding workflow
- Tightened ``metabolism_chain()`` traversal so carcinogen-specific chains do
  not absorb unrelated unlabeled branches through shared enzymes
- ``filter_knowledge_graph()`` now returns detached model copies rather than
  aliasing the source graph objects
- Restored a clean strict ``mypy`` pass for the shipped source tree

0.0.3 (2026-03-21)
------------------

Release metadata synchronized for the current PyPI/GitHub publication.

- Updated package and app version identifiers to ``0.0.3``

0.0.2 (2026-03-19)
------------------

Current development release.

- Added graph ingestion modes: ``exploratory`` and ``strict``
- Added canonical grounding metadata, record origin tracking, and custom predicates
- Added validated/exploratory graph visibility filtering for preview, data, export, and persistence
- Added visibility-aware JSON, HTML, JS, and GEXF export helpers
- Added revision visibility tracking and SQLite schema migration support
- Added clean repository shutdown via context-manager support and fixed prior SQLite resource warnings
- Added external extraction backend abstraction and mode-aware seeded graph preparation

0.0.1 (2026-03-17)
------------------

Initial release.

- Pydantic v2 data models for 7 node types and 10 edge types
- NetworkX MultiDiGraph engine with load/merge/validate
- Literature extraction workflow with structured-output validation
- Streamlit app with manual entry, text extraction, and gene panel loading
- D3.js force-directed graph viewer with dark theme
- Export to JSON, ``graph-data.js`` (D3 viewer), and GEXF (Gephi)
- Curated Tier 1 (13 genes) and Tier 2 (23 genes) reference panels
- Referenced activity-score tables for 18 genes, including evidence metadata
- Test coverage across models, engine, exporter, storage, and reference data
- CI/CD: ruff linting, pytest matrix (3.10–3.12), PyPI publish workflow
