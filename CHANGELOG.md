# Changelog

All notable changes to ExposoGraph will be documented in this file.

## [Unreleased]

### Added

- **Production-readiness scaffolding** — added a root `Makefile`, an in-tree
  GitHub Actions workflow at `.github/workflows/ci.yml`, and a new
  `docs/production-readiness.rst` page so the repository now encodes its
  required vs. advisory quality gates instead of relying on undocumented local
  commands.
- **Ruff backlog triage** — fixed the non-`E501` lint failures directly
  (`__version__` export, lambda cleanup, unused locals/imports, test import
  ordering, example-script cleanup) and scoped figure/UI/test modules whose
  long authored geometry or layout literals are intentionally clearer without
  wrapping. The remaining Ruff backlog is now line-length-only.
- **Connected seven previously orphan carcinogens** in the showcase graph —
  added procarcinogen → reactive-intermediate `ACTIVATES` edges for Benzidine,
  NDMA, Vinyl Chloride, and 17β-Estradiol, plus full metabolism chains for
  DMBA (CYP1A1/1B1 → 3,4-epoxide → EPHX1 → diol → diol-epoxide → DNA adduct)
  and MeIQx (CYP1A2 → N-OH-MeIQx → NAT1 → acetoxy ester → C8-dG). Updated
  `reference_graphs/full_legends_graph.js` and regenerated
  `ExposoGraph/map/graph-data.js` and `ExposoGraph/temp_kg/graph-data.js`.
- **Added five auxiliary Phase I CYP nodes to the bundled reference graph** —
  `CYP2B6`, `CYP2C9`, `CYP2C19`, `CYP2D6`, and `CYP2F1` now appear in the
  bundled viewer as lightweight CYP450 pathway members, keeping the map aligned
  with expanded interaction-model coverage without introducing the full kinetic
  matrix as extra viewer edges. These nodes are present in the current bundled
  reference viewer graph (**212 nodes / 313 edges**).
- **Added an interaction-source catalog and explicit expansion backlog** —
  `interaction_parameters.json` and `parameter_provenance.json` now record
  BRENDA, Rendic & Guengerich 2012, PharmGKB/ClinPGx, IARC Vol. 100A-F,
  ATSDR, and PubMed as prioritized curation sources, and expose a structured
  `not yet in JSON` backlog showing that 27 requested pairs are already
  parameterized and 43 remain pending. The backlog now also includes a
  green/yellow/red scientific-validity triage for the requested Tier A / Tier B
  Ki-curation list, separating direct competition-pair candidates from
  provisional and model-mismatched entries. `interaction_engine.py` now
  exposes `get_interaction_source_catalog()` and
  `get_interaction_expansion_backlog()`.
- **Standalone D3.js viewer overhaul** (`ExposoGraph/map/index.html`):
  - Shape-coded node-type legend glyphs matching the on-canvas shapes
    (diamond / hexagon / rounded rect / circle) with the same fill-opacity
    ladder used by the simulation.
  - Responsive sidebar width with `clamp(280px, 22vw, 420px)` plus a
    draggable 5 px splitter (min 240 px, max 640 px; double-click to reset).
  - Live, selection-driven header counts for *nodes*, *edges*, *node types*,
    and *edge types* — wired into search, type highlight, carcinogen-class
    filter, tissue filter, and reset paths via a new `updateHeaderStats()`
    helper. Force simulation re-centers on container resize.
  - Wrappable panel headers and legend labels so long strings
    (e.g. "Carcinogens (IARC-classified)", "Filter by Carcinogen Class") no
    longer clip at narrow sidebar widths.
- New procarcinogen → intermediate entries in `_EDGE_SPECS`
  (`example_graphs.py`) keep the Python builders in sync with the
  hand-curated `reference_graphs/full_legends_graph.js` source of truth.

### Changed

- `README.md` and `docs/deployment.rst` now reflect the staged release model:
  regression tests and docs are blocking CI jobs, while coverage, Ruff, and
  strict mypy run as advisory jobs until those backlogs are reduced.
- `build_full_legends_graph()` now returns **107 nodes / 124 edges**
  (previously 96 / 102), and `build_full_legends_graph(include_heavy_metals=True)`
  currently resolves to the same **212 / 313** bundled footprint as
  `build_reference_graph()`. `build_full_legends_architecture_summary()`
  remains 107 / 124 base-showcase summary.
- Single source of truth for `IARCGroup` — the duplicate enum in
  `expanded_metals.py` now imports from `db_clients.iarc`, and `__all__` in
  `ExposoGraph/__init__.py` was deduplicated.
- `pyproject.toml` — added a `per-file-ignores` block for `E501` on the
  long-string data-literal modules (`reference_data.py`, `expanded_metals.py`,
  `wave2_classes.py`, `oxidative_stress.py`, `population_genomics.py`,
  `tissue_subgraphs.py`, `exposure_engine.py`) so line-length warnings no
  longer swamp genuine lint findings.

### Fixed

- Removed the dead "Risk Modules (v2.0)" sidebar panel (CSS + HTML) from
  `map/index.html`; the live counters and filter panels are now the canonical
  sidebar UI.
- Updated 453 hard-coded node / edge / type-count assertions across
  `tests/test_example_graphs.py`, `tests/test_reference_graph.py`, and
  `tests/test_figure_architecture.py` to match the new graph sizes. All tests
  green.

## [0.0.5] - 2026-04-20

### Added

- **Quantitative flux engine** (`flux_engine.py`) — Michaelis–Menten and Hill-equation
  activation / detoxification flux modeling across 7 carcinogen classes (PAH, Aflatoxin,
  Aldehyde, Nitrosamine, NDMA, HCA, Benzene) with genotype activity scores and curated
  or GTEx v8 tissue weights, plus a `flux_cli` wrapper and `sensitivity_analysis`
- **Exposure integration engine** (`exposure_engine.py`) — NHANES/EPA IRIS/OSHA-derived
  exposure profiles, lifetime excess cancer risk (LECR) scoring with de-minimis and
  action-threshold flags, and an `exposure_cli` wrapper
- **Multi-carcinogen interaction engine** (`interaction_engine.py`) — co-exposure
  synergy matrix including GSH depletion, enzyme induction, and competitive inhibition
  effects; `interaction_cli` wrapper; critical-interaction detection
- **Tissue-specific subgraphs** (`tissue_subgraphs.py`) — Liver, Lung, Kidney, Breast,
  Prostate, Bladder, Colon, Brain views with quantitative GTEx weights and
  cross-carcinogen tissue metabolism profiles
- **Population genomics** (`population_genomics.py`) — 1000 Genomes and gnomAD v4.1
  allele-frequency lookups across EUR/EAS/SAS/AFR/AMR/OCE populations,
  population-aware panel builders, and All of Us BigQuery SQL templates
- **Oxidative stress module** (`oxidative_stress.py`) — ROS marker and antioxidant-
  response panels for heavy-metal and persistent-organic exposures, plus integration
  helpers for the core metabolism graph
- **TK/TD modeling** (`tk_td_modeling.py`) — 1- and 2-compartment pharmacokinetic
  simulations, steady-state and cumulative-dose derivations, and `scipy.integrate`-
  based ODE solving
- **Cross-species scaling** (`cross_species.py`) — allometric scaling (default
  exponent 0.75), species physiology comparisons, and study-duration extrapolation
- **Expanded heavy-metal catalogue** (`expanded_metals.py`) — metal-specific profiles,
  IARC-group filtering, cancer-site associations, and oxidative-stress integration
- **Wave 2 carcinogen classes** (`wave2_classes.py`) — Aldehydes, Dioxins/AhR ligands,
  Dietary N-Nitroso compounds, and Chlorinated solvents with per-class enzyme panels
- **Unified high-level API** (`unified_api.py`) — `PatientRiskProfile` and
  `patient_risk_query()` combining flux, exposure, interaction, and tissue reports
- **Synthetic population simulation** (`population_simulation/`) — ancestry-stratified
  cohort generator, validation framework, phenotype extraction for All of Us, and
  batch simulation runner
- **14 IARC carcinogen class coverage** and an expanded reference graph shipping
  **162 nodes / 225 edges** across all bundled showcase artifacts and viewer bundles
- **Smoke test suite** for all 13 new modules (flux, exposure, interaction, tissue,
  population genomics, oxidative stress, TK/TD, cross-species, expanded metals, wave2
  classes, unified API, population simulation, CLI entrypoints)

### Changed

- `numpy` (≥1.26) and `scipy` (≥1.11) are now required runtime dependencies; they
  were previously optional or absent
- `paper_architecture_overrides()` now reflects the expanded reference graph
  (128 nodes / 161 edges with the heavy-metal overlay) used by the architecture
  infographic
- Reference gene panels retuned to 23 Tier 2 genes and 36 total panel entries,
  with the DNA Repair sub-panel narrowed from 7 to 5 enzymes

### Removed

- `render_genotype_comparison_figure` has been retired from
  `exemplar_pathways_figure.py`; its callers are now covered by the unified-profile
  comparison helpers in `unified_api.py`
- The `viewer-genotype-feedback` Dash component has been removed from the advanced
  viewer layout in favour of the new patient-profile workflow

## [0.0.4] - 2026-04-15

### Fixed

- **Parallel edge preservation** — `GraphEngine` now preserves multiple edges
  that share the same `(source, type, target)` triple instead of silently
  overwriting earlier evidence records
- **KEGG fixed-width parsing** — the KEGG client now correctly parses
  multi-line `GENE` sections with numeric gene IDs and multi-line `PATHWAY`
  sections from `get/{id}` records, so seeded graphs retain the expected gene
  symbols and pathway memberships
- **`metabolism_chain()` scope leakage** — carcinogen-specific pathway
  traversal no longer absorbs unrelated unlabeled branches merely because they
  share an upstream enzyme node
- **Filtered graph aliasing** — `filter_knowledge_graph()` now returns
  detached model copies instead of reusing the original `Node` and `Edge`
  objects

### Changed

- Strict `mypy` checks are back in sync with the shipped source tree
- Documentation now reflects the current `0.0.4` package version and the
  clarified semantics for graph filtering, KEGG seeding, and metabolism-chain
  traversal

## [0.0.3] - 2026-03-21

### Changed

- Version bump to `0.0.3` for the current release candidate and synchronized
  package, app, and documentation metadata ahead of PyPI/GitHub publication

## [0.0.2] - 2026-03-19

### Added

- **Graph analysis module** (`graph_analysis.py`) — shortest path, all shortest paths,
  degree/betweenness centrality, metabolism chain traversal, pathway subgraph extraction,
  and variant impact scoring with activity-score integration
- **Multi-LLM support** — pluggable backend architecture (`llm_backend.py`) with OpenAI
  (structured output + JSON-mode fallback, exponential backoff retry) and Ollama
  (`/api/chat`) backends; `LLMProvider` enum in `config.py`; token/cost tracking via
  `UsageRecord` dataclass
- **Public database clients** (`db_clients/`) — KEGG REST API client (pathway + gene
  lookups), CTD batch query client (chemical-gene interactions with organism filtering),
  bundled IARC carcinogen classification data with monograph volume references
- **Seeder module** (`seeder.py`) — converts KEGG pathways and CTD interactions directly
  into `KnowledgeGraph` objects with provenance tracking and heuristic edge-type inference
- **Streamlit UI updates** — LLM provider selector (OpenAI/Ollama) in extraction tab with
  Ollama-specific URL and model inputs; token usage display on successful extraction
- **Pre-built example** — `examples/bap_graph.json` (20 nodes, 20 edges covering full
  BaP activation, detoxification, adduct formation, and repair pathways) and
  `examples/build_bap_graph.py` demo script
- **GEXF export** (`exporter.py`) — Gephi-compatible graph export with automatic
  JSON-serialization of complex node/edge attributes
- **Comprehensive test suite** — 226 tests across 14 test modules; `test_integration.py`
  for end-to-end pipeline validation; `test_llm_backend.py`, `test_db_clients.py`,
  `test_seeder.py`, `test_graph_analysis.py` for new modules
- **CI/CD** — GitHub Actions workflow with pytest-cov (85% coverage gate, currently 96%),
  ruff linting, mypy strict type checking
- **Sphinx documentation** scaffolding with ReadTheDocs configuration

### Fixed

- **Fabricated PMID** — replaced non-existent PMID 41024270 (OGG1) with verified
  PMID 25588927 (Zhou et al. meta-analysis)
- **PMID-title mismatches** — corrected titles for PMID 29194389 (CYP2A6, Tanner &
  Tyndale 2017) and PMID 23665933 (CYP3A4, Okubo et al. 2013)
- **ClinPGx URL pattern** — `clinpgx.org/gene/{symbol}` does not resolve; replaced with
  PharmGKB accession ID lookup (`clinpgx.org/gene/{accession_id}`) for all 28 gene panels
- **IARC references** — added specific IARC Monograph volume numbers and years to all
  30 classification entries (previously generic "Group N" only)
- **Missing PubMed references** — added literature citations to 5 activity score metadata
  entries (CYP2D6, CYP2C9, CYP2C19, UGT1A1, XPC) that previously had none
- **GEXF export crash** — NetworkX GEXF writer failed on list/dict node attributes
  (provenance, curation); fixed by JSON-serializing non-scalar values before write

### Changed

- `extract_graph()` now delegates to pluggable `LLMBackend` protocol; new
  `extract_graph_with_usage()` returns both the graph and a `UsageRecord`
- `__init__.py` public API expanded with all new module exports
- `pyproject.toml` — added `[ollama]`, `[db]`, and `[docs]` optional dependency groups;
  pytest-cov configuration; mypy overrides for new dependencies

## [0.0.1] - 2026-03-17

### Added

- Initial release: Pydantic v2 models, NetworkX graph engine, LLM extraction (OpenAI),
  JSON/D3.js/HTML export, Streamlit UI, reference gene panels and activity scores
