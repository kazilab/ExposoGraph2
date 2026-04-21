# ExposoGraph 2.0
<!-- PyPI version badge -->
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://ExposoGraph.streamlit.app)
[![PyPI version](https://img.shields.io/pypi/v/ExposoGraph.svg)](https://pypi.org/project/ExposoGraph/)
[![Documentation Status](https://readthedocs.org/projects/ExposoGraph/badge/?version=latest)](https://ExposoGraph.readthedocs.io/en/latest/?badge=latest)
<!-- [![ResearchSquare](https://img.shields.io/badge/ResearchSquare-rs--9202489%2Fv1-00A0E0.svg)](https://www.researchsquare.com/article/rs-9202489/v1) -->
<!-- [![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.03.22.713456-b31b1b.svg)](https://doi.org/10.64898/2026.03.22.713456) -->
<!-- PyPI version badge -->
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-kazilab%2FExposoGraph-181717?logo=github&logoColor=white)](https://github.com/kazilab/ExposoGraph)
[![@KaziLab.se](https://img.shields.io/website?url=https://www.kazilab.se/)](https://www.kazilab.se/)
<!-- PyPI version badge -->

Build, curate, and export carcinogen metabolism knowledge graphs using LLM-powered extraction and manual entry, then run quantitative multi-carcinogen, tissue-aware risk analysis on the resulting graph.

**ExposoGraph 2.0** expands coverage to **14 IARC carcinogen classes**, ships a bundled reference knowledge graph of **212 nodes / 313 edges**, adds seven quantitative risk modules: enzyme flux modeling, exposure integration, multi-carcinogen interactions, tissue-specific subgraphs, population genomics, oxidative stress integration, toxicokinetic/toxicodynamic modeling, cross-species scaling, and a synthetic population simulator.

Developed by: **Data analysis team @ KaziLab**
Contact: **exposograph@kazilab.se**
Copyright: **KaziLab**

## Acknowledgement

Parts of this documentation and code were created with assistance from ChatGPT Codex and Claude Code.

## Features

### Graph authoring and analysis
- **Multi-LLM Extraction** — Describe a carcinogen metabolism pathway in plain English; OpenAI (GPT-4o) or local Ollama models extract structured nodes and edges automatically
- **Mode-Aware Ingestion** — Use `exploratory` mode to keep provisional entities or `strict` mode to keep only canonically grounded content
- **Manual Entry** — Add and annotate nodes and edges with full provenance and curation tracking
- **Graph Analysis** — Shortest path, centrality, metabolism chain traversal, pathway subgraph, variant impact scoring
- **Public DB Integration** — KEGG pathway lookups, CTD chemical-gene interactions, IARC carcinogen classifications
- **Interactive Preview** — Color-coded Streamlit AGraph visualization with hover metadata, search/filter controls, zoom, and downloadable Plotly HTML exports
- **Advanced Graph Viewer** — Dash Cytoscape viewer with sidebar search, legends, carcinogen filters, detail panel, image export, and saved layout JSON
- **Standalone D3.js Viewer** — Zero-dependency `map/index.html` for the bundled reference `graph-data.js` payload, with shape-coded nodes, a draggable sidebar resizer, a runtime bundle subtitle, and live header counts that react to search, node-type, carcinogen-class, and tissue filters
- **Validated vs Exploratory Views** — Filter the current graph to `all`, `validated only`, or `exploratory only` in the UI and export pipeline
- **Rich Annotations** — Structured provenance records, source manifests, curated KEGG pathway coverage, PubMed IDs, tissue context, pharmacogenomic variants, activity scores
- **Multiple Export Formats** — Standalone Plotly HTML, parseable app HTML, JSON, D3.js viewer (`graph-data.js`), GEXF (Gephi)
- **Viewer Data Contract** — Export a Cytoscape-ready JSON bundle and saved preset layout for richer web-style exploration without maintaining custom JavaScript
- **Validation** — Referential integrity checks at model level, dangling edge detection, carcinogen context validation
- **Persistent Storage** — SQLite-backed revision history with explicit export visibility tracking and atomic operations

### Quantitative risk modules (2.0)
- **Flux Engine** — Michaelis–Menten / Hill-equation activation vs. detoxification flux ratios per carcinogen class, weighted by genotype activity scores and GTEx v8 tissue expression
- **Exposure Integration** — NHANES / EPA IRIS / OSHA-derived exposure profiles produce lifetime excess cancer risk (LECR) scores with de-minimis and action-threshold flags
- **Multi-Carcinogen Interactions** — Synergy matrix across co-exposures with GSH depletion, enzyme induction, and competitive inhibition effects
- **Tissue-Specific Subgraphs** — Liver, Lung, Kidney, Breast, Prostate, Bladder, Colon, Brain views with quantitative GTEx expression weights
- **Population Genomics** — 1000 Genomes and gnomAD v4.1 allele-frequency lookups across EUR/EAS/SAS/AFR/AMR/OCE populations
- **Oxidative Stress** — ROS marker and antioxidant-response annotations for heavy-metal and persistent-organic exposures
- **TK/TD Modeling** — 1- and 2-compartment pharmacokinetic simulations with allometric cross-species scaling (rat → human)
- **Population Simulation** — Synthetic cohort generator with ancestry-stratified genotype sampling and lifetime risk distributions
- **Unified API** — `patient_risk_query()` combines flux, exposure, interaction, and tissue reports in a single call

## Quick Start

### Try Without an API Key

A pre-built Benzo[a]pyrene metabolism graph is included:

```bash
pip install -e .
python examples/build_bap_graph.py
```

This loads `examples/bap_graph.json` (20 nodes, 20 edges covering the full BaP → BPDE → DNA adduct pathway), runs graph analysis, and exports to HTML and JSON.

For a no-API-key demonstration of strict vs exploratory handling and
validated-only exports:

```bash
python examples/mode_visibility_demo.py
```

<details>
<summary>Sample output</summary>

```
Loaded graph: 20 nodes, 20 edges

Shortest path CYP1A1 → BPDE-dG: CYP1A1 → BPDE → BPDE_dG

Top-5 nodes by degree centrality:
  CYP1A1        0.263
  BPDE          0.211
  CYP1B1        0.158
  BPDE_dG       0.158
  BPDE_GSH      0.158

BaP metabolism chain: 16 nodes, 13 edges
  Activation edges:    4
  Detoxification edges: 3
  Adduct edges:        1
  Repair edges:        3

Variant impact score for CYP1A1:
  Activity score:        1.0
  Downstream adducts:    1
  Impact score:          1.00
```

</details>

### Sample JSON

```json
{
  "nodes": [
    {"id": "BaP", "label": "Benzo[a]pyrene", "type": "Carcinogen", "group": "PAH", "iarc": "Group 1"},
    {"id": "CYP1A1", "label": "CYP1A1", "type": "Enzyme", "phase": "I", "role": "Activation"},
    {"id": "BPDE", "label": "BPDE", "type": "Metabolite", "reactivity": "High"},
    {"id": "BPDE_dG", "label": "BPDE-N2-dG", "type": "DNA_Adduct"}
  ],
  "edges": [
    {"source": "CYP1A1", "target": "BPDE", "type": "ACTIVATES", "carcinogen": "BaP"},
    {"source": "BPDE", "target": "BPDE_dG", "type": "FORMS_ADDUCT", "carcinogen": "BaP"}
  ]
}
```

### Streamlit App

```bash
pip install -e ".[streamlit]"
streamlit run ExposoGraph/app.py
```

App mode defaults to `stateless`, which disables server-side saves and
is appropriate for public web deployment. To enable local revision history
and file saves on your own machine:

```bash
export ExposoGraph_MODE=local
streamlit run ExposoGraph/app.py
```

### Jupyter

```bash
pip install -e ".[notebook]"
jupyter lab
```

No notebook file is currently bundled in this repository. Use the installed
package from your own notebook, or start from the runnable examples in
`examples/`.

### Standalone D3.js Viewer (`map/index.html`)

A zero-install HTML viewer is bundled at `ExposoGraph/map/index.html` with its
graph payload in `ExposoGraph/map/graph-data.js` (**212 nodes / 313 edges** by
default — the current bundled reference graph). This shipped payload is a
curated graph export, not a rendering of the quantitative interaction engine.
Open the HTML file
directly in a browser or serve the `map/` folder statically:

```bash
python -m http.server --directory ExposoGraph/map 8000
# Open http://localhost:8000/
```

Viewer features:

- **Shape-coded node glyphs** — the sidebar *Node Types* legend renders each
  entry with its on-canvas shape (diamond = Carcinogen, hexagon = DNA adduct,
  rounded rectangle = Pathway / Tissue, circle = Enzyme / Gene / Metabolite).
  The shipped reference bundle currently contains `Carcinogen`, `Enzyme`,
  `Metabolite`, `DNA_Adduct`, and `Pathway` nodes.
- **Resizable sidebar** — the width auto-scales with the viewport via
  `clamp(280px, 22vw, 420px)` and a vertical splitter between the sidebar and
  the graph can be dragged to any width between 240 px and 640 px. Double-click
  the splitter to reset to the responsive default.
- **Runtime bundle subtitle** — the header reports the currently loaded bundle
  size and whether the heavy-metal overlay is present, without crowding the
  clickable sidebar controls.
- **Live header counts** — the four counters in the top-right (*nodes*,
  *edges*, *node types*, *edge types*) update dynamically as you search, click
  a node-type legend entry, pick a carcinogen class, or apply a GTEx tissue
  filter — they always reflect the currently-visible subset.
- **Interactive filters** — search by label / id / group / variant, toggle
  node types, filter by carcinogen class, and slice by tissue + GTEx
  expression threshold. The shipped bundle embeds GTEx weights on selected
  nodes; it does not currently ship standalone `Tissue` nodes.

To regenerate the bundled data from a Python graph:

```python
from ExposoGraph import build_reference_engine, to_graph_data_js

engine = build_reference_engine()
to_graph_data_js(engine, "ExposoGraph/map/graph-data.js")
```

### Advanced Viewer

```bash
pip install -e ".[viewer]"
```

```python
from ExposoGraph import (
    GraphVisibility,
    ViewerLayoutMode,
    launch_dash_viewer,
    write_cytoscape_bundle,
)

write_cytoscape_bundle(
    engine,
    "exports/graph_cytoscape.json",
    visibility=GraphVisibility.ALL,
    layout_mode=ViewerLayoutMode.COSE,
)

launch_dash_viewer(
    engine,
    visibility=GraphVisibility.ALL,
    layout_mode=ViewerLayoutMode.COSE,
    port=8050,
)
```

### Python Library

```bash
pip install -e .
```

```python
from ExposoGraph import (
    GraphEngine,
    GraphMode,
    GraphVisibility,
    centrality,
    extract_graph,
    metabolism_chain,
    to_json,
)

# LLM extraction (requires OpenAI API key)
# exploratory: keep unmatched or custom content
kg = extract_graph(
    "Benzo[a]pyrene is activated by CYP1A1...",
    mode=GraphMode.EXPLORATORY,
)
engine = GraphEngine()
engine.merge(kg, mode=GraphMode.EXPLORATORY)

# Analysis
scores = centrality(engine, method="degree")
chain = metabolism_chain(engine, "BaP")

# Visibility-aware export
to_json(engine, "graph_validated.json", visibility=GraphVisibility.VALIDATED_ONLY)
```

#### Using Ollama (Local LLM)

```python
from ExposoGraph import GraphMode, extract_graph
from ExposoGraph.llm_backend import OllamaBackend

backend = OllamaBackend(base_url="http://localhost:11434")
kg = extract_graph(
    "BaP is activated by CYP1A1...",
    backend=backend,
    model="llama3.1",
    mode=GraphMode.EXPLORATORY,
)
```

#### Public Database Integration

```python
from ExposoGraph import GraphMode
from ExposoGraph.db_clients import IARCClassifier
from ExposoGraph.seeder import seed_from_ctd, seed_from_kegg_pathway

# Seed from KEGG pathway
kg = seed_from_kegg_pathway("hsa05204", mode=GraphMode.STRICT)

# Seed from CTD
kg = seed_from_ctd("Benzo(a)pyrene", mode=GraphMode.EXPLORATORY)

# IARC classification lookup
clf = IARCClassifier()
clf.classify("Benzo[a]pyrene")  # → IARCGroup.GROUP_1
```

#### Behavior Notes

- `GraphEngine` preserves parallel edges, so multiple evidence records can coexist even when they share the same source, edge type, and target.
- `metabolism_chain()` stays scoped to the carcinogen-linked metabolism subgraph and does not pull in unrelated unlabeled branches through a shared enzyme.
- `filter_knowledge_graph()` returns detached copies of nodes and edges, so mutating a filtered result does not mutate the original graph.
- The KEGG client and KEGG-based seeding support fixed-width multi-line KEGG records, including numeric `GENE` rows and continued `PATHWAY` sections.

#### Reference Curation Metadata

```python
from ExposoGraph import CURATION_SOURCE_MANIFEST, REFERENCE_KEGG_PATHWAYS

primary_sources = CURATION_SOURCE_MANIFEST["primary_sources"]
kegg_ids = [entry["pathway_id"] for entry in REFERENCE_KEGG_PATHWAYS]
```

#### Canonical Reference Graph

```python
from ExposoGraph import (
    build_reference_architecture_summary,
    build_reference_engine,
    build_reference_graph,
    write_reference_exports,
)

graph = build_reference_graph()
engine = build_reference_engine()
summary = build_reference_architecture_summary()
artifacts = write_reference_exports("exports/reference")

len(graph.nodes)            # 212
len(graph.edges)            # 313
engine.validate()           # []
summary.node_count          # 212
summary.edge_count          # 313
artifacts["graph_data_js"]  # bundled viewer export path
```

#### Summary

The legacy `build_full_legends_*` showcase API remains available for the
paper-aligned **107-node / 124-edge** base example graph. The current shipped
viewer payload now matches the bundled reference graph instead; if you call
`build_full_legends_graph(include_heavy_metals=True)`, it currently resolves to
the same **212 / 313** bundled footprint as `build_reference_graph()`.

```python
from ExposoGraph import build_full_legends_architecture_summary

summary = build_full_legends_architecture_summary()

summary.node_count          # 107
summary.edge_count          # 124
summary.node_type_counts    # {'Carcinogen': 15, 'Enzyme': 41, 'Metabolite': 33,
                            #  'DNA_Adduct': 12, 'Pathway': 6}
summary.edge_type_counts    # {'ACTIVATES': 42, 'DETOXIFIES': 24, 'TRANSPORTS': 7,
                            #  'FORMS_ADDUCT': 16, 'REPAIRS': 11, 'PATHWAY': 24}
summary.carcinogen_classes  # grouped class inventories for section 2.2 rewrites
```

#### Optional Androgen Module

```python
from ExposoGraph import build_androgen_module_graph, build_full_legends_graph

androgen_only = build_androgen_module_graph()
showcase_with_androgen = build_full_legends_graph(include_androgen_module=True)
```

### Quantitative Risk Quick Start

```python
from ExposoGraph import (
    compute_pathway_flux,
    compute_lifetime_cancer_risk,
    compute_interaction_matrix,
    patient_risk_query,
    generate_synthetic_cohort,
)

# 1. Per-pathway activation/detox flux ratio
flux = compute_pathway_flux("PAH", {"CYP1A1": "NM", "GSTM1": "NM"}, tissue="Lung")
print(flux.net_ratio, flux.risk_classification)

# 2. Lifetime excess cancer risk against a published slope factor
risk = compute_lifetime_cancer_risk(
    "PAH", {"CYP1A1": "NM", "GSTM1": "NM"}, daily_dose_mg_kg=0.001,
)
print(risk.lecr, risk.exceeds_action_threshold)

# 3. Multi-carcinogen interaction matrix with synergy factors
interactions = compute_interaction_matrix({"PAH": 1.0, "HCA": 0.5}, tissue="Liver")
print(interactions.interaction_factor, interactions.summary)

# 4. Unified patient profile combining all of the above
profile = patient_risk_query(
    {"CYP1A1": "NM", "GSTM1": "NM", "NAT2": "NM"},
    tissue="Liver",
)

# 5. Synthetic cohort with ancestry-stratified genotypes
cohort = generate_synthetic_cohort(n=1000, seed=42)
```

Three command-line entry points wrap the same engines for batch use:

```bash
python -m ExposoGraph.flux_cli --help
python -m ExposoGraph.exposure_cli --help
python -m ExposoGraph.interaction_cli --help
```

#### Interaction Provenance and Backlog

```python
from ExposoGraph import (
    assumed_ki_pairs,
    get_interaction_expansion_backlog,
    get_interaction_source_catalog,
    get_parameter_provenance,
)

catalog = get_interaction_source_catalog()
backlog = get_interaction_expansion_backlog()
provenance = get_parameter_provenance()

catalog[0]["source"]  # "BRENDA enzyme database"
catalog[1]["source"]  # Rendic & Guengerich 2012 review
backlog["remaining_pairs_to_parameterize"]  # 43
backlog["scientific_validity_triage"]["red_exclude_or_rename_pairs"]
len(assumed_ki_pairs())  # competitive pairs still using Ki = Km
```

The interaction model keeps numeric `Km_uM`, `Vmax_relative`, and any explicit
`Ki_uM` values in `ExposoGraph/data/interaction_parameters.json`, while
`ExposoGraph/data/parameter_provenance.json` stores source citations,
confidence grades, `ki_status`, and the structured expansion backlog.

## Configuration

Set your OpenAI API key as an environment variable:

```bash
export OPENAI_API_KEY="sk-..."
```

Or enter it in the Streamlit sidebar when running the app.

For Streamlit Cloud deployment, add the key to `.streamlit/secrets.toml`:

```toml
OPENAI_API_KEY = "sk-..."
```

## Graph Modes and Visibility

Two separate controls now shape how data moves through the system:

- **Graph mode** controls ingestion behavior:
  - `exploratory` keeps unmatched or custom entities and marks them as provisional
  - `strict` keeps only canonically grounded nodes and edges
- **Graph visibility** controls viewing and export behavior:
  - `all`
  - `validated_only`
  - `exploratory_only`

Typical pattern:

```python
from ExposoGraph import (
    GraphEngine,
    GraphMode,
    GraphRepository,
    GraphVisibility,
    ViewerLayoutMode,
    launch_dash_viewer,
    write_cytoscape_bundle,
    extract_graph,
    to_interactive_html,
    to_plotly_html,
)

engine = GraphEngine()
kg = extract_graph(
    "BaP induces CYP1A1 and forms BPDE adducts",
    mode=GraphMode.STRICT,
)
engine.merge(kg, mode=GraphMode.STRICT)

to_interactive_html(
    engine,
    "validated_graph.html",
    visibility=GraphVisibility.VALIDATED_ONLY,
)

to_plotly_html(
    engine,
    "validated_graph_plotly.html",
    visibility=GraphVisibility.VALIDATED_ONLY,
)

write_cytoscape_bundle(
    engine,
    "validated_graph_cytoscape.json",
    visibility=GraphVisibility.VALIDATED_ONLY,
    layout_mode=ViewerLayoutMode.PRESET,
)

with GraphRepository("data/ExposoGraph.sqlite3") as repo:
    repo.save_engine(
        graph_key="bap_validated",
        graph_name="BaP Validated",
        engine=engine,
        visibility=GraphVisibility.VALIDATED_ONLY,
    )
```

## Project Structure

```
ExposoGraph/
├── __init__.py              # Public API exports
├── _version.py              # Single source of truth for __version__
├── app.py                   # Streamlit UI orchestrator
├── branding.py              # Version and metadata
├── config.py                # App modes, graph modes, and visibility enums
├── engine.py                # NetworkX-backed graph engine
├── exporter.py              # JSON, D3.js, HTML, GEXF export
├── graph_filters.py         # Validated/exploratory graph filtering helpers
├── graph_analysis.py        # Shortest path, centrality, metabolism chains
├── grounding.py             # Canonical grounding and strict-mode preparation
├── llm_backend.py           # Pluggable LLM backends (OpenAI, Ollama)
├── llm_extractor.py         # LLM prompt + extraction pipeline
├── models.py                # Pydantic data models (Node, Edge, KnowledgeGraph)
├── reference_data.py        # Gene panels and activity scores
├── seeder.py                # DB-to-KnowledgeGraph conversion
├── storage.py               # SQLite revision history
├── flux_engine.py           # Measured kinetics plus data-driven proxy flux modeling
├── flux_cli.py              # CLI wrapper around the flux engine
├── exposure_engine.py       # NHANES/EPA IRIS exposure integration + LECR
├── exposure_cli.py          # CLI wrapper around exposure integration
├── interaction_engine.py    # Multi-carcinogen synergy & interaction matrix
├── interaction_cli.py       # CLI wrapper around the interaction engine
├── tissue_subgraphs.py      # GTEx-weighted tissue-specific subgraphs
├── population_genomics.py   # 1000G / gnomAD allele-frequency lookups
├── oxidative_stress.py      # ROS markers and antioxidant response panels
├── tk_td_modeling.py        # 1- and 2-compartment TK/TD models
├── cross_species.py         # Allometric scaling and species comparisons
├── expanded_metals.py       # Expanded heavy-metal catalogue (Wave 2)
├── wave2_classes.py         # Wave 2 carcinogen class profiles
├── unified_api.py           # PatientRiskProfile / patient_risk_query
├── population_simulation/   # Synthetic cohort generator + validation
├── data/                    # Bundled JSON parameter files
│   ├── kinetic_parameters.json
│   ├── exposure_database.json
│   ├── interaction_parameters.json
│   ├── proxy_flux_parameters.json
│   ├── proxy_flux_provenance.json
│   └── tissue_expression_data.json
├── db_clients/
│   ├── kegg.py              # KEGG REST API client
│   ├── ctd.py               # CTD chemical-gene interaction client
│   └── iarc.py              # Bundled IARC classification data
├── ui_extract.py            # Tab: LLM extraction
├── ui_manual.py             # Tab: manual node/edge entry
├── ui_preview.py            # Tab: interactive graph preview
└── ui_data.py               # Tab: raw data view
examples/
├── bap_graph.json       # Pre-built BaP metabolism graph (no API key needed)
├── build_bap_graph.py   # Demo script: load → analyze → export
└── mode_visibility_demo.py  # Demo script: strict ingestion + filtered export/save
tests/
├── test_integration.py  # End-to-end pipeline test
├── test_engine.py
├── test_models.py
├── test_exporter.py
├── test_graph_analysis.py
├── test_llm_backend.py
├── test_llm_extractor.py
├── test_db_clients.py
├── test_seeder.py
├── test_reference_data.py
├── test_config.py
└── test_storage.py
```

`kinetic_parameters.json` remains the source for literature-backed Michaelis-Menten and Hill fits. `proxy_flux_parameters.json` and `proxy_flux_provenance.json` carry the receptor-mediated and semi-quantitative proxy blocks used for classes that do not yet have full measured kinetic calibration.

## Node & Edge Types

**Nodes:** Carcinogen, Enzyme, Gene, Metabolite, DNA_Adduct, Pathway, Tissue

**Edges:** ACTIVATES, DETOXIFIES, TRANSPORTS, FORMS_ADDUCT, REPAIRS, PATHWAY, EXPRESSED_IN, INDUCES, INHIBITS, ENCODES, CUSTOM

## Development

```bash
pip install -e ".[all]"
make test                        # pytest --no-cov
make docs                        # sphinx dummy build
make lint
make typecheck
make test-cov                    # pytest with the configured 85% coverage gate
```

The repository now includes a staged GitHub Actions workflow in
`.github/workflows/ci.yml`. Required CI jobs run the regression suite and docs
build. Coverage, Ruff, and strict mypy are present as advisory jobs until the
existing backlog is reduced.

As of **April 21, 2026**, the local readiness audit is:

- `python -m pytest --no-cov` — passed (`456` tests)
- `python -m sphinx -b dummy docs docs/_build/dummy` — passed
- `python -m pytest` — failed the `85%` coverage gate (`64.96%`)
- `ruff check .` — `162` findings, all currently `E501` line-length cases
- `python -m mypy ExposoGraph` — failing across several strict-mode modules

See [docs/production-readiness.rst](docs/production-readiness.rst) for the
current release blockers and staged adoption plan.

### Optional dependency groups

| Group | Install | Provides |
|-------|---------|----------|
| `llm` | `pip install -e ".[llm]"` | OpenAI API support |
| `ollama` | `pip install -e ".[ollama]"` | Ollama local LLM support |
| `db` | `pip install -e ".[db]"` | KEGG/CTD HTTP clients |
| `streamlit` | `pip install -e ".[streamlit]"` | Streamlit web app |
| `notebook` | `pip install -e ".[notebook]"` | Jupyter + Plotly/Matplotlib |
| `dev` | `pip install -e ".[dev]"` | pytest, ruff, mypy |
| `docs` | `pip install -e ".[docs]"` | Sphinx + Furo |
| `all` | `pip install -e ".[all]"` | Everything |

## License

MIT
