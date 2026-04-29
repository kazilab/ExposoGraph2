Quick Start
===========

Installation
------------

From PyPI (once published):

.. code-block:: bash

   pip install ExposoGraph

From source (development):

.. code-block:: bash

   git clone https://github.com/kazilab/ExposoGraph.git
   cd ExposoGraph
   pip install -e ".[all]"

Optional dependency groups:

- ``streamlit`` — Streamlit UI and agraph visualization
- ``viewer`` — Dash Cytoscape advanced graph viewer
- ``notebook`` — Jupyter, Plotly, and Matplotlib
- ``dev`` — pytest, ruff, mypy
- ``docs`` — Sphinx, MyST, and Furo for documentation builds
- ``all`` — everything above

Streamlit App
-------------

.. code-block:: bash

   pip install -e ".[streamlit]"
   streamlit run ExposoGraph/app.py

App mode defaults to ``stateless``. To set it explicitly:

.. code-block:: bash

   export ExposoGraph_MODE=stateless

Set your OpenAI API key in the sidebar, or via environment variable:

.. code-block:: bash

   export OPENAI_API_KEY="sk-..."

For local persistence and revision history, switch to local mode:

.. code-block:: bash

   export ExposoGraph_MODE=local
   streamlit run ExposoGraph/app.py

Jupyter
-------

.. code-block:: bash

   pip install -e ".[notebook]"
   jupyter lab

No notebook file is currently bundled in the repository. Use the installed
package from your own notebook, or start from the runnable examples in
``examples/``.

Standalone D3.js Viewer
-----------------------

A zero-install HTML viewer ships alongside the package at
``ExposoGraph/map/index.html``, with its graph payload in
``ExposoGraph/map/graph-data.js`` (**214 nodes / 321 edges** by default — the
current bundled reference graph). This shipped payload is a curated graph
export, not the quantitative interaction engine. Open the HTML file in a
browser, or serve it statically:

.. code-block:: bash

   python -m http.server --directory ExposoGraph/map 8000
   # Open http://localhost:8000/

Viewer features:

- **Shape-coded node legend** — the *Node Types* legend renders each entry
  with the same glyph the force-directed canvas uses (diamond = Carcinogen,
  hexagon = DNA adduct, rounded rectangle = Pathway / Tissue, circle =
  Enzyme / Gene / Metabolite). The default shipped bundle currently contains
  ``Carcinogen``, ``Enzyme``, ``Metabolite``, ``DNA_Adduct``, and ``Pathway``
  nodes.
- **Resizable sidebar** — the width auto-scales with the viewport via
  ``clamp(280px, 22vw, 420px)`` and a vertical splitter between the sidebar
  and the graph can be dragged between 240 px and 640 px. Double-click the
  splitter to reset to the responsive default.
- **Runtime bundle subtitle** — the header reports the active bundle size and
  whether the heavy-metal overlay is present, without taking up sidebar space.
- **Live header counts** — the top-right counters for *nodes*, *edges*,
  *node types*, and *edge types* update dynamically as you search or apply a
  legend / carcinogen-class / tissue filter, always reflecting the currently
  highlighted subset.
- **Interactive filters** — search by label, id, group, or variant; toggle
  node types; filter by carcinogen class; and slice by tissue with a GTEx
  expression-weight threshold. The shipped bundle embeds GTEx weights on
  selected nodes; it does not currently ship standalone ``Tissue`` nodes.

To regenerate the bundled data from a Python graph:

.. code-block:: python

   from ExposoGraph import build_reference_engine, to_graph_data_js

   engine = build_reference_engine()
   to_graph_data_js(engine, "ExposoGraph/map/graph-data.js")

Advanced Viewer
---------------

.. code-block:: bash

   pip install -e ".[viewer]"

.. code-block:: python

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

Python Library
--------------

.. code-block:: python

   from ExposoGraph import (
       GraphEngine,
       GraphMode,
       GraphVisibility,
       extract_graph,
       to_json,
   )

   # LLM-powered extraction
   # exploratory keeps provisional nodes and edges
   kg = extract_graph(
       "Benzo[a]pyrene is activated by CYP1A1...",
       mode=GraphMode.EXPLORATORY,
   )
   engine = GraphEngine()
   engine.merge(kg, mode=GraphMode.EXPLORATORY)

   print(engine.node_count, "nodes")
   to_json(engine, "validated_only.json", visibility=GraphVisibility.VALIDATED_ONLY)

Interaction Provenance Helpers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from ExposoGraph import (
       assumed_ki_pairs,
       get_interaction_expansion_backlog,
       get_interaction_source_catalog,
       get_parameter_provenance,
   )

   catalog = get_interaction_source_catalog()
   backlog = get_interaction_expansion_backlog()
   provenance = get_parameter_provenance()

   print(catalog[0]["source"])  # BRENDA
   print(backlog["remaining_pairs_to_parameterize"])  # 43
   print(len(assumed_ki_pairs()))

The interaction stack keeps numeric ``Km_uM`` / ``Vmax_relative`` values in
``interaction_parameters.json``, source citations / confidence grades /
``ki_status`` in ``parameter_provenance.json``, and the structured expansion
backlog plus green/yellow/red scientific-validity triage in the interaction
metadata block.

Graph Modes
^^^^^^^^^^^

ExposoGraph uses two ingestion modes:

- ``exploratory`` keeps unmatched and custom content, annotated as provisional
- ``strict`` keeps only canonically grounded nodes and edges

.. code-block:: python

   from ExposoGraph import GraphEngine, GraphMode, extract_graph

   strict_kg = extract_graph(
       "BaP activates CYP1A1 and forms BPDE adducts",
       mode=GraphMode.STRICT,
   )

   engine = GraphEngine()
   warnings = engine.merge(strict_kg, mode=GraphMode.STRICT)
   print(warnings)

Loading Reference Gene Panels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from ExposoGraph import (
       GraphEngine,
       build_full_panel,
       get_activity_score_references,
       get_activity_scores,
   )

   # Load all 36 Tier 1 + Tier 2 genes
   kg = build_full_panel()
   engine = GraphEngine()
   engine.load(kg)

   # Look up activity scores for a gene
   scores = get_activity_scores("CYP2D6")
   for s in scores:
       print(f"  {s['allele']}: {s['value']} — {s['phenotype']}")

   refs = get_activity_score_references("CYP2D6")
   for ref in refs or []:
       print(f"  {ref['source_db']}: {ref.get('pmid') or ref.get('record_id')}")

Exporting
^^^^^^^^^

.. code-block:: python

   from ExposoGraph import (
       GraphVisibility,
       ViewerLayoutMode,
       to_gexf,
       to_graph_data_js,
       to_interactive_html,
       to_json,
       to_plotly_html,
       write_cytoscape_bundle,
   )

   # Standalone parseable app HTML
   to_interactive_html(
       engine,
       "exports/graph.html",
       visibility=GraphVisibility.ALL,
   )

   # Standalone Plotly HTML
   to_plotly_html(
       engine,
       "exports/graph_plotly.html",
       visibility=GraphVisibility.ALL,
   )

   # Validated-only HTML
   to_interactive_html(
       engine,
       "exports/graph_validated.html",
       visibility=GraphVisibility.VALIDATED_ONLY,
   )

   # Cytoscape-ready JSON bundle
   write_cytoscape_bundle(
       engine,
       "exports/graph_cytoscape.json",
       visibility=GraphVisibility.ALL,
       layout_mode=ViewerLayoutMode.COSE,
   )

   # D3.js viewer format
   to_graph_data_js(
       engine,
       "exports/graph-data.js",
       visibility=GraphVisibility.ALL,
   )

   # Plain JSON
   to_json(
       engine,
       "output.json",
       visibility=GraphVisibility.EXPLORATORY_ONLY,
   )

   # GEXF (Gephi)
   to_gexf(
       engine,
       "output.gexf",
       visibility=GraphVisibility.VALIDATED_ONLY,
   )

Filtered Revisions
^^^^^^^^^^^^^^^^^^

In ``local`` app mode, SQLite revision saves can persist either the full graph
or the current visibility slice.

.. code-block:: python

   from ExposoGraph import GraphRepository, GraphVisibility

   with GraphRepository("data/ExposoGraph.sqlite3") as repo:
       repo.save_engine(
           graph_key="bap_demo",
           graph_name="BaP Demo",
           engine=engine,
           visibility=GraphVisibility.VALIDATED_ONLY,
           note="Validated subset only",
       )

See also ``examples/mode_visibility_demo.py`` for a runnable no-API-key
example that demonstrates strict vs exploratory merge behavior and
visibility-aware export.

Biomarker Mapping Validation
----------------------------

Use the scaffold checker to keep ``ExposoGraph/data/biomarker_mapping.json``
valid, traceable, and forward update-compatible.
The preserved ``ExposoGraph/data/biomarker_mapping_old.json`` snapshot is kept
for comparison while the new JSON is rebuilt from the YAML source registry.

.. code-block:: bash

   # Validate only (no write)
   make check-biomarker-mapping

   # Rebuild the JSON from the split YAML source manifest and compare against the old snapshot
   make build-biomarker-mapping
   make compare-biomarker-mapping

   # Same check via installed console script
   exposograph-check-biomarker-mapping --mapping ExposoGraph/data/biomarker_mapping.json

   # Rebuild via installed console script
   exposograph-build-biomarker-mapping \
     --source ExposoGraph/_biomarker_scaffold/data/registries/biomarkers_master.yaml \
     --out ExposoGraph/data/biomarker_mapping.json \
     --old ExposoGraph/data/biomarker_mapping_old.json

   # Normalize and write missing trace/update fields
   python -m ExposoGraph._biomarker_scaffold.scripts.registries.check_mapping \
     --mapping ExposoGraph/data/biomarker_mapping.json --fix --write

   # Rebuild from the scaffold source manifest and compare with the old snapshot
   python -m ExposoGraph._biomarker_scaffold.scripts.registries.build_mapping \
     --source ExposoGraph/_biomarker_scaffold/data/registries/biomarkers_master.yaml \
     --out ExposoGraph/data/biomarker_mapping.json \
     --old ExposoGraph/data/biomarker_mapping_old.json
