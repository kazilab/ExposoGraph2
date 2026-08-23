# Reference Map viewer refactor: single source of truth + Python-side filtering

Status: approved, in progress on `refactor_ui` (forked from `qa`)

## Motivation

The "Reference Map" tab (`ui_map_viewer.py`) currently bundles a hand-maintained
`ExposoGraph/map/graph-data.js` and does all graph filtering in JavaScript
inside `ExposoGraph/map/index.html`. This has caused real drift: `graph-data.js`
(278 nodes / 461 edges, legacy `ACTIVATES`/`DETOXIFIES`/... edge schema, no
`Substrate`/`Receptor` nodes) is materially behind `ExposoGraph/map/graph-data.json`
(327 nodes / 563 edges, current directional edge schema, includes `Substrate`/
`Receptor`). Filtering logic embedded in `index.html` is also not reusable
outside the browser, and has real bugs (filters reset each other instead of
intersecting; the carcinogen filter only expands one hop instead of following
full metabolic paths; the tissue filter hides everything instead of just
de-emphasizing under-expressed enzymes).

This doc records the agreed plan before making any code changes.

## Key findings

- `engine.py` already has `load_reference_graph()`, which loads the canonical
  `graph-data.json` (plus tissue-expression and interaction-parameter overlays).
  It is not currently called anywhere in the running app.
- `_app_shared.start_engine()` is dead/broken: it calls a nonexistent
  `engine.load_base_data()` against a nonexistent `map/graph_data.json`, and
  `app.py` imports `start_engine` but never calls it. Every tab besides the
  Reference Map is silently running against an empty `GraphEngine` today.
- `engine.py` already has several unused-but-ready building blocks added in
  prior commits (`34a3afe`, `d9e0217`): `get_tissues`, `get_carcinogen_groups`,
  `carcinogens_by_group`, `filtered_subgraph`, `enzymes_by_tissue_threshold`,
  and `paths_from_carcinogen`/`paths_to_carcinogen` (full maximal simple
  directed path enumeration). None have call sites outside `engine.py`/tests.
- There are two independent D3 rendering paths: `ui_map_viewer.py`
  (static bundle of `map/index.html` + `map/graph-data.js`, no engine
  involved) and `ui_d3_viewer.py` (`exporter.to_interactive_html_string`,
  engine-driven, different template). The mechanism the second path uses --
  `exporter._inline_graph_data_script()` swapping the `<script src="./graph-data.js">`
  tag in a template for an inline `<script>const GRAPH_DATA = {...}</script>` --
  is what we reuse for the map tab, fed a Python-filtered subgraph instead of
  the raw engine dump.
- Current `index.html` bugs: a single global `activeView` means switching
  filter kind (search/type/carcinogen/tissue) always resets the others first
  -- no intersection is possible. `applyCarcinogenFilter` only expands one hop
  (plus a hardcoded second hop for `FORMS_ADDUCT`/`REPAIRS`), not full directed
  paths. `applyTissueFilter` fades out everything not touching a qualifying
  enzyme rather than leaving the rest of the graph untouched.

## Decisions (confirmed)

1. **Ask #1 (single source of truth):** wire the engine's existing
   `load_reference_graph()` into app startup; retire `graph-data.js` as an
   input (it may still be regenerated as an export artifact via
   `exporter.to_graph_data_js`, never read as input again).
2. **Asks #2/#3 (Python-side filtering, engine-owned):** new `GraphEngine`
   methods, operating on plain node/edge dicts so they're reusable outside
   Streamlit (and, later, by `flux_engine.py` -- not part of this change).
3. **Ask #4 (keep the force-directed layout):** `index.html` is trimmed, not
   rewritten -- the D3 simulation, zoom/pan, drag, tooltip, and detail panel
   stay as-is. Only the filter-computation functions and the controls that
   drove them are removed, replaced by Streamlit widgets driving Python
   filtering.
4. **Ask #5 (multi-select intersection):** node-type filter and carcinogen
   class/group filter intersect (AND). Node type filter is a strict
   membership filter -- unselected/excluded types are absent from the
   rendered subgraph, not dimmed.
5. **Ask #6 (carcinogen + tissue filter semantics):**
   - Carcinogen filtering stays at the class/group level (not per-carcinogen).
     For each selected group, resolve its Carcinogen node ids, take the union
     of `paths_from_carcinogen` over all of them, and union the resulting
     nodes/edges across groups.
   - Tissue filtering never removes anything. It is a non-destructive overlay
     applied *after* the node-type/carcinogen intersection: Enzyme nodes (and
     only their directly incident edges) in the already-shown subgraph get a
     `_dimmed` annotation when their `tissue_weights[tissue]` is below the
     threshold. A dimmed enzyme that has already been excluded by the
     node-type filter is simply absent, not dimmed -- dimming only ever
     applies to enzymes already present in the shown subgraph.
6. **Ask #7 (hide Substrates):** `Substrate` nodes are excluded by default in
   the new node-type subgraph helper, so every caller gets this for free; they
   remain untouched in the underlying `GraphEngine`/knowledge graph.

## New `GraphEngine` methods (engine.py)

- `carcinogen_group_paths_subgraph(groups)` -- resolves groups to Carcinogen
  node ids via the existing `carcinogens_by_group`, unions
  `paths_from_carcinogen` over each, returns `{"nodes": [...], "edges": [...]}`.
- `subgraph_by_node_types(node_types, *, exclude_types=("Substrate",))` --
  generalizes the existing single-type `nodes_by_type`/`subgraph_by_node_type`
  to multi-select, with the Substrate exclusion baked in as the default.
- `map_viewer_subgraph(node_types=None, carcinogen_groups=None)` -- intersects
  the two axes above (`None`/empty axis imposes no restriction, matching the
  identity-element convention `filtered_subgraph` already uses); always drops
  Substrates.
- `dim_by_tissue_threshold(subgraph, tissue, threshold)` -- annotates a given
  subgraph's Enzyme nodes/incident edges with `_dimmed` based on
  `tissue_weights`, leaving every other node/edge untouched. Returns a new
  subgraph dict (non-mutating).

Existing `graph_filtering.py`/`graph_filters.py` (`GraphFilterCriteria`,
`filter_graph_by_criteria`) are left as-is -- they serve a distinct, tested
concern (validated/exploratory export visibility) with different (strict-AND,
non-path-based) semantics, and aren't the mechanism ask #3 is pointed at.

## Exporter change (exporter.py)

Add a helper that turns an arbitrary `{"nodes", "edges"}` dict (rather than a
whole `GraphEngine`) into the inlined-`GRAPH_DATA` HTML string, reusing the
existing `_inline_graph_data_script`/template-loading machinery. `ui_map_viewer.py`
uses this against the existing `map/index.html` template.

## `ui_map_viewer.py` / `index.html`

`ui_map_viewer.py` becomes: Streamlit multiselects for node type and
carcinogen group, a select + slider for tissue/threshold; calls the new
engine methods in sequence; renders via the new exporter helper.

`index.html` keeps: D3 force simulation (`forceLink`/`forceManyBody`/
`forceCenter`/`forceCollide`/`forceX`/`forceY`), `ticked()`, zoom/pan, drag,
tooltip, detail panel, and (kept as client-side, non-removing) free-text
search highlight -- search isn't one of the seven asks and doesn't remove
data, so it stays as lightweight JS for keystroke responsiveness.

`index.html` loses: `applyTypeHighlight`/`highlightType`,
`applyCarcinogenFilter`/`filterByCarcinogen`, `applyTissueFilter`/
`filterByTissue`, `activeView`/`restoreActiveView`, and the sidebar controls
that drove them. A small addition renders `_dimmed` nodes/edges at reduced
opacity using the existing styling constants.

## Files touched

- `ExposoGraph/engine.py` -- new methods above.
- `ExposoGraph/_app_shared.py` -- fix `start_engine`.
- `ExposoGraph/app.py` -- call `start_engine(engine)`.
- `ExposoGraph/exporter.py` -- new dict-based render helper.
- `ExposoGraph/ui_map_viewer.py` -- rewritten.
- `ExposoGraph/map/index.html` -- trimmed.
- `ExposoGraph/map/graph-data.js` -- no longer read as input.
- `tests/test_engine.py` -- coverage for the new methods.

## Sequencing

1. Fix `start_engine`/`app.py` bootstrap (ask #1) -- isolated, verifiable on
   its own, unblocks the other tabs too.
2. Add the new `GraphEngine` methods with unit tests (asks #3, #6, #7).
3. Add the exporter dict-based render helper.
4. Rewrite `ui_map_viewer.py` + trim `index.html` together (asks #2, #4, #5),
   since they must land together to keep the tab working.
