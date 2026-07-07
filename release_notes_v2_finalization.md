# ExposoGraph 2.0 Finalization Release Notes Draft

This draft summarizes the v2 finalization work for human review before release
publication.

## Implemented In v2

- Module 3 remains the simple individual-carcinogen flux workflow through
  `compute_pathway_flux`.
- Module 5 remains the advanced multi-carcinogen mechanism-resolved interaction
  workflow through `compute_interaction_matrix`.
- Module 5 exposes inhibition-mode-aware handling for competitive, pure
  non-competitive, uncompetitive, and mixed inhibition where implemented.
- Module 5 uses centralized Ki/Km/IC50 resolution with visible warnings,
  provenance, uncertainty/status metadata, and review-required behavior.
- Module 5 reports activation-scaled GSH contribution where available, explicit
  fallbacks where needed, and adjusted-risk fields with diagnostic-output
  boundaries.
- Module 5 exposes Shapley/eight-state mechanism attribution, while retaining
  pairwise synergy values as descriptive heatmap output.
- Figure 3 has a local source-of-truth regeneration command:
  `python tools/generate_figure3.py`.
- Package-root public helpers expose the reference graph, reference engine,
  Module 3 workflow, Module 5 workflow, and unified patient API surfaces.

## Validated In v2

- Module 3 public contract and JSON-safe output behavior.
- Module 5 public contract, inhibition parameter resolution, GSH coupling,
  adjusted-risk output, model-card fields, and biological-output transparency.
- Module 3 and Module 5 coexistence without deprecation or replacement.
- KG/provider local data paths and package-root graph helper access.
- API/CLI/UI user pathways where applicable.
- Figure 3 source-of-truth generation from the integrated Module 5 outputs.
- Local platform integration through `tools/local_platform_integration_check.py`.

## Deferred To v3 Or Later

- HIF1A/VEGF gene-level knowledge-graph expansion.
- Cobalt-HIF-VEGF pathway modeling.
- Expanded edge ontology design.
- Population modeling.
- Population-genomics expansion.
- All of Us adapter expansion.
- Broad text-extraction adapter expansion.
- Full Cytoscape/Dash restoration beyond currently supported optional helpers.
- Full migration of every runtime JSON parameter into the knowledge graph.
- New external data-source integrations not required for v2 finalization.

## Accepted v2 Model Boundaries

- Module 3 is not deprecated in v2.0.
- Module 5 does not replace Module 3 in v2.0.
- Schema-compatible heavy-metal content is included where present in the
  accepted baseline, but HIF1A/VEGF pathway expansion is deferred.
- Diagnostic biological outputs support inspection and transparency; they are
  not authoritative adjusted-risk multipliers.
- Pairwise synergy heatmaps are descriptive. Shapley/eight-state attribution is
  the mechanism-attribution interpretation where implemented.
- GSH v2 behavior is a bounded local model with explicit fallbacks, warnings,
  and provenance.
