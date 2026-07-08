Module 5 Mechanism-Resolved Model
=================================

Module 5 is the advanced multi-carcinogen mechanism-resolved interaction
workflow in ExposoGraph 2.0. It combines exposure-scaled baseline risk,
induction, inhibition, GSH redox status, adjusted-risk reporting, and
mechanism attribution into transparent local outputs.

Inhibition Modes
----------------

Where implemented by the accepted v2 code path, Module 5 uses a centralized
reversible-inhibition resolver for:

- competitive inhibition;
- pure non-competitive inhibition;
- uncompetitive inhibition;
- mixed inhibition.

The model does not silently convert every incomplete kinetic record into a
quantified effect. Missing or mode-incomplete kinetic evidence is surfaced
through warnings, provenance, status fields, and review-required flags.

Ki, Km, IC50, Warning, And Provenance Behavior
----------------------------------------------

Module 5 uses centralized Ki/Km/IC50 handling. Curated Ki values are preferred.
Permitted local or proxy values are marked with their source kind, resolution
method, uncertainty, warnings, and metadata. IC50 conversion requires an
explicit inhibition mode and assay context. A single IC50 is not enough to
resolve mixed inhibition without both required arms.

The v2 output contract keeps warning and review information visible in
compatibility payloads, mechanism-resolved risk records, model-card fields, and
unified API biological-output integration.

GSH Coupling And Adjusted Risk
------------------------------

Module 5 uses activation-scaled GSH contribution where upstream activation
information is available. If direct upstream activation is absent, the model can
use explicit fallback inputs when supplied, or neutral scaling with a warning.

The implemented adjusted-risk structure is:

.. code-block:: text

   baseline_relative_risk = BASELINE_RISK_SCORES[carcinogen] * exposure_multiplier
   final_mechanism_multiplier = induction_multiplier * inhibition_burden_multiplier * matrix_gsh_penalty
   adjusted_relative_risk = baseline_relative_risk * final_mechanism_multiplier

Diagnostic biological-output records are not authoritative adjusted risk. They
are included for inspection and transparency, while authoritative adjusted risk
comes from the matrix-level mechanism-resolved path.

Synergy And Mechanism Attribution
---------------------------------

Module 5 reports Shapley/eight-state mechanism attribution where implemented.
That attribution is the authoritative mechanism-attribution interpretation for
induction, competition, and GSH state contributions.

Pairwise ``synergy_matrix`` output is retained as descriptive heatmap and
viewer-compatible information. Residual fields are numerical reconstruction
checks, not an additional unexplained mechanism.

Figure 3 Source Of Truth
------------------------

Figure 3 is a Module 5 figure and uses the current local interaction and
synergy/decomposition outputs. The accepted source-of-truth regeneration
command is:

.. code-block:: bash

   python tools/generate_figure3.py

The command writes the notebook and CSV/PNG/PDF/SVG outputs under
``Figures_Notebook/``. Figure regeneration was completed before WP12; WP12 does
not rerun it.

Heavy-Metal And v3 Boundary
---------------------------

Schema-compatible heavy-metal content is included where present in the accepted
v2 baseline. ExposoGraph 2.0 does not claim new HIF1A/VEGF gene-level
treatment, cobalt-HIF-VEGF pathway modeling, expanded edge ontology, population
modeling, population-genomics expansion, All of Us adapter expansion, broad
text-extraction adapter expansion, full Cytoscape/Dash restoration beyond
currently supported optional helpers, full migration of every runtime JSON
parameter into the knowledge graph, or new external data-source integrations.
Those items are deferred to v3 or later unless already explicitly supported by
the accepted v2 package.

Accepted v2 Model Boundaries
----------------------------

- Module 3 remains the simple individual-carcinogen workflow.
- Module 5 remains the advanced multi-carcinogen mechanism-resolved workflow.
- Module 5 warnings, provenance, model-card fields, and review-required flags
  are part of the user-facing contract.
- GSH v2 behavior is a bounded local model with visible fallbacks and warnings.
- Pairwise synergy heatmaps are descriptive; Shapley/eight-state attribution is
  the mechanism-attribution interpretation where implemented.
