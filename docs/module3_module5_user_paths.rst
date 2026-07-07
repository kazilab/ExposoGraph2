Module 3 And Module 5 User Paths
================================

ExposoGraph 2.0 keeps two neighboring quantitative workflows.

Module 3 is the simple individual-carcinogen workflow. It answers one
single-carcinogen activation and detoxification flux question at a time.

Module 5 is the advanced multi-carcinogen mechanism-resolved interaction
workflow. It answers a co-exposure question and reports mechanism-aware
interaction, warning, provenance, and transparency fields.

Module 3 remains in v2.0 and is not deprecated. Module 5 remains in v2.0 and
does not replace Module 3.

Module 3 Simple Workflow
------------------------

Use Module 3 when the analysis question is a single carcinogen class in one
tissue context.

Canonical entry points:

- ``ExposoGraph.compute_pathway_flux(...)``
- ``ExposoGraph.flux_engine.compute_pathway_flux(...)``
- ``python -m ExposoGraph.flux_cli --help``

The default workflow requires a carcinogen class, a genotype map, and a
tissue. It returns a ``PathwayFluxResult`` with activation enzymes, detox
enzymes, total activation, total detox, net activation/detoxification ratio,
susceptibility score, risk classification, warnings, and tissue/parameter
metadata. The default Module 3 path does not require an interaction context or
co-substrate exposure profile.

Runnable local example:

.. code-block:: bash

   python examples/module3_simple_workflow.py

Module 5 Advanced Workflow
--------------------------

Use Module 5 when the analysis question involves co-exposure or
multi-carcinogen interaction.

Canonical entry points:

- ``ExposoGraph.compute_interaction_matrix(...)``
- ``ExposoGraph.interaction_engine.compute_interaction_matrix(...)``
- ``python -m ExposoGraph.interaction_cli --help``

Module 5 accepts a multi-carcinogen exposure profile plus optional genotype,
tissue, lifestyle, mechanism toggles, and perturbation inputs. It reports
``InteractionMatrixResult`` fields such as independent risks,
interaction-adjusted risks, synergy matrix, GSH status, induction effects,
competitive effects, total interaction risk, interaction factor, mechanism
attribution, and mechanism-resolved risks.

Runnable local example:

.. code-block:: bash

   python examples/module5_advanced_interaction_workflow.py

Coexistence Contract
--------------------

The accepted v2.0 relationship is coexistence:

- Module 3 is retained as the direct single-carcinogen flux workflow.
- Module 5 is retained as the advanced co-exposure workflow.
- Module 3 is not described as deprecated.
- Module 5 is not described as a replacement for Module 3.
- The unified patient API may orchestrate both workflows, but orchestration
  does not collapse them into one model.

Unified API Labels
------------------

``patient_risk_query(...)`` exposes both surfaces when interaction analysis is
enabled:

- ``profile.flux_profile`` for Module 3-style flux profiling.
- ``profile.interactions`` for Module 5 interaction analysis.
- ``profile.workflow_labels["module3_simple"]`` for the stable Module 3 user
  label.
- ``profile.workflow_labels["module5_advanced"]`` for the stable Module 5
  user label.

The labels are active workflow labels, not deprecation notices.

Figure 3 Usage Note
-------------------

Figure 3 is a Module 5 figure. It is regenerated from the local source-of-truth
command recorded in the repository:

.. code-block:: bash

   python tools/generate_figure3.py

WP12 documents the command but does not regenerate the figure outputs.
