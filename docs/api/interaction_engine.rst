``ExposoGraph.interaction_engine``
==================================

Quantitative multi-carcinogen interaction modeling plus provenance helpers for
the curated competitive-substrate panel.

The interaction stack is split across two JSON layers:

- ``interaction_parameters.json`` stores the numeric ``Km_uM``,
  ``Vmax_relative``, and explicit ``Ki_uM`` values used by the competitive
  inhibition model
- ``parameter_provenance.json`` stores source citations, confidence grades,
  and per-pair ``ki_status`` metadata

The interaction metadata layer also exposes a prioritized source catalog
(``BRENDA``, Rendic & Guengerich 2012, PharmGKB/ClinPGx, IARC, ATSDR, and
PubMed) together with the structured expansion backlog and its
green/yellow/red scientific-validity triage.

Public Data Accessors
---------------------

.. autofunction:: ExposoGraph.interaction_engine.get_parameter_provenance

.. autofunction:: ExposoGraph.interaction_engine.get_interaction_source_catalog

.. autofunction:: ExposoGraph.interaction_engine.get_interaction_expansion_backlog

.. autofunction:: ExposoGraph.interaction_engine.assumed_ki_pairs

Core Modeling Functions
-----------------------

.. autofunction:: ExposoGraph.interaction_engine.compute_interaction_matrix

.. autofunction:: ExposoGraph.interaction_engine.decompose_synergy

.. autofunction:: ExposoGraph.interaction_engine.monte_carlo_synergy_ci

.. autofunction:: ExposoGraph.interaction_engine.identify_critical_interactions

Primary Result Types
--------------------

.. autoclass:: ExposoGraph.interaction_engine.InteractionMatrixResult
   :members:

.. autoclass:: ExposoGraph.interaction_engine.CriticalInteraction
   :members:

.. autoclass:: ExposoGraph.interaction_engine.SynergyDecomposition
   :members:

.. autoclass:: ExposoGraph.interaction_engine.SynergyConfidenceInterval
   :members:
