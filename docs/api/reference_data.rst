``ExposoGraph.reference_data``
==============================

Curated gene panels and activity scores from the CarcinoGenomic Platform.

Gene Panel Data
---------------

.. data:: ExposoGraph.reference_data.TIER1_GENES

   List of 13 Tier 1 core carcinogen-metabolizing enzyme dictionaries.

.. data:: ExposoGraph.reference_data.TIER2_GENES

   List of 23 Tier 2 extended gene panel dictionaries.

.. data:: ExposoGraph.reference_data.CURATION_SOURCE_MANIFEST

   Structured curation-source manifest describing the manuscript-aligned
   primary sources and the implementation-level supporting sources used by
   ExposoGraph.

.. data:: ExposoGraph.reference_data.REFERENCE_KEGG_PATHWAYS

   Curated KEGG pathway catalog used to track the reference pathway context
   for carcinogen metabolism and DNA-damage interpretation.

.. data:: ExposoGraph.reference_data.ACTIVITY_SCORES

   Dictionary mapping gene IDs to lists of per-allele activity score entries.
   Each entry has ``allele``, ``value``, ``phenotype``, and ``confidence`` keys.

.. data:: ExposoGraph.reference_data.ACTIVITY_SCORE_METADATA

   Dictionary mapping gene IDs to evidence metadata for each activity score
   table, including ``evidence_basis``, ``note``, and supporting
   ``references``.

Panel Builders
--------------

.. autofunction:: ExposoGraph.reference_data.build_tier1_panel

.. autofunction:: ExposoGraph.reference_data.build_tier2_panel

.. autofunction:: ExposoGraph.reference_data.build_full_panel

Lookups
-------

.. autofunction:: ExposoGraph.reference_data.get_activity_scores

.. autofunction:: ExposoGraph.reference_data.get_activity_score_metadata

.. autofunction:: ExposoGraph.reference_data.get_activity_score_references
