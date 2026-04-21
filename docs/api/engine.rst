``ExposoGraph.engine``
======================

NetworkX-backed graph engine for building and querying the knowledge graph.

``GraphEngine.load(...)`` and ``GraphEngine.merge(...)`` both accept a
``mode`` argument:

- ``exploratory`` keeps provisional content after grounding
- ``strict`` merges only canonically grounded content and returns warnings for dropped items

The engine is backed by a NetworkX ``MultiDiGraph`` and preserves parallel
edges. Distinct evidence edges with the same source, predicate, and target are
kept as separate graph edges instead of replacing earlier records.

.. autoclass:: ExposoGraph.engine.GraphEngine
   :members:
   :undoc-members:
   :show-inheritance:
