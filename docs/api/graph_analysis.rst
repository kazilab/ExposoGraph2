``ExposoGraph.graph_analysis``
==============================

Domain-aware graph analysis helpers built on top of :class:`ExposoGraph.engine.GraphEngine`.

Notable behavior:

- ``metabolism_chain()`` follows carcinogen-linked metabolism edges without
  pulling in unrelated unlabeled branches that merely share an upstream enzyme
- ``variant_impact_score()`` combines representative activity scores with
  downstream adduct and repair topology

.. automodule:: ExposoGraph.graph_analysis
   :members:
   :undoc-members:
