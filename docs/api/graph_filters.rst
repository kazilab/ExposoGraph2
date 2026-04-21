``ExposoGraph.graph_filters``
=============================

Helpers for deriving validated-only or exploratory-only graph views without
mutating the stored graph.

The filtering helpers return detached ``Node`` and ``Edge`` copies, so editing
the filtered result does not mutate the source ``KnowledgeGraph`` or the
backing engine state.

.. automodule:: ExposoGraph.graph_filters
   :members:
   :undoc-members:
