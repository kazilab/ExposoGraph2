``ExposoGraph.storage``
=======================

SQLite-backed graph revision storage for local mode.

Revision summaries carry a ``visibility`` field so saved revisions can track
whether they contain the full graph or a filtered validated/exploratory slice.
``GraphRepository`` also supports context-manager usage for clean SQLite
connection shutdown.

.. autoclass:: ExposoGraph.storage.GraphRevisionSummary
   :members:
   :undoc-members:

.. autoclass:: ExposoGraph.storage.GraphRevision
   :members:
   :undoc-members:

.. autoclass:: ExposoGraph.storage.GraphRepository
   :members:
   :undoc-members:
