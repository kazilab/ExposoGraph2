``ExposoGraph.example_graphs``
==============================

Canonical reference-graph builders plus legacy helpers.

The canonical reference builders expose the current fully bundled reference
graph (**214 nodes / 321 edges**). The ``build_full_legends_*`` showcase API
now reads the same bundled full-legends payload, so its default footprint is
also **214 nodes / 321 edges**.

The shipped D3 viewer payload in ``ExposoGraph/map/graph-data.js`` currently
matches the bundled reference graph footprint. It remains a curated graph
export rather than a rendering of the quantitative interaction-engine
competition matrix.

Data Structures
---------------

.. autoclass:: ExposoGraph.example_graphs.ArchitectureInventoryGroup
   :members:

.. autoclass:: ExposoGraph.example_graphs.ArchitectureSummary
   :members:

Canonical Reference Builders
----------------------------

.. autofunction:: ExposoGraph.example_graphs.build_reference_graph

.. autofunction:: ExposoGraph.example_graphs.build_reference_engine

.. autofunction:: ExposoGraph.example_graphs.build_reference_architecture_summary

.. autofunction:: ExposoGraph.example_graphs.write_reference_exports

Legacy Showcase Builders
------------------------

.. autofunction:: ExposoGraph.example_graphs.build_androgen_module_graph

.. autofunction:: ExposoGraph.example_graphs.build_androgen_module_engine

.. autofunction:: ExposoGraph.example_graphs.build_full_legends_graph

.. autofunction:: ExposoGraph.example_graphs.build_full_legends_engine

.. autofunction:: ExposoGraph.example_graphs.build_full_legends_architecture_summary

.. autofunction:: ExposoGraph.example_graphs.write_full_legends_exports
