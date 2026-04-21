``ExposoGraph.example_graphs``
==============================

Canonical reference-graph builders plus legacy helpers.

The canonical reference builders expose the current fully bundled reference
graph (**212 nodes / 313 edges**). The older ``build_full_legends_*`` showcase
API remains available for **107-node / 124-edge**
base example graph.

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
