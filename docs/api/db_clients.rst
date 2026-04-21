``ExposoGraph.db_clients``
==========================

Public database clients for KEGG, CTD, and the bundled IARC reference catalog.

KEGG notes
----------

The KEGG client parses the fixed-width KEGG REST record format used by
``get/{id}`` responses. This includes:

- multi-line ``GENE`` sections where the first token may be a numeric KEGG gene ID
- multi-line ``PATHWAY`` sections in gene records

Those parsing rules are important for seeding because they preserve pathway
member gene symbols and per-gene pathway memberships from live KEGG records.

.. automodule:: ExposoGraph.db_clients.kegg
   :members:
   :undoc-members:

.. automodule:: ExposoGraph.db_clients.ctd
   :members:
   :undoc-members:

.. automodule:: ExposoGraph.db_clients.iarc
   :members:
   :undoc-members:
