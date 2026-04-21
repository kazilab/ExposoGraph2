``ExposoGraph.seeder``
======================

Database seeding helpers for KEGG, CTD, and IARC-backed enrichment.

KEGG and CTD seeders accept a ``mode`` argument so seeded graphs can be
prepared in ``exploratory`` or ``strict`` mode before merge.

The KEGG-backed path uses the fixed-width KEGG REST record parser from
``ExposoGraph.db_clients.kegg``. Multi-line ``GENE`` and ``PATHWAY`` sections
are supported so pathway seeding keeps the expected gene symbols and
memberships from live KEGG records.

.. automodule:: ExposoGraph.seeder
   :members:
   :undoc-members:
