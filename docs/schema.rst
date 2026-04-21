Schema Reference
================

ExposoGraph uses a typed schema for all nodes and edges in the knowledge graph.
Types are defined as Python enums and enforced by Pydantic models.
The core ontology remains fixed, while matching/provenance metadata captures
whether a record is canonical, alias-matched, unmatched, or custom.

Node Types
----------

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Type
     - Description
     - Key Fields
   * - ``Carcinogen``
     - Chemical carcinogenic agents (PAHs, HCAs, nitrosamines, etc.)
     - ``group``, ``iarc``
   * - ``Enzyme``
     - Metabolizing, transport, and repair proteins
     - ``phase`` (I/II/III only), ``role``, ``group``, ``tier``
   * - ``Gene``
     - Gene loci — for pharmacogenomic variants or expression context
     - ``variant``, ``phenotype``, ``activity_score``
   * - ``Metabolite``
     - Intermediate and terminal metabolites
     - ``reactivity``
   * - ``DNA_Adduct``
     - Covalent DNA lesions formed by reactive metabolites
     -
   * - ``Pathway``
     - Biological/KEGG pathways
     -
   * - ``Tissue``
     - Anatomical tissues or organs
     -

Edge Types
----------

.. list-table::
   :header-rows: 1
   :widths: 22 55 23

   * - Type
     - Description
     - Typical Source → Target
   * - ``ACTIVATES``
     - Enzyme converts procarcinogen to reactive metabolite
     - Enzyme → Metabolite
   * - ``DETOXIFIES``
     - Enzyme conjugates/inactivates a metabolite
     - Enzyme → Metabolite
   * - ``TRANSPORTS``
     - Efflux transporter moves conjugate out of cell
     - Enzyme → Metabolite
   * - ``FORMS_ADDUCT``
     - Reactive metabolite covalently modifies DNA
     - Metabolite → DNA_Adduct
   * - ``REPAIRS``
     - DNA repair enzyme removes a lesion
     - Enzyme → DNA_Adduct
   * - ``PATHWAY``
     - Node belongs to a biological pathway
     - Node → Pathway
   * - ``EXPRESSED_IN``
     - Gene or enzyme is expressed in a tissue
     - Gene/Enzyme → Tissue
   * - ``INDUCES``
     - Substance induces enzyme expression/activity
     - Carcinogen → Enzyme
   * - ``INHIBITS``
     - Substance inhibits enzyme expression/activity
     - Carcinogen → Enzyme
   * - ``ENCODES``
     - Gene encodes an enzyme
     - Gene → Enzyme
   * - ``CUSTOM``
     - User-defined exploratory predicate
     - Any validated or provisional node pair

Annotation Fields
-----------------

All nodes support these optional annotation fields:

- ``source_db`` — Provenance databases (NCBI Gene, GTEx, ClinPGx, CTD, IARC, KEGG, etc.)
- ``evidence`` — Brief evidence note
- ``pmid`` — PubMed ID
- ``tissue`` — Relevant tissue context

For repair proteins, ``group`` is the recommended place to store the repair
class (for example ``DNA Repair (BER)``, ``DNA Repair (NER)``, or
``DNA Repair (Direct Reversal)``). ``phase`` is reserved for Phase I/II/III
metabolism and transport labels.

Edges also support ``carcinogen`` (the parent carcinogen context node ID)
and ``label`` (short description of the reaction).

Structured Provenance
---------------------

Nodes and edges may also carry a ``provenance`` list with one or more records.
Each record can store:

- ``source_db`` — Source database or catalog
- ``record_id`` — Stable database identifier or accession when available
- ``evidence`` — Evidence summary for the specific record
- ``pmid`` — PubMed ID
- ``tissue`` — Tissue-specific context
- ``citation`` — Human-readable citation text
- ``url`` — Link to the source when available

Legacy top-level fields such as ``source_db`` and ``pmid`` remain supported.
When present, they are normalized into a single provenance record for backward
compatibility.

Grounding and Match Metadata
----------------------------

Nodes and edges also support a lightweight grounding layer:

- ``origin`` — where the record came from: ``imported``, ``seeded``, ``user``, or ``llm``
- ``match_status`` — grounding state: ``unknown``, ``canonical``, ``alias``, ``unmatched``, or ``custom``
- ``canonical_id`` / ``canonical_label`` / ``canonical_namespace`` — canonical mapping for grounded nodes
- ``canonical_predicate`` / ``canonical_namespace`` — canonical mapping for grounded edges
- ``custom_type`` / ``custom_predicate`` — exploratory labels for intentionally user-defined content

``CUSTOM`` nodes must include ``custom_type`` and ``CUSTOM`` edges must include
``custom_predicate``.

This metadata is orthogonal to curation. A node may be canonical but still in
``Draft`` curation status, or custom but manually reviewed.

Curation Fields
---------------

Nodes and edges may include a ``curation`` object with review metadata:

- ``status`` — ``Draft``, ``In Review``, ``Reviewed``, ``Approved``, or ``Rejected``
- ``confidence`` — ``Low``, ``Medium``, or ``High``
- ``curator`` — Person who created or updated the record
- ``reviewed_by`` — Reviewer identity
- ``reviewed_at`` — Review timestamp or date string
- ``notes`` — Free-text curation rationale

Validation
----------

The :class:`~ExposoGraph.engine.GraphEngine` enforces referential integrity and
mode-aware merge behavior:

- Every edge ``source`` and ``target`` must reference an existing node
- If ``carcinogen`` is set on an edge, that node must also exist
- :meth:`~ExposoGraph.engine.GraphEngine.load` and
  :meth:`~ExposoGraph.engine.GraphEngine.merge` accept ``mode="exploratory"``
  or ``mode="strict"``
- ``strict`` mode drops non-canonical nodes and edges and returns warning strings
