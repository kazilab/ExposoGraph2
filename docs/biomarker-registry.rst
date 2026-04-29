Biomarker Registry Workflow
===========================

The biomarker mapping registry is split into source-specific YAML files under
``ExposoGraph/_biomarker_scaffold/data/registries/``. The file
``biomarkers_master.yaml`` is now a manifest that points at the source files;
it is no longer the single row store.

Current split
-------------

.. list-table::
   :header-rows: 1
   :widths: 28 18 18 36

   * - File
     - Phase
     - Tier
     - Role
   * - ``source_nhanes.yaml``
     - ``direct_measurement``
     - ``1``
     - Direct NHANES biomarker rows
   * - ``source_literature.yaml``
     - ``validated_proxy``
     - ``2``
     - Literature-backed proxy and model rows
   * - ``source_brenda.yaml``
     - ``enzyme_kinetics``
     - ``1``
     - Reserved for BRENDA Km / kinetics rows
   * - ``source_pubchem.yaml``
     - ``chemical_identity``
     - ``0``
     - Reserved for chemical identity metadata
   * - ``source_comptox.yaml``
     - ``chemical_identity``
     - ``0``
     - Reserved for chemical identity metadata

The manifest keeps the source ordering and update policy together with the
split files. The rebuilt JSON also records the source breakdown in its
top-level metadata and per-entry trace fields. Source paths are stored
relative to the repository root so the snapshot stays portable across
machines.

Phase and tier labels
---------------------

- ``registry_phase`` describes where a row belongs in the curation workflow.
- ``registry_tier`` is the maintenance confidence bucket for the source file.
- ``source_note`` should explain the scope of the file in plain language.

For this registry:

- ``direct_measurement`` means direct NHANES measurement or biomonitoring data.
- ``validated_proxy`` means literature-backed or scenario proxy rows that are
  still curated and referenced.
- ``enzyme_kinetics`` means direct Km / enzyme curation rows when those are
  added later.
- ``chemical_identity`` means identifier-only support files with no biomarker
  rows yet.

How updates work
----------------

1. Edit the source YAML that owns the biomarker row.
2. Keep each ``biomarker`` + ``lifestyle_factor`` combination in exactly one
   source file.
3. If the row changes source family, update the manifest entry and the source
   file metadata together.
4. Add new source families by creating a new ``source_*.yaml`` file and adding
   it to ``biomarkers_master.yaml``.
5. Rebuild the JSON snapshot with ``make build-biomarker-mapping``.
6. Compare the new snapshot against the preserved old file with
   ``make compare-biomarker-mapping``.
7. Validate and test with ``make check-biomarker-mapping`` and
   ``python -m pytest --no-cov tests/test_biomarker_mapping.py
   tests/test_biomarker_mapping_build.py -q``.

The preserved comparison file is ``ExposoGraph/data/biomarker_mapping_old.json``.
The rebuilt snapshot is ``ExposoGraph/data/biomarker_mapping.json``.
