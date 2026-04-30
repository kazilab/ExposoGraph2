Production Readiness
====================

As of **April 29, 2026**, ExposoGraph clears every configured quality gate and
is suitable for internal research use, and reproducible
local analysis. All gates now run as blocking jobs in CI; there is no longer
a split between required and advisory checks.

The repository ships an in-tree GitHub Actions workflow at
``.github/workflows/ci.yml`` plus a root ``Makefile`` so the quality gates are
defined locally instead of living only in documentation.

Current Gate Status
-------------------

Every gate below is green locally and blocking in CI:

- ``python -m pytest --no-cov`` — **472 tests** passing
- ``python -m sphinx -b dummy docs docs/_build/dummy`` — clean build
- ``python -m pytest`` — **71.09%** coverage vs. a **70%** floor
  (see the configured ``[tool.coverage.run] omit`` list for modules that are
  deliberately excluded)
- ``ruff check .`` — zero findings
- ``python -m mypy ExposoGraph`` — strict mode clean across 68 source files

Local Commands
--------------

Use the root ``Makefile`` to run the same checks locally:

.. code-block:: bash

   make test          # pytest --no-cov
   make test-cov      # pytest with the configured coverage floor
   make docs          # Sphinx dummy build
   make lint          # ruff check .
   make typecheck     # mypy ExposoGraph
   make ci            # every gate above, in one pass

The ``make ci-advisory`` alias is preserved for existing workflows and now runs
the same checks as ``make ci``.

Scope Decisions
---------------

The coverage floor and lint configuration encode some honest trade-offs. These
are called out here so reviewers can challenge them explicitly:

- **Coverage floor is 70%, not 85%.** The floor matches the currently tested
  surface (71.09% at the time of writing). The largest remaining gaps live in
  quantitative modules with heavy scientific content — ``exposure_engine``,
  ``tissue_subgraphs``, ``tk_td_modeling``, ``cross_species``,
  ``oxidative_stress``, ``expanded_metals``, and parts of
  ``population_simulation``. Raising the floor is a good follow-up once those
  modules grow focused tests.
- **Coverage omit list.** Streamlit/Dash UI entry points (``app.py``,
  ``_app_shared.py``, ``ui_*.py``), thin CLI launchers
  (``exposure_cli.py``, ``flux_cli.py``, ``interaction_cli.py``), and the
  population-simulation analysis/batch scripts
  (``allofus_adapter.py``, ``batch_runner.py``, ``phenotype_extractor.py``,
  ``population_analysis.py``, ``validation_framework.py``) are excluded from
  coverage. They are exercised by hand against real cohorts rather than in
  the unit-test suite.
- **Per-file Ruff ignores.** ``E501`` is suppressed in three groups of files:
  data-literal modules (reference panels, IARC tables, allele frequencies),
  scientific-formula modules (``interaction_engine``, ``exposure_engine``,
  ``flux_engine``), and orchestration / embedded-content modules
  (``unified_api``, ``exporter`` for HTML/CSS/JS, ``llm_extractor`` for prompt
  literals, ``batch_runner`` for argparse epilogs), plus the existing figure
  and UI modules. Wrapping those lines obscures the underlying equations,
  payloads, notebook figure geometry, or example commands.
  ``ExposoGraph/app.py`` keeps its ``E402`` exemption for the Streamlit
  ``sys.path`` bootstrap.

Remaining Release Blockers
--------------------------

Release-level follow-ups that remain open:

1. Revisit package metadata such as the current ``Development Status :: 3 -
   Alpha`` classifier once there is a published release cadence.
2. Add a release/publish workflow only after the GitHub repository settings
   and PyPI trusted-publisher configuration are in place.
3. Grow tests in the quantitative modules listed above so the coverage floor
   can be raised in a future iteration.
