Deployment
==========

Streamlit Cloud
---------------

1. Push the repository to GitHub.

2. Connect the repo on `Streamlit Cloud <https://share.streamlit.io>`_.

3. Set the **main file path** to ``ExposoGraph/app.py``.

4. Add any site-specific service credentials in **Secrets management**.

5. Deploy. The app reads secrets via ``st.secrets`` and falls back to
   environment variables for local development.

6. For privacy-safe public deployment, keep the default stateless mode or
   set:

   .. code-block:: toml

      ExposoGraph_MODE = "stateless"

   In stateless mode the app does not write user graphs, revisions, or HTML
   snapshots to server-side storage. Users must download the interactive HTML
   output to keep their work.

Local Development
-----------------

.. code-block:: bash

   export ExposoGraph_MODE=local
   pip install -e ".[streamlit]"
   streamlit run ExposoGraph/app.py

Or copy the example secrets file:

.. code-block:: bash

   cp .streamlit/secrets.toml.example .streamlit/secrets.toml
   # Edit .streamlit/secrets.toml with local service settings

Continuous Integration
----------------------

The repository now includes a GitHub Actions workflow at
``.github/workflows/ci.yml``. It runs the same required quality gates as the
local ``Makefile``:

- regression suite without coverage (``make test``)
- Sphinx dummy docs build (``make docs``)
- biomarker mapping validation (``make check-biomarker-mapping``)
- coverage-gated test run (``make test-cov``)
- Ruff linting (``make lint``)
- strict mypy type checking (``make typecheck``)

Run the same commands locally from the repository root:

.. code-block:: bash

   make ci
   make ci-advisory

``make ci-advisory`` is retained as a compatibility alias and currently runs
the same gates as ``make ci``.

PyPI Publishing
---------------

The project uses `hatch <https://hatch.pypa.io/>`_ as its build backend.
This checkout does not yet ship an automated publish workflow; package builds
are still created locally until GitHub repository settings and
`trusted publishers <https://docs.pypi.org/trusted-publishers/>`_ are wired up.

To build locally:

.. code-block:: bash

   pip install build
   python -m build
   # Outputs dist/exposograph-<version>-py3-none-any.whl
