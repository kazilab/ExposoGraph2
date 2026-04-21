"""Sphinx configuration for ExposoGraph documentation."""

import importlib.util
import os
import sys

# Add project root so autodoc can find the package
sys.path.insert(0, os.path.abspath(".."))

from ExposoGraph import APP_NAME, COPYRIGHT_HOLDER, DEVELOPED_BY, __version__

# -- Project information -----------------------------------------------------
project = APP_NAME
author = DEVELOPED_BY
copyright = f"2026, {COPYRIGHT_HOLDER}"
release = __version__

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Accept both .rst and .md
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- Options for HTML output -------------------------------------------------
html_theme = "furo" if importlib.util.find_spec("furo") else "alabaster"
html_static_path = ["_static"]
html_title = APP_NAME

# -- Extension configuration -------------------------------------------------

# autodoc
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# napoleon (Google-style docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
    "pydantic": ("https://docs.pydantic.dev/latest/", None),
}

# myst-parser
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
