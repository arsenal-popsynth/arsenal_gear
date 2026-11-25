# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath("../arsenal_gear"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "arsenal_gear"
project_copyright = "2024, BW Keller"
author = "BW Keller"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "myst_parser",
    "sphinxcontrib.mermaid",
]

autosummary_generate = True
autodoc_typehints_format = "short"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for MyST -------------------------------------------------
# https://myst-parser.readthedocs.io/en/latest/sphinx/intro.html

myst_enable_extensions = ["dollarmath", "amsmath"]

# -- Options for sphinxcontrib-mermaid ---------------------------------------
# https://github.com/mgaitan/sphinxcontrib-mermaid
# The default is a fixed 500px, which shrinks larger diagrams until their text
# is hard to read; "auto" scales them to the page width instead.

mermaid_height = "auto"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "classic"
html_static_path = ["_static"]
html_sidebars = {
    "**": ["globaltoc.html", "relations.html", "sourcelink.html", "searchbox.html"]
}


import arsenal_gear

version = ".".join(arsenal_gear.__version__.split(".")[:2])

release = arsenal_gear.__version__
