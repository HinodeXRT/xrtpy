# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import os

from datetime import datetime
from sphinx.application import Sphinx

# -- Project information -----------------------------------------------------

project = "xrtpy"
author = "Joy Velasquez, Nick Murphy, and Jonathan Slavin"
copyright = f"2021–{datetime.utcnow().year}, {author}"

# The full version, including alpha/beta/rc tags
# from xrtpy import __version__
# release = __version__

release = "0.2.0dev"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_automodapi.automodapi",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.graphviz",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "nbsphinx",
    "sphinx_copybutton",
    "sphinx_gallery.load_style",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_changelog",
    "sphinxcontrib.bibtex",
    "hoverxref.extension",
]

bibtex_bibfiles = ["bibliography.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "author_year"
bibtex_cite_id = "{key}"

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    ".DS_Store",
    "_build",
    "_links.rst",
    "_substitutions.rst",
    "Thumbs.db",
]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ".rst"

# The root toctree document.
root_doc = "index"

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "ndcube": ("https://docs.sunpy.org/projects/ndcube/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
    "sunpy": ("https://docs.sunpy.org/en/stable/", None),
}

hoverxref_intersphinx = [
    "astropy",
    "ndcube",
    "numpy",
    "python",
    "scipy",
    "sunpy",
]

autoclass_content = "both"
autodoc_typehints_format = "short"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

linkcheck_allowed_redirects = {
    r"https://doi\.org/.+": r"https://.+",  # DOI links are more persistent
    r"https://docs.+\.org": r"https://docs.+\.org/en/.+",
    r"https://docs.+\.io": r"https://docs.+\.io/en/.+",
    r"https://docs.+\.com": r"https://docs.+\.com/en/.+",
    r"https://docs.+\.dev": r"https://docs.+\.dev/en/.+",
    r"https://.+\.readthedocs\.io": r"https://.+\.readthedocs\.io/en/.+",
    r"https://.+/github\.io": r"https://.+/github\.io/en/.+",
    r"https://pip\.pypa\.io": r"https://pip\.pypa\.io/en/.+",
    r"https://www.python.org/dev/peps/pep.+": "https://peps.python.org/pep.+",
}

linkcheck_anchors = True
linkcheck_anchors_ignore = [
    "/room",
    r".+openastronomy.+",
    "L[0-9].+",
    "!forum/plasmapy",
]

# Use a code highlighting style that meets the WCAG AA contrast standard
pygments_style = "default"

# adapted from sphinx-hoverxref conf.py
if os.environ.get("READTHEDOCS"):
    # Building on Read the Docs
    hoverxref_api_host = "https://readthedocs.org"

    if os.environ.get("PROXIED_API_ENDPOINT"):
        # Use the proxied API endpoint
        # - A RTD thing to avoid a CSRF block when docs are using a
        #   custom domain
        hoverxref_api_host = "/_"

hoverxref_tooltip_maxwidth = 600  # RTD main window is 696px
hoverxref_auto_ref = True
hoverxref_mathjax = True
hoverxref_sphinxtabs = True

# hoverxref has to be applied to these
hoverxref_domains = ["py", "cite"]
hoverxref_roles = ["confval", "term"]

hoverxref_role_types = {
    # roles with cite domain
    "p": "tooltip",
    "t": "tooltip",
    #
    # roles with py domain
    "attr": "tooltip",
    "class": "tooltip",
    "const": "tooltip",
    "data": "tooltip",
    "exc": "tooltip",
    "func": "tooltip",
    "meth": "tooltip",
    "mod": "tooltip",
    "obj": "tooltip",
    #
    # roles with std domain
    "confval": "tooltip",
    "hoverxref": "tooltip",
    "ref": "tooltip",
    "term": "tooltip",
}

# Specify patterns to ignore when doing a nitpicky documentation build.

python_role = "py:.*"

nitpick_ignore_regex = [
    (python_role, "and"),
    (python_role, "array .*"),
    (python_role, "array_like"),
    (python_role, "callable"),
    (python_role, "function"),
    (python_role, ".*integer.*"),
    (python_role, "iterable"),
    (python_role, "key"),
    (python_role, "keyword-only"),
    (python_role, ".* object"),
    (python_role, "optional"),
    (python_role, "or"),
    (python_role, ".*real number.*"),
    (python_role, ".*Unit.*"),
]

# This is added to the end of RST files — a good place to put substitutions to
# be used globally.
rst_epilog = ""
for epilog_file in ["_links.rst", "_substitutions.rst"]:
    with open(epilog_file) as file:
        rst_epilog += file.read()


def setup(app: Sphinx) -> None:
    app.add_config_value("revision", "", True)
    app.add_css_file("css/admonition_color_contrast.css")
    app.add_css_file("css/plasmapy.css", priority=600)
