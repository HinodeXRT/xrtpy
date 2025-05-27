"""
Configuration file for the Sphinx documentation builder.
"""

# -- stdlib imports ------------------------------------------------------------
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path

from packaging.version import Version

# -- Read the Docs Specific Configuration --------------------------------------
# This needs to be done before xrtpy is imported
on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if on_rtd:
    os.environ["SUNPY_CONFIGDIR"] = "/home/docs/"
    os.environ["HOME"] = "/home/docs/"
    os.environ["LANG"] = "C"
    os.environ["LC_ALL"] = "C"
    os.environ["PARFIVE_HIDE_PROGRESS"] = "True"

# -- Imports -------------------------------------------------------------------
from astropy.utils.exceptions import AstropyDeprecationWarning
from matplotlib import MatplotlibDeprecationWarning
from sunpy.util.exceptions import (
    SunpyDeprecationWarning,
    SunpyPendingDeprecationWarning,
)

# -- Project information -------------------------------------------------------
project = "XRTpy"
author = "Joy Velasquez, Nick Murphy, and Jonathan Slavin"
copyright = f"2021-{datetime.now(tz=timezone.utc).year}, {author}"

# The full version, including alpha/beta/rc tags
from xrtpy import __version__

_version_ = Version(__version__)
# NOTE: Avoid "post" appearing in version string in rendered docs
if _version_.is_postrelease:
    version = release = f"{_version_.major}.{_version_.minor}.{_version_.micro}"
else:
    version = release = str(_version_)
is_development = _version_.is_devrelease
# -- General configuration ---------------------------------------------------

# We want to make sure all the following warnings fail the build
warnings.filterwarnings("error", category=SunpyDeprecationWarning)
warnings.filterwarnings("error", category=SunpyPendingDeprecationWarning)
warnings.filterwarnings("error", category=MatplotlibDeprecationWarning)
warnings.filterwarnings("error", category=AstropyDeprecationWarning)

# For the linkcheck
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
linkcheck_anchors_ignore = []

# sphinxext-opengraph
ogp_image = "https://raw.githubusercontent.com/HinodeXRT/xrtpy/main/docs/_static/images/xrtpy_logo.png"
ogp_use_first_image = True
ogp_description_length = 160
ogp_custom_meta_tags = [
    '<meta property="og:ignore_canonical" content="true" />',
]

# Suppress warnings about overriding directives as we overload some of the
# doctest extensions.
suppress_warnings = [
    "app.add_directive",
]

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx_copybutton",
    "sphinx_gallery.gen_gallery",
    "sphinx_issues",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinxcontrib.bibtex",
    "sphinxext.opengraph",
]

# Set automodapi to generate files inside the generated directory
automodapi_toctreedirnm = "generated/api"

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
# html_extra_path = ['robots.txt']

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "_links.rst",
    "_substitutions.rst",
]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# The reST default role (used for this markup: `text`) to use for all
# documents. Set to the "smart" one.
default_role = "obj"

# Disable having a separate return type row
napoleon_use_rtype = False

# Disable google style docstrings
napoleon_google_docstring = False

# Disable the use of param, which prevents a distinct "Other Parameters" section
napoleon_use_param = False

# Enable nitpicky mode, which forces links to be non-broken
nitpicky = True
# This is not used. See docs/nitpick-exceptions file for the actual listing.
nitpick_ignore = []
with Path("nitpick-exceptions").open() as f:
    for line in f:
        if line.strip() == "" or line.startswith("#"):
            continue
        dtype, target = line.split(None, 1)
        target = target.strip()
        nitpick_ignore.append((dtype, target))

bibtex_bibfiles = ["bibliography.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "author_year"
bibtex_cite_id = "{key}"

# This is added to the end of RST files — a good place to put substitutions to be used globally.
rst_epilog = ""
for epilog_file in ["_links.rst", "_substitutions.rst"]:
    with Path(epilog_file).open() as file:
        rst_epilog += file.read()

# Configure sphinx-issues
issues_github_path = "HinodeXRT/xrtpy"

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://numpy.org/doc/stable/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/numpy.inv"),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/scipy.inv"),
    ),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "ndcube": ("https://docs.sunpy.org/projects/ndcube/en/stable/", None),
    "sunpy": ("https://docs.sunpy.org/en/stable/", None),
}

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "logo": {
        "text": "XRTpy",
    },
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/HinodeXRT/xrtpy/",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/xrtpy/",
            "icon": "fa-brands fa-python",
        },
    ],
}
html_context = {
    "github_user": "HinodeXRT",
    "github_repo": "xrtpy",
    "github_version": "main",
    "doc_path": "docs",
}
html_logo = "_static/images/XRTpy_logo.png"
html_sidebars = {
    # Sidebar removal
    "about_xrt*": [],
    "install*": [],
    "getting_started*": [],
    "bibliography*": [],
    "glossary*": [],
    "feedback_communication*": [],
    "contributing*": [],
    "code_of_conduct*": [],
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

# Render inheritance diagrams in SVG
graphviz_output_format = "svg"

graphviz_dot_args = [
    "-Nfontsize=10",
    "-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Efontsize=10",
    "-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Gfontsize=10",
    "-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
]

# -- Sphinx Gallery ------------------------------------------------------------

sphinx_gallery_conf = {
    "backreferences_dir": Path("generated") / "modules",
    "filename_pattern": "^((?!skip_).)*$",
    "examples_dirs": Path("..") / "examples",
    "within_subsection_order": "ExampleTitleSortKey",
    "gallery_dirs": Path("generated") / "gallery",
    "matplotlib_animations": True,
    "default_thumb_file": "https://raw.githubusercontent.com/HinodeXRT/xrtpy/main/docs/_static/images/XRTpy_logo.png",
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("xrtpy"),
    "only_warn_on_example_error": True,
}

# -- Options for sphinx-copybutton ---------------------------------------------

# Python Repl + continuation, Bash, ipython and qtconsole + continuation, jupyter-console + continuation
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
