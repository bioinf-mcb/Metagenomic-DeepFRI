# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import configparser
import datetime
import os
import re

import semantic_version
import sphinx_bootstrap_theme

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

docssrc_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(docssrc_dir)

# -- Project information -----------------------------------------------------
import mDeepFRI  # noqa: E402

project = mDeepFRI.__name__
author = re.match('(.*) <.*>', mDeepFRI.__author__).group(1)
year = datetime.date.today().year
copyright = '{}, {}'.format("2022" if year == 2022 else "2020-{}".format(year),
                            author)

# extract the semantic version
semver = semantic_version.Version.coerce(mDeepFRI.__version__)
version = str(semver.truncate(level="patch"))
release = str(semver)

# extract the project URLs from ``setup.cfg``
cfgparser = configparser.ConfigParser()
cfgparser.read(os.path.join(project_dir, "setup.cfg"))
project_urls = dict(
    map(str.strip, line.split(" = ", 1))
    for line in cfgparser.get("metadata", "project_urls").splitlines()
    if line.strip())

# -- Sphinx Setup ------------------------------------------------------------


def setup(app):
    # Add custom stylesheet
    app.add_css_file("css/main.css")


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    # "sphinx.ext.imgconverter",
    "sphinx.ext.napoleon",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.extlinks",
    "sphinxcontrib.jquery",
    "sphinx_bootstrap_theme",
    "nbsphinx",
    "recommonmark",
    "sphinx_click.ext",
    "IPython.sphinxext.ipython_console_highlighting",
]

bibtex_bibfiles = ['_static/bibtex/references.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['build']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# https://bootswatch.com/
html_theme = "bootstrap"

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

html_theme_options = {
    # Bootswatch (http://bootswatch.com/) theme.
    "bootswatch_theme":
    "flatly",
    # Choose Bootstrap version.
    "bootstrap_version":
    "3",
    # Tab name for entire site. (Default: "Site")
    "navbar_site_name":
    "Documentation",
    # HTML navbar class (Default: "navbar") to attach to <div> element.
    # For black navbar, do "navbar navbar-inverse"
    "navbar_class":
    "navbar",
    # Render the next and previous page links in navbar. (Default: true)
    "navbar_sidebarrel":
    True,
    # Render the current pages TOC in the navbar. (Default: true)
    "navbar_pagenav":
    False,
    # A list of tuples containing pages or urls to link to.
    "navbar_links":
    [("GitHub", cfgparser.get("metadata", "url").strip(), True)] +
    [(k, v, True) for k, v in project_urls.items() if k in {"Zenodo", "PyPI"}],
    "admonition_use_panel":
    True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ['.rst', '.md']

# The master toctree document.
master_doc = 'index'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

html_sidebars = {
    "*": ["localtoc.html"],
    "api/*": ["localtoc.html"],
}

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = mDeepFRI.__name__

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).

# -- Extension configuration -------------------------------------------------

# -- Options for imgmath extension -------------------------------------------

imgmath_image_format = "svg"

# -- Options for napoleon extension ------------------------------------------

napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True
napoleon_include_private_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_rtype = False

# -- Options for autodoc extension -------------------------------------------

autoclass_content = "class"
autodoc_member_order = 'groupwise'
autosummary_generate = []

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "biopython": ("https://biopython.org/docs/latest/api/", None),
}

# -- Options for recommonmark extension --------------------------------------

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

# -- Options for nbsphinx extension ------------------------------------------

nbsphinx_execute = 'auto'
nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 96}",
]
