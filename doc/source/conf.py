# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import os

from sphinx_gallery.sorting import ExplicitOrder

# -- Project information -----------------------------------------------------

project = 'pfsspy'
copyright = '2018-2022 pfsspy contributors'
author = 'pfsspy contributors'

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx_automodapi.automodapi',
    'sphinx_gallery.gen_gallery',
]
napoleon_google_docstring = False

sphinx_gallery_conf = {
    'ignore_pattern': '.*helpers.py',
    'examples_dirs': '../../examples',
    'gallery_dirs': 'auto_examples',
    'subsection_order': ExplicitOrder(['../../examples/using_pfsspy',
                                       '../../examples/finding_data',
                                       '../../examples/utils',
                                       '../../examples/pfsspy_info',
                                       '../../examples/testing']),
    'reference_url': {'sphinx_gallery': None}
}

# The master toctree document.
master_doc = 'index'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

default_role = 'py:obj'

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'pfsspydoc'

os.environ["JSOC_EMAIL"] = 'jsoc@sunpy.org'

# -- Options for intersphinx extension
# ------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/3/': None,
                       'https://docs.astropy.org/en/stable': None,
                       'https://docs.sunpy.org/en/stable': None,
                       'https://numpy.org/doc/stable': None,
                       'https://streamtracer.readthedocs.io/en/stable': None,
                       'https://docs.scipy.org/doc/scipy/': None,
                       'https://reproject.readthedocs.io/en/stable/': None}

nitpicky = True
