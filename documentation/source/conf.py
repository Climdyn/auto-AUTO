# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'auto-AUTO'
copyright = '2025, Jonathan Demaeyer and Oisín Hamilton'
author = 'Jonathan Demaeyer and Oisin Hamilton'
release = 'v0.5.2'
version = release

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    #    'nbsphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosectionlabel',
    # 'sphinxcontrib.bibtex'
]

templates_path = ['_templates']
exclude_patterns = ['files/notebooks/*', '**.ipynb_checkpoints']

# -- Debugging ---------------------------------------------------------------

# nitpicky = True
#
# nitpick_ignore = [('py:class', 'optional'),
#                   ('py:class', 'iterable'),
#                   ('py:class', '2-tuple'),
#                   ('py:class', 'Sympy expression'),
#                   ('py:class', 'callable')]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Options for LaTeX output -------------------------------------------------

latex_elements = {
    'printindex': r'\def\twocolumn[#1]{#1}\printindex',
}


# -- Extension configuration -------------------------------------------------
# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       'numpy': ('https://numpy.org/doc/stable/', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/', None),
                       'matplotlib': ('https://matplotlib.org/stable/', None),
                       }

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# mathjax config

mathjax_path = r"https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/latest.js?config=TeX-MML-AM_CHTML"

# autosection label option

autosectionlabel_prefix_document = True

# bibtex config

# bibtex_bibfiles = ['files/model/ref.bib']

# -- Global links definition -----------------

rst_epilog = """
.. role:: raw-html(raw)
   :format: html

.. |AUTO| replace:: :raw-html:`<a href="https://github.com/auto-07p/auto-07p">AUTO</a>`
.. |Matplotlib| replace:: :raw-html:`<a href="https://matplotlib.org/">Matplotlib</a>`

"""
