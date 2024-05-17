# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import shutil

# Path
# Based on https://github.com/pyscaffold/pyscaffold
__location__ = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(__location__, "../src"))

# API doc
try:  # for Sphinx >= 1.7
    from sphinx.ext import apidoc
except ImportError:
    from sphinx import apidoc

output_dir = os.path.join(__location__, "api")
module_dir = os.path.join(__location__, "../src/eeepy")
try:
    shutil.rmtree(output_dir)
except FileNotFoundError:
    pass

try:
    import sphinx

    cmd_line = f"sphinx-apidoc --implicit-namespaces -f -o {output_dir} {module_dir} -T"

    args = cmd_line.split(" ")
    if tuple(sphinx.__version__.split(".")) >= ("1", "7"):
        # This is a rudimentary parse_version to avoid external dependencies
        args = args[1:]

    apidoc.main(args)
except Exception as e:
    print(f"Running `sphinx-apidoc` failed!\n{e}")


try:
    from eeepy import __version__ as version
except ImportError:
    version = ""

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'eeepy'
copyright = '2023, ecyht2'
author = 'ecyht2'
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Options for autodoc
autodoc_typehints = 'both'
