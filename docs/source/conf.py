# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../src/ensemblify'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ensemblify'
copyright = '2025, Nuno P. Fernandes'
author = 'Nuno P. Fernandes'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "autoapi.extension",
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo' # 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Options for autodoc -----------------------------------------------------

autodoc_mock_imports = ["pyrosetta"] # ,"ensemblify", "MDAnalysis", "mdtraj", "numpy", "pandas", "pyarrow", "scikit-learn", "scipy", "tqdm", "biopython", "plotly", "pyyaml", "ray",
autodoc_typehints = 'description'

# -- Options for napoleon -----------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_type_aliases = None
napoleon_attr_annotations = True

# -- Options for autoapi -----------------------------------------------------

autoapi_dirs = ['../../src/ensemblify']
autoapi_type = "python"
# autoapi_template_dir = "_templates/autoapi"
autoapi_options = [
    "members",
    "show-module-summary",
]
autoapi_ignore = ['*third_party*','*cli*']

autoapi_keep_files = True
autodoc_typehints = "signature"
suppress_warnings = ['autoapi.python_import_resolution']