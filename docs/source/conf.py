# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ensemblify'
copyright = '2025, Nuno P. Fernandes'
author = 'Nuno P. Fernandes'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon", # allow for the Google and Numpy docstring formats
    "autoapi.extension", # automatically generate API Reference
    "sphinx_copybutton", # add a copy button to code blocks
    "myst_parser", # support for markdown files
    "sphinx_tabs.tabs", # support for tabs in markdown
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = []
suppress_warnings = ['autoapi.python_import_resolution','myst.xref_missing']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_title = "Ensemblify"
html_logo = "../assets/logo.png"
html_theme = 'furo'
html_theme_options = {
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/npfernandes/ensemblify",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}

# These folders are copied to the documentation's HTML output
html_static_path = ['_static']

# These paths are either relative to html_static_path or fully qualified paths (eg. https://...)
html_css_files = [
    'css/custom.css',
]

# # -- Options for napoleon -----------------------------------------------------
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
autoapi_template_dir = "_templates/autoapi"
autoapi_options = [
    "members",
    "show-module-summary",
]
autoapi_ignore = ['*third_party*','*cli*']
autoapi_add_toctree_entry = False
autoapi_keep_files = False
autodoc_typehints = "description"

# -- Options for sphinx_copybutton -------------------------------------------------
copybutton_prompt_text = r"\$ |\(ensemblify_env\) \$ |>>> "
copybutton_prompt_is_regexp = True
copybutton_copy_empty_lines = False

# -- Options for myst --------------------------------------------------------
myst_enable_extensions = [
    "html_image",
    "colon_fence",
]
myst_heading_anchors = 6

# -- Options for sphinx_tabs ------------------------------------------------
sphinx_tabs_disable_tab_closing = True
sphinx_tabs_disable_css_loading = True
