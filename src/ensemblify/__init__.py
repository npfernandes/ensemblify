"""
Ensemblify - a Python library for generating and analyzing ensembles of protein structures.
===========================================================================================
:Author(s): Nuno P. Fernandes
:Year: 2024
:Copyright: GNU Public License v3

Main Features
-------------

- Generate protein conformational ensembles by changing flexible regions.
- Convert generated ensembles to trajectory file format.
- Calculate theoretical SAXS curves from generated ensembles.
- Analyze generated ensemble's structural properties through interactive plots.
- Reweight generated ensembles using experimental data.
- Re-analyze structural properties of generated ensembles using weights calculated from 
  experimental data, compare to non-reweighted structural properties.

How to access documentation
----------------------------

Documentation is available in two forms:

- Docstrings provided with the code
- An online ReadtheDocs, available at https://ensemblify.readthedocs.io/en/latest/index.html.

Use the built-in `help` function to view a function or module's docstring:

>>> import ensemblify as ey
>>> help(ey)
>>> help(ey.generation)
>>> help(ey.generate_ensemble)

Available subpackages
---------------------

`generation`
    Generate an ensemble of structures.
`conversion`
    Convert ensembles to trajectory files and from these calculate SAXS curves.
`analysis`
    Calculate and plot data describing your ensemble.
`reweighting`
    Reweight a generated ensemble using experimental data.
`utils`
    Auxilliary functions used by other modules.

Utilities
---------

`test()`
    Run Ensemblify unit tests.
`show_config()`
    View Ensemblify's current general configuration.
`update_config()`
    Update Ensemblify's current general configuration.
`clash_checking.check_steric_clashes()`
    Check an already generated ensemble for steric clashes, reporting any found.
`pipeline.ensemblify_pipeline()`
    Function to use all of Ensemblify's functionalities sequentially.

Citation
--------

When using Ensemblify in published work, please cite

    PUB

"""
from ensemblify.config import show_config, update_config
from ensemblify.utils import df_from_pdb, df_to_pdb, extract_pdb_info

__version__ = '1.0.0'

__all__ = ['df_from_pdb',
           'df_to_pdb',
           'extract_pdb_info',
           'show_config',
           'update_config']
