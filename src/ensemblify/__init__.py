"""
Ensemblify - a Python library for generating and analyzing ensembles of protein structures.
=====

Main Features
-------------
    - Generate conformational ensembles of flexible regions of protein structures.
    - Convert the generated ensembles to trajectory file format.
    - Calculate theoretical SAXS curves from generated ensembles.
    - Analyze generated ensemble's structural properties through interactive plots.
    - Reweigh generated ensembles using experimental data.
    - Re-analyze structural properties of generated ensembles using weights calculated
    from experimental data, compare to non-reweighted structural properties.

How to use the documentation
----------------------------
Documentation is available in two forms: docstrings provided with the code, and an online
ReadtheDocs, available from `<https://ensemblifyreadthedocs.org>`.

The docstring examples assume that `ensemblify` has been imported as ``ey``::

  >>> import ensemblify as ey

Code snippets are indicated by three greater-than signs::

  >>> x = 42
  >>> x = x + 1

Use the built-in ``help`` function to view a function's docstring::

  >>> help(ey.generate_ensemble)

Available subpackages
---------------------
generation
    Generate an ensemble of structures.
conversion
    Convert ensembles to trajectory files and from these calculate SAXS curves.
analysis
    Calculate and plot data describing your ensemble.
reweighting
    Reweigh a generated ensemble using experimental data.
utils
    Generally useful auxilliary functions for the other modules.

Utilities
---------
test
    Run ensemblify unittests.
show_config
    View Ensemblify's current general configuration.
update_config
    Update Ensemblify's current general configuration.
pipeline
    Function to use all of Ensemblify's functionalities sequentially.
"""
from ensemblify.analysis import analyze_trajectory
from ensemblify.conversion import ensemble2traj,traj2saxs
from ensemblify.generation import generate_ensemble
from ensemblify.reweighting import reweigh_ensemble
from ensemblify.pipeline import ensemblify_pipeline
from ensemblify.config import show_config,update_config
from ensemblify.utils import df_from_pdb, df_to_pdb, extract_pdb_info, kde

__all__ = ['analyze_trajectory',
           'ensemble2traj',
           'traj2saxs',
           'generate_ensemble',
           'reweigh_ensemble',
           'ensemblify_pipeline',
           'show_config',
           'update_config',
           'df_from_pdb',
           'df_to_pdb',
           'extract_pdb_info',
           'kde']
