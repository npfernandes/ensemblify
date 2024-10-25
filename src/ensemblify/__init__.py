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
Documentation is available in two forms: docstrings provided with the code, and an online
ReadtheDocs, available from https://ensemblify.readthedocs.io/en/latest/index.html.

Use the built-in `help` function to view a function or module's docstring:

```
>>> import ensemblify as ey
>>> help(ey)
>>> help(ey.generation)
>>> help(ey.generate_ensemble)
```

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
`check_steric_clashes()`
    Check an already generated ensemble for steric clashes, reporting any found.
`pipeline()`
    Function to use all of Ensemblify's functionalities sequentially.

Citation
--------

When using Ensemblify in published work, please cite

    PUB

"""

from ensemblify.analysis import analyze_trajectory
from ensemblify.clash_checking import check_steric_clashes
from ensemblify.config import show_config, update_config
from ensemblify.conversion import ensemble2traj, traj2saxs
from ensemblify.generation import generate_ensemble
from ensemblify.pipeline import ensemblify_pipeline
from ensemblify.reweighting import reweight_ensemble
from ensemblify.utils import df_from_pdb, df_to_pdb, extract_pdb_info

__version__ = '1.0.0'

__all__ = ['analyze_trajectory',
           'ensemble2traj',
           'traj2saxs',
           'generate_ensemble',
           'reweight_ensemble',
           'check_steric_clashes',
           'df_from_pdb',
           'df_to_pdb',
           'extract_pdb_info',
           'ensemblify_pipeline',
           'show_config',
           'update_config']
