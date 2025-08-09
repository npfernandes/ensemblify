"""
Reweighting - ``ensemblify.reweighting``
========================================

:Author(s): Nuno P. Fernandes
:Year: 2024
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains functions for reweighting a generated ensemble of protein conformations using
experimental data, outputting an interactive graphical dashboard.

Example applications
--------------------

Reweight an ensemble using experimental SAXS data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``ensemblify.generation.reweight_ensemble`` function can be used to reweight a conformational
ensemble using an experimental SAXS data file.

For example, we can reweight the ensemble of Histatin5 (Hst5), an intrinsically disordered
peptide with 24 aminoacid residues, using experimental SAXS data of Hst5.

Assuming the path to the required trajectory, topology and experimental data files are assigned
to the variables HST5_TRAJ_PATH, HST5_TOP_PATH and HST5_EXP_SAXS_DATA, respectively, and your
desired trajectory ID is 'Hst5', you should run:

>>> import ensemblify.reweighting as er
>>> er.reweight_ensemble(HST5_TRAJ_PATH, HST5_TOP_PATH, 'Hst5', HST5_EXP_SAXS_DATA)

Available Functions
-------------------

``reweight_ensemble``
    Apply Bayesian Maximum Entropy (BME) reweighting to a conformational ensemble, given
    experimental SAXS data.

"""

from ensemblify.reweighting.ensemble import reweight_ensemble

__all__ = ['reweight_ensemble']
