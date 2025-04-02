"""
Reweighting - ``ensemblify.reweighting``
========================================

:Author(s): Nuno P. Fernandes
:Year: 2024
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains functions for reweighting a generated ensemble of protein conformations using
experimental SAXS data, outputting an interactive graphical dashboard.

Example applications
--------------------

Reweight an ensemble
~~~~~~~~~~~~~~~~~~~~
The ``ensemblify.generation.reweight_ensemble`` function can be used to reweight a conformational
ensemble using an experimental SAXS data file.

For example, we can reweight the ensemble of Histatin5 (Hst5), an intrinsically disordered
protein (IDP) with 24 aminoacid residues, using experimental SAXS data of Hst5.

Assuming the path to the required Hst5 .xtc trajectory file and path to the required Hst5
.dat experimental SAXS data file are assigned to variables named,
respectively, HST5_TRAJECTORY and HST5_SAXS, you should run:

>>> import ensemblify.reweighting as er
>>> er.reweight_ensemble(HST5_TRAJECTORY,HST5_SAXS)

Available Functions
-------------------

- ``reweight_ensemble``

      Apply Bayesian Maximum Entropy (BME) reweighting to a conformational ensemble, given
      experimental SAXS data.

"""

from ensemblify.reweighting.ensemble import reweight_ensemble

__all__ = ['reweight_ensemble']
