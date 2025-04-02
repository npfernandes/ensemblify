"""
Generation - ``ensemblify.generation``
======================================

:Author(s): Nuno P. Fernandes
:Year: 2024
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains functions for generating an ensemble of protein conformations from an input
structure using an input parameters file.

Example applications
--------------------

Generate an ensemble
~~~~~~~~~~~~~~~~~~~~

The ``ensemblify.generation.generate_ensemble`` function can be used to generate a conformational
ensemble from a .yaml parameters file.

For example, we can generate an ensemble for Histatin5, an intrinsically disordered protein (IDP)
with 24 aminoacid residues.

Assuming the path to the required input parameters file is assigned to a variable named
HST5_PARAMS, you should run:

>>> import ensemblify.generation as eg
>>> eg.generate_ensemble(HST5_PARAMS)

Available Functions
-------------------

``generate_ensemble``
    Generate an ensemble of conformations given a parameters file.

"""

from ensemblify.generation.ensemble import generate_ensemble

__all__ = ['generate_ensemble']
