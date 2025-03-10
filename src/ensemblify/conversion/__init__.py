"""
Conversion - `ensemblify.conversion`
====================================

:Author(s): Nuno P. Fernandes
:Year: 2024
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains functions for converting an ensemble of protein conformations from a set of
.pdb files into a single compressed .xtc trajectory file and for calculating a theoretical SAXS
curve from a .xtc trajectory file.

Example applications
--------------------

Convert an ensemble to trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `ensemblify.conversion.ensemble2traj` function can be used to create a .xtc trajectory file
from an ensemble of .pdb files.

For example, we can create a trajectory from the set of .pdb structures of Histatin5, an
intrinsically disordered protein (IDP) with 24 aminoacid residues.

Assuming that the path to the directory where Hst5 .pdb structures are stored and the path to
the directory where the created .xtc trajectory file will be stored are assigned to variables
named, respectively, HST5_ENSEMBLE_DIR and HST5_TRAJECTORY_DIR, you should run:

>>> import ensemblify.conversion as ec
>>> ec.ensemble2traj(HST5_ENSEMBLE_DIR,HST5_TRAJECTORY_DIR,'Hst5')

Calculate a theoretical SAXS curve from a trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `ensemblify.conversion.traj2saxs` function can be used to back-calculate an average SAXS
curve from a .xtc trajectory file.

For example, we can calculate a SAXS curve from a .xtc trajectory file of Histatin5, an
intrinsically disordered protein (IDP) with 24 aminoacid residues.

Assuming the path to the required Hst5 .xtc trajectory file is assigned to a variable named
HST5_TRAJECTORY, you should run:

>>> import ensemblify.conversion as ec
>>> ec.traj2saxs(HST5_TRAJECTORY)

Available Functions
-------------------

- `ensemble2traj`

      Create a .xtc trajectory file from an ensemble of .pdb files.

- `traj2saxs`

      Calculate a theoretical SAXS curve from a trajectory file using PEPSI-SAXS.

"""

from ensemblify.conversion.ensemble2trajectory import ensemble2traj
from ensemblify.conversion.trajectory2saxs import traj2saxs

__all__ = ['ensemble2traj','traj2saxs']
