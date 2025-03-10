ensemblify.conversion.trajectory2saxs
=====================================

.. py:module:: ensemblify.conversion.trajectory2saxs

.. autoapi-nested-parse::

   Calculate a set of SAXS curves (.dat) from a trajectory file (.xtc).



Functions
---------

.. autoapisummary::

   ensemblify.conversion.trajectory2saxs.traj2saxs


Module Contents
---------------

.. py:function:: traj2saxs(trajectory: str, topology: str, trajectory_id: str, exp_saxs_file: str) -> str

   Calculate a set of theoretical SAXS curves from a trajectory file using PEPSI-SAXS.

   Calculation is done in chunks distributed across available processor cores.
   A Universe object is created with the given trajectory and topology, which
   allows for writing of temporary individual .pdb files from trajectory frames
   that are then used for SAXS curve calculation at every frame.

   :param trajectory: Path to trajectory file used in Universe object creation.
   :type trajectory: :py:class:`str`
   :param topology: Path to topology file used in Universe object creation.
   :type topology: :py:class:`str`
   :param trajectory_id: Prefix identifier for created files.
   :type trajectory_id: :py:class:`str`
   :param exp_saxs_file: Path to the experimental SAXS data for this protein, used in PEPSI-SAXS
                         for SAXS curve calculation.
   :type exp_saxs_file: :py:class:`str`

   :returns:     Path to the file containing the set of calculated SAXS curves, one for every frame of
                 the trajectory.
   :rtype: str

   Adapted from:
       https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb


