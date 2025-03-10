ensemblify.conversion.ensemble2trajectory
=========================================

.. py:module:: ensemblify.conversion.ensemble2trajectory

.. autoapi-nested-parse::

   Create a trajectory (.xtc) file from an ensemble of .pdb files.



Functions
---------

.. autoapisummary::

   ensemblify.conversion.ensemble2trajectory.ensemble2traj


Module Contents
---------------

.. py:function:: ensemble2traj(ensemble_dir: str | None = os.getcwd(), trajectory_dir: str | None = os.getcwd(), trajectory_id: str | None = '', trajectory_size: int | None = 10000) -> tuple[str, str]

   Create a trajectory (.xtc) file from an ensemble of .pdb files.

   Additionally, one of the .pdb files used to create this trajectory is kept as
   a .pdb topology file.

   :param ensemble_dir: Path to directory where all the .pdb files are stored. Defaults to current working
                        directory.
   :type ensemble_dir: :py:class:`str`
   :param trajectory_dir: Path to directory where trajectory .xtc file will be created. Will be created if it
                          does not exist. Defaults to current working directory.
   :type trajectory_dir: :py:class:`str`
   :param trajectory_id: Prefix identifier for any created files.
   :type trajectory_id: :py:class:`str`
   :param trajectory_size: Number of randomly sampled .pdb files to use for trajectory creation.
   :type trajectory_size: :py:class:`int`

   :returns:

                 trajectory_path (str):
                     Path to created trajectory .xtc file.
                 topology_path (str):
                     Path to created topology .pdb file.
   :rtype: tuple[str,str]


