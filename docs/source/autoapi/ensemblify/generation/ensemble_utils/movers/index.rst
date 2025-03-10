ensemblify.generation.ensemble_utils.movers
===========================================

.. py:module:: ensemblify.generation.ensemble_utils.movers

.. autoapi-nested-parse::

   Custom Mover classes created from the PyRosetta Mover class.



Classes
-------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.movers.SetRandomDihedralsMover


Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.movers.setup_mover


Module Contents
---------------

.. py:class:: SetRandomDihedralsMover(databases: dict[str, dict[str, pandas.DataFrame]], variance: float, log_file: str)



   Custom PyRosetta Mover object that sets random dihedral angles in target residues, taken from
   a given database.

   Inherits from pyrosetta.rosetta.protocols.moves.Mover.

   .. attribute:: databases

      All the available databases to sample from. Mapping of database_ids to
      databases nested dicts, that map residue 1lettercodes to dihedral
      angle values dataframes.

      :type: :py:class:`dict`

   .. attribute:: variance

      New dihedral angle values inserted into sampling regions are sampled from a Gaussian
      distribution centred on the value found in database and percentage variance equal to
      this value.

      :type: :py:class:`float`

   .. attribute:: log_file

      Path to .log file for warnings or error messages related to sampling.

      :type: :py:class:`str`


   .. py:method:: get_name()

      Return the name of this mover.



   .. py:method:: apply(pose: pyrosetta.rosetta.core.pose.Pose, target_resnum: int, database_id: str, secondary_structure: str | None, sampling_mode: str)

      Apply the mover to a Pose, on the given residue number (PyRosetta Pose numbering).

      :param pose: Pose on which the mover will be applied.
      :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
      :param target_resnum: Residue number where to apply the mover (PyRosetta Pose numbering).
      :type target_resnum: :py:class:`int`
      :param database_id: Identifier for which database to sample from.
      :type database_id: :py:class:`str`
      :param secondary_structure: Which secondary structure element to force in this target region.
                                  If None, sample database without restraint.
      :type secondary_structure: :py:class:`str, optional`
      :param sampling_mode: Whether to sample the database considering neighbouring residues ('TRIPEPTIDE')
                            or not ('SINGLERESIDUE').
      :type sampling_mode: :py:class:`str`



.. py:function:: setup_mover(mover_id: str, databases: dict[str, dict[str, pandas.DataFrame]], variance: float, log_file: str) -> pyrosetta.rosetta.protocols.moves.Mover

   Create custom PyRosetta Mover.

   Setup Mover object given a database to sample from and a mover id.

   :param mover_id: Identifier for which CustomMover to create.
   :type mover_id: :py:class:`str`
   :param databases: Mapping of database_ids to databases nested dicts, that map residue 1lettercodes
                     to dihedral angle values dataframes.
   :type databases: :py:class:`dict[str,dict[str,pd.DataFrame]]`
   :param variance: New dihedral angle values inserted into sampling regions are sampled from a Gaussian
                    distribution centred on the value found in database and percentage variance equal to
                    this value.
   :type variance: :py:class:`float`
   :param log_file: Path to .log file for warnings or error messages related to sampling.
   :type log_file: :py:class:`str`

   :returns:     Custom PyRosetta Mover object.
   :rtype: pyrosetta.rosetta.protocols.moves.Mover


