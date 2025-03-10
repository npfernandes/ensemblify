ensemblify.generation.ensemble_utils.samplers
=============================================

.. py:module:: ensemblify.generation.ensemble_utils.samplers

.. autoapi-nested-parse::

   Custom Sampler classes to generate protein conformations.



Classes
-------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.samplers.MonteCarloSampler


Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.samplers.setup_samplers


Module Contents
---------------

.. py:class:: MonteCarloSampler(scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction, databases: dict, mover_id: str, smp_params: dict[str, int], variance: float, log_file: str)

   Custom MonteCarlo sampler for sampling dihedral angles.

   .. attribute:: scorefxn

      PyRosetta score function to be used for evaluating Pose objects
      during sampling.

      :type: :py:class:`pyrosetta.rosetta.core.scoring.ScoreFunction`

   .. attribute:: databases

      All the available databases to sample from. Mapping of database_ids to
      databases nested dicts, that map residue 1lettercodes to dihedral
      angle values dataframes.

      :type: :py:class:`dict`

   .. attribute:: mover

      Custom PyRosetta Mover used to apply dihedral angle changes to a Pose.

      :type: :py:class:`pyrosetta.rosetta.protocols.moves.Mover`

   .. attribute:: params

      
      
      Hyperparameters for this sampler (temperature and maximum loops):
          temperature (int):
              A measure of how probable it is to accept Pose objects with a worse score than
              the current one after applying our Mover, according to the acceptance criterion.
          maximum loops (int):
              The maximum amount of attempts without accepting a Move before moving on to
              the next residue to sample.

      :type: :py:class:`dict`

   .. attribute:: log_file

      Path to .log file for warnings or error messages related to sampling.

      :type: :py:class:`str`


   .. py:method:: apply(pose: pyrosetta.rosetta.core.pose.Pose, target: list[int], chain: str, database_id: str, ss_bias: tuple[tuple[str, tuple[int, int], str], Ellipsis] | None, sampling_mode: str)

      Perform MC sampling on the given pose, in the given target residue range.

      :param pose: Pose to be modified during sampling.
      :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
      :param target: Residue range on which sampling will be applied.
      :type target: :py:class:`list[int]`
      :param chain: Letter identifier for the current chain being sampled.
      :type chain: :py:class:`str`
      :param database_id: Identifier for which database to sample from.
      :type database_id: :py:class:`str`
      :param ss_bias: Information about types of secondary structure biases, including which chain and
                      residue numbers they should be applied on.
      :type ss_bias: :py:class:`tuple[tuple[str,tuple[int,int],str],...] | None`
      :param sampling_mode: Whether to sample the database considering neighbouring residues ('TRIPEPTIDE')
                            or not ('SINGLERESIDUE').
      :type sampling_mode: :py:class:`str`



.. py:function:: setup_samplers(sampler_params: dict[str, dict[str, int]], variance: float, scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction, databases: dict[str, dict[str, pandas.DataFrame]], log_file: str) -> dict[str, MonteCarloSampler]

   Create all Sampler objects to be used during sampling.

   Create a dictionary with all the samplers that will be used during sampling,
   given a list of sampler_ids and certain parameters.

   :param sampler_params: Parameters for each sampler to setup.
   :type sampler_params: :py:class:`dict[str,dict[str,int]`
   :param scorefxn: PyRosetta score function, with desired constraints already added.
   :type scorefxn: :py:class:`pyrosetta.rosetta.core.scoring.ScoreFunction`
   :param databases: Mapping of database_ids to databases nested dicts, that map residue 1lettercodes
                     to dihedral angle values dataframes.
   :type databases: :py:class:`dict[str,dict[str,pd.DataFrame]]`

   :returns:     Mapping of sampler_ids to sampler objects to use during sampling.
   :rtype: dict[str,MonteCarloSampler]


