ensemblify.generation.ensemble_utils.sampling_utils
===================================================

.. py:module:: ensemblify.generation.ensemble_utils.sampling_utils

.. autoapi-nested-parse::

   Auxiliary functions for setting up the conformational sampling process.



Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.sampling_utils.setup_sampling_logging
   ensemblify.generation.ensemble_utils.sampling_utils.setup_sampling_parameters
   ensemblify.generation.ensemble_utils.sampling_utils.setup_sampling_initial_pose
   ensemblify.generation.ensemble_utils.sampling_utils.sample_pdb


Module Contents
---------------

.. py:function:: setup_sampling_logging(sampling_log: str) -> tuple[logging.Logger, str, str]

   Setup logging handlers and files for PyRosetta sampling and Ray.

   :param sampling_log: Path to sampling .log file.
   :type sampling_log: :py:class:`str`

   :returns:

                 logger (logging.Logger):
                     The Logger object associated with the sampling .log file
                 ray_log (str):
                     Filepath to .log file with Ray log messages.
                 pyrosetta_log (str):
                     Filepath to .log file with PyRosetta log messages.
   :rtype: tuple[logging.Logger,str,str]


.. py:function:: setup_sampling_parameters(parameters_file: str) -> dict

   Update the parameters dictionary before sampling.

   In the 'targets' parameter, change a target from e.g. [1,54] to (range(1,55)).
   If using an AlphaFold model as a starting structure, keep in the sampling ranges only regions
   of at least a certain contiguous size where each residue's pLDTT is below a threshold.
   Targets, secondary structure biases and contacts are also updated to tuples instead of lists.

   :param parameters_file: Path to parameters file following the Ensemblify template.
   :type parameters_file: :py:class:`str`

   :returns:     The updated params dictionary.
   :rtype: dict


.. py:function:: setup_sampling_initial_pose(params: dict, sampling_log: str) -> pyrosetta.rosetta.core.pose.Pose

   Create initial Pose for sampling and apply any required constraints.

   Constraints applied to the pose are stored in a constraints.cst file stored in the same
   directory as the sampling_log file.

   :param params: Parameters following the Ensemblify template.
   :type params: :py:class:`dict`
   :param sampling_log: Path to the sampling .log file.
   :type sampling_log: :py:class:`str`

   :returns:     Object to be used as the starting structure for the sampling process.
   :rtype: pyrosetta.rosetta.core.pose.Pose


.. py:function:: sample_pdb(ppose: pyrosetta.distributed.packed_pose.core.PackedPose, databases: dict[str, dict[str, pandas.DataFrame]], targets: dict[str, tuple[tuple[str, tuple[int, Ellipsis], str, str]]], output_path: str, job_name: str, decoy_num: str = '', log_file: str = os.path.join(os.getcwd(), 'pyrosetta.log'), ss_bias: tuple[tuple[tuple[str, tuple[int, int], str], Ellipsis], int] | None = None, variance: float | None = 0.1, sampler_params: dict[str, dict[str, int]] = {'MC': {'temperature': 200, 'max_loops': 200}}, scorefxn_id: str = 'score0', scorefxn_weight: float = 1.0, minimizer_id: str = 'dfpmin_armijo_nonmonotone', minimizer_tolerance: float = 0.001, minimizer_maxiters: int = 5000, minimizer_finalcycles: int = 5, cst_weight: int = 1, cstviolation_threshold: float = 0.015, cstviolation_maxres: int = 20) -> str | None

   Sample dihedral angles from a database into target regions of a given structure.

   :param ppose: Reference to the initial structure.
   :type ppose: :py:class:`pyrosetta.distributed.packed_pose.core.PackedPose`
   :param databases: Reference to the databases dictionary.
   :type databases: :py:class:`dict[str,dict[str,pd.DataFrame]]`
   :param targets: Dictionary detailing the target regions for sampling in each chain.
   :type targets: :py:class:`dict[str,tuple[tuple[str,tuple[int,...],str,str]]]`
   :param output_path: Path to directory where sampled structures will be written to.
   :type output_path: :py:class:`str`
   :param job_name: Prefix identifier for generated structures.
   :type job_name: :py:class:`str`
   :param decoy_num: Identifier to differentiate between different decoys of the same batch in a
                     multiprocessing context. Defaults to ''.
   :type decoy_num: :py:class:`str`
   :param log_file: Path to the PyRosetta .log file. Defaults to 'pyrosetta.log' in current working
                    directory.
   :type log_file: :py:class:`str`
   :param ss_bias: Secondary Structure Bias with the desired percentage of total structures to respect
                   this bias. Defaults to None.
   :type ss_bias: :py:class:`tuple[tuple[tuple[str,tuple[int,int],str],...],int] | None`
   :param variance: New dihedral angle values inserted into sampling regions are sampled from a Gaussian
                    distribution centred on the value found in database and percentage variance equal to
                    this value. Defaults to 0.10 (10%).
   :type variance: :py:class:`float`
   :param sampler_params: Parameters for the used sampler, assumes MonteCarloSampler is used. Defaults to
                          {'MC':{'temperature':200,'max_loops':200}}.
   :type sampler_params: :py:class:`dict[str,dict[str,int]]`
   :param scorefxn_id: PyRosetta ScoreFunction identifier. Must pertain to a .wst weights file present in
                       /.../pyrosetta/database/scoring/weights/ . Defaults to 'score0'.
   :type scorefxn_id: :py:class:`str`
   :param scorefxn_weight: Weight for the repulsive Van der Waals term in the ScoreFunction. Will only have an
                           effect if the ScoreFunction has a repulsive Van der Waals term. Defaults to 1.0.
   :type scorefxn_weight: :py:class:`float`
   :param minimizer_id: PyRosetta minimization algorithm identifier used in MinMover.
                        Defaults to 'dfpmin_armijo_nonmonotone'.
   :type minimizer_id: :py:class:`str`
   :param minimizer_tolerance: Tolerance value for the PyRosetta MinMover object. Defaults to 0.001.
   :type minimizer_tolerance: :py:class:`float`
   :param minimizer_maxiters: Maximum iterations value for the PyRosetta MinMover object. Defaults to 5000.
   :type minimizer_maxiters: :py:class:`int`
   :param minimizer_finalcycles: Number of times to apply the MinMover to our final structure. Defaults to 5.
   :type minimizer_finalcycles: :py:class:`int`
   :param cst_weight: Weight of the AtomPairConstraint term in the ScoreFunction. Defaults to 1.
   :type cst_weight: :py:class:`int`
   :param cstviolation_threshold: Any residue with AtomPairConstraint score term value above this threshold is considered
                                  in violation of the applied constraints. Defaults to 0.015.
   :type cstviolation_threshold: :py:class:`float`
   :param cstviolation_maxres: Number of residues allowed to be above the constraint violation threshold.
                               Defaults to 20.
   :type cstviolation_maxres: :py:class:`int`

   :returns:     Path to the sampled .pdb structure. Only written and returned if the
                 sampled structure is valid (does not violate the applied constraints).
                 Otherwise, return None.
   :rtype: str | None


