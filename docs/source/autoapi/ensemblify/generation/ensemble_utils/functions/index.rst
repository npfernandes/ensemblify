ensemblify.generation.ensemble_utils.functions
==============================================

.. py:module:: ensemblify.generation.ensemble_utils.functions

.. autoapi-nested-parse::

   Auxiliary functions for sampling.



Functions
---------

.. autoapisummary::

   ensemblify.generation.ensemble_utils.functions.add_intrachain_constraints
   ensemblify.generation.ensemble_utils.functions.add_contacts_constraints
   ensemblify.generation.ensemble_utils.functions.get_targets_from_plddt
   ensemblify.generation.ensemble_utils.functions.setup_pose
   ensemblify.generation.ensemble_utils.functions.setup_minmover
   ensemblify.generation.ensemble_utils.functions.derive_constraint_targets
   ensemblify.generation.ensemble_utils.functions.apply_pae_constraints
   ensemblify.generation.ensemble_utils.functions.apply_constraints
   ensemblify.generation.ensemble_utils.functions.setup_fold_tree


Module Contents
---------------

.. py:function:: add_intrachain_constraints(pose: pyrosetta.rosetta.core.pose.Pose, constraint_targets: tuple[tuple[int, int], Ellipsis], constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet, constraint_function: pyrosetta.rosetta.core.scoring.func.HarmonicFunc | pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc, stdev: float, tolerance: float | None = None) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet

   Add constraints between non-sampled residues in a chain to a constraint set.

   Create all AtomPairConstraints between residues of a Pose object
   present in constraint_targets and add them to a ConstraintSet.

   :param pose: Target Pose object for constraints.
   :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
   :param constraint_targets: Residues between which AtomPairConstraints will be applied.
   :type constraint_targets: :py:class:`tuple[tuple[int,int],...]`
   :param constraint_set: Set of constraints to later be applied to Pose.
   :type constraint_set: :py:class:`pyrosetta.rosetta.core.scoring.constraints.ConstraintSet`
   :param constraint_function: Function to use for each added constraint.
   :type constraint_function: :py:class:`pyrosetta.rosetta.core.scoring.func.HarmonicFunc | pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc`
   :param stdev: Standard deviation value to use in constraints.
   :type stdev: :py:class:`float`
   :param tolerance: Tolerance value to use in constraints (if applicable). Defaults to None.
   :type tolerance: :py:class:`float, optional`

   :returns:     Updated constraint set, with intrachain constraints.
   :rtype: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet


.. py:function:: add_contacts_constraints(pose: pyrosetta.rosetta.core.pose.Pose, contacts: tuple[tuple[tuple[str, tuple[int, int]], tuple[str, tuple[int, int]]], Ellipsis] | None, constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet, constraint_function: pyrosetta.rosetta.core.scoring.func.HarmonicFunc | pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc, stdev: float, tolerance: float | None = None) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet

   Add constraints between residues of different chains/domains that must remain in contact.

   Create all contact constraints (dimerization sites, intrachain folding) targeting residues
   of a Pose object and add them to a ConstraintSet.

   :param pose: Target Pose object for constraints.
   :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
   :param contacts: Residue ranges where two regions are interacting.
   :type contacts: :py:class:`tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...]`
   :param constraint_set: Set of constraints to later be applied to Pose.
   :type constraint_set: :py:class:`pyrosetta.rosetta.core.scoring.constraints.ConstraintSet`
   :param constraint_function: Function to use for each added constraint.
   :type constraint_function: :py:class:`pyrosetta.rosetta.core.scoring.func.HarmonicFunc | pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc`
   :param stdev: Standard deviation value to use in constraints.
   :type stdev: :py:class:`float`
   :param tolerance: Tolerance value to use in constraints (if applicable). Defaults to None.
   :type tolerance: :py:class:`float, optional`

   :returns:     Updated constraint set, with contact constraints.
   :rtype: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet


.. py:function:: get_targets_from_plddt(parameters: dict) -> dict[str, list[int]]

   Get, for each chain, lists of residues with pLDDT value below the threshold.

   The input structure defined in the parameters dictionary must be an AlphaFold model,
   i.e. have the pLDDT value for each residue in the .pdb B-Factor column.

   :param parameters: Dictionary following Ensemblify parameters template.
   :type parameters: :py:class:`dict`

   :returns:     Mapping of each chain to the residue numbers contained in it pertaining
                 to sampled residues with pLDDT below the threshold. For example:

                 {'A': [[234,235,236,237],[536,537,538,539]], 'B': [[124,125,126,127,128,129]] },

                 when the contiguous_res parameter is equal to 4 residues.
   :rtype: dict[str,list[int]]


.. py:function:: setup_pose(input_structure: str) -> pyrosetta.rosetta.core.pose.Pose

   Initialize a Pose object from a sequence, a .txt file containing the sequence or a PDB file.

   The created Pose object is then changed to 'centroid' configuration.

   :param input_structure: Filepath to the input .pdb structure, .txt with sequence or the actual sequence string.
   :type input_structure: :py:class:`str`

   :returns:     Our initial Pose for sampling.
   :rtype: pyrosetta.rosetta.core.pose.Pose


.. py:function:: setup_minmover(scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction, min_id: str, tolerance: float, max_iters: int, dofs: tuple[str, str] = ('bb', 'chi')) -> pyrosetta.rosetta.protocols.minimization_packing.MinMover

   Setup the MoveMap and MinMover for last minimization steps in the sampling process.

   :param scorefxn: Score function used during sampling to evaluate our Pose conformations.
   :type scorefxn: :py:class:`pyrosetta.rosetta.core.scoring.ScoreFunction`
   :param min_id: Identifier for the PyRosetta minimization algorithm.
   :type min_id: :py:class:`str`
   :param tolerance: Value for the MinMover tolerance.
   :type tolerance: :py:class:`float`
   :param max_iters: Maximum iterations of the MinMover.
   :type max_iters: :py:class:`int`
   :param dofs: Defines what angles to set as flexible during minimization.
                Defaults to backbone and sidechain, i.e. ('bb','chi').
   :type dofs: :py:class:`tuple[str,str], optional`

   :returns:     PyRosetta MinMover for last minimization steps in the sampling process.
   :rtype: pyrosetta.rosetta.protocols.minimization_packing.MinMover


.. py:function:: derive_constraint_targets(pose: pyrosetta.rosetta.core.pose.Pose, sampling_targets: dict[str, tuple[tuple[str, tuple[int, Ellipsis], str, str], Ellipsis]]) -> tuple[tuple[int, int], Ellipsis]

   Derive the list of residues to keep constrained based on sampling targets.

   Given a Pose and the target residue ranges for sampling, mark all non-sampled
   residues as constraint targets.
   In the case of a multichain input structure, assumes chains are properly labeled.

   :param pose: Initial Pose object for sampling.
   :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
   :param sampling_targets: Dictionary detailing the target regions for sampling in each chain.
   :type sampling_targets: :py:class:`dict[str,tuple[tuple[str,tuple[int,...],str,str],...]`

   :returns:     All the residue number pairs representing regions on which to apply constraints.
   :rtype: tuple[tuple[int,int],...]


.. py:function:: apply_pae_constraints(pose: pyrosetta.rosetta.core.pose.Pose, pae_filepath: str, plddt_targets: dict, cutoff: float = 10.0, flatten_cutoff: float = 10.0, flatten_value: float = 10.0, weight: float = 1.0, tolerance: float | None = None, adjacency_threshold: int = 8, plddt_scaling_factor: float = 1.0)

   Apply constraints to a Pose created from an AlphaFold structure
   based on the Predicted Aligned Error (PAE) matrix.

   The strength of the constraints scales with the value of the predicted aligned error for that
   residue pair.
   After scaling with PAE value, applied constraints are made weaker when between a residue with
   high pLDDT and a residue with low pLDDT.
   This function should be used instead of the apply_constraints function, they should not be used
   in tandem.

   Adapted from:
       https://github.com/matteoferla/pyrosetta-help

   :param pose: The Pose constraints will be applied to.
   :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
   :param pae_filepath: Path to the pae matrix .json file.
   :type pae_filepath: :py:class:`str`
   :param plddt_targets: Mapping of each chain to the residue numbers that
                         will be sampled (pdb numbering), all with plddt below a threshold.
   :type plddt_targets: :py:class:`dict`
   :param cutoff: Only consider PAE values below this number (low error).
   :type cutoff: :py:class:`float`
   :param flatten_cutoff: Any PAE values below this value will be changed to match flatten value.
   :type flatten_cutoff: :py:class:`float`
   :param flatten_value: Any PAE values below flatten cutoff will be changed to match this value.
   :type flatten_value: :py:class:`float`
   :param weight: Along with the error value, determines the strength of the applied AtomPairConstraints
   :type weight: :py:class:`float`
   :param tolerance: Defines the tolerance of the FlatHarmonicFunction of AtomPairConstraints
                     created from the PAE matrix. Defaults to None.
   :type tolerance: :py:class:`float, optional`
   :param adjacency_threshold: How far away two residues need to be to consider their PAE value. Neighbours are skipped
                               as PAE is best used for determining between domain or between chain confidence.
   :type adjacency_threshold: :py:class:`int`
   :param plddt_scaling_factor: Any constraints setup between residues where one of them has a low pLDDT and another a
                                high pLDDT will be scaled by multiplying its weight by this factor. The higher this
                                value the weaker those constraints will be.
   :type plddt_scaling_factor: :py:class:`float`


.. py:function:: apply_constraints(pose: pyrosetta.rosetta.core.pose.Pose, cst_targets: tuple[tuple[int, int], Ellipsis], contacts: tuple[tuple[tuple[str, tuple[int, int]], tuple[str, tuple[int, int]]], Ellipsis] | None, stdev: float = 10.0, tolerance: float | None = 0.001)

   Apply constraints to non-sampled regions of a Pose object.

   Apply all appropriate intra-region constraints based on cst_targets and all interregion
   constraints based on given restraints.
   not be sampled.
   This function is incompatible with the apply_pae_constraints function, and they should not be
   used in tandem.

   :param pose: Target Pose object for constraints.
   :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
   :param cst_targets: Residues between which AtomPairConstraints will be applied.
   :type cst_targets: :py:class:`tuple[tuple[int,int],...]`
   :param contacts: Residue ranges where two chains are interacting.
   :type contacts: :py:class:`tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...]`
   :param stdev: Standard deviation value to use in constraints.
   :type stdev: :py:class:`float`
   :param tolerance: Tolerance value to use in constraints (if applicable). Defaults to None.
   :type tolerance: :py:class:`float, optional`


.. py:function:: setup_fold_tree(pose: pyrosetta.rosetta.core.pose.Pose, constraint_targets: tuple[tuple[int, int], Ellipsis], contacts: tuple[tuple[tuple[str, tuple[int, int]], tuple[str, tuple[int, int]]], Ellipsis] | None)

   Change a Pose's FoldTree in order to minimize "lever arm" effects during sampling.

   Perform slight alterations to the given Pose's FoldTree to minimize "lever arm" effects
   that might result in movement in constrained regions. These changes are based on which
   residues are going to be constrained. First the most "central" residue of the constrained
   residues in each chain is found, and the FoldTree is changed to a tree that has this
   residue as a parent in that chain: it starts from this residue and goes in both the
   N-terminal and C-terminal direction of the protein chain.

   Resulting FoldTree (2 chains example):

       Chain 1:  1 <----------- "central" res ------------> chain.size()

       Chain 2:  1 <----------- "central" res ------------> chain.size()

       Jump between the two "central" residues.

   If there are inter-chain contacts, after deriving the optimal "central" residues they are
   updated so that every central residue is part of a region that is in contact with another
   chain. This avoids cases where, in multi-chain, proteins, certain folded domains would not
   move relative to other chains which would inadvertedly bias conformational sampling.

   :param pose: Pose object whose FoldTree will be updated.
   :type pose: :py:class:`pyrosetta.rosetta.core.pose.Pose`
   :param constraint_targets: Residues between which AtomPairConstraints will be applied.
   :type constraint_targets: :py:class:`tuple[tuple[int,int],...]`
   :param contacts: Residue ranges where two chains are interacting.
   :type contacts: :py:class:`tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...]`

   Reference:
       See https://docs.rosettacommons.org/demos/latest/tutorials/fold_tree/fold_tree
       for more information about the Rosetta FoldTree.


