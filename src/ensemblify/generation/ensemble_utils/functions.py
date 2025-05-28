"""Auxiliary functions for sampling."""

# IMPORTS
## Standard Library Imports
import json
from copy import deepcopy

## Third Party Imports
import numpy as np
import pyrosetta
import pyrosetta.distributed.io as io

## Local Imports
from ensemblify.utils import df_from_pdb

# FUNCTIONS
def add_intrachain_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    constraint_targets: tuple[tuple[int,int],...],
    constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet,
    stdev: float | None = 10.0,
    tolerance: float | None = None,
) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
    """
    Add constraints between desired residues of the same domain to a constraint set.
    
    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Target Pose object for constraints. Only used to extract residue ID and coordinates.
        constraint_targets (tuple[tuple[int,int],...]):
            Residues between which AtomPairConstraints will be added.
        constraint_set (pyrosetta.rosetta.core.scoring.constraints.ConstraintSet):
            Set of constraints to later be applied to Pose.
        stdev (float, optional):
            Standard deviation value to use in constraints. Defaults to 10.0.
        tolerance (float, optional):
            Tolerance value to use in constraints. Defaults to None.
    
    Returns:
        pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
            Updated ConstraintSet object, with added intrachain constraints.
    """
    # Initialize working ConstraintSet object
    working_cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
    working_cs.detached_copy(constraint_set)

    # Iterate over target residue ranges and add constraints
    for target in constraint_targets:
        target_range = list(range(target[0],target[1]+1))
        for i,res_id_1 in enumerate(target_range):
            for res_id_2 in target_range[i+1:]:
                a1 = pyrosetta.rosetta.core.id.AtomID(pose.residue(res_id_1).atom_index('CA'),
                                                      res_id_1)
                a2 = pyrosetta.rosetta.core.id.AtomID(pose.residue(res_id_2).atom_index('CA'),
                                                      res_id_2)
                a1_xyz = pose.residue(a1.rsd()).xyz(a1.atomno())
                a2_xyz = pose.residue(a2.rsd()).xyz(a2.atomno())
                d = (a1_xyz - a2_xyz).norm()

                if tolerance:
                    apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
                        a1,
                        a2,
                        pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc(x0_in=d,
                                                                             sd_in=stdev,
                                                                             tol_in=tolerance)
                    )

                else:
                    apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
                        a1,
                        a2,
                        pyrosetta.rosetta.core.scoring.func.HarmonicFunc(x0_in=d,
                                                                         sd_in=stdev)
                    )

                working_cs.add_constraint(apc)

    return working_cs


def add_contacts_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    contacts: tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...] | None,
    constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet,
    stdev: float | None = 10.0,
    tolerance: float | None = None,
    ) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
    """
    Add constraints between residue regions whose relative position must be conserved to a
    constraint set.

    This would include, for example, dimerization sites or long range intrachain folding.
    The word 'contacts' is used here for convenience, these could be any two regions whose
    relative position should be conserved, even if they are far apart in the Pose.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Target Pose object for constraints.
        contacts (tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...]):
            Residue ranges of regions whose relative position should be conserved.
        constraint_set (pyrosetta.rosetta.core.scoring.constraints.ConstraintSet):
            Set of constraints to later be applied to Pose.
        stdev (float, optional):
            Standard deviation value to use in constraints. Defaults to 10.0.
        tolerance (float, optional):
            Tolerance value to use in constraints (if applicable). Defaults to None.
    
    Returns:
        pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
            Updated ConstraintSet, with added relative position constraints.
    """
    # Initialize working ConstraintSet object
    working_cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
    working_cs.detached_copy(constraint_set)

    if contacts is None:
        return working_cs

    # Iterate over residue ranges and add constraints
    for contact in contacts:
        # contact: ( ('X', (x1,x2) ) , ( 'Y', (y1,y2) ) )

        # Chain X
        chain_x = contact[0][0]

        # Contact residue range for X
        x_start = contact[0][1][0]
        x_end = contact[0][1][1]
        inter_range_x = [ pyrosetta.rosetta.core.pose.pdb_to_pose(pose, res_id, chain_x)
                            for res_id in range(x_start, x_end+1) ]

        # Chain Y
        chain_y = contact[1][0]

        # Contact residue range for Y
        y_start = contact[1][1][0]
        y_end = contact[1][1][1]
        inter_range_y = [ pyrosetta.rosetta.core.pose.pdb_to_pose(pose, res_id, chain_y)
                            for res_id in range(y_start,y_end+1) ]

        # Apply inter-region constraints between region (x1,x2) and (y1,y2)
        for res_id_x in inter_range_x:
            for res_id_y in inter_range_y:
                a1 = pyrosetta.rosetta.core.id.AtomID(pose.residue(res_id_x).atom_index('CA'),
                                                      res_id_x)
                a2 = pyrosetta.rosetta.core.id.AtomID(pose.residue(res_id_y).atom_index('CA'),
                                                      res_id_y)
                a1_xyz = pose.residue(a1.rsd()).xyz(a1.atomno())
                a2_xyz = pose.residue(a2.rsd()).xyz(a2.atomno())
                d = (a1_xyz - a2_xyz).norm()

                if tolerance:
                    apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
                        a1,
                        a2,
                        pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc(x0_in=d,
                                                                            sd_in=stdev,
                                                                            tol_in=tolerance)
                    )

                else:
                    apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
                        a1,
                        a2,
                        pyrosetta.rosetta.core.scoring.func.HarmonicFunc(x0_in=d,
                                                                        sd_in=stdev)
                    )

                working_cs.add_constraint(apc)

    return working_cs


def get_targets_from_plddt(parameters: dict) -> dict[str,list[int]]:
    """Get, for each chain, lists of residues with pLDDT value below the threshold.

    The input structure defined in the parameters dictionary must be an AlphaFold model,
    i.e. have the pLDDT value for each residue in the .pdb B-Factor column.

    Args:
        parameters (dict):
            Dictionary following Ensemblify parameters template.
    
    Returns:
        dict[str,list[int]]:
            Mapping of each chain to the residue numbers contained in it pertaining
            to sampled residues with pLDDT below the threshold. For example:

            {'A': [[234,235,236,237],[536,537,538,539]], 'B': [[124,125,126,127,128,129]] },

            when the contiguous_res parameter is equal to 4 residues.
    """
    # Get unfiltered sampling residues for each chain
    targets = parameters['targets']
    chains_residues = {}
    for chain,regions in targets.items():
        sampling_res = set()
        for region in regions:
            residues = region[1]
            for i in range(residues[0],residues[1]+1):
                sampling_res.add(i)
        chains_residues[chain] = sampling_res

    # Get b-factors (pLDDT) for residues in each chain
    input_af_model = parameters['sequence']
    af_model_df = df_from_pdb(input_af_model)

    # Check if B-Factor column was properly read from file
    assert not af_model_df['B-Factor'].isin(['0.0']).all(), ('B-Factor column read from file is '
                                                             'empty! Make sure your input .pdb '
                                                             'file is properly formatted.')

    chain_resn_bfact = af_model_df[['ChainID','ResidueNumber','B-Factor']]

    # Get low confidence residues (pLDTT < threshold) in each chain
    chains_low_confidence_resn = {}
    for i in range(chain_resn_bfact.shape[0]):
        chain, resn, bfact = chain_resn_bfact.iloc[i]
        try:
            chains_low_confidence_resn[chain]
        except KeyError:
            chains_low_confidence_resn[chain] = []
        if (resn in chains_residues[chain] and
            resn not in chains_low_confidence_resn[chain] and
            bfact < parameters['plddt_params']['threshold']):
            chains_low_confidence_resn[chain].append(resn)

    # Get which residues participate in contiguous regions of at least a certain size
    lcr_sets = {}
    for chain,resrange in chains_low_confidence_resn.items():
        lcr_sets[chain] = set(resrange)

    bfact_sampling_targets = {}
    for chain,lcresn in chains_low_confidence_resn.items():
        bfact_sampling_targets[chain] = []
        curr_streak = []
        for resn in lcresn:
            if not curr_streak:
                curr_streak.append(resn)
            if resn+1 in lcr_sets[chain]:
                curr_streak.append(resn+1)
            else:
                if len(curr_streak) >= parameters['plddt_params']['contiguous_res']:
                    bfact_sampling_targets[chain].append(curr_streak)
                curr_streak = []

    return bfact_sampling_targets


def setup_pose(
    input_structure: str,
    make_centroid: bool = False,
    ) -> pyrosetta.rosetta.core.pose.Pose:
    """Initialize a Pose object from a sequence, a .txt file containing the sequence or a PDB file.
     
    If desired, the created Pose object is then changed to 'centroid' configuration.

    Args:
        input_structure (str):
            Filepath to the input .pdb structure, .txt with sequence or the actual sequence string.
        make_centroid (str, optional):
            Whether to convert the Pose side-chains to centroid configuration. Defaults to False.

    Returns:
        pyrosetta.rosetta.core.pose.Pose:
            A PyRosetta Pose object initialized from the input structure/sequence.
    """
    pose = None
    if input_structure.endswith('.pdb'):
        # Returns PackedPose so we need to convert to Pose
        pose = io.to_pose(io.pose_from_file(input_structure))

    elif input_structure.endswith('.txt'):
        with open(input_structure,'r',encoding='utf-8') as input_sequence:
            # Returns PackedPose so we need to convert to Pose
            pose = io.to_pose(io.pose_from_sequence(input_sequence.read().strip()))

    else:
        # Returns PackedPose so we need to convert to Pose
        pose = io.to_pose(io.pose_from_sequence(input_structure))

    assert pose is not None, 'Invalid input structure/sequence!'

    # Swap to centroid configuration if requested
    if make_centroid:
        pyrosetta.rosetta.protocols.simple_moves.SwitchResidueTypeSetMover('centroid').apply(pose)

    return pose


def setup_minmover(
    scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction,
    min_id: str,
    tolerance: float,
    max_iters: int | None = None,
    dofs: tuple[str,str] = ('bb','chi'),
    ) -> pyrosetta.rosetta.protocols.minimization_packing.MinMover:
    """Setup a PyRosetta MinMover object given the necessary parameters.

    Args:
        scorefxn (pyrosetta.rosetta.core.scoring.ScoreFunction):
            Score function that will be used during Pose minimization.
        min_id (str):
            Identifier for the used PyRosetta minimization algorithm.
        tolerance (float):
            Value for the MinMover tolerance.
        max_iters (int):
            Maximum iterations of the MinMover. Defaults to None, meaning
            the MinMover object's default value.
        dofs (tuple[str,str], optional):
            Defines what angles to set as flexible during minimization.
            Defaults to backbone and sidechain, i.e. ('bb','chi').

    Returns:
        pyrosetta.rosetta.protocols.minimization_packing.MinMover:
            A PyRosetta MinMover object setup with desired parameters.
    """
    # Setup the MoveMap object with desired degrees of freedom
    mmap = pyrosetta.rosetta.core.kinematics.MoveMap()
    if 'bb' in dofs:
        mmap.set_bb(True) # we want to modify backbone torsion angles (phi psi)
    if 'chi' in dofs:
        mmap.set_chi(True) # and chi torsion angles (side chains)
    
    # Setup the MinMover object with required paremeters
    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover(mmap,
                                                                          scorefxn,
                                                                          min_id,
                                                                          tolerance,
                                                                          True) # whether to use neighbour list
    if max_iters is not None:
        min_mover.max_iter(max_iters)

    return min_mover

def derive_constraint_targets(
    pose: pyrosetta.rosetta.core.pose.Pose,
    sampling_targets: dict[str,tuple[tuple[str,tuple[int,...],str,str],...]],
    ) -> tuple[tuple[int,int],...]:
    """Derive the list of residues to keep constrained based on sampling targets.
    
    Given a Pose and the target residue ranges for sampling, mark all non-sampled
    residues as constraint targets.
    In the case of a multichain input structure, assumes chains are properly labeled.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Initial Pose object for sampling.
        sampling_targets (dict[str,tuple[tuple[str,tuple[int,...],str,str],...]):
            Dictionary detailing the target regions for sampling in each chain.

    Returns:
        tuple[tuple[int,int],...]:
            All the residue number pairs representing regions on which to apply constraints.
    """
    targets = deepcopy(sampling_targets)

    # Map chain letters (str) to chain ids (int)
    if pose.num_chains() == 1:
        chain_letter = pyrosetta.rosetta.core.pose.get_chain_from_chain_id(1,pose)
        if chain_letter == ' ':
            # In single-chain PDBs chain letter can be empty so we force it to be 'A' here
            chain_letter_ids = { 'A': 1}
        else:
            chain_letter_ids = { chain_letter: 1}

    else: # if there is more than 1 chain we assume they are properly labeled
        chain_letter_ids = {}
        i = 1
        while i < pose.num_chains() + 1:
            chain_letter_ids[pyrosetta.rosetta.core.pose.get_chain_from_chain_id(i,pose)] = i
            i += 1

    constraint_targets = []
    target_chains = list(targets.keys())

    # Constrain all the residues of chains with no sampled regions
    for chain_letter,chain_id in chain_letter_ids.items():
        if chain_letter not in target_chains:
            chain_start = pose.chain_begin(chain_id)
            chain_end = pose.chain_end(chain_id)
            constraint_targets.append([chain_start,chain_end])

    # Constrain non-sampled regions of chains with sampled regions
    for i,chain in enumerate(target_chains):
        # Get chain start and end residue numbers
        chain_start = pose.chain_begin(chain_letter_ids[chain])
        chain_end = pose.chain_end(chain_letter_ids[chain])

        # Get tuple of target tuples for this chain
        chain_targets = targets[chain]

        # Get starting and ending residues of target regions for sampling
        start_residues = []
        end_residues = []
        for target in chain_targets: # e.g. target = ('MC' , (X1,Y1), 'all', 'TRIPEPTIDE')
            # Get start and end residue numbers for target region
            start_res = pyrosetta.rosetta.core.pose.pdb_to_pose(pose, target[1][0], chain)
            end_res = pyrosetta.rosetta.core.pose.pdb_to_pose(pose, target[1][-1], chain)
            start_res_sampled = False
            end_res_sampled = False

            # Chain start and end edge cases are deal with in if block below
            if start_res == chain_start:
                start_res_sampled = True
            else:
                start_residues.append(start_res - 1)
            if end_res == chain_end:
                end_res_sampled = True
            else:
                end_residues.append(end_res + 1)

        # If we sample the first res but not the last
        if len(start_residues) < len(end_residues):
            start_residues.append(chain_end)

        # If we sample the last res but not the first
        elif len(start_residues) > len(end_residues):
            end_residues = [chain_start] + end_residues

        # Account for edge cases
        else:
            if not start_res_sampled:
                start_residues.append(chain_end)
            if not end_res_sampled:
                end_residues = [chain_start] + end_residues

        # Assign constraint targets
        while len(end_residues) > 0 and len(start_residues) > 0:
            start_res = start_residues.pop(0)
            end_res = end_residues.pop(0)
            # We want to constrain the non-sampled regions so we grab end_res first
            constraint_targets.append((end_res,start_res))

    return tuple(constraint_targets)


def apply_pae_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    pae_filepath: str,
    plddt_targets: dict,
    cutoff: float = 10.0,
    flatten_cutoff: float = 10.0,
    flatten_value: float = 10.0,
    weight: float = 1.0,
    tolerance: float | None = None,
    adjacency_threshold: int = 8,
    plddt_scaling_factor: float = 1.0):
    """Apply constraints to a Pose created from an AlphaFold structure
    based on the Predicted Aligned Error (PAE) matrix.
    
    The strength of the constraints scales with the value of the predicted aligned error for that
    residue pair.
    After scaling with PAE value, applied constraints are made weaker when between a residue with
    high pLDDT and a residue with low pLDDT.
    This function should be used instead of the apply_constraints function, they should not be used
    in tandem.
    
    Adapted from:
        https://github.com/matteoferla/pyrosetta-help
    
    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            The Pose constraints will be applied to.
        pae_filepath (str):
            Path to the pae matrix .json file.
        plddt_targets (dict):
            Mapping of each chain to the residue numbers that
            will be sampled (pdb numbering), all with plddt below a threshold.
        cutoff (float):
            Only consider PAE values below this number (low error).
        flatten_cutoff (float):
            Any PAE values below this value will be changed to match flatten value.
        flatten_value (float):
            Any PAE values below flatten cutoff will be changed to match this value.
        weight (float):
            Along with the error value, determines the strength of the applied AtomPairConstraints
        tolerance (float, optional):
            Defines the tolerance of the FlatHarmonicFunction of AtomPairConstraints
            created from the PAE matrix. Defaults to None.
        adjacency_threshold (int):
            How far away two residues need to be to consider their PAE value. Neighbours are skipped
            as PAE is best used for determining between domain or between chain confidence.
        plddt_scaling_factor (float):
            Any constraints setup between residues where one of them has a low pLDDT and another a
            high pLDDT will be scaled by multiplying its weight by this factor. The higher this
            value the weaker those constraints will be.
    """
    # Get PAE Matrix
    with open(pae_filepath,'r',encoding='utf-8-sig') as f:
        pae_content = json.load(f)

    # Check if PAE Matrix comes from UniProt accession download, correct its type if so
    if isinstance(pae_content,list):
        pae_content = pae_content[0]

    # Read PAE Matrix
    try:
        pae_matrix = np.array(pae_content['predicted_aligned_error'])
    except KeyError:
        try:
            pae_matrix = np.array(pae_content['pae'])
        except KeyError as e:
            raise AssertionError('PAE matrix in the .json file must have key '
                                 '"predicted_aligned_error" or "pae".') from e

    # Check which residues have pLDTT below threshold
    low_plddt_res = []
    for chain_id,targets in plddt_targets.items():
        for target in targets:
            res_range = target[1]
            for res in res_range:
                low_plddt_res.append(pose.pdb_info().pdb2pose(chain_id,res))

    tmp_pose = pyrosetta.rosetta.core.pose.Pose()
    tmp_pose.detached_copy(pose)

    # Create constraint set
    working_cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet(pose.constraint_set())

    # Apply constraints based on pae
    for r1_idx, r2_idx in np.argwhere(pae_matrix < cutoff):
        # Between two residues with low plddt, pae value is
        # often high so its already discarded here

        if abs(r1_idx - r2_idx) < adjacency_threshold:
            # This will also discard constraints between residues and themselves
            continue

        elif r1_idx <= r2_idx:
            # Do not add redundant constraints (more constraints mean longer energy calculations)
            continue

        elif (r1_idx+1 not in low_plddt_res and r2_idx+1 in low_plddt_res or \
              r1_idx+1 in low_plddt_res and r2_idx+1 not in low_plddt_res):
            # Add weaker constraints between residues when one of them has low plddt and
            # the other has high plddt

            error = pae_matrix[r1_idx, r2_idx]
            if error < flatten_cutoff:
                error = flatten_value

            ca1 = pyrosetta.rosetta.core.id.AtomID(pose.residue(r1_idx + 1).atom_index('CA'),
                                                   r1_idx + 1)
            ca2 = pyrosetta.rosetta.core.id.AtomID(pose.residue(r2_idx + 1).atom_index('CA'),
                                                   r2_idx + 1)
            ca1_xyz = pose.residue(ca1.rsd()).xyz(ca1.atomno())
            ca2_xyz = pose.residue(ca2.rsd()).xyz(ca2.atomno())
            d = (ca1_xyz - ca2_xyz).norm()

            if tolerance is not None:
                # Higher sd_in -> Lower apc value -> Lower Pose score => Weaker cst
                fun = pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc(
                    x0_in=d,
                    sd_in=error*weight*plddt_scaling_factor,
                    tol_in=tolerance
                    )
            else:
                fun = pyrosetta.rosetta.core.scoring.func.HarmonicFunc(
                    x0_in=d,
                    sd_in=error*weight*plddt_scaling_factor
                    )

            apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(ca1,ca2,fun)
            working_cs.add_constraint(apc)

        elif r1_idx+1 not in low_plddt_res and r2_idx+1 not in low_plddt_res:
            # Add stronger (not weakened) constraints between res when both of them have high plddt
            # This if clause also includes pairs of non-sampled residues, if they have low PAE
            error = pae_matrix[r1_idx, r2_idx]
            if error < flatten_cutoff:
                error = flatten_value

            ca1 = pyrosetta.rosetta.core.id.AtomID(pose.residue(r1_idx + 1).atom_index('CA'),
                                                   r1_idx + 1)
            ca2 = pyrosetta.rosetta.core.id.AtomID(pose.residue(r2_idx + 1).atom_index('CA'),
                                                   r2_idx + 1)
            ca1_xyz = pose.residue(ca1.rsd()).xyz(ca1.atomno())
            ca2_xyz = pose.residue(ca2.rsd()).xyz(ca2.atomno())
            d = (ca1_xyz - ca2_xyz).norm()

            if tolerance is not None:
                # Lower sd_in -> Higher apc value -> Higher Pose score => Stronger cst
                fun = pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc(x0_in=d,
                                                                           sd_in=error*weight,
                                                                           tol_in=tolerance)
            else:
                fun = pyrosetta.rosetta.core.scoring.func.HarmonicFunc(x0_in=d,
                                                                       sd_in=error*weight)

            apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(ca1,ca2,fun)
            working_cs.add_constraint(apc)

    # Apply constraint set to pose
    setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    setup.constraint_set(working_cs)
    setup.apply(tmp_pose)

    # Update working pose
    pose.assign(tmp_pose)


def apply_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    cst_targets: tuple[tuple[int,int],...],
    contacts: tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...] | None = None,
    stdev: float | None = 10.0,
    tolerance: float | None = 0.001):
    """Apply energy constraints to desired regions of a Pose object.

    This function is incompatible with the apply_pae_constraints function, and they should not be
    used on the same Pose.
    The word 'contacts' is used here for convenience, these could be any two regions whose
    relative position should be conserved, even if they are far apart in the Pose.
    

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Target Pose object for constraints. It is modified in place.
        cst_targets (tuple[tuple[int,int],...]):
            Residue ranges defining regions that make up folded protein domains.
            AtomPairConstraints will be applied between constituting residues.
        contacts (tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...]):
            Pairs of residue ranges defining regions whose relative position should be conserved.
            AtomPairConstraints will be applied between residues belonging to different regions.
            Defaults to None.
        stdev (float, optional):
            Standard deviation value to use in constraints. Defaults to 10.0.
        tolerance (float, optional):
            Tolerance value to use in constraints. Defaults to 0.001.
    """
    # Initialize temporary Pose object
    tmp_pose = pyrosetta.rosetta.core.pose.Pose()
    tmp_pose.detached_copy(pose)

    # Initialize working ConstraintSet object
    cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet(pose.constraint_set())

    # Add to ConstraintSet
    ## Conserve the structure of folded domains made up of contiguous residues
    cs_intra = add_intrachain_constraints(pose=tmp_pose,
                                          constraint_targets=cst_targets,
                                          constraint_set=cs,
                                          stdev=stdev,
                                          tolerance=tolerance)

    ## Additional constraints between regions that are not contiguous in the Pose,
    # but whose relative position should be conserved
    cs_intra_contacts = add_contacts_constraints(pose=tmp_pose,
                                                 contacts=contacts,
                                                 constraint_set=cs_intra,
                                                 stdev=stdev,
                                                 tolerance=tolerance)

    # Update ConstraintSet of temporary Pose
    setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    setup.constraint_set(cs_intra_contacts)
    setup.apply(tmp_pose)

    # Overwrite Pose with final ConstraintSet
    pose.assign(tmp_pose)


def setup_fold_tree(
    pose: pyrosetta.rosetta.core.pose.Pose,
    constraint_targets: tuple[tuple[int,int],...],
    contacts: tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...] | None ):
    """Change a Pose's FoldTree in order to minimize "lever arm" effects during sampling.
        
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

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Pose object whose FoldTree will be updated.
        constraint_targets (tuple[tuple[int,int],...]):
            Residues between which AtomPairConstraints will be applied.
        contacts (tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...]):
            Residue ranges where two chains are interacting.
    
    Reference:
        See https://docs.rosettacommons.org/demos/latest/tutorials/fold_tree/fold_tree
        for more information about the Rosetta FoldTree.
    """
    # Prepare working pose
    ft_pose = pyrosetta.rosetta.core.pose.Pose()
    ft_pose.detached_copy(pose)

    #central_residues = [227,483,739] # N246TRIMER (sets the fold tree manually)
    #central_residues = [227,483] # N246DIMER (sets the fold tree manually)
    #central_residues = [227] # N246MONOMER
    #central_residues = [52,133,214] # NLINKERTRIMER
    #central_residues = [227,656,1085,1514] # NFullLengthTetramer (sets the fold tree manually)
    # Uncomment what is below (so that the FoldTree is set automatically)

    # Calculate the optimal central residues
    central_residues = []
    for i in range(1,ft_pose.num_chains()+1):
        chain_start = ft_pose.chain_begin(i)
        chain_end = ft_pose.chain_end(i)
        ideal_central_res = (chain_start + chain_end) // 2

        # Start from any value in constrained res range
        central_res = constraint_targets[0][0]

        # Get how far the current res is from the ideal mid-point
        minimum_distance = abs(central_res - ideal_central_res)

        for res_range in constraint_targets:
            for res in range(res_range[0],res_range[1]+1):
                if ft_pose.chain(res) == i:
                    dist = abs(res - ideal_central_res) # distance between residue numbers
                    if dist < minimum_distance:
                        minimum_distance = dist
                        central_res = res
        central_residues.append(central_res)

    # If there are contacts, use them to update central residues if needed
    # This is important in multi-chain proteins, to avoid cases where the central residue of a
    # chain is not in a contacted region, which would make it so that folded domain that central
    # residue belongs to would not move relative to other chains during sampling
    if contacts is not None:
        chains_contact_regions = {}
        for contact in contacts:
            # contact: ( ('X', (x1,x2) ) , ( 'Y', (y1,y2) ) )

            # Chain X
            chain_x = contact[0][0]

            # Contact residue range for X
            x_start = contact[0][1][0]
            x_end = contact[0][1][1]
            inter_range_x = set([pyrosetta.rosetta.core.pose.pdb_to_pose(
                                    ft_pose,
                                    res_id,
                                    chain_x)
                                for res_id in range(x_start, x_end+1) ])
            try:
                if inter_range_x not in chains_contact_regions[chain_x]:
                    chains_contact_regions[chain_x].append(inter_range_x)
            except KeyError:
                chains_contact_regions[chain_x] = [inter_range_x]

            # Chain Y
            chain_y = contact[1][0]

            # Contact residue range for Y
            y_start = contact[1][1][0]
            y_end = contact[1][1][1]
            inter_range_y = [pyrosetta.rosetta.core.pose.pdb_to_pose(
                                ft_pose,
                                res_id,
                                chain_y)
                            for res_id in range(y_start,y_end+1) ]
            try:
                if inter_range_y not in chains_contact_regions[chain_y]:
                    chains_contact_regions[chain_y].append(inter_range_y)
            except KeyError:
                chains_contact_regions[chain_y] = [inter_range_y]

        for chain,regions in chains_contact_regions.items():
            regions_lists = []
            for region in regions:
                regions_lists.append(sorted([x for x in region]))
            chains_contact_regions[chain] = regions_lists

        updated_central_residues = []
        for cen_res in central_residues:
            # Get chain of current central res
            res_chain = ft_pose.pdb_info().chain(cen_res)

            # Get regions of that chain that are in contact
            chain_contact_regions = chains_contact_regions[res_chain]

            # See if the central residue is already inside a region that is in contact
            in_contact_region = False
            for contact_region in chain_contact_regions:
                if cen_res in contact_region:
                    in_contact_region = True
                    break

            # If not, replace it with the nearest residue that is inside a region in contact
            if not in_contact_region:
                distances = {}
                for region in chain_contact_regions:
                    distances[region[0]] = abs(cen_res - region[0])
                    distances[region[1]] = abs(cen_res - region[1])
                min_distance = min(distances.values())
                for resnum, distance in distances.items():
                    if distance == min_distance:
                        updated_central_residues.append(resnum)
            else:
                updated_central_residues.append(cen_res)

        # Update central residues with new residue numbers
        central_residues = updated_central_residues

    # Update FoldTree using new 'central residues'
    ft = ft_pose.fold_tree()

    # Split tree
    for cen_res in central_residues:
        ft.split_existing_edge_at_residue(cen_res)

    # Get chain starting residues
    chain_starts = []
    for chain_num in range(1,ft_pose.num_chains()+1):
        chain_starts.append(ft_pose.chain_begin(chain_num))

    # Delete old jumps
    jump_counter = 1
    for i in chain_starts[1:]:
        ft.delete_unordered_edge(chain_starts[0],i,jump_counter)
        jump_counter += 1

    # Set new jumps between chain's new parents
    jump_counter = 1
    for i in central_residues[1:]:
        ft.add_edge(central_residues[0],i,jump_counter)
        jump_counter += 1

    # Reorder tree to flow to N and C terminals
    for cen_res in central_residues:
        ft.reorder(cen_res)

    # Check validity before continuing
    assert ft.check_fold_tree(), print('Invalid FoldTree setup:\n',ft)
    pose.fold_tree(ft)


def _prep_target(seq_len: int,target: list[int]) -> list[int]:
    """In target, replace seq_len with seq_len -1 if it is present.
    
    Check the given target residue range. If the last residue of the
    chain is included replace it with the second to last to avoid
    attempting to sample a fragment with only two residues.

    Args:
        seq_len (int):
            Length of the chain of the protein being sampled.
        target (list[int]):
            Range of target residues for MC sampler.

    Returns:
        list[int]:
            Updated range of target residues for MC sampler.
    """
    target_new = deepcopy(target)

    if target_new[-1] == seq_len:
        try:
            target_new[-2]
        except IndexError:
            # If the last target range is only the last res, replace it with second to last
            target_new = [seq_len - 1]
        else:
            # If last target range is more than only last res, remove last element
            target_new = target_new[:-1]

    return target_new
