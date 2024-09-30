"""Auxiliary functions for sampling."""

# IMPORTS
## Standard Library Imports
from copy import deepcopy
from typing import Optional,Union
import json

## Third Party Imports
import numpy as np
import pyrosetta
import pyrosetta.distributed.io as io
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pose import Pose, pdb_to_pose, get_chain_from_chain_id
from pyrosetta.rosetta.core.scoring import constraints,func
from pyrosetta.rosetta.protocols.constraint_movers import ConstraintSetMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover

## Local Imports
from ensemblify.utils import df_from_pdb


# FUNCTIONS
def add_intrachain_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    constraint_targets: tuple[tuple[int,int],...],
    constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet,
    constraint_function: Union[pyrosetta.rosetta.core.scoring.func.HarmonicFunc,pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc],
    stdev: float,
    tolerance: Optional[float] = None,
) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
    """
    Add constraints between non-sampled residues in a chain to a constraint set.
    
    Create all AtomPairConstraints between residues of a Pose object
    present in constraint_targets and add them to a ConstraintSet.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            target Pose object for constraints.
        constraint_targets:
            residues between which AtomPairConstraints will be applied.
        constraint_set (pyrosetta.rosetta.core.scoring.constraints.ConstraintSet):
            set of constraints to later be applied to Pose.
        constraint_function (pyrosetta.rosetta.core.scoring.func.HarmonicFunc | pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc):
            function to use for each added constraint.
        stdev:
            standard deviation value to use in constraints.
        tolerance:
            tolerance value to use in constraints (if applicable). Defaults to None.
    
    Returns:
        working_cs (pyrosetta.rosetta.core.scoring.constraints.ConstraintSet):
            updated constraint set, with intrachain constraints.
    """
    AtomPairConstraint = constraints.AtomPairConstraint
    working_cs = constraints.ConstraintSet()
    working_cs.detached_copy(constraint_set)
    for target in constraint_targets:
        target_range = list(range(target[0],target[1]+1))
        for i,res_id_1 in enumerate(target_range):
            for res_id_2 in target_range[i+1:]:
                a1 = AtomID(pose.residue(res_id_1).atom_index('CA'), res_id_1)
                a2 = AtomID(pose.residue(res_id_2).atom_index('CA'), res_id_2)
                a1_xyz = pose.residue(a1.rsd()).xyz(a1.atomno())
                a2_xyz = pose.residue(a2.rsd()).xyz(a2.atomno())
                d = (a1_xyz - a2_xyz).norm()

                if tolerance:
                    apc = AtomPairConstraint(a1,a2, constraint_function(x0_in=d,
                                                                        sd_in=stdev,
                                                                        tol_in=tolerance))
                else:
                    apc = AtomPairConstraint(a1,a2, constraint_function(x0_in=d,
                                                                        sd_in=stdev))

                working_cs.add_constraint(apc)

    return working_cs


def add_contacts_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    contacts: tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...] | None,
    constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet,
    constraint_function: Union[pyrosetta.rosetta.core.scoring.func.HarmonicFunc,pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc],
    stdev: float,
    tolerance: Optional[float] = None,
    ) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
    """
    Add constraints between residues of different chains/domains that must remain in contact.

    Create all contact constraints (dimerization sites, intrachain folding) targeting residues
    of a Pose object and add them to a ConstraintSet.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            target Pose object for constraints.
        contacts:
            residue ranges where two regions are interacting.
        constraint_set (pyrosetta.rosetta.core.scoring.constraints.ConstraintSet):
            set of constraints to later be applied to Pose.
        constraint_function (pyrosetta.rosetta.core.scoring.func.HarmonicFunc | pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc):
            function to use for each added constraint.
        stdev:
            standard deviation value to use in constraints.
        tolerance:
            tolerance value to use in constraints (if applicable). Defaults to None.
    
    Returns:
        working_cs (pyrosetta.rosetta.core.scoring.constraints.ConstraintSet):
            updated constraint set, with contact constraints.
    """
    AtomPairConstraint = constraints.AtomPairConstraint # for interactions

    working_cs = constraints.ConstraintSet()
    working_cs.detached_copy(constraint_set)

    if contacts is not None:
        for contact in contacts:
            # contact: ( ('X', (x1,x2) ) , ( 'Y', (y1,y2) ) )

            # Chain X
            chain_x = contact[0][0]

            # Contact residue range for X
            x_start = contact[0][1][0]
            x_end = contact[0][1][1]
            inter_range_x = [ pdb_to_pose(pose, res_id, chain_x)
                              for res_id in range(x_start, x_end+1) ]

            # Chain Y
            chain_y = contact[1][0]

            # Contact residue range for Y
            y_start = contact[1][1][0]
            y_end = contact[1][1][1]
            inter_range_y = [ pdb_to_pose(pose, res_id, chain_y)
                              for res_id in range(y_start,y_end+1) ]

            # Apply inter-region constraints between X and Y
            for res_id_x in inter_range_x:
                for res_id_y in inter_range_y:
                    a1 = AtomID(pose.residue(res_id_x).atom_index('CA'), res_id_x)
                    a2 = AtomID(pose.residue(res_id_y).atom_index('CA'), res_id_y)
                    a1_xyz = pose.residue(a1.rsd()).xyz(a1.atomno())
                    a2_xyz = pose.residue(a2.rsd()).xyz(a2.atomno())
                    d = (a1_xyz - a2_xyz).norm()

                    if tolerance:
                        apc = AtomPairConstraint(a1,a2, constraint_function(x0_in=d,
                                                                            sd_in=stdev,
                                                                            tol_in=tolerance))
                    else:
                        apc = AtomPairConstraint(a1,a2, constraint_function(x0_in=d,
                                                                            sd_in=stdev))

                    working_cs.add_constraint(apc)

    #setup constraints for other experimental info here similarly

    return working_cs


def get_targets_from_plddt(parameters: dict) -> dict[str,list[int]]:
    """Get, for each chain, lists of residues with pLDDT value below the threshold.

    The input structure defined in the parameters dictionary must be an AlphaFold model,
    i.e. have the pLDDT value for each residue in the .pdb B-Factor column.

    Args:
        parameters:
            dictionary following Ensemblify parameters template.
    
    Returns:
        bfact_sampling_targets:
            mapping of each chain to the residue numbers contained in it pertaining
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


def setup_pose(input_structure: str) -> pyrosetta.rosetta.core.pose.Pose:
    """Initialize a Pose object from a sequence or a PDB file.
     
    The created Pose object is changed to 'centroid' configuration.

    Args:
        input_structure:
            filepath to the input .pdb structure or .txt sequence.

    Returns:
        initial_pose (pyrosetta.rosetta.core.pose.Pose):
            our initial Pose for sampling.
    """
    # Returns PackedPose so we need to convert to Pose
    initial_pose = io.to_pose(io.pose_from_file(input_structure))

    # Swap to centroid conformation
    SwitchResidueTypeSetMover('centroid').apply(initial_pose)

    return initial_pose


def setup_minmover(
    scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction,
    min_id: str,
    tolerance: float,
    max_iters: int,
    dofs: tuple[str,str] = ('bb','chi'),
    ) -> pyrosetta.rosetta.protocols.minimization_packing.MinMover:
    """Setup the MoveMap and MinMover for last minimization steps in the sampling process.

    Args:
        scorefxn (pyrosetta.rosetta.core.scoring.ScoreFunction):
            score function used during sampling to evaluate our Pose conformations.
        min_id:
            identifier for the PyRosetta minimization algorithm.
        tolerance:
            value for the MinMover tolerance.
        max_iters:
            maximum iterations of the MinMover.
        dofs:
            defines what angles to set as flexible during minimization.
            Defaults to backbone and sidechain.

    Returns:
        min_mover (pyrosetta.rosetta.protocols.minimization_packing.MinMover):
            PyRosetta MinMover for last minimization steps in the sampling process.
    """
    # Setup the MoveMap and MinMover for last minimization step
    mmap = MoveMap()
    for dof in dofs:
        if dof == 'bb':
            mmap.set_bb(True) # we want to modify backbone torsion angles (phi psi)
        elif dof == 'chi':
            mmap.set_chi(True) # and chi torsion angles (side chains)

    min_mover = MinMover(mmap,scorefxn,min_id,tolerance,True)
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
            initial Pose object for sampling.
        sampling_targets:
            dictionary detailing the target regions for sampling in each chain.

    Returns:
        constraint_targets:
            all the residue numbers pairs representing regions on which to apply constraints.
    """
    targets = deepcopy(sampling_targets)

    # Map chain letters (str) to chain ids (int)
    if pose.num_chains() == 1:
        chain_letter = get_chain_from_chain_id(1,pose)
        if chain_letter == ' ':
            # In single-chain PDBs chain letter can be empty so we force it to be 'A' here
            chain_letter_ids = { 'A': 1}
        else:
            chain_letter_ids = { chain_letter: 1}

    else: # if there is more than 1 chain we assume they are properly labeled
        chain_letter_ids = {}
        i = 1
        while i < pose.num_chains() + 1:
            chain_letter_ids[get_chain_from_chain_id(i,pose)] = i
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
            start_res = pdb_to_pose(pose, target[1][0], chain)
            end_res = pdb_to_pose(pose, target[1][-1], chain)
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

    cst_targets = tuple(constraint_targets)

    return cst_targets


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
            the Pose constraints will be applied to.
        pae_filepath:
            path to the pae matrix .json file.
        plddt_targets:
            mapping of each chain to the residue numbers that
            will be sampled (pdb numbering), all with plddt below a threshold.
        cutoff:
            only consider PAE values below this number (low error).
        flatten_cutoff:
            any PAE values below this value will be changed to match flatten value.
        flatten_value:
            any PAE values below flatten cutoff will be changed to match this value.
        weight:
            along with the error value, determines the strength of the applied AtomPairConstraints
        tolerance:
            if given, defines the tolerance of the FlatHarmonicFunction of AtomPairConstraints
            created from the PAE matrix. Defaults to None.
        adjacency_threshold:
            how far away two residues need to be to consider their PAE value. Neighbours are skipped
            as PAE is best used for determining between domain or between chain confidence.
        plddt_scaling_factor:
            any constraints setup between residues where one of them has a low pLDDT and another a high
            pLDDT will be scaled by multiplying its weight by this factor. The higher this value the weaker
            those constraints will be.
    """

    # Get error matrix
    with open(pae_filepath,'r',encoding='utf-8-sig') as f:
        pae_content = json.load(f)

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

    cst_pose = Pose()
    cst_pose.detached_copy(pose)

    # Create constraint set
    working_cs = constraints.ConstraintSet(pose.constraint_set())

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

            FlatHarmonicFunc = func.FlatHarmonicFunc
            HarmonicFunc = func.HarmonicFunc
            AtomPairConstraint = constraints.AtomPairConstraint

            ca1 = AtomID(pose.residue(r1_idx + 1).atom_index('CA'), r1_idx + 1)
            ca2 = AtomID(pose.residue(r2_idx + 1).atom_index('CA'), r2_idx + 1)
            ca1_xyz = pose.residue(ca1.rsd()).xyz(ca1.atomno())
            ca2_xyz = pose.residue(ca2.rsd()).xyz(ca2.atomno())
            d = (ca1_xyz - ca2_xyz).norm()

            if tolerance is not None:
                # Higher sd_in -> Lower apc value -> Lower Pose score => Weaker cst
                fun = FlatHarmonicFunc(x0_in=d,
                                       sd_in=error*weight*plddt_scaling_factor,
                                       tol_in=tolerance)
            else:
                fun = HarmonicFunc(x0_in=d,
                                   sd_in=error*weight)

            apc = AtomPairConstraint(ca1,ca2,fun)
            working_cs.add_constraint(apc)

        elif r1_idx+1 not in low_plddt_res and r2_idx+1 not in low_plddt_res:
            # Add stronger (not weakened) constraints between res when both of them have high plddt
            error = pae_matrix[r1_idx, r2_idx]
            if error < flatten_cutoff:
                error = flatten_value

            FlatHarmonicFunc = func.FlatHarmonicFunc
            HarmonicFunc = func.HarmonicFunc
            AtomPairConstraint = constraints.AtomPairConstraint

            ca1 = AtomID(pose.residue(r1_idx + 1).atom_index('CA'), r1_idx + 1)
            ca2 = AtomID(pose.residue(r2_idx + 1).atom_index('CA'), r2_idx + 1)
            ca1_xyz = pose.residue(ca1.rsd()).xyz(ca1.atomno())
            ca2_xyz = pose.residue(ca2.rsd()).xyz(ca2.atomno())
            d = (ca1_xyz - ca2_xyz).norm()

            if tolerance is not None:
                # Lower sd_in -> Higher apc value -> Higher Pose score => Stronger cst
                fun = FlatHarmonicFunc(x0_in=d,
                                    sd_in=error*weight,
                                    tol_in=tolerance)
            else:
                fun = HarmonicFunc(x0_in=d,
                                   sd_in=error*weight)

            apc = AtomPairConstraint(ca1,ca2,fun)
            working_cs.add_constraint(apc)

    # Apply constraint set to pose
    setup = ConstraintSetMover()
    setup.constraint_set(working_cs)
    setup.apply(cst_pose)

    # Update working pose
    pose.assign(cst_pose)


def apply_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    cst_targets: tuple[tuple[int,int],...],
    contacts: tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...] | None,
    stdev: float = 10.0,
    tolerance: float | None = 0.001):
    """Apply constraints to non-sampled regions of a Pose object.

    Apply all appropriate intra-region constraints based on cst_targets and all interregion
    constraints based on given restraints.
    not be sampled.
    This function is incompatible with the apply_pae_constraints function, and they should not be
    used in tandem.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            target Pose object for constraints.
        cst_targets:
            residues between which AtomPairConstraints will be applied.
        contacts:
            residue ranges where two chains are interacting.
        stdev:
            standard deviation value to use in constraints.
        tolerance:
            tolerance value to use in constraints (if applicable). Defaults to None.
    """
    cst_pose = Pose()
    cst_pose.detached_copy(pose)

    # Setup constraint function
    if tolerance is not None:
        cst_func = func.FlatHarmonicFunc
    else:
        cst_func = func.HarmonicFunc

    # Create constraint set
    cs = constraints.ConstraintSet(pose.constraint_set())

    # Add to constraint set
    ## Constraints for non sampled regions (conserve secondary structure)
    cs_intra = add_intrachain_constraints(pose=cst_pose,
                                          constraint_targets=cst_targets,
                                          constraint_set=cs,
                                          constraint_function=cst_func,
                                          stdev=stdev,
                                          tolerance=tolerance)

    ## Additional constraints to direct sampling (conserve tertiary and quaternary structure)
    cs_intra_contacts = add_contacts_constraints(pose=cst_pose,
                                                 contacts=contacts,
                                                 constraint_set=cs_intra,
                                                 constraint_function=cst_func,
                                                 stdev=stdev,
                                                 tolerance=tolerance)

    # Apply constraint set to pose
    setup = ConstraintSetMover()
    setup.constraint_set(cs_intra_contacts)
    setup.apply(cst_pose)

    # Update working pose
    pose.assign(cst_pose)


def setup_fold_tree(
    pose: pyrosetta.rosetta.core.pose.Pose,
    constraint_targets: tuple[tuple[int,int],...]):
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

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Pose object whose FoldTree will be updated.
        constraint_targets:
            residues between which AtomPairConstraints will be applied.
    """
    # Prepare working pose
    ft_pose = Pose()
    ft_pose.detached_copy(pose)

    #central_residues = [227,483,739] # N246TRIMER (sets the fold tree manually)
    #central_residues = [227,483] # N246DIMER (sets the fold tree manually)
    #central_residues = [227] # N246MONOMER

    #central_residues = [52,133,214] # NLINKERTRIMER

    #central_residues = [227,656,1085,1514] # NFullLengthTetramer (sets the fold tree manually)

    # Uncomment what is below (so that the FoldTree is set automatically)
    central_residues = []
    for i in range(1,ft_pose.num_chains()+1):
        chain_start = ft_pose.chain_begin(i)
        chain_end = ft_pose.chain_end(i)
        ideal_central_res = chain_start + chain_end // 2

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


def prep_target(seq_len: int,target: list[int]) -> list[int]:
    """In target, replace seq_len with seq_len -1 if it is present.
    
    Check the given target residue range. If the last residue of the
    chain is included replace it with the second to last to avoid
    attempting to sample a fragment with only two residues.

    Args:
        seq_len:
            length of the chain of the protein being sampled.
        target:
            range of target residues for MC sampler.

    Returns:
        target_new:
            updated range of target residues for MC sampler.
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
