"""Auxiliary functions for applying energy restraints (PyRosetta Constraints) to Pose objects."""

# IMPORTS
## Standard Library Imports
import json

## Third Party Imports
import numpy as np
import pyrosetta

# FUNCTIONS
def make_atom_pair_constraint(
    pose: pyrosetta.rosetta.core.pose.Pose,
    res_num_1: int,
    res_num_2: int,
    stdev: float | None = 10.0,
    tolerance: float | None = None
    ) -> pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint:
    """Create an AtomPairConstraint between the alpha carbons of two residues in a Pose.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Target Pose object for constraints. Used to extract residue ID and coordinates.
        res_num_1 (int):
            The first residue number.
        res_num_2 (int):
            The second residue number.
        stdev (float, optional):
            Standard deviation value to use in constraints. Defaults to 10.0.
        tolerance (float, optional):
            Tolerance value to use in constraints. Defaults to None.

    Returns:
        pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint:
            The created AtomPairConstraint object.
    """
    ca1 = pyrosetta.rosetta.core.id.AtomID(pose.residue(res_num_1).atom_index('CA'), res_num_1)
    ca2 = pyrosetta.rosetta.core.id.AtomID(pose.residue(res_num_2).atom_index('CA'), res_num_2)
    ca1_xyz = pose.residue(ca1.rsd()).xyz(ca1.atomno())
    ca2_xyz = pose.residue(ca2.rsd()).xyz(ca2.atomno())
    d = (ca1_xyz - ca2_xyz).norm()

    if tolerance:
        apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
            ca1,
            ca2,
            pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc(x0_in=d,
                                                                 sd_in=stdev,
                                                                 tol_in=tolerance)
        )
    else:
        apc = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
            ca1,
            ca2,
            pyrosetta.rosetta.core.scoring.func.HarmonicFunc(x0_in=d,
                                                             sd_in=stdev)
        )

    return apc

    
def add_intrachain_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    constraint_targets: tuple[tuple[int,int],...],
    constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet,
    stdev: float | None = 10.0,
    tolerance: float | None = None,
    ) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
    """Add constraints between desired residues of the same domain to a constraint set.
    
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
        for i, res_num_1 in enumerate(target_range):
            for res_num_2 in target_range[i+1:]:
                apc = make_atom_pair_constraint(pose=pose,
                                                res_num_1=res_num_1,
                                                res_num_2=res_num_2,
                                                stdev=stdev,
                                                tolerance=tolerance)

                working_cs.add_constraint(apc)

    return working_cs


def add_contacts_constraints(
    pose: pyrosetta.rosetta.core.pose.Pose,
    contacts: tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...] | None,
    constraint_set: pyrosetta.rosetta.core.scoring.constraints.ConstraintSet,
    stdev: float | None = 10.0,
    tolerance: float | None = None,
    ) -> pyrosetta.rosetta.core.scoring.constraints.ConstraintSet:
    """Add constraints between residue regions whose relative position must remain conserved.

    This would include, for example, dimerization sites or long range intrachain folding.
    The word 'contacts' is used here for convenience, these could be any two regions whose
    relative position should be conserved, even if they are far apart in the Pose.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            Target Pose object for constraints.
        contacts (tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...]):
            Residue ranges of regions whose relative position should be conserved. If None,
            no constraints will be added.
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
        for res_num_x in inter_range_x:
            for res_num_y in inter_range_y:
                apc = make_atom_pair_constraint(pose=pose,
                                                res_num_1=res_num_x,
                                                res_num_2=res_num_y,
                                                stdev=stdev,
                                                tolerance=tolerance)

                working_cs.add_constraint(apc)

    return working_cs


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
    """Apply energy constraints to a AlphaFold Pose based on its PAE matrix.
    
    The strength of the constraints scales with the value of the predicted aligned error for that
    residue pair.
    After scaling with PAE value, applied constraints are made weaker when between a residue with
    high pLDDT and a residue with low pLDDT.
    This function is incompatible with the apply_constraints function, and they should not be
    used on the same Pose.
    
    Adapted from:
        https://github.com/matteoferla/pyrosetta-help/blob/main/pyrosetta_help/alphafold/constraints.py
    
    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            The Pose constraints will be applied to.
        pae_filepath (str):
            Path to the PAE matrix .json file.
        plddt_targets (dict):
            Mapping of each chain to the residue numbers that will be sampled (pdb numbering),
            already filtered for residues with pLDDT below a threshold.
        cutoff (float):
            Only consider PAE values below this number (low error).
        flatten_cutoff (float):
            Any PAE values below this value will be changed to match flatten_value.
        flatten_value (float):
            Any PAE values below flatten_cutoff will be changed to match this value.
        weight (float):
            Along with the error value, determines the strength of the applied AtomPairConstraints.
        tolerance (float, optional):
            Defines the tolerance of the FlatHarmonicFunction of AtomPairConstraints created from
            the PAE matrix. Defaults to None.
        adjacency_threshold (int):
            How far away two residues need to be to consider their PAE value. Neighbours are
            skipped as PAE is best used for determining between domain or between chain confidence.
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

    # Initialize temporary Pose object
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

            # Higher stdev -> Lower apc value increase rate -> Lower Pose score => Weaker cst
            apc = make_atom_pair_constraint(pose=tmp_pose,
                                            res_num_1=r1_idx + 1,
                                            res_num_2=r2_idx + 1,
                                            stdev=error*weight*plddt_scaling_factor,
                                            tolerance=tolerance)
            
            working_cs.add_constraint(apc)

        elif r1_idx+1 not in low_plddt_res and r2_idx+1 not in low_plddt_res:
            # Add stronger (not weakened) constraints between res when both of them have high plddt
            # This if clause also includes pairs of non-sampled residues, if they have low PAE
            # (e.g. folded proteins domains whose structure should remain conserved)

            error = pae_matrix[r1_idx, r2_idx]

            if error < flatten_cutoff:
                error = flatten_value
            
            # Lower stdev -> Higher apc value increase rate -> Higher Pose score => Stronger cst
            apc = make_atom_pair_constraint(pose=tmp_pose,
                                            res_num_1=r1_idx + 1,
                                            res_num_2=r2_idx + 1,
                                            stdev=error*weight,
                                            tolerance=tolerance)
            
            working_cs.add_constraint(apc)

    # Apply constraint set to pose
    setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    setup.constraint_set(working_cs)
    setup.apply(tmp_pose)

    # Update working pose
    pose.assign(tmp_pose)
    