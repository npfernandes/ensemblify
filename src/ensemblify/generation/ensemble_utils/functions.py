"""Auxiliary functions for sampling."""

# IMPORTS
## Standard Library Imports
from copy import deepcopy

## Third Party Imports
import pyrosetta

## Local Imports
from ensemblify.utils import df_from_pdb

# FUNCTIONS
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
        contacts (tuple[tuple[tuple[str,tuple[int,int]],tuple[str,tuple[int,int]]],...], optional):
            Residue ranges where two chains are interacting.
    
    Reference:
        See https://docs.rosettacommons.org/demos/latest/tutorials/fold_tree/fold_tree
        for more information about the Rosetta FoldTree.
    """
    # Prepare working pose
    ft_pose = pyrosetta.rosetta.core.pose.Pose()
    ft_pose.detached_copy(pose)

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
