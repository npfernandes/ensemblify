"""Auxiliary functions for sampling."""

# IMPORTS
## Standard Library Imports
from collections import deque  
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
 
    Upate the given Pose's FoldTree to minimize "lever arm" effects that might result in movement
    in constrained regions during sampling. The goal is to try to maximmize the amount of
    constrained residues that are 'upstream' in relation to sampled residues.

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

    num_chains = ft_pose.num_chains()

    # Calculate optimal central residue per chain (constrained res closest to chain midpoint)
    central_residues = {}
    for i in range(1, num_chains + 1):
        chain_start = ft_pose.chain_begin(i)
        chain_end = ft_pose.chain_end(i)
        ideal_mid = (chain_start + chain_end) // 2
        best_res, min_dist = ideal_mid, float('inf')
        for first, last in constraint_targets:
            for res in range(first, last + 1):
                if ft_pose.chain(res) == i:
                    dist = abs(res - ideal_mid) # distance between residue numbers
                    if dist < min_dist:
                        min_dist, best_res = dist, res
        central_residues[i] = best_res

    # Build inter-chain contact graph
    # contact_graph[(c_lo, c_hi)] = (res_in_c_lo, res_in_c_hi)
    contact_graph = {}
    if contacts is not None:
        for contact in contacts:
            chain_x_id, (x_first, x_last) = contact[0]
            chain_y_id, (y_first, y_last) = contact[1]
            x_mid = (x_first + x_last) // 2
            y_mid = (y_first + y_last) // 2
            x_res = pyrosetta.rosetta.core.pose.pdb_to_pose(ft_pose, x_mid, chain_x_id)  
            y_res = pyrosetta.rosetta.core.pose.pdb_to_pose(ft_pose, y_mid, chain_y_id)
            cx, cy = ft_pose.chain(x_res), ft_pose.chain(y_res)
            if cx == cy:  
                continue
            key = (min(cx, cy), max(cx, cy))
            if key not in contact_graph:
                contact_graph[key] = (x_res, y_res) if cx < cy else (y_res, x_res)

    # Root chain = chain with most sampled residues
    sampled_count = {i: 0 for i in range(1, num_chains + 1)}
    for lo, hi in constraint_targets:  
        for res in range(lo, hi + 1):
            sampled_count[ft_pose.chain(res)] += 1
    root_chain = max(sampled_count, key=sampled_count.get)

    # Breadth first search spanning tree of contact graph from root_chain
    adjacency = {i: [] for i in range(1, num_chains + 1)}
    for (ci, cj), (ri, rj) in contact_graph.items():
        adjacency[ci].append((cj, ri, rj))  
        adjacency[cj].append((ci, rj, ri))
    visited = {root_chain}
    queue = deque([root_chain])
    spanning_edges = [] # (upstream_res, downstream_res)
    while queue:
        current = queue.popleft()
        for neighbor, res_in_current, res_in_neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(neighbor) 
                queue.append(neighbor)
                spanning_edges.append((res_in_current, res_in_neighbor))

    # Connect isolated chains (no contacts) directly to root
    for chain_num in range(1, num_chains + 1):
        if chain_num not in visited:
            spanning_edges.append((central_residues[root_chain], central_residues[chain_num]))
            visited.add(chain_num)

    # Assign cutpoints (chain boundaries between jump residues)
    # Process edges with fewest valid options first to avoid conflicts.
    def valid_cutpoints(r1, r2, used):
        first, last = min(r1, r2), max(r1, r2)
        return [ft_pose.chain_end(c) for c in range(1, num_chains)
                if first < ft_pose.chain_end(c) < last and ft_pose.chain_end(c) not in used]

    used_cuts = set()
    # Sort by number of valid options (most constrained first)
    ordered = sorted(spanning_edges, key=lambda e: len(valid_cutpoints(e[0], e[1], set())))
    jump_specs = []  # (r1, r2, cutpoint)
    for r1, r2 in ordered:
        options = valid_cutpoints(r1, r2, used_cuts)
        if not options:
            raise ValueError(
                f'No valid cutpoint between residues {r1} (chain {ft_pose.chain(r1)}) and {r2} '
                f'(chain {ft_pose.chain(r2)}). Please check that contacts span distinct chains '
                'and chain order is correct.'
            )
        cut = options[0]
        used_cuts.add(cut)
        jump_specs.append((r1, r2, cut))

    # Build new FoldTree from scratch
    ft = pyrosetta.rosetta.core.kinematics.FoldTree()
    ft.simple_tree(ft_pose.size())
    for r1, r2, cut in jump_specs:
        ft.new_jump(r1, r2, cut)

    # Root at the sampled chain's central residue, upstream of sampled regions
    ft.reorder(central_residues[root_chain])

    # Check validity before continuing
    assert ft.check_fold_tree(), f'Invalid FoldTree setup:\n{ft}'

    # Finally, assign new FoldTree and exit
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
