"""Auxiliary functions for combining PDB structures with disordered tails/linkers."""

# IMPORTS
## Standard Library Imports
import os
import sys

## Third Party Imports
import pyrosetta

## Local Imports
from ensemblify.modelling._alignment import (
    setup_sequences_structures,
    perform_global_alignment,
    parse_alignment,
    extract_indels_mutations
)
from ensemblify.modelling.objects import setup_pose, setup_minmover
from ensemblify.modelling.constraints import apply_constraints
from ensemblify.modelling.pdb_processing import process_pdb_structure

# FUNCTIONS
def apply_basic_idealize(
    pose: pyrosetta.rosetta.core.pose.Pose,
    start_res: int,
    end_res: int,
    sfxn: pyrosetta.rosetta.core.scoring.ScoreFunction,
    fast: bool = False):
    """Idealize bond geometry of a range of residues in a Pose object.
    
    Uses the basic_idealize protocol from PyRosetta. The provided Pose object is modified in place.
    In short, the Mover stochastically picks a move-able position (between the provived start and
    end points), forces that position into ideal bond geometry, and tries to use minimization to
    bring the coordinates back to very near their starting points.

    Args:
        pose (pyrosetta.rosetta.core.pose.Pose):
            The Pose object to be idealized.
        start_res (int):
            The starting residue index (1-indexed) for the idealization range.
        end_res (int):
            The ending residue index (1-indexed) for the idealization range.
        sfxn (pyrosetta.rosetta.scoring.ScoreFunction):
            The ScoreFunction to be used for the idealization.
        fast (bool):
            If True, uses a faster version of the idealization protocol.
            Defaults to False, which uses the standard idealization protocol.
    """
    # Initialize required Rosetta data structure
    idealize_vector = pyrosetta.rosetta.utility.vector1_unsigned_long()

    # Add start and end residues to the idealization vector
    idealize_vector.append(start_res)
    idealize_vector.append(end_res)

    # Idealize the specified range of residues in the Pose object
    pyrosetta.rosetta.protocols.idealize.basic_idealize(pose,
                                                        idealize_vector,
                                                        sfxn,
                                                        int(fast))


def pack_sidechains(
    score_function: pyrosetta.rosetta.core.scoring.ScoreFunction,
    pose: pyrosetta.rosetta.core.pose.Pose):
    """Add side-chains to a Pose using a PackRotamersMover with the given score function.
    
    Changes provided Pose object in place.

    Args:
        score_function (pyrosetta.rosetta.core.scoring.ScoreFunction):
            The score function to use for packing.
        pose (pyrosetta.rosetta.core.pose.Pose):
            The Pose object to pack side-chains into.
    """
    # Initialize PyRosetta standard packer task
    packer_task = pyrosetta.standard_packer_task(pose)
    packer_task.restrict_to_repacking()
    
    # Initialize PackRotamersMover with ScoreFunction and pack side-chains
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(score_function,
                                                                                packer_task)
    packer.apply(pose)


def _setup_fusion_objects(
    sequences: list[str],
    pdb_structures: list[str],
    ) -> list[tuple[str, tuple[int, int]]]:
    """Setup fusion objects by aligning sequences with PDB structures and updating the list of
    sequences with the corresponding PDB filepaths where applicable.

    Args:
        sequences (list[str]):
            List of sequences to be fused, where each sequence corresponds to a segment of the
            full-length structure.
        pdb_structures (list[str]):
            List of .pdb file paths corresponding to the folded domains in the full length
            structure.

    Returns:
        list[tuple[str, tuple[int, int]]]:
            List of fusion objects, where sequences corresponding to the provided PDB structures
            are replaced with their corresponding PDB filepaths and their residue number range in
            the full length sequence.
    """
    # Setup folded sequences for reference
    sequences2pdb = { setup_pose(pdb).sequence() : pdb for pdb in pdb_structures }

    # Initialize fusion objects list with sequences
    fusion_objects = [x for x in sequences]

    # Update fusion objects with corresponding PDB structure filepaths
    fl_sequence = ''.join(sequences)
    
    for sequence_idx,sequence in enumerate(sequences):
        # Perform the alignment
        alignment = perform_global_alignment(target=fl_sequence,
                                             query=sequence)

        # Parse the alignment
        aligned_target, aligned_query, alignment_encoding = parse_alignment(alignment=alignment)
        
        # Extract inserts, deletions and mutations from aligned sequences
        indels, _ = extract_indels_mutations(aligned_target=aligned_target,
                                             aligned_query=aligned_query,
                                             alignment_encoding=alignment_encoding)
        
        # Interpret indels, update fusion objects and register residue numbers
        ## Here we will assume only insertions are present (the expected case)

        # If there is only one insertion, our query sequence is one of the terminal regions
        if len(indels) == 1:
            indel = indels[0]
            if indel[1] == indel[2] == 1:
                # If insert is on 1st position of query, our query is the C-terminal region
                res_range = (len(indel[0]) + 1, len(indel[0]) + len(sequence))
            elif indel[1] == indel[2] == len(sequence):
                # If insert is on last position of query, our query is the N-terminal region
                res_range = (1, len(sequence))
        elif len(indels) == 2:
            # If there are two indels, we know our query sequence is a linker
            res_range = (len(indels[0][0]) + 1,
                         len(indels[0][0]) + len(sequence))

        # If current sequence is a folded domain, replace it with the corresponding PDB file
        # and register its residue number range in the full length sequence
        if sequence in sequences2pdb:
            fusion_objects[sequence_idx] = (sequences2pdb[sequence], res_range)
        else:
            # If current sequence is a disordered region, simply add its residue number range
            # in the full length sequence
            fusion_objects[sequence_idx] = (sequence, res_range)

    return fusion_objects


def _setup_first_fusion_piece(
    fusion_object: tuple[str, tuple[int, int]],
    scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction,
    minmover: pyrosetta.rosetta.protocols.minimization_packing.MinMover,
    ) -> pyrosetta.rosetta.core.pose.Pose:
    """
    Initialize the first fusion piece, which is the N-terminal fusion object.
    
    If the fusion object is a PDB file, it will be loaded as a Pose.
    If it is a sequence, a Pose will be created from that sequence, side-chains are repacked
    and the structure is minimized.
    
    Args:
        fusion_object (tuple[str, tuple[int, int]]):
            A tuple containing the fusion object as a PDB file or sequence string and its residue
            number range in the full length sequence.
        scorefxn (pyrosetta.rosetta.core.scoring.ScoreFunction):
            The score function for packing side-chains, if applicable.
        minmover (pyrosetta.rosetta.protocols.minimization_packing.MinMover):
            The MinMover for minimization, if applicable.

    Returns:
        pyrosetta.rosetta.core.pose.Pose:
            The initialized Pose object for the first fusion piece.
    """
    # Initialize Pose from the first fusion object
    fusion_piece = setup_pose(fusion_object[0])

    # If its a N-terminal tail, add side-chains and minimize
    if not fusion_object[0].endswith('.pdb'):
        # Add sidechains to the linker structure
        pack_sidechains(score_function=scorefxn,
                        pose=fusion_piece)
            
        # Minimize the linker structure
        minmover.apply(fusion_piece)

    return fusion_piece


def _set_active_dof(
    movemap: pyrosetta.rosetta.core.kinematics.MoveMap,
    cst_targets: tuple[tuple[int, int], ...],
    ):
    """Set the active degrees of freedom (dof) in a MoveMap.

    Changes the provided MoveMap object in place.

    Args:
        movemap (pyrosetta.rosetta.core.kinematics.MoveMap):
            The MoveMap object to modify.
        cst_targets (tuple[tuple[int, int], ...]):
            A tuple of tuples, where each inner tuple contains the start and end residue indices
            for the folded regions in the Pose that will not be considered flexible in the MoveMap.
    """
    for start_res, end_res in cst_targets:
        # Set backbone and chi angles for residues in the folded regions
        for i in range(start_res, end_res + 1):
            movemap.set_bb(i, True)
            movemap.set_chi(i, True)


def create_fusion_pose(
    sequences: list[str],
    pdb_structures: list[str],
    output_name: str,
    output_dir: str | None = None,
    fusion_verbose: bool = True,
    pyrosetta_verbose: bool = False,
    ) -> pyrosetta.rosetta.core.pose.Pose:
    """Create a full-length fused Pose object from provided sequences and PDB structures.

    Args:
        sequences (list[str]):
            List of sequences to be fused, where each sequence corresponds to a segment of the
            full-length structure.
        pdb_structures (list[str]):
            List of .pdb file paths corresponding to the folded domains in the full length
            structure.
        output_name (str):
            Name for the output fused .pdb file (without extension).
        output_dir (str | None):
            Directory where the output .pdb file will be saved. Defaults to current working
            directory.
        fusion_verbose (bool):
            Whether to print verbose output during the fusion process. Defaults to True.
        pyrosetta_verbose (bool):
            Whether to initialize PyRosetta with verbose output. Defaults to False.

    Returns:
        pyrosetta.rosetta.core.pose.Pose:
            The final fused Pose object containing all sequences and PDB structures.
    """
    # Setup output directory
    if output_dir is None:
        output_dir = os.getcwd()
    elif not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Initialize PyRosetta
    if pyrosetta_verbose:
        pyrosetta.init()
    else:
        pyrosetta.init(options='-mute all',
                       silent=True)

    # Setup PyRosetta QoL variables and objects
    rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue
    rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue
    sff = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory

    ## ScoreFunction for packing the side-chains of the inserted linker segments
    sfxn_sc = sff.create_score_function('ref2015_cart')

    ## ScoreFunction for minimization steps applied to fusion sites
    sfxn_idealize = sff.create_score_function('ref2015')
    
    ## MinMover object for linker minimization
    minmover = setup_minmover(scorefxn=sfxn_sc,
                              min_id='lbfgs_armijo_nonmonotone',
                              tolerance=0.0001)
    
    # Initialize final Pose object that will hold the full-length structure
    FLpose = pyrosetta.rosetta.core.pose.Pose()

    # Update sequences list with corresponding PDB structures (where applicable)
    fusion_objects = _setup_fusion_objects(sequences=sequences,
                                           pdb_structures=pdb_structures)

    if fusion_verbose:
        print('Fusing provided PDB structures and sequences into a full-length Pose...')

    # Initialize our first fusion piece, the N-terminal fusion object
    if fusion_verbose:
        print(f'Setting up {fusion_objects[0][0]} as the first fusion piece (N-terminal)...')
    Npose = _setup_first_fusion_piece(fusion_object=fusion_objects[0],
                                      scorefxn=sfxn_sc,
                                      minmover=minmover)
    
    # Setup fusion tasks based on remaining fusion objects
    ## ( Object to add into current Pose, Type of object to add {Linker,Domain} )
    fusion_tasks = [
        (obj[0], setup_pose(obj[0]), 'Folded Domain') if obj[0].endswith('.pdb')
        else (obj[0], setup_pose(obj[0]),'Disordered Region') for obj in fusion_objects[1:]
    ]

    # Perform fusion tasks sequentially   
    for fusion_object, fusion_pose, fusion_type in fusion_tasks:
        if fusion_verbose:
            print(f'Fusing {fusion_object} as a {fusion_type}...')
            sys.stdout.flush()

        # Initialize fusion object to add to Npose C-terminal
        Cpose = fusion_pose

        # Define the start and end residues for the fusion site,
        # whose bond geometry will later be idealized
        start_res = Npose.total_residue()
        if fusion_type == 'Disordered Region':
            # If we are fusing a linker, all of it is idealized
            end_res = start_res + Cpose.total_residue()
        else:
            # If adding folded domain, only the residues in the fusion site are idealized
            end_res = start_res + 1

        # Remove terminal tag from C-terminal residue of Npose and
        # N-terminal residue of Cpose (the two that will be attached)
        rm_upper(Npose.conformation(), start_res)
        rm_lower(Cpose.conformation(), 1)

        if fusion_type == 'Disordered Region':
            # Add sidechains to the linker structure
            pack_sidechains(score_function=sfxn_sc,
                            pose=Cpose)
            
            # Minimize the linker structure
            minmover.apply(Cpose)

        # Update Npose by adding Cpose to its C-terminal
        Npose = pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(Npose,
                                                                           Cpose,
                                                                           Npose.total_residue(),
                                                                           Npose.total_residue())

        # Update Npose FoldTree after fusion
        ft = pyrosetta.rosetta.core.kinematics.FoldTree()
        ft.simple_tree(Npose.total_residue())
        Npose.fold_tree(ft)

        # Idealize the structure of the residues making up the fusion site
        apply_basic_idealize(pose=Npose,
                             start_res=start_res,
                             end_res=end_res,
                             sfxn=sfxn_idealize)
        
        # Finally, update the FLPose with the fused NPose+Cpose
        FLpose.assign(Npose)

    if fusion_verbose:
        print('Relaxing flexible regions of the fused Pose...')
        sys.stdout.flush()

    # Setup residue numbers of folded regions (constraint targets) based on full length sequence
    cst_targets = tuple([fusion_object[1] for fusion_object in fusion_objects
                         if fusion_object[0].endswith('.pdb')])

    # Add constraints to folded domains of FLPose before FastRelax of flexible regions
    apply_constraints(pose=FLpose,
                      cst_targets=cst_targets)

    # Add weight to AtomPairConstraint in the ScoreFunction to be used in FastRelax
    sfxn_idealize.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint, 10.0)

    # Setup MoveMap for FastRelax
    movemap = pyrosetta.rosetta.core.kinematics.MoveMap()
    _set_active_dof(movemap=movemap,
                    cst_targets=cst_targets)

    # Setup FastRelax protocol using ScoreFunction and MoveMap
    fastrelax = pyrosetta.rosetta.protocols.relax.FastRelax()
    fastrelax.min_type('lbfgs_armijo_nonmonotone')
    fastrelax.set_scorefxn(sfxn_idealize)
    fastrelax.set_movemap(movemap)

    # Perform FastRelax on flexible regions of the fused Pose
    fastrelax.apply(FLpose)
    
    # Dump the final relaxed Pose object to a PDB file
    print('Pose fusion complete.')
    FLpose.dump_pdb(f'{os.path.join(output_dir,output_name)}.pdb')

    return FLpose


def fuse_structures(
    input_fastas: list[str],
    input_pdbs: list[str],
    output_name: str,
    output_dir: str | None = None,
    ) -> tuple[str,str]:
    """Fuse sequences and PDB structures into a full-length PDB structure.

    Args:
        input_fastas (list[str]):
            Path(s) to FASTA file(s) containing the sequences of all protein domains
            (folded + disordered) to be fused, in order from N- to C-terminal.
        input_pdbs (list[str]):
            Path(s) to PDB file(s) containing the structures of folded protein domains
            to be fused, in order from N- to C-terminal.
        output_name (str):
            Name for the output fused PDB file (without extension).
        output_dir (str):
            Directory where the output .pdb file will be saved. Defaults to current working
            directory.
    Returns:
        tuple[str,str]:
            str:
                Path to the steric clash report from applying Pulchra to the fused PDB structure.
            str:
                Path to the fused PDB structure.
    """
    # Setup output directory
    if output_dir is None:
        output_dir = os.getcwd()
    elif not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Setup sequences and PDB structures from input files
    sequences, \
    pdb_structures = setup_sequences_structures(input_fastas=input_fastas,
                                                input_pdbs=input_pdbs)

    # Create the fused Pose object from sequences and PDB structures
    ## This step outputs a .pdb file with the fused structure
    _ = create_fusion_pose(sequences=sequences,
                           pdb_structures=pdb_structures,
                           output_name=output_name,
                           output_dir=output_dir)
    
    # Define the path for the output fused PDB file
    fused_pdb_path = f'{os.path.join(output_dir,output_name)}.pdb'

    # Process the fused PDB to optimize its structure and avoid steric clashes
    print('Processing fused Pose...')
    processed_fused_pdb_log, \
    processed_fused_pdb = process_pdb_structure(pdb=fused_pdb_path)
    print('Processing complete.\nYour full-length structure is ready.')

    return processed_fused_pdb_log, processed_fused_pdb
