"""Auxiliary functions for creating PyRosetta objects."""

# IMPORTS
## Third Party Imports
import pyrosetta
import pyrosetta.distributed.io as io

# FUNCTIONS
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
            Maximum iterations of the MinMover. Defaults to None, meaning the MinMover object's
            default value.
        dofs (tuple[str,str], optional):
            Defines what angles to set as flexible during minimization. Defaults to backbone and
            sidechain, i.e. ('bb','chi').

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
                                                                          True) # use neighbor list
    if max_iters is not None:
        min_mover.max_iter(max_iters)

    return min_mover
