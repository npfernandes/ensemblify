"""Auxiliary functions for setting up the conformational sampling process."""

# IMPORTS
## Standard Library Imports
import logging
import logging.config
import os
import re
import sys
from copy import deepcopy
from timeit import default_timer as timer

## Third Party Imports
import pandas as pd
import pyrosetta
import pyrosetta.distributed
import ray
import yaml
from tqdm import tqdm

## Local Imports
from ensemblify.generation.ensemble_utils.functions import (
    derive_constraint_targets,
    get_targets_from_plddt,
    _prep_target,
    setup_fold_tree,
)
from ensemblify.generation.ensemble_utils.samplers import setup_samplers
from ensemblify.modelling.constraints import apply_constraints, apply_pae_constraints
from ensemblify.modelling.objects import setup_minmover, setup_pose

# FUNCTIONS
def setup_sampling_logging(sampling_log: str) -> tuple[logging.Logger,str,str]:
    """Setup logging handlers and files for PyRosetta sampling and Ray.
    
    Args:
        sampling_log (str):
            Path to sampling .log file.
    
    Returns:
        tuple[logging.Logger,str,str]:
            logger (logging.Logger): 
                The Logger object associated with the sampling .log file
            ray_log (str):
                Filepath to .log file with Ray log messages.
            pyrosetta_log (str):
                Filepath to .log file with PyRosetta log messages.
    """
    ray_log = os.path.join(os.path.split(sampling_log)[0],'ray.log')
    pyrosetta_log = os.path.join(os.path.split(sampling_log)[0],'pyrosetta.log')
    logging_schema = {
            'version': 1, # Always 1
            'formatters': {
                'default': {
                    'format': '[%(asctime)s] %(name)s - %(process)s - %(levelname)s - %(message)s',
                    'datefmt': '%Y-%m-%d %H:%M:%S'
                },
            },
            # Handlers use the formatter names declared above
            'handlers': {
                'samplingfile': {
                    'class': 'logging.FileHandler',
                    'formatter': 'default',
                    'level': 'INFO',
                    'filename': sampling_log
                },
                'pyrosettafile': {
                    'class': 'logging.FileHandler',
                    'formatter': 'default',
                    'level': 'INFO',
                    'filename': pyrosetta_log
                },
                'rayfile': {
                    'class': 'logging.FileHandler',
                    'formatter': 'default',
                    'level': 'INFO',
                    'filename': ray_log
                },
            },
            # Loggers use the handler names declared above
            'loggers' : {
                '__main__': {
                    'handlers': ['samplingfile'], # Use a list even if one handler is used
                    'level': 'INFO',
                    'propagate': False
                },
                'pyrosetta': {
                    'handlers' : ['pyrosettafile'],
                    'level' : 'INFO',
                    'propagate' : False
                },
                'ray': {
                    'handlers' : ['rayfile'],
                    'level' : 'INFO',
                    'propagate' : False
                },
                },
            # Standalone kwarg for the root logger
            'root' : {
                'handlers': ['samplingfile'],
                'level': 'INFO'
            }
        }
    logging.config.dictConfig(logging_schema)
    logger = logging.getLogger(__name__)
    return logger, ray_log, pyrosetta_log


def _setup_ray_worker_logging():
    logger = logging.getLogger('ray')
    logger.setLevel(logging.ERROR)


def _remove_ansi(file: str):
    """Replace a file with a copy of it without ANSI characters.

    Args:
        file (str):
            Path to text file to change.
    """
    with open(file,'r') as f:
        clean_text = re.sub(r'\033\[[0-9;]*[mGKHF]','', f.read())
    with open(file,'w') as t:
        t.write(clean_text)


def _check_contained_in(
    container: list[int],
    elements: list[list[int]],
) -> list[list[int]]:
    """Check what elements, if any, have all their numbers present in container.

    Args:
        container (list[int]):
            List of integer numbers.
        elements (list[list[int]]):
            List of lists of integer numbers.

    Returns:
        list[list[int]]:
            What lists of integer numbers from elements have all their numbers present
            inside container.
    """
    elements_in_container = []
    for element in elements:
        if all([x in container for x in element]):
            elements_in_container.append(element)
    return elements_in_container


def setup_sampling_parameters(parameters_file: str) -> dict:
    """ 
    Update the parameters dictionary before sampling.
    
    In the 'targets' parameter, change a target from e.g. [1,54] to (range(1,55)).
    If using an AlphaFold model as a starting structure, keep in the sampling ranges only regions
    of at least a certain contiguous size where each residue's pLDTT is below a threshold.
    Targets, secondary structure biases and contacts are also updated to tuples instead of lists.

    Args:
        parameters_file (str):
            Path to parameters file following the Ensemblify template.

    Returns:
        dict:
            The updated params dictionary.
    """
    with open(parameters_file,'r',encoding='utf-8-sig') as document:
        prepared_parameters = yaml.safe_load(document)

    if prepared_parameters['alphafold']:
        # Get list of target residues considering pldTT threshold
        chains_targets = get_targets_from_plddt(parameters=prepared_parameters)

        # Update targets to a continuous range according to pldTT threshold
        updated_targets_all = {}
        for chain in prepared_parameters['targets']:
            updated_targets = []
            for target in prepared_parameters['targets'][chain]:
                sampler,region,database,mode = target
                updated_regions = _check_contained_in(container=[x for x in range(region[0],
                                                                                 region[1]+1)],
                                                     elements=chains_targets[chain])
                for updated_region in updated_regions:
                    updated_targets.append((sampler,updated_region,database,mode))
            updated_targets_all[chain] = tuple(updated_targets)
        prepared_parameters['targets'] = updated_targets_all
    else:
        # Update targets to a continuous range not using pldTT info
        updated_targets_all = {}
        for chain in prepared_parameters['targets']:
            updated_targets = []
            for target in prepared_parameters['targets'][chain]:
                updated_target = deepcopy(target)
                updated_target[1] = tuple(range(target[1][0],target[1][1]+1))
                updated_targets.append(tuple(updated_target))
            updated_targets_all[chain] = tuple(updated_targets)
        prepared_parameters['targets'] = updated_targets_all

    # Change ss_bias content to tuples
    ss_bias = prepared_parameters['restraints']['ss_bias']
    if ss_bias is not None:
        ss_biases = prepared_parameters['restraints']['ss_bias'][0]
        for i,bias in enumerate(ss_biases):
            for j,bias_range in enumerate(bias):
                if isinstance(bias_range,list):
                    prepared_parameters['restraints']['ss_bias'][0][i][j] = tuple(bias_range)
            prepared_parameters['restraints']['ss_bias'][0][i] = tuple(bias)

        prepared_parameters['restraints']['ss_bias'][0] = tuple(ss_bias[0])
        prepared_parameters['restraints']['ss_bias'] = tuple(ss_bias)

    # Change contacts content to tuples
    contacts = prepared_parameters['restraints']['contacts']
    if contacts is not None:
        for i,contact in enumerate(contacts):
            for j, contact_region in enumerate(contact):
                for k,res_range in enumerate(contact_region):
                    if isinstance(res_range,list):
                        prepared_parameters['restraints']['contacts'][i][j][k] = tuple(res_range)
                prepared_parameters['restraints']['contacts'][i][j] = tuple(contact_region)
            prepared_parameters['restraints']['contacts'][i] = tuple(contact)
        prepared_parameters['restraints']['contacts'] = tuple(contacts)

    return prepared_parameters


def setup_sampling_initial_pose(
    params: dict,
    sampling_log: str,
) -> pyrosetta.rosetta.core.pose.Pose:
    """Create initial Pose for sampling and apply any required constraints.

    Constraints applied to the pose are stored in a constraints.cst file stored in the same
    directory as the sampling_log file.

    Args:
        params (dict):
            Parameters following the Ensemblify template.
        sampling_log (str):
            Path to the sampling .log file.
    
    Returns:
        pyrosetta.rosetta.core.pose.Pose:
            Object to be used as the starting structure for the sampling process.
    """
    logger = logging.getLogger(__name__)

    logger.info('Setting up sampling initial structure...')

    # Initialize PyRosetta
    pyrosetta.distributed.init()
    logger.info('PyRosetta was initialized sucessfully.')

    # Create pose
    logger.info('Creating Pose from input structure/sequence...')
    initial_pose = setup_pose(input_structure=params['sequence'],
                              make_centroid=True)

    if params['pae'] != 'None':
        logger.info('Applying pae constraints...')
        apply_pae_constraints(pose=initial_pose,
                              pae_filepath=params['pae'],
                              plddt_targets=params['targets'],
                              cutoff=params['pae_params']['cutoff'],
                              flatten_cutoff=params['pae_params']['flatten_cutoff'],
                              flatten_value=params['pae_params']['flatten_value'],
                              weight=params['pae_params']['weight'],
                              tolerance=params['pae_params']['tolerance'],
                              adjacency_threshold=params['pae_params']['adjacency_threshold'],
                              plddt_scaling_factor=params['pae_params']['plddt_scaling_factor'])

    # Derive constraint targets from sampling targets
    logger.info('Deriving regions to constrain (if any)...')
    c_targets = derive_constraint_targets(pose=initial_pose,
                                          sampling_targets=params['targets'])

    # Apply constraints if needed
    if c_targets:
        logger.info('Applying constraints...')
        apply_constraints(pose=initial_pose,
                          cst_targets=c_targets,
                          contacts=params['restraints']['contacts'],
                          stdev=params['constraints']['stdev'],
                          tolerance=params['constraints']['tolerance'])

        # Adjust fold_tree to improve sampling yield
        logger.info('Setting up FoldTree...')
        setup_fold_tree(pose=initial_pose,
                        constraint_targets=c_targets,
                        contacts=params['restraints']['contacts'])

    # Save applied constraints to file
    cs = initial_pose.constraint_set()
    if cs.has_constraints():
        logger.info('Saving applied constraints to file...')
        pyrosetta.rosetta.core.scoring.constraints.ConstraintIO.write_constraints(
            os.path.join(os.path.split(sampling_log)[0],'constraints.cst'),
            cs,
            initial_pose
            )

    # Save initial pose information to log file
    logger.info(f'Initial Pose:\n{initial_pose}')
    if c_targets:
        logger.info(f'Intra-chain constrained regions (Pose numbering):\n{c_targets}')
    logger.info('Initial structure has been setup sucessfully.')

    return initial_pose


def _get_dbs_mem_size(databases: dict[str,pd.DataFrame]) -> int:
    """Get a rough estimate of the size of a databases dictionary, in bytes.

    Store the size of the databases dictionary and for each database, the database ID,
    the database dictionary, the aminoacid ID and DataFrame.

    Args:
        databases (dict[str,pd.DataFrame]):
            Mapping of database IDs to Ensemblify databases.

    Returns:
        int:
            Total memory size of each elements of the databases, in bytes.
    """
    total_mem = sys.getsizeof(databases)
    for db_id,db in databases.items():
        total_mem += sys.getsizeof(db_id) + sys.getsizeof(db)
        for aa,df in db.items():
            total_mem += sys.getsizeof(aa) + df.memory_usage(index=True,deep=True).sum()
    return total_mem


@ray.remote(num_returns=1,
            num_cpus=1)
def sample_pdb(
    ppose: pyrosetta.distributed.packed_pose.core.PackedPose,
    databases: dict[str,dict[str,pd.DataFrame]],
    targets: dict[str,tuple[tuple[str,tuple[int,...],str,str]]],
    output_path: str,
    job_name: str,
    decoy_num: str = '',
    log_file: str = None,
    ss_bias: tuple[tuple[tuple[str,tuple[int,int],str],...],int] | None = None,
    variance: float | None = 0.10,
    sampler_params: dict[str,dict[str,int]] = {'MC':{'temperature':200,'max_loops':200}},
    scorefxn_id: str = 'score0',
    scorefxn_weight: float = 1.0,
    minimizer_id: str = 'dfpmin_armijo_nonmonotone',
    minimizer_tolerance: float = 0.001,
    minimizer_maxiters: int = 5000,
    minimizer_finalcycles: int = 5,
    cst_weight: int = 1,
    cstviolation_threshold: float = 0.015,
    cstviolation_maxres: int = 20,
    ) -> str | None:
    """Sample dihedral angles from a database into target regions of a given structure.

    Args:
        ppose (pyrosetta.distributed.packed_pose.core.PackedPose):
            Reference to the initial structure.
        databases (dict[str,dict[str,pd.DataFrame]]):
            Reference to the databases dictionary.
        targets (dict[str,tuple[tuple[str,tuple[int,...],str,str]]]):
            Dictionary detailing the target regions for sampling in each chain.
        output_path (str):
            Path to directory where sampled structures will be written to.
        job_name (str):
            Prefix identifier for generated structures.
        decoy_num (str, optional):
            Identifier to differentiate between different decoys of the same batch in a
            multiprocessing context. Defaults to ''.
        log_file (str, optional):
            Path to the PyRosetta .log file. Defaults to 'pyrosetta.log' in current working
            directory.
        ss_bias (tuple[tuple[tuple[str,tuple[int,int],str],...],int], optional):
            Secondary Structure Bias with the desired percentage of total structures to respect
            this bias. Defaults to None.
        variance (float, optional):
            New dihedral angle values inserted into sampling regions are sampled from a Gaussian
            distribution centered on the value found in database and percentage variance equal to
            this value. Defaults to 0.10 (10%).
        sampler_params (dict[str,dict[str,int]], optional):
            Parameters for the used sampler, assumes MonteCarloSampler is used. Defaults to
            {'MC':{'temperature':200,'max_loops':200}}.
        scorefxn_id (str, optional):
            PyRosetta ScoreFunction identifier. Must pertain to a .wst weights file present in
            /.../pyrosetta/database/scoring/weights/ . Defaults to 'score0'.
        scorefxn_weight (float, optional):
            Weight for the repulsive Van der Waals term in the ScoreFunction. Will only have an
            effect if the ScoreFunction has a repulsive Van der Waals term. Defaults to 1.0.
        minimizer_id (str, optional):
            PyRosetta minimization algorithm identifier used in MinMover.
            Defaults to 'dfpmin_armijo_nonmonotone'.
        minimizer_tolerance (float, optional):
            Tolerance value for the PyRosetta MinMover object. Defaults to 0.001.
        minimizer_maxiters (int, optional):
            Maximum iterations value for the PyRosetta MinMover object. Defaults to 5000.
        minimizer_finalcycles (int, optional):
            Number of times to apply the MinMover to our final structure. Defaults to 5.
        cst_weight (int, optional):
            Weight of the AtomPairConstraint term in the ScoreFunction. Defaults to 1.
        cstviolation_threshold (float, optional):
            Any residue with AtomPairConstraint score term value above this threshold is considered
            in violation of the applied constraints. Defaults to 0.015.
        cstviolation_maxres (int, optional):
            Number of residues allowed to be above the constraint violation threshold.
            Defaults to 20.

    Returns:
        str | None:
            Path to the sampled .pdb structure. Only written and returned if the
            sampled structure is valid (does not violate the applied constraints).
            Otherwise, return None.
    """
    import pyrosetta
    import pyrosetta.distributed
    import pyrosetta.distributed.io as io

    pyrosetta.distributed.maybe_init()

    # Setup log file
    if log_file is None:
        log_file = os.path.join(os.getcwd(),'pyrosetta.log')

    # Setup score function
    scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_id)

    # Set VdW weight in score function
    vdw_cst = pyrosetta.rosetta.core.scoring.ScoreType.vdw
    scorefxn.set_weight(vdw_cst,scorefxn_weight)

    # Set AtomPairConstraint weight in score function
    ap_cst = pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint
    scorefxn.set_weight(ap_cst,cst_weight)

    # Setup minimization mover
    min_mover = setup_minmover(scorefxn=scorefxn,
                               min_id=minimizer_id,
                               tolerance=minimizer_tolerance,
                               max_iters=minimizer_maxiters)
    # Setup sampler
    samplers = setup_samplers(sampler_params=sampler_params,
                              variance=variance,
                              scorefxn=scorefxn,
                              databases=databases,
                              log_file=log_file)

    # Grab initial pose (already setup with constraints)
    initial_pose = io.to_pose(ppose)

    # Setup working pose
    working_pose = pyrosetta.rosetta.core.pose.Pose()
    working_pose.detached_copy(initial_pose)

    # Do sampling
    tqdm.write(f'[DECOY {decoy_num}] Sampling ...')
    start = timer()
    for chain in targets:
        for sampler,target_dirty,database_id,sampling_mode in targets[chain]:

            # Fix residue numbering for multiple chain or negative resnum situations
            target = [ pyrosetta.rosetta.core.pose.pdb_to_pose(working_pose, res_num, chain)
                       for res_num in target_dirty ]

            assert sampler == 'MC', 'Sampler ID must be \'MC\'!'

            if working_pose.num_chains() > 1:
                chain_end = pyrosetta.rosetta.core.pose.chain_end_res(
                        working_pose,
                        pyrosetta.rosetta.core.pose.get_chain_id_from_chain(chain,
                                                                            working_pose))
            else:
                chain_end = working_pose.size()

            # So we don't sample fragments with length 2
            target_prep = _prep_target(seq_len=chain_end,
                                       target=target)

            if ss_bias is not None:
                ssb = ss_bias[0]
            else:
                ssb = None

            samplers['MC'].apply(pose=working_pose,
                                 target=target_prep,
                                 chain=chain,
                                 database_id=database_id,
                                 ss_bias=ssb,
                                 sampling_mode=sampling_mode)

            # add if clauses and logic for future samplers here
    end = timer()
    tqdm.write(f'[DECOY {decoy_num}] Sampling time: {end-start}')

    # Perform energy minimization
    tqdm.write(f'[DECOY {decoy_num}] Minimizing ...')
    start = timer()
    for _ in range(minimizer_finalcycles):
        pose_score = scorefxn(working_pose)
        if float(pose_score) < 0:
            break
        min_mover.apply(working_pose)
    end = timer()
    tqdm.write(f'[DECOY {decoy_num}] Minimization time: {end-start}')

    # Check if constraints are respected
    residue_total_energies = working_pose.energies().residue_total_energies_array()
    respectful = True
    counter = 0
    for e in residue_total_energies['atom_pair_constraint']:
        if e > cstviolation_threshold:
            counter += 1
    if counter > cstviolation_maxres:
        respectful = False

    # If constraints are respected then output pose as .pdb
    if respectful:
        # Setup pdb file name
        output_filename = os.path.join(output_path,f'{decoy_num}_{job_name}.pdb')

        # Output pose
        io.to_pose(working_pose).dump_pdb(output_filename)
        return output_filename

    return None
