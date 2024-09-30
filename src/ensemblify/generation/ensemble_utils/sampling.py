"""Function to sampling conformational ensembles from an input structure."""

# IMPORTS
## Standard Library Imports
import glob
import os

## Third Party Imports
import pyrosetta.distributed.io as io
import ray
from tqdm import tqdm

## Local Imports
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.generation.ensemble_utils.inputs_processing import register_input_clashes
from ensemblify.utils import HashableDict
from ensemblify.generation.ensemble_utils.sampling_utils import (setup_sampling_logging,
                                                                   setup_ray_worker_logging,
                                                                   setup_sampling_parameters,
                                                                   setup_sampling_initial_pose,
                                                                   sample_pdb)
from ensemblify.generation.ensemble_utils.movers_utils import setup_databases
from ensemblify.generation.ensemble_utils.pdb_processing import process_pdb

# FUNCTIONS
def run_sampling(
    input_parameters: str,
    input_clashes_file: str | None,
    valid_pdbs_dir: str,
    sampling_log: str):
    """Perform conformational sampling according to input parameters.

    The Pulchra clashes file for the input structure must be provided, along with the path to
    the directory where sampled structures will be stored and a sampling log file.
    Additional log files will be created in the same directory as the provided log file.

    Args:
        input_parameters:
            path to parameters file following the Ensemblify template.
        input_clashes_file:
            path to Pulchra log file for the input structure.
        valid_pdbs_dir:
            path to the directory where sampled structures will be stored.
        sampling_log:
            path to the .log file for the sampling process.
    """
    # Setup logging
    logger, RAY_LOG, PYROSETTA_LOG = setup_sampling_logging(sampling_log=sampling_log)

    # Setup parameters file
    PARAMETERS = setup_sampling_parameters(parameters_file=input_parameters)

    # Save setup sampling parameters
    sampling_params_path = os.path.join(os.path.split(input_parameters)[0],
                                        '3_parameters_sampling.yaml')
    with open(sampling_params_path,'w',encoding='utf-8') as f:
        for key in PARAMETERS:
            f.write(f'{key}:\n')
            if isinstance(PARAMETERS[key],dict):
                for nested_key in PARAMETERS[key]:
                    f.write(f'  {nested_key}: ')
                    if PARAMETERS[key][nested_key] is not None:
                        f.write(f'{PARAMETERS[key][nested_key]}\n')
                    else:
                        f.write('\n')
            else:
                f.write(f'  {PARAMETERS[key]}\n')

    # Setup input objects
    clashes_input = register_input_clashes(input_clashes_file=input_clashes_file)
    initial_pose =  setup_sampling_initial_pose(params=PARAMETERS,
                                                sampling_log=sampling_log)
    databases = HashableDict(setup_databases(PARAMETERS['databases']))

    # Setup sampling constants and variables
    logger.info('Setting up sampling...')

    ## Assign constants
    VALID_PDBS_DIR = valid_pdbs_dir
    SAMPLING_TARGETS = HashableDict(PARAMETERS['targets'])
    OUTPUT_PATH = PARAMETERS['output_path']

    if PARAMETERS['faspr_path'] is None:
        FASPR_PATH = GLOBAL_CONFIG['FASPR_PATH']
    else:
        FASPR_PATH = PARAMETERS['faspr_path']

    if PARAMETERS['pulchra_path'] is None:
        PULCHRA_PATH = GLOBAL_CONFIG['PULCHRA_PATH']
    else:
        PULCHRA_PATH = PARAMETERS['pulchra_path']

    JOB_NAME = PARAMETERS['job_name']
    SS_BIAS = PARAMETERS['restraints']['ss_bias']
    SAMPLER_PARAMS = HashableDict(PARAMETERS['sampler_params'])
    SCOREFXN_ID = PARAMETERS['scorefxn']['id']
    SCOREFXN_WEIGHT = PARAMETERS['scorefxn']['weight']
    MINIMIZER_ID = PARAMETERS['minimizer']['id']
    MINIMIZER_TOL = PARAMETERS['minimizer']['tolerance']
    MINMIZER_MAXITERS = PARAMETERS['minimizer']['max_iters']
    MINMIZER_FINALCYCLES = PARAMETERS['minimizer']['finalcycles']
    CONSTRAINTS_WEIGHT = PARAMETERS['constraints']['weight']
    CONSTRAINTS_VIOL_THRESHOLD = PARAMETERS['constraints_violation']['threshold']
    CONSTRAINTS_VIOL_MAXRES = PARAMETERS['constraints_violation']['maxres']

    ## Setup sampling cutpoints
    GOAL_ENSEMBLE_SIZE = int(PARAMETERS['size'])
    if SS_BIAS is not None:
        SS_BIAS_CUTPOINT = round(GOAL_ENSEMBLE_SIZE * (SS_BIAS[1] / 100))
    else:
        SS_BIAS_CUTPOINT = None

    # Setup progress bar
    pbar = tqdm(total=GOAL_ENSEMBLE_SIZE,
                desc='Ensemblifying... ',
                unit='valid_pdb',
                miniters=1,
                dynamic_ncols=True)

    logger.info('Sampling has been setup successfully.')

    # Setup Ray
    ray.init(num_cpus=PARAMETERS['core_amount'],
             object_store_memory=1000000000,
             log_to_driver=False,
             runtime_env={'worker_process_setup_hook': setup_ray_worker_logging})

    # Put input clashes information in Ray object store memory
    if clashes_input != []:
        stored_input_clashes = ray.put(clashes_input)
    else:
        stored_input_clashes = None

    # Put databases in Ray object store memory
    stored_databases = ray.put(databases)

    # Put packed initial pose in Ray object store memory
    stored_initial_pose = ray.put(io.to_packed(initial_pose))

    # Setup counters
    disrespectful_pdbs = 0
    clashed_pdbs = 0
    valid_pdbs = 0

    # Setup sampling aux variables
    unfinished_obj_refs_batch = []
    MAX_PENDING_TASKS = PARAMETERS['core_amount'] # Max number of tasks to allow in the system at once
    decoy_num = 0

    # Sample with secondary structure bias (if applicable)
    if SS_BIAS_CUTPOINT is not None:
        logger.info('Starting sampling with secondary structure bias...')
        while valid_pdbs < SS_BIAS_CUTPOINT or len(os.listdir(VALID_PDBS_DIR)) < SS_BIAS_CUTPOINT:

            # Submit new tasks if we haven't reached the target count
            while len(unfinished_obj_refs_batch) < MAX_PENDING_TASKS:
                sampling_future = sample_pdb.remote(stored_initial_pose,
                                                    stored_databases,
                                                    SAMPLING_TARGETS,
                                                    OUTPUT_PATH,
                                                    JOB_NAME,
                                                    str(decoy_num),
                                                    PYROSETTA_LOG,
                                                    SS_BIAS,
                                                    SAMPLER_PARAMS,
                                                    SCOREFXN_ID,
                                                    SCOREFXN_WEIGHT,
                                                    MINIMIZER_ID,
                                                    MINIMIZER_TOL,
                                                    MINMIZER_MAXITERS,
                                                    MINMIZER_FINALCYCLES,
                                                    CONSTRAINTS_WEIGHT,
                                                    CONSTRAINTS_VIOL_THRESHOLD,
                                                    CONSTRAINTS_VIOL_MAXRES,)

                processing_future = process_pdb.remote(sampling_future,
                                                       FASPR_PATH,
                                                       PULCHRA_PATH,
                                                       stored_input_clashes,
                                                       SAMPLING_TARGETS,
                                                       VALID_PDBS_DIR,
                                                       SS_BIAS_CUTPOINT,)

                unfinished_obj_refs_batch.append(processing_future)
                decoy_num += 1

            # Returns the ObjectRefs that are ready
            finished_tasks, unfinished_obj_refs_batch = ray.wait(unfinished_obj_refs_batch,
                                                                num_returns=1)

            for finished_task in finished_tasks:
                try:
                    result = ray.get(finished_task)
                except ray.exceptions.RayTaskError as e:
                    with open(RAY_LOG,'a',encoding='utf-8') as f:
                        f.write(str(e))
                    continue
                else:
                    # Process result
                    if isinstance(result,type(None)):
                        disrespectful_pdbs += 1
                    elif isinstance(result,bool):
                        if result:
                            clashed_pdbs += 1
                        else:
                            valid_pdbs += 1
                            if pbar.n < SS_BIAS_CUTPOINT:
                                pbar.update(1)

    # If we wanted 100% with secondary structure bias, exit
    if len(os.listdir(VALID_PDBS_DIR)) == SS_BIAS_CUTPOINT == GOAL_ENSEMBLE_SIZE:
        logger.info('Desired number of structurally biased structures has been reached.')
        logger.info('Cancelling any remainings tasks and shutting down ray client...')
        ray.shutdown()

        # Close progress bar (required)
        pbar.set_description('Ensemblified!')
        pbar.close()

        logger.info('Sampling has finished sucessfully.')

        logger.info('Cleaning up any remaining intermediate .pdb files...')
        for pdb in glob.glob(os.path.join(OUTPUT_PATH,'*.pdb')):
            try:
                os.remove(pdb)
            except FileNotFoundError:
                continue

        final_log_msg = (f'There are {len(os.listdir(VALID_PDBS_DIR))} valid pdbs, '
                        f'{clashed_pdbs+disrespectful_pdbs} were discarded '
                        f'( {clashed_pdbs} clashed | {disrespectful_pdbs} violated constraints).')

        logger.info(final_log_msg)
        print(final_log_msg)
        return None # exit
    else:
        # Do not cleanup purposefully (any leftover will be overwritten later anyway)
        logger.info(f'There are {len(os.listdir(VALID_PDBS_DIR))} valid pdbs, '
                    f'{clashed_pdbs+disrespectful_pdbs} were discarded '
                    f'( {clashed_pdbs} clashed | {disrespectful_pdbs} violated constraints).')

    # Make sure valid pdbs starts from cutpoint
    if SS_BIAS_CUTPOINT is not None:
        valid_pdbs = SS_BIAS_CUTPOINT

    # Remove secondary structure bias
    SS_BIAS = None

    logger.info('Starting sampling without secondary structure bias...')
    while valid_pdbs < GOAL_ENSEMBLE_SIZE or len(os.listdir(VALID_PDBS_DIR)) < GOAL_ENSEMBLE_SIZE:

        # Submit new tasks if we haven't reached the target count
        while len(unfinished_obj_refs_batch) < MAX_PENDING_TASKS:
            sampling_future = sample_pdb.remote(stored_initial_pose,
                                                stored_databases,
                                                SAMPLING_TARGETS,
                                                OUTPUT_PATH,
                                                JOB_NAME,
                                                decoy_num,
                                                PYROSETTA_LOG,
                                                SS_BIAS,
                                                SAMPLER_PARAMS,
                                                SCOREFXN_ID,
                                                SCOREFXN_WEIGHT,
                                                MINIMIZER_ID,
                                                MINIMIZER_TOL,
                                                MINMIZER_MAXITERS,
                                                MINMIZER_FINALCYCLES,
                                                CONSTRAINTS_WEIGHT,
                                                CONSTRAINTS_VIOL_THRESHOLD,
                                                CONSTRAINTS_VIOL_MAXRES,)

            processing_future = process_pdb.remote(sampling_future,
                                                   FASPR_PATH,
                                                   PULCHRA_PATH,
                                                   stored_input_clashes,
                                                   SAMPLING_TARGETS,
                                                   VALID_PDBS_DIR,
                                                   GOAL_ENSEMBLE_SIZE,)

            unfinished_obj_refs_batch.append(processing_future)
            decoy_num += 1

        # Returns the ObjectRefs that are ready
        finished_tasks, unfinished_obj_refs_batch = ray.wait(unfinished_obj_refs_batch,
                                                             num_returns=1)

        for finished_task in finished_tasks:
            try:
                result = ray.get(finished_task)
            except ray.exceptions.RayTaskError as e:
                with open(RAY_LOG,'a',encoding='utf-8') as f:
                    f.write(str(e))
                continue
            else:
                # Process result
                if isinstance(result,type(None)):
                    disrespectful_pdbs += 1
                elif isinstance(result,bool):
                    if result:
                        clashed_pdbs += 1
                    else:
                        valid_pdbs += 1
                        if pbar.n < GOAL_ENSEMBLE_SIZE:
                            pbar.update(1)

    logger.info('Desired number of structurally non biased structures has been reached.')

    logger.info('Cancelling any remainings tasks and shutting down ray client...')
    ray.shutdown()

    # Close progress bar (required)
    pbar.set_description('Ensemblified!')
    pbar.close()

    logger.info('Sampling has finished sucessfully.')

    logger.info('Cleaning up any remaining intermediate .pdb files...')
    for pdb in glob.glob(os.path.join(OUTPUT_PATH,'*.pdb')):
        try:
            os.remove(pdb)
        except FileNotFoundError:
            continue

    final_log_msg = (f'There are {len(os.listdir(VALID_PDBS_DIR))} valid pdbs, '
                     f'{clashed_pdbs+disrespectful_pdbs} were discarded '
                     f'( {clashed_pdbs} clashed | {disrespectful_pdbs} violated constraints).')

    logger.info(final_log_msg)
    print(final_log_msg)
