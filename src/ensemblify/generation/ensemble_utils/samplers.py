"""Custom Sampler classes to generate protein conformations."""

# IMPORTS
## Standard Library Imports
import random

## Third Party Imports
import numpy as np
import pandas as pd
import pyrosetta
from pyrosetta.rosetta.core.pose import Pose

## Local Imports
from ensemblify.generation.ensemble_utils.movers import setup_mover

# CLASSES
class MonteCarloSampler():
    """Custom MonteCarlo sampler for sampling dihedral angles.
    
    Attributes:
        scorefxn (pyrosetta.rosetta.core.scoring.ScoreFunction):
            PyRosetta score function to be used for evaluating Pose objects
            during sampling.
        databases (dict):
            All the available databases to sample from. Mapping of database_ids to
            databases nested dicts, that map residue 1lettercodes to dihedral
            angle values dataframes.
        mover (pyrosetta.rosetta.protocols.moves.Mover):
            Custom PyRosetta Mover used to apply dihedral angle changes to a Pose.
        params (dict):
            Hyperparameters for this sampler (temperature and maximum loops):
                temperature (int):
                    A measure of how probable it is to accept Pose objects with a worse score than
                    the current one after applying our Mover, according to the acceptance criterion.
                maximum loops (int):
                    The maximum amount of attempts without accepting a Move before moving on to
                    the next residue to sample.
        log_file (str):
            Path to .log file for warnings or error messages related to sampling.
    """

    def __init__(self,
        scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction,
        databases: dict,
        mover_id: str,
        smp_params: dict[str,int],
        variance: float,
        log_file: str):
        """Initializes the instance based on the given parameters.
        
        Args:
            scorefxn (pyrosetta.rosetta.core.scoring.ScoreFunction):
                PyRosetta score function, with desired constraints already added.
            databases (dict):
                Mapping of database_ids to databases nested dicts, that map residue 1lettercodes
                to dihedral angle values dataframes.
            mover_id (str):
                Identifier for which mover to use in this sampler object.
            smp_params (dict[str,int]):
                Parameters for the instantiated sampler.
            variance (float):
                New dihedral angle values inserted into sampling regions are sampled from a
                Gaussian distribution centered on the value found in database and percentage
                variance equal to this value.
            log_file (str):
                Path to .log file for warnings or error messages related to sampling.
        """
        self.scorefxn = scorefxn
        self.databases = databases
        self.params = smp_params  #smp_params is a dictionary
        self.mover = setup_mover(mover_id,databases,variance,log_file)
        self.log_file = log_file

    def apply(self,
        pose: pyrosetta.rosetta.core.pose.Pose,
        target: list[int],
        chain: str,
        database_id: str,
        ss_bias: tuple[tuple[str,tuple[int,int],str],...] | None,
        sampling_mode: str):
        """Perform MC sampling on the given pose, in the given target residue range.

        Args:
            pose (pyrosetta.rosetta.core.pose.Pose):
                Pose to be modified during sampling.
            target (list[int]):
                Residue range on which sampling will be applied.
            chain (str):
                Letter identifier for the current chain being sampled.
            database_id (str):
                Identifier for which database to sample from.
            ss_bias (tuple[tuple[str,tuple[int,int],str],...] | None):
                Information about types of secondary structure biases, including which chain and
                residue numbers they should be applied on.
            sampling_mode (str):
                Whether to sample the database considering neighbouring residues ('TRIPEPTIDE')
                or not ('SINGLERESIDUE').

        """
        # Shuffling step (Increase diversity between sampled structures)
        target_shuffled = random.sample(target, len(target))

        # Calculate initial score
        current_score = self.scorefxn(pose)

        # Setup working Pose
        working_pose = Pose()
        working_pose.assign(pose)

        # Sample each residue
        for target_res in target_shuffled:

            # Looping variables
            tryagain = True
            loops = 0

            while tryagain and loops < self.params['max_loops']:

                # Check for secondary structure bias
                secondary_structure = None
                if ss_bias is not None:
                    for chain_region_ss in ss_bias:
                        c, r, ss = chain_region_ss
                        if c == chain and target_res in range(r[0],r[1]+1):
                            secondary_structure = ss

                # Sample dihedrals, following sampling mode and respecting any ss bias
                self.mover.apply(working_pose,
                                 target_res,
                                 database_id,
                                 secondary_structure,
                                 sampling_mode)

                # Score Pose created by applying mover
                new_score = self.scorefxn(working_pose)

                # If it was a good move accept, else acceptance criterion
                if new_score <= current_score:
                    pose.assign(working_pose)
                    tryagain = False
                    current_score = new_score
                else:
                    delta = new_score - current_score
                    prob = min(1,np.exp(-delta/self.params['temperature']))
                    n = random.random()
                    # Accept according to criterion, or undo move if rejected
                    if n < prob:
                        pose.assign(working_pose)
                        tryagain = False
                        current_score = new_score
                    else:
                        working_pose.assign(pose)
                loops += 1

            # Log MC failure for this res
            if loops == self.params['max_loops']:
                print(f'Max loops reached: residue {target_res}\n')
                with open(self.log_file,'a',encoding='utf-8') as f:
                    f.write(f'Max loops reached: residue {target_res}\n')


# FUNCTIONS
def setup_samplers(
    sampler_params: dict[str,dict[str,int]],
    variance: float,
    scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction,
    databases: dict[str,dict[str,pd.DataFrame]],
    log_file: str,
    ) -> dict[str,MonteCarloSampler]:
    """Create all Sampler objects to be used during sampling.
    
    Create a dictionary with all the samplers that will be used during sampling,
    given a list of sampler_ids and certain parameters.

    Args:
        sampler_params (dict[str,dict[str,int]):
            Parameters for each sampler to setup.
        scorefxn (pyrosetta.rosetta.core.scoring.ScoreFunction):
            PyRosetta score function, with desired constraints already added.
        databases (dict[str,dict[str,pd.DataFrame]]):
            Mapping of database_ids to databases nested dicts, that map residue 1lettercodes
            to dihedral angle values dataframes.

    Returns:
        dict[str,MonteCarloSampler]:
            Mapping of sampler_ids to sampler objects to use during sampling.
    """
    sampler_ids = list(sampler_params.keys())
    samplers = {}
    for smp_id in sampler_ids:
        if smp_id == 'MC':
            samplers[smp_id] = MonteCarloSampler(scorefxn=scorefxn,
                                                 databases=databases,
                                                 mover_id='set_random_dihedrals',
                                                 smp_params=sampler_params['MC'],
                                                 variance=variance,
                                                 log_file=log_file)
        # add more samplers here

    return samplers
