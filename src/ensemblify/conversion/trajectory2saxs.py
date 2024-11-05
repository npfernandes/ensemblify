"""Calculate a SAXS curve (.dat) from a trajectory file (.xtc)."""

# IMPORTS
## Standard Library Imports
import os
from concurrent.futures import ProcessPoolExecutor

## Third Party Imports
import MDAnalysis as mda
import numpy as np
from tqdm import tqdm

## Local Imports
from ensemblify.conversion.conversion_utils import calc_saxs_data

# FUNCTIONS
def traj2saxs(
    trajectory: str,
    topology: str,
    trajectory_id: str,
    exp_saxs_file: str,
    ) -> str:
    """Calculate a theoretical SAXS curve from a trajectory file using PEPSI-SAXS.
    
    Calculation is done in chunks distributed across availabe processor cores.
    A Universe object is created with the given trajectory and topology, which
    allows for writing of temporary individual .pdb files from trajectory frames
    that are then used for SAXS curve calculation at every frame.

    Args:
        trajectory:
            path to trajectory file used in Universe object creation.
        topology:
            path to topology file used in Universe object creation.
        trajectory_id:
            prefix identifier for created files.
        exp_saxs_file:
            path to the experimental SAXS data for this protein, used in PEPSI-SAXS
            for SAXS curve calculation.

    Returns:
        calc_saxs_file:
            path to the calculated SAXS curve file, with the data calculated for
            every frame of the trajectory.

    Adapted from:
        https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
    """
    # Setup the Universe
    u = mda.Universe(topology, # topology
                     trajectory) # trajectory

    # Setup multiprocessing arguments
    frame_indices = range(len(u.trajectory))
    trajectory_ids = [trajectory_id] * len(frame_indices)
    universes = [u] * len(frame_indices)
    exp_saxs_files = [exp_saxs_file] * len(frame_indices)
    calc_saxs_log = os.path.join(os.path.split(exp_saxs_file)[0],f'{trajectory_id}_calc_saxs.log')
    calc_saxs_logs = [calc_saxs_log] * len(frame_indices)

    # Calculate SAXS data in chunks
    with ProcessPoolExecutor() as ppe:
        calc_saxs_chunks = list(tqdm(ppe.map(calc_saxs_data,
                                             trajectory_ids,
                                             universes,
                                             frame_indices,
                                             exp_saxs_files,
                                             calc_saxs_logs),
                                     total=len(frame_indices),
                                     desc=f'Calculating {trajectory_id} SAXS data... '))

    # Build the full SAXS data file
    ## Join the calculated chunks
    all_calc_saxs = np.vstack(calc_saxs_chunks)
    ## Make a column with the frame indices
    col0 = np.arange(1,len(all_calc_saxs) + 1).reshape(len(all_calc_saxs), 1)
    ## Join indices with data
    all_calc_saxs = np.hstack((col0,all_calc_saxs))

    # Save calculated SAXS data
    calc_saxs_file = os.path.join(os.path.split(exp_saxs_file)[0],f'{trajectory_id}_calc_saxs.dat')
    np.savetxt(calc_saxs_file,
               all_calc_saxs,
               encoding='utf-8')

    return calc_saxs_file
