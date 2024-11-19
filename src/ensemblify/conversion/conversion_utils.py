"""Auxiliary functions for the conversion module."""

# IMPORTS
## Standard Library Imports
import glob
import os
import subprocess
import warnings

## Third Party Imports
import MDAnalysis as mda
import numpy as np

## Local Imports
from ensemblify.config import GLOBAL_CONFIG

# FUNCTIONS
def move_topol_pdb(
    job_name: str,
    origin_path: str,
    destination_path:str,
    ) -> str:
    """Move a .pdb file from origin to destination directory.
    
    Given a path to a directory containing .pdb files, moves a single .pdb
    from that directory to a destination directory, to serve as a topology
    file for future trajectory analysis.

    Args:
        job_name:
            prefix identifier for moved .pdb file.
        origin_path:
            directory where topology file of interest is located.
        destination_path:
            directory where topology file will be moved to.
    
    Returns:
        topology_path:
            path to the moved topology file.
    """
    topology_path = os.path.join(destination_path,f'{job_name}_top.pdb')
    for pdb in glob.glob(os.path.join(origin_path,'*.pdb')):
        with open(pdb,'r',encoding='utf-8-sig') as f, open(topology_path,'w',encoding='utf-8') as t:
            t.write(f.read())
            break

    return topology_path


def join_pdbs(
    pdbs_dir: str,
    job_name: str,
    ensemble_dir: str,
    n_models: int,
    ) -> str:
    """Join a randomly sampled number of .pdb files in a directory into a single multimodel
    .pdb file.

    Args:
        job_name:
            prefix identifier for created multimodel .pdb file.
        ensemble_dir:
            path to directory where ensemble pdb will be created.
        pdbs_dir:
            path to directory where numbered .pdb files are stored.
        n_models:
            number of .pdb files to randomly sample from the specified directory.
    
    Returns:
        ensemble_path:
            path to created multimodel ensemble .pdb file.
    """
    ensemble_path = os.path.join(ensemble_dir,f'{job_name}_ensemble.pdb')
    with open(ensemble_path,'x',encoding='utf-8') as output:
        model = 1
        for pdb in glob.glob(os.path.join(pdbs_dir,'*.pdb')):
            if model <= n_models:
                output.write(f'MODEL {model}\n')
                output.write(f'REMARK {pdb}\n')
                with open(pdb, 'r',encoding='utf-8-sig') as pdb_file:
                    content = pdb_file.readlines()
                    for line in content:
                        if line.startswith(('ATOM', 'HETA', 'TER')):
                            output.write(line)
                output.write('ENDMDL\n')
                model += 1
            else:
                break
    return ensemble_path


def calc_saxs_data(
    trajectory_id: str,
    universe: mda.Universe,
    frame_index: int,
    exp_saxs_file: str,
    calc_saxs_log: str,
    ) -> np.ndarray:
    """Calculate a theoretical SAXS curve for a frame of a MDAnalysis
    Universe object using PEPSI-SAXS.

    Calculation is done for a single temporary .pdb file created  from
    the current frame of the Universe object.

    Args:
        trajectory_id:
            prefix identifier for created files.
        universe:
            Universe object containing the trajectory being analyzed.
        frame_index:
            current frame being calculated.
        exp_saxs_file:
            path to .dat file with experimental SAXS data for the current
            protein to be used by PEPSI-SAXS.
        calc_saxs_log:
            path to .log file for the SAXS curve calculation of each frame.

    Returns:
        calc_saxs:
            ndarray with the values for the SAXS curve calculated from this
            frame of the trajectory in the Universe object.
    """
    # Setup tmp files
    frame_file = os.path.join(os.path.split(exp_saxs_file)[0],
                              f'{trajectory_id}_tmp_frame_{frame_index}.pdb')
    output_file = os.path.join(os.path.split(exp_saxs_file)[0],
                               f'{trajectory_id}_tmp_saxs_{frame_index}.dat')

    # Set the Universe to point to this frame
    universe.trajectory[frame_index]  # no need to assign variable

    # Save the current frame to a tmp file
    with warnings.catch_warnings():
        # Suppress UserWarnings related to Unit cell dimensions and 'formalcharges'
        warnings.filterwarnings('ignore', category=UserWarning)
        with mda.Writer(frame_file, universe.atoms.n_atoms) as W:
            W.write(universe.atoms)

    # Calculate SAXS data for this frame
    pepsi_saxs_path = GLOBAL_CONFIG['PEPSI_SAXS_PATH']
    assert pepsi_saxs_path is not None, 'Pepsi-SAXS installation not found!'

    pepsi_comm = f'{pepsi_saxs_path} {frame_file} {exp_saxs_file} -o {output_file} -cst -x'
    subprocess.run(pepsi_comm.split(),
                   stdout=open(calc_saxs_log,'a',encoding='utf-8'),
                   stderr=subprocess.STDOUT,
                   check=True)

    calc_saxs = np.loadtxt(output_file)[..., 3] # saxs data for frame

    # Clean up tmp files
    os.remove(frame_file)
    os.remove(output_file)

    return calc_saxs
