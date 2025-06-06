"""Auxiliary functions for the conversion module."""

# IMPORTS
## Standard Library Imports
import glob
import os
import random
import subprocess
import warnings

## Third Party Imports
import MDAnalysis as mda
import numpy as np

## Local Imports
from ensemblify.config import GLOBAL_CONFIG

# FUNCTIONS
def move_topology_pdb(
    topology_name: str,
    origin_dir: str,
    destination_dir: str,
    ) -> str:
    """Move a .pdb file from origin to destination directory, to act as a topology file.
    
    Given a path to a directory containing .pdb files, moves a single .pdb from that directory
    to a destination directory, to serve as a topology file for future trajectory analysis.

    Args:
        topology_name (str):
            Prefix identifier for moved .pdb file.
        origin_dir (str):
            Directory where file of interest is located.
        destination_dir (str):
            Directory where file will be moved to.
    
    Returns:
        str:
            Path to the moved topology file.
    """
    topology_path = os.path.join(destination_dir,f'{topology_name}_top.pdb')
    for pdb in glob.glob(os.path.join(origin_dir,'*.pdb')):
        with open(pdb,'r',encoding='utf-8-sig') as f, open(topology_path,'w',encoding='utf-8') as t:
            t.write(f.read())
            break

    return topology_path


def _sample_without_topology(
    pdbs_dir: str,
    topology_path: str | None = None,
    n_models: int | None = None,
    ) -> list[str]:
    """Sample .pdb files from a directory, excluding a topology file.
    
    Args:
        pdbs_dir (str):
            Path to directory where .pdb files are stored.
        topology_path (str):
            Path to a topology .pdb file to be ignored when sampling .pdb files from
            pdbs_dir. Defaults to None.
        n_models (int):
            Number of .pdb files to randomly sample from the provided directory.
            Defaults to all .pdb files in the directory.

    Returns:
        list[str]:
            List of paths to sampled .pdb files."""
    
    # Grab total .pdb files
    total_pdbs = glob.glob(os.path.join(pdbs_dir,'*.pdb'))
    
    # Remove topology_path from the list of .pdb files, if applicable
    if topology_path is not None:
        try:
            total_pdbs.remove(topology_path)
        except ValueError:
            pass
    
    # Assign number of models to sample if not given
    if n_models is None:
        n_models = len(total_pdbs)
    
    # Sample desired number of .pdb files (without replacement)
    sampled_pdbs = random.sample(total_pdbs, n_models)
    
    return sampled_pdbs


def join_pdbs(
    pdbs_dir: str,
    multimodel_name: str,
    multimodel_dir: str,
    n_models: int | None = None,
    topology_path: str | None = None,
    ) -> str:
    """Join a randomly sampled number of .pdb files in a directory into a single multimodel
    .pdb file.

    Args:
        pdbs_dir (str):
            Path to directory where numbered .pdb files are stored.
        multimodel_name (str):
            Prefix identifier for created multimodel .pdb file.
        multimodel_dir (str):
            Path to directory where ensemble pdb will be created.
        n_models (int):
            Number of .pdb files to randomly sample from the specified directory.
            If None, all .pdb files in the directory will be used.
        topology_path (str):
            Path to a topology .pdb file to be ignored when sampling .pdb files from
            multimodel_dir. Required if pdbs_dir matches multimodel_dir. Defaults to None.
    
    Returns:
        str:
            Path to created multimodel ensemble .pdb file.
    """  
    # Sample desired number of .pdb files (without replacement)
    sampled_pdbs2join = _sample_without_topology(pdbs_dir=pdbs_dir,
                                                 topology_path=topology_path,
                                                 n_models=n_models)

    # Setup multimodel .pdb filepath
    ensemble_path = os.path.join(multimodel_dir,f'{multimodel_name}_ensemble.pdb')

    # Write the sampled .pdb files to the multimodel .pdb file
    with open(ensemble_path,'x',encoding='utf-8') as output:
        model = 1
        for pdb in sampled_pdbs2join:
            output.write(f'MODEL {model}\n')
            output.write(f'REMARK {pdb}\n')
            with open(pdb, 'r',encoding='utf-8-sig') as pdb_file:
                content = pdb_file.readlines()
                for line in content:
                    if line.startswith(('ATOM', 'HETA', 'TER')):
                        output.write(line)
            output.write('ENDMDL\n')
            model += 1
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
        trajectory_id (str):
            Prefix identifier for created files.
        universe (mda.Universe):
            Universe object containing the trajectory being analyzed.
        frame_index (int):
            Current frame being calculated.
        exp_saxs_file (str):
            Path to .dat file with experimental SAXS data for the current
            protein to be used by PEPSI-SAXS.
        calc_saxs_log (str):
            Path to .log file for the SAXS curve calculation of each frame.

    Returns:
        np.ndarray:
            Numpy ndarray with the values for the SAXS curve calculated from this
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
        with mda.Writer(frame_file, universe.atoms.n_atoms) as w:
            w.write(universe.atoms)

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
