"""Create a trajectory (.xtc) file from an ensemble of .pdb files."""

# IMPORTS
## Standard Library Imports
import os
import subprocess

## Third Party Imports
from tqdm import tqdm

## Local Imports
from ensemblify.conversion.conversion_utils import join_pdbs, move_topol_pdb

# FUNCTIONS
def ensemble2traj(
    job_name: str,
    ensemble_dir: str,
    trajectory_dir: str,
    trajectory_size: int | None = 10000,
    ) -> tuple[str,str]:
    """Create a trajectory (.xtc) file from an ensemble of .pdb files.
    
    Additionally, one of the .pdb files used to create this trajectory is kept as
    a .pdb topology file.

    Args:
        job_name:
            prefix identifier for any created files.
        ensemble_dir:
            path to directory where all the .pdb files are stored.
        trajectory_dir:
            path to directory where trajectory .xtc file will be created.
        trajectory_size:
            number of randomly sampled .pdb files to use for trajectory creation.
    
    Returns:
        A tuple (trajectory_path,topology_path) where:
            trajectory_path:
                path to created trajectory .xtc file.
            topology_path:
                path to created topology .pdb file.
    """

    # Create trajectory directory if non existent
    if not os.path.isdir(trajectory_dir):
        os.mkdir(trajectory_dir)

    # Setup trajectory path and creation log file
    trajectory_creation_log = os.path.join(trajectory_dir,'ensemble_to_xtc.log')
    trajectory_path = os.path.join(trajectory_dir,f'{job_name}_trajectory.xtc')

    # Initialize pbar
    with tqdm(total=4,unit='step') as pbar:

        # Keep one .pdb to serve as topology file for later analysis
        pbar.set_description(f'Moving {job_name} topology pdb... ')
        topology_path = move_topol_pdb(job_name=job_name,
                                       origin_path=ensemble_dir,
                                       destination_path=trajectory_dir)
        pbar.update(1)

        # Join pdbs into a single multimodel .pdb file
        pbar.set_description(f'Joining {job_name} pdbs... ')
        ensemble_pdb_path = join_pdbs(pdbs_dir=ensemble_dir,
                                      job_name=job_name,
                                      ensemble_dir=trajectory_dir,
                                      n_models=trajectory_size)
        pbar.update(1)

        # From a multimodel .pdb file, create a .xtc trajectory file
        pbar.set_description(f'Creating {job_name} trajectory... ')
        subprocess.run(['gmx', 'trjconv', '-f', ensemble_pdb_path, '-o', f'{trajectory_path}'],
                       stdout=open(trajectory_creation_log,'a',encoding='utf-8'),
                       stderr=subprocess.STDOUT,
                       check=True)
        pbar.update(1)

        # Remove created multi_model ensemble
        if os.path.isfile(ensemble_pdb_path):
            os.remove(ensemble_pdb_path)
        pbar.update(1)

        pbar.set_description(f'{job_name} trajectory creation completed! ')

    return trajectory_path,topology_path
