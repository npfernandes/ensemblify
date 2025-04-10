"""Create a trajectory (.xtc) file from an ensemble of .pdb files."""

# IMPORTS
## Standard Library Imports
import os
import subprocess

## Third Party Imports
from tqdm import tqdm

## Local Imports
from ensemblify.conversion.conversion_utils import join_pdbs, move_topology_pdb

# CONSTANTS
MOVE_PDB_MSG = 'Moving{}topology .pdb... '
JOIN_PDBS_MSG = 'Joining{}.pdbs... '
CREATE_TRJ_MSG = 'Creating{}trajectory... '
RMV_MM_PDB_MSG = 'Removing{}multimodel pdb... '
FINAL_MSG = '{}Trajectory creation complete! '

# FUNCTIONS
def ensemble2traj(
    ensemble_dir: str | None = os.getcwd(),
    trajectory_dir: str | None = os.getcwd(),
    trajectory_id: str | None = '',
    trajectory_size: int | None = None,
    ) -> tuple[str,str]:
    """Create a trajectory (.xtc) file from an ensemble of .pdb files.
    
    Additionally, one of the .pdb files used to create this trajectory is kept as
    a .pdb topology file.

    Args:
        ensemble_dir (str):
            Path to directory where all the .pdb files are stored. Defaults to current working
            directory.
        trajectory_dir (str):
            Path to directory where trajectory .xtc file will be created. Will be created if it
            does not exist. Defaults to current working directory.
        trajectory_id (str):
            Prefix identifier for any created files.
        trajectory_size (int):
            Number of randomly sampled .pdb files to use for trajectory creation.
            Defaults to all .pdb files in the ensemble directory.
    
    Returns:
        tuple[str,str]:
            trajectory_path (str):
                Path to created trajectory .xtc file.
            topology_path (str):
                Path to created topology .pdb file.
    """
    # Create trajectory directory if non existent
    if not os.path.isdir(trajectory_dir):
        os.mkdir(trajectory_dir)

    # Setup trajectory name
    if trajectory_id:
        trajectory_id_msg = ' ' + trajectory_id + ' '
    else:
        trajectory_id_msg = ' '

    # Setup trajectory path and creation log file
    trajectory_creation_log = os.path.join(trajectory_dir,'ensemble_to_xtc.log')
    trajectory_path = os.path.join(trajectory_dir,f'{trajectory_id}_trajectory.xtc')

    # Initialize pbar
    with tqdm(total=4,unit='step') as pbar:

        # Keep one .pdb to serve as topology file for later analysis
        pbar.set_description(MOVE_PDB_MSG.format(trajectory_id_msg))
        topology_path = move_topology_pdb(topology_name=trajectory_id,
                                          origin_dir=ensemble_dir,
                                          destination_dir=trajectory_dir)
        pbar.update(1)

        # Join pdbs into a single multimodel .pdb file
        pbar.set_description(JOIN_PDBS_MSG.format(trajectory_id_msg))
        ensemble_pdb_path = join_pdbs(pdbs_dir=ensemble_dir,
                                      multimodel_name=trajectory_id,
                                      multimodel_dir=trajectory_dir,
                                      n_models=trajectory_size,
                                      topology_path=topology_path)
        pbar.update(1)

        # From a multimodel .pdb file, create a .xtc trajectory file
        pbar.set_description(CREATE_TRJ_MSG.format(trajectory_id_msg))
        try:
            subprocess.run(['gmx', 'trjconv',
                            '-f', ensemble_pdb_path, '-o', f'{trajectory_path}'],
                            stdout=open(trajectory_creation_log,'a',encoding='utf-8'),
                            stderr=subprocess.STDOUT,
                            check=True)
        except FileNotFoundError:
            subprocess.run([f'{os.environ.get("GMXBIN")}/gmx', 'trjconv',
                            '-f', ensemble_pdb_path, '-o', f'{trajectory_path}'],
                            stdout=open(trajectory_creation_log,'a',encoding='utf-8'),
                            stderr=subprocess.STDOUT,
                            check=True)
        pbar.update(1)

        # Remove created multi_model ensemble
        pbar.set_description(RMV_MM_PDB_MSG.format(trajectory_id_msg))
        if os.path.isfile(ensemble_pdb_path):
            os.remove(ensemble_pdb_path)
        pbar.update(1)

        pbar.set_description(FINAL_MSG.format(trajectory_id_msg))

    return trajectory_path,topology_path
