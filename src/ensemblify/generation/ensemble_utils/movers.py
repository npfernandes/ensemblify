"""Custom Mover classes created from the PyRosetta Mover class."""

# IMPORTS
## Standard Library Imports
import math

## Third Party Imports
import numpy as np
import pandas as pd
import pyrosetta
from pyrosetta.rosetta.protocols.moves import Mover

## Local Imports
from ensemblify.generation.ensemble_utils.movers_utils import get_ss_bounds
from ensemblify.config import GLOBAL_CONFIG

# CLASSES
class SetRandomDihedralsMover(Mover):
    """
    Custom PyRosetta Mover object that sets random dihedral angles in
    target residues, taken from a given database.

    Inherits from pyrosetta.rosetta.protocols.moves.Mover.

    Attributes:
        databases (dict):
            all the available databases to sample from. Mapping of database_ids to
            databases nested dicts, that map residue 1lettercodes to dihedral
            angle values dataframes.
        log_file (str):
            path to .log file for warnings or error messages related to sampling.
    """

    def __init__(self, databases: dict[str,dict[str,pd.DataFrame]], log_file: str):
        """Initializes the instance from a dictionary of database(s) to sample from.

        Args:
            databases:
                mapping of database_ids to databases nested dicts, that map residue 1lettercodes
                to dihedral angle values dataframes.
            log_file:
                path to .log file for warnings or error messages related to sampling.
        """
        Mover.__init__(self)
        self.databases = databases
        self.log_file = log_file

    def get_name(self):
        """Return the name of this mover."""
        return self.__class__.__name__

    def apply(self,
        pose: pyrosetta.rosetta.core.pose.Pose,
        target_resnum: int,
        database_id: str,
        secondary_structure: str | None,
        sampling_mode: str):
        """Apply the mover to a Pose, on the given residue number (PyRosetta Pose numbering).

        Args:
            pose (pyrosetta.rosetta.core.pose.Pose):
                Pose on which the mover will be applied.
            target_resnum:
                residue number where to apply the mover (PyRosetta Pose numbering).
            database_id:
                identifier for which database to sample from.
            secondary_structure:
                which secondary structure element to force in this target region.
                If None, sample database without restraint.
            sampling_mode:
                whether to sample the database considering neighbouring residues ('TRIPEPTIDE')
                or not ('SINGLERESIDUE').
        """
        # Slice sequence given target res number
        fragment = pose.sequence()[target_resnum-2:target_resnum+1]

        if fragment: # avoid attempting to sample empty fragments (e.g. if trying to slice [-1:2])

            # Get desired database
            database = self.databases[database_id]

            # Get residue
            if len(fragment) == 3:
                residue = fragment[1]
            else:
                # single residue databases
                residue = fragment

            # Get angles for this res from database
            res_dihedrals = database[residue]

            # Respect sampling mode
            if sampling_mode == 'TRIPEPTIDE':
                # Filter res db to look only at curr frag
                all_dihedrals = res_dihedrals[res_dihedrals['FRAG'].str.match(fragment)]

            elif sampling_mode == 'SINGLERESIDUE':
                # Do not filter db by fragment
                all_dihedrals = res_dihedrals

            # Respect secondary structure bias
            if secondary_structure is not None:
                # Filter db to include only dihedrals for the desired secondary structure
                phi_bounds, psi_bounds = get_ss_bounds(secondary_structure)
                phi_mask = all_dihedrals['PHI2'].map(lambda x: math.radians(phi_bounds[0]) <= x <= math.radians(phi_bounds[1]))
                psi_mask = all_dihedrals['PSI2'].map(lambda x: math.radians(psi_bounds[0]) <= x <= math.radians(psi_bounds[1]))
                dihedrals = all_dihedrals[phi_mask & psi_mask]
            else:
                dihedrals = all_dihedrals

            # Check if the previous filtering step resulted in at least one valid row to sample from
            if dihedrals.shape[0] == 0:
                if secondary_structure == 'alpha_helix':
                    # Set canonical alpha helix values with omega 180
                    phi, psi = GLOBAL_CONFIG['ALPHA_HELIX_CANON']
                    omg = 180
                elif secondary_structure == 'beta_strand':
                    # Set canonical beta strand values with omega 180
                    phi, psi = GLOBAL_CONFIG['BETA_STRAND_CANON']
                    omg = 180

                pose.set_phi(target_resnum,phi)
                pose.set_psi(target_resnum,psi)
                pose.set_omega(target_resnum,omg)

                no_rows_found_msg = ('No rows respecting filters found in database, canonical '
                                     f'values for {secondary_structure} set on residue '
                                     f'{target_resnum}!\n')
                print(no_rows_found_msg, end='')
                with open(self.log_file,'a',encoding='utf-8') as f:
                    f.write(no_rows_found_msg)

            else:
                # Choose row from which to take dihedral values
                choice_angle= np.random.randint(dihedrals.shape[0])

                # Reset indexing (if we filtered the database by fragment it became non-sequential)
                dihedrals_indexed = dihedrals.reset_index(drop=True)

                # Get the dihedral values and convert from radians to degrees
                phi = dihedrals_indexed['PHI2'][choice_angle]
                psi = dihedrals_indexed['PSI2'][choice_angle]
                omg = dihedrals_indexed['OMG2'][choice_angle]

                phi = math.degrees(phi)
                psi = math.degrees(psi)
                omg = math.degrees(omg)

                # Perform the mover operation (change Phi/Psi/Omega angles)
                # sample from a normal distribution to add diversity
                pose.set_phi(target_resnum,np.random.normal(phi,abs(phi*0.10)))
                pose.set_psi(target_resnum,np.random.normal(psi,abs(psi*0.10)))
                pose.set_omega(target_resnum,np.random.normal(omg,abs(omg*0.10)))


# FUNCTIONS
def setup_mover(
    mover_id: str,
    databases: dict[str,dict[str,pd.DataFrame]],
    log_file: str,
    ) -> pyrosetta.rosetta.protocols.moves.Mover:
    """Create custom PyRosetta Mover.
    
    Setup Mover object given a database to sample from and a mover id.

    Args:
        mover_id:
            identifier for which CustomMover to create.
        databases:
            mapping of database_ids to databases nested dicts, that map residue 1lettercodes
            to dihedral angle values dataframes.
        log_file:
            path to .log file for warnings or error messages related to sampling.

    Returns:
        my_mover (pyrosetta.rosetta.protocols.moves.Mover):
            Custom PyRosetta Mover object.
    """
    if mover_id == 'set_random_dihedrals':
        my_mover = SetRandomDihedralsMover(databases,log_file)
    return my_mover
