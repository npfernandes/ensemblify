"""Auxiliary functions to read, manipulate and write .pdb files."""

# IMPORTS
## Standard Library Imports
import os

## Third Party Imports
import pandas as pd
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

# FUNCTIONS
def df_from_pdb(pdb: str) -> pd.DataFrame:
    """Convert the information in a .pdb file into a pandas DataFrame using BioPDB.

    Args:
        pdb (str):
            Path to the .pdb file.

    Returns:
        pd.DataFrame:
            The given pdb's information in DataFrame format.
    """

    # Create a PDBParser object to parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb)

    # Create empty lists to store atom information
    atom_data = []

    # Iterate through the structure to extract atom information
    for model in structure:
        for chain in model:
            i = 1
            for residue in chain:
                for atom in residue:
                    atom_info = {
                        'AtomNumber': atom.get_serial_number() ,
                        'AtomName': atom.get_id(),
                        'ResidueName': residue.get_resname(),
                        'ChainID': chain.get_id(),
                        'ResidueNumber': residue.id[1],
                        'X': atom.get_coord()[0],
                        'Y': atom.get_coord()[1],
                        'Z': atom.get_coord()[2],
                        'Occupancy': atom.get_occupancy(),
                        'B-Factor': atom.get_bfactor(),
                        'Element': atom.element
                    }
                    atom_data.append(atom_info)
                i+= 1

    # Create a DataFrame from the extracted atom data
    df = pd.DataFrame(atom_data)

    return df


def df_to_pdb(df: pd.DataFrame, output_pdb_filename: str):
    """Write content of a DataFrame containing PDB file info as a .pdb file using BioPDB.

    Args:
        df (pd.DataFrame):
            DataFrame containing PDB information.
        output_pdb_filename (str):
            Path to the output .pdb.
    """
    # Create an instance of PDBIO to write the PDB file
    pdb_io = PDBIO()

    # Create a Structure object to hold the PDB data
    structure = Structure('protein')

    # Create a Model
    model = Model(0)

    # Iterate through each Atom in the df and add it to the correct Residue,Chain
    current_chain = None
    current_residue = None
    for _, row in df.iterrows():
        chain_id = row['ChainID']
        residue_id = (' ', row['ResidueNumber'], ' ')

        # Check if we're still in the same chain
        if current_chain is None or current_chain.id != chain_id:
            current_chain = Chain(chain_id)
            model.add(current_chain) # add the Chain to the Model
            current_residue = None

        # Check if we're still in the same residue
        if current_residue is None or current_residue.id != residue_id:
            current_residue = Residue(residue_id, row['ResidueName'],'')
            current_chain.add(current_residue) # add the Residue to the Chain

        # Create an Atom object
        atom = Atom(name=row['AtomName'], coord=(row['X'], row['Y'], row['Z']),
                    bfactor=row['B-Factor'], occupancy=row['Occupancy'],
                    altloc=' ', fullname=row['AtomName'],
                    serial_number= ['AtomNumber'], element=row['Element'])

        # Add the Atom to the Residue
        current_residue.add(atom)

    # Add the Model to the Structure
    structure.add(model)

    # Write the Structure to a PDB file
    pdb_io.set_structure(structure)
    pdb_io.save(output_pdb_filename)


def extract_pdb_info(pdb: str) -> dict[int,tuple[str,int,int]]:
    """Extract from a .pdb file info about number of chains, chain letters, starting residue
    numbers and chain size.

    Args:
        topology (str):
            Path to .pdb topology file.
    
    Returns:
        dict[int,tuple[str,int,int]]:
            Mapping of chain numbers to their letter, starting residue number and chain size.
    """
    # Setup info dict
    ch_res = {}

    # Read pdb into DataFrame
    df = df_from_pdb(pdb)
    nrows = df.shape[0]

    # If ChainID is fully empty, replace it with 'A'
    if df['ChainID'].isin([' ']).all():
        df['ChainID'] = 'A'

    # Setup current chain, rescount for chainsize and chain number id
    curr_chain_id = df.loc[nrows-1,'ChainID']
    rescount = 1
    chain_number = 1

    # Iterate from the bottom of the .pdb, skipping last res (because we look at i+1)
    for i in range(nrows-2,-1,-1):
        res_num = df.loc[i,'ResidueNumber']
        chain_id = df.loc[i,'ChainID']
        if chain_id != curr_chain_id:
            ch_res[chain_number] = (curr_chain_id,df.loc[i+1,'ResidueNumber'],rescount)
            chain_number += 1
            curr_chain_id = chain_id
            rescount = 1
        elif res_num != df.loc[i+1,'ResidueNumber']:
            rescount += 1

    # Add the first chain to the dict (as there is no change in ChainID while iterating
    # through to the first residue of the .pdb, this chain would never be added)
    ch_res[chain_number] = (df.loc[0,'ChainID'], df.loc[0,'ResidueNumber'],rescount)

    return ch_res


def cleanup_pdbs(pdbs: list[str]):
    """Delete all .pdb files in the given list.

    Args:
        pdbs (list[str]):
            Paths to .pdb files to delete.
    """
    for leftover_pdb in pdbs:
        try:
            os.remove(leftover_pdb)
        except FileNotFoundError:
            continue
