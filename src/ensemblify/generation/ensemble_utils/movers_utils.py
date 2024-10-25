"""Auxiliary functions for Custom PyRosetta Movers and to read database files into memory."""

# IMPORTS
## Standard Library Imports
import os

## Third Party Imports
import pandas as pd

## Local Imports
from ensemblify.config import GLOBAL_CONFIG

# CONSTANTS
ALLOWED_DATABASE_FORMATS = ['.csv','.parquet','.pkl']
ALLOWED_SECONDARY_STRUCTURE = ['alpha_helix','beta_strand']
DATABASE_OPTIMIZED_COL_DTYPES = { GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['OMG1'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['OMG2'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['OMG3'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PHI1'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PHI2'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PHI3'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PSI1'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PSI2'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PSI3'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['FRAG'] : 'category'}

# FUNCTIONS
def read_database(database_path: str) -> pd.DataFrame:
    """Read a database file into a pandas DataFrame.
    
    If possible, read only the desired set of columns into a pandas.DataFrame
    (depends on database file format).

    Args:
        database_path:
            filepath to database file.

    Returns:
        database: database as a pandas.DataFrame.
    """
    assert os.path.splitext(database_path)[1] \
          in ALLOWED_DATABASE_FORMATS, f'Database format must be in {ALLOWED_DATABASE_FORMATS}'

    if database_path.endswith('.csv'): # pick columns, convert dtypes
        database = pd.read_csv(database_path,
                               usecols=list(GLOBAL_CONFIG['USED_DATABASE_COLNAMES'].values()),
                               dtype=DATABASE_OPTIMIZED_COL_DTYPES)

    elif database_path.endswith('.parquet'): # pick columns
        database = pd.read_parquet(database_path,
                                   columns=list(GLOBAL_CONFIG['USED_DATABASE_COLNAMES'].values()))

    elif database_path.endswith('.pkl'): # no customization
        database = pd.read_pickle(database_path)

    # .upper() all column names
    database.rename(columns=str.upper, inplace=True)

    return database


def trim_database(database: pd.DataFrame,columns_to_keep: list[str]):
    """Removes columns in a database whose names are not in a given list.

    Args:
        database:
            target database in DataFrame format.
        columns_to_keep:
            columns names to keep in the database.
    """
    # Identify columns to drop
    columns_to_drop = list(set(database.columns) - set(columns_to_keep))

    # Drop columns not present in the list, inplace to save memory
    database.drop(columns=columns_to_drop,
                  inplace=True)


def optimize_database(database: pd.DataFrame) -> dict[str,pd.DataFrame]:
    """Reduce a database's memory usage to what is strictly necessary.

    Datatypes of database's columns are optimized and database is broken
    into 20 pieces, one for each aminoacid residue.

    Args:
        database:
            unoptimized dihedral angle database.

    Returns:
        res_angles:
            mapping of aminoacid 1lettercode to their corresponding
            dihedral angle values in the optimized database.
    """

    # Optimize our database, chaging datatypes to ones that are appropriate but use less memory
    optimized_db = database.astype(DATABASE_OPTIMIZED_COL_DTYPES)

    # Freeup memory
    del database

    # Now we break apart our full database into a dataFrame for each aminoacid residue and
    # create the new database dictionary
    res_angles = {
        'A' : optimized_db[optimized_db['FRAG'].str.contains('.+A.+')],
        'R' : optimized_db[optimized_db['FRAG'].str.contains('.+R.+')],
        'N' : optimized_db[optimized_db['FRAG'].str.contains('.+N.+')],
        'D' : optimized_db[optimized_db['FRAG'].str.contains('.+D.+')],
        'C' : optimized_db[optimized_db['FRAG'].str.contains('.+C.+')],
        'Q' : optimized_db[optimized_db['FRAG'].str.contains('.+Q.+')],
        'E' : optimized_db[optimized_db['FRAG'].str.contains('.+E.+')],
        'G' : optimized_db[optimized_db['FRAG'].str.contains('.+G.+')],
        'H' : optimized_db[optimized_db['FRAG'].str.contains('.+H.+')],
        'I' : optimized_db[optimized_db['FRAG'].str.contains('.+I.+')],
        'L' : optimized_db[optimized_db['FRAG'].str.contains('.+L.+')],
        'K' : optimized_db[optimized_db['FRAG'].str.contains('.+K.+')],
        'M' : optimized_db[optimized_db['FRAG'].str.contains('.+M.+')],
        'F' : optimized_db[optimized_db['FRAG'].str.contains('.+F.+')],
        'P' : optimized_db[optimized_db['FRAG'].str.contains('.+P.+')],
        'S' : optimized_db[optimized_db['FRAG'].str.contains('.+S.+')],
        'T' : optimized_db[optimized_db['FRAG'].str.contains('.+T.+')],
        'V' : optimized_db[optimized_db['FRAG'].str.contains('.+V.+')],
        'W' : optimized_db[optimized_db['FRAG'].str.contains('.+W.+')],
        'Y' : optimized_db[optimized_db['FRAG'].str.contains('.+Y.+')],
    }

    return res_angles


def setup_databases(databases_paths: dict[str,str]) -> dict[str,dict[str,pd.DataFrame]]:
    """Setup the databases the movers can access during sampling of dihedral angles.

    Databases are read into memory, trimmed into only necessary columns and optimized
    to use the least amount of memory possible.

    Args:
        databases_paths:
            mapping of database_ids to filepaths where the specified databases are stored.

    Returns:
        databases:
            mapping of database_ids to mappings of aminoacid 1lettercode to their
            corresponding dihedral angle values in the optimized database.
    """
    databases = {}
    for db_id in databases_paths:

        # Load database into DataFrame
        original_database= read_database(database_path=databases_paths[db_id])

        # Keep only columns need for sampling
        trim_database(database=original_database,
                      columns_to_keep=GLOBAL_CONFIG['USED_DATABASE_COLNAMES'].values())

        # Optimize db datatypes for memory usage
        optimized_database = optimize_database(database=original_database)

        # Free up memory
        del original_database

        # Update databases with optimized version
        databases[db_id] = optimized_database

    return databases


def get_ss_bounds(secondary_structure: str) -> tuple[tuple[int,int],tuple[int,int]]:
    """Return the allowed range for the phi and psi angle values of a given secondary structure.

    Args:
        secondary_structure:
            identifier for a protein secondary structure.

    Returns:
        A tuple (phi_bounds,psi_bounds) where:
            phi_bounds: tuple with the lower and upper bounds for phi dihedral angle values
            for the secondary structure in question.
            psi_bounds: tuple with the lower and upper bounds for psi dihedral angle values
            for the secondary structure in question.
    """
    assert secondary_structure in ALLOWED_SECONDARY_STRUCTURE, ('Desired secondary structure '
                                                                'must be in '
                                                                f'{ALLOWED_SECONDARY_STRUCTURE}')

    if secondary_structure == 'alpha_helix':
        canonical_helix = GLOBAL_CONFIG['ALPHA_HELIX_CANON']
        phi_bounds = (canonical_helix[0] - 7, canonical_helix[0] + 7) #  ~ canonical_value ± 7°
        psi_bounds = (canonical_helix[1] - 7, canonical_helix[1] + 7)  # ~ canonical_value ± 7°
    elif secondary_structure == 'beta_strand':
        canonical_strand = GLOBAL_CONFIG['BETA_STRAND_CANON']
        phi_bounds = (canonical_strand[0] - 7, canonical_strand[0] + 7) # ~ canonical value ± 7°
        psi_bounds = (canonical_strand[1] - 7, canonical_strand[1] + 7) # ~ canonical value ± 7°
    return phi_bounds, psi_bounds
