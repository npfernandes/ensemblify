"""
Auxiliary functions for custom Movers created from the PyRosetta Mover class.
"""

# IMPORTS
## Third Party Imports
import pandas as pd

## Local Imports
from ensemblify.config import GLOBAL_CONFIG

# CONSTANTS
DATABASE_OPTIMIZED_COL_DTYPES = { GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['OM1'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['OM2'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['OM3'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PH1'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PH2'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PH3'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PS1'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PS2'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['PS3'] : 'float32',
                                  GLOBAL_CONFIG['USED_DATABASE_COLNAMES']['FRG'] : 'category'}

# FUNCTIONS
def read_database(database_path: str) -> pd.DataFrame:
    """Read a database file into a pandas DataFrame.
    
    If possible, read only the desired set of columns into a pandas.DataFrame
    (depends on database file format).
    Currently supported database file extensions: .csv, .h5, .pkl.

    Args:
        database_path:
            filepath to database file.

    Returns:
        database: database as a pandas.DataFrame.
    """
    if database_path.endswith('.csv'): # pick columns, convert dtypes
        database = pd.read_csv(database_path,
                               usecols=list(GLOBAL_CONFIG['USED_DATABASE_COLNAMES'].values()),
                               dtype=DATABASE_OPTIMIZED_COL_DTYPES,
                               encoding='utf-8-sig')

    elif database_path.endswith('.h5'): # pick columns
        database = pd.read_hdf(database_path,
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

    # Now we break apart our full database into a dataFrame for each aminoacid residue
    A_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+A.+')]
    R_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+R.+')]
    N_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+N.+')]
    D_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+D.+')]
    C_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+C.+')]
    Q_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+Q.+')]
    E_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+E.+')]
    G_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+G.+')]
    H_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+H.+')]
    I_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+I.+')]
    L_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+L.+')]
    K_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+K.+')]
    M_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+M.+')]
    F_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+F.+')]
    P_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+P.+')]
    S_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+S.+')]
    T_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+T.+')]
    V_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+V.+')]
    W_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+W.+')]
    Y_angles_optimized = optimized_db[optimized_db['FRAGMENT'].str.contains('.+Y.+')]

    # Make sure we free up memory
    del optimized_db

    # Create the new database dictionary
    res_angles = {
        'A' : A_angles_optimized,
        'R' : R_angles_optimized,
        'N' : N_angles_optimized,
        'D' : D_angles_optimized,
        'C' : C_angles_optimized,
        'Q' : Q_angles_optimized,
        'E' : E_angles_optimized,
        'G' : G_angles_optimized,
        'H' : H_angles_optimized,
        'I' : I_angles_optimized,
        'L' : L_angles_optimized,
        'K' : K_angles_optimized,
        'M' : M_angles_optimized,
        'F' : F_angles_optimized,
        'P' : P_angles_optimized,
        'S' : S_angles_optimized,
        'T' : T_angles_optimized,
        'V' : V_angles_optimized,
        'W' : W_angles_optimized,
        'Y' : Y_angles_optimized,
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
    if secondary_structure == 'alpha_helix':
        canonical_helix = GLOBAL_CONFIG['ALPHA_HELIX_CANON']
        phi_bounds = (canonical_helix[0] - 7, canonical_helix[0] + 7) #  ~ canonical_value ± 7°
        psi_bounds = (canonical_helix[1] - 7, canonical_helix[1] + 7)  # ~ canonical_value ± 7°
    elif secondary_structure == 'beta_strand':
        canonical_strand = GLOBAL_CONFIG['BETA_STRAND_CANON']
        phi_bounds = (canonical_strand[0] - 7, canonical_strand[0] + 7) # ~ canonical value ± 7°
        psi_bounds = (canonical_strand[1] - 7, canonical_strand[1] + 7) # ~ canonical value ± 7°
    return phi_bounds, psi_bounds

