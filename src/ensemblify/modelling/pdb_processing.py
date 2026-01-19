"""Auxiliary functions for processing .pdb files."""

# IMPORTS
## Standard Library Imports
import logging
import os
import subprocess

## Local Imports
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.utils import df_from_pdb, df_to_pdb, extract_pdb_info

# FUNCTIONS
def _setup_logger(pdb: str, log_file: str) -> logging.Logger:
    """Setup a Logger object for this pdb, with output to log_file.

    Args:
        pdb (str):
            Path to .pdb file, will be the name of the logger.
        log_file (str):
            Filepath to log file.

    Returns:
        logging.Logger:
            Logger object.
    """
    # Create logger object named pdb
    logger = logging.getLogger(pdb)

    # Setup logger level and formatter
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Assign logging directory
    ## Look for logs directory where .pdb file is stored
    LOGS_DIR = os.path.join(os.path.dirname(pdb),'logs')
    if not os.path.isdir(LOGS_DIR):
        ## Look for an input_logs directory (we might be processing the input structure)
        LOGS_DIR = os.path.join(os.path.dirname(pdb),'input_logs')
    if not os.path.isdir(LOGS_DIR):
        ## If processing files outside the pipeline simply create a .log file where .pdb is
        LOGS_DIR = os.path.dirname(pdb)

    # Setup handler and add it to Logger
    file_handler = logging.FileHandler(os.path.join(LOGS_DIR,log_file))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def apply_faspr_single(faspr_path: str,pdb: str) -> str | None:
    """Apply FASPR to a .pdb file and log the outcome.

    Args:
        faspr_path (str):
            Path to FASPR executable.
        pdb (str):
            Path to .pdb file.
    
    Returns:
        str | None:
            Path to .pdb file from FASPR output, with filename equal to input .pdb with the suffix
            '_faspr' added. None when this function raises an exception.

    Raises:
        subprocess.CalledProcessError:
            If FASPR was not applied succesfully.
    """
    # Setup logging and output
    logger = _setup_logger(pdb,'faspr.log')
    logger.info('Applying FASPR')
    faspr_output_filepath = f'{os.path.splitext(pdb)[0]}_faspr.pdb'

    try:
        # Run FASPR, raise CalledProcessError if return code is not 0
        subprocess.run([faspr_path, '-i', pdb, '-o', faspr_output_filepath],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE,
                       check=True)
        logger.info('Done FASPR')
        return faspr_output_filepath

    except subprocess.CalledProcessError as e:
        logger.error(f'FASPR failed with return code {e.returncode}')
        logger.error(f'Stderr: {e.stderr.decode()}')
        raise e


def apply_rewrite_single(pdb: str) -> str:
    """Convert a .pdb file into single chain with sequential numbering.
    
    Necessary when a multichain .pdb will be input into PULCHRA, as it does not support multiple
    chain structures. The output modified version has the _rewrite suffix added to its name.
    
    Args:
        pdb (str):
            Path to input .pdb file for conversion.
    
    Returns
        str:
            Path to modified .pdb. Filename is the same as the input, with _rewrite suffix added.
    """
    # Get information in .pdb as DataFrame
    df = df_from_pdb(pdb)
    df_after_faspr = df.copy(deep=True)

    # Sequential resnumbers
    new_res_nums = []
    offset = 0
    curr_chain_id = df_after_faspr.loc[0,'ChainID'] # start with 'A'
    curr_chain_size = 0
    for i, resnum in enumerate(df_after_faspr['ResidueNumber']):
        if df_after_faspr.loc[i,'ChainID'] != curr_chain_id:
            offset += curr_chain_size
            curr_chain_id = df_after_faspr.loc[i,'ChainID']
            curr_chain_size = 0
        new_res_nums.append(resnum + offset)
        if i > 0:
            if resnum != df_after_faspr.loc[i-1,'ResidueNumber']:
                curr_chain_size += 1
        else:
            curr_chain_size += 1
    df_after_faspr['ResidueNumber'] = new_res_nums

    # Change all chains to A
    df_after_faspr['ChainID'] = 'A'

    rewrite_filename = pdb[:-4] + '_rewrite.pdb'
    df_to_pdb(df_after_faspr,rewrite_filename)
    return rewrite_filename


def apply_pulchra_single(pulchra_path: str,pdb: str) -> tuple[str,str] | tuple[None,None]:
    """Apply PULCHRA to a .pdb file and log the outcome.

    Args:
        pulchra_path (str):
            Path to PULCHRA executable.
        pdb (str):
            Path to .pdb file.

    Returns:
        tuple[str,str] | tuple[None,None]:
            str | None:
                Path to PULCHRA output structure. Same filename as input .pdb with added .rebuilt
                suffix. None when this function raises an exception.
            str | None:
                PULCHRA stdout used for later clash checking. None when this function raises an
                exception.
        
    Raises:
        subprocess.CalledProcessError:
            If PULCHRA was not applied sucessfully.
    """
    # Setup logging and output
    logger = _setup_logger(pdb,'pulchra.log')
    logger.info('Applying PULCHRA')
    rebuilt_filename = f'{os.path.splitext(pdb)[0]}.rebuilt.pdb'

    try:
        # Run PULCHRA, raise CalledProcessError if return code is not 0
        result = subprocess.run([pulchra_path, '-q', '-e', '-v', pdb],
                                capture_output=True,
                                check=True)

        logger.info('Done PULCHRA')
        # Capture PULCHRA stdout
        pulchra_output = result.stdout.decode('utf-8', errors='replace').strip()

        return rebuilt_filename, pulchra_output

    except subprocess.CalledProcessError as e:
        logger.error(f'PULCHRA failed with return code {e.returncode}')
        logger.error(f'Stderr: {e.stderr.decode()}')
        raise e


def apply_restore_single(pdb: str, reference_pdb: str) -> str:
    """Restore chain, residue number and B-Factor info to .pdb from reference .pdb.
        
    Restore chain, residue numbering and B-Factor information to a post-Pulchra .pdb
    file, following the information in a reference .pdb file (either the first output
    of sampling process or the sampling input .pdb).

    Args:
        pdb (str):
            Path to the PULCHRA output .pdb structure (ending in .rebuilt suffix).
        reference_pdb (str):
            Path to .pdb file to use as reference for restoring the chains and residue numbering.
    
    Returns:
        str:
            Path to the .pdb structure with restored chain and residue numbering.
            Filename matches that of input, with _restored suffix added.
    """

    # Get information from reference pdb
    pdb_info = extract_pdb_info(reference_pdb)

    # Assign each chain a residue number range
    ordered_chains = [pdb_info[x][0] for x in sorted(pdb_info.keys(),reverse=True)]
    chain_res_ranges = {}
    for i,chain_number in enumerate(sorted(pdb_info.keys(),reverse=True)):
        chain_letter, chain_start, chain_size = pdb_info[chain_number]
        try:
            prev_chain_range = chain_res_ranges[ordered_chains[i-1]]
            chain_res_ranges[chain_letter] = list(range(prev_chain_range[-1] + 1 ,
                                                 prev_chain_range[-1] + chain_size + 1 ))
        except (KeyError,IndexError):
            chain_res_ranges[chain_letter] = list(range(chain_start,chain_start + chain_size))

    # Grab our rebuilt pdb file as DataFrame
    df = df_from_pdb(pdb)
    df_copy = df.copy(deep=True)

    # Restore chain IDs
    for i,res_num in enumerate(df_copy['ResidueNumber']):
        for chain_letter, chain_residues in chain_res_ranges.items():
            if res_num in chain_residues:
                df_copy.loc[i,'ChainID'] = chain_letter

    # Restore residue numbers
    restored_res_nums = []
    offset = 0
    offset_history = []
    curr_chain = df_copy.loc[0,'ChainID']
    for i,res_num in enumerate(df_copy['ResidueNumber']):
        if df_copy.loc[i,'ChainID'] != curr_chain:
            # When we change chain subtract an offset from current residue number
            curr_chain = df_copy.loc[i,'ChainID']
            offset = len(set(restored_res_nums))
            offset_history.append(offset)
        restored_res_nums.append(res_num - sum(offset_history))
    df_copy['ResidueNumber'] = restored_res_nums

    # Update occupancy to 1
    df_copy['Occupancy'] = 1.0

    # Restore B-Factor
    ref_df = df_from_pdb(reference_pdb)
    for i,res_num in enumerate(df_copy['ResidueNumber']):
        res_num_entries = ref_df.loc[ref_df['ResidueNumber'] == res_num]
        ca_entry = res_num_entries[res_num_entries['AtomName'] == 'CA']
        ca_b_factor = ca_entry['B-Factor']

        res_num_entries_restored = df_copy.loc[df_copy['ResidueNumber'] == res_num]
        df_copy.loc[res_num_entries_restored.index,'B-Factor'] = ca_b_factor.values[0]

    # Save restored .pdb file
    restored_filename = pdb[:-4] + '_restored.pdb'
    df_to_pdb(df_copy,restored_filename)

    return restored_filename


def process_pdb_structure(
    pdb: str,
    faspr_path: str | None = None,
    pulchra_path: str | None = None,
    ) -> tuple[str,str]:
    """Apply FASPR and PULCHRA to PDB structure.
    
    FASPR is used for side-chain repacking and PULCHRA is used for backbone optimization and steric
    clash checking.
    The processed PDB file and a PULCHRA clashes report file are created in the same directory as
    the provided file.

    Args:
        fused_pdb (str):
            Path to the .pdb file with the PDB structure to process.
        faspr_path (str):
            Path to the FASPR executable. Defaults to the path set in Ensemblify global configuration.
        pulchra_path (str):
            Path to the PULCHRA executable. Defaults to the path set in Ensemblify global configuration. 

    Returns:
        tuple[str,str]:
            str:
                Path to the .log file with the output of applying PULCHRA to the provided structure.
            str:
                Path to the .pdb file resulting from applying FASPR and PULCHRA to the provided structure.
    """
    # Setup FASPR and PULCHRA paths
    if faspr_path is None:
        faspr_path = GLOBAL_CONFIG['FASPR_PATH']

    if pulchra_path is None:
        pulchra_path = GLOBAL_CONFIG['PULCHRA_PATH']

    clashes_file = os.path.join(os.path.dirname(pdb),
                                f'{os.path.splitext(os.path.basename(pdb))[0]}_clashes.txt')

    if os.path.isfile(clashes_file):
        os.remove(clashes_file)

    # Process input pdb
    FASPR_PDB = apply_faspr_single(faspr_path=faspr_path,
                                   pdb=pdb)

    FASPR_REWRITE_PDB = apply_rewrite_single(pdb=FASPR_PDB)

    REBUILT_PDB, clashes_pulchra_buffer = apply_pulchra_single(pulchra_path=pulchra_path,
                                                               pdb=FASPR_REWRITE_PDB)

    # Write Pulchra output to log file
    with open(clashes_file, 'a', errors='ignore',encoding='utf-8') as log_file:
        log_file.write(f'{clashes_pulchra_buffer}\n')
        log_file.flush()

    REBUILT_RESTORED_PDB = apply_restore_single(pdb=REBUILT_PDB,
                                                reference_pdb=pdb)

    # Rename processed input_pdb from ..._rebuilt_restored.pdb to ..._processed.pdb
    processed_pdb = f'{REBUILT_RESTORED_PDB[:-35]}_processed.pdb'
    os.rename(REBUILT_RESTORED_PDB,processed_pdb)

    # cleanup intermediate files
    os.remove(FASPR_PDB)
    os.remove(FASPR_REWRITE_PDB)
    os.remove(REBUILT_PDB)

    return clashes_file, processed_pdb
