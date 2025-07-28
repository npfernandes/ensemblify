"""Auxiliary functions for processing .pdb files from sampling output."""

# IMPORTS
## Standard Library Imports
import io
import logging
import os
import re
import subprocess

## Third Party Imports
import ray

## Local Imports
from ensemblify.utils import cleanup_pdbs, df_from_pdb, df_to_pdb, extract_pdb_info

# FUNCTIONS
def check_clashes(
    sampled_pdb: str,
    pulchra_output_buffer: str,
    sampling_targets: dict[str,tuple[tuple[str,tuple[int,...],str,str]]],
    input_clashes: list[tuple[str,str]] | None,
    ) -> bool:
    """Check if there are recorded steric clashes in given PULCHRA output.

    Clashes present in input structure are not considered.
    Clashes are only considered when at least one residue belongs to a sampled region.

    Args:
        sampled_pdb (str):
            Filepath to .pdb file output from conformational sampling.
        pulchra_output_buffer (str):
            Stdout from applying PULCHRA to the sampled .pdb structure.
        sampling_targets (dict[str,tuple[tuple[str,tuple[int,...],str,str]]]):
            Mapping of chain identifiers to sampled residue numbers.
        input_clashes (list[tuple[str,str]], optional):
            Clashes present in the sampling input structure, that will be ignored if
            present in the given PULCHRA output.
    Returns:
        bool:
            True if the given PULCHRA output mentions clashes not already present in sampling
            input structure, False otherwise.
    """

    # Setup regular expressions
    regex_res_num = re.compile(r'([A-Z]{3}\[-?[0-9]+\])') # to find e.g. LYS[309]
    regex_num = re.compile(r'(-?[0-9]+)') # to find e.g. 309

    # Get info regarding chains and res nums
    pdb_info = extract_pdb_info(sampled_pdb) # (chain_letter, start_res, chain_size)
    ordered_chains_letters_sizes = [ (pdb_info[x][0],
                                      pdb_info[x][2]) for x in sorted(list(pdb_info.keys()))]

    # Get chain offsets according to size of previous chains
    chain_offsets = {}
    offset = 0
    for chain_letter,chain_size in ordered_chains_letters_sizes:
        chain_offsets[chain_letter] = offset
        offset += chain_size

    # Get sampled residues
    sampled_residues = set()
    for chain_letter,chain_size in ordered_chains_letters_sizes:
        offset = chain_offsets[chain_letter]
        all_target_res = [x[1] for x in sampling_targets[chain_letter]]

        # Get sampled residues
        for target_res in all_target_res:
            sampled_residues.update(range(target_res[0] + offset,target_res[-1] + 1 + offset))

    # Get PULCHRA output content
    with io.StringIO(pulchra_output_buffer) as pulchra_output_stream:
        clashes_file = pulchra_output_stream.readlines()

    # Find the clashes in sampled .pdb
    clashes_sampled_pdb = []
    for line in clashes_file:
        if line.startswith('STERIC CONFLICT'):

            # Get the 2 residues participating in the clash
            clash = tuple(re.findall(regex_res_num,line))

            # Check if this clash has not been recorded yet, in both 'directions'
            if clash not in clashes_sampled_pdb and clash[::-1] not in clashes_sampled_pdb:
                res1 = int(re.findall(regex_num,clash[0])[0])
                res2 = int(re.findall(regex_num,clash[1])[0])
                # Check if both residue numbers are part of sampled regions
                if res1 in sampled_residues or res2 in sampled_residues:
                    if input_clashes is not None:
                        # Check if clash is not present in input clashes, in both 'directions'
                        if clash not in input_clashes and clash[::-1] not in input_clashes:
                            clashes_sampled_pdb.append(clash)
                    else:
                        clashes_sampled_pdb.append(clash)

    clashed = len(clashes_sampled_pdb) > 0

    return clashed


def setup_logger(pdb: str, log_file: str) -> logging.Logger:
    """Setup a Logger object for this pdb, with output to log_file.

    Args:
        pdb (str):
            Path to .pdb file, will be the name of the logger.
        log_file (str):
            Filepath to log file.

    Returns:
        logger (logging.Logger):
            Logger object.
    """
    # Create logger object named pdb
    logger = logging.getLogger(pdb)

    # Setup logger level and formatter
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Assign logging directory
    ## Look for logs directory where .pdb file is stored
    LOGS_DIR = os.path.join(os.path.split(pdb)[0],'logs')
    if not os.path.isdir(LOGS_DIR):
        ## Look for an input_logs directory (we might be processing the input structure)
        LOGS_DIR = os.path.join(os.path.split(pdb)[0],'input_logs')
    if not os.path.isdir(LOGS_DIR):
        ## If processing files outside the pipeline simply create a .log file where .pdb is
        LOGS_DIR = os.path.split(pdb)[0]

    # Setup handler and add it to Logger
    file_handler = logging.FileHandler(os.path.join(LOGS_DIR,log_file))
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def apply_faspr_single(faspr_path: str, pdb: str) -> str | None:
    """Apply FASPR to a .pdb file. Log outcome.

    Args:
        faspr_path (str):
            Path to FASPR executable or its alias.
        pdb (str):
            Path to .pdb file.
    
    Returns:
        str | None:
            Path to .pdb file from FASPR output, with filename equal to
            input .pdb with the suffix '_faspr' added.
        
    Raises:
        subprocess.CalledProcessError if FASPR was not applied succesfully.
    """
    # Setup logging and output
    logger = setup_logger(pdb,'faspr.log')
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
    
    Necessary when a multichain .pdb will be input into Pulchra, as it does
    not support multiple chain structures.
    The output modified version has the _rewrite suffix added to its name.
    
    Args:
        pdb (str):
            Path to input .pdb file for conversion.
    
    Returns
        str:
            Path to modified .pdb. Filename is the same as the input, with
            _rewrite suffix added.
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


def apply_pulchra_single(pulchra_path: str, pdb: str) -> tuple[str,str] | tuple[None,None]:
    """Apply PULCHRA to a .pdb file. Log outcome.

    Args:
        pulchra_path (str):
            Path to PULCHRA executable or its alias.
        pdb (str):
            Path to .pdb file.

    Returns:
        tuple[str,str] | tuple[None,None]:
            rebuilt_filename (str | None):
                Path to PULCHRA output structure. Same filename as input .pdb
                with added .rebuilt suffix.
            pulchra_output (str | None):
                PULCHRA stdout used for later clash checking.
        
    Raises:
        subprocess.CalledProcessError if PULCHRA was not applied sucessfully.
    """
    # Setup logging and output
    logger = setup_logger(pdb,'pulchra.log')
    logger.info('Applying PULCHRA')
    rebuilt_filename = pdb[:-4] + '.rebuilt.pdb'

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
    """Restore chain, residue number and B-Factor info to pdb from reference pdb.
        
    Restore chain, residue numbering and B-Factor information to a post-Pulchra .pdb
    file, following the information in a reference .pdb file (either the first output
    of sampling process or the sampling input .pdb).

    Args:
        pdb (str):
            Path to the PULCHRA output .pdb structure (ending in .rebuilt suffix).
        reference_pdb (str):
            Path to .pdb file to use as reference for restoring the chains and
            residue numbering.
    
    Returns:
        str:
            Path to the .pdb structure with restored chain and residue numbering.
            Filename matches that of input, with _restores suffix added.
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


@ray.remote(num_returns=1,
            num_cpus=1)
def process_pdb(
    sampled_pdb: str | None,
    faspr_path: str,
    pulchra_path: str,
    input_clashes: list[tuple[str,str]],
    sampling_targets: dict[str,tuple[tuple[str,tuple[int,...],str,str]]],
    valid_pdbs_dir: str,
    goal_ensemble_size: int,
    ) -> bool | None:
    """Repack the side-chains and check for steric clashes in a sampled .pdb structure.
    
    Side-chain repacking is done by passing the structure through FASPR.
    The resulting .pdb file is then rewritten into a single chain with sequential residue numbering
    before being passed into PULCHRA, as it does not support multi-chain structures.
    Clash checking is done by passing the structure through PULCHRA and checking its output.
    If no clashes are present, the resulting .pdb file has its chain and residue numbering
    information restored to its original status.

    Args:
        sampled_pdb (str | None):
            Sampled .pdb structure, unprocessed. If None, processing is cancelled.
        faspr_path (str):
            Path to FASPR executable or its alias.
        pulchra_path (str):
            Path to PULCHRA executable or its alias.
        input_clashes (list[tuple[str,str]]):
            List of clashes present in the input structure that, if present, will be ignored.
        sampling_targets (dict[str,tuple[tuple[str,tuple[int,...],str,str]]]):
            Mapping of chain letters to target regions for sampling.
        valid_pds_dir (str):
            Path to directory where valid structures will be output.
        goal_ensemble_size (int):
            If the number of structures in valid pdbs directory is ever greater than this value
            do not write any more structures into the directory.
    
    Returns:
        bool | None:
            True if the sampled .pdb structure has steric clashes, False otherwise.
            None if an error occured.
    """

    if isinstance(sampled_pdb,type(None)):
        return None

    # If sampled_pdb does not exist, cancel processing and return it as clashed (invalid)
    if not os.path.isfile(sampled_pdb):
        return True

    # Try to apply pdb processing
    try:
        # Apply FASPR
        faspr_filename = apply_faspr_single(faspr_path=faspr_path,
                                            pdb=sampled_pdb)
    except subprocess.CalledProcessError:
        # If FASPR failed cancel processing and return it as clashed (invalid)
        cleanup_pdbs([sampled_pdb])
        return True

    # Rewrite .pdb into single chain and sequential numbering
    rewrite_filename = apply_rewrite_single(pdb=faspr_filename)

    try:
        # Apply PULCHRA
        rebuilt_filename, pulchra_output_buffer = apply_pulchra_single(pulchra_path=pulchra_path,
                                                                       pdb=rewrite_filename)
    except subprocess.CalledProcessError:
        # If PULCHRA failed cancel processing and return it as clashed (invalid)
        cleanup_pdbs([sampled_pdb,faspr_filename,rewrite_filename])
        return True

    # Check for steric clashes
    clashed = check_clashes(sampled_pdb=sampled_pdb,
                            pulchra_output_buffer=pulchra_output_buffer,
                            sampling_targets=sampling_targets,
                            input_clashes=input_clashes)

    if not clashed:
        # Restore multichain and residue numbering info
        restored_pdb = apply_restore_single(pdb=rebuilt_filename,
                                            reference_pdb=sampled_pdb)

        # Save final valid .pdb file to valid pdbs directory
        i = len(os.listdir(valid_pdbs_dir)) + 1
        if i <= goal_ensemble_size:
            with open(restored_pdb,'r',encoding='utf-8-sig') as f:
                pdb_data = f.read()
            valid_pdb_path = os.path.join(valid_pdbs_dir,f'{i}.pdb')
            with open(valid_pdb_path,'w',encoding='utf-8') as t:
                t.write(pdb_data)

        # Remove restored .pdb
        cleanup_pdbs([restored_pdb])

    # Cleanup temporary .pdb files
    cleanup_pdbs([sampled_pdb,faspr_filename,rewrite_filename,rebuilt_filename])

    return clashed

