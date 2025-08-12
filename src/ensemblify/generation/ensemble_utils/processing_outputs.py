"""Auxiliary functions for processing files from sampling output."""

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
from ensemblify.utils import cleanup_pdbs, extract_pdb_info

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

