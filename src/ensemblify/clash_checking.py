"""Simple module to re-check previously generated ensembles for steric clashes."""

# IMPORTS
## Standard Library Imports
import glob
import io
import os
import re
from concurrent.futures import ProcessPoolExecutor
from subprocess import CalledProcessError

## Third Party Imports
import yaml
from tqdm import tqdm

## Local Imports
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.generation.ensemble_utils.inputs_processing import process_input_pdb, register_input_clashes
from ensemblify.generation.ensemble_utils.pdb_processing import apply_pulchra_single, apply_rewrite_single, cleanup_pdbs
from ensemblify.utils import extract_pdb_info

# FUNCTIONS
def process_pulchra_output(
    sampled_pdb: str,
    pulchra_output_buffer: str,
    sampling_targets: dict[str,tuple[tuple[str,tuple[int,...],str,str]]] | None = None,
    input_clashes: list[tuple[str,str]] | None = None,
    ) -> list[str]:
    """Check if there are recorded steric clashes in given PULCHRA output.

    Clashes present in input structure (if provided) are not ignored.
    Clashes are only considered when at least one residue belongs to a sampled region (if those
    regions are provided).

    Args:
        sampled_pdb:
            filepath to .pdb file output from conformational sampling.
        pulchra_output_buffer:
            stdout from applying PULCHRA to the sampled .pdb structure.
        sampling_targets:
            mapping of chain identifiers to sampled residue numbers.
        input_clashes:
            clashes present in the sampling input structure, that will be ignored if
            present in the given PULCHRA output.
    Returns:
        clashes_sampled_pdb:
            list of clashes present in sampled .pdb file.
    """

    # Setup regular expressions
    regex_res_num = re.compile(r'([A-Z]{3}\[-?[0-9]+\][A-Z]{1,3}[0-9]*)') # LYS[309]OXT, LYS[309]01
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
    if sampling_targets is not None:
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
                if len(sampled_residues) > 0:
                    # Check if both residue numbers are part of sampled regions
                    if res1 in sampled_residues or res2 in sampled_residues:
                        if input_clashes is not None:
                            # Check if clash is not present in input clashes, in both 'directions'
                            if clash not in input_clashes and clash[::-1] not in input_clashes:
                                clashes_sampled_pdb.append(clash)
                        else:
                            clashes_sampled_pdb.append(clash)
                else:
                    if input_clashes is not None:
                        # Check if clash is not present in input clashes, in both 'directions'
                        if clash not in input_clashes and clash[::-1] not in input_clashes:
                            clashes_sampled_pdb.append(clash)
                    else:
                        clashes_sampled_pdb.append(clash)

    return clashes_sampled_pdb


def check_report_pdb_clashes(
    pdb2check: str,
    sampling_targets: dict[str,tuple[tuple[str,tuple[int,...],str,str]]] | None = None,
    input_clashes: list[tuple[str,str]] | None = None,
    ) -> tuple[str,list[str] | None]:
    """Check for steric clashes in a .pdb file, optionally considering sampling targets and clashes
    in the input structure.

    A steric clash is reported when the distance between any two non bonded atoms is less than two
    angstrom.

    Args:
        pdb2check:
            path to .pdb file to check for steric clashes.
        sampling_targets:
            mapping of chains to sampled regions following Ensemblify parameters style. If
            provided, clashes are only checked for in these regions. Defaults to None.
        input_clashes: 
            list of clashes detected in the ensemble generation input structure. If provided,
            clashes detailed here will be ignored if found in the .pdb to check (only in sampled
            regions, if sampling targets is provided). Defaults to None.
    Returns:
        A tuple[pdb2check,steric_clashes] where:
            pdb2check:
                the path to the sampled .pdb to check for clashes.
            steric_clashes:
                list of clashes present in sampled .pdb file or None, if PULCHRA erred.
    """
    # Rewrite .pdb into single chain and sequential numbering (for PULCHRA compatibility)
    rewrite_filename = apply_rewrite_single(pdb=pdb2check)

    try:
        # Apply PULCHRA
        rebuilt_filename,\
        pulchra_output_buffer = apply_pulchra_single(pulchra_path=GLOBAL_CONFIG['PULCHRA_PATH'],
                                                     pdb=rewrite_filename)
    except CalledProcessError:
        # If PULCHRA failed cleanup and return None
        cleanup_pdbs([rewrite_filename])
        return pdb2check, None

    # Process PULCHRA output, get a list of steric clashes
    steric_clashes = process_pulchra_output(sampled_pdb=pdb2check,
                                            pulchra_output_buffer=pulchra_output_buffer,
                                            sampling_targets=sampling_targets,
                                            input_clashes=input_clashes)

    # Cleanup temporary .pdb files
    cleanup_pdbs([rewrite_filename,rebuilt_filename])

    return pdb2check, steric_clashes


def check_steric_clashes(
    ensemble_dir: str,
    sampling_targets: str | dict[str,tuple[tuple[str,tuple[int,...],str,str]]] | None = None,
    input_structure: str | None = None,
) -> tuple[str,str]:
    """Check a generated ensemble for steric clashes, outputting clash reports.

    A directory is created inside the ensemble directory where clash reports (simple and detailed)
    will be stored, as well as any files output by processing the input structure (if provided).

    Args:
        ensemble_dir:
            path to directory where ensemble .pdb structures are stored.
        sampling_targets:
            mapping of chains to sampled regions following Ensemblify parameters style or path to
            this mapping in .yaml format. If a path is provided, it can either be solely the
            mapping or a full Ensemblify parameters file. If provided, clashes are only checked
            for in these regions. Defaults to None.
        input_structure:
            path to input structure used to generate the ensemble. If provided, steric clashes
            present in this structure (only in sampled regions, if sampling targets is provided)
            are ignored if they are detected in any of the sampled structures. Defaults to None.

    Returns:
        A tuple[clash_report,clash_report_detailed] where:
            clash_report:
                path to file with simplified ensemble clash report, i.e. total number of clashed
                structures and how many clashes were detected in each structure.
            clash_report_detailed:
                path to file with detailed ensemble clash report, i.e. how many clashes were
                detected in each structure and the atoms involved in the detected clash.
    """
    # Setup sampling targets
    if isinstance(sampling_targets,str):
        with open(sampling_targets,'r',encoding='utf-8-sig') as smp_targets_file:
            content = yaml.safe_load(smp_targets_file)
            try:
                smp_targets = content['targets']
            except KeyError:
                smp_targets = content
    else:
        smp_targets = sampling_targets

    # Create clash checking directory
    clash_checking_directory = os.path.join(ensemble_dir,'clash_checking')
    if not os.path.isdir(clash_checking_directory):
        os.mkdir(clash_checking_directory)

    # Grab .pdb files from ensemble directory
    pdbs_2_check = glob.glob(os.path.join(ensemble_dir,'*.pdb'))
    ensemble_size = len(pdbs_2_check)

    # Setup clash report files
    clash_report = os.path.join(clash_checking_directory,'clash_report.txt')
    clash_report_detailed = os.path.join(clash_checking_directory,'clash_report_detailed.txt')

    # Setup multiprocessing variables
    sampling_targets_all = [smp_targets] * ensemble_size
    if input_structure is not None:
        input_processing_directory = os.path.join(clash_checking_directory,'input_processing')
        if not os.path.isdir(input_processing_directory):
            os.mkdir(input_processing_directory)
        input_clashes_file, _  = process_input_pdb(faspr_path=GLOBAL_CONFIG['FASPR_PATH'],
                                                   pulchra_path=GLOBAL_CONFIG['PULCHRA_PATH'],
                                                   inputs_dir=input_processing_directory,
                                                   input_pdb=input_structure)
        input_clashes = register_input_clashes(input_clashes_file)
    else:
        input_clashes = None
    input_clashes_all = [input_clashes] * ensemble_size

    # Check ensemble for steric clashes, using multiprocessing
    with ProcessPoolExecutor() as ppe:
        results = list(tqdm(ppe.map(check_report_pdb_clashes,
                                    pdbs_2_check,
                                    sampling_targets_all,
                                    input_clashes_all),
                            desc='Checking ensemble for steric clashes... ',
                            total=ensemble_size)) # (pdb2check, steric_clashes)

    clashed_pdbs = 0
    erred_pdbs = []
    for pdb2check, steric_clashes in results:
        if steric_clashes is not None:
            nclashes = len(steric_clashes)
            if nclashes > 0:
                clashed_pdbs += 1
            with open(clash_report,'a',encoding='utf-8') as simplified_report:
                simplified_report.write(f'{pdb2check} : {nclashes}\n')
            with open(clash_report_detailed,'a',encoding='utf-8') as detailed_report:
                detailed_report.write(f'{pdb2check} : {nclashes}\n')
                for clash in steric_clashes:
                    detailed_report.write(f'{str(clash)}\n')
        else:
            erred_pdbs.append(pdb2check)

    if erred_pdbs:
        print('The following .pdb files could not be processed:')
        for erred_pdb in erred_pdbs:
            print(erred_pdb)

    with open(clash_report,'r+',encoding='utf-8') as simplified_report:
        content = simplified_report.read()
        simplified_report.seek(0,0)
        simplified_report.write(f'Total number of clashed pdbs: {clashed_pdbs}\n' + content)

    print(('Ensemble steric clash check has finished. Please check the reports created in the '
           'clash_checking directory.'))

    return clash_report, clash_report_detailed
