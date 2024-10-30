"""Auxiliary functions to process input parameters file and the input structure."""

# IMPORTS
## Standard Library Imports
import copy
import os
import re
import requests
import shutil

## Third Party Imports
import yaml

## Local Imports
from ensemblify.generation.ensemble_utils.pdb_processing import (
    apply_faspr_single, apply_pulchra_single, apply_rewrite_single,
    apply_restore_single
)
from ensemblify.utils import df_from_pdb, df_to_pdb

# CONSTANTS
VALID_PARAMS_TYPES = {
    'job_name' : str,
    'sequence' : str,
    'alphafold' : bool,
    'pae': str,
    'size' : int,
    'databases' : dict,
    'targets' : dict,
    'restraints' : dict,
    'core_amount' : int,
    'output_path' : str,
    'faspr_path' : str,
    'pulchra_path' : str,
    'scorefxn': dict,
    'minimizer' : dict,
    'sampler_params' : dict,
    'constraints' : dict,
    'constraints_violation': dict,
    'plddt_params': dict,
    'pae_params': dict
}

ADVANCED_PARAMS_DEFAULTS = {
    'alphafold': False,
    'pae': None,
    'restraints': {'ss_bias': None,
                    'contacts': None},
    'core_amount': os.cpu_count() - 1, 
    'output_path': os.getcwd(),
    'scorefxn': {'id': 'score0',
                    'weight': 1.0},
    'minimizer': {'id': 'dfpmin_armijo_nonmonotone',
                    'tolerance': 0.001,
                    'max_iters': 5000,
                    'finalcycles': 5},
    'constraints': {'weight': 1.0,
                    'stdev': 10.0,
                    'tolerance': 0.001},
    'constraints_violation': {'threshold': 0.015,
                                'maxres': 20},
    'plddt_params': {'threshold': 70,
                        'contiguous_res': 4},
    'pae_params': {'cutoff': 10.0,
                    'flatten_cutoff': 10.0,
                    'flatten_value': 10.0,
                    'weight': 1.0,
                    'tolerance': 0.001,
                    'adjacency_threshold': 8,
                    'plddt_scaling_factor': 1.0}
}

# FUNCTIONS
def process_input_pdb(
    faspr_path: str,
    pulchra_path: str,
    inputs_dir: str,
    input_pdb: str,
    ) -> tuple[str,str]:
    """Process input structure prior to sampling.
    
    Apply FASPR and Pulchra to the .pdb of our sampling input structure, outputting the resulting
    .pdb file. Allows us to start the sampling from a pdb file with less clashes, while taking note
    of the clashes still present in the input pdb so those are ignored later.

    Args:
        faspr_path:
            path to the FASPR executable (or its alias).
        pulchra_path:
            path to the PULCHRA executable (or its alias).
        inputs_dir:
            path to directory where files and .logs resulting from input processing will be stored.
        input_pdb:
            path to the .pdb input structure to process.

    Returns:
        A tuple (clashes_file,processed_pdb) where:
            clashes_file:
                .log file with the output of applying PULCHRA to our input structure
            processed_pdb:
                is the .pdb file resulting from applying FASPR and PULCHRA to our input structure.
    """

    # Save input pdb into input directory, if chainID is empty replace with 'A'
    pdb = os.path.join(inputs_dir,os.path.split(input_pdb)[1])
    input_pdb_df = df_from_pdb(input_pdb)
    if input_pdb_df['ChainID'][0] == ' ':
        input_pdb_df['ChainID'] = 'A'
    df_to_pdb(input_pdb_df,pdb)

    # Setup logs directory and log files
    LOGS_DIR = os.path.join(os.path.split(pdb)[0],'input_logs')
    if os.path.isdir(LOGS_DIR):
        shutil.rmtree(LOGS_DIR) # remove
    os.mkdir(LOGS_DIR) # make

    FASPR_LOG = os.path.join(LOGS_DIR, f'faspr_{os.path.split(pdb)[1][:-4]}_input.log')
    PULCHRA_LOG = os.path.join(LOGS_DIR, f'pulchra_{os.path.split(pdb)[1][:-4]}_input.log')
    clashes_file = os.path.join(LOGS_DIR, os.path.split(pdb)[1][:-4] + '_input_clashes.txt')

    if os.path.isfile(FASPR_LOG):
        os.remove(FASPR_LOG)
    if os.path.isfile(PULCHRA_LOG):
        os.remove(PULCHRA_LOG)
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


def register_input_clashes(input_clashes_file: str | None) -> list[tuple[str,str]]:
    """Register clashes in input structure to be ignored later.

    Args:
        input_clashes_file:
            path to the PULCHRA output file for the input structure.

    Returns:
        clashes_input:
            list of clashes in the input structure. Can be empty list if input_clashes_file
            is None. For example:
            [ ('ARG[277]', 'GLY[287]'),('ARG[706]', 'GLY[716]'), ... ].
    """
    # Get steric clashes present in input .pdb
    clashes_input = []

    if input_clashes_file is not None: # input structure could have been .txt
        with open(input_clashes_file,'r',errors='replace',encoding='utf-8-sig') as in_clashes:
            lines = in_clashes.readlines()

        for line in lines:
            if line.startswith('STERIC CONFLICT'):
                # Get a list containing the 2 residues participating in the clash
                clash_input = re.findall(re.compile(r'([A-Z]{3}\[-?[0-9]+\])'),line)
                if (clash_input not in clashes_input and
                    clash_input[::-1] not in clashes_input): # check for both directions
                    clashes_input.append(clash_input)

        clashes_input = [tuple(clash) for clash in clashes_input]

    return clashes_input


def get_protein_info(uniprot_accession: str) -> dict:
    """Get information about a protein from the AlphaFold Database using a given UniProt accession.

    Args:
        uniprot_accession:
            UniProt accession to use in request for AlphaFold's Database API.

    Returns:
        protein_info:
            information about the protein identified by the given UniProt accession, including
            links to its .pdb structure and .json PAE matrix.
    
    Adapted from:
        https://github.com/PDBeurope/afdb-notebooks/blob/main/AFDB_API.ipynb
    """
    api_endpoint = 'https://alphafold.ebi.ac.uk/api/prediction/'
    url = f'{api_endpoint}{uniprot_accession}'  # Construct the URL for API

    try:
        # Use a timeout to handle potential connection issues
        response = requests.get(url, timeout=10)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            protein_info = response.json()[0]
            return protein_info
        else:
            # Raise an exception for better error handling
            response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f'Error: {e}')


def download_from_url(url: str) -> str:
    """Download the contents of a url.

    Args:
        url:
            url path.
    
    Returns:
        content:
            content of url.
    
    Raises:
        RequestException with HTTPError if request timed out or encountered an issue.
    """
    # Download file
    try:
        # Use a timeout to handle potential connection issues
        response = requests.get(url, timeout=10)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            content = response.content
            return content
        else:
            # Raise an exception for better error handling
            response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f'Error: {e}')


def setup_ensemble_gen_params(input_params: dict, inputs_dir: str) -> tuple[str,str | None]:
    """Update sampling input parameters, store steric clashes present in input structure.
    
        If a UniProt accession is given in either the sequence of pae fields, replace it with the
    corresponding downloaded .pdb or .json file.

    - 'sequence' field is updated with the path to the processed input structure. If a UniProt
    accession is provided, replace it with the corresponding downloaded .pdb file.
    - 'pae' field is updated with the path to the downloaded .json PAE matrix file, if a UniProt
    accession is provided.
    - 'output_path' field is updated with the ensemble directory inside the created directory named 'job_name'.
    - File with the updated parameters is saved to the inputs directory.
    - Input structure is processed and any steric clashes present after processing are stored in a file so they can later be ignored on non-sampled regions of structures resulting from the sampling process.

    Args:
        input_params:
            parameters following the Ensemblify template.
        inputs_dir:
            path to directory where files and .logs resulting from input processing will be stored.

    Returns:
        A tuple (processed_parameters_path, input_clashes) where:
            processed_parameters_path:
                path to file where updated parameters are stored.
            input_clashes:
                path to file with PULCHRA output from processing input .pdb.
    """
    input_params_processed = copy.deepcopy(input_params)

    # Check if there are any UniProt Accessions in user inputs
    input_sequence = input_params_processed['sequence']
    input_pae_matrix = input_params_processed['pae']
    uniprotid_pattern = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|'
                                   r'[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
    uniprotid_match_sequence = re.search(uniprotid_pattern,input_sequence)
    uniprotid_match_paematrix = re.search(uniprotid_pattern,input_pae_matrix)

    # Update input .pdb structure if match found
    if uniprotid_match_sequence is not None:
        # Get protein info from UniProt accession
        protein_info = get_protein_info(uniprotid_match_sequence.group(0))

        # Extract PDB content
        pdb_url = protein_info.get('pdbUrl')
        pdb_content = download_from_url(pdb_url)

        # Write .pdb file
        input_pdb_path = os.path.join(input_params_processed['output_path'],
                                      input_params_processed['job_name'],
                                      'ensemble',
                                      'inputs',
                                      f'{uniprotid_match_sequence.group(0)}_structure.pdb')
        with open(input_pdb_path,'wb') as t:
            t.write(pdb_content)

        # Update sequence with .pdb
        input_params_processed['sequence'] = input_pdb_path

    # Update input PAE matrix if match found
    if uniprotid_match_paematrix is not None:
        try:
            # Extract PAE URL
            pae_url = protein_info.get('paeDocUrl')

        except NameError:
            # Get protein info from UniProt accession
            protein_info = get_protein_info(uniprotid_match_paematrix.group(0))

            # Extract PAE URL
            pae_url = protein_info.get('paeDocUrl')

        finally:
            # Extract PAE content
            pae_content = download_from_url(pae_url)

            # Write PAE .json file
            pae_matrix_path = os.path.join(input_params_processed['output_path'],
                                           input_params_processed['job_name'],
                                           'ensemble',
                                           'inputs',
                                           f'{uniprotid_match_paematrix.group(0)}_pae_matrix.json')
            with open(pae_matrix_path,'wb') as t:
                t.write(pae_content)

            # Update sequence with PAE
            input_params_processed['pae'] = pae_matrix_path

    # Process input .pdb through FASPR and PULCHRA
    input_sequence = input_params_processed['sequence']
    if input_sequence.endswith('.pdb'):
        # Apply FASPR and PULCHRA to input pdb
        input_clashes,\
        new_input_pdb = process_input_pdb(faspr_path=input_params_processed['faspr_path'],
                                          pulchra_path=input_params_processed['pulchra_path'],
                                          inputs_dir=inputs_dir,
                                          input_pdb=input_sequence)

        # Update input pdb with filepath to processed pdb
        input_params_processed['sequence'] = f'{new_input_pdb}'
    else:
        input_clashes = None # If given a .txt skip this step

    # Update output directory to ensemble folder
    OUTPUT_PATH = input_params_processed['output_path']
    JOB_NAME = input_params_processed['job_name']
    ENSEMBLE_DIR = os.path.join(OUTPUT_PATH,JOB_NAME,'ensemble')
    input_params_processed['output_path'] = ENSEMBLE_DIR

    # Save processed parameters in inputs directory
    processed_parameters_path = os.path.join(inputs_dir,'2_parameters_ensemble_gen.yaml')
    with open(processed_parameters_path,'w',encoding='utf-8') as f:
        for key in input_params_processed:
            f.write(f'{key}:\n')
            if isinstance(input_params_processed[key],dict):
                for nested_key in input_params_processed[key]:
                    f.write(f'  {nested_key}: ')
                    if input_params_processed[key][nested_key] is not None:
                        f.write(f'{input_params_processed[key][nested_key]}\n')
                    else:
                        f.write('\n')
            else:
                f.write(f'  {input_params_processed[key]}\n')

    return processed_parameters_path, input_clashes


def read_input_parameters(parameter_path: str) -> dict:
    """Read input parameters and assert the validity of its contents.

    Args:
        parameter_path:
            filepath to the parameter .yaml file.

    Returns:
        params:
            dictionary with all the parameters validated for correct python types.
    """
    with open(parameter_path,'r',encoding='utf-8-sig') as parameters_file:
        params = yaml.safe_load(parameters_file)

    for key in params:
        # Check parameter types
        if params[key] is None:
            try:
                params[key] = ADVANCED_PARAMS_DEFAULTS[key]
            except KeyError as e:
                raise AssertionError(f'{key} cannot be empty!') from e
        elif not isinstance(params[key],VALID_PARAMS_TYPES[key]):
            raise AssertionError(f'{key} must be of type {VALID_PARAMS_TYPES[key]} !')

        # Always leave one core free
        elif key =='core_amount' and params[key] >= os.cpu_count():
            params[key] = ADVANCED_PARAMS_DEFAULTS[key]

        # Take care of nested parameters
        elif key == 'databases':
            dbs = params[key]
            for db_id in dbs:
                assert isinstance(dbs[db_id],str), 'Database Path must be of type str !'
                assert (dbs[db_id].endswith('.pkl') or
                        dbs[db_id].endswith('.csv') or
                        dbs[db_id].endswith('.parquet')), ('Database must be a .pkl, '
                                                           '.csv or .parquet file!')

        elif key == 'restraints':
            restraints = params[key]
            for restraint_id in restraints:
                if restraints[restraint_id] is not None:
                    assert isinstance(restraints[restraint_id],list), ('Restraints must be of '
                                                                       'type list !')

        elif key == 'scorefxn':
            scorefxn_params = params[key]
            assert isinstance(scorefxn_params,dict), ('Score function parameters must be of '
                                                      'type dict!')
            for scorefxn_param in scorefxn_params:
                if scorefxn_param == 'id':
                    assert isinstance(scorefxn_params[scorefxn_param],str), ('Score function id '
                                                                             'must be of type str!')
                elif scorefxn_param == 'weight':
                    assert isinstance(scorefxn_params[scorefxn_param],float), ('Score function '
                                                                               'weight must be '
                                                                               'of type float!')

        elif key == 'minimizer':
            minimizer_params = params[key]
            assert isinstance(minimizer_params,dict), ('Minimizer parameters must '
                                                       'be of type dict!')
            for minimizer_param in minimizer_params:
                if minimizer_param == 'id':
                    assert isinstance(minimizer_params[minimizer_param],str), ('Minimizer id must '
                                                                               'be of type str!')
                elif minimizer_param == 'tolerance':
                    assert isinstance(minimizer_params[minimizer_param],float), ('Minimizer '
                                                                                 'tolerance must '
                                                                                 'be of type '
                                                                                 'float!')
                elif minimizer_param == 'max_iters':
                    assert isinstance(minimizer_params[minimizer_param],int), ('Minimizer maximum '
                                                                               'iterations must '
                                                                               'be of type int!')
                elif minimizer_param == 'finalcycles':
                    assert isinstance(minimizer_params[minimizer_param],int), ('Minimizer final '
                                                                               'cycles must be of '
                                                                               'type int!')

        elif key == 'sampler_params':
            sampler_params = params[key]
            for sampler_id in sampler_params:
                assert isinstance(sampler_params[sampler_id],dict), ('Sampler Parameters must be '
                                                                     'of type dict !')
                for smp_param_id in sampler_params[sampler_id]:
                    assert isinstance(sampler_params[sampler_id][smp_param_id],int), ('Sampler '
                                                                                      'Parameters '
                                                                                      'must be of '
                                                                                      'type int !')

        elif key == 'constraints':
            constraints_params = params[key]
            assert isinstance(constraints_params,dict), ('Constraints parameters must be of '
                                                         'type dict!')
            for constraints_param in constraints_params:
                if constraints_param == 'weight':
                    assert isinstance(constraints_params[constraints_param],float), ('Constraints '
                                                                                     'weight must '
                                                                                     'be of type '
                                                                                     'float!')
                elif constraints_param == 'stdev':
                    assert isinstance(constraints_params[constraints_param],float), ('Constraints '
                                                                                     'stdev must '
                                                                                     'be of type '
                                                                                     'float!')
                elif constraints_param == 'tolerance':
                    assert isinstance(constraints_params[constraints_param],float), ('Constraints '
                                                                                     'tolerance '
                                                                                     'must be of '
                                                                                     'type float!')

        elif key == 'constraint_violation':
            constraint_violation_params = params[key]
            assert isinstance(scorefxn_params,dict), ('Constraints violation parameters must be '
                                                      'of type dict!')
            for constraint_violation_param in constraint_violation_params:
                if constraint_violation_param == 'threshold':
                    assert isinstance(constraint_violation_params[constraint_violation_param],
                                      float), ('Constraints violation threshold must be of type '
                                               'float!')
                elif constraint_violation_param == 'maxres':
                    assert isinstance(constraint_violation_params[constraint_violation_param],
                                      int), ('Constraints violation maximum residues must be of '
                                             'type int!')

        elif key == 'plddt_params':
            plddt_params = params[key]
            assert isinstance(plddt_params,dict), ('plDDT parameters must be of type dict!')
            for plddt_param in plddt_params:
                if plddt_param == 'threshold':
                    assert isinstance(plddt_params[plddt_param],int), ('plDDT threshold must be '
                                                                       'of type int!')
                elif plddt_param == 'contiguous_res':
                    assert isinstance(plddt_params[plddt_param],int), ('plDDT contigous residues '
                                                                       'must be of type int!')

        elif key == 'pae_params':
            pae_params = params[key]
            assert isinstance(pae_params,dict), ('PAE parameters must be of type dict!')
            for pae_param in pae_params:
                if pae_param == 'cutoff':
                    assert isinstance(pae_params[pae_param],float), ('PAE cutoff must be of type '
                                                                     'float!')
                elif pae_param == 'flatten_cutoff':
                    assert isinstance(pae_params[pae_param],float), ('PAE flatten cutoff must be '
                                                                     'of type float!')
                elif pae_param == 'weight':
                    assert isinstance(pae_params[pae_param],float), ('PAE weight must be of type '
                                                                     'float!')
                elif pae_param == 'tolerance':
                    assert isinstance(pae_params[pae_param],float), ('PAE tolerance must be of '
                                                                     'type float!')
                elif pae_param == 'adjacency_threshold':
                    assert isinstance(pae_params[pae_param],int), ('PAE adjacency threshold must '
                                                                   'be of type int!')
                elif pae_param == 'plddt_scaling_factor':
                    assert isinstance(pae_params[pae_param],float), ('PAE pLDDT scaling factor '
                                                                     'must be of type float!')

    return params
