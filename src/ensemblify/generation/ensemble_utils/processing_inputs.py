"""Auxiliary functions to process input parameters file and input structure."""

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
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.modelling import fuse_structures
from ensemblify.modelling.pdb_processing import (
    apply_faspr_single,
    apply_pulchra_single,
    apply_rewrite_single,
    apply_restore_single
)
from ensemblify.utils import cif_to_pdb, df_from_pdb, df_to_pdb

# CONSTANTS
VALID_PARAMS_TYPES = {
    'job_name': str,
    'sequence': str | dict,
    'alphafold': bool,
    'pae': str,
    'size': int,
    'databases': dict,
    'targets': dict,
    'restraints': dict,
    'core_amount': int,
    'output_path': str,
    'faspr_path': str,
    'pulchra_path': str,
    'scorefxn': dict,
    'minimizer': dict,
    'variability': dict,
    'sampler_params': dict,
    'constraints': dict,
    'constraints_violation': dict,
    'plddt_params': dict,
    'pae_params': dict
}

ADVANCED_PARAMS_DEFAULTS = {
    'alphafold': False,
    'pae': None,
    'restraints': {'ss_bias': None,
                   'contacts': None},
    'output_path': os.getcwd(),
    'faspr_path': None,
    'pulchra_path': None,
    'core_amount': os.cpu_count() - 1, 
    'scorefxn': {'id': 'score0',
                 'weight': 1.0},
    'minimizer': {'id': 'dfpmin_armijo_nonmonotone',
                  'tolerance': 0.001,
                  'max_iters': 5000,
                  'finalcycles': 5},
    'variability': {'variance': 0.10},
    'sampler_params': {'MC': {'temperature': 200,
                              'max_loops': 200}},
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
                   'plddt_scaling_factor': 30.0}
}

DISORDER_DATA_KEYS = [
    'curated-disorder-merge',
    'prediction-disorder-mobidb_lite',
]

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
        faspr_path (str):
            Path to the FASPR executable.
        pulchra_path (str):
            Path to the PULCHRA executable.
        inputs_dir (str):
            Path to directory where files and .logs resulting from input processing will be stored.
        input_pdb (str):
            Path to the .pdb input structure to process.

    Returns:
        tuple[str,str]:
            clashes_file (str):
                Path to the .log file with the output of applying PULCHRA to our input structure.
            processed_pdb (str):
                Path to the .pdb file resulting from applying FASPR and PULCHRA to our input structure.
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

    # Setup FASPR and PULCHRA paths
    if faspr_path is None:
        faspr_path = GLOBAL_CONFIG['FASPR_PATH']

    if pulchra_path is None:
        pulchra_path = GLOBAL_CONFIG['PULCHRA_PATH']

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
        input_clashes_file (str, optional):
            Path to the PULCHRA output file for the input structure. If None, no clashes are
            registered.

    Returns:
        list[tuple[str,str]]:
            List of clashes in the input structure. Can be empty list if input_clashes_file
            is None. For example:

            [ ('ARG[277]', 'GLY[287]'),('ARG[706]', 'GLY[716]'), ... ]

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
        uniprot_accession (str):
            UniProt accession to use in request for AlphaFold's Database API.

    Returns:
        dict:
            Information about the protein identified by the given UniProt accession, including
            links to its .pdb structure and .json PAE matrix.

    Raises:
        RequestException with HTTPError if request timed out or encountered an issue.

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


def _download_from_url(url: str) -> str:
    """Download the contents of a URL.

    Args:
        url (str):
            URL path.
    
    Returns:
        str:
            Content of URL.
    
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
    
    If a UniProt accession is given in either the sequence or PAE fields, replace it with the
    corresponding downloaded .pdb or .json file.

    - 'sequence' field is updated with the path to the processed input structure. If a UniProt
      accession is provided, replace it with the corresponding downloaded .pdb file.
    - 'pae' field is updated with the path to the downloaded .json PAE matrix file, if a UniProt
      accession is provided.
    - 'output_path' field is updated with the ensemble directory inside the created directory named
      'job_name'.
    - File with the updated parameters is saved to the inputs directory.
    - Input structure is processed and any steric clashes present after processing are stored in a
      file so they can later be ignored on non-sampled regions of structures resulting from the
      sampling process.


    Args:
        input_params (dict):
            Parameters following the Ensemblify template.
        inputs_dir (str):
            Path to directory where files and .logs resulting from input processing will be stored.

    Returns:
        tuple[str,str] | None:
            processed_parameters_path (str):
                Path to file where updated parameters are stored.
            input_clashes (str):
                Path to file with PULCHRA output from processing input .pdb.
    """
    input_params_processed = copy.deepcopy(input_params)
    input_sequence = input_params_processed['sequence']
    input_pae_matrix = input_params_processed['pae']

    # Check if there are any UniProt Accessions in user inputs
    uniprotid_pattern = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|'
                                   r'[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
    uniprotid_match_sequence = re.search(uniprotid_pattern,input_sequence)
    if input_pae_matrix is not None:
        uniprotid_match_paematrix = re.search(uniprotid_pattern,input_pae_matrix)
    else:
        uniprotid_match_paematrix = None

    # Update input .pdb structure if match found
    if uniprotid_match_sequence is not None:
        # Get protein info from UniProt accession
        protein_info = get_protein_info(uniprotid_match_sequence.group(0))

        # Extract PDB content
        pdb_url = protein_info.get('pdbUrl')
        pdb_content = _download_from_url(pdb_url)

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
            protein_info = get_alphafold_data(uniprotid_match_paematrix.group(0))

            # Extract PAE URL
            pae_url = protein_info.get('paeDocUrl')

        finally:
            # Extract PAE content
            pae_content = _download_from_url(pae_url)

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

    # Convert input PDBx/mmCIF to PDB if applicable
    input_sequence = input_params_processed['sequence']
    if input_sequence.endswith('.cif'):
        pdb_path = os.path.splitext(input_sequence)[0] + '.pdb'
        cif_to_pdb(cif_path=input_sequence,
                   pdb_path=pdb_path)
        input_params_processed['sequence'] = pdb_path

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
        input_clashes = None # If given a .txt or sequence string skip this step

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
        parameter_path (str):
            Path to the parameter .yaml file.

    Returns:
        dict:
            Dictionary with all the parameters validated for correct Python types.
    """
    with open(parameter_path,'r',encoding='utf-8-sig') as parameters_file:
        params = yaml.safe_load(parameters_file)

    # Fill in empty parameters
    for key,default_value in ADVANCED_PARAMS_DEFAULTS.items():
        try:
            params[key]
        except KeyError:
            params[key] = default_value

    # Check all parameters
    for key in params:
        # Check empty parameters
        if params[key] is None:
            try:
                params[key] = ADVANCED_PARAMS_DEFAULTS[key]
            except KeyError as e:
                raise AssertionError(f'{key} cannot be empty!') from e
        
        # Check parameter types
        elif not isinstance(params[key],VALID_PARAMS_TYPES[key]):
            raise AssertionError(f'{key} must be of type {VALID_PARAMS_TYPES[key]} !')

        # Apply modelling module if applicable
        elif key == 'sequence' and isinstance(params[key],dict):
            job_name = params['job_name']
            output_name = f'{job_name}_fused'
            output_dir = os.path.join(params['output_path'],
                                      f'{job_name}_fusion')

            # Fuse structure and update 'sequence' accordingly
            _, params[key] = fuse_structures(input_fastas=params['sequence']['fastas'],
                                             input_pdbs=params['sequence']['pdbs'],
                                             output_name=output_name,
                                             output_dir=output_dir)

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

        elif key == 'variability':
            variability_params = params[key]
            assert isinstance(variability_params,dict), ('Variability parameters must be of '
                                                           'type dict!')
            for variability_param in variability_params:
                if variability_param == 'variance':
                    assert isinstance(variability_params[variability_param],float), ('Variance '
                                                                                     'must be of '
                                                                                     'type '
                                                                                     'float!')

        elif key == 'sampler_params':
            sampler_params = params[key]
            assert isinstance(sampler_params,dict), ('Samplers parameters must be of type dict!')
            for sampler_id in sampler_params:
                assert isinstance(sampler_params[sampler_id],dict), (f'{sampler_id} sampler '
                                                                     'parameters must be of '
                                                                     'type dict !')
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


def get_mobidb_data(uniprot_accession: str) -> dict:
    """Get disorder data for a protein from MobiDB using a given UniProt accession.

    Args:
        uniprot_accession (str):
            UniProt accession to use in request for MobiDB's API.

    Returns:
        dict:
            Intrinsic disorder data for the protein identified by the given UniProt accession,
            including experimentally confirmed and predicted disordered regions.
            Returns None if no disorder data is found or if the request fails.

    Raises:
        RequestException with HTTPError if request timed out or encountered an issue.
    """
    api_endpoint = 'https://mobidb.org/api/download_page'
    params = {
        'acc': uniprot_accession,
        'limit': 10,
        'format': 'json',
    }

    try:
        # Use a timeout to handle potential connection issues
        response = requests.get(api_endpoint, params=params, timeout=10)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            data = response.json()
            entries = data.get('data', [])
            if not entries:
                raise RuntimeError(f'Error: No MobiDB entry found for accession {uniprot_accession}')

            entry = entries[0]
            disorder_data = {}
            for disorder_key in DISORDER_DATA_KEYS:
                if disorder_key in entry:
                    disorder_data[disorder_key] = entry[disorder_key]
            if not disorder_data:
                raise RuntimeError(f'Error: No disorder annotations found in MobiDB for {uniprot_accession}')

            return disorder_data

        else:
            # Raise an exception for better error handling
            response.raise_for_status()

    except requests.exceptions.RequestException as e:
        print(f'Error: {e}')


def format_mobidb_report(
    uniprot_accession: str,
    disorder_data: dict,
    ) -> str:
    """Format a human-readable report of disordered data from MobiDB.

    Args:
        uniprot_accession (str):
            UniProt accession used to extract the MobiDB data.
        disorder_data (dict):
            Extracted data from MobiDB.

    Returns:
        str:
            Formatted multi-line report string.
    """
    lines = [
        f'\nMobiDB report for {uniprot_accession} (Most to least reliable)',
        '=' * 68,
    ]

    for data_type, region_data in disorder_data.items():
        regions = region_data.get('regions', [])
        fraction = region_data.get('content_fraction', 0.0)
        count = region_data.get('content_count', 0)

        lines.append(f'Evidence Type       : {data_type}')
        lines.append(f'Disordered residues : {count} ({fraction:.1%} of sequence)')
        lines.append(f'Regions             : {len(regions)}')

        for i, (start, end) in enumerate(regions, start=1):
            lines.append(f'    [{i}] residues {start}-{end}  (length {end - start + 1})')

        lines.append('')

    return '\n'.join(lines)


def get_ted_data(uniprot_accession: str) -> dict:
    """Get structured domain annotations for a protein from The Encyclopedia of Domains (TED).

    Args:
        uniprot_accession (str):
            UniProt accession used to query TED.

    Returns:
        dict:
            Dictionary containing parsed domain annotations.
            Returns None if request fails or no domains are found.

    Raises:
        RequestException with HTTPError if request timed out or encountered an issue.
    """
    api_endpoint = f'https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_accession}'

    try:
        response = requests.get(api_endpoint, timeout=10)

        if response.status_code == 200:
            data = response.json()
            entries = data.get('data', [])

            if not entries:
                raise RuntimeError(
                    f'Error: No TED domain entries found for accession {uniprot_accession}'
                )

            domain_data = []

            for entry in entries:
                chopping = entry.get('chopping', '')
                segments = []
                if chopping:
                    for part in chopping.split(','):
                        if '_' in part:
                            domains = part.split('_')
                            for domain in domains:
                                start, end = domain.split('-')
                                segments.append((int(start), int(end)))
                        elif '-' in part:
                            start, end = part.split('-')
                            segments.append((int(start), int(end)))

                domain_data.append({
                    'ted_id': entry.get('ted_id'),
                    'consensus_level': entry.get('consensus_level'),
                    'segments': segments,
                    'nres_domain': entry.get('nres_domain'),
                    'plddt': entry.get('plddt'),
                    'cath_label': entry.get('cath_label'),
                })

            if not domain_data:
                raise RuntimeError(
                    f'Error: No valid domain annotations found for {uniprot_accession}'
                )

            return {'domains': domain_data}

        else:
            response.raise_for_status()

    except requests.exceptions.RequestException as e:
        print(f'Error: {e}')

def format_ted_report(uniprot_accession: str, domain_data: dict) -> str:
    """Format a human-readable report of TED domain annotations.

    Args:
        uniprot_accession (str):
            UniProt accession used to extract the TED data.
        disorder_data (dict):
            Extracted data from TED.

    Returns:
        str:
            Formatted multi-line report string.
    """

    lines = [
        f'\nTED report for {uniprot_accession}',
        '=' * 68,
    ]

    domains = domain_data.get('domains', [])

    for i, domain in enumerate(domains, start=1):
        segments = domain.get('segments', [])

        lines.append(f'Domain #{i}')
        lines.append(f'TED ID            : {domain.get("ted_id")}')
        lines.append(f'Consensus level   : {domain.get("consensus_level")}')
        lines.append(f'CATH annotation   : {domain.get("cath_label")}')
        lines.append(f'Residues          : {domain.get("nres_domain")}')
        lines.append(f'Mean pLDDT        : {domain.get("plddt"):.2f}')
        lines.append(f'Segments          : {len(segments)}')

        for j, (start, end) in enumerate(segments, start=1):
            lines.append(
                f'    [{j}] residues {start}-{end} (length {end - start + 1})'
            )

        lines.append('')

    return '\n'.join(lines)


def get_alphafold_plddt(uniprot_accession: str) -> list:
    """Retrieve per-residue pLDDT scores from AlphaFold."""

    info = get_protein_info(uniprot_accession)
    if info is None:
        return None

    pdb_url = info.get('pdbUrl')
    if pdb_url is None:
        raise RuntimeError(f'Error: No PDB URL found for {uniprot_accession}')

    pdb_content = _download_from_url(pdb_url)
    
    # Extract pLDDT from PDB
    seen_residues = set()
    plddt_scores = []
    for line in pdb_content.decode().splitlines():
        if line.startswith("ATOM"):
            res_id = int(line[22:26].strip())

            # Only take first atom per residue
            if res_id not in seen_residues:
                seen_residues.add(res_id)
                plddt = float(line[60:66].strip())
                plddt_scores.append(plddt)

    return plddt_scores


def extract_plddt_regions(plddt_scores: list, plddt_range: tuple = (70,101)) -> list:
    """Convert pLDDT scores into contiguous regions based on a given pLDDT range."""

    regions = []
    start = None

    for i, score in enumerate(plddt_scores, start=1):
        if score >= plddt_range[0] and score < plddt_range[1]:
            if start is None:
                start = i
        else:
            if start is not None:
                regions.append((start, i - 1))
                start = None

    if start is not None:
        regions.append((start, len(plddt_scores)))

    return regions


def get_alphafold_data(
    uniprot_accession: str,
    detail: str = 'low+very_low',
    ) -> dict:
    """Convert AlphaFold pLDDT into disorder-like region format."""

    plddt_scores = get_alphafold_plddt(uniprot_accession)
    if not plddt_scores:
        raise RuntimeError(f'Error: Could not retrieve pLDDT scores for {uniprot_accession}')

    # Assign confidence levels
    if detail == 'all':
        confidence_levels = [
            ('very_high', '(pLDDT >= 90)', (90,101)),
            ('high', '(70 <= pLDDT < 90)', (70,90)),
            ('low', '(50 <= pLDDT < 70)', (50,70)),
            ('very_low', '(pLDDT < 50)', (0,50)),       
        ]
    elif detail == 'low+very_low':
        confidence_levels = [
            ('low+very_low', '(pLDDT < 70)', (0,71)),      
        ]

    confidence_regions_data = {}
    for level_id, level_metadata, level_range in confidence_levels:
        regions = extract_plddt_regions(plddt_scores, plddt_range=level_range)
        total_res = len(plddt_scores)
        confident_res = sum(end - start + 1 for start, end in regions)
        confidence_regions_data[f'alphafold_{level_id}_confidence {level_metadata}'] = {
            'regions': regions,
            'content_count': confident_res,
            'content_fraction': (confident_res / total_res if total_res else 0.0)
        }
    return confidence_regions_data


def format_alphafold_report(uniprot_accession: str, afdb_data: dict) -> str:
    """Format a human-readable report of data from AlphaFold Protein Structure Database.

    Args:
        uniprot_accession (str):
            UniProt accession used to extract the AFDB data.
        afdb_data (dict):
            Extracted AFDB data.

    Returns:
        str:
            Formatted multi-line report string.
    """
    lines = [
        f'\nAFDB confidence regions report for {uniprot_accession}',
        '=' * 68,
    ]

    for confidence_level, region_data in afdb_data.items():
        regions = region_data.get('regions', [])
        fraction = region_data.get('content_fraction', 0.0)
        count = region_data.get('content_count', 0)

        lines.append(f'Confidence level    : {confidence_level}')
        lines.append(f'Disordered residues : {count} ({fraction:.1%} of sequence)')
        lines.append(f'Regions             : {len(regions)}')

        for i, (start, end) in enumerate(regions, start=1):
            lines.append(f'    [{i}] residues {start}-{end}  (length {end - start + 1})')

        lines.append('')

    return '\n'.join(lines)


def suggest_targets(user_input: str):
    """Query external databases for information about a given UniProt accession and report
    on potential targets for ensemble generation.
    
    Args:
        user_input (str):
            User input to be checked for UniProt accessions and queried against external databases.
    
    Raises:
        AssertionError if no valid UniProt accession is found in the user input.
    """

    # Attempt to extract UniProt accession
    uniprotid_pattern = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|'
                                    r'[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
    uniprotid_match_targets = re.search(uniprotid_pattern,user_input)

    if uniprotid_match_targets is not None:

        # Extract UniProt accession
        accession = uniprotid_match_targets.group(0)

        # MobiDB disordered regions
        try:
            # Get disorder data from UniProt accession
            disorder_data = get_mobidb_data(accession)
            if disorder_data is not None:
                # Report on disordered regions
                report = format_mobidb_report(accession, disorder_data)
                print(report)
            else:
                raise RuntimeError(
                    f'Error: Could not retrieve MobiDB disorder information for {accession}.'
                )
        except RuntimeError as e:
            print(e)

        # TED structured domains
        try:
            # Get TED domains data from UniProt accession
            domain_data = get_ted_data(accession)
            if domain_data is not None:
                # Report on structured domains
                report = format_ted_report(accession, domain_data)
                print(report)
            else:
                raise RuntimeError(
                    f'Error: Could not retrieve TED domain annotations for {accession}.'
                )
        except RuntimeError as e:
            print(e)

        # AlphaFold Protein Structure Database confidence levels
        try:
            # Get AlphaFold data from UniProt accession
            afdb_data = get_alphafold_data(accession)
            if afdb_data is not None:
                # Report on AF results
                report = format_alphafold_report(accession, afdb_data)
                print(report)
            else:
                raise RuntimeError(
                    f'Error: Could not retrieve AlphaFold prediction data for {accession}.'
                )
        except RuntimeError as e:
            print(e)

    else:
        raise AssertionError(f'Error: {user_input} must be a valid UniProt accession!')
