"""Auxiliary functions for reweighting ensembles."""

# IMPORTS
## Standard Library Imports
import contextlib
import glob
import os
import subprocess
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed

## Third Party Imports
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.linear_model import LinearRegression

## Local Imports
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.reweighting.third_party import simple_BME

# CONSTANTS
ERROR_READ_MSG = 'Error in reading {} from file.'
SUCCESSFUL_READ_MSG = '{} has been read from file.'
PROVIDED_MSG = '{} data has been provided.'
NOT_PROVIDED_MSG = 'No {} data was provided.'
ALLOWED_DATA_MSG_TAGS = ('cmatrix','dmatrix','ss_freq','structural_metrics')
DATA_NAMES = {'cmatrix': 'contact matrix',
              'dmatrix': 'distance matrix',
              'ss_freq': 'secondary structure assignment frequency matrix',
              'structural_metrics': 'structural metrics distributions'}
FOUND_EXP_SAXS_MSG = 'Found processed experimental SAXS data.'
FOUND_EXP_CALC_SAXS_MSG = 'Found processed experimental SAXS data and calculated SAXS data.'
FOUND_RW_DATA_MSG = 'Found calculated BME reweighting data with theta values: {}'
EXP2FIT = {'SAXS': 'scale+offset',
           'RDC': 'scale',
           'CS': None,
           'JCOUPLINGS': None}
FIT_TYPES = ['scale','scale+offset']

# FUNCTIONS
def _derive_traj_id_from_file(filename: str) -> str:
    """Attempt to extract a trajectory identifier from a experimental data filename.

    Args:
        filename (str):
            Name of file to extract trajectory identifier from.

    Returns:
        str:
            The text before the first experimental type identifier, or the whole filename if
            no experimental types are present.
    """
    # Split the filename into components and handle the file extension
    filename_components, _ = os.path.splitext(filename)
    filename_components = filename_components.split('_')

    # Find the earliest occurrence of any experimental type
    earliest_exp_index = len(filename_components)  # Default to end of list
    
    for exp_type in EXP2FIT.keys():
        try:
            exp_index = filename_components.index(exp_type)
            earliest_exp_index = min(earliest_exp_index, exp_index)
        except ValueError:
            # exp_type not found in filename_components
            continue
    
    # Use components before the earliest experimental type
    trajectory_id = '_'.join(filename_components[:earliest_exp_index])
    
    return trajectory_id


def process_exp_saxs_data(experimental_data_path: str) -> str:
    """Check formatting and units in input experimental SAXS data file.

    If values for q are in Ångstrom, they are converted to nanometer.
    Any q-values above 5nm^(-1) are removed, as SAXS calculations are not reliable in that
    range.

    Args:
        experimental_data_path (str):
            Path to experimental SAXS data file.
    Returns:
        str:
            Path to experimental SAXS data file with any applied changes.
    
    Adapted from:
        https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
    """
    # Read experimental data into array
    exp_saxs_input = np.loadtxt(experimental_data_path)
    assert np.shape(exp_saxs_input)[1] == 3, ('Unable to read file. Make sure the file only '
                                              'contains 3 columns (q,I,sigma) '
                                              'and #commented lines.')

    # Convert from A to nm if required
    if exp_saxs_input[...,0][-1] <  1:
        print('q is in Å units. Converting to nm.')
        exp_saxs_input[...,0] = exp_saxs_input[...,0]*10

    # Remove values from the SAXS profile related to noise or aggregation (> 5 nm^(-1))
    n_unreliable = (exp_saxs_input[...,0] >= 5).sum()
    if n_unreliable > 0:
        print(f'Found {n_unreliable} q-values above 5 nm^(-1). SAXS calculations are not reliable '
               'in that region of the spectrum. Those datapoints will be removed.')
        exp_saxs_input = exp_saxs_input[(exp_saxs_input[...,0] < 5)]

    # Save processed data
    trajectory_id = _derive_traj_id_from_file(filename=os.path.basename(experimental_data_path))
    processed_exp_saxs = os.path.join(os.path.dirname(experimental_data_path),
                                      f'{trajectory_id}_SAXS_exp_input_processed.dat')
    np.savetxt(processed_exp_saxs,
               exp_saxs_input,
               header=' DATA=SAXS')
    
    return processed_exp_saxs


def correct_exp_saxs_error(experimental_data_path: str) -> str:
    """Correct experimental error of input experimental data file using BIFT.

    Bayesian Indirect Fourier Transformation (BIFT) can identify whether the
    experimental error in small-angle scattering data is over- or
    underestimated. The error values are then scaled accordingly.

    Reference:
        Larsen, A.H. and Pedersen, M.C. (2021), Experimental noise in small-angle scattering
        can be assessed using the Bayesian indirect Fourier transformation. J. Appl. Cryst.,
        54: 1281-1289. https://doi.org/10.1107/S1600576721006877

    Args:
        experimental_data_path (str):
            Path to experimental SAXS data file.
    Returns:
        str:
            Path to experimental SAXS data file with corrected errors.

    Adapted from:
        https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
    """
    # Assign working directory
    working_dir, experimental_data_file = os.path.split(experimental_data_path)

    # Prepare input file for BIFT
    input_file = os.path.join(working_dir,'inputfile.dat')
    with open(input_file,'w',encoding='utf-8') as f:
        f.write(f'{experimental_data_file}' + ('\n' * 18))

    # Run BIFT
    try:
        subprocess.run(f'{GLOBAL_CONFIG["BIFT_PATH"]} < {os.path.basename(input_file)}',
                       shell=True,
                       cwd=working_dir,
                       stdout=open(os.path.join(working_dir,'bift.log'),'w',encoding='utf-8'),
                       stderr=subprocess.STDOUT,
                       check=True)
    except subprocess.CalledProcessError as e:
        raise e

    # Save experimental data with corrected errors
    trajectory_id = _derive_traj_id_from_file(filename=experimental_data_file)
    corrected_exp_saxs = os.path.join(working_dir,
                                      f'{trajectory_id}_SAXS_exp.dat')
    np.savetxt(corrected_exp_saxs,
               np.loadtxt(os.path.join(working_dir,'rescale.dat')),
               header=' DATA=SAXS')

    # Save scale factor
    scale_factor = np.loadtxt(os.path.join(working_dir,'scale_factor.dat'))[0,1]

    # Save used parameters in log file
    with (open(os.path.join(working_dir,'parameters.dat'),'r',encoding='utf-8-sig') as f,
          open(os.path.join(working_dir,'bift.log'),'a',encoding='utf-8') as t):
        t.write('------------------------ Input Parameters ------------------------\n')
        t.write(f.read())

    # Remove unnecessary bift output files
    for tmp_file in ['data.dat','dummy.dat','fit.dat','gs_pr.dat','gx_pr.dat','inputfile.dat',
                     'parameters.dat','pr.dat','rescale.dat','scale_factor.dat','st_pr.dat']:
        os.remove(os.path.join(working_dir,tmp_file))

    print('Experimental errors on SAXS intensities have been corrected with BIFT using scale '
          f'factor {scale_factor}.')

    return corrected_exp_saxs


def bme(
    theta: int,
    exp_files: str | list[str],
    calc_files: str | list[str],
    output_dir: str,
    exp_types: str | None | list[str|None],
    ) -> tuple[int,tuple[float,float,float],np.ndarray]:
    """Apply the Bayesian Maximum Entropy (BME) algorithm.
     
    Uses the provided value for the theta parameter, and possibly applies iterative BME
    according to the type of experimental data used.

    Reference:
        Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
        Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
        doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

    Args:
        theta (int):
            Value for the theta parameter to be used in BME algorithm.
        exp_files (str | list[str]):
            Path to .dat file with experimental SAXS curve.
        calc_file (str | list[str]):
            Path to .dat file with SAXS curves calculated for each conformer of an ensemble.
        output_dir (str):
            Path to directory where all the files resulting from the reweighting procedure will be
            stored.
        exp_types (str | None | list[str|None]):
            Type(s) of experimental data provided.
            If a list is provided, it must follow the same order as the exp_files list.

    Returns:
        tuple[int, tuple[float, float, float], np.ndarray]:
            int:
                Value for the theta parameter used in BME algorithm (same as input).
            tuple[float, float, float]:
                float:
                    The value for the chisquare of fitting the ensemble with uniform
                    weights to the experimental data.
                float:
                    The value for the chisquare of fitting the reweighted ensemble to
                    the experimental data.
                float:
                    The fraction of effective frames being used in the reweighted ensemble.
            np.ndarray:
                An array containing the new weights of the ensemble, one for each frame.

    Adapted from:
        https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
    """
    # Change current working directory
    old_cd = os.getcwd()
    os.chdir(output_dir)

    # Setup paths for experimental and calculated data files
    if isinstance(exp_files, str):
        exps = [exp_files]
    elif isinstance(exp_files, list):
        exps = [x for x in exp_files]

    if isinstance(calc_files, str):
        calcs = [calc_files]
    elif isinstance(calc_files, list):
        calcs = [x for x in calc_files]

    if isinstance(exp_types, str):
        types = [exp_types]
    elif isinstance(exp_types, list):
        types = [x for x in exp_types]

    # Create reweight object
    rew = simple_BME.SimpleReweight(f'ibme_t{theta}')

    # Load experimental and calculated data files, scale and offset if necessary
    for exp, calc, exp_type in zip(exps,calcs,types):
        rew.load(exp_file=exp,
                 calc_file=calc,
                 exp_type=exp_type)

    # Do reweighting
    with contextlib.redirect_stdout(open(os.devnull, 'w',encoding='utf-8')):
        rew.ibme(theta=theta,
                 iterations=25,
                 ftol=0.001)

    # Restore working directory
    os.chdir(old_cd)

    weights = rew.get_ibme_weights()[-1] # get the final weights
    stats = rew.get_ibme_stats()[-1] # get the final stats

    return theta, stats, weights


def bme_ensemble_reweighting(
    exp_data: str | list[str],
    exp_type: str | list[str],
    calc_data: str | list[str],
    thetas: list[int],
    output_dir: str,
    ) -> tuple[np.ndarray,np.ndarray]:
    """Apply Bayesian/Maximum Entropy (BME) reweighting on calculated+experimental data.

    The algorithm is applied using different theta values and the results for each value are stored.

    Reference:
        Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
        Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
        doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

    Args:
        exp_data (str | list[str]):
            Path to .dat file(s) with experimental data.
        exp_type (str | list[str]):
            Type(s) of experimental data. If a list is provided, it must follow the same order as
            the exp_data list.
        calc_saxs_file (str):
            Path to .dat file(s) with experimental data calculated from a conformational ensemble.
        thetas (list[int]):
            Values of theta to try when applying BME.
        output_dir (str):
            Path to directory where output files from reweighting protocol will be stored.

    Returns:
        tuple[np.ndarray,np.ndarray]:
            stats (np.ndarray):
                An array where each row corresponds to a different theta value with columns
                (chi2_before,chi2_after,phi) where:
                    chi2_before:
                        The value for the chisquare of fitting the ensemble with uniform
                        weights to the experimental data.
                    chi2_after:
                        The value for the chisquare of fitting the reweighted ensemble to
                        the experimental data.
                    phi:
                        The fraction of effective frames being used in the reweighted ensemble.
            weights (np.ndarray):
                An array where each row corresponds to a different theta value with columns
                containing the set of weights of the ensemble, one for each frame.
    """
    # Check inputs
    if isinstance(exp_data, str):
        exp_data_input = [exp_data]
    elif isinstance(exp_data, list):
        exp_data_input = [x for x in exp_data]

    if isinstance(exp_type, str):
        exp_type_input = [exp_type]
    elif isinstance(exp_type, list):
        exp_type_input = [x for x in exp_type]

    if isinstance(calc_data, str):
        calc_data_input = [calc_data]
    elif isinstance(calc_data, list):
        calc_data_input = [x for x in calc_data]

    # Setup necessary BME inputs
    exp_files_BME = [os.path.abspath(x) for x in exp_data_input]
    exp_types_BME = [ x for x in exp_type_input]
    calc_files_BME = [os.path.abspath(x) for x in calc_data_input]

    # Parallelize the computation across the theta values
    results = []
    with ProcessPoolExecutor() as ppe:
        futures = [ppe.submit(bme,
                              theta,
                              exp_files_BME,
                              calc_files_BME,
                              output_dir,
                              exp_types_BME) for theta in thetas]
    
        for future in tqdm(as_completed(futures),
                           desc='Reweighting ensemble... ',
                           total=len(thetas)):
            
            results.append(future.result())

    # Extract the results
    results.sort()  # sort by theta (first element in tuples)
    thetas, stats, weights = zip(*results)

    # Convert to numpy arrays for easier indexing
    stats = np.array(stats)

    return stats, weights


def average_saxs_profiles(
    exp_saxs_file: str,
    calc_saxs_file: str,
    rw_calc_saxs_file: str,
    weights: np.ndarray,
    ) -> tuple[float,float]:
    """Average the SAXS intensities for uniform and reweighted calculated SAXS data.
    The uniform data is then scaled and offset by linear regression fitting to experimental data.

    Args:
        exp_saxs_file (str):
            Path to .dat file with experimental SAXS data.
        calc_saxs_file (str):
            Path to .dat file with SAXS data calculated from a conformational ensemble.
        rw_calc_saxs_file (str):
            Path to .dat file with SAXS data calculated from a conformational ensemble considering
            the weights (from iBME) for each frame.
        weights (np.ndarray):
            Array resulting from iBME with weights for each data point. Defaults to uniform weights.

    Returns:
        tuple[float,float]:
            i_prior (float):
                an array of SAXS intensities averaged over all the frames of a SAXS data file
                calculated from a conformational ensemble with uniform weights.
            i_post (float):
                an array of SAXS intensities averaged over all the frames of a SAXS data file
                calculated from a conformational ensemble with the provided set of weights.
    """
    # Average the uniform and reweighted calculated saxs intensities
    i_prior = np.average(np.loadtxt(calc_saxs_file)[...,1:],
                         axis=0)
    i_post = np.average(np.loadtxt(rw_calc_saxs_file)[...,1:],
                        axis=0,
                        weights=weights)

    # Get experimental intensities and errors
    _, i_exp, err = np.loadtxt(exp_saxs_file,
                               unpack=True)

    # Perform ordinary least squares Linear Regression, fitting the calculated SAXS data to the
    # experimental SAXS data, taking into account the experimental error of each data point
    model = LinearRegression()
    model.fit(X=i_prior.reshape(-1,1),
              y=i_exp,
              sample_weight=1/(err**2))

    # Adjust uniform saxs profile by linear fitting of scale and offset from experimental data
    a = model.coef_[0]
    b = model.intercept_
    i_prior = a * i_prior + b

    return i_prior,i_post


def attempt_read_calculated_data(
    data: pd.DataFrame | str | None,
    data_msg_tag: str,
    calc_fn: Callable,
    *args,
    **kwargs,
    ) -> pd.DataFrame:
    """Attempt to read data from file, else calculate it using provided function.

    If data is given directly as a DataFrame, it is simply returned. Otherwise, it
    is either read from file or calculated using the provided function and arguments.

    Args:
        data (pd.DataFrame | str | None):
            A DataFrame with the desired data, the path to the data in .csv format or None.
        data_msg_tag (str):
            String identifier for which data we are working with so prints to console are
            correct.
        calc_fn (Callable):
            An object with a __call__ method, e.g. a function to be used in calculating the
            data if it is not provided.

    Returns:
        pd.DataFrame:
            Desired data in DataFrame format.
    """
    assert data_msg_tag in ALLOWED_DATA_MSG_TAGS, ('Data message tag must be in '
                                                   f'{ALLOWED_DATA_MSG_TAGS} !')
    data_name = DATA_NAMES[data_msg_tag]

    # Attempt read of data
    if isinstance(data,str):
        assert data.endswith('.csv'), ('Calculated data must be in provided .csv format!')
        try:
            data_df = pd.read_csv(data,index_col=0)
        except:
            print(ERROR_READ_MSG.format(data_name))
            data_df = calc_fn(*args,**kwargs)
        else:
            print(SUCCESSFUL_READ_MSG.format(data_name.capitalize()))
    elif isinstance(data,pd.DataFrame):
        print(PROVIDED_MSG,format(data_name.capitalize()))
        data_df = data
    else:
        print(NOT_PROVIDED_MSG.format(data_name))
        data_df = calc_fn(*args,**kwargs)
    return data_df


def attempt_read_reweighting_data(
    reweighting_output_directory: str,
    trajectory_id: str,
    exp_type: str,
    ) -> tuple[str | None, str | None, np.ndarray | None,
               np.ndarray | None, np.ndarray | None]:
    """Attempt to read reweighting data from output directory, returning None if not found.

    Args:
        reweighting_output_directory (str):
            Directory where data should be searched.
        trajectory_id (str):
            Prefix for filenames to look for in directory.
    Returns:
        tuple[str | None, str | None, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
            str | None:
                The experimental data (if found) or None (if not found).
            str | None:
                The calculated data (if found) or None (if not found).
            np.ndarray | None:
                The array of theta values (if found) or None (if not found).
            np.ndarray | None:
                The fitting statistics (if found) or None (if not found).
            np.ndarray | None:
                The set of weights (if found) or None (if not found).
    """
    # TODO Add logic for other exp data
    # Check for experimental SAXS data file
    exp_saxs_file = os.path.join(reweighting_output_directory,f'{trajectory_id}_exp_saxs.dat')
    if not os.path.isfile(exp_saxs_file):
        exp_saxs_file = None
        return None, None, None, None, None

    # Check for calculated SAXS data file
    calc_saxs_file = os.path.join(reweighting_output_directory,f'{trajectory_id}_calc_saxs.dat')
    if not os.path.isfile(calc_saxs_file):
        calc_saxs_file = None
        print(FOUND_EXP_SAXS_MSG)
        return exp_saxs_file, None, None, None, None

    # Check for BME reweighting results
    bme_results_dir = os.path.join(reweighting_output_directory,'bme_reweighting_results')
    if not os.path.isdir(bme_results_dir):
        print(FOUND_EXP_CALC_SAXS_MSG)
        return exp_saxs_file, calc_saxs_file, None, None, None

    ## Check which theta values are present (if any)
    theta_values = []
    for theta_log in glob.glob(os.path.join(bme_results_dir,'ibme_t*.log')):
        theta_log_prefix = os.path.split(theta_log)[1].split('_')[1]
        if '.log' not in theta_log_prefix:
            theta_value = int(theta_log_prefix[1:])
            if theta_value not in theta_values:
                theta_values.append(theta_value)

    if not theta_values:
        print(FOUND_EXP_CALC_SAXS_MSG)
        return exp_saxs_file, calc_saxs_file, None, None, None

    ## Check which weights/stats are present (if any)
    ## weights : f'ibme_t{THETA_VALUE}.weights.dat'
    ## stats : f'ibme_t{THETA_VALUE}_ibme_{ITERATION_NUMBER}.log' with the highest ITERATION_NUMBER

    all_weights = []
    all_stats = []

    for theta in sorted(theta_values):
        # Get weights
        weights_files = glob.glob(os.path.join(bme_results_dir,f'ibme_t{theta}_*.weights.dat'))
        try:
            weights = np.loadtxt(weights_files[0],
                                 usecols=1)
            all_weights.append(weights)
        except (IndexError, FileNotFoundError):
            print(FOUND_EXP_CALC_SAXS_MSG)
            return exp_saxs_file, calc_saxs_file, None, None, None

        # Get stats
        stats_files_sorted = sorted(glob.glob(os.path.join(bme_results_dir,
                                                           f'ibme_t{theta}_ibme_*.log')),
                                    key=lambda x : int(os.path.split(x)[1].split('_')[-1][:-4]))
        try:
            with open(stats_files_sorted[-1],'r',encoding='utf-8-sig') as stats_file:
                lines = stats_file.readlines()
            chi2_before = float(lines[2].strip().split(' ')[-1])
            chi2_after = float(lines[5].strip().split(' ')[-1])
            phi = float(lines[-1].strip().split(' ')[-1])
            stats = (chi2_before,chi2_after,phi)
            all_stats.append(stats)
        except (IndexError, FileNotFoundError):
            print(FOUND_EXP_CALC_SAXS_MSG)
            return exp_saxs_file, calc_saxs_file, None, None, None

    if len(all_weights) != len(theta_values) or len(all_stats) != len(theta_values):
        print(FOUND_EXP_CALC_SAXS_MSG)
        return exp_saxs_file, calc_saxs_file, None, None, None

    weights = np.array(all_weights)
    stats = np.array(all_stats)
    thetas_array = np.array(sorted(theta_values))

    print(FOUND_RW_DATA_MSG.format(sorted(theta_values)))

    return exp_saxs_file, calc_saxs_file, thetas_array, stats, weights
