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
import sklearn
from tqdm import tqdm

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

# FUNCTIONS
def derive_traj_id_from_file(filename: str) -> str:
    """Attempt to extract a trajectory identifier from a experimental data filename.

    Args:
        filename:
            name of file to extract trajectory identifier from.

    Returns:
        trajectory_id:
            the text before the first underscore of filename, or the whole filename if
            no underscores are present.
    """
    trajectory_id = filename
    filename_components = filename.split('_')
    no_exp = True
    for i,s in enumerate(filename_components):
        if s == 'exp':
            no_exp = False
            trajectory_id = '_'.join([x for x in filename_components[:i] ])
    if no_exp:
        filename_components[-1] = filename_components[-1][:-4]
        trajectory_id = '_'.join([x for x in filename_components])
    return trajectory_id


def process_exp_data(experimental_data_path: str) -> str:
    """Check formatting and units in input experimental data file.

    If values for q are in Ångstrom, they are converted to nanometer.
    Any q-values above 5nm^(-1) are removed, as SAXS calculations are not reliable in that
    range.

    Args:
        experimental_data_path:
            path to experimental SAXS data file.
    Returns:
        processed_exp_saxs:
            path to experimental SAXS data file with any applied changes.
    
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
    trajectory_id = derive_traj_id_from_file(filename=os.path.split(experimental_data_path)[1])
    processed_exp_saxs = os.path.join(os.path.split(experimental_data_path)[0],
                                      f'{trajectory_id}_exp_saxs_input_processed.dat')
    np.savetxt(processed_exp_saxs,exp_saxs_input)

    return processed_exp_saxs


def correct_exp_error(experimental_data_path: str) -> str:
    """Correct experimental error of input experimental data file using BIFT.

    Bayesian Indirect Fourier Transformation (BIFT) can identify whether the
    experimental error in small-angle scattering data is over- or
    underestimated. The error values are then scaled accordingly.

    Reference:
        Larsen, A.H. and Pedersen, M.C. (2021), Experimental noise in small-angle scattering
        can be assessed using the Bayesian indirect Fourier transformation. J. Appl. Cryst.,
        54: 1281-1289. https://doi.org/10.1107/S1600576721006877

    Args:
        experimental_data_path:
            path to experimental SAXS data file.
    Returns:
        corrected_exp_saxs:
            path to experimental SAXS data file with corrected errors.

    Adapted from:
        https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
    """
    # Deduce working directory
    working_dir, experimental_data_file = os.path.split(experimental_data_path)

    # Prepare input file for BIFT
    input_file = os.path.join(working_dir,'inputfile.dat')
    with open(input_file,'w',encoding='utf-8') as f:
        f.write(f'{experimental_data_file}\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')

    # Run BIFT
    try:
        subprocess.run(f'{GLOBAL_CONFIG["BIFT_PATH"]} < {os.path.split(input_file)[1]}',
                       shell=True,
                       cwd=working_dir,
                       stdout=open(os.path.join(working_dir,'bift.log'),'w',encoding='utf-8'),
                       stderr=subprocess.STDOUT,
                       check=True)
    except subprocess.CalledProcessError as e:
        raise e

    # Save rescaled experimental data
    trajectory_id = derive_traj_id_from_file(filename=experimental_data_file)
    corrected_exp_saxs = os.path.join(working_dir,
                                      f'{trajectory_id}_exp_saxs.dat')
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

    print( 'Experimental errors on SAXS intensities have been corrected with BIFT using scale '
          f'factor {scale_factor}.')

    return corrected_exp_saxs


def ibme(
    theta: int,
    exp_file: str,
    calc_file: str,
    output_dir: str,
    ) -> tuple[int,tuple[float,float,float],np.ndarray]:
    """Apply the Iterative Bayesian Maximum Entropy (BME) algorithm on calculated SAXS data,
    given a value for the theta parameter.

    The used algorithm is explained in:
        https://github.com/KULL-Centre/BME/blob/main/notebook/example_04.ipynb

    Reference:
        Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
        Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
        doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

    Args:
        theta:
            value for the theta parameter to be used in BME algorithm.
        exp_file:
            path to .dat file with experimental SAXS curve.
        calc_file:
            path to .dat file with SAXS curve calculated from an ensemble.
        output_dir:
            path to directory where all the files resulting from the reweighting procedure will be
            stored.

    Returns:
        A tuple (theta, stats, weights) where:
            theta:
                value for the theta parameter used in BME algorithm (same as input).
            stats:
                a tuple (chi2_before,chi2_after,phi) where:
                    chi2_before:
                        the value for the chisquare of fitting the ensemble with uniform
                        weights to the experimental data.
                    chi2_after:
                        the value for the chisquare of fitting the reweighted ensemble to
                        the experimental data.
                    phi:
                        the fraction of effective frames being used in the reweighted ensemble.
            weights:
                an array containing the new weights of the ensemble, one for each frame.

    Adapted from:
        https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
    """
    # Change current working directory
    old_cd = os.getcwd()
    os.chdir(output_dir)

    # Create reweight object
    rew = simple_BME.SimpleReweight(f'ibme_t{theta}')

    # Load files
    rew.load(exp_file=exp_file,
             calc_file=calc_file)

    # Do reweighting
    with contextlib.redirect_stdout(open(os.devnull, 'w',encoding='utf-8')):
        rew.ibme(theta=theta,
                 iterations=25,
                 ftol=0.001)

    # Restore working directory
    os.chdir(old_cd)

    weights = rew.get_ibme_weights()[-1] # get the final weights
    stats = rew.get_ibme_stats()[-1] # get the final stats

    return theta,stats,weights


def bme_ensemble_reweighting(
    exp_saxs_file: str,
    calc_saxs_file: str,
    thetas: list[int],
    output_dir: str,
    ) -> tuple[np.ndarray,np.ndarray]:
    """Perform Bayesian Maximum Entropy (BME) reweighting on calculated SAXS data based on
    experimental SAXS data of the same protein.

    Applies the iterative BME algorithm, explained in:
        https://github.com/KULL-Centre/BME/blob/main/notebook/example_04.ipynb
    The algorithm is applied using different theta values and the results for each value are stored.

    Reference:
        Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
        Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
        doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

    Args:
        exp_saxs_file:
            path to .dat file with experimental SAXS data.
        calc_saxs_file:
            path to .dat file with SAXS data calculated from a conformational ensemble.
        thetas:
            values of theta to try when applying iBME.
        output_dir:
            path to directory where output files from reweighting protocol will be stored.

    Returns:
        A tuple (stats,weights) where:
            stats:
                an array where each row corresponds to a different theta value with columns
                (chi2_before,chi2_after,phi) where:
                    chi2_before:
                        the value for the chisquare of fitting the ensemble with uniform
                        weights to the experimental data.
                    chi2_after:
                        the value for the chisquare of fitting the reweighted ensemble to
                        the experimental data.
                    phi:
                        the fraction of effective frames being used in the reweighted ensemble.
            weights:
                an array where each row corresponds to a different theta value with column
                contains the set of weights of
                the ensemble, one for each frame.
    """
    # Setup variables for parallel computing
    exp_file = os.path.abspath(exp_saxs_file)
    calc_file = os.path.abspath(calc_saxs_file)

    # Parallelize the computation across the theta values
    results = []
    with ProcessPoolExecutor() as ppe:
        futures = [ppe.submit(ibme, theta, exp_file, calc_file, output_dir) for theta in thetas]
        for future in tqdm(as_completed(futures),total=len(thetas),desc='Reweighting ensemble... '):
            results.append(future.result())

    # Extract the results
    results.sort()  # sort by theta (first element)
    thetas, stats, weights = zip(*results)

    # Convert to numpy arrays for easier indexing
    stats = np.array(stats)

    return stats,weights


def average_saxs_profiles(
    exp_saxs_file: str,
    calc_saxs_file: str,
    rw_calc_saxs_file: str,
    weights: np.ndarray,
    ) -> tuple[float,float]:
    """Average the SAXS intensities for uniform and reweighted calculated SAXS data.
    The uniform data is then scaled and offset by linear regression fitting to experimental data.

    Args:
        exp_saxs_file:
            path to .dat file with experimental SAXS data.
        calc_saxs_file:
            path to .dat file with SAXS data calculated from a conformational ensemble.
        rw_calc_saxs_file:
            path to .dat file with SAXS data calculated from a conformational ensemble considering
            the weights (from iBME) for each frame.
        weights:
            array resulting from iBME with weights for each data point. Defaults to uniform weights.

    Returns:
        A tuple i_prior,i_post where:
            i_prior:
                an array of SAXS intensities averaged over all the frames of a SAXS data file
                calculated from a conformational ensemble with uniform weights.
            i_post:
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
    model = sklearn.linear_model.LinearRegression()
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
        data:
            a DataFrame with the desired data, the path to the data in .csv format or None.
        data_msg_tag:
            string identifier for which data we are working with so prints to console are
            correct.
        calc_fn:
            an object with a __call__ method, e.g. a function to be used in calculating the
            data if it is not provided.

    Returns:
        pd.DataFrame:
            desired data in DataFrame format.
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
    ) -> tuple[str|None, str|None, np.ndarray|None,
               np.ndarray|None, np.ndarray|None]:
    """Attempt to read reweighting data from output directory, returning None if not found.

    Args:
        reweighting_output_directory:
            directory where data should be searched.
        trajectory_id:
            prefix for filenames to look for in directory.
    Returns:
        A tuple exp_saxs_file, calc_saxs_file, thetas_array, stats, weights where each variable
        is either the corresponding data (if found) or None (if not found).
    """
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
    ## weights = 'ibme_t{THETA_VALUE}.weights.dat'
    ## stats = 'ibme_t{THETA_VALUE}_ibme_{ITERATION_NUMBER}.log' with the highest ITERATION_NUMBER

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
