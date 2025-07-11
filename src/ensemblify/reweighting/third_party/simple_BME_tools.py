"""
The code below was adapted from https://github.com/KULL-Centre/BME.
See the THIRD_PARTY_NOTICE_BME.txt file for more details.

Reference:
    S. Bottaro , T. Bengsten and K. Lindorff-Larsen, "Integrating Molecular Simulation and
    Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach," pp. 219-240, Feb.
    2020. In: Z. Gáspári, (eds) *Structural Bioinformatics*, *Methods in Molecular Biology*,
    vol. 2112, Humana, New York, NY. (https://doi.org/10.1007/978-1-0716-0270-6_15)
"""

# IMPORTS
## Third Party Imports
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# CONSTANTS
EXP_TYPES = ['JCOUPLINGS','CS','SAXS','RDC']
FIT_TYPES = ['scale', 'scale+offset']

# FUNCTIONS
def parse(exp_file: str, calc_file: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Parse provided files, do sanity checks and read data.
    
    Args:
        exp_file (str):
            Path to the experimental data file.
        calc_file (str):
            Path to the calculated data file.
    
    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, str]:
            np.ndarray:
                Labels of the experimental data.
            np.ndarray:
                Experimental data in the format {value, error}.
            np.ndarray:
                Calculated data in the format {value}.
            str:
                Log message summarizing the data read and any checks performed.
    """

    # Setup log msg
    log = ''

    # Check experimental data file format
    with open(exp_file,'r',encoding='utf-8-sig') as fh:
        first = fh.readline()
    assert first[0] == '#', (f'Error. First line of exp file {exp_file} must be in the format '
                             f'# DATA={EXP_TYPES}')

    # Check experimental data type
    data_string = (first.split('DATA=')[-1].split()[0]).strip()
    assert data_string in EXP_TYPES , (f'Error. DATA in {exp_file} must be one of the following: '
                                       f'{EXP_TYPES} ')
    log += f'# Reading {data_string} data \n'

    # Read experimental data from file into DataFrame
    df_exp = pd.read_csv(exp_file,
                         sep='\s+',
                         header=None,
                         comment='#')

    # Sanity check and column renaming
    assert df_exp.shape[1] == 3, ('Error. Experimental datafile must be in the format '
                                  'LABEL VALUE ERROR')
    df_exp = df_exp.rename(columns={0: 'label',
                                    1: 'val',
                                    2: 'sigma'})

    # Read calculated data from file into DataFrame
    df_calc = pd.read_csv(calc_file,
                          sep='\s+',
                          header=None,
                          comment='#')
    # Drop frame numbering
    df_calc = df_calc.drop(columns=[0])

    # Check that the number of experimental data matches the calculated data
    assert df_calc.shape[1] == df_exp.shape[0], ('Error: Number of experimental data in '
                                                 f'{exp_file} ({df_exp.shape[0]}) must match '
                                                 f'the calculated data in {calc_file} '
                                                 f'({df_calc.shape[1]})')

    log += f'# Reading {df_exp.shape[0]} experimental data from {exp_file} \n'
    log += f'# Reading {df_calc.shape[0]} calculated samples from {calc_file} \n'

    # Extract data labels into separate array
    labels = df_exp['label'].values

    # Extract experimental data values and respective errors into separate array
    # The experimental data is in the format {value, error}
    exp = np.array(df_exp[['val','sigma']].values)

    # Extract calculated data values into separate array
    # The calculated data is in the format {value}
    calc = np.array(df_calc.values)

    return labels, exp, calc, log


def subsample(
    labels: np.ndarray,
    exp: np.ndarray,
    calc: np.ndarray,
    use_samples: list | None = None,
    use_data: list | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Reduce data samples to the specified subsets. If no subset is specified,
    the full unchanged dataset is returned.

    Args:
        labels (np.ndarray):
            Labels of the experimental data.
        exp (np.ndarray):
            Experimental data in the format {I, sigma, bound}.
        calc (np.ndarray):
            Calculated data in the format {I, sigma, bound}.
        use_samples (list | None):
            List of indices of the samples to use from the calculated data.
        use_data (list | None):
            List of indices of the experimental data to use.
    
    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, str]:
            np.ndarray:
                Labels of the experimental data.
            np.ndarray:
                Experimental data in the format {I, sigma, bound}.
            np.ndarray:
                Calculated data in the format {I, sigma, bound}.
            str:
                Log message summarizing the subset of data used, if applicable.
    """
    # Setup log msg
    log = ''

    # Check if any subset is specified
    if use_samples is not None:
        calc = calc[use_samples,:]
        log += f'# Using a subset of samples ({calc.shape[0]}) \n'
    if use_data is not None:
        labels = labels[use_data]
        exp = exp[use_data,:]
        calc = calc[:,use_data]
        log += f'# Using a subset of datapoints ({exp.shape[0]}) \n'

    return labels,exp,calc,log


def check_data(
    labels: np.ndarray,
    exp: np.ndarray,
    calc: np.ndarray,
    sample_weights: np.ndarray,
    ) -> str:
    """Initial sanity check between calculated and experimental data.
    
    Args:
        labels (np.ndarray):
            Labels of the experimental data.
        exp (np.ndarray):
            Experimental data in the format {value, error}.
        calc (np.ndarray):
            Calculated averages for each data point in the sample, in the format {value}.
        sample_weights (np.ndarray):
            Weights for each sample in the calculated data.
    
    Returns:
        str:
            Log message summarizing the comparison between calculated and experimental data.
    """
    # Setup log msg
    log = ''

    # Calculate the average value for each calculated data point in the sample,
    # weighted by sample weights
    calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                      axis=0)

    # Calculate the difference between calculated averages and experimental values
    diff = calc_avg - exp[:,0]
    
    # Compare calculated and experimental data using reduced chi2, RMSD, and value differences
    chi2 = np.average((diff / exp[:,1]) ** 2)
    rmsd = np.sqrt(np.average(diff**2))
    viol = np.where(diff > 1.)[0]
    log += f'CHI2: {chi2:.5f} \n'
    log += f'RMSD: {rmsd:.5f} \n'
    log += f'VIOLATIONS (Data points with difference > 1.0): {len(viol)} \n'

    # Check for values in calculated data that are outside the experimental range
    m_min = np.min(calc,axis=0)
    m_max = np.max(calc,axis=0)

    diff_min = (m_min - exp[:,0]) / exp[:,1]
    ii_min = np.where(diff_min > 1.)[0]
    if len(ii_min) > 0:
        log += '##### WARNING ########## \n'
        log += '# The minimum value of the following data is higher than expt range: \n'
        log += f'# {"label":14s} {"exp_avg":8s} {"exp_sigma":8s} {"min_calc":8s} \n'
        for j in ii_min:
            log += f'# {labels[j]:15f} {exp[j,0]:8.4f} {exp[j,1]:8.4f} {m_min[j]:8.4f} \n'

    diff_max = (exp[:,0] - m_max) / exp[:,1]
    ii_max = np.where(diff_max > 1.)[0]
    if len(ii_max) > 0:
        log += '##### WARNING ########## \n'
        log += '# The maximum value of the following data is higher than expt range: \n'
        log += f'# {"label":14s} {"exp_avg":8s} {"exp_sigma":8s} {"max_calc":8s} \n'
        for j in ii_max:
            log += f'# {labels[j]:15f} {exp[j,0]:8.4f} {exp[j,1]:8.4f} {m_max[j]:8.4f} \n'

    return log


def normalize(
    exp: np.ndarray,
    calc: np.ndarray,
    sample_weights: np.ndarray,
    normalization_type: str = 'z-score',
    ) -> tuple[str,float,float]:
    """Normalize experimental and calculated data.
    
    Modifies the data arrays in place.
    
    Args:
        exp (np.ndarray):
            Experimental data in the format {value, error}.
        calc (np.ndarray):
            Calculated averages for each data point in the sample, in the format {value}.
        sample_weights (np.ndarray):
            Weights for each sample in the calculated data.
        normalization_type (str):
            Type of standardization to apply. Currently only 'z-score' and 'min-max' is supported.

    Returns:
        tuple[str, float, float]:
            str:
                Log message indicating standardization type.
            float:
                Average value used in the standardization.
            float:
                Standard deviation value used in the standardization.

    Raises:
        ValueError:
            If an unknown standardization type is provided.
    """
    # Setup log msg
    log = ''

    if normalization_type == 'z-score':
        # Calculate the average value for each calculated data point in the sample,
        # weighted by sample weights
        calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                        axis=0)

        # Calculate the variance of calculated values around their weighted average
        calc_var = np.average((calc_avg - calc)**2,
                            weights=sample_weights,
                            axis=0)

        # Combine calculated and experimental uncertainties into a single standard deviation
        calc_std = np.average(np.array([np.sqrt(calc_var),
                                        exp[:,1]]),
                            axis=0)

        # Perform Z-score normalization
        exp[:,0] -= calc_avg # center experimental values
        exp[:,0] /= calc_std # scale by combined uncertainty
        calc -= calc_avg # center calculated values  
        calc /= calc_std # scale by combined uncertainty
        exp[:,1] /= calc_std # scale experimental errors
        log += '# Z-score normalization \n'
    
    elif normalization_type == 'min-max':
        # Get min and max values of calculated data
        mmin = calc.min(axis=0)
        mmax = calc.max(axis=0)
        delta = mmax-mmin

        # Peform Min-Max normalization
        exp[:,0] -= mmin
        exp[:,0] /= delta
        calc -= mmin
        calc /= delta
        exp[:,1] /= delta
        log += '# MinMax normalization \n'

    else:
        raise ValueError(f'Unknown standardization type: {normalization_type}')

    return log, calc_avg, calc_std
    

def split_array_by_groups(array, groups_dict, type=None):
    """
    Splits a NumPy array into multiple arrays based on named groups of indices.

    Type of splitting can be specified as either 'row' or 'column'.
    If 'row', the array is split into groups of rows.
    If 'column', the array is split into groups of columns.
    
    Args:
        array: NumPy array to split.
        groups_dict: Dictionary with group names as keys and column indices as values

    Returns:
        Dictionary with group names as keys and split arrays as values
    """
    if type == 'column':
        return {name: array[:, indices[0]:indices[1]+1] for name, indices in groups_dict.items()}
    elif type == 'row':
        return {name: array[indices[0]:indices[1]+1, :] for name, indices in groups_dict.items()}
    else:
        raise ValueError("Type must be 'row' or 'column'.")
    

def fit_and_scale(
    exp: np.ndarray,
    calc: np.ndarray,
    calc_weights: np.ndarray | None = None,
    fit_type: str = None,
    ) -> tuple[np.ndarray,str]:
    """Fit a linear regression model between calculated and experimental data and scale calculated
    data accordingly.

    Args:
        exp (np.ndarray):
            Experimental data in the format {value, error}.
        calc (np.ndarray):
            Calculated data. Will be averaged and possibly scaled.
        calc_weights (np.ndarray | None):
            Weights for each row in the calculated data, used in averaging.
            If None, uniform weights are assumed.
        fit_type (str):
            Type of data scaling to be applied after fitting. Can be 'scale', 'scale+offset'
            or None. If None, no data scaling is applied and input data is not changed.

    Returns:
        np.ndarray:
            The (possibly) rescaled calculated data (if scale_type!=None).
        str:
            Log message summarizing the fit and scaling process.
    """
    # Check input
    assert fit_type in FIT_TYPES, f'Fit type must be in {FIT_TYPES} or None'

    # Setup log msg
    log = ''

    log += f'# Using scaling: {fit_type} \n'

    # Check if scaling is required.
    if fit_type is None:
        return calc, log

    # Setup calculated data weights
    if calc_weights is None:
        # If no weights are provided, assume uniform weights
        calc_weights = np.ones(calc.shape[0])/calc.shape[0]

    # Calculate the average value for each calculated data point in the sample,
    # weighted by sample weights
    calc_avg = np.sum(calc*calc_weights[:,np.newaxis],
                      axis=0)
    
    # Extract experimental values
    exp_avg = exp[:,0]

    # Perform appropriate fit
    if fit_type == 'scale':
        fit_intercept = False
    else:
        fit_intercept = True

    # Since the magnitude of experimental noise is not constant, the experimental data
    # are heteroskedastic. To improve our maximum likelihood estimation in the LinearRegression
    # process below, we use the inverse of the variance of experimental noise at each point as
    # the experimental sample weights.
    inv_var = 1. / exp[:,1]**2

    # Perform linear regression fitting towards experimental data
    model = LinearRegression(fit_intercept=fit_intercept)
    model.fit(X=calc_avg.reshape(-1,1),
              y=exp_avg.reshape(-1,1),
              sample_weight=inv_var)
    
    # Get the slope and intercept of the fitted model
    slope = model.coef_[0]
    intercept = model.intercept_
    calc = slope * calc + intercept # scale calculated data

    log += f'# Slope={slope}; Offset={intercept}\n'
    
    return calc, log


def calc_chi(exp: np.ndarray, calc: np.ndarray, sample_weights: np.ndarray) -> float:
    """Calculate the reduced chi2 value between experimental and calculated data.
    
    Args:
        exp (np.ndarray):
            Experimental data in the format {value, error}.
        calc (np.ndarray):
            Calculated averages for each data point in the sample, in the format {value}.
        sample_weights (np.ndarray):
            Weights for each sample in the calculated data.

    Returns:
        float:
            Reduced chi2 value.
    """
    # Calculate the average value for each calculated data point in the sample,
    # weighted by sample weights
    calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                      axis=0)

    # Calculate the difference between calculated averages and experimental values
    diff = calc_avg-exp[:,0]

    # Calculate the reduced chi2 value using experimental errors
    chi2 = np.average((diff/exp[:,1])**2)

    return chi2


def srel(w0: np.ndarray, w1: np.ndarray) -> float:
    """Calculate relative Shannon entropy between two sets of weight values.
    
    Args:
        w0 (np.ndarray):
            Reference weights.
        w1 (np.ndarray):
            Weights to compare against the reference.
    
    Returns:
        float:
            Relative Shannon entropy between the two sets of weights.
    """
    # Filter out very small values to avoid numerical issues
    idxs = np.where(w1 > 1.0E-50)

    return np.sum(w1[idxs] * np.log(w1[idxs] / w0[idxs]))
