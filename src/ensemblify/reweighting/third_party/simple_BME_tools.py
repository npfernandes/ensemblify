# IMPORTS
## Third Party Imports
import numpy as np
import pandas as pd

# CONSTANTS
EXP_TYPES = ['SAXS'] # ['NOE','JCOUPLINGS','CS','SAXS','RDC']
BOUND_TYPES = ['UPPER','LOWER']

# FUNCTIONS
def srel(w0: np.ndarray, w1: np.ndarray) -> float:
    """Calculate relative entropy."""
    idxs = np.where(w1 > 1.0E-50)
    return np.sum(w1[idxs] * np.log(w1[idxs] / w0[idxs]))


def parse(exp_file: str, calc_file: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Parse files."""

    # Setup log msg
    log = []

    # Check experimental data file format
    with open(exp_file,'r',encoding='utf-8-sig') as fh:
        first = fh.readline()
    assert first[0] == '#', (f'Error. First line of exp file {exp_file} must be in the format '
                              '# DATA=[{EXP_TYPES}] [BOUND=UPPER/LOWER]')

    # Check experimental data type
    data_string = (first.split('DATA=')[-1].split()[0]).strip()
    assert data_string in EXP_TYPES , (f'Error. DATA in {exp_file} must be one of the following: '
                                       f'{EXP_TYPES} ')

    log.append(f'# Reading {data_string} data \n')

    # If it is not an average but a boundary it can be specified
    bound_string = None
    if len(first.split('BOUND=')) == 2:
        bound_string = (first.split('BOUND=')[-1].split()[0]).strip()
        assert bound_string in BOUND_TYPES , (f'Error. {bound_string} is not known. BOUND in '
                                              f'{exp_file} must be one of the following: '
                                              f'{BOUND_TYPES} ')
        log.append(f'# {bound_string}-bound data \n')

    # Read experimental data
    df_exp = pd.read_csv(exp_file,
                         sep='\s+',
                         header=None,
                         comment='#')

    assert df_exp.shape[1]==3, ('Error. Experimental datafile must be in the format '
                                'LABEL VALUE ERROR')
    df_exp = df_exp.rename(columns={0: 'label', 1: 'val',2:'sigma'})

    # Read calculated datafile
    df_calc = pd.read_csv(calc_file,
                          sep='\s+',
                          header=None,
                          comment='#')
    # Drop frame
    df_calc = df_calc.drop(columns=[0])

    assert df_calc.shape[1] == df_exp.shape[0], ('Error: Number of experimental data in '
                                                 f'{exp_file} ({df_exp.shape[0]}) must match '
                                                 f'the calculated data in {calc_file} '
                                                 f'({df_calc.shape[1]})')

    log.append(f'# Reading {df_exp.shape[0]} experimental data from {exp_file} \n')
    log.append(f'# Reading {df_calc.shape[0]} calculated samples from {calc_file} \n')

    df_exp = df_exp.rename(columns={'val': 'avg'})
    df_exp = df_exp.rename(columns={'sigma':'sigma2'})

    # Define bounds
    df_exp['bound'] = 0
    if bound_string == 'UPPER':
        df_exp['bound'] = 1.0
    if bound_string == 'LOWER':
        df_exp['bound'] = -1.0

    labels = df_exp['label'].values
    exp = np.array(df_exp[['avg','sigma2','bound']].values)
    calc = np.array(df_calc.values)

    # Create log msg
    log_msg = ''.join(log)

    return labels, exp, calc, log_msg


def subsample(
    label: np.ndarray,
    exp: np.ndarray,
    calc: np.ndarray,
    use_samples: list | None = None,
    use_data: list | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Possibly remove some of the samples."""

    log = []
    if use_samples is not None:
        calc = calc[use_samples,:]
        log.append(f'# Using a subset of samples ({calc.shape[0]}) \n')
    if use_data is not None:
        label = label[use_data]
        exp = exp[use_data,:]
        calc = calc[:,use_data]
        log.append(f'# Using a subset of datapoints ({exp.shape[0]}) \n')

    # Create log msg
    log_msg = ''.join(log)

    return label,exp,calc,log_msg


def calc_chi(exp: np.ndarray, calc: np.ndarray, sample_weights: np.ndarray) -> float:
    """Calculate chi2."""

    calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                      axis=0)
    diff = calc_avg-exp[:,0]
    ii = np.where(((diff<0) & (exp[:,2]<0)) | ((diff>0) & (exp[:,2]>0)))[0]
    ff = [1 if (exp[j,2] == 0 or j in ii) else 0 for j in range(exp.shape[0])]
    diff *= ff # to_zero
    return  np.average((diff/exp[:,1])**2)


def check_data(
    labels: np.ndarray,
    exp: np.ndarray,
    calc: np.ndarray,
    sample_weights: np.ndarray,
    ) -> str:
    """Sanity check."""

    calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                      axis=0)

    log = []

    diff = calc_avg-exp[:,0]
    ii = np.where(((diff<0) & (exp[:,2]<0)) | ((diff>0) & (exp[:,2]>0)))[0]
    ff = [1 if exp[j,2] == 0 else 0 for j in range(exp.shape[0])]
    nviol = 0
    if len(ii) > 0:
        log.append(f'# The ensemble violates the following {len(ii)} boundaries: \n')
        log.append(f'# {"label":14s} {"exp_avg":8s} {"calc_avg":8s} \n')
        for j in ii:
            log.append(f'# {labels[j]:14s} {exp[j,0]:8.4f} {calc_avg[j]:8.4f} \n')
            ff[j] = 1
            nviol += 1
    diff *= ff # to_zero
    chi2 = np.average((diff / exp[:,1]) ** 2)
    rmsd = np.sqrt(np.average(diff**2))
    ii_out = np.where(diff > 1.)[0]
    nviol += len(ii_out)
    log.append(f'CHI2: {chi2:.5f} \n')
    log.append(f'RMSD: {rmsd:.5f} \n')
    log.append(f'VIOLATIONS: {nviol} \n')

    m_min = np.min(calc,axis=0)
    m_max = np.max(calc,axis=0)

    diff_min = ff * (m_min - exp[:,0]) / exp[:,1]
    ii_min = np.where(diff_min > 1.)[0]
    if len(ii_min) > 0:
        log.append('##### WARNING ########## \n')
        log.append('# The minimum value of the following data is higher than expt range: \n')
        log.append(f'# {"label":14s} {"exp_avg":8s} {"exp_sigma":8s} {"min_calc":8s} \n')
        for j in ii_min:
            log.append(f'# {labels[j]:15f} {exp[j,0]:8.4f} {exp[j,1]:8.4f} {m_min[j]:8.4f} \n')

    diff_max = ff * (exp[:,0] - m_max) / exp[:,1]
    ii_max = np.where(diff_max > 1.)[0]
    if len(ii_max) > 0:
        log.append('##### WARNING ########## \n')
        log.append('# The maximum value of the following data is higher than expt range: \n')
        log.append(f'# {"label":14s} {"exp_avg":8s} {"exp_sigma":8s} {"max_calc":8s} \n')
        for j in ii_max:
            log.append(f'# {labels[j]:15f} {exp[j,0]:8.4f} {exp[j,1]:8.4f} {m_max[j]:8.4f} \n')

    # Create log msg
    log_msg = ''.join(log)

    return log_msg


def standardize(
    exp: np.ndarray,
    calc: np.ndarray,
    sample_weights: np.ndarray,
    ) -> tuple[str,float,float]:
    """Standardize dataset. Modify array in place."""

    log = []

    calc_avg = np.sum(calc * sample_weights[:,np.newaxis],
                      axis=0)

    calc_var = np.average((calc_avg - calc)**2,
                          weights=sample_weights,
                          axis=0)

    #v2 = np.sqrt(np.average(np.array([calc_var,exp[:,1]]),axis=0)) # std
    calc_std = np.average(np.array([np.sqrt(calc_var),
                                    exp[:,1]]),
                          axis=0)

    exp[:,0] -= calc_avg
    exp[:,0] /= calc_std
    calc -= calc_avg
    calc /= calc_std
    exp[:,1] /= calc_std

    log.append('# Z-score normalization \n')

    # Create log msg
    log_msg = ''.join(log)

    return log_msg,calc_avg,calc_std
