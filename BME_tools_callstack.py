# IMPORTS
## Third Party Imports
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# CONSTANTS
EXP_TYPES = ['NOE','JCOUPLINGS','CS','SAXS','RDC']
BOUND_TYPES = ['UPPER','LOWER']
AVERAGING_TYPES = ['linear','power_3','power_6','power_4']
FIT_TYPES = ['scale','scale+offset','no']

# FUNCTIONS
def srel(
    w0,
    w1,
    ):
    """Calculate relative entropy.

    Args:
        w0:

        w1:


    Returns:
        rel_ent:

    """
    idxs = np.where(w1 > 1.0E-50)
    return np.sum(w1[idxs] * np.log(w1[idxs] / w0[idxs]))


def parse(
    exp_file,
    calc_file,
    averaging='auto',
    ):
    """Parse files.

    Args:
        exp_file:

        calc_file:

        averaging:
            . Defaults to 'auto'.

    Returns:
        A tuple X where:
            :

            :

            :

            :

    """
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
    bound_string=None
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

    # Determine averaging
    if averaging == 'auto':
        if data_string == 'NOE':
            averaging = 'power_6'
        else:
            averaging = 'linear'
    else:
        assert averaging in AVERAGING_TYPES, f'averaging type must be in {AVERAGING_TYPES} '
    log.append(f'# Using {averaging} averaging \n')

    if averaging.split('_')[0] == 'power':
        noe_power = int(averaging.split('_')[-1])
        df_exp['avg'] = np.power(df_exp['val'],
                                 -noe_power)

        df_exp['sigma2'] = noe_power * df_exp['avg']*df_exp['sigma'] / df_exp['val']
        df_calc = np.power(df_calc,
                           -noe_power)

        # if bound constraints, swap lower and upper
        if bound_string == 'LOWER':
            bound_string = 'UPPER'
        elif bound_string == 'UPPER':
            bound_string = 'LOWER'
    else:
        df_exp = df_exp.rename(columns={'val': 'avg'})
        df_exp = df_exp.rename(columns={'sigma':'sigma2'})

    # Define bounds
    df_exp['bound'] = 0
    if bound_string == 'UPPER': df_exp['bound'] = 1.0
    if bound_string == 'LOWER': df_exp['bound'] = -1.0
    #df_exp["tag"] = exp_file

    labels = df_exp['label'].values
    exp = np.array(df_exp[['avg','sigma2','bound']].values)
    calc = np.array(df_calc.values)

    # Create log msg
    log_msg = ''.join(log)

    return labels,exp,calc,log_msg,averaging


def fit_and_scale(
    exp,
    calc,
    sample_weights,
    fit):
    """Perform linear regression.

    Args:
        exp:

        calc:

        sample_weights:

        fit:


    Returns:
        :
    """
    assert fit in FIT_TYPES, f'fit type must be in {FIT_TYPES} '

    calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                      axis=0)

    log = [f'# Using {fit} scaling \n']

    if fit == 'no':
        # Create log msg
        log_msg = ''.join(log)
        return calc_avg,log_msg
    else:
        exp_avg = exp[:,0]

        if fit == 'scale':
            fit_intercept = False
        else:
            fit_intercept = True

        oversigma = (1./exp[:,1]**2)

        reg = LinearRegression(fit_intercept=fit_intercept).fit(calc_avg.reshape(-1,1),
                                                                exp_avg.reshape(-1,1),
                                                                sample_weight=oversigma)

        r_value = reg.score(calc_avg.reshape(-1,1),
                            exp_avg.reshape(-1,1),
                            sample_weight=oversigma)

        slope,intercept = reg.coef_[0],reg.intercept_

        calc *=slope
        calc +=intercept

        calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                          axis=0).reshape(-1,1)

        log.append(f'# Slope={slope:8.4f}; Offset={intercept:8.4f}; r2={r_value:8.4f} \n')

        # Create log msg
        log_msg = ''.join(log)

        return calc_avg,log_msg


def subsample(
    label,
    exp,
    calc,
    use_samples,
    use_data):
    """Remove some of the samples.

    Args:
        label:

        exp:

        calc:

        use_samples:

        use_data:


    Returns:
        :
    """

    log = []
    if len(use_samples) != 0:
        calc = calc[use_samples,:]
        log.append(f'# Using a subset of samples ({calc.shape[0]}) \n')
    if len(use_data) != 0:
        label = label[use_data]
        exp = exp[use_data,:]
        calc = calc[:,use_data]
        log.append(f'# Using a subset of datapoints ({exp.shape[0]}) \n')

    # Create log msg
    log_msg = ''.join(log)

    return label,exp,calc,log_msg


def calc_chi(
    exp,
    calc,
    sample_weights):
    """Calculate chi2.

    Args:
        exp:

        calc:

        sample_weights:


    Returns:
        :
    """
    calc_avg = np.sum(calc*sample_weights[:,np.newaxis],
                      axis=0)
    diff = calc_avg-exp[:,0]
    ii = np.where(((diff<0) & (exp[:,2]<0)) | ((diff>0) & (exp[:,2]>0)))[0]
    ff = [1 if (exp[j,2] == 0 or j in ii) else 0 for j in range(exp.shape[0])]
    diff *= ff # to_zero
    return  np.average((diff/exp[:,1])**2)


def check_data(
    label,
    exp,
    calc,
    sample_weights):
    """Sanity check.

    Args:
        label:

        exp:

        calc:

        sample_weights:


    Returns:
        :
    """

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
            log.append(f'# {label[j]:14s} {exp[j,0]:8.4f} {calc_avg[j]:8.4f} \n')
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
            log.append(f'# {label[j]:15f} {exp[j,0]:8.4f} {exp[j,1]:8.4f} {m_min[j]:8.4f} \n')

    diff_max = ff * (exp[:,0] - m_max) / exp[:,1]
    ii_max = np.where(diff_max > 1.)[0]
    if len(ii_max) > 0:
        log.append('##### WARNING ########## \n')
        log.append('# The maximum value of the following data is higher than expt range: \n')
        log.append(f'# {"label":14s} {"exp_avg":8s} {"exp_sigma":8s} {"max_calc":8s} \n')
        for j in ii_max:
            log.append(f'# {label[j]:15f} {exp[j,0]:8.4f} {exp[j,1]:8.4f} {m_max[j]:8.4f} \n')

    # Create log msg
    log_msg = ''.join(log)

    return log_msg


def calc_diffs(
    exp_avg,
    exp_sigma,
    bounds,
    calc_avg):
    """Calculate Chi2, RMSD, number of violations.

    Args:
        exp_avg:

        exp_sigma:

        bounds:

        calc_avg:


    Returns:
        :
    """

    diff = (calc_avg - exp_avg) / exp_sigma
    # boundaries violations
    ii = np.where(((diff < -1.) & (bounds < 0)) | ((diff > 1.) & (bounds > 0)))[0]
    ff = [1 if (bounds[j]==0 or j in ii) else 0 for j in range(exp_avg.shape[0])]
    # other violations
    viol = np.where(((diff**2 > 1) & (bounds == 0)))[0]
    #to_zero = np.zeros(len(diff))
    viol_idx = np.array(list(ii) + list(viol))
    nviol = viol_idx.shape[0]
    diff *= ff # to_zero
    chi2 = np.average((diff)**2)
    rmsd = np.sqrt(np.average((diff * exp_sigma)**2))
    return chi2,rmsd,nviol,viol_idx


def calc_stats(
    label,
    exp,
    calc,
    sample_weights0,
    sample_weights1,
    averaging,
    fit,
    outfile):
    """Calculate statistics.

    Args:
        label:

        exp:

        calc:

        sample_weights0:

        sample_weights1:

        averaging:

        fit:

        outfile:


    Returns:
        : 
    """
    assert averaging in AVERAGING_TYPES, f'averaging type must be in {AVERAGING_TYPES} '

    calc_avg0 = np.sum(calc * sample_weights0[:,np.newaxis],
                       axis=0)
    calc_avg1 = np.sum(calc * sample_weights1[:,np.newaxis],
                       axis=0)

    exp_avg = exp[:,0]
    oversigma = 1. / exp[:,1]**2
    #oversigma = np.ones(len(exp))
    intercept0, intercept1 = 0., 0.
    slope0, slope1 = 1., 1.

    if fit == 'scale':
        slope0 = np.dot(exp_avg,calc_avg0) / np.dot(calc_avg0,calc_avg0)
        slope1 = np.dot(exp_avg,calc_avg1) / np.dot(calc_avg1,calc_avg1)
    elif fit == 'scale+offset':
        #slope0, intercept0, r_value0, p_value0, std_err0 = linregress(calc_avg0[:,0],exp[:,0])
        #slope1, intercept1, r_value1, p_value1, std_err1 = linregress(calc_avg1[:,0],exp[:,0])
        reg0 = LinearRegression(fit_intercept=True).fit(calc_avg0.reshape(-1,1),
                                                        exp_avg.reshape(-1,1),
                                                        sample_weight=oversigma)
        slope0, intercept0 = reg0.coef_[0], reg0.intercept_

        reg1 = LinearRegression(fit_intercept=True).fit(calc_avg1.reshape(-1,1),
                                                        exp_avg.reshape(-1,1),
                                                        sample_weight=oversigma)
        slope1, intercept1 = reg1.coef_[0], reg1.intercept_

    log = [f'# Slope={slope0:8.4f}/{slope1:8.4f}; Offset={intercept0:8.4f}/{intercept1:8.4f} \n']

    if averaging.split('_')[0] == 'power':
        noe_power = int(averaging.split('_')[-1])

        calc_avg0 = np.power(np.sum((slope0 * calc + intercept0) * sample_weights0[:,np.newaxis],
                                    axis=0),
                             -1/noe_power)

        calc_avg1 = np.power(np.sum((slope1 * calc + intercept1) * sample_weights1[:,np.newaxis],
                                    axis=0),
                             -1/noe_power)

        exp_avg = np.power(exp[:,0],
                           -1/noe_power)

        exp_sigma = (exp_avg * exp[:,1]) / (noe_power * exp[:,0])
        bounds = -exp[:,2]
    else:
        calc_avg0 = np.sum((slope0 * calc + intercept0) * sample_weights0[:,np.newaxis],
                           axis=0)
        calc_avg1 = np.sum((slope1 * calc + intercept1) * sample_weights1[:,np.newaxis],
                           axis=0)

        exp_avg = exp[:,0]
        exp_sigma = exp[:,1]
        bounds = exp[:,2]

    log.append(f'# {averaging} averaging \n')

    chi2_0, rmsd_0, nviol_0, viol_idx_0 = calc_diffs(exp_avg,
                                                     exp_sigma,
                                                     bounds,
                                                     calc_avg0)

    chi2_1, rmsd_1, nviol_1, viol_idx_1 = calc_diffs(exp_avg,
                                                     exp_sigma,
                                                     bounds,
                                                     calc_avg1)

    if outfile is not None:
        violation0 = [1 if k in viol_idx_0 else 0 for k in range(len(calc_avg0))]
        violation1 = [1 if k in viol_idx_1 else 0 for k in range(len(calc_avg1))]
        violation = [f'{i1}{i2}' for i1,i2 in zip(violation0,violation1)]

        df = pd.DataFrame({'label': label,
                           'exp_avg': exp_avg,
                           'exp_sigma': exp_sigma,
                           'calc_avg': calc_avg0,
                           'calc_avg_rew': calc_avg1,
                           'violation': violation})

        with open(outfile, 'w', encoding='utf-8') as fh:
            fh.write(f'# {" ".join(list(df.columns))} \n')

        df['label'] = df['label'].map(lambda x: f'{x:-15}')
        df.to_csv(outfile,
                  index=False,
                  sep='\t',
                  header=False,
                  float_format='%8.3e',
                  mode='a')

        with open(outfile, 'a', encoding='utf-8') as fh:
            fh.write(f'# CHI2:       {chi2_0:8.4f} {chi2_1:8.4f} \n')
            fh.write(f'# RMSD:       {rmsd_0:8.4f} {rmsd_1:8.4f} \n')
            fh.write(f'# Violations: {nviol_0:8d}  {nviol_1:8d} \n')
            fh.write(f'# ndata:      {exp.shape[0]:8d} {calc.shape[0]:8d} \n')

    stats = [chi2_0, rmsd_0, nviol_0, chi2_1, rmsd_1, nviol_1]

    # Create log msg
    log_msg = ''.join(log)

    return stats,log_msg


def standardize(
    exp,
    calc,
    sample_weights,
    normalize='zscore'):
    """Standardize dataset. Modify array in place.

    Args:
        exp:

        calc:

        sample_weights:

        normalize:
            . Defaults to 'zscore'.

    Returns:
        :
    """

    log = []

    if normalize == 'zscore':
        v1 = np.sum(calc * sample_weights[:,np.newaxis],
                    axis=0)
        calc_var = np.average((v1 - calc)**2,
                              weights=sample_weights,
                              axis=0)
        #v2 = np.sqrt(np.average(np.array([calc_var,exp[:,1]]),axis=0)) # std
        v2 = np.average(np.array([np.sqrt(calc_var),
                                  exp[:,1]]),
                        axis=0)
        exp[:,0] -= v1
        exp[:,0] /= v2
        calc -= v1
        calc /= v2
        exp[:,1] /= v2

        log.append('# Z-score normalization \n')

    elif normalize == 'minmax':
        mmin = calc.min(axis=0)
        mmax = calc.max(axis=0)
        delta = mmax-mmin
        exp[:,0] -= mmin
        exp[:,0] /= delta
        calc -= mmin
        calc /= delta
        exp[:,1] /= delta
        log.append('# MinMax normalization \n')

    # Create log msg
    log_msg = ''.join(log)

    return log_msg,v1,v2
