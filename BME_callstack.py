# IMPORTS
## Standard Library Imports
import contextlib
import sys
import os
import time

## Third Party Imports
import numpy as np
import pandas as pd
from scipy.optimize import minimize #version >0.13.0 for dogleg optimizer
from scipy.special import logsumexp
from sklearn.linear_model import LinearRegression

## Local Imports
import BME_tools_callstack as bt

# FUNCTIONS
def progress(
    count,
    total,
    suffix=''):

    total -= 1
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write(f'[{bar}] {percents}% ...{suffix}\r')
    sys.stdout.flush()


class Reweight:
    """Reweight class.

    Returns:
        :
    """
    def __init__(
        self,
        name,
        w0=[]):
        """Initialize.

        Args:
            name:

            w0:
                . Defaults to [].
        """

        self.name = name

        if len(w0) != 0:
            self.w0 = w0/np.sum(w0)
        else:
            self.w0= []

        self.w_opt = []
        self.lambdas = []

        self.logfile = f'{name}.log'

        self.labels = []
        self.experiment =  []
        self.calculated =  []
        self.standardized = False

    def _write_log(
        self,
        msg):

        with open(self.logfile,'w',encoding='utf-8') as log_fh:
            log_fh.write(msg)

    def read_file(
        self,
        exp_file,
        calc_file,
        averaging='auto',
        fit='no',
        use_samples=[],
        use_data=[]):

        # read file
        log = ''
        label,exp,calc,log,averaging = bt.parse(exp_file,
                                                calc_file,
                                                averaging=averaging)
        self._write_log(log)

        # remove datapoints if use_samples or use_data is not empty
        label,exp, calc, log = bt.subsample(label,
                                            exp,
                                            calc,
                                            use_samples,
                                            use_data)
        self._write_log(log)

        if len(self.w0) == 0:
            self.w0 = np.ones(calc.shape[0]) / calc.shape[0]
            self._write_log(f'Initialized uniform weights {calc.shape[0]}\n')

        # fit/scale
        _, log = bt.fit_and_scale(exp,
                                  calc,
                                  self.w0,
                                  fit=fit)
        self._write_log(log)

        # do sanity checks
        log  = bt.check_data(label,
                             exp,
                             calc,
                             self.w0)
        self._write_log(log)

        return label,exp,calc

    def load(
        self,
        exp_file,
        calc_file,
        averaging='auto',
        fit='no',
        use_samples=[],
        use_data=[],
        weight=1):

        label, exp, calc = self.read_file(exp_file,
                                          calc_file,
                                          averaging=averaging,
                                          fit=fit,
                                          use_samples=use_samples,
                                          use_data=use_data)

        if len(self.experiment) == 0:
            self.experiment = exp
            self.calculated = calc
            self.labels = label
            self.weights = np.ones(exp.shape[0]) * weight
        else:
            self.experiment = np.vstack([self.experiment,exp])
            self.calculated = np.hstack([self.calculated,calc])
            self.labels = np.hstack([self.labels,label])
            self.weights = np.hstack([self.weights,np.ones(exp.shape[0]) * weight])
        # note to self: implement weight

    def load_array(
        self,
        label,
        exp,
        calc,
        weight=1):
        """Add data from external array.

        Args:
            label:

            exp:

            calc:

            weight:
                . Defaults to 1.
        """
        if len(self.experiment) == 0:
            self.experiment = exp
            self.calculated = calc
            self.labels = label
            if len(self.w0) == 0:
                self.w0 = np.ones(calc.shape[0]) / calc.shape[0]
            self.weights = np.ones(exp.shape[0]) * weight
        else:
            self.experiment = np.vstack([self.experiment,exp])
            self.calculated = np.hstack([self.calculated,calc])
            self.labels = np.hstack([self.labels,label])
            self.weights = np.hstack([self.weights,np.ones(exp.shape[0]) * weight])

    def predict_array(
        self,
        label,
        exp,
        calc,
        outfile=None,
        averaging='linear',
        fit='no'):
        stats, _ = bt.calc_stats(label,
                                 exp,
                                 calc,
                                 self.w0,
                                 self.w_opt,
                                 averaging=averaging,
                                 outfile=outfile,
                                 fit=fit)
        return stats

    def get_lambdas(self):
        return self.lambdas

    def get_iterations(self):
        return self.niter

    def get_nsamples(self):
        return self.calculated.shape[0]

    def get_ndata(self):
        return self.experiment.shape[0]

    def get_labels(self):
        return self.labels

    def get_experiment(self):
        return np.copy(self.experiment)

    def get_calculated(self):
        return np.copy(self.calculated)

    def get_name(self):
        return self.name

    def get_weights(self):
        return np.copy(self.w_opt)

    def get_w0(self):
        return np.copy(self.w0)

    def set_lambdas(self,lambda0):
        if len(self.lambdas) == 0:
            self.lambdas = lambda0
        else:
            print('# Overriding lambdas is not possible')
            sys.exit(1)

    def fit(
        self,
        theta,
        lambdas_init=True):
        """Optimize.

        Args:
            theta:

            lambdas_init:
                . Defaults to True.

        Returns:
            : 
        """

        if not self.standardized:
            bt.standardize(self.experiment,
                           self.calculated,
                           self.w0,
                           normalize='zscore')
            self.standardized = True

        def maxent(lambdas):
            # weights
            arg = -np.sum(lambdas[np.newaxis,:] * self.calculated,axis=1) - tmax + np.log(self.w0)

            logz = logsumexp(arg)
            ww = np.exp(arg - logz)

            avg = np.sum(ww[:,np.newaxis] * self.calculated,
                         axis=0)

            # gaussian integral
            eps2 = 0.5 * np.sum((lambdas * lambdas) * theta_sigma2)

            # experimental value
            sum1 = np.dot(lambdas,self.experiment[:,0])
            fun = sum1 + eps2 + logz

            # gradient
            jac = self.experiment[:,0] + lambdas * theta_sigma2 - avg

            # divide by theta to avoid numerical problems
            return  fun/theta, jac/theta

        if lambdas_init:
            lambdas = np.zeros(self.experiment.shape[0],
                               dtype=np.longdouble)

            self._write_log('Lagrange multipliers initialized from zero\n')
        else:
            assert len(self.lambdas) == self.experiment.shape[0]
            lambdas = np.copy(self.lambdas)
            self._write_log('Warm start\n')

        bounds = []
        for j in range(self.experiment.shape[0]):
            if self.experiment[j,2] == 0:
                bounds.append([None,None])
            elif self.experiment[j,2] == -1:
                bounds.append([None,0.0])
            else:
                bounds.append([0.0,None])

        opt = {'maxiter': 50000,
               'disp': False}

        tmax = np.log((sys.float_info.max) / 5.)

        theta_sigma2 = theta * self.weights * self.experiment[:,1]**2

        chi2_before  = bt.calc_chi(self.experiment,
                                   self.calculated,
                                   self.w0)

        self._write_log((f'Optimizing {self.experiment.shape[0]} data and '
                         f'{self.calculated.shape[0]} samples. Theta={theta} \n'))
        self._write_log(f'CHI2 before optimization: {chi2_before:8.4f} \n')

        mini_method = 'L-BFGS-B'
        start_time = time.time()

        result = minimize(maxent,
                          lambdas,
                          options=opt,
                          method=mini_method,
                          jac=True,
                          bounds=bounds)

        self._write_log(f'Execution time: {(time.time() - start_time):.2f} seconds\n')

        if result.success:
            self._write_log((f'Minimization using {mini_method} successful '
                             f'(iterations:{result.nit})\n'))
            arg = -np.sum(result.x[np.newaxis,:] * self.calculated,axis=1) - tmax
            w_opt = self.w0 * np.exp(arg)
            w_opt /= np.sum(w_opt)
            self.lambdas = np.copy(result.x)
            self.w_opt = np.copy(w_opt)
            self.niter = result.nit
            chi2_after = bt.calc_chi(self.experiment,
                                     self.calculated,
                                     w_opt)
            phi = np.exp(-bt.srel(self.w0,
                                  w_opt))

            self._write_log(f'CHI2 after optimization: {chi2_after:8.4f} \n')
            self._write_log(f'Fraction of effective frames: {phi:8.4f} \n')

            return chi2_before, chi2_after, phi

        else:
            self._write_log(f'Minimization using {mini_method} failed\n')
            self._write_log(f'Message: {result.message}\n')
            self.niter = -1
            return np.NaN, np.NaN, np.NaN

    def ibme(
        self,
        theta,
        ftol=0.01,
        iterations=50,
        lr_weights=True,
        offset=True):

        current_weights = self.get_w0()
        w0 = self.get_w0()
        name = self.get_name()
        labels = self.get_labels()
        exp = self.get_experiment()
        calc = self.get_calculated()

        self.ibme_weights = []
        self.ibme_stats = []

        if lr_weights:
            inv_var = 1. / exp[:,1]**2
        else:
            inv_var = np.ones(len(exp))

        log = []
        rr_old = np.NaN

        for it in range(iterations):

            calc_avg = np.sum(calc * current_weights[:,np.newaxis],
                              axis=0)

            model = LinearRegression(fit_intercept=offset)
            model.fit(calc_avg.reshape(-1,1),
                      exp[:,0],
                      inv_var)

            alpha = model.coef_[0] # scale factor
            beta = model.intercept_
            calc = alpha * calc + beta

            r1 = Reweight(f'{name}_ibme_{it}',
                          w0=np.copy(w0))
            r1.load_array(labels,
                          np.copy(exp),
                          np.copy(calc))
            rr = r1.fit(theta=theta)

            if it == 0:
                chi2_0 = rr[0]
                calc_0 = np.copy(calc)

            current_weights = np.copy(r1.get_weights())

            diff = abs(rr_old - rr[1])
            rr_old = rr[1]

            log.append((f'Iteration:{it:3d} scale: {alpha:7.4f} offset: {beta:7.4f} '
                        f'chi2: {rr[1]:7.4f} diff: {diff:7.4e}\n'))
            self.ibme_weights.append(current_weights)
            self.ibme_stats.append(rr)

            if diff < ftol:
                line = (f'Iterative procedure converged below tolerance {diff:.2e} '
                        f'after {it} iterations\n')
                print(line,end='')
                log.append(line)
                break

        self._write_log(''.join(log)+ '\n')

        n1 = f'{self.name}_{it}.calc.dat'
        n2 = f'{self.name}_{it}.weights.dat'

        df = pd.DataFrame(calc)
        df.to_csv(n1,
                  sep=' ',
                  header=False,
                  float_format='%8.4e')

        df = pd.DataFrame(current_weights)
        df.to_csv(n2,
                  sep=' ',
                  header=False,
                  float_format='%8.4e')

        phi = np.exp(-bt.srel(w0,
                              current_weights))
        self.w_opt = current_weights

        return chi2_0,rr[1],phi,calc_0,calc

    def get_ibme_weights(self):

        try:
            return self.ibme_weights
        except AttributeError:
            print('# iBME weights not available. Call iBME first')
            sys.exit(1)

    def get_ibme_stats(self):

        try:
            return self.ibme_stats
        except AttributeError:
            print('# iBME stats not available. Call iBME first')
            sys.exit(1)


def myibme(
    theta: int,
    exp_file: str,
    calc_file: str,
    output_dir: str,
    ) -> tuple[int,tuple[float,float,float],np.ndarray]:
    """Perform the Iterative Bayesian Maximum Entropy (BME) algorithm on calculated SAXS data,
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

    """
    # Change current working directory
    old_cd = os.getcwd()
    os.chdir(output_dir)

    # Create reweight object
    rew = Reweight(f'ibme_t{theta}')

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

if __name__ == '__main__':
    THETA = 100
    EXP_FILE = '/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/REWEIGHTING/Hst5/Hst5_exp_saxs.dat'
    CALC_FILE = '/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/REWEIGHTING/Hst5/Hst5_calc_saxs.dat'
    OUTPUT_DIR = '/home/tiagogomes/reweighting_callstack'
    myibme(theta=THETA,exp_file=EXP_FILE,calc_file=CALC_FILE,output_dir=OUTPUT_DIR)
