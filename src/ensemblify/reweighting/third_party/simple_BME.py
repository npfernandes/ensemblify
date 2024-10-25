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
## Standard Library Imports
import sys
import time

## Third Party Imports
import numpy as np
import pandas as pd
from scipy.optimize import minimize #version >0.13.0 for dogleg optimizer
from scipy.special import logsumexp
from sklearn.linear_model import LinearRegression

## Local Imports
import ensemblify.reweighting.third_party.simple_BME_tools as bt

# CLASS
class SimpleReweight:
    """SimpleReweight class."""
    def __init__(self, name: str, w0: list | None = None):
        """Initialize SimpleReweight instance.

        Args:
            name:
                name of current instance, will appear in created files.
            w0:
                initial weights. Defaults to None which will lead to uniform weight initialization.
        """

        self.name = name
        self.w0 = w0

        self.w_opt = None
        self.lambdas = None

        self.logfile = f'{name}.log'

        self.labels = None
        self.experiment =  None
        self.calculated =  None
        self.standardized = False

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

    def set_lambdas(self, lambda0: np.ndarray):
        if self.lambdas is None:
            self.lambdas = lambda0
        else:
            print('# Overriding lambdas is not possible')
            sys.exit(1)

    def _write_log(self, msg: str):
        """Write message to log file."""
        with open(self.logfile,'w',encoding='utf-8') as log_fh:
            log_fh.write(msg)

    def read_file(
        self,
        exp_file: str,
        calc_file: str,
        use_samples: list | None = None,
        use_data: list | None = None,
        ) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
        """Read experimental and calculated data files.

        Args:
            exp_file:
                path to file with experimental data.
            calc_file:
                path to file with calculated calculated data.
            use_samples:
                Use only this subset of calculated data indices. Defaults to None, and all
                samples are used.
            use_data:
                Use only this subset of experimental data indices. Defaults to None, and all
                samples are used.

        Returns:
            label,exp,calc:
                data labels, experimental data, calculated data.
        """

        # read file
        log = ''
        labels, exp, calc, log = bt.parse(exp_file,
                                         calc_file)
        self._write_log(log)

        # remove datapoints if use_samples or use_data is not empty
        labels, exp, calc, log = bt.subsample(labels,
                                              exp,
                                              calc,
                                              use_samples,
                                              use_data)
        self._write_log(log)

        # Initialize uniform weights
        if self.w0 is None:
            self.w0 = np.ones(calc.shape[0]) / calc.shape[0]
            self._write_log(f'Initialized uniform weights {calc.shape[0]}\n')
        else:
            self._write_log('Warm start\n')

        # do sanity checks
        log  = bt.check_data(labels,
                             exp,
                             calc,
                             self.w0)
        self._write_log(log)

        return labels,exp,calc

    def load(
        self,
        exp_file: str,
        calc_file: str,
        use_samples: list | None = None,
        use_data: list | None = None,
        weight: int = 1,
        ) -> None:
        """Load data from files into class attributes.

        Args:
            exp_file:
                path to file with experimental data.
            calc_file:
                path to file with calculated data.
            use_samples:
                Use only this subset of calculated data indices. Defaults to None, and all
                samples are used.
            use_data:
                Use only this subset of experimental data indices. Defaults to None, and all
                samples are used.
            weight:
                value to multiply all weights by. Defaults to 1.
        """

        labels, exp, calc = self.read_file(exp_file,
                                          calc_file,
                                          use_samples=use_samples,
                                          use_data=use_data)

        if self.experiment is None:
            self.experiment = exp
            self.calculated = calc
            self.labels = labels
            self.weights = np.ones(exp.shape[0]) * weight
        else:
            self.experiment = np.vstack([self.experiment,exp])
            self.calculated = np.hstack([self.calculated,calc])
            self.labels = np.hstack([self.labels,labels])
            self.weights = np.hstack([self.weights,np.ones(exp.shape[0]) * weight])
        # note to self: implement weight

    def load_array(
        self,
        labels: np.ndarray,
        exp: np.ndarray,
        calc: np.ndarray,
        weight: int = 1,
        ) -> None:
        """Load data from external array into class attributes.

        Args:
            labels:
                array of data labels.
            exp:
                array of experimental data values.
            calc:
                array of calculated data values.
            weight:
                value to multiply all weights by. Defaults to 1.
        """
        if self.experiment is None:
            self.experiment = exp
            self.calculated = calc
            self.labels = labels
            if self.w0 is None:
                self.w0 = np.ones(calc.shape[0]) / calc.shape[0]
            self.weights = np.ones(exp.shape[0]) * weight
        else:
            self.experiment = np.vstack([self.experiment,exp])
            self.calculated = np.hstack([self.calculated,calc])
            self.labels = np.hstack([self.labels,labels])
            self.weights = np.hstack([self.weights,np.ones(exp.shape[0]) * weight])

    def fit(self, theta: int) -> tuple[float, float, float]:
        """Optimize.

        Args:
            theta:
                theta value to use in reweighting.

        Returns:
            chi2_before,chi2_after,phi:
                chisquare value before and after fitting and fraction of effective frames.
        """

        if not self.standardized:
            bt.standardize(self.experiment,
                           self.calculated,
                           self.w0)
            self.standardized = True

        def maxent(lambdas: np.ndarray) -> tuple[float,float]:
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

        lambdas = np.zeros(self.experiment.shape[0],
                           dtype=np.longdouble)
        self._write_log('Lagrange multipliers initialized from zero\n')

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
        theta: int,
        ftol: float = 0.01,
        iterations: int = 50,
        offset: bool = True,
        ) -> tuple[float | float | float | np.ndarray | np.ndarray]:
        """_summary_

        Args:
            theta:
                theta value to use in iBME iterations.
            ftol:
                tolerance for minimization procedure. Defaults to 0.01.
            iterations:
                number of iBME iterations to perform. Defaults to 50.
            offset:
                whether to offset calculated data at each step when fitting it to experimental
                data. Defaults to True.

        Returns:
            chi2_0,rr[1],phi,calc_0,calc:
                Initial chisquare value, final chisquare value, fraction of effective frames,
                initial calculated data, final calculated data.
        """

        current_weights = self.get_w0()
        w0 = self.get_w0()
        name = self.get_name()
        labels = self.get_labels()
        exp = self.get_experiment()
        calc = self.get_calculated()

        self.ibme_weights = []
        self.ibme_stats = []

        # Setup inverse variance (weights)
        inv_var = 1. / exp[:,1]**2

        log = []
        rr_old = np.NaN

        for it in range(iterations):

            calc_avg = np.sum(calc * current_weights[:,np.newaxis],
                              axis=0)

            model = LinearRegression(fit_intercept=offset)
            model.fit(X=calc_avg.reshape(-1,1),
                      y=exp[:,0],
                      sample_weight=inv_var)

            alpha = model.coef_[0] # scale factor
            beta = model.intercept_
            calc = alpha * calc + beta

            r1 = SimpleReweight(name=f'{name}_ibme_{it}',
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
