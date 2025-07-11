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
from collections.abc import Callable

## Third Party Imports
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from sklearn.linear_model import LinearRegression

## Local Imports
import ensemblify.reweighting.third_party.simple_BME_tools as bt

# CONSTANTS
TMAX = np.log((sys.float_info.max) / 5.)
"""Maximum safe argument for exponential functions. Used to prevent overflow when computing 
exp(arg) by ensuring the result stays below sys.float_info.max with a safety factor of 5."""

EXP2FIT = {'SAXS': 'scale+offset',
           'RDC': 'scale',
           'CS': None,
           'JCOUPLINGS': None}
"""Mapping of experimental data types to data scaling methods."""

# CLASS
class SimpleReweight:
    """The SimpleReweight class for performing Bayesian/Maximum Entropy reweighting.

    It supports plain averaged experimental data like Chemical Shifts (CS) and 3J couplings
    (JCOUPLINGS).
    It supports data that must be scaled or scaled+offset, like Residual Dipolar Couplings
    (RDCs) and Small Angle X-Ray Scattering (SAXS) respectively.
    It does NOT support non linearly averaged data or data that comes in the form of lower and
    upper bounds like nuclear Overhauser effects (NOEs).

    Attributes:
        name (str):
            Identifier for the instance, used in created file names.
        w0 (np.ndarray | None):
            Initial weights. If None, uniform weights are initialized.
        w_opt (np.ndarray | None):
            Optimized weights after fitting to experimental data.
        lambdas (np.ndarray | None):
            Lagrange multipliers used in the optimization.
        log_fd (TextIOWrapper[_WrappedBuffer]):
            Open file descriptor for logging output.
        labels (np.ndarray | None):
            Labels for the experimental data.
        experiment (np.ndarray | None):
            Experimental data values.
        calculated (np.ndarray | None):
            Back-calculated experimental data values.
        normalized (bool):
            Indicates if the data has been normalized.
        exp_types (dict):
            Mapping of data point indices to their experimental data type. For example, for a set
            of 10 data points:

                { 'SAXS': (0,4), 'RDC': (5,10) }
                
            indicates that the first five data points correspond to SAXS data and the final five
            data points correspond to RDC data.
    """
    def __init__(self, name: str, w0: np.ndarray | None = None):
        """Initialize SimpleReweight instance.

        Args:
            name (str):
                Identifier for the instance, used in created file names.
            w0 (np.ndarray | None):
                Initial weights. If None, uniform weights are initialized.
        """
        # Setup input attributes
        self.name = name
        self.w0 = w0

        # Setup other attributes
        self.w_opt = None 
        self.lambdas = None
        self.log_fd = open(f'{name}.log','w',encoding='utf-8')
        self.labels = None
        self.experiment =  None
        self.calculated =  None
        self.normalized = False
        self.exp_types = {}

    def get_lambdas(self) -> np.ndarray | None:
        return self.lambdas

    def get_iterations(self) -> int:
        return self.niter

    def get_nsamples(self) -> int:
        return self.calculated.shape[0]

    def get_ndata(self) -> int:
        return self.experiment.shape[0]

    def get_labels(self) -> np.ndarray | None:
        return self.labels

    def get_experiment(self) -> np.ndarray | None:
        return np.copy(self.experiment)

    def get_calculated(self) -> np.ndarray | None:
        return np.copy(self.calculated)

    def get_name(self) -> str:
        return self.name

    def get_weights(self) -> np.ndarray | None:
        return np.copy(self.w_opt)

    def get_w0(self) -> np.ndarray | None:
        return np.copy(self.w0)

    def add_exp_type_idxs(self, exp_type: str, idxs: tuple[int, int]):
        """Add a new experimental data type to the exp_types dictionary.

        Args:
            exp_type (str):
                Type of experimental data to add.
            idxs (tuple[int, int]):
                Indices of the experimental data points corresponding to this type.
        """
        assert exp_type in bt.EXP_TYPES , (f'Error. Experimental type must be one of the '
                                           f'following: {bt.EXP_TYPES} ')
        self.exp_types[exp_type] = idxs

    def set_lambdas(self, lambda0: np.ndarray):
        if self.lambdas is None:
            self.lambdas = lambda0
        else:
            print('# Overriding lambdas is not possible')
            sys.exit(1)

    def log_message(self, message: str) -> None:
        """Helper method to write and flush log messages."""
        self.log_fd.write(message)
        self.log_fd.flush()

    def __del__(self):
        """Close log file when object is destroyed."""
        if hasattr(self, 'log_fd') and not self.log_fd.closed:
            self.log_fd.close()

    def read_file(
        self,
        exp_file: str,
        calc_file: str,
        exp_type: str,
        use_samples: list[int] | None = None,
        use_data: list[int] | None = None,
        ) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
        """Extract experimental and calculated data from files.

        Calculated data is fit to experimental data, then scaled and offset if applicable.

        Args:
            exp_file (str):
                Path to file with experimental data.
            calc_file (str):
                Path to file with back-calculated experimental data.
            exp_type (str):
                Type of experimental data, used to determine how to fit the calculated data.
            use_samples (list[int] | None):
                Use only this subset of calculated data indices. Defaults to None, and all
                samples are used.
            use_data (list[int] | None):
                Use only this subset of experimental data indices. Defaults to None, and all
                samples are used.

        Returns:
            tuple[np.ndarray,np.ndarray,np.ndarray]:
                np.ndarray:
                    experimental data labels.
                np.ndarray:
                    experimental data values.
                np.ndarray:
                    calculated data values.
        """
        # Read file
        labels, exp, calc, log = bt.parse(exp_file=exp_file,
                                          calc_file=calc_file)
        self.log_message(log)

        # Use only a subset of datapoints if use_samples or use_data is specified
        labels, exp, calc, log = bt.subsample(labels=labels,
                                              exp=exp,
                                              calc=calc,
                                              use_samples=use_samples,
                                              use_data=use_data)
        self.log_message(log)

        # Initialize weights, by default uniform
        if self.w0 is None:
            self.w0 = np.ones(calc.shape[0]) / calc.shape[0]
            self.log_message(f'Initialized uniform weights {calc.shape[0]}\n')
        else:
            self.log_message('Using provided starting weights\n')

        # Fit calculated data to experimental data, scale appropriately
        assert exp_type in bt.EXP_TYPES , (f'Error. DATA in {exp_file} must be one of the '
                                           f'following: {bt.EXP_TYPES} ')
        fit_type = EXP2FIT[exp_type]
        scaled_calc, log = bt.fit_and_scale(exp=exp,
                                            calc=calc,
                                            calc_weights=self.w0,
                                            fit_type=fit_type)
        self.log_message(log)

        # Do sanity checks and comparisons between experimental and calculated data
        log = bt.check_data(labels=labels,
                            exp=exp,
                            calc=calc,
                            sample_weights=self.w0)
        self.log_message(log)

        return labels, exp, calc

    def load(
        self,
        exp_file: str,
        calc_file: str,
        exp_type: str | None = None,
        use_samples: list[int] | None = None,
        use_data: list[int] | None = None,
        ) -> None:
        """Load data from files into class attributes.

        Args:
            exp_file (str):
                Path to file with experimental data.
            calc_file (str):
                Path to file with calculated data.
            exp_type (str):
                Type of experimental data, used to determine how to fit the calculated data.
            use_samples (list[int] | None):
                Use only this subset of calculated data indices. Defaults to None, and all
                samples are used.
            use_data (list[int] | None):
                Use only this subset of experimental data indices. Defaults to None, and all
                samples are used.
        """
        # Check if experimental data type is provided
        if exp_type is None:
            # Read first line of experimental data file to get type
            with open(exp_file, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
                if 'DATA=' in first_line:
                    for equality in first_line.split(' '):
                        if 'DATA=' in equality:
                            exp_type = equality.split('=')[1].upper()
                            break
                else:
                    raise ValueError('Experimental data type not specified in the first line of '
                                    'the experimental data file. Please add it or specify it '
                                    'using the exp_type argument.')

        assert exp_type in bt.EXP_TYPES , (f'Error. DATA in {exp_file} must be one of the '
                                           f'following: {bt.EXP_TYPES} ')
        
        # Read files and parse data
        labels, exp, calc = self.read_file(exp_file=exp_file,
                                           calc_file=calc_file,
                                           exp_type=exp_type,
                                           use_samples=use_samples,
                                           use_data=use_data)

        # If its the first time loading data, initialize class attributes
        if self.experiment is None:
            self.experiment = exp
            self.calculated = calc
            self.labels = labels
            self.add_exp_type_idxs(exp_type, (0, exp.shape[0] - 1))
        else:
            # If data is already loaded, append new data to existing arrays
            self.experiment = np.vstack([self.experiment,exp])
            self.calculated = np.hstack([self.calculated,calc])
            self.labels = np.hstack([self.labels,labels])
            self.add_exp_type_idxs(exp_type, (self.experiment.shape[0] - exp.shape[0],
                                              self.experiment.shape[0] - 1)
                                    )

    def load_array(
        self,
        labels: np.ndarray,
        exp: np.ndarray,
        calc: np.ndarray):
        """Load data from external array directly into class attributes.

        Args:
            labels (np.ndarray):
                array of data labels.
            exp (np.ndarray):
                array of experimental data values.
            calc (np.ndarray):
                array of calculated data values.
        """
        # If its the first time loading data, initialize class attributes
        if self.experiment is None:
            self.experiment = exp
            self.calculated = calc
            self.labels = labels
            if self.w0 is None:
                self.w0 = np.ones(calc.shape[0]) / calc.shape[0]
        else:
            # If data is already loaded, append new data to existing arrays
            self.experiment = np.vstack([self.experiment,exp])
            self.calculated = np.hstack([self.calculated,calc])
            self.labels = np.hstack([self.labels,labels])


    def define_maxent(self,
                      theta: int) -> Callable:
        """Define the maximum entropy function for optimization.
        Args:
            theta (int):
                value for the theta parameter to use in defining the maximum entropy function.
        
        Returns:
            function:
                A function that takes Lagrange multipliers and returns the maximum entropy value
                and its gradient.
        """
        def maxent(lambdas: np.ndarray) -> tuple[float,float]:
            """Calculate the maximum entropy value and its gradient.
            
            Args:
                lambdas (np.ndarray):
                    Lagrange multipliers for the optimization.
                    
            Returns:
                tuple[float,float]:
                    float:
                        Maximum entropy value.
                    float:
                        Gradient of the maximum entropy value with respect to the Lagrange
                        multipliers.
            
            Reference:
                S. Bottaro , T. Bengsten and K. Lindorff-Larsen, "Integrating Molecular Simulation
                and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach," pp.
                219-240, Feb. 2020. In: Z. Gáspári, (eds) *Structural Bioinformatics*, *Methods in
                Molecular Biology*, vol. 2112, Humana, New York, NY.
                (https://doi.org/10.1007/978-1-0716-0270-6_15)
            """
            # Compute weights
            arg = - np.sum(lambdas[np.newaxis,:]*self.calculated,axis=1) + np.log(self.w0) - TMAX
            log_z = logsumexp(arg)
            ww = np.exp(arg - log_z)

            # Compute weighted averages of calculated data
            avg = np.sum(ww[:,np.newaxis] * self.calculated,
                         axis=0)

            # Scale experimental uncertainties using theta
            theta_sigma2 = theta * self.experiment[:,1]**2
            
            # Compute gaussian integral
            eps2 = 0.5 * np.sum((lambdas * lambdas) * theta_sigma2)

            # Experimental value
            sum1 = np.dot(lambdas,self.experiment[:,0])

            # Corresponds to function described in Eq. 6 in the reference paper
            fun = log_z + sum1 + eps2

            # Compute gradient
            jac = self.experiment[:,0] + lambdas * theta_sigma2 - avg

            # Divide by theta to avoid numerical problems
            return  fun/theta, jac/theta

        return maxent
    
    def fit(self, theta: int) -> tuple[float | None, float | None, float | None]:
        """Optimize set of weights using BME reweighting.

        Args:
            theta (int):
                value for the theta parameter to use in reweighting.

        Returns:
            tuple[float | None, float | None, float | None]:
                float | None:
                    chisquare value before optimization. None if optimization fails.
                float | None:
                    chisquare value after optimization. None if optimization fails.
                float | None:
                    fraction of effective frames. None if optimization fails.
        """
        # Normalize data if not already done
        if not self.normalized:
            bt.normalize(exp=self.experiment,
                         calc=self.calculated,
                         sample_weights=self.w0)
            self.normalized = True

        # Calculate initial reduced chisquare value
        chi2_before  = bt.calc_chi(exp=self.experiment,
                                   calc=self.calculated,
                                   sample_weights=self.w0)
        # Log initial values
        self.log_message((f'Optimizing {self.experiment.shape[0]} data and '
                           f'{self.calculated.shape[0]} samples. Theta={theta} \n'))
        self.log_message(f'CHI2 before optimization: {chi2_before:8.4f} \n')
        self.log_fd.flush()

        # Setup optimization parameters
        ## Initialize Lagrange multipliers
        lambdas = np.zeros(self.experiment.shape[0],
                           dtype=np.longdouble)
        self.log_message('Lagrange multipliers initialized from zero\n')

        ## Define the maximum entropy function using chosen theta
        maxent_func = self.define_maxent(theta=theta)

        ## Define minimization options
        opt = {'maxiter': 50000,
               'disp': False}
        mini_method = 'L-BFGS-B'

        # Perform optimization
        start_time = time.time()
        result = minimize(fun=maxent_func,
                          x0=lambdas,
                          options=opt,
                          method=mini_method,
                          jac=True)
        self.log_message(f'Execution time: {(time.time() - start_time):.2f} seconds\n')

        # Check if the optimization was successful
        if result.success:
            self.log_message((f'Minimization using {mini_method} successful '
                               f'(iterations:{result.nit})\n'))

            # Calculate optimized set of weights
            arg = -np.sum(result.x[np.newaxis,:] * self.calculated,axis=1) - TMAX
            w_opt = self.w0 * np.exp(arg)
            w_opt /= np.sum(w_opt)

            # Register final set of lagrange multipliers, optimized weights and number of iters
            self.lambdas = np.copy(result.x)
            self.w_opt = np.copy(w_opt)
            self.niter = result.nit

            # Calculate final reduced chisquare value and fraction of effective frames
            chi2_after = bt.calc_chi(self.experiment,
                                     self.calculated,
                                     w_opt)
            phi = np.exp(-bt.srel(self.w0,
                                  w_opt))

            self.log_message(f'CHI2 after optimization: {chi2_after:8.4f} \n')
            self.log_message(f'Fraction of effective frames: {phi:8.4f} \n')
            return chi2_before, chi2_after, phi

        else:
            # Abort
            self.log_message(f'Minimization using {mini_method} failed\n')
            self.log_message(f'Message: {result.message}\n')
            self.niter = -1
            return np.NaN, np.NaN, np.NaN

    def ibme(
        self,
        theta: int,
        ftol: float = 0.01,
        iterations: int = 50,
        ) -> tuple[float, float, float, np.ndarray, np.ndarray]:
        """Iterative BME.

        Args:
            theta (int):
                theta value to use in iBME iterations.
            ftol (float):
                tolerance for minimization procedure. Defaults to 0.01.
            iterations (int):
                number of iBME iterations to perform. Defaults to 50.
            offset (bool):
                whether to offset calculated data at each step when fitting it to experimental
                data. Defaults to True.

        Returns:
            tuple[float, float, float, np.ndarray, np.ndarray]:
                float:
                    Initial chi-square value before optimization.
                float:
                    Final chi-square value after optimization.
                float:
                    Fraction of effective frames after optimization.
                np.ndarray:
                    Initial calculated data before optimization.
                np.ndarray:
                    Rescaled calculated data after optimization.
        """
        # Here we assume that the data and weights are already loaded
        current_weights = self.get_w0()
        w0 = self.get_w0()
        name = self.get_name()
        labels = self.get_labels()
        exp = self.get_experiment()
        calc = self.get_calculated()

        # Setup list to track iBME weights and stats
        self.ibme_weights = []
        self.ibme_stats = []

        # Setup log msg
        log = ''

        # Initialize fit statistics of previous iteration
        fit_stats_old = np.NaN

        for it in range(iterations):
            # Split calculated and experimental data by experimental data type
            groups_calc = bt.split_array_by_groups(array=calc,
                                                   groups_dict=self.exp_types,
                                                   type='column')
            groups_exp = bt.split_array_by_groups(array=exp,
                                                  groups_dict=self.exp_types,
                                                  type='row')

            rescaled_calcs = []
            for exp_type in self.exp_types.keys():
                rescaled_calc, _ = bt.fit_and_scale(exp=groups_exp[exp_type],
                                                    calc=groups_calc[exp_type],
                                                    calc_weights=current_weights,
                                                    fit_type=EXP2FIT[exp_type])
                rescaled_calcs.append(rescaled_calc)

            # Concatenate rescaled calculated data for each experimental data type
            calc = np.hstack(rescaled_calcs)

            # Create SimpleReweight instance, load exp and calc data and fit calc data to exp
            r1 = SimpleReweight(name=f'{name}_ibme_{it}',
                                w0=np.copy(w0))
            r1.load_array(labels,
                          np.copy(exp),
                          np.copy(calc))
            fit_stats = r1.fit(theta=theta)

            # If at first iteration, store initial chi2 and calculated data
            if it == 0:
                chi2_0 = fit_stats[0]
                calc_0 = np.copy(calc)

            # Extract weights after reweighting
            current_weights = np.copy(r1.get_weights())

            # Calculate and log difference in fit statistics
            diff = abs(fit_stats_old - fit_stats[1])
            log += (f'Iteration:{it:3d} chi2: {fit_stats[1]:7.4f} diff: {diff:7.4e}\n')
            
            # Store old fit statistics and update tracking lists
            fit_stats_old = fit_stats[1]
            self.ibme_weights.append(current_weights)
            self.ibme_stats.append(fit_stats)

            # Check end condition
            if diff < ftol:
                line = (f'Iterative procedure converged below tolerance {diff:.2e} '
                        f'after {it} iterations\n')
                print(line,end='')
                log += line
                break
        
        # Save the final calculated data and weights to files
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

        # Write log msgs
        chi2_msg = f'Done. Initial chi2: {chi2_0:8.4f} Final chi2: {fit_stats[1]:8.4f}'
        print(chi2_msg)
        log += chi2_msg + '\n'
        
        outfiles_msg = f'Done. Writing output files {n1} {n2}'        
        print(outfiles_msg)
        log += outfiles_msg + '\n'

        self.log_message(log+'\n')

        # Extract final fraction of effective frames and set of optimized weights
        phi = np.exp(-bt.srel(w0,
                              current_weights))
        self.w_opt = current_weights

        # Return initial chi2, final chi2, fraction eff frames, initial calc data, final calc data
        return chi2_0, fit_stats[1], phi, calc_0, calc

    def get_ibme_weights(self) -> list[np.ndarray]:

        try:
            return self.ibme_weights
        except AttributeError:
            print('# iBME weights not available. Call iBME first')
            sys.exit(1)

    def get_ibme_stats(self) -> list[tuple[float, float, float, np.ndarray, np.ndarray]]:

        try:
            return self.ibme_stats
        except AttributeError:
            print('# iBME stats not available. Call iBME first')
            sys.exit(1)
