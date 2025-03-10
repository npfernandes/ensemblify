ensemblify.reweighting.data
===========================

.. py:module:: ensemblify.reweighting.data

.. autoapi-nested-parse::

   Auxiliary functions for reweighting ensembles.



Functions
---------

.. autoapisummary::

   ensemblify.reweighting.data.process_exp_data
   ensemblify.reweighting.data.correct_exp_error
   ensemblify.reweighting.data.ibme
   ensemblify.reweighting.data.bme_ensemble_reweighting
   ensemblify.reweighting.data.average_saxs_profiles
   ensemblify.reweighting.data.attempt_read_calculated_data
   ensemblify.reweighting.data.attempt_read_reweighting_data


Module Contents
---------------

.. py:function:: process_exp_data(experimental_data_path: str) -> str

   Check formatting and units in input experimental data file.

   If values for q are in Ångstrom, they are converted to nanometer.
   Any q-values above 5nm^(-1) are removed, as SAXS calculations are not reliable in that
   range.

   :param experimental_data_path: Path to experimental SAXS data file.
   :type experimental_data_path: :py:class:`str`

   :returns:     Path to experimental SAXS data file with any applied changes.
   :rtype: str

   Adapted from:
       https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb


.. py:function:: correct_exp_error(experimental_data_path: str) -> str

   Correct experimental error of input experimental data file using BIFT.

   Bayesian Indirect Fourier Transformation (BIFT) can identify whether the
   experimental error in small-angle scattering data is over- or
   underestimated. The error values are then scaled accordingly.

   Reference:
       Larsen, A.H. and Pedersen, M.C. (2021), Experimental noise in small-angle scattering
       can be assessed using the Bayesian indirect Fourier transformation. J. Appl. Cryst.,
       54: 1281-1289. https://doi.org/10.1107/S1600576721006877

   :param experimental_data_path: Path to experimental SAXS data file.
   :type experimental_data_path: :py:class:`str`

   :returns:     Path to experimental SAXS data file with corrected errors.
   :rtype: str

   Adapted from:
       https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb


.. py:function:: ibme(theta: int, exp_file: str, calc_file: str, output_dir: str) -> tuple[int, tuple[float, float, float], numpy.ndarray]

   Apply the Iterative Bayesian Maximum Entropy (BME) algorithm on calculated SAXS data,
   given a value for the theta parameter.

   The used algorithm is explained in:
       https://github.com/KULL-Centre/BME/blob/main/notebook/example_04.ipynb

   Reference:
       Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
       Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
       doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

   :param theta: Value for the theta parameter to be used in BME algorithm.
   :type theta: :py:class:`int`
   :param exp_file: Path to .dat file with experimental SAXS curve.
   :type exp_file: :py:class:`str`
   :param calc_file: Path to .dat file with SAXS curve calculated from an ensemble.
   :type calc_file: :py:class:`str`
   :param output_dir: Path to directory where all the files resulting from the reweighting procedure will be
                      stored.
   :type output_dir: :py:class:`str`

   :returns:

                 theta (int):
                     Value for the theta parameter used in BME algorithm (same as input).
                 stats (tuple[float,float,float]):
                     chi2_before (float):
                         The value for the chisquare of fitting the ensemble with uniform
                         weights to the experimental data.
                     chi2_after (float):
                         The value for the chisquare of fitting the reweighted ensemble to
                         the experimental data.
                     phi (float):
                         The fraction of effective frames being used in the reweighted ensemble.
                 weights (np.ndarray):
                     An array containing the new weights of the ensemble, one for each frame.
   :rtype: tuple[int,tuple[float,float,float],np.ndarray]

   Adapted from:
       https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb


.. py:function:: bme_ensemble_reweighting(exp_saxs_file: str, calc_saxs_file: str, thetas: list[int], output_dir: str) -> tuple[numpy.ndarray, numpy.ndarray]

   Perform Bayesian Maximum Entropy (BME) reweighting on calculated SAXS data based on
   experimental SAXS data of the same protein.

   Applies the iterative BME algorithm, explained in:
       https://github.com/KULL-Centre/BME/blob/main/notebook/example_04.ipynb

   The algorithm is applied using different theta values and the results for each value are stored.

   Reference:
       Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
       Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
       doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

   :param exp_saxs_file: Path to .dat file with experimental SAXS data.
   :type exp_saxs_file: :py:class:`str`
   :param calc_saxs_file: Path to .dat file with SAXS data calculated from a conformational ensemble.
   :type calc_saxs_file: :py:class:`str`
   :param thetas: Values of theta to try when applying iBME.
   :type thetas: :py:class:`list[int]`
   :param output_dir: Path to directory where output files from reweighting protocol will be stored.
   :type output_dir: :py:class:`str`

   :returns:

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
   :rtype: tuple[np.ndarray,np.ndarray]


.. py:function:: average_saxs_profiles(exp_saxs_file: str, calc_saxs_file: str, rw_calc_saxs_file: str, weights: numpy.ndarray) -> tuple[float, float]

   Average the SAXS intensities for uniform and reweighted calculated SAXS data.
   The uniform data is then scaled and offset by linear regression fitting to experimental data.

   :param exp_saxs_file: Path to .dat file with experimental SAXS data.
   :type exp_saxs_file: :py:class:`str`
   :param calc_saxs_file: Path to .dat file with SAXS data calculated from a conformational ensemble.
   :type calc_saxs_file: :py:class:`str`
   :param rw_calc_saxs_file: Path to .dat file with SAXS data calculated from a conformational ensemble considering
                             the weights (from iBME) for each frame.
   :type rw_calc_saxs_file: :py:class:`str`
   :param weights: Array resulting from iBME with weights for each data point. Defaults to uniform weights.
   :type weights: :py:class:`np.ndarray`

   :returns:

                 i_prior (float):
                     an array of SAXS intensities averaged over all the frames of a SAXS data file
                     calculated from a conformational ensemble with uniform weights.
                 i_post (float):
                     an array of SAXS intensities averaged over all the frames of a SAXS data file
                     calculated from a conformational ensemble with the provided set of weights.
   :rtype: tuple[float,float]


.. py:function:: attempt_read_calculated_data(data: pandas.DataFrame | str | None, data_msg_tag: str, calc_fn: collections.abc.Callable, *args, **kwargs) -> pandas.DataFrame

   Attempt to read data from file, else calculate it using provided function.

   If data is given directly as a DataFrame, it is simply returned. Otherwise, it
   is either read from file or calculated using the provided function and arguments.

   :param data: A DataFrame with the desired data, the path to the data in .csv format or None.
   :type data: :py:class:`pd.DataFrame | str | None`
   :param data_msg_tag: String identifier for which data we are working with so prints to console are
                        correct.
   :type data_msg_tag: :py:class:`str`
   :param calc_fn: An object with a __call__ method, e.g. a function to be used in calculating the
                   data if it is not provided.
   :type calc_fn: :py:class:`Callable`

   :returns:     Desired data in DataFrame format.
   :rtype: pd.DataFrame


.. py:function:: attempt_read_reweighting_data(reweighting_output_directory: str, trajectory_id: str) -> tuple[str | None, str | None, numpy.ndarray | None, numpy.ndarray | None, numpy.ndarray | None]

   Attempt to read reweighting data from output directory, returning None if not found.

   :param reweighting_output_directory: Directory where data should be searched.
   :type reweighting_output_directory: :py:class:`str`
   :param trajectory_id: Prefix for filenames to look for in directory.
   :type trajectory_id: :py:class:`str`

   :returns:

                 exp_saxs_file (str | None):
                     The corresponding data (if found) or None (if not found).
                 calc_saxs_file (str | None):
                     The corresponding data (if found) or None (if not found).
                 thetas_array (np.ndarray | None):
                     The corresponding data (if found) or None (if not found).
                 stats (np.ndarray | None):
                     The corresponding data (if found) or None (if not found).
                 weights (np.ndarray | None):
                     The corresponding data (if found) or None (if not found).
   :rtype: tuple[str|None, str|None, np.ndarray|None, np.ndarray|None, np.ndarray|None]


