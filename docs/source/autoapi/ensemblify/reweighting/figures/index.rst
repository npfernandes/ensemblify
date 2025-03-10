ensemblify.reweighting.figures
==============================

.. py:module:: ensemblify.reweighting.figures

.. autoapi-nested-parse::

   Auxiliary functions for creating reweighting figures.



Functions
---------

.. autoapisummary::

   ensemblify.reweighting.figures.create_effective_frames_fit_fig
   ensemblify.reweighting.figures.create_reweighting_fits_fig
   ensemblify.reweighting.figures.create_reweighting_metrics_fig


Module Contents
---------------

.. py:function:: create_effective_frames_fit_fig(stats: numpy.ndarray, thetas: numpy.ndarray, choices: int | list[int] | None = None, title_text: str | None = None, colors: list[str] = None) -> plotly.graph_objects.Figure

   Create a Figure plotting the fraction of effective frames vs the chisquare value, resulting
   from applying BME using different theta values.

   The fraction of effective frames of an ensemble after reweighting is plotted agaisnt the
   chisquare value of the fitting of the data calculated from the reweighted ensemble to the
   experimental data.
   Each data point results from the application of the Bayesian Maximum Entropy (BME) algorithm to
   the calculated+experimental data using different values for the theta parameter.

   :param stats: An array where each row corresponds to a different theta value with columns
                 (chi2_before,chi2_after,phi) where:
                     chi2_before:
                         The value for the chisquare of fitting the ensemble with uniform
                         weights to the experimental data.
                     chi2_after:
                         The value for the chisquare of fitting the reweighted ensemble to
                         the experimental data.
                     phi:
                         The fraction of effective frames being used in the reweighted ensemble.
   :type stats: :py:class:`np.ndarray`
   :param thetas: Array of values for the theta parameter used when applying BME algorithm.
   :type thetas: :py:class:`np.ndarray`
   :param choices: Theta value(s) chosen for reweighting ensemble, corresponding data points will be
                   highlighted in the created Figure. Defaults to None.
   :type choices: :py:class:`int | list[int] | None`
   :param title_text: Title for the created Figure. Defaults to None.
   :type title_text: :py:class:`str | None`
   :param colors: Hexcodes for the colors to use for highlighting theta values. Defaults to ['#E69F00',
                  '#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].
   :type colors: :py:class:`list[str]`

   :returns:     The created plot, optionally with data points corresponding to highlighted theta
                 values in different colors.
   :rtype: go.Figure


.. py:function:: create_reweighting_fits_fig(q: numpy.ndarray, i_exp: numpy.ndarray, err: numpy.ndarray, i_prior: numpy.ndarray, i_posts: numpy.ndarray | list[numpy.ndarray], title_text: str | None = None, colors: list[str] = None) -> plotly.graph_objects.Figure

   Create a multiplot Figure showcasing the differences between uniform and reweighted
   calculated SAXS data, when fit to experimental data.

   :param q: An array with momentum transfer values, common to all SAXS curves being deal with.
   :type q: :py:class:`np.ndarray`
   :param i_exp: An array with experimentally measured SAXS intensities.
   :type i_exp: :py:class:`np.ndarray`
   :param err: An array with the experimental error of the provided experimental SAXS intensities.
   :type err: :py:class:`np.ndarray`
   :param i_prior: An array of SAXS intensities averaged over all the frames of a SAXS data file
                   calculated from a conformational ensemble with uniform weights.
   :type i_prior: :py:class:`np.ndarray`
   :param i_posts: An array or list of arrays of SAXS intensities averaged over all the frames of a SAXS
                   data file calculated from a conformational ensemble with the provided set of weights.
   :type i_posts: :py:class:`np.ndarray | list[np.ndarray]`
   :param title_text: A title for the created multiplot Figure. Defaults to None.
   :type title_text: :py:class:`str, optional`
   :param colors: Color to attribute to the plotted prior and posterior traces, in order of input.
                  Defaults to ['#E69F00', '#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].
   :type colors: :py:class:`list[str], optional`

   :returns:

                 A multiplot Figure containing four plots:
                     - the fitting of i_prior and i_post(s) to the experimental SAXS data i_exp.
                     - the previous plot in log scale.
                     - Kraty plot for i_prior and i_post fitted to experimental data.
                     - residuals between i_prior/i_post(s) and i_exp.
   :rtype: go.Figure


.. py:function:: create_reweighting_metrics_fig(metrics: pandas.DataFrame, weight_sets: numpy.ndarray | list[numpy.ndarray], title_text: str | None = None, colors: list[str] = None) -> plotly.graph_objects.Figure

   Create a Figure with probability distribution plots for calculated structural metrics, using
   uniform or unequal weights.

   :param metrics: A DataFrame with the calculated structural metrics, one row per frame in the
                   conformational ensemble.
   :type metrics: :py:class:`pd.DataFrame`
   :param weight_sets: An array or list of arrays containing the weights for calculating the probability
                       distributions of each structural metric, for each set of weights.
   :type weight_sets: :py:class:`np.ndarray | list[np.ndarray]`
   :param title_text: Title for the created Figure. Defaults to None.
   :type title_text: :py:class:`str | None`
   :param colors: Hexcodes for colors to use for the traces relative to each i_post, in order of input.
                  Defaults to ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].
   :type colors: :py:class:`list[str]`

   :returns:     A Figure plotting the structural metrics distributions for uniformly and unequally
                 weighted conformational ensembles.
   :rtype: go.Figure


