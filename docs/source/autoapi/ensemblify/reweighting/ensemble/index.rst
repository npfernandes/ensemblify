ensemblify.reweighting.ensemble
===============================

.. py:module:: ensemblify.reweighting.ensemble

.. autoapi-nested-parse::

   Reweight a conformational ensemble using experimental data.



Functions
---------

.. autoapisummary::

   ensemblify.reweighting.ensemble.reweight_ensemble


Module Contents
---------------

.. py:function:: reweight_ensemble(trajectory: str, topology: str, trajectory_id: str, exp_saxs_data: str, output_dir: str | None = os.getcwd(), thetas: list[int] | None = None, calculated_cmatrix: pandas.DataFrame | str | None = None, calculated_dmatrix: pandas.DataFrame | str | None = None, calculated_ss_frequency: pandas.DataFrame | str | None = None, calculated_metrics_data: pandas.DataFrame | str | None = None, compare_rg: bool | None = True, compare_dmax: bool | None = True, compare_eed: bool | None = True, compare_cmdist: bool | None = None)

   Apply Bayesian Maximum Entropy (BME) reweighting to a conformational ensemble, given
   experimental SAXS data.

   Inside the output directory, a directory named trajectory_id will be created and that is where
   all output files will be stored.
   If calculated metrics data is provided, it will not be recalculated.
   If data for the center mass distance is to be taken from the given calculated metrics data,
   the compare_cmdist mapping must be provided. The identifiers of this mapping will be matched
   to column names of the given DataFrame, if present.

   :param trajectory: Path to .xtc trajectory file where conformational ensemble is stored.
   :type trajectory: :py:class:`str`
   :param topology: Path to .pdb topology file corresponding to any one frame of the ensemble.
   :type topology: :py:class:`str`
   :param trajectory_id: Prefix trajectory identifier to be added to plotted traces and output files.
   :type trajectory_id: :py:class:`str`
   :param exp_saxs_data: Path to .dat file with experimental SAXS data.
   :type exp_saxs_data: :py:class:`str`
   :param output_dir: Path to output directory. Is created if it does not exist. Defaults to current working
                      directory. After output_dir is setup, a directory named trajectory_id is created inside
                      it, where the interactive .html plots and reweighting output files will be stored.
   :type output_dir: :py:class:`str`
   :param thetas: List of values to try as the theta parameter in BME. The ensemble will be reweighted
                  each time using a different theta value. The effect of different theta values can be
                  analyzed in the created effective frames figure.
   :type thetas: :py:class:`list[int]`
   :param calculated_cmatrix: DataFrame with the calculated average contact matrix for the current trajectory or
                              path to this file in .csv format. Defaults to None, and this data is calculated anew.
   :type calculated_cmatrix: :py:class:`pd.DataFrame | str, optional`
   :param calculated_dmatrix: DataFrame with the calculated average distance matrix for the current trajectory or
                              path to this file in .csv format. Defaults to None, and this data is calculated anew.
   :type calculated_dmatrix: :py:class:`pd.DataFrame | str, optional`
   :param calculated_ss_frequency: DataFrame with the calculated secondary structure assignment frequency matrix for the
                                   current trajectory or path to this file in .csv format. Defaults to None, and this data
                                   is calculated anew.
   :type calculated_ss_frequency: :py:class:`pd.DataFrame | str, optional`
   :param calculated_metrics_data: DataFrame with calculated structural metrics (columns) for each frame of the trajectory
                                   (rows) or path to this DataFrame in .csv format. Defaults to None, and this data is
                                   calculated anew.
   :type calculated_metrics_data: :py:class:`pd.DataFrame | str, optional`
   :param compare_rg: Whether to calculate/consider the radius of gyration when comparing structural metrics
                      between uniform and reweighted conformational ensembles. Defaults to True.
   :type compare_rg: :py:class:`bool, optional`
   :param compare_dmax: Whether to calculate/consider the maximum distance between any two alpha carbons when
                        comparing structural metrics between uniform and reweighted conformational ensembles.
                        Defaults to True.
   :type compare_dmax: :py:class:`bool, optional`
   :param compare_eed: Whether to calculate/consider the distance from the N to C terminal when comparing
                       structural metrics between uniform and reweighted conformational ensembles. Defaults
                       to True.
   :type compare_eed: :py:class:`bool, optional`
   :param compare_cmdist: Mapping of identifiers to tuples with two selection strings for creating MDAnalysis
                          AtomGroups, whose center mass distance will be calculated. For example:

                              {'inter_domain' : ('resid 1:30', 'resid 110:140')}

                          If None, no center mass distances are calculated or compared.
                          See https://userguide.mdanalysis.org/stable/selections.html for more information about
                          MDAnalysis selections.
                          Defaults to None.
   :type compare_cmdist: :py:class:`bool, optional`


