ensemblify.analysis.figures
===========================

.. py:module:: ensemblify.analysis.figures

.. autoapi-nested-parse::

   Auxiliary functions for creating analysis figures.



Functions
---------

.. autoapisummary::

   ensemblify.analysis.figures.create_ramachandran_figure
   ensemblify.analysis.figures.create_contact_map_fig
   ensemblify.analysis.figures.create_distance_matrix_fig
   ensemblify.analysis.figures.create_ss_frequency_figure
   ensemblify.analysis.figures.create_metrics_traces
   ensemblify.analysis.figures.create_metrics_fig
   ensemblify.analysis.figures.create_single_metrics_fig_directly
   ensemblify.analysis.figures.create_analysis_figures


Module Contents
---------------

.. py:function:: create_ramachandran_figure(dihedrals_matrix: pandas.DataFrame | str, trajectory_id: str | None = None, output_path: str | None = None) -> plotly.graph_objects.Figure

   Create a ramachandran plot Figure from a calculated dihedral angles matrix.

   :param dihedrals_matrix: Calculated dihedral angles matrix DataFrame or path to calculated matrix in .csv format.
                            If difference is True, this should be the difference dihedral angles matrix between the
                            uniformly weighted and the reweighted dihedral angles matrix.
   :type dihedrals_matrix: :py:class:`pd.DataFrame | str`
   :param trajectory_id: Used on Figure title and prefix for saved ramachandran plot filename. Defaults to None.
   :type trajectory_id: :py:class:`str, optional`
   :param output_path: Path to output .html file or output directory where created Figure will be stored.
                       If directory, written file is named 'ramachandran_plot.html', optionally with
                       trajectory_id prefix. Defaults to None.
   :type output_path: :py:class:`str, optional`

   :returns:     Ploty Figure object displaying a ramachandran plot.
   :rtype: go.Figures


.. py:function:: create_contact_map_fig(contact_matrix: pandas.DataFrame | str, topology: str, trajectory_id: str | None = None, output_path: str | None = None, reweighted: bool | None = False, difference: bool | None = False) -> plotly.graph_objects.Figure

   Create a contact map Figure from a calculated contact matrix.

   The topology provides information about number of chains, their chain letters and
   residue numbers.

   :param contact_matrix: Calculated contact matrix DataFrame or path to calculated matrix in .csv format.
                          If difference is True, this should be the difference contact matrix between the
                          uniformly weighted and the reweighted contact matrix.
   :type contact_matrix: :py:class:`pd.DataFrame | str`
   :param topology: Path to topology .pdb file.
   :type topology: :py:class:`str`
   :param trajectory_id: Used on Figure title and prefix for saved contact map filename. Defaults to None.
   :type trajectory_id: :py:class:`str, optional`
   :param output_path: Path to output .html file or output directory where created Figure will be stored.
                       If directory, written file is named 'contact_map.html', optionally with
                       trajectory_id prefix. Defaults to None.
   :type output_path: :py:class:`str, optional`
   :param reweighted: Boolean stating whether we are creating a reweighted contact map figure or a default
                      one. Defaults to False.
   :type reweighted: :py:class:`bool, optional`
   :param difference: Boolean stating whether we are creating a difference contact map figure or a default
                      one. Defaults to False.
   :type difference: :py:class:`bool, optional`

   :returns:     Ploty Figure object displaying a contact map.
   :rtype: go.Figure


.. py:function:: create_distance_matrix_fig(distance_matrix: pandas.DataFrame | str, topology: str, trajectory_id: str | None = None, output_path: str | None = None, max_colorbar: int | None = None, min_colorbar: int | None = None, reweighted: bool | None = False, difference: bool | None = False) -> plotly.graph_objects.Figure

   Create a distance matrix Figure from a calculated distance matrix.

   The topology provides information about number of chains, their chain letters and
   residue numbers.

   :param distance_matrix: Calculated distance matrix DataFrame or path to calculated matrix in .csv format.
                           If difference is True, this should be the difference distance matrix between the
                           uniformly weighted and the reweighted distance matrix.
   :type distance_matrix: :py:class:`pd.DataFrame | str`
   :param topology: Path to topology .pdb file.
   :type topology: :py:class:`str`
   :param trajectory_id: Used on Figure title and prefix for saved distance matrix filename. Defaults to None.
   :type trajectory_id: :py:class:`str, optional`
   :param output_path: Path to output .html file or output directory where created Figure will be stored.
                       If directory, written file is named 'distance_matrix.html', optionally with
                       trajectory_id prefix. Defaults to None.
   :type output_path: :py:class:`str, optional`
   :param max_colorbar: Maximum limit for the distance colorbar. Defaults to None, in which case it is
                        derived from the data.
   :type max_colorbar: :py:class:`int, optional`
   :param min_colorbar: Minimum limit for the distance colorbar. Defaults to None, in which case it is
                        derived from the data.
   :type min_colorbar: :py:class:`int, optional`
   :param reweighted: Boolean stating whether we are creating a reweighted distance matrix figure or a
                      default one. Defaults to False.
   :type reweighted: :py:class:`bool, optional`
   :param difference: Boolean stating whether we are creating a difference distance matrix figure or a
                      default one. Defaults to False.
   :type difference: :py:class:`bool, optional`

   :returns:     Ploty Figure object displaying a distance matrix.
   :rtype: go.Figure


.. py:function:: create_ss_frequency_figure(ss_frequency: pandas.DataFrame | str, topology: str, trajectory_id: str | None = None, output_path: str | None = None, reweighted: bool = False, difference: bool = False) -> plotly.graph_objects.Figure

   Create a secondary structure frequency Figure from a secondary structure assignment
   frequency matrix.

   The topology provides information about number of chains, their chain letters and
   residue numbers.

   :param ss_frequency: Calculated secondary structure assignment frequency matrix DataFrame or path to
                        calculated matrix in .csv format.
   :type ss_frequency: :py:class:`pd.DataFrame | str`
   :param topology: Path to topology .pdb file.
   :type topology: :py:class:`str`
   :param trajectory_id: Used on Figure title and prefix for saved ss_frequency filename. Defaults to None.
   :type trajectory_id: :py:class:`str, optional`
   :param output_path: Path to output .html file or output directory where created Figure will be stored.
                       If directory, written file is named 'ss_frequency.html', optionally with
                       trajectory_id prefix. Defaults to None.
   :type output_path: :py:class:`str, optional`
   :param reweighted: Boolean stating whether we are creating a reweighted secondary structure frequency
                      figure or a default one. Defaults to False.
   :type reweighted: :py:class:`bool, optional`
   :param difference: Boolean stating whether we are creating a difference secondary structure frequency
                      figure or a default one. Defaults to False.
   :type difference: :py:class:`bool, optional`

   :returns:     Stacked line plot with the secondary structure frequencies of each secondary
                 structure type for each residue in the structure.
   :rtype: go.Figure


.. py:function:: create_metrics_traces(metrics: pandas.DataFrame | str, trajectory_id: str, color: str = '#1f77b4') -> tuple[list[plotly.graph_objects.Box], list[plotly.graph_objects.Histogram], list[plotly.graph_objects.Scatter], list[float], list[float]]

   Create Ploty Box, Histogram and Scatter (KDE) traces to be used in the creation of
   a Structural Metrics Figure.

   :param metrics: DataFrame where columns are the desired structural metrics and rows are the frames
                   of the trajectory or path to that DataFrame in .csv format.
   :type metrics: :py:class:`pd.DataFrame | str`
   :param trajectory_id: prefix identifier for trace names.
   :type trajectory_id: :py:class:`str`
   :param color: hex code for the color the created traces will be. Defaults to '#636EFA', or light
                 blue.
   :type color: :py:class:`str, optional`

   :returns:

                 box_traces (list[go.Box]):
                     A list of the boxplot traces, one for each structural metric.
                 hist_traces (list[go.Histogram]):
                     A list of the histogram traces, one for each structural metric.
                 scatter_traces (list[go.Scatter]):
                     A list of the scatter Kernel Density Estimate (KDE) traces, one
                     for each structural metric.
                 avg_values (list[float]):
                     A list of the values of the mean for each metric.
                 avg_stderr (list[float]):
                     A list of the values of the standard error of the mean for each metric.
   :rtype: tuple[list[go.Box], list[go.Histogram], list[go.Scatter], list[float], list[float]]


.. py:function:: create_metrics_fig(trajectory_ids: list[str], total_box_traces: dict[str, list[plotly.graph_objects.Box]], total_hist_traces: dict[str, list[plotly.graph_objects.Histogram]], total_scatter_traces: dict[str, list[plotly.graph_objects.Scatter]], total_avg_values: dict[str, list[float]], total_avg_stderr_values: dict[str, list[float]], output_path: str | None = None) -> plotly.graph_objects.Figure

   Create a Structural Metrics Figure from previously created Box, Histogram and Scatter traces.

   :param trajectory_ids: List of prefix identifiers that must match the prefix identifiers used for naming the
                          created traces.
   :type trajectory_ids: :py:class:`list[str]`
   :param total_box_traces: Mapping of trajectory_ids to a list of created Box traces.
   :type total_box_traces: :py:class:`dict[str,list[go.Box]]`
   :param total_hist_traces: Mapping of trajectory_ids to a list of created Histogram traces.
   :type total_hist_traces: :py:class:`dict[str,list[go.Histogram]]`
   :param total_scatter_traces: Mapping of trajectory_ids to a list of created Scatter traces.
   :type total_scatter_traces: :py:class:`dict[str,list[go.Scatter]]`
   :param total_avg_values: Mapping of trajectory_ids to a list of mean values.
   :type total_avg_values: :py:class:`dict[str,list[float]]`
   :param output_path: Path to output .html file or output directory where the created Figure will be stored.
                       If directory, written file is named 'structural_metrics.html'. Defaults to None.
   :type output_path: :py:class:`str, optional`

   :returns:     Structural Metrics Figure that allows for comparison between all the created traces.
   :rtype: go.Figure


.. py:function:: create_single_metrics_fig_directly(metrics: pandas.DataFrame | str, trajectory_id: str | None = None) -> plotly.graph_objects.Figure

   Create a Structural Metrics Figure for a single trajectory directly from its trajectory
   id and calculated metrics data.

   :param metrics: DataFrame where columns are the desired structural metrics and rows are the frames
                   of the trajectory or path to that DataFrame in .csv format.
   :type metrics: :py:class:`pd.DataFrame | str`
   :param trajectory_id: Prefix identifier for trace names.
   :type trajectory_id: :py:class:`str, optional`

   :returns:     Structural Metrics Figure with the provided calculated metrics data.
   :rtype: go.Figure


.. py:function:: create_analysis_figures(analysis_data: dict[str, list[pandas.DataFrame]] | None, topologies: list[str], trajectory_ids: list[str], output_directory: str | None = os.getcwd(), color_palette: list[str] | None = None) -> dict[str, list[plotly.graph_objects.Figure]]

   Create interactive figures given analysis data for one or more pairs of trajectory,topology
     files.

   :param analysis_data: Mapping of data identifiers to lists of DataFrames with the calculated analysis data,
                         one element for each given trajectory,topology,trajectory_id trio.
   :type analysis_data: :py:class:`dict[str,list[pd.DataFrame]], optional`
   :param topologies: List of paths to .pdb topology files.
   :type topologies: :py:class:`list[str]`
   :param trajectory_ids: Prefix trajectory identifiers to distinguish between calculated data files.
   :type trajectory_ids: :py:class:`list[str]`
   :param output_directory: Path to directory where created Figures will be stored. Defaults to current
                            working directory.
   :type output_directory: :py:class:`str, optional`
   :param color_palette: List of color hexcodes, to associate one with each trajectory when creating the
                         Structural Metrics interactive figure.
   :type color_palette: :py:class:`list[str], optional`

   :returns:     Mapping of figure identifiers to lists of the created Figures, one for each trajectory
                 outlined in the given analysis data. For example:

                 data = {'ContactMaps' : [ContactMap1,ContactMap2,ContactMap3],
                 'DistanceMatrices' : [DistanceMatrix1,DistanceMatrix2,DistanceMatrix3],
                 'SecondaryStructureFrequencies' : [SSFrequency1,SSFrequency2,SSFrequency3],
                 'StructuralMetrics' : [StructuralMetrics1,StructuralMetrics2,StructuralMetrics3]}
   :rtype: dict[str,list[go.Figure]]


