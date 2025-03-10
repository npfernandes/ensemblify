ensemblify.analysis.trajectory
==============================

.. py:module:: ensemblify.analysis.trajectory

.. autoapi-nested-parse::

   Analyze a trajectory,topology pair, outputting a graphical dashboard.



Functions
---------

.. autoapisummary::

   ensemblify.analysis.trajectory.analyze_trajectory


Module Contents
---------------

.. py:function:: analyze_trajectory(trajectories: list[str] | str, topologies: list[str] | str, trajectory_ids: list[str] | str, output_directory: str | None = os.getcwd(), ramachandran_data: bool | None = True, distancematrices: bool | None = True, contactmatrices: bool | None = True, ssfrequencies: bool | None = True, rg: bool | None = True, dmax: bool | None = True, eed: bool | None = True, cm_dist: dict[str, tuple[str, str]] | None = None, color_palette: list[str] | None = None) -> dict[str, pandas.DataFrame]

   Calculate structural data and create interactive figures for given trajectory and
   topology files.

   :param trajectories: List of paths to .xtc trajectory files or string with the path to a single .xtc
                        trajectory file.
   :type trajectories: :py:class:`list[str] | str`
   :param topologies: List of paths to .pdb topology files or string with the path to a single .pdb
                      topology file.
   :type topologies: :py:class:`list[str] | str`
   :param trajectory_ids: List of prefix trajectory identifiers to distinguish between calculated data
                          files or string with a single prefix trajectory identifier.
   :type trajectory_ids: :py:class:`list[str] | str`
   :param output_directory: Path to directory where calculated data and created figures will be stored.
                            If it does not exist, it is created. Defaults to current working directory.
   :type output_directory: :py:class:`str`
   :param ramachandran_data: Whether to calculate a dihedral angles matrix for each trajectory,topology
                             file pair.
   :type ramachandran_data: :py:class:`bool`
   :param distancematrices: Whether to calculate an alpha carbon distance matrix for each trajectory,topology
                            file pair and create the corresponding distance matrix interactive figure.
   :type distancematrices: :py:class:`bool`
   :param contactmatrices: Whether to calculate a contact frequency matrix for each trajectory,topology
                           file pair and create the corresponding contact map interactive figure.
   :type contactmatrices: :py:class:`bool`
   :param ssfrequencies: Whether to calculate a secondary structure assignment frequency matrix for each
                         trajectory, topology file pair and create the corresponding secondary structure
                         frequency interactive figure.
   :type ssfrequencies: :py:class:`bool`
   :param rg: Whether to calculate and plot the radius of gyration for each trajectory,topology
              file pair.
   :type rg: :py:class:`bool`
   :param dmax: Whether to calculate and plot the maximum distance between any two alpha carbons
                for each trajectory,topology file pair.
   :type dmax: :py:class:`bool`
   :param eed: Whether to calculate and plot the distance between the N and C terminal for each
               trajectory, topology file pair.
   :type eed: :py:class:`bool`
   :param cm_dist: Mapping of identifiers to tuples with two selection strings for creating MDAnalysis
                   AtomGroups, whose center mass distance will be calculated and plotted. For example:

                       {'inter_domain' : ('resid 1:30', 'resid 110:140')}

                   If None, no center mass distances are calculated.
                   See https://userguide.mdanalysis.org/stable/selections.html for more information about
                   MDAnalysis selections.
   :type cm_dist: :py:class:`dict[str,tuple[str,str]]`
   :param color_palette: List of color hexcodes, to associate one with each trajectory in the created Structural
                         Metrics interactive dashboard.
   :type color_palette: :py:class:`list[str]`

   :returns:     Mapping of data identifiers to DataFrames containing the calculated analysis data for
                 each frame of each given trajectory. For convenience, this is returned as a variable
                 as well as saved to output directory.
   :rtype: dict[str,pd.DataFrame]


