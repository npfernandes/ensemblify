ensemblify.analysis.data
========================

.. py:module:: ensemblify.analysis.data

.. autoapi-nested-parse::

   Auxiliary functions for calculating structural properties data.



Functions
---------

.. autoapisummary::

   ensemblify.analysis.data.calculate_ramachandran_data
   ensemblify.analysis.data.calculate_contact_matrix_frame
   ensemblify.analysis.data.calculate_contact_matrix
   ensemblify.analysis.data.calculate_distance_matrix_frame
   ensemblify.analysis.data.calculate_distance_matrix
   ensemblify.analysis.data.calculate_ss_assignment
   ensemblify.analysis.data.calculate_ss_frequency
   ensemblify.analysis.data.calc_rg
   ensemblify.analysis.data.calc_eed
   ensemblify.analysis.data.calc_dmax
   ensemblify.analysis.data.calc_cm_dist
   ensemblify.analysis.data.calculate_metrics_data
   ensemblify.analysis.data.calculate_analysis_data


Module Contents
---------------

.. py:function:: calculate_ramachandran_data(trajectory: str, topology: str, output_path: str | None = os.getcwd()) -> pandas.DataFrame

   Calculate a dihedral angles matrix from trajectory and topology files.

   Phi and Psi dihedral angle values are calculated for each residue in each trajectory frame.
   Optionally saves the matrix to output directory in .csv format, defaulting to current working
   directory.

   :param trajectory: Path to .xtc trajectory file.
   :type trajectory: :py:class:`str`
   :param topology: Path to .pdb topology file.
   :type topology: :py:class:`str`
   :param output_path: Path to output .csv file or output directory. If directory, written file is named
                       'ramachandran_data.csv'. Defaults to current working directory.
   :type output_path: :py:class:`str, optional`

   :returns:     DataFrame with Phi and Psi values of each residue for each frame of the trajectory.
   :rtype: pd.DataFrame


.. py:function:: calculate_contact_matrix_frame(u: MDAnalysis.Universe, frame_idx: int, frame_weight: float) -> numpy.ndarray

   Calculates a contact matrix for a frame of a trajectory.

   :param u: `MDAnalysis.Universe` object containing the trajectory being analyzed.
   :type u: :py:class:`mda.Universe`
   :param frame_idx: Number of the frame to be analyzed.
   :type frame_idx: :py:class:`int`
   :param frame_weight: Contacts found in this frame will be assigned this value in the
                        resulting matrix instead of the default value of 1. In a uniformly
                        weighted matrix, this value will be of 1 / number of trajectory frames.
   :type frame_weight: :py:class:`float`

   :returns:     Contact matrix for the current frame.
   :rtype: np.ndarray


.. py:function:: calculate_contact_matrix(trajectory: str, topology: str, weights: numpy.ndarray | None = None, output_path: str | None = os.getcwd()) -> pandas.DataFrame

   Calculate a contact frequency matrix from a trajectory and topology files.

   The contact frequency of a residue pair is calculated from the number of times they are in
   contact over all frames in the trajectory.
   Optionally saves the matrix to output directory in .csv format.
   Uses multiprocessing whenever possible.

   :param trajectory: Path to .xtc trajectory file.
   :type trajectory: :py:class:`str`
   :param topology: Path to .pdb topology file.
   :type topology: :py:class:`str`
   :param weights: Array of weights to be used when calculating the contact matrix. If None, uniform
                   weights are used.
   :type weights: :py:class:`np.ndarray, optional`
   :param output_path: Path to output .csv file or output directory. If directory, written file is named
                       'contact_matrix.csv'. Defaults to current working directory.
   :type output_path: :py:class:`str, optional`

   :returns:     DataFrame with the frequency of each residue contact in the trajectory.
   :rtype: pd.DataFrame


.. py:function:: calculate_distance_matrix_frame(u: MDAnalysis.Universe, frame_idx: int, frame_weight: float) -> numpy.ndarray

   Calculates a distance matrix for the alpha carbons of a trajectory frame.

   :param u: `MDAnalysis.Universe` object containing the trajectory being analyzed.
   :type u: :py:class:`mda.Universe`
   :param frame_idx: Number of the frame to be analyzed.
   :type frame_idx: :py:class:`int`
   :param frame_weight: Distances calculated for this frame will be multiplied by this value
                        in the resulting frame matrix. In a uniformly weighted matrix, calculated
                        distances will be multiplied by 1 / number of trajectory frames.
   :type frame_weight: :py:class:`float`

   :returns:     Distance matrix for the current frame.
   :rtype: np.ndarray


.. py:function:: calculate_distance_matrix(trajectory: str, topology: str, weights: numpy.ndarray | None = None, output_path: str | None = os.getcwd()) -> pandas.DataFrame

   Calculate an alpha carbon average distance matrix from a trajectory and topology files.

   The distances between different pairs of alpha carbons pair is calculated for each trajectory
   frame and the values are then averaged to create the final distance matrix.

   Optionally save the matrix to output directory in .csv format.
   Uses multiprocessing whenever possible.

   :param trajectory: Path to .xtc trajectory file.
   :type trajectory: :py:class:`str`
   :param topology: Path to .pdb topology file.
   :type topology: :py:class:`str`
   :param weights: Array of weights to be used when calculating the distance matrix. If None, uniform
                   weights are used.
   :type weights: :py:class:`np.ndarray, optional`
   :param output_path: Path to output .csv file or output directory. If directory, written file is named
                       'distance_matrix.csv'. Defaults to current working directory.
   :type output_path: :py:class:`str, optional`

   :returns:     DataFrame with the average distance between each pair of alpha carbons in the
                 trajectory.
   :rtype: pd.DataFrame


.. py:function:: calculate_ss_assignment(trajectory: str, topology: str, output_path: str | None = None) -> pandas.DataFrame

   Calculate a secondary structure assignment matrix from a trajectory and topology files.

   For each residue in each frame of the trajectory, calculate it's secondary structure
   assignment using DSSP. The simplified DSSP codes used here are:

       'H' : Helix. Either of the 'H', 'G', or 'I' codes.

       'E' : Strand. Either of the 'E', or 'B' codes.

       'C' : Coil. Either of the 'T', 'S' or ' ' codes.

   Optionally save the resulting matrix to output directory in .csv format.

   :param trajectory: Path to .xtc trajectory file.
   :type trajectory: :py:class:`str`
   :param topology: Path to .pdb topology file.
   :type topology: :py:class:`str`
   :param output_path: Path to output .csv file or output directory. If directory, written file is named
                       'ss_assignment.csv'. Defaults to None, and no file is written.
   :type output_path: :py:class:`str, optional`

   :returns:     DataFrame holding the secondary structure assignment matrix.
   :rtype: pd.DataFrame


.. py:function:: calculate_ss_frequency(trajectory: str, topology: str, weights: numpy.ndarray | None = None, output_path: str | None = os.getcwd()) -> pandas.DataFrame

   Calculate secondary structure assignment frequencies from a trajectory and topology files.

   :param trajectory: Path to .xtc trajectory file.
   :type trajectory: :py:class:`str`
   :param topology: Path to .pdb topology file.
   :type topology: :py:class:`str`
   :param weights: Optional array of weight values to be used in secondary structure
                   assignment reweighting. If None, uniform weights are used.
   :type weights: :py:class:`np.ndarray, optional`
   :param output_path: Path to output .csv file or output directory. If directory, written file is named
                       'ss_frequency.csv'. Defaults to current working directory.
   :type output_path: :py:class:`str, optional`

   :returns:     Secondary structure frequencies matrix for trajectory being analyzed.
   :rtype: pd.DataFrame


.. py:function:: calc_rg(u: MDAnalysis.Universe) -> float

   Calculate the radius of gyration of the current frame.

   :param u: Universe pointing to the current frame.
   :type u: :py:class:`mda.Universe`

   :returns:     Radius of gyration of the protein in the current frame.
   :rtype: float


.. py:function:: calc_eed(u: MDAnalysis.Universe) -> float

   Calculate the distance from the N to the C terminal in the current frame.

   :param u: Universe pointing to the current frame.
   :type u: :py:class:`mda.Universe`

   :returns:     End-to-end distance of the protein in the current frame.
   :rtype: float


.. py:function:: calc_dmax(u: MDAnalysis.Universe) -> float

   Calculate the maximum of the distances between any two alpha carbons in the current frame.

   :param u: Universe pointing to the current frame.
   :type u: :py:class:`mda.Universe`

   :returns:     Maximum of the distances between any two alpha carbons of the protein in the current
                 frame.
   :rtype: float


.. py:function:: calc_cm_dist(u: MDAnalysis.Universe, sel1: str, sel2: str) -> float

   Calculate the distance between the center of mass of two atom selections in current frame.

   :param u: Universe pointing to the current frame.
   :type u: :py:class:`mda.Universe`
   :param sel1: MDAnalysis selection string for selecting an AtomGroup whose center of mass will be
                calculated.
   :type sel1: :py:class:`str`
   :param sel2: MDAnalysis selection string for selecting an AtomGroup whose center of mass will be
                calculated.
   :type sel2: :py:class:`str`

   :returns:     Center of mass distance between AtomGroups selected by sel1 and sel2.
   :rtype: float


.. py:function:: calculate_metrics_data(trajectory: str, topology: str, rg: bool | None = True, dmax: bool | None = True, eed: bool | None = True, cm_dist: dict[str, tuple[str, str]] | None = None, output_path: str | None = os.getcwd()) -> pandas.DataFrame

   Calculate structural metrics for each frame of a trajectory.

   :param trajectory: Path to .xtc trajectory file.
   :type trajectory: :py:class:`str`
   :param topology: Path to .pdb topology file.
   :type topology: :py:class:`str`
   :param rg: Whether to calculate the radius of gyration of the protein.
   :type rg: :py:class:`bool, optional`
   :param dmax: Whether to calculate the maximum distance between any two alpha carbons in the protein.
   :type dmax: :py:class:`bool, optional`
   :param eed: Whether to calculate the distance from the N to C terminal of the protein.
   :type eed: :py:class:`bool, optional`
   :param cm_dist: Mapping of identifiers to tuples with two selection strings for creating MDAnalysis
                   AtomGroups, whose center mass distance will be calculated. For example:

                       {'inter_domain' : ('resid 1:30', 'resid 110:140')}

                   If None, no center mass distances are calculated.
                   See https://userguide.mdanalysis.org/stable/selections.html for more information about
                   MDAnalysis selections.
   :type cm_dist: :py:class:`dict[str,tuple[str,str]], optional`
   :param output_path: Path to output .csv file or output directory. If directory, written file is named
                       'structural_metrics.csv'. Defaults to current working directory.
   :type output_path: :py:class:`str, optional`

   :returns:     DataFrame where columns are the desired structural metrics and rows are the frames
                 of the trajectory.
   :rtype: pd.DataFrame


.. py:function:: calculate_analysis_data(trajectories: list[str], topologies: list[str], trajectory_ids: list[str], output_directory: str | None = os.getcwd(), ramachandran_data: bool = True, distancematrices: bool = True, contactmatrices: bool = True, ssfrequencies: bool = True, rg: bool = True, dmax: bool = True, eed: bool = True, cm_dist: dict[str, tuple[str, str]] | None = None) -> dict[str, list[pandas.DataFrame]]

   Calculate  structural data for each given pair of trajectory,topology files.

   :param trajectories: List of paths to .xtc trajectory files.
   :type trajectories: :py:class:`list[str]`
   :param topologies: List of paths to .pdb topology files.
   :type topologies: :py:class:`list[str]`
   :param trajectory_ids: Prefix trajectory identifiers to distinguish between calculated data files.
   :type trajectory_ids: :py:class:`list[str]`
   :param output_directory: Path to directory where calculated data will be stored. Defaults to current
                            working directory.
   :type output_directory: :py:class:`str, optional`
   :param ramachandran_data: Whether to calculate a dihedral angles matrix for each trajectory,topology
                             file pair.
   :type ramachandran_data: :py:class:`bool`
   :param distancematrices: Whether to calculate an alpha carbon distance matrix for each trajectory,topology
                            file pair.
   :type distancematrices: :py:class:`bool`
   :param contactmatrices: Whether to calculate a contact frequency matrix for each trajectory,topology
                           file pair.
   :type contactmatrices: :py:class:`bool`
   :param ssfrequencies: Whether to calculate a secondary structure assignment frequency matrix for each
                         trajectory, topology file pair.
   :type ssfrequencies: :py:class:`bool`
   :param rg: Whether to calculate the radius of gyration for each trajectory,topology file pair.
   :type rg: :py:class:`bool`
   :param dmax: Whether to calculate the maximum distance between any two alpha carbons for each
                trajectory,topology file pair.
   :type dmax: :py:class:`bool`
   :param eed: Whether to calculate the distance between the N and C terminal for each trajectory,
               topology file pair.
   :type eed: :py:class:`bool`
   :param cm_dist: Mapping of identifiers to tuples with two selection strings for creating MDAnalysis
                   AtomGroups, whose center mass distance will be calculated. If None, no center mass
                   distances are calculated. See https://userguide.mdanalysis.org/stable/selections.html
                   for more information about MDAnalysis selections. For example:

                   {'inter_domain' : ('resid 1:30', 'resid 110:140')}
   :type cm_dist: :py:class:`dict[str,tuple[str,str]], optional`

   :returns:     Mapping of data identifiers to lists of DataFrames with the calculated analysis data,
                 one element for each given trajectory,topology,trajectory_id trio. For example:

                 data = {
                 'DistanceMatrices' : [DistanceMatrix1,DistanceMatrix2,DistanceMatrix3],
                 'ContactMatrices' : [ContactMatrix1,ContactMatrix2,ContactMatrix3],
                 'SecondaryStructureFrequencies' : [SSFrequency1,SSFrequency2,SSFrequency3],
                 'StructuralMetrics' : [StructuralMetrics1,StructuralMetrics2, StructuralMetrics3]}
   :rtype: dict[str,list[pd.DataFrame]]


