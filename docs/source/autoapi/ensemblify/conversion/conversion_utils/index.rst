ensemblify.conversion.conversion_utils
======================================

.. py:module:: ensemblify.conversion.conversion_utils

.. autoapi-nested-parse::

   Auxiliary functions for the conversion module.



Functions
---------

.. autoapisummary::

   ensemblify.conversion.conversion_utils.move_topol_pdb
   ensemblify.conversion.conversion_utils.join_pdbs
   ensemblify.conversion.conversion_utils.calc_saxs_data


Module Contents
---------------

.. py:function:: move_topol_pdb(job_name: str, origin_path: str, destination_path: str) -> str

   Move a .pdb file from origin to destination directory.

   Given a path to a directory containing .pdb files, moves a single .pdb
   from that directory to a destination directory, to serve as a topology
   file for future trajectory analysis.

   :param job_name: Prefix identifier for moved .pdb file.
   :type job_name: :py:class:`str`
   :param origin_path: Directory where topology file of interest is located.
   :type origin_path: :py:class:`str`
   :param destination_path: Directory where topology file will be moved to.
   :type destination_path: :py:class:`str`

   :returns:     Path to the moved topology file.
   :rtype: str


.. py:function:: join_pdbs(pdbs_dir: str, job_name: str, ensemble_dir: str, n_models: int) -> str

   Join a randomly sampled number of .pdb files in a directory into a single multimodel
   .pdb file.

   :param job_name: Prefix identifier for created multimodel .pdb file.
   :type job_name: :py:class:`str`
   :param ensemble_dir: Path to directory where ensemble pdb will be created.
   :type ensemble_dir: :py:class:`str`
   :param pdbs_dir: Path to directory where numbered .pdb files are stored.
   :type pdbs_dir: :py:class:`str`
   :param n_models: Number of .pdb files to randomly sample from the specified directory.
   :type n_models: :py:class:`int`

   :returns:     Path to created multimodel ensemble .pdb file.
   :rtype: str


.. py:function:: calc_saxs_data(trajectory_id: str, universe: MDAnalysis.Universe, frame_index: int, exp_saxs_file: str, calc_saxs_log: str) -> numpy.ndarray

   Calculate a theoretical SAXS curve for a frame of a MDAnalysis
   Universe object using PEPSI-SAXS.

   Calculation is done for a single temporary .pdb file created  from
   the current frame of the Universe object.

   :param trajectory_id: Prefix identifier for created files.
   :type trajectory_id: :py:class:`str`
   :param universe: Universe object containing the trajectory being analyzed.
   :type universe: :py:class:`mda.Universe`
   :param frame_index: Current frame being calculated.
   :type frame_index: :py:class:`int`
   :param exp_saxs_file: Path to .dat file with experimental SAXS data for the current
                         protein to be used by PEPSI-SAXS.
   :type exp_saxs_file: :py:class:`str`
   :param calc_saxs_log: Path to .log file for the SAXS curve calculation of each frame.
   :type calc_saxs_log: :py:class:`str`

   :returns:     Numpy ndarray with the values for the SAXS curve calculated from this
                 frame of the trajectory in the Universe object.
   :rtype: np.ndarray


