ensemblify.analysis
===================

.. py:module:: ensemblify.analysis

.. autoapi-nested-parse::

   Analysis - `ensemblify.analysis`
   ================================

   :Author(s): Nuno P. Fernandes
   :Year: 2024
   :Copyright: GNU Public License v3

   .. versionadded:: 1.0.0

   This module contains functions for performing structural analysis on an ensemble of protein
   conformations, outputting an interactive graphical dashboard.

   Example applications
   --------------------

   Analyze an ensemble in trajectory format
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   The `ensemblify.analysis.analyze_trajectory` function can be used to create an interactive
   analysis graphical dashboard from a conformational ensemble in .xtc trajectory format.

   For example, we can analyze an ensemble of Histatin5, an intrinsically disordered protein
   (IDP) with 24 aminoacid residues.

   Assuming that the path to the required Hst5 .xtc trajectory
   file is assigned to a variable named HST5_TRAJECTORY_PATH, you should run:

   >>> import ensemblify.analysis as ea
   >>> ea.analyze_trajectory(HST5_TRAJECTORY_PATH)


   Available Functions
   -------------------

   - `analyze_trajectory`

         Calculate structural data and create interactive figures for given trajectory and topology
         files.

   - `calculate_analysis_data`

         Calculate  structural data for each given pair of trajectory,topology files.

   - `create_analysis_figures`

         Create interactive figures given analysis data for one or more pairs of trajectory,topology
         files.

   - `calculate_ramachandran_data`

         Calculate a dihedral angles matrix from trajectory and topology files.

   - `calculate_contact_matrix`

         Calculate a contact frequency matrix from a trajectory and topology files.

   - `calculate_distance_matrix`

         Calculate an alpha carbon average distance matrix from a trajectory and topology files.

   - `calculate_ss_assignment`

         Calculate a secondary structure assignment matrix from a trajectory and topology files.

   - `calculate_ss_frequency`

         Calculate secondary structure assignment frequencies from a trajectory and topology files.

   - `calculate_metrics_data`

         Calculate structural metrics for each frame of a trajectory.

   - `create_ramachandran_figure`

         Create a ramachandran plot Figure from a calculated dihedral angles matrix.

   - `create_contact_map_fig`

         Create a contact map Figure from a calculated contact matrix.

   - `create_distance_matrix_fig`

         Create a distance matrix Figure from a calculated distance matrix.

   - `create_ss_frequency_figure`

         Create a secondary structure frequency Figure from a secondary structure assignment frequency
         matrix.

   - `create_metrics_traces`

         Create Ploty Box, Histogram and Scatter (KDE) traces from calculated structural metrics data
         to be used in the creation of a Structural Metrics Figure.

   - `create_metrics_fig`

         Create a Structural Metrics Figure from previously created Box, Histogram and Scatter traces.

   - `create_single_metrics_fig_directly`

         Create a Structural Metrics Figure for a single trajectory, directly from its calculated
         metrics data and trajectory ID.



Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/ensemblify/analysis/colors/index
   /autoapi/ensemblify/analysis/data/index
   /autoapi/ensemblify/analysis/figures/index
   /autoapi/ensemblify/analysis/trajectory/index


