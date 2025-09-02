"""
Analysis - ``ensemblify.analysis``
==================================

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
The ``ensemblify.analysis.analyze_trajectory`` function can be used to create an interactive
analysis graphical dashboard from a conformational ensemble in .xtc trajectory format.

For example, we can analyze an ensemble of Histatin5 (Hst5), an intrinsically disordered peptide
with 24 aminoacid residues.

Assuming that the path to the required trajectory and topology files are assigned to the variables
HST5_TRAJ_PATH and HST5_TOP_PATH, respectively, and that your desired trajectory ID is 'Hst5', you
can run:

>>> import ensemblify.analysis as ea
>>> ea.analyze_trajectory(HST5_TRAJ_PATH, HST5_TOP_PATH, 'Hst5')


Available Functions
-------------------

``analyze_trajectory``
    Calculate structural data and create interactive figures for given trajectory and topology
    files.
``calculate_analysis_data``
    Calculate  structural data for each given pair of trajectory,topology files.
``create_analysis_figures``
    Create interactive figures given analysis data for one or more pairs of trajectory,topology
    files.
``calculate_ramachandran_data``
    Calculate a dihedral angles matrix from trajectory and topology files.
``calculate_contact_matrix``
    Calculate a contact frequency matrix from a trajectory and topology files.
``calculate_distance_matrix``
    Calculate an alpha carbon average distance matrix from a trajectory and topology files.
``calculate_ss_assignment``
    Calculate a secondary structure assignment matrix from a trajectory and topology files.
``calculate_ss_frequency``
    Calculate secondary structure assignment frequencies from a trajectory and topology files.
``calculate_metrics_data``
    Calculate structural metrics for each frame of a trajectory.
``create_ramachandran_figure``
    Create a ramachandran plot Figure from a calculated dihedral angles matrix.
``create_contact_map_fig``
    Create a contact map Figure from a calculated contact matrix.
``create_distance_matrix_fig``
    Create a distance matrix Figure from a calculated distance matrix.
``create_ss_frequency_figure``
    Create a secondary structure frequency Figure from a secondary structure assignment frequency
    matrix.
``create_metrics_traces``
    Create Plotly Box, Histogram and Scatter (KDE) traces from calculated structural metrics data
    to be used in the creation of a Structural Metrics Figure.
``create_metrics_fig``
    Create a Structural Metrics Figure from previously created Box, Histogram and Scatter traces.
``create_single_metrics_fig_directly``
    Create a Structural Metrics Figure for a single trajectory, directly from its calculated
    metrics data and trajectory ID.
"""

from ensemblify.analysis.trajectory import analyze_trajectory
from ensemblify.analysis.data import (
    calculate_analysis_data, calculate_contact_matrix, calculate_distance_matrix,
    calculate_metrics_data, calculate_ramachandran_data, calculate_ss_assignment,
    calculate_ss_frequency,
)
from ensemblify.analysis.figures import (
    create_analysis_figures, create_contact_map_fig,
    create_distance_matrix_fig, create_metrics_fig, create_metrics_traces,
    create_ramachandran_figure, create_single_metrics_fig_directly,
    create_ss_frequency_figure
)

__all__ = ['analyze_trajectory',
           'calculate_analysis_data',
           'calculate_contact_matrix',
           'calculate_distance_matrix',
           'calculate_metrics_data',
           'calculate_ramachandran_data',
           'calculate_ss_assignment',
           'calculate_ss_frequency',
           'create_analysis_figures',
           'create_contact_map_fig',
           'create_distance_matrix_fig',
           'create_metrics_fig',
           'create_metrics_traces',
           'create_ramachandran_figure',
           'create_single_metrics_fig_directly',
           'create_ss_frequency_figure']
