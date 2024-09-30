"""A one-line summary of the module or program, terminated by a period.

Leave one blank line.  The rest of this docstring should contain an
overall description of the module or program.  Optionally, it may also
contain a brief description of exported classes and functions and/or usage
examples.

Typical usage example:

  foo = ClassFoo()
  bar = foo.FunctionBar()
"""

from ensemblify.analysis.trajectory import analyze_trajectory
from ensemblify.analysis.trajectory_utils import (calculate_analysis_data,create_analysis_figures,
                                                  calculate_ramachandran_data,calculate_distance_matrix,
                                                  calculate_contact_matrix,calculate_ss_assignment,
                                                  calculate_ss_frequency,calculate_metrics_data,
                                                  create_contact_map_fig,create_ss_frequency_figure,
                                                  create_metrics_fig)

__all__ = ['analyze_trajectory',
           'calculate_analysis_data',
           'create_analysis_figures',
           'calculate_ramachandran_data',
           'calculate_distance_matrix',
           'calculate_contact_matrix',
           'calculate_ss_assignment',
           'calculate_ss_frequency',
           'calculate_metrics_data',
           'create_contact_map_fig',
           'create_ss_frequency_figure',
           'create_metrics_fig']
