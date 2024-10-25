"""Analyze a trajectory,topology pair, outputting a graphical dashboard."""

# IMPORTS
## Standard Library Imports
import os

## Third Party Imports
import pandas as pd
from plotly.offline import get_plotlyjs

## Local Imports
from ensemblify.analysis.trajectory_utils import calculate_analysis_data, create_analysis_figures
from ensemblify.config import GLOBAL_CONFIG

# FUNCTIONS
def analyze_trajectory(
    trajectories: list[str] | str,
    topologies: list[str] | str,
    trajectory_ids: list[str] | str,
    output_directory: str = None,
    ramachandran_data: bool = True,
    distancematrices: bool = True,
    contactmatrices: bool = True,
    ssfrequencies: bool = True,
    rg: bool = True,
    dmax: bool = True,
    eed: bool = True,
    cm_dist: dict | None = None,
    color_palette: list[str] | None = None,
    ) -> dict[str,pd.DataFrame]:
    """Calculate structural data and create interactive figures for given trajectory and
    topology files.

    Args:
        trajectories:
            list of paths to .xtc trajectory files or string with the path to a single .xtc
            trajectory file.
        topologies:
            list of paths to .pdb topology files or string with the path to a single .pdb
            topology file.
        trajectory_ids:
            list of prefix trajectory identifiers to distinguish between calculated data
            files or string with a single prefix trajectory identifier.
        output_directory:
            path to directory where calculated data and created figures will be stored.
            If it does not exist, it is created. Defaults to current working directory.
        ramachandran_data:
            whether to calculate a dihedral angles matrix for each trajectory,topology 
            file pair.
        distancematrices:
            whether to calculate an alpha carbon distance matrix for each trajectory,topology
            file pair and create the corresponding distance matrix interactive figure.
        contactmatrices:
            whether to calculate a contact frequency matrix for each trajectory,topology
            file pair and create the corresponding contact map interactive figure.
        ssfrequencies:
            whether to calculate a secondary structure assignment frequency matrix for each
            trajectory, topology file pair and create the corresponding secondary structure
            frequency interactive figure.
        rg:
            whether to calculate and plot the radius of gyration for each trajectory,topology
            file pair.
        dmax:
            whether to calculate and plot the maximum distance between any two alpha carbons
            for each trajectory,topology file pair.
        eed:
            whether to calculate and plot the distance between the N and C terminal for each
            trajectory, topology file pair.
        cm_dist:
            mapping of identifiers to tuples with two selection strings for creating MDAnalysis
            AtomGroups, whose center mass distance will be calculated and plotted. For example:
                {'inter_domain' : ('resid 1:30', 'resid 110:140')}
            If None, no center mass distances are calculated.
            See https://userguide.mdanalysis.org/stable/selections.html for more information about
            MDAnalysis selections.
        color_palette:
            list of color hexcodes, to associate one with each trajectory in the created Structural
            Metrics interactive dashboard.

    Returns:
        analysis_data:
            mapping of data identifiers to DataFrames containing the calculated analysis data for
            each frame of each given trajectory. For convenience, this is returned as a variable
            as well as saved to output directory.
    """
    # Setup trajectory,topology,trajectory_id
    if isinstance(trajectories,str):
        trajectories = list(trajectories)
    if isinstance(topologies,str):
        topologies = list(topologies)
    if isinstance(trajectory_ids,str):
        trajectory_ids = list(trajectory_ids)

    # Setup output directory
    if output_directory is None:
        output_directory = os.getcwd()
    elif not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Setup color palette
    if color_palette is None:
        color_palette = ['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                         '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52']

    # Calculate analysis data
    analysis_data = calculate_analysis_data(trajectories=trajectories,
                                            topologies=topologies,
                                            trajectory_ids=trajectory_ids,
                                            output_directory=output_directory,
                                            ramachandran_data=ramachandran_data,
                                            distancematrices=distancematrices,
                                            contactmatrices=contactmatrices,
                                            ssfrequencies=ssfrequencies,
                                            rg=rg,
                                            dmax=dmax,
                                            eed=eed,
                                            cm_dist=cm_dist)

    # Create Figures
    figures = create_analysis_figures(analysis_data=analysis_data,
                                      topologies=topologies,
                                      trajectory_ids=trajectory_ids,
                                      output_directory=output_directory,
                                      color_palette=color_palette)

    ################################
    ### Offline Interactive Html ###
    ################################
    print(f'Building {trajectory_ids} analysis dashboard...')
    # Build divs for html
    # DISTANCE MATRICES
    distance_matrices = figures['DistanceMatrices']
    if distance_matrices:
        distance_matrices_htmls = []
        for i,cmap_fig in enumerate(distance_matrices):
            distance_matrices_htmls.append(cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                            full_html=False,
                                                            include_plotlyjs=False,
                                                            div_id=f'distance_matrix_{i}'))
        distance_matrices_div = '\n'.join(distance_matrices_htmls)
    else:
        distance_matrices_div = ''

    # CONTACT MAPS
    contact_maps = figures['ContactMaps']
    if contact_maps:
        contact_maps_htmls = []
        for i,cmap_fig in enumerate(contact_maps):
            contact_maps_htmls.append(cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                       full_html=False,
                                                       include_plotlyjs=False,
                                                       div_id=f'contact_map_{i}'))
        contact_maps_div = '\n'.join(contact_maps_htmls)
    else:
        contact_maps_div = ''

    # SECONDARY STRUCTURE FREQUENCIES
    ss_freqs = figures['SecondaryStructureFrequencies']
    if ss_freqs:
        ss_freqs_htmls = []
        for i,ss_freqs_fig in enumerate(ss_freqs):
            ss_freqs_htmls.append(ss_freqs_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                       full_html=False,
                                                       include_plotlyjs=False,
                                                       div_id=f'ss_frequency_{i}'))
        ss_assigns_div = '\n'.join(ss_freqs_htmls)
    else:
        ss_assigns_div = ''

    # STRUCTURAL METRICS
    metrics_fig = figures['StructuralMetrics']
    if metrics_fig is not None:
        metrics_div = metrics_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                          full_html=False,
                                          include_plotlyjs=False,
                                          div_id='metrics')
    else:
        metrics_div = ''

    # Build html
    dashboard_html = f'''
    <!DOCTYPE html>
    <html>
        <head>
            <script type="text/javascript">{get_plotlyjs()}</script>
            <style>
                .flex-container {{
                    display: flex;
                    flex-wrap: wrap;  /* wrap or nowrap */
                }}
                .metricsdiv {{
                    width: auto;
                }}
            </style>
        </head>
        <body>
            <div class="flex-container">
                {contact_maps_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {distance_matrices_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {ss_assigns_div}
            </div>
            <div class="metricsdiv">
            {metrics_div}
            </div>
        </body>
    </html>
    '''

    # Save dashboard
    with open(os.path.join(output_directory,'analysis_dashboard.html'),'w',encoding='utf-8') as f:
        f.write(dashboard_html)

    print('Ensemble analysis calculation has finished. Please consult the interactive '
          'analysis_dashboard.html figure.')

    # For convenience, return the calculated analysis data.
    return analysis_data
