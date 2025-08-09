"""Analyze a trajectory,topology pair, outputting a graphical dashboard."""

# IMPORTS
## Standard Library Imports
import os

## Third Party Imports
import pandas as pd
from plotly.offline import get_plotlyjs

## Local Imports
from ensemblify.analysis.data import calculate_analysis_data
from ensemblify.analysis.figures import create_analysis_figures
from ensemblify.analysis.colors import DEFAULT_TRACE_COLOR_PALETTE
from ensemblify.config import GLOBAL_CONFIG

# FUNCTIONS
def analyze_trajectory(
    trajectories: list[str] | str,
    topologies: list[str] | str,
    trajectory_ids: list[str] | str,
    output_directory: str | None = None,
    ramachandran_data: bool | None = True,
    distancematrices: bool | None = True,
    contactmatrices: bool | None = True,
    ssfrequencies: bool | None = True,
    rg: bool | None = True,
    dmax: bool | None = True,
    eed: bool | None = True,
    cm_dist: dict[str,tuple[str,str]] | None = None,
    color_palette: list[str] | None = None,
    ) -> dict[str,pd.DataFrame]:
    """Calculate structural data and create figures for trajectory,topology file pairs.

    Args:
        trajectories (list[str] | str):
            List of paths to .xtc trajectory files or string with the path to a single .xtc
            trajectory file.
        topologies (list[str] | str):
            List of paths to .pdb topology files or string with the path to a single .pdb
            topology file.
        trajectory_ids (list[str] | str):
            List of prefix trajectory identifiers to distinguish between calculated data
            files or string with a single prefix trajectory identifier.
        output_directory (str, optional):
            Path to directory where calculated data and created figures will be stored.
            If it does not exist, it is created. Defaults to current working directory.
        ramachandran_data (bool, optional):
            Whether to calculate a dihedral angles matrix for each trajectory,topology
            file pair.
        distancematrices (bool, optional):
            Whether to calculate an alpha carbon distance matrix for each trajectory,topology
            file pair and create the corresponding distance matrix interactive figure.
        contactmatrices (bool, optional):
            Whether to calculate a contact frequency matrix for each trajectory,topology
            file pair and create the corresponding contact map interactive figure.
        ssfrequencies (bool, optional):
            Whether to calculate a secondary structure assignment frequency matrix for each
            trajectory,topology file pair and create the corresponding secondary structure
            frequency interactive figure.
        rg (bool, optional):
            Whether to calculate and plot a probability distribution for the radius of gyration
            for each trajectory,topology file pair.
        dmax (bool, optional):
            Whether to calculate and plot a probability distribution for the maximum distance
            between any two alpha carbons for each trajectory,topology file pair.
        eed (bool, optional):
            Whether to calculate and plot a probability distribution for the distance between
            the N- and C-terminal for each trajectory,topology file pair.
        cm_dist (dict[str,tuple[str,str]], optional):
            Mapping of arbitrary string identifiers to tuples containing two selection strings
            for creating MDAnalysis AtomGroups. A probability distribution for the center mass
            distance between the two AtomGroups will be calculated and plotted. For example, to
            calculate the distance between the centers of mass of two domains, one comprising
            residues 1-30 and the other comprising residues 110-140:

                {'inter_domain' : ('resid 1:30', 'resid 110:140')}

            If None, no center mass distances are calculated.
            See https://userguide.mdanalysis.org/stable/selections.html for more information about
            MDAnalysis selection strings.
        color_palette (list[str], optional):
            List of color hexcodes, to associate one with each trajectory,topology pair in the
            created figures. If None, the default color palette is used.

    Returns:
        dict[str,pd.DataFrame]:
            Mapping of data identifiers to DataFrames containing the calculated analysis data for
            each frame of each given trajectory. For convenience, this is returned as a variable
            as well as saved to output directory.
    """
    # Setup trajectory,topology,trajectory_id
    if isinstance(trajectories,str):
        trajectories = [trajectories]
    if isinstance(topologies,str):
        topologies = [topologies]
    if isinstance(trajectory_ids,str):
        trajectory_ids = [trajectory_ids]

    # Setup output directory
    if output_directory is None:
        output_directory = os.getcwd()
    else:
        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

    # Setup color palette
    if color_palette is None:
        color_palette = DEFAULT_TRACE_COLOR_PALETTE

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
    # Build HTML divs

    # DISTANCE MATRICES
    distance_matrices = figures['DistanceMatrices']
    if distance_matrices:
        distance_matrices_htmls = []
        for i,cmap_fig in enumerate(distance_matrices):
            distance_matrices_htmls.append(
                cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                 full_html=False,
                                 include_plotlyjs=False,
                                 div_id=f'distance_matrix_{i}')
                )
        distance_matrices_div = f'{"&nbsp;" * 6}'.join(distance_matrices_htmls)
    else:
        distance_matrices_div = ''

    # CONTACT MAPS
    contact_maps = figures['ContactMaps']
    if contact_maps:
        contact_maps_htmls = []
        for i,cmap_fig in enumerate(contact_maps):
            contact_maps_htmls.append(
                cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                 full_html=False,
                                 include_plotlyjs=False,
                                 div_id=f'contact_map_{i}')
                )
        contact_maps_div = f'{"&nbsp;" * 6}'.join(contact_maps_htmls)
    else:
        contact_maps_div = ''

    # SECONDARY STRUCTURE FREQUENCIES
    ss_freqs = figures['SecondaryStructureFrequencies']
    if ss_freqs:
        ss_freqs_htmls = []
        for i,ss_freqs_fig in enumerate(ss_freqs):
            ss_freqs_htmls.append(
                ss_freqs_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                     full_html=False,
                                     include_plotlyjs=False,
                                     div_id=f'ss_frequency_{i}')
                )
        ss_assigns_div = f'{"&nbsp;" * 6}'.join(ss_freqs_htmls)
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

    # Build html dashboard
    multi_traj_option = f'''
        <div class="flex-container" id="cmaps_container">
        {contact_maps_div}
        </div>
        <div class="flex-container" id="dmatrices_container">
        {distance_matrices_div}
        </div>
        <div class="flex-container" id="ssfreqs_container">
        {ss_assigns_div}
        </div>
    '''
    single_traj_option = f'''
        <div class="flex-container">
        {contact_maps_div}
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        {distance_matrices_div}
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        {ss_assigns_div}
        </div>
    '''

    if len(contact_maps) > 1 or len(distance_matrices) > 1 or len(ss_freqs) > 1:
        figures_div = multi_traj_option
    else:
        figures_div = single_traj_option

    dashboard_html = f'''
    <!DOCTYPE html>
    <html>
        <head>
            <script type="text/javascript">{get_plotlyjs()}</script>
            <style>
                @import url("https://fonts.googleapis.com/css2?family=Jersey+25&display=swap");
        
                header {{
                    padding-top:5px;
                    padding-bottom: 5px;
                    font-family: "Jersey 25", sans-serif;
                    font-size: 140px;
                    text-align: center;
                    background: radial-gradient(circle, rgba(0,38,66,1) 25%, rgba(48,77,109,1) 100%);
                    color: white;
                    display: inline-block;
                    width: 100%;
                    border-radius: 10px;
                }}
            
                .header-container {{
                    text-align: center;
                }}
                
                .flex-container {{
                    display: flex;
                    justify-content: center;
                    align-items: center;
                    flex-wrap: wrap;
                    column-gap: 10px;
                    row-gap: 20px;
                }}
                .metricsdiv {{
                    display:flex;
                    width: auto;
                    justify-content: center;
                    align-items: center;
                }}
            </style>
        </head>
        <body>
            <div class="header-container">
            <header>Ensemblify Analysis</header>
            </div>
            <br>
            {figures_div}
            <br><br>
            <div class="metricsdiv">
            {metrics_div}
            </div>
            <br>
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
