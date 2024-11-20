"""Auxiliary functions for creating analysis figures."""

# IMPORTS
## Standard Library Imports
import importlib.resources
import math
import os

## Third Party Imports
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

## Local Imports
from ensemblify.analysis.colors import PARULA_COLORSCALE, DEFAULT_TRACE_COLOR_PALETTE
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.utils import extract_pdb_info, get_array_extremum, kde, round_to_nearest_multiple

# CONSTANTS
DEFAULT_AXIS = {
    'showgrid': False,
    'color': 'black',
    'showline': True,
    'mirror': True,
    'linewidth': 4,
    'linecolor': 'black',
    'ticks': 'outside',
    'ticklen': 10,
    'tickwidth': 4,
}

DEFAULT_LAYOUT = {
    'plot_bgcolor': '#FFFFFF',
    'paper_bgcolor': '#FFFFFF',
    'font': {
        'family': 'Helvetica',
        'color': 'black',
        'size': 30,
    },
    'modebar_remove': ['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
    'xaxis': DEFAULT_AXIS,
    'yaxis': DEFAULT_AXIS,
}

# FUNCTIONS
def get_figure_layout_elements(
    topology: str,
    num_res: int,
    ssfreq: bool | None = False,
    ) -> tuple[list[str],list[str],list[dict],list[int]]:
    """Create tick labels and values and lines to divide chain regions for a Figure.

    Args:
        topology:
            topology file associated with the trajectory being analyzed.
        num_res:
            total number of residues in the analyzed system (equivalent to axis size in
            a Figure).
        ssfreq:
            Whether we are dealing with a secondary structure assignment frequency Figure or not.
            Defaults to False.

    Returns:
        A tuple x_labels,y_labels,chain_dividers,tickvals where:
            x_labels:
                all the possible tick labels for the x axis.
            y_labels:
                all the possible tick labels for the y axis. Only relevant when not building a
                secondary structure assignment frequency Figure.
            chain_dividers:
                dictionaries to create Plotly shapes, i.e. lines to add to the Figure that separate
                between different protein chains.
            tickvals:
                where to place the ticks in the x axis.
    """
    # Extract info regarding chains and resnums
    top_info = extract_pdb_info(topology)
    resranges = {}
    chain_letters = []

    # Start from the last chain as the .pdb was also parsed from the last res
    for chain_number in range(len(top_info.keys()),0,-1):
        chain_letter, starting_res, chain_size = top_info[chain_number]
        resranges[chain_letter] = [ x for x in range(starting_res, starting_res + chain_size)]
        chain_letters.append(chain_letter)

    # Create tick labels that respect chain id
    if len(chain_letters) > 1:
        x_labels = [f'{chain_letter}{resnum}' for chain_letter in chain_letters
                    for resnum in resranges[chain_letter]]

        y_labels = [f'{chain_letter}{resnum}' for chain_letter in chain_letters
                    for resnum in resranges[chain_letter]]
    else:
        x_labels = [f'{resnum}' for chain_letter in chain_letters
                    for resnum in resranges[chain_letter]]

        y_labels = [f'{resnum}' for chain_letter in chain_letters
                    for resnum in resranges[chain_letter]]

    # Setup chain dividers lines
    chain_dividers = []
    chain_ends = [] # to be used in tickvals
    chain_starts = [] # to be used in tickvals
    cumulative_residues = 0

    for chain_letter in chain_letters[:-1]:
        if ssfreq:
            chain_starts.append(cumulative_residues + 1)
        else:
            chain_starts.append(cumulative_residues)
        chain_size = len(resranges[chain_letter])
        chain_end = cumulative_residues + chain_size
        chain_ends.append(chain_end)

        if ssfreq:
            chain_dividers.append(dict(type='line',
                                       xref='x',
                                       x0=chain_end+0.5,
                                       x1=chain_end+0.5,
                                       y0=0,
                                       y1=1,
                                       line=dict(color='black',
                                                 width=2)))
        else:
            chain_dividers.append(dict(type='line',
                                       xref='x',
                                       x0=chain_end-0.5,
                                       x1=chain_end-0.5,
                                       y0=0-0.5,
                                       y1=num_res+0.5,
                                       line=dict(color='black',
                                                 width=2)))
            chain_dividers.append(dict(type='line',
                                       yref='y',
                                       y0=chain_end-0.5,
                                       y1=chain_end-0.5,
                                       x0=0-0.5,
                                       x1=num_res+0.5,
                                       line=dict(color='black',
                                                 width=2)))
        cumulative_residues += chain_size
    if ssfreq:
        chain_starts.append(num_res - len(resranges[chain_letters[-1]]) + 1)
    else:
        chain_starts.append(num_res - len(resranges[chain_letters[-1]]))
    chain_ends.append(num_res)

    # Setup tick values
    tickvals = []
    chain_counter = 0
    curr_val = chain_starts[chain_counter]
    tick_step = num_res // len(chain_letters) // 4 # always 5 ticks per axis

    while curr_val <= num_res:
        try:
            chain_end = chain_ends[chain_counter]
        except IndexError:
            tickvals.append(curr_val)
        else:
            if chain_end - curr_val <= tick_step:
                chain_counter += 1
                try:
                    curr_val = chain_starts[chain_counter]
                    tickvals.append(curr_val)
                except IndexError:
                    if chain_ends[-1] - curr_val <= 3:
                        if ssfreq:
                            tickvals.append(chain_ends[-1])
                        else:
                            tickvals.append(chain_ends[-1] - 1)
                    else:
                        tickvals.append(curr_val)
                        tickvals.append(chain_ends[-1])
            else:
                tickvals.append(curr_val)
        curr_val += tick_step

    return x_labels,y_labels,chain_dividers,tickvals


def create_ramachandran_figure(
    dihedrals_matrix: pd.DataFrame | str,
    trajectory_id: str | None = None,
    output_path: str | None = None,
    ) -> go.Figure:
    """Create a ramachandran plot Figure from a calculated dihedral angles matrix.

    Args:
        dihedrals_matrix:
            calculated dihedral angles matrix DataFrame or path to calculated matrix in .csv format.
            If difference is True, this should be the difference dihedral angles matrix between the
            uniformly weighted and the reweighted dihedral angles matrix.
        trajectory_id:
            used on Figure title and prefix for saved ramachandran plot filename. Defaults to None.
        output_path:
            path to output .html file or output directory where created Figure will be stored.
            If directory, written file is named 'ramachandran_plot.html', optionally with
            trajectory_id prefix. Defaults to None.

    Returns:
        rama_fig:
            Ploty Figure object displaying a ramachandran plot.
    """
    if isinstance(dihedrals_matrix,str):
        assert dihedrals_matrix.endswith('.csv'), ('Dihedral angles matrix file must '
                                                   'be in .csv format!')
        dihedrals_matrix = pd.read_csv(dihedrals_matrix,index_col=0)

    # Create Ramachandran Plot figure
    rama_fig = go.Figure(layout=DEFAULT_LAYOUT)

    # Add Ramachandran Reference Contours
    rama_ref_data = np.load(
        importlib.resources.files('ensemblify.analysis')
        .joinpath('rama_ref_data.npy')
        .open('rb')
        ).flatten()
    # Ramachandran Regions Reference:
    # https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/analysis/data/rama_ref_data.npy

    X, Y = np.meshgrid(np.arange(-180, 180, 4),
                       np.arange(-180, 180, 4))

    rama_fig.add_trace(go.Contour(x=X.flatten(),
                                  y=Y.flatten(),
                                  z=rama_ref_data,
                                  name='Show/Hide Allowed+MarginallyAllowed regions countours',
                                  line_width=3,
                                  contours=dict(start=1,
                                                end=17,
                                                size=16,
                                                showlabels=False,
                                                coloring='lines'),
                                  colorscale=['#FF7124',
                                              '#FF000D'],
                                  hoverinfo='none',
                                  showscale=False,
                                  showlegend=True))

    rama_fig.add_trace(go.Scattergl(x=dihedrals_matrix['Phi'],
                                    y=dihedrals_matrix['Psi'],
                                    mode='markers',
                                    marker=dict(color='black',
                                                size=0.6),
                                    name='dihedrals',
                                    showlegend=False))

    # Create quadrant dividers
    shapes = []
    shapes.append(dict(type='line',
                       xref='x',
                       x0=0,
                       x1=0,
                       y0=-180,
                       y1=180,
                       line=dict(color='black',
                                 width=2)))
    shapes.append(dict(type='line',
                       yref='y',
                       y0=0,
                       y1=0,
                       x0=-180,
                       x1=180,
                       line=dict(color='black',
                                 width=2)))

    # Update Figure Layout
    if trajectory_id is not None:
        rama_title = f'{trajectory_id} Ramachandran Plot'
    else:
        rama_title = 'Ramachandran Plot'

    rama_title = 'Ramachandran Plot'

    rama_fig.update_layout(width=900,
                           height=900,
                           title=dict(text=rama_title,
                                      xref='paper',
                                      x=0.5),
                           xaxis=dict(title='\u03A6', # Phi
                                      range=[-180,180],
                                      tick0=-180,
                                      dtick=60,
                                      ticksuffix='\u00B0'),
                           yaxis=dict(title='\u03A8', # Psi
                                      range=[-180,180],
                                      tick0=-180,
                                      dtick=60,
                                      ticksuffix='\u00B0',
                                      title_standoff=5),
                            legend=dict(orientation='h',
                                        xref='paper',
                                        xanchor='center',
                                        x=0.5,
                                        y=1.065,
                                        font_size=18),
                           shapes=shapes)

    if output_path is not None:
        # Save ramachandran plot
        if os.path.isdir(output_path):
            if trajectory_id is not None:
                output_filename = f'{trajectory_id}_ramachandran_plot.html'
            else:
                output_filename = 'ramachandran_plot.html'
            rama_fig.write_html(os.path.join(output_path,output_filename),
                                config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

        elif output_path.endswith('.html'):
            rama_fig.write_html(output_path,
                                config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

        else:
            print(('Ramachandran plot was not saved to disk, '
                    'output path must be a directory or .html filepath!'))

    return rama_fig


def create_contact_map_fig(
    contact_matrix: pd.DataFrame | str,
    topology: str,
    trajectory_id: str | None = None,
    output_path: str | None = None,
    reweighted: bool | None = False,
    difference: bool | None = False,
    ) -> go.Figure:
    """Create a contact map Figure from a calculated contact matrix.

    The topology provides information about number of chains, their chain letters and
    residue numbers.

    Args:
        contact_matrix:
            calculated contact matrix DataFrame or path to calculated matrix in .csv format.
            If difference is True, this should be the difference contact matrix between the
            uniformly weighted and the reweighted contact matrix.
        topology:
            path to topology .pdb file.
        trajectory_id:
            used on Figure title and prefix for saved contact map filename. Defaults to None.
        output_path:
            path to output .html file or output directory where created Figure will be stored.
            If directory, written file is named 'contact_map.html', optionally with
            trajectory_id prefix. Defaults to None.
        reweighted:
            boolean stating whether we are creating a reweighted contact map figure or a default
            one. Defaults to False.
        difference:
            boolean stating whether we are creating a difference contact map figure or a default
            one. Defaults to False.

    Returns:
        cmap_fig:
            Ploty Figure object displaying a contact map.
    """
    assert not(reweighted and difference), ('Contact Map Figure can\'t simultaneously be '
                                            'difference and reweighted!')

    if isinstance(contact_matrix,str):
        assert contact_matrix.endswith('.csv'), 'Contact matrix file must be in .csv format!'
        contact_matrix = pd.read_csv(contact_matrix,index_col=0)

    # Get x and y labels, chain dividers, and tickvals needed later for Figure creation
    x_labels,\
    y_labels,\
    chain_dividers,\
    tickvals = get_figure_layout_elements(topology=topology,
                                          num_res=len(contact_matrix.columns))

    # Create Contact Map Figure
    cmap_fig = go.Figure(layout=DEFAULT_LAYOUT)

    # Add our data, setup Figure title
    if not difference:
        contact_matrix.replace(0,np.nan,inplace=True)
        cmap_fig.add_trace(go.Heatmap(z=contact_matrix,
                                      zmin=0,
                                      zmax=1,
                                      x=x_labels,
                                      y=y_labels,
                                      transpose=True,
                                      hoverongaps=False,
                                      colorscale=px.colors.sequential.Reds))

        if reweighted:
            if trajectory_id is not None:
                cmap_title = f'{trajectory_id} Reweighted Contact Map'
            else:
                cmap_title = 'Reweighted Contact Map'
        else:
            if trajectory_id is not None:
                cmap_title = f'{trajectory_id} Contact Map'
            else:
                cmap_title = 'Contact Map'
    else:
        # Create hovertext
        hovertext = []
        for yi, yy in enumerate(x_labels):
            hovertext.append([])
            for xi, xx in enumerate(x_labels):
                value = round(contact_matrix.iat[yi,xi],3)
                if value != 0.0:
                    text = f'x: {xx}<br />y: {yy}<br />z: {value}'
                else:
                    text = f'x: {xx}<br />y: {yy}<br />z: {contact_matrix.iat[yi,xi]}'
                hovertext[-1].append(text)

        contact_matrix.replace(0,np.nan,inplace=True)
        cmap_fig.add_trace(go.Heatmap(z=contact_matrix,
                                      zmin=-1,
                                      zmax=1,
                                      x=x_labels,
                                      y=y_labels,
                                      hoverongaps=False,
                                      hoverinfo='text',
                                      hovertext=hovertext,
                                      colorscale=PARULA_COLORSCALE,
                                      zmid=0))

        if trajectory_id is not None:
            cmap_title = f'{trajectory_id} Difference Contact Map'
        else:
            cmap_title = 'Difference Contact Map'

    # Update Figure Layout
    cmap_fig.update_layout(width=900,
                           height=900,
                           title=dict(text=cmap_title,
                                      x=0.5,
                                      subtitle=dict(text=('Frequency of contacts between '
                                                          'any atoms of a residue pair'),
                                                    font=dict(color='gray',
                                                              size=24,))),
                           margin_t=150, # to fit subtitle
                           xaxis=dict(title='Residue',
                                      tickvals=tickvals,
                                      constrain='domain'),
                           yaxis=dict(title='Residue',
                                      title_standoff=5,
                                      tickvals=tickvals,
                                      constrain='domain',
                                      scaleanchor='x'),
                           shapes=chain_dividers)

    if output_path is not None:
        # Save contact map
        if os.path.isdir(output_path):
            if trajectory_id is not None:
                if reweighted:
                    output_filename = f'{trajectory_id}_contact_map_reweighted.html'
                elif difference:
                    output_filename = f'{trajectory_id}_contact_map_difference.html'
                else:
                    output_filename = f'{trajectory_id}_contact_map.html'
            elif reweighted:
                output_filename = 'contact_map_reweighted.html'
            elif difference:
                output_filename = 'contact_map_difference.html'
            else:
                output_filename = 'contact_map.html'
            cmap_fig.write_html(os.path.join(output_path,output_filename),
                                config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

        elif output_path.endswith('.html'):
            cmap_fig.write_html(output_path,
                                config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

        else:
            print(('Contact map was not saved to disk, '
                   'output path must be a directory or .html filepath!'))

    return cmap_fig


def create_distance_matrix_fig(
    distance_matrix: pd.DataFrame | str,
    topology: str,
    trajectory_id: str | None = None,
    output_path: str | None = None,
    max_colorbar: int | None = None,
    min_colorbar: int | None = None,
    reweighted: bool | None = False,
    difference: bool | None = False,
    ) -> go.Figure:
    """Create a distance matrix Figure from a calculated distance matrix.

    The topology provides information about number of chains, their chain letters and
    residue numbers.

    Args:
        distance_matrix:
            calculated distance matrix DataFrame or path to calculated matrix in .csv format.
            If difference is True, this should be the difference distance matrix between the
            uniformly weighted and the reweighted distance matrix.
        topology:
            path to topology .pdb file.
        trajectory_id:
            used on Figure title and prefix for saved distance matrix filename. Defaults to None.
        output_path:
            path to output .html file or output directory where created Figure will be stored.
            If directory, written file is named 'distance_matrix.html', optionally with
            trajectory_id prefix. Defaults to None.
        max_colorbar:
            maximum limit for the distance colorbar. Defaults to None, in which case it is
            derived from the data.
        min_colorbar:
            minimum limit for the distance colorbar. Defaults to None, in which case it is
            derived from the data.
        reweighted:
            boolean stating whether we are creating a reweighted distance matrix figure or a
            default one. Defaults to False.
        difference:
            boolean stating whether we are creating a difference distance matrix figure or a
            default one. Defaults to False.

    Returns:
        dmatrix_fig:
            Ploty Figure object displaying a distance matrix.
    """
    assert not(reweighted and difference), ('Distance Matrix Figure can\'t simultaneously be '
                                            'difference and reweighted!')
    if isinstance(distance_matrix,str):
        assert distance_matrix.endswith('.csv'), 'Distance matrix file must be in .csv format!'
        distance_matrix = pd.read_csv(distance_matrix,index_col=0)

    # Get x and y labels, chain dividers, and tickvals needed later for Figure creation
    x_labels,\
    y_labels,\
    chain_dividers,\
    tickvals = get_figure_layout_elements(topology=topology,
                                          num_res=len(distance_matrix.columns))

    # Create Distance Matrix Figure
    dmatrix_fig = go.Figure(layout=DEFAULT_LAYOUT)

    # Add our data
    if not difference:
        distance_matrix.replace(0,np.nan,inplace=True)

        if max_colorbar is None:
            max_colorbar = math.ceil(np.max(distance_matrix))

        dmatrix_fig.add_trace(go.Heatmap(z=distance_matrix,
                                         x=x_labels,
                                         y=y_labels,
                                         zmin=0,
                                         zmax=max_colorbar,
                                         colorbar_title='&#197;',
                                         colorbar_ticklabeloverflow='allow',
                                         transpose=True,
                                         hoverongaps=False,
                                         colorscale=px.colors.sequential.Reds))

        if reweighted:
            if trajectory_id is not None:
                dmatrix_title = f'{trajectory_id} Reweighted Distance Matrix'
            else:
                dmatrix_title = 'Reweighted Distance Matrix'
        else:
            if trajectory_id is not None:
                dmatrix_title = f'{trajectory_id} Distance Matrix'
            else:
                dmatrix_title = 'Distance Matrix'
    else:
        # Create hovertext
        hovertext = []
        for yi, yy in enumerate(x_labels):
            hovertext.append([])
            for xi, xx in enumerate(x_labels):
                value = round(distance_matrix.iat[yi,xi],3)
                if value != 0.0:
                    text = f'x: {xx}<br />y: {yy}<br />z: {value}'
                else:
                    text = f'x: {xx}<br />y: {yy}<br />z: {distance_matrix.iat[yi,xi]}'
                hovertext[-1].append(text)

        distance_matrix.replace(0,np.nan,inplace=True)
        dmatrix_fig.add_trace(go.Heatmap(z=distance_matrix,
                                         x=x_labels,
                                         y=y_labels,
                                         hoverongaps=False,
                                         colorbar_title='&#197;',
                                         hoverinfo='text',
                                         hovertext=hovertext,
                                         colorscale=PARULA_COLORSCALE,
                                         zmin=min_colorbar,
                                         zmid=0,
                                         zmax=max_colorbar))

        if trajectory_id is not None:
            dmatrix_title = f'{trajectory_id} Difference Distance Matrix'
        else:
            dmatrix_title = 'Difference Distance Matrix'

    # Update Figure Layout
    dmatrix_fig.update_layout(width=900,
                              height=900,
                              title=dict(text=dmatrix_title,
                                         x=0.5,
                                         subtitle=dict(text=('Average distance between alpha '
                                                             'carbons of a residue pair'),
                                                       font=dict(color='gray',
                                                                 size=24,))),
                              margin_t=150, # to fit subtitle
                              xaxis=dict(title='Residue',
                                         tickvals=tickvals,
                                         constrain='domain',),
                              yaxis=dict(title='Residue',
                                         tickvals=tickvals,
                                         title_standoff=5,
                                         scaleanchor='x',
                                         constrain='domain',),
                              shapes=chain_dividers)

    if output_path is not None:
        # Save distance matrix
        if os.path.isdir(output_path):
            if trajectory_id is not None:
                if reweighted:
                    output_filename = f'{trajectory_id}_distance_matrix_reweighted.html'
                elif difference:
                    output_filename = f'{trajectory_id}_distance_matrix_difference.html'
                else:
                    output_filename = f'{trajectory_id}_distance_matrix.html'
            elif reweighted:
                output_filename = 'distance_matrix_reweighted.html'
            elif difference:
                output_filename = 'distance_matrix_difference.html'
            else:
                output_filename = 'distance_matrix.html'
            dmatrix_fig.write_html(os.path.join(output_path,output_filename),
                                   config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

        elif output_path.endswith('.html'):
            dmatrix_fig.write_html(output_path,
                                config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

        else:
            print(('Distance Matrix was not saved to disk, '
                   'output path must be a directory or .html filepath!'))

    return dmatrix_fig


def create_ss_frequency_figure(
    ss_frequency: pd.DataFrame | str,
    topology: str,
    trajectory_id: str | None = None,
    output_path: str | None = None,
    reweighted: bool = False,
    difference: bool = False,
    ) -> go.Figure:
    """Create a secondary structure frequency Figure from a secondary structure assignment
    frequency matrix.

    The topology provides information about number of chains, their chain letters and
    residue numbers.

    Args:
        ss_frequency:
            calculated secondary structure assignment frequency matrix DataFrame or path to
            calculated matrix in .csv format.
        topology:
            path to topology .pdb file.
        trajectory_id:
            used on Figure title and prefix for saved ss_frequency filename. Defaults to None.
        output_path:
            path to output .html file or output directory where created Figure will be stored.
            If directory, written file is named 'ss_frequency.html', optionally with
            trajectory_id prefix. Defaults to None.
        reweighted:
            boolean stating whether we are creating a reweighted secondary structure frequency
            figure or a default one. Defaults to False.
        difference:
            boolean stating whether we are creating a difference secondary structure frequency
            figure or a default one. Defaults to False.

    Returns:
        ss_freq_fig:
            a stacked line plot with the secondary structure frequencies of each secondary
            structure type for each residue in the structure.
    """
    assert not(reweighted and difference), ('Secondary Structure Frequency Figure can\'t '
                                            'simultaneously be difference and reweighted!')

    if isinstance(ss_frequency,str):
        assert ss_frequency.endswith('.csv'), ('Secondary structure assignment frequency '
                                                'matrix must be in .csv format!')
        ss_frequency = pd.read_csv(ss_frequency,index_col=0)

    # Get x and y labels, chain dividers, and tickvals needed later for Figure creation
    x_labels,\
    _,\
    chain_dividers,\
    tickvals = get_figure_layout_elements(topology=topology,
                                          num_res=len(ss_frequency.columns),
                                          ssfreq=True)

    # Create Figure
    ss_freq_fig = go.Figure(layout=DEFAULT_LAYOUT)

    # Adding traces for each secondary structure type
    colors = ['#1f77b4',  # Blue
              '#2ca02c',  # Green
              '#d62728'  # Red
              ]

    if difference:
        for structure,color in zip(ss_frequency.index,colors):
            # Create hovertext
            hovertext = [f'x: {x_label}<br />y: {round(ss_frequency.loc[structure].iloc[i],5)}'
                          for i,x_label in enumerate(x_labels)]
            ss_freq_fig.add_trace(go.Scatter(x=list(range(1,len(ss_frequency.columns)+1)),
                                             y=ss_frequency.loc[structure],
                                             mode='lines',
                                             marker_color=color,
                                             line_width=4,
                                             name=structure,
                                             hoverinfo='text',
                                             hovertext=hovertext))

            # Setup Figure Layout
            if trajectory_id is not None:
                ss_freq_title = f'{trajectory_id} Difference Sec. Struct. Frequencies'
            else:
                ss_freq_title = 'Difference Secondary Structure Frequencies'
            range_y = [-1,1]
    else:
        for structure,color in zip(ss_frequency.index,colors):
            # Create hovertext
            hovertext = [f'x: {x_label}<br />y: {round(ss_frequency.loc[structure].iloc[i],5)}'
                          for i,x_label in enumerate(x_labels)]
            ss_freq_fig.add_trace(go.Scatter(x=list(range(1,len(ss_frequency.columns)+1)),
                                             y=ss_frequency.loc[structure],
                                             mode='lines',
                                             stackgroup='one', # remove for non stacked plot
                                             marker_color=color,
                                             line_width=0,
                                             name=structure,
                                             hoverinfo='text',
                                             hovertext=hovertext))

        # Setup Figure Layout
        if reweighted:
            if trajectory_id is not None:
                ss_freq_title = f'{trajectory_id} Reweighted Sec. Struct. Frequencies'
            else:
                ss_freq_title = 'Reweighted Secondary Structure Frequencies'
        else:
            if trajectory_id is not None:
                ss_freq_title = f'{trajectory_id} Secondary Structure Frequencies'
            else:
                ss_freq_title = 'Secondary Structure Frequencies'
        range_y = [0,1]

    # Setup tick text
    x_text = []
    for x_val in tickvals:
        x_t = x_labels[x_val-1]
        x_text.append(x_t)

    # Update Figure Layout
    ss_freq_fig.update_layout(width=1000,
                              height=750,
                              title=dict(text=ss_freq_title,
                                         x=0.5,
                                         subtitle=dict(text=('Frequency of each secondary '
                                                             'structure assignment code for '
                                                             'each residue '),
                                                       font=dict(color='gray',
                                                                 size=24))),
                              margin_t=150, # to fit subtitle
                              xaxis=dict(title='Residue',
                                         ticks='outside',
                                         tickvals=tickvals,
                                         ticktext=x_text),
                              yaxis=dict(title='Frequency',
                                         range=range_y),
                              shapes=chain_dividers)

    if output_path is not None:
        # Save Secondary Structure frequency
        if os.path.isdir(output_path):
            if trajectory_id is not None:
                if reweighted:
                    output_filename = f'{trajectory_id}_ss_frequency_reweighted.html'
                elif difference:
                    output_filename = f'{trajectory_id}_ss_frequency_difference.html'
                else:
                    output_filename = f'{trajectory_id}_ss_frequency.html'
            elif reweighted:
                output_filename = 'ss_frequency_reweighted.html'
            elif difference:
                output_filename = 'ss_frequency_difference.html'
            else:
                output_filename = 'ss_frequency.html'
            ss_freq_fig.write_html(os.path.join(output_path,output_filename),
                                   config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

        elif output_path.endswith('.html'):
            ss_freq_fig.write_html(output_path,
                                   config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])
        else:
            print(('Secondary structure frequency graph was not saved to disk, '
                   'output path must be a directory or .html filepath!'))

    return ss_freq_fig


def create_metrics_traces(
    metrics: pd.DataFrame | str,
    trajectory_id: str,
    color: str = '#636EFA',
    ) -> tuple[list[go.Box], list[go.Histogram], list[go.Scatter]]:
    """Create Ploty Box, Histogram and Scatter (KDE) traces to be used in the creation of
    a Structural Metrics Figure.

    Args:
        metrics:
            DataFrame where columns are the desired structural metrics and rows are the frames
            of the trajectory or path to that DataFrame in .csv format.
        trajectory_id:
            prefix identifier for trace names.
        color:
            hex code for the color the created traces will be. Defaults to '#636EFA', or light
            blue.

    Returns:
        A tuple (box_traces,hist_traces,scatter_traces) where:
            box_traces:
                a list of the boxplot traces, one for each structural metric.
            hist_traces:
                a list of the histogram traces, one for each structural metric.
            scatter_traces:
                a list of the scatter Kernel Density Estimate (KDE) traces, one
                for each structural metric.
            avg_values:
                a list of the values of the mean for each metric.
    """
    if isinstance(metrics,str):
        assert metrics.endswith('.csv'), ('Structural metrics matrix '
                                          'must be in .csv format!')
        metrics = pd.read_csv(metrics,index_col=0)

    # Create Box traces
    box_traces = []
    for col_name in metrics.columns:
        box_trace = go.Box(x=metrics[col_name],
                           name=f'{trajectory_id}_{col_name}',
                           orientation='h',
                           boxmean=True,
                           marker_color=color,
                           showlegend=False)
        box_traces.append(box_trace)

    # Create Histogram and Scatter (KDE) traces
    hist_traces = []
    scatter_traces = []
    avg_values = []
    avg_stderr_values = []

    for col_name in metrics.columns:
        data_array = np.array(metrics[col_name])

        hist_bin_range = (math.floor(np.min(data_array)), math.ceil(np.max(data_array)))
        hist_bin_edges = np.histogram_bin_edges(a=data_array,
                                                range=hist_bin_range,
                                                bins='fd')
        hist_bin_size = hist_bin_edges[1] - hist_bin_edges[0]

        hist_trace = go.Histogram(x=data_array,
                                  name=f'{trajectory_id}_{col_name}',
                                  xbins=dict(start=hist_bin_range[0],
                                             end=hist_bin_range[1],
                                             size=hist_bin_size),
                                  histnorm='probability',
                                  marker=dict(color=color),
                                  opacity=0.7,
                                  visible=False)
        hist_traces.append(hist_trace)

        kde_x, kde_y, avg, avg_stderr = kde(data=data_array)
        avg_values.append(avg)
        avg_stderr_values.append(avg_stderr)

        scatter_trace = go.Scatter(x=kde_x,
                                   y=kde_y,
                                   mode='lines',
                                   name=f'{trajectory_id}_{col_name}',
                                   marker_color=color,
                                   line=dict(width=4),
                                   legend='legend',
                                   visible=True)
        scatter_traces.append(scatter_trace)

    return box_traces,hist_traces,scatter_traces,avg_values,avg_stderr_values


def create_metrics_fig(
    trajectory_ids: list[str],
    total_box_traces: dict[str,list[go.Box]],
    total_hist_traces: dict[str,list[go.Histogram]],
    total_scatter_traces: dict[str,list[go.Scatter]],
    total_avg_values: dict[str,list[float]],
    total_avg_stderr_values: dict[str,list[float]],
    output_path: str | None = None,
    ) -> go.Figure:
    """Create a Structural Metrics Figure from previously created Box, Histogram and Scatter traces.

    Args:
        trajectory_ids:
            list of prefix identifiers that must match the prefix identifiers used for naming the
            created traces.
        total_box_traces:
            mapping of trajectory_ids to a list of created Box traces.
        total_hist_traces:
            mapping of trajectory_ids to a list of created Histogram traces.
        total_scatter_traces:
            mapping of trajectory_ids to a list of created Scatter traces.
        total_avg_values:
            mapping of trajectory_ids to a list of mean values.
        output_path:
            path to output .html file or output directory where the created Figure will be stored.
            If directory, written file is named 'structural_metrics.html'. Defaults to None.

    Returns:
        metrics_fig:
            structural metrics dashboard for comparison between all the created traces.
    """
    # Get dimensions of dashboard
    nrows = len(trajectory_ids) + 2 # last plot occupies 2 rows
    nmetrics = len(list(total_box_traces.values())[0]) #ncolumns of last row

    # Setup x_axis titles and column titles
    x_axis_titles = {}
    col_titles = []
    for box_trace in total_box_traces[trajectory_ids[0]]:
        # Trace names are f'{trajectory_id}_{metric_name}'
        metricname = '_'.join(box_trace['name'].split('_')[trajectory_ids[0].count('_')+1:])
        cm_dist_count = 0
        if  metricname == 'rg':
            col_titles.append('Radius of gyration (<i>R<sub>g</sub></i>)')
            x_axis_titles['rg'] = '<i>R<sub>g</sub></i>'
        elif metricname == 'dmax':
            col_titles.append('Maximum distance (<i>D<sub>max</sub></i>)')
            x_axis_titles['dmax'] = '<i>D<sub>max</sub></i>'
        elif metricname == 'eed':
            col_titles.append('End-to-end Distance (<i>D<sub>ee</sub></i>)')
            x_axis_titles['eed'] = '<i>D<sub>ee</sub></i>'
        else:
            cm_dist_count += 1
            col_titles.append(f'{metricname} (<i>D<sub>cm{cm_dist_count}</sub></i>)')
            x_axis_titles[metricname] = f'<i>D<sub>cm{cm_dist_count}</sub></i>'

    # Setup Figure specs
    specs = [[{}] * nmetrics] * (nrows)
    specs[-2] = [{'rowspan': 2}] * nmetrics
    specs[-1] = [None] * nmetrics

    # Create metrics figure
    row_titles = trajectory_ids+['']
    metrics_fig = make_subplots(rows=nrows,
                                cols=nmetrics,
                                column_titles=col_titles,
                                row_titles=row_titles,
                                horizontal_spacing=0.255/nmetrics,
                                vertical_spacing=0.45/nrows,
                                specs=specs)

    # Add the Box traces (all rows except last 2)
    for rownum in range(1,nrows-1):
        colnum = 1
        trajectory_id = trajectory_ids[rownum-1]
        for box_trace in total_box_traces[trajectory_id]:

            # Add trace
            metrics_fig.add_trace(box_trace,
                                  row=rownum,
                                  col=colnum)

            # Update axes
            metrics_fig.update_yaxes(**DEFAULT_AXIS,
                                     showticklabels=False, # Remove trace names from boxplot y axis
                                     row=rownum,
                                     col=colnum)

            metrics_fig.update_xaxes(**DEFAULT_AXIS,
                                     row=rownum,
                                     col=colnum)

            colnum += 1

    # Add the Histogram and Scatter traces for last row
    # Store min and max values for each column to have all x_axis with the same dimensions
    min_max_values = {}
    for trajectory_id in trajectory_ids:
        for colnum in range(1,nmetrics+1): # On last row, go through all columns

            # Add traces
            hist_trace = total_hist_traces[trajectory_id][colnum-1]
            scatter_trace = total_scatter_traces[trajectory_id][colnum-1]
            metrics_fig.add_trace(hist_trace,
                                  row=nrows-1,
                                  col=colnum)
            metrics_fig.add_trace(scatter_trace,
                                  row=nrows-1,
                                  col=colnum)

            # Add mean dashed lines
            mean_value = total_avg_values[trajectory_id][colnum-1]
            mean_value_stderr = total_avg_stderr_values[trajectory_id][colnum-1]

            metrics_fig.add_shape(dict(name=scatter_trace.name,
                                       type='line',
                                       xref='x',
                                       x0=mean_value,
                                       x1=mean_value,
                                       y0=0,
                                       y1=np.interp(mean_value,
                                                    scatter_trace.x,
                                                    scatter_trace.y),
                                       line=dict(dash='dot',
                                                 color=hist_trace.marker.color,
                                                 width=4)),
                                       legend='legend',
                                       row=nrows-1,
                                       col=colnum)

            # Allows for hovering the dashed line to get mean value
            metrics_fig.add_trace(go.Scatter(x=[mean_value],
                                             y=[y for y in
                                                np.arange(0,
                                                          np.interp(mean_value,
                                                                    scatter_trace.x,
                                                                    scatter_trace.y) + 0.001,
                                                          0.001)],
                                             mode='markers',
                                             marker_color=hist_trace.marker.color,
                                             hovertext=(f'Avg: {round(mean_value,2)} &plusmn; '
                                                        f'{round(mean_value_stderr,2)}'),
                                             hoverinfo='text',
                                             hoverlabel_bgcolor=hist_trace.marker.color,
                                             fill='toself',
                                             name=(f'Avg: {round(mean_value,2)} &plusmn; '
                                                   f'{round(mean_value_stderr,2)}'),
                                             opacity=0,
                                             showlegend=False),
                                  row=nrows-1,
                                  col=colnum)

            # Store min and max values
            # Set starting comparison values if first trace in column
            try:
                min_max_values[colnum]
            except KeyError:
                min_max_values[colnum] = (math.inf,0)

            # Set max value if greater
            if max(scatter_trace.x) > min_max_values[colnum][1]:
                max_val = max(scatter_trace.x)
            else:
                max_val = min_max_values[colnum][1]

            # Set min value if lesser
            if min(scatter_trace.x) < min_max_values[colnum][0]:
                min_val = min(scatter_trace.x)
            else:
                min_val = min_max_values[colnum][0]

            # Update min and max values in dict
            min_max_values[colnum] = (min_val,max_val)

            # Update axes
            metric_id = scatter_trace['name'].split('_')[trajectory_id.count('_')+1:] # e.g. rg
            metric_name = x_axis_titles['_'.join(metric_id)]

            metrics_fig.update_xaxes(**DEFAULT_AXIS,
                                     title_text=f'{metric_name} (&#197;)', # angstrom symbol
                                     row=nrows-1,
                                     col=colnum)

            metrics_fig.update_yaxes(**DEFAULT_AXIS,
                                     title_text=f'KDE ({metric_name})',
                                     rangemode='tozero',
                                     row=nrows-1,
                                     col=colnum)

    # Equalize x axis range
    for rownum in range(1,nrows):
        for colnum, (min_val,max_val) in min_max_values.items():
            if min_val-round(max_val*0.05) > 0:
                x_axis_min = min_val-round(max_val*0.05)
            else:
                x_axis_min = 0
            metrics_fig.update_xaxes(range=[x_axis_min,max_val + round(max_val*0.05)],
                                     row=rownum,
                                     col=colnum)

    # Make note of which traces are Histograms for the button
    hst_idxs = []
    for i,trace in enumerate(metrics_fig.data):
        if isinstance(trace,go.Histogram):
            hst_idxs.append(i)

    # Save shapes list for button
    mean_line_shapes = metrics_fig.layout.shapes

    # Create buttons
    toggle_histograms_button = dict(type='buttons',
                                    buttons=[dict(method='restyle',
                                                  label='Toggle Histograms',
                                                  visible=True,
                                                  args=[{'visible':False,
                                                         'showlegend':False}, # Traces
                                                        {}, # Layout
                                                        hst_idxs], # Target trace indexes
                                                  args2=[{'visible':True,
                                                          'showlegend':False}, # Traces
                                                         {}, # Layout
                                                         hst_idxs])], # Target trace indexes
                                    showactive=False,
                                    pad=dict(l=0,
                                             r=0,
                                             t=5,
                                             b=5),
                                    bgcolor='#FFFFFF',
                                    font=dict(size=20,
                                              family='Helvetica',
                                              color='black'),
                                    xanchor='left',
                                    x=1.015,
                                    yanchor='bottom',
                                    y=1.01)

    toggle_mean_lines_button = dict(type='buttons',
                                    buttons=[dict(method='relayout',
                                                  label='Toggle Mean Lines',
                                                  visible=True,
                                                  args=['shapes', mean_line_shapes],
                                                  args2=['shapes', [] ])],
                                    showactive=False,
                                    pad=dict(l=0,
                                             r=0,
                                             t=5,
                                             b=5),
                                    bgcolor='#FFFFFF',
                                    font=dict(size=20,
                                              family='Helvetica',
                                              color='black'),
                                    xanchor='left',
                                    x=1.015,
                                    yanchor='top',
                                    y=1)
    # Update Figure Layout
    metrics_fig.update_layout(**DEFAULT_LAYOUT)

    metrics_fig.update_layout(height=240 * nrows,
                              width=620 * nmetrics,
                              legend=dict(font_size=20,
                                          title='KDE Plots',
                                          yanchor='bottom',
                                          y=0),
                              updatemenus=[toggle_histograms_button,toggle_mean_lines_button],
                              margin=dict(t=80,
                                          l=80,
                                          r=0,
                                          b=80),
                                          )

    metrics_fig.update_xaxes(title_standoff=30)

    metrics_fig.update_yaxes(ticks='',
                             title_standoff=15)

    # Fix yaxis title moving to the left when we try to
    # showticklabels=False by making its ticklabels invisible
    # see https://github.com/plotly/plotly.js/issues/6552
    metrics_fig.update_yaxes(tickfont=dict(color='rgba(0,0,0,0)',
                                           size=1))


    # Update column titles position
    for annotation in metrics_fig.layout.annotations[:nmetrics]:
        annotation['yshift'] = 10

    # Update column titles size
    metrics_fig.update_annotations(font_size=34)

    # Reduce size of row titles in subplots
    metrics_fig.for_each_annotation(lambda a: a.update(font_size=20,
                                                       x=0.9825)
                                    if a.text in row_titles
                                    else ())

    # Save Structural Metrics figure
    if output_path is not None:
        if os.path.isdir(output_path):
            metrics_fig.write_html(os.path.join(output_path,'structural_metrics.html'),
                                   config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])
        elif output_path.endswith('.html'):
            metrics_fig.write_html(output_path,
                                   config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])
        else:
            print(('Structural Metrics dashboard was not saved to disk, '
                   'output path must be a directory or .html filepath!'))

    return metrics_fig

def create_analysis_figures(
    analysis_data: dict[str,list[pd.DataFrame]] | None,
    topologies: list[str],
    trajectory_ids: list[str],
    output_directory: str | None = os.getcwd(),
    color_palette: list[str] | None = None,
    ) -> dict[str,list[go.Figure]]:
    """Create interactive figures given analysis data for one or more pairs of trajectory,topology
      files.

    Args:
        analysis_data:
            mapping of data identifiers to lists of DataFrames with the calculated analysis data,
            one element for each given trajectory,topology,trajectory_id trio.
        topologies:
            list of paths to .pdb topology files.
        trajectory_ids:
            prefix trajectory identifiers to distinguish between calculated data files.
        output_directory:
            path to directory where created Figures will be stored. Defaults to current
            working directory.
        color_palette:
            list of color hexcodes, to associate one with each trajectory when creating the
            Structural Metrics interactive figure.
    Returns:
        figures:
            mapping of figure identifiers to lists of the created Figures, one for each trajectory
            outlined in the given analysis data. For example:
            data = {'ContactMaps' : [ContactMap1,ContactMap2,ContactMap3],
                    'DistanceMatrices' : [DistanceMatrix1,DistanceMatrix2,DistanceMatrix3],
                    'SecondaryStructureFrequencies' : [SSFrequency1,SSFrequency2,SSFrequency3],
                    'StructuralMetrics' : [StructuralMetrics1,StructuralMetrics2,StructuralMetrics3] }
    """
    # Setup color palette
    if color_palette is None:
        color_palette = DEFAULT_TRACE_COLOR_PALETTE

    # If data is not available check output directory for it
    if analysis_data is None:
        analysis_data = {'ContactMatrices' : [],
                         'DistanceMatrices' : [],
                         'SecondaryStructureFrequencies' : [],
                         'StructuralMetrics' : [] }
        data_ids_2_data = {'contact_matrix.csv': 'ContactMatrices',
                           'distance_matrix.csv': 'DistanceMatrices',
                           'ss_frequency.csv': 'SecondaryStructureFrequencies',
                           'structural_metrics.csv': 'StructuralMetrics'}
        for data_id, data_name in data_ids_2_data.items():
            for trajectory_id in trajectory_ids:
                try:
                    filename = os.path.join(output_directory,f'{trajectory_id}_{data_id}')
                    calculated_data = pd.read_csv(filename,index_col=0)
                    analysis_data[data_name].append(calculated_data)
                except FileNotFoundError:
                    continue
                else:
                    print(f'Found calculated {trajectory_id}_{data_id}')

    # Create figures
    figures = {'ContactMaps' : [],
               'DistanceMatrices' : [],
               'SecondaryStructureFrequencies' : [],
               'StructuralMetrics' : None }

    ## Get maximum of DistanceMatrices colorbar
    if analysis_data['DistanceMatrices']:
        max_data = get_array_extremum(analysis_data['DistanceMatrices'])
        max_colorbar = round_to_nearest_multiple(max_data,5)
    else:
        max_colorbar = None

    for i,(trajectory_id,topology,color) in enumerate(zip(trajectory_ids,topologies,color_palette)):
        print(f'Creating {trajectory_id} analysis figures...')

        try:
            cmatrix = analysis_data['ContactMatrices'][i]
        except (KeyError,IndexError):
            pass
        else:
            cmap_fig_out = os.path.join(output_directory,
                                        f'{trajectory_id}_contact_map.html')

            cmap_fig = create_contact_map_fig(contact_matrix=cmatrix,
                                              trajectory_id=trajectory_id,
                                              topology=topology,
                                              output_path=cmap_fig_out)

            figures['ContactMaps'].append(cmap_fig)

        try:
            dmatrix = analysis_data['DistanceMatrices'][i]
        except (KeyError,IndexError):
            pass
        else:
            dmatrix_fig_out = os.path.join(output_directory,
                                           f'{trajectory_id}_distance_matrix.html')

            dmatrix_fig = create_distance_matrix_fig(distance_matrix=dmatrix,
                                                     trajectory_id=trajectory_id,
                                                     topology=topology,
                                                     output_path=dmatrix_fig_out,
                                                     max_colorbar=max_colorbar)

            figures['DistanceMatrices'].append(dmatrix_fig)

        try:
            ssfreq = analysis_data['SecondaryStructureFrequencies'][i]
        except (KeyError,IndexError):
            pass
        else:
            ssfreq_fig_out = os.path.join(output_directory,
                                          f'{trajectory_id}_ss_frequency.html')

            ssfreq = create_ss_frequency_figure(ss_frequency=ssfreq,
                                                topology=topology,
                                                trajectory_id=trajectory_id,
                                                output_path=ssfreq_fig_out)

            figures['SecondaryStructureFrequencies'].append(ssfreq)

    total_box_traces = {}
    total_hist_traces = {}
    total_scatter_traces = {}
    total_avg_values = {}
    total_avg_stderr_values = {}

    for i,(trajectory_id,topology,color) in enumerate(zip(trajectory_ids,topologies,color_palette)):
        try:
            metrics = analysis_data['StructuralMetrics'][i]
        except (KeyError,IndexError):
            pass
        else:
            box_traces,\
            hist_traces,\
            scatter_traces,\
            avg_values,\
            avg_stderr_values = create_metrics_traces(metrics=metrics,
                                                      trajectory_id=trajectory_id,
                                                      color=color)

            total_box_traces[trajectory_id] = box_traces
            total_hist_traces[trajectory_id] = hist_traces
            total_scatter_traces[trajectory_id] = scatter_traces
            total_avg_values[trajectory_id] = avg_values
            total_avg_stderr_values[trajectory_id] = avg_stderr_values

    if total_box_traces:
        metrics_fig_out = os.path.join(output_directory,
                                       'structural_metrics.html')

        metrics_fig = create_metrics_fig(trajectory_ids=trajectory_ids,
                                         total_box_traces=total_box_traces,
                                         total_hist_traces=total_hist_traces,
                                         total_scatter_traces=total_scatter_traces,
                                         total_avg_values=total_avg_values,
                                         total_avg_stderr_values=total_avg_stderr_values,
                                         output_path=metrics_fig_out)

        figures['StructuralMetrics'] = metrics_fig

    return figures
