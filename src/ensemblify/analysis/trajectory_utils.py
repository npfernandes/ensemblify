"""Auxiliary functions for the analysis module."""

# IMPORTS
## Standard Library Imports
import importlib.resources
import math
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from functools import reduce

## Third Party Imports
import MDAnalysis as mda
import mdtraj
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy
from MDAnalysis.analysis import dihedrals
from plotly.subplots import make_subplots
from tqdm import tqdm

## Local Imports
from ensemblify.analysis.colors import PARULA_COLORSCALE
from ensemblify.analysis.third_party.simple_mdreader import SimpleMDreader
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.utils import extract_pdb_info, get_array_extremum, kde, round_to_nearest_multiple

# FUNCTIONS
def calc_rg(u: mda.Universe) -> float:
    """Calculate the radius of gyration of the current frame.
    
    Args:
        u:
            Universe pointing to the current frame.
    
    Returns:
        rg:
            radius of gyration of the protein in the current frame.
    """
    protein = u.select_atoms('protein')
    rg = protein.radius_of_gyration()
    return rg


def calc_eed(u: mda.Universe) -> float:
    """Calculate the distance from the N to the C terminal in the current frame.
    
    Args:
        u:
            Universe pointing to the current frame.
    
    Returns:
        eed:
            end-to-end distance of the protein in the current frame.

    """
    nterm = u.select_atoms('protein and name N')[0]
    cterm = u.select_atoms('protein and name C')[-1]
    eed = np.linalg.norm(cterm.position - nterm.position)
    return eed


def calc_dmax(u: mda.Universe) -> float:
    """Calculate the maximum of the distances between any two alpha carbons in the current frame.

    Args:
        u:
            Universe pointing to the current frame.
    
    Returns:
        dmax:
            Maximum of the distances between any two alpha carbons of the protein in the current
            frame.

    """
    ca_selection = u.select_atoms('protein and name CA')
    ca_coordinates = ca_selection.positions #expose numpy array of coords
    distance_matrix_pool = scipy.spatial.distance.cdist(ca_coordinates, ca_coordinates)
    maximum_distance_pool = distance_matrix_pool.max()
    dmax = np.linalg.norm(maximum_distance_pool)
    return dmax


def calc_cm_dist(
    u: mda.Universe,
    sel1: str,
    sel2: str,
    ) -> float:
    """Calculate the distance between the center of mass of two atom selections
    in the current frame.

    Args:
        u:
            Universe pointing to the current frame.
        sel1:
            MDAnalysis selection string for selecting an AtomGroup whose center of mass will be
            calculated.
        sel2:
            MDAnalysis selection string for selecting an AtomGroup whose center of mass will be
            calculated.
    
    Returns:
        cm_dist:
            Center of mass distance between AtomGroups selected by sel1 and sel2.

    """
    cm1 = u.select_atoms(sel1).center_of_mass()
    cm2 = u.select_atoms(sel2).center_of_mass()
    cm_dist = np.linalg.norm(cm1 - cm2)
    return cm_dist


def calculate_ramachandran_data(
    trajectory: str,
    topology: str,
    output_path: str | None = os.getcwd(),
    ) -> pd.DataFrame:
    """Calculate a dihedral angles matrix from trajectory and topology files.

    Phi and Psi dihedral angle values are calculated for each residue in each trajectory frame.
    Optionally saves the matrix to output directory in .csv format, defaulting to current working
    directory.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'ramachandran_data.csv'. Defaults to current working directory.

    Returns:
        dihedrals_matrix:
            DataFrame with Phi and Psi values of each residue for each frame of the trajectory.
    """
    # Create Universe
    u = mda.Universe(topology, trajectory)
    protein = u.select_atoms('protein')

    # Calculate dihedral angles of each conformation
    with warnings.catch_warnings():
        # Suppress warning for first and last res not having angles
        warnings.filterwarnings('ignore', category=UserWarning)
        rama = dihedrals.Ramachandran(protein).run()
        rama_angles = rama.results.angles

    # Get our phi psi angle values
    rama_xs = []
    rama_ys = []
    for frame_dihedrals_matrix in rama_angles:
        for phi,psi in frame_dihedrals_matrix:
            rama_xs.append(phi)
            rama_ys.append(psi)

    # Create DataFrame with data
    dihedrals_matrix = pd.DataFrame({'Phi':rama_xs,
                                     'Psi':rama_ys})

    # Save dihedrals matrix
    if os.path.isdir(output_path):
        dihedrals_matrix_output_filename = 'ramachandran_data.csv'
        dihedrals_matrix.to_csv(os.path.join(output_path,dihedrals_matrix_output_filename))
    elif output_path.endswith('.csv'):
        dihedrals_matrix.to_csv(output_path)
    else:
        print(('Ramachandran data was not saved to disk, '
               'output path must be a directory or .csv filepath!'))

    return dihedrals_matrix


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
    rama_fig = go.Figure()

    # Add Ramachandran Reference Contours
    rama_ref_data = np.load(importlib.resources.files('ensemblify.analysis').joinpath('rama_ref_data.npy').open('rb')).flatten()
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
                           plot_bgcolor='#FFFFFF',
                           font=dict(family='Helvetica',
                                     color='black',
                                     size=30),
                           modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                           title=dict(text=rama_title,
                                      xref='paper',
                                      x=0.5),
                           xaxis=dict(title='\u03A6', # Phi
                                      range=[-180,180],
                                      tick0=-180,
                                      dtick=60,
                                      ticks='outside',
                                      ticksuffix='\u00B0',
                                      ticklen=10,
                                      tickwidth=4,
                                      showgrid=False),
                           yaxis=dict(title='\u03A8', # Psi
                                      range=[-180,180],
                                      tick0=-180,
                                      dtick=60,
                                      ticks='outside',
                                      ticksuffix='\u00B0',
                                      ticklen=10,
                                      tickwidth=4,
                                      showgrid=False,
                                      title_standoff=5),
                            legend=dict(orientation='h',
                                        xref='paper',
                                        xanchor='center',
                                        x=0.5,
                                        y=1.065,
                                        font_size=18),
                           shapes=shapes)

    rama_fig.update_xaxes(showline=True,
                          linewidth=4,
                          linecolor='black',
                          color='black',
                          mirror=True)

    rama_fig.update_yaxes(showline=True,
                          linewidth=4,
                          linecolor='black',
                          color='black',
                          mirror=True)

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


def calculate_contact_matrix_frame(
    u: mda.Universe,
    frame_idx: int,
    frame_weight: float,
    ) -> np.ndarray:
    """Calculates a contact matrix for a frame of a trajectory.

    Args:
        u:
            `MDAnalysis.Universe` object containing the trajectory being analyzed.
        frame_idx:
            number of the frame to be analyzed.
        frame_weight:
            contacts found in this frame will be assigned this value in the
            resulting matrix instead of the default value of 1. In a uniformly
            weighted matrix, this value will be of 1 / number of trajectory frames.

    Returns:
        np.ndarray:
            contact matrix for the current frame.
    """
    # Point universe to frame of interest
    u.trajectory[frame_idx]

    # Create results contact matrix
    contact_matrix = np.array([[0.0] * len(u.residues)] *len(u.residues))

    # For each residue, iterate over all other residues
    for res1 in u.residues:

        # Select current residue's atoms
        current_res_atom_selection = res1.atoms

        # Expose coordinates np.array
        current_res_atom_coordinates = current_res_atom_selection.positions

        # Only calculate distances once for each pair, ignoring neighbours
        for res2 in u.residues[res1.resindex + 3:]:

            # Select current residue's atoms
            target_res_atom_selection = res2.atoms

            # Expose coordinates np.array
            target_res_atom_coordinates = target_res_atom_selection.positions

            # Calculate distance matrix
            distance_matrix = scipy.spatial.distance.cdist(current_res_atom_coordinates,
                                                            target_res_atom_coordinates,
                                                            'euclidean')

            if np.argwhere(distance_matrix < 4.5).shape[0] > 0:
                # Add contacts on both halves of matrix
                contact_matrix[res1.resindex,res2.resindex] = 1.0
                contact_matrix[res2.resindex,res1.resindex] = 1.0

    # Reweigh matrix
    contact_matrix *= frame_weight

    return contact_matrix


def calculate_contact_matrix(
    trajectory: str,
    topology: str,
    weights: np.ndarray | None = None,
    output_path: str | None = os.getcwd(),
    ) -> pd.DataFrame:
    """Calculate a contact frequency matrix from a trajectory and topology files.
    
    The contact frequency of a residue pair is calculated from the number of times they are in
    contact over all frames in the trajectory.
    Optionally saves the matrix to output directory in .csv format.
    Uses multiprocessing whenever possible.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        weights:
            array of weights to be used when calculating the contact matrix. If None, uniform
            weights are used.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'contact_matrix.csv'. Defaults to current working directory.

    Returns:
        contact_matrix:
            DataFrame with the frequency of each residue contact in the trajectory.
    """
    # Setup Universe object
    u = mda.Universe(topology,trajectory)
    trajectory_size = len(u.trajectory)

    # Setup multiprocessing variables
    frame_idxs = np.array(range(trajectory_size))
    universes = [u] * trajectory_size

    if weights is None:
        # Setup multiprocessing variables
        weights = np.array([1/trajectory_size] * trajectory_size )

        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            contact_matrix_array = reduce(lambda x,y: np.add(x,y),
                                          tqdm(ppe.map(calculate_contact_matrix_frame,
                                                       universes,
                                                       frame_idxs,
                                                       weights),
                                               desc='Calculating contact matrix...',
                                               total=trajectory_size))
    else:
        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            contact_matrix_array = reduce(lambda x,y: np.add(x,y),
                                          tqdm(ppe.map(calculate_contact_matrix_frame,
                                                       universes,
                                                       frame_idxs,
                                                       weights),
                                               desc='Calculating reweighted contact matrix...',
                                               total=trajectory_size))

    # Convert calculated averaged matrix to DataFrame
    contact_matrix = pd.DataFrame(data=contact_matrix_array,
                                  index=list(range(1,contact_matrix_array.shape[0]+1)),
                                  columns=[str(x) for x in range(1,contact_matrix_array.shape[0]+1)])

    # Save contact matrix
    if os.path.isdir(output_path):
        if weights is None:
            contact_matrix_output_filename = 'contact_matrix_csv'
        else:
            contact_matrix_output_filename = 'contact_matrix_reweighted.csv'
        contact_matrix.to_csv(os.path.join(output_path,contact_matrix_output_filename))
    elif output_path.endswith('.csv'):
        contact_matrix.to_csv(output_path)
    else:
        print(('Contact matrix was not saved to disk, '
               'output path must be a directory or .csv filepath!'))

    return contact_matrix


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

    # Create Contact Map Figure
    cmap_fig = go.Figure()

    # Add our data
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

    # Setup chain dividers lines
    num_res = len(contact_matrix.columns)
    chain_ends = [] # to be used in tickvals
    chain_begins = [] # to be used in tickvals
    shapes = []
    cumulative_residues = 0

    for chain_letter in chain_letters[:-1]:
        chain_begins.append(cumulative_residues)
        chain_size = len(resranges[chain_letter])
        chain_end = cumulative_residues + chain_size
        chain_ends.append(chain_end)

        shapes.append(dict(type='line',
                           xref='x',
                           x0=chain_end-0.5,
                           x1=chain_end-0.5,
                           y0=0-0.5,
                           y1=num_res+0.5,
                           line=dict(color='black',
                                     width=2)))
        shapes.append(dict(type='line',
                           yref='y',
                           y0=chain_end-0.5,
                           y1=chain_end-0.5,
                           x0=0-0.5,
                           x1=num_res+0.5,
                           line=dict(color='black',
                                     width=2)))
        cumulative_residues += chain_size
    chain_begins.append(num_res - len(resranges[chain_letters[-1]]))
    chain_ends.append(num_res)

    # Setup tick values
    tickvals = []
    curr_chain = 0
    curr_val = chain_begins[curr_chain]
    tick_step = num_res // len(chain_letters) // 4 # always 5 ticks per axis

    while curr_val <= num_res:
        try:
            chain_end = chain_ends[curr_chain]
        except IndexError:
            tickvals.append(curr_val)
        else:
            if chain_end - curr_val <= tick_step:
                curr_chain += 1
                try:
                    curr_val = chain_begins[curr_chain]
                    tickvals.append(curr_val)
                except IndexError:
                    if chain_ends[-1] - curr_val <= 3:
                        tickvals.append(chain_ends[-1] - 1)
                    else:
                        tickvals.append(curr_val)
                        tickvals.append(chain_ends[-1])
            else:
                tickvals.append(curr_val)
        curr_val += tick_step

    # Update Figure Layout
    if difference:
        if trajectory_id is not None:
            cmap_title = f'{trajectory_id} Difference Contact Map'
        else:
            cmap_title = 'Difference Contact Map'
    elif reweighted:
        if trajectory_id is not None:
            cmap_title = f'{trajectory_id} Reweighted Contact Map'
        else:
            cmap_title = 'Reweighted Contact Map'
    else:
        if trajectory_id is not None:
            cmap_title = f'{trajectory_id} Contact Map'
        else:
            cmap_title = 'Contact Map'

    cmap_fig.update_layout(width=900,
                           height=900,
                           plot_bgcolor='#FFFFFF',
                           font=dict(family='Helvetica',
                                     color='black',
                                     size=30),
                           modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                           title=dict(text=cmap_title,
                                      x=0.5,
                                      subtitle=dict(text=('Frequency of contacts between '
                                                          'any atoms of a residue pair'),
                                                    font=dict(color='gray',
                                                              size=24,))),
                           margin_t=150, # to fit subtitle
                           xaxis=dict(title='Residue',
                                      tickvals=tickvals,
                                      ticks='outside',
                                      ticklen=10,
                                      tickwidth=4,
                                      showgrid=False,
                                      constrain='domain',),
                           yaxis=dict(title='Residue',
                                      tickvals=tickvals,
                                      ticks='outside',
                                      ticklen=10,
                                      tickwidth=4,
                                      showgrid=False,
                                      title_standoff=5,
                                      scaleanchor='x',
                                      constrain='domain',),
                           shapes=shapes)

    cmap_fig.update_xaxes(showline=True,
                          linewidth=4,
                          linecolor='black',
                          color='black',
                          mirror=True)

    cmap_fig.update_yaxes(showline=True,
                          linewidth=4,
                          linecolor='black',
                          color='black',
                          mirror=True)

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


def calculate_distance_matrix_frame(
    u: mda.Universe,
    frame_idx: int,
    frame_weight: float,
    ) -> np.ndarray:
    """Calculates a distance matrix for the alpha carbons of a trajectory frame.

    Args:
        u:
            `MDAnalysis.Universe` object containing the trajectory being analyzed.
        frame_idx:
            number of the frame to be analyzed.
        frame_weight:
            distances calculated for this frame will be multiplied by this value
            in the resulting frame matrix. In a uniformly weighted matrix, calculated
            distances will be multiplied by 1 / number of trajectory frames.

    Returns:
            np.ndarray: distance matrix for the current frame.
    """
    # Point universe to frame of interest
    u.trajectory[frame_idx]

    # Select alpha carbons
    ca_selection = u.select_atoms('protein and name CA')

    # Expose coordinates np.array
    ca_coordinates = ca_selection.positions

    # Calculate distance matrix
    distance_matrix = scipy.spatial.distance.cdist(ca_coordinates,ca_coordinates,'euclidean')

    # Ignore neighbours
    for ca1_idx, ca2_idx in np.argwhere(distance_matrix):
        if abs(ca1_idx - ca2_idx) <= 2:
            distance_matrix[ca1_idx,ca2_idx] = 0.0

    # Reweigh matrix
    distance_matrix *= frame_weight

    return distance_matrix


def calculate_distance_matrix(
    trajectory: str,
    topology: str,
    weights: np.ndarray | None = None,
    output_path: str | None = os.getcwd(),
    ) -> pd.DataFrame:
    """Calculate an alpha carbon average distance matrix from a trajectory and topology files.
    
    The distances between different pairs of alpha carbons pair is calculated for each trajectory
    frame and the values are then averaged to create the final distance matrix. 
    
    Optionally save the matrix to output directory in .csv format.
    Uses multiprocessing whenever possible.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        weights:
            array of weights to be used when calculating the distance matrix. If None, uniform
            weights are used.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'distance_matrix.csv'. Defaults to current working directory.

    Returns:
        distance_matrix:
            DataFrame with the average distance between each pair of alpha carbons in the
            trajectory.
    """
    # Setup Universe object
    u = mda.Universe(topology,trajectory)
    trajectory_size = len(u.trajectory)

    # Setup multiprocessing variables
    frame_idxs = np.array(range(trajectory_size))
    universes = [u] * trajectory_size

    if weights is None:
        weights = np.array([1/trajectory_size] * trajectory_size )

        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            distance_matrix_array = reduce(lambda x,y: np.add(x,y),
                                           tqdm(ppe.map(calculate_distance_matrix_frame,
                                                        universes,
                                                        frame_idxs,
                                                        weights),
                                                desc='Calculating distance matrix... ',
                                                total=trajectory_size))
    else:
        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            distance_matrix_array = reduce(lambda x,y: np.add(x,y),
                                           tqdm(ppe.map(calculate_distance_matrix_frame,
                                                        universes,
                                                        frame_idxs,
                                                        weights),
                                                desc='Calculating reweighted distance matrix... ',
                                                total=trajectory_size))

    # Convert calculated averaged matrix to DataFrame
    distance_matrix = pd.DataFrame(data=distance_matrix_array,
                                   index=list(range(1,distance_matrix_array.shape[0]+1)),
                                   columns=[str(x) for x in range(1,distance_matrix_array.shape[0]+1)])
    
    # Save distance matrix
    if os.path.isdir(output_path):
        if weights is None:
            distance_matrix_output_filename = 'distance_matrix.csv'
        else:
            distance_matrix_output_filename = 'distance_matrix_reweighted.csv'
        distance_matrix.to_csv(os.path.join(output_path,distance_matrix_output_filename))
    elif output_path.endswith('.csv'):
        distance_matrix.to_csv(output_path)
    else:
        print(('Distance matrix was not saved to disk, '
               'output path must be a directory or .csv filepath!'))

    return distance_matrix


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

    # Create Distance Matrix Figure
    dmatrix_fig = go.Figure()

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

    # Setup chain dividers lines
    num_res = len(distance_matrix.columns)
    chain_ends = [] # to be used in tickvals
    chain_begins = [] # to be used in tickvals
    shapes = []
    cumulative_residues = 0

    for chain_letter in chain_letters[:-1]:
        chain_begins.append(cumulative_residues)
        chain_size = len(resranges[chain_letter])
        chain_end = cumulative_residues + chain_size
        chain_ends.append(chain_end)

        shapes.append(dict(type='line',
                           xref='x',
                           x0=chain_end-0.5,
                           x1=chain_end-0.5,
                           y0=0-0.5,
                           y1=num_res+0.5,
                           line=dict(color='black',
                                     width=2)))
        shapes.append(dict(type='line',
                           yref='y',
                           y0=chain_end-0.5,
                           y1=chain_end-0.5,
                           x0=0-0.5,
                           x1=num_res+0.5,
                           line=dict(color='black',
                                     width=2)))
        cumulative_residues += chain_size
    chain_begins.append(num_res - len(resranges[chain_letters[-1]]))
    chain_ends.append(num_res)

    # Setup tick values
    tickvals = []
    curr_chain = 0
    curr_val = chain_begins[curr_chain]
    tick_step = num_res // len(chain_letters) // 4 # always 5 ticks per axis

    while curr_val <= num_res:
        try:
            chain_end = chain_ends[curr_chain]
        except IndexError:
            tickvals.append(curr_val)
        else:
            if chain_end - curr_val <= tick_step:
                curr_chain += 1
                try:
                    curr_val = chain_begins[curr_chain]
                    tickvals.append(curr_val)
                except IndexError:
                    if chain_ends[-1] - curr_val <= 3:
                        tickvals.append(chain_ends[-1] - 1)
                    else:
                        tickvals.append(curr_val)
                        tickvals.append(chain_ends[-1])
            else:
                tickvals.append(curr_val)
        curr_val += tick_step

    # Update Figure Layout
    if difference:
        if trajectory_id is not None:
            dmatrix_title = f'{trajectory_id} Difference Distance Matrix'
        else:
            dmatrix_title = 'Difference Distance Matrix'
    elif reweighted:
        if trajectory_id is not None:
            dmatrix_title = f'{trajectory_id} Reweighted Distance Matrix'
        else:
            dmatrix_title = 'Reweighted Distance Matrix'
    else:
        if trajectory_id is not None:
            dmatrix_title = f'{trajectory_id} Distance Matrix'
        else:
            dmatrix_title = 'Distance Matrix'

    dmatrix_fig.update_layout(width=900,
                              height=900,
                              plot_bgcolor='#FFFFFF',
                              font=dict(family='Helvetica',
                                        color='black',
                                        size=30),
                              modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                              title=dict(text=dmatrix_title,
                                         x=0.5,
                                         subtitle=dict(text=('Average distance between alpha '
                                                             'carbons of a residue pair'),
                                                       font=dict(color='gray',
                                                                 size=24,))),
                              margin_t=150, # to fit subtitle
                              xaxis=dict(title='Residue',
                                         tickvals=tickvals,
                                         ticks='outside',
                                         ticklen=10,
                                         tickwidth=4,
                                         showgrid=False,
                                         constrain='domain',),
                              yaxis=dict(title='Residue',
                                         tickvals=tickvals,
                                         ticks='outside',
                                         ticklen=10,
                                         tickwidth=4,
                                         showgrid=False,
                                         title_standoff=5,
                                         scaleanchor='x',
                                         constrain='domain',),
                              shapes=shapes)

    dmatrix_fig.update_xaxes(showline=True,
                             linewidth=4,
                             linecolor='black',
                             color='black',
                             mirror=True)

    dmatrix_fig.update_yaxes(showline=True,
                             linewidth=4,
                             linecolor='black',
                             color='black',
                             mirror=True)

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


def calculate_ss_assignment(
    trajectory: str,
    topology: str,
    output_path: str | None = None,
    ) -> pd.DataFrame:
    """Calculate a secondary structure assignment matrix from a trajectory and topology files.
    
    For each residue in each frame of the trajectory, calculate it's secondary structure
    assignment using DSSP. The simplified DSSP codes used here are:
        'H' : Helix. Either of the 'H', 'G', or 'I' codes.
        'E' : Strand. Either of the 'E', or 'B' codes.
        'C' : Coil. Either of the 'T', 'S' or ' ' codes.
    Optionally save the resulting matrix to output directory in .csv format.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'ss_assignment.csv'. Defaults to None, and no file is written.

    Returns:
        ss_assign:
            DataFrame holding the secondary structure assignment matrix.
    """
    # Load trajectory
    traj_md = mdtraj.load(trajectory,top=topology)

    # Calculate DSSP
    dssp_mdt = mdtraj.compute_dssp(traj_md, simplified=True)
    ss_assign = pd.DataFrame(dssp_mdt)

    # Rename columns to correct residue numbering
    top_info = extract_pdb_info(topology)
    resranges = {}
    chain_letters = []
    for chain_number in range(len(top_info.keys()),0,-1):
        chain_letter, starting_res, chain_size = top_info[chain_number]
        resranges[chain_letter] = [ x for x in range(starting_res, starting_res + chain_size)]
        chain_letters.append(chain_letter)

    if len(chain_letters) > 1:
        full_column_names = [f'{chain_letter}{resnum}' for chain_letter in chain_letters
                    for resnum in resranges[chain_letter]]
    else:
        full_column_names = [f'{resnum}' for chain_letter in chain_letters
                    for resnum in resranges[chain_letter]]

    ss_assign.columns = full_column_names

    # Save ss assignment
    if output_path is not None:
        if os.path.isdir(output_path):
            ss_assign.to_csv(os.path.join(output_path,'ss_assignment.csv'))
        elif output_path.endswith('.csv'):
            ss_assign.to_csv(output_path)
        else:
            print(('Secondary structure assignment matrix was not saved to disk, '
                   'output path must be a directory or .csv filepath!'))

    return ss_assign


def calculate_ss_frequency(
    trajectory: str,
    topology: str,
    weights: np.ndarray | None = None,
    output_path: str | None = os.getcwd(),
    ) -> pd.DataFrame:
    """Calculate secondary structure assignment frequencies from a trajectory and topology files.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        weights:
            optional array of weight values to be used in secondary structure
            assignment reweighting. If None, uniform weights are used.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'ss_frequency.csv'. Defaults to current working directory.

    Returns:
        frequency:
            secondary structure frequencies matrix for trajectory being analyzed.
    """
    # Calculate ss assignment
    ss_assignment = calculate_ss_assignment(trajectory=trajectory,
                                            topology=topology)

    if weights is None:
        # Count the frequency of each secondary structure element
        frequency = ss_assignment.apply(lambda x: pd.Series(x).value_counts())
        frequency = frequency.fillna(0)
        frequency = frequency / ss_assignment.shape[0]
    else:
        print('Calculating reweighted secondary structure assignment frequency matrix...')
        # Iterate over each column and compute reweighted frequency
        reweighted_freqs = {'C': [],
                            'E': [],
                            'H': []}
        for label in ss_assignment:
            c_weighted_sum = ((ss_assignment[label] == 'C') * weights).sum()
            e_weighted_sum = ((ss_assignment[label] == 'E') * weights).sum()
            h_weighted_sum = ((ss_assignment[label] == 'H') * weights).sum()

            reweighted_freqs['C'].append(c_weighted_sum)
            reweighted_freqs['E'].append(e_weighted_sum)
            reweighted_freqs['H'].append(h_weighted_sum)

        # Get frequency DataFrame
        frequency = pd.DataFrame(reweighted_freqs,
                                 index=ss_assignment.columns).T

    # Save ss assignment frequency
    if os.path.isdir(output_path):
        if weights is None:
            ss_frequency_output_filename = 'ss_frequency.csv'
        else:
            ss_frequency_output_filename = 'ss_frequency_reweighted.csv'
        frequency.to_csv(os.path.join(output_path,ss_frequency_output_filename))
    elif output_path.endswith('.csv'):
        frequency.to_csv(output_path)
    else:
        print(('Secondary structure assignment frequency matrix was not saved to disk, '
               'output path must be a directory or .csv filepath!'))

    return frequency


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
    else:
        x_labels = [f'{resnum}' for chain_letter in chain_letters
                    for resnum in resranges[chain_letter]]

    # Create Figure
    ss_freq_fig = go.Figure()

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

    # Setup chain dividers lines
    num_res = len(ss_frequency.columns)
    chain_ends = [] # to be used in tickvals
    chain_begins = [] # to be used in tickvals
    shapes = []
    cumulative_residues = 0

    for chain_letter in chain_letters[:-1]:
        chain_begins.append(cumulative_residues+1)
        chain_size = len(resranges[chain_letter])
        chain_end = cumulative_residues + chain_size
        chain_ends.append(chain_end)

        shapes.append(dict(type='line',
                           xref='x',
                           x0=chain_end,
                           x1=chain_end,
                           y0=0,
                           y1=1,
                           line=dict(color='black',
                                     width=2)))

        cumulative_residues += chain_size
    chain_begins.append(num_res - len(resranges[chain_letters[-1]]) + 1)
    chain_ends.append(num_res)

    # Setup axis tick values
    tickvals = []
    chain_counter = 0
    curr_val = chain_begins[chain_counter]
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
                    curr_val = chain_begins[chain_counter]
                    tickvals.append(curr_val)
                except IndexError:
                    if chain_ends[-1] - curr_val <= 3:
                        tickvals.append(chain_ends[-1])
                    else:
                        tickvals.append(curr_val)
                        tickvals.append(chain_ends[-1])
            else:
                tickvals.append(curr_val)
        curr_val += tick_step

    # Setup tick text
    x_text = []
    for x_val in tickvals:
        x_t = x_labels[x_val-1]
        x_text.append(x_t)

    # Update Figure Layout
    if difference:
        if trajectory_id is not None:
            ss_freq_title = f'{trajectory_id} Difference Sec. Struct. Frequencies'
        else:
            ss_freq_title = 'Difference Secondary Structure Frequencies'
    elif reweighted:
        if trajectory_id is not None:
            ss_freq_title = f'{trajectory_id} Reweighted Sec. Struct. Frequencies'
        else:
            ss_freq_title = 'Reweighted Secondary Structure Frequencies'
    else:
        if trajectory_id is not None:
            ss_freq_title = f'{trajectory_id} Secondary Structure Frequencies'
        else:
            ss_freq_title = 'Secondary Structure Frequencies'

    if difference:
        range_y = [-1,1]
    else:
        range_y = [0,1]

    ss_freq_fig.update_layout(width=1000,
                              height=750,
                              font=dict(family='Helvetica',
                                        color='black',
                                        size=30),
                              plot_bgcolor='#FFFFFF',
                              paper_bgcolor='#FFFFFF',
                              modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                              title=dict(text=ss_freq_title,
                                         x=0.5,
                                         subtitle=dict(text=('Frequency of each secondary '
                                                             'structure assignment code for '
                                                             'each residue '),
                                                       font=dict(color='gray',
                                                                 size=24,))),
                              margin_t=150, # to fit subtitle
                              xaxis=dict(title='Residue',
                                         ticks='outside',
                                         tickvals=tickvals,
                                         ticktext=x_text,
                                         ticklen=10,
                                         tickwidth=4,
                                         showgrid=False),
                              yaxis=dict(title='Frequency',
                                         ticks='outside',
                                         ticklen=10,
                                         tickwidth=4,
                                         range=range_y,
                                         showgrid=False),
                              shapes=shapes)

    ss_freq_fig.update_xaxes(showline=True,
                             linewidth=4,
                             linecolor='black',
                             color='black',
                             mirror=True)

    ss_freq_fig.update_yaxes(showline=True,
                             linewidth=4,
                             linecolor='black',
                             color='black',
                             mirror=True)

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


def calculate_metrics_data(
    trajectory: str,
    topology: str,
    rg: bool | None = True,
    dmax: bool | None = True,
    eed: bool | None = True,
    cm_dist: dict[str,tuple[str,str]] | None = None,
    output_path: str | None = os.getcwd(),
    ) -> pd.DataFrame:
    """Calculate structural metrics for each frame of a trajectory.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        rg:
            whether to calculate the radius of gyration of the protein.
        dmax:
            whether to calculate the maximum distance between any two alpha carbons in the protein.
        eed:
            whether to calculate the distance from the N to C terminal of the protein.
        cm_dist:
            mapping of identifiers to tuples with two selection strings for creating MDAnalysis
            AtomGroups, whose center mass distance will be calculated. For example:
                {'inter_domain' : ('resid 1:30', 'resid 110:140')}
            If None, no center mass distances are calculated.
            See https://userguide.mdanalysis.org/stable/selections.html for more information about
            MDAnalysis selections.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'structural_metrics.csv'. Defaults to current working directory.

    Returns:
        traj_analysis:
            DataFrame where columns are the desired structural metrics and rows are the frames
            of the trajectory.
    """
    # Initialize MDReader
    u = SimpleMDreader(trajectory=trajectory,
                       topology=topology)

    # Calculate trajectory metrics
    results = []
    if rg:
        print('Calculating rg...')
        rgs = np.array(u.do_in_parallel(calc_rg,u))
        results.append(('rg',rgs))
    if eed:
        print('Calculating eed...')
        eeds = np.array(u.do_in_parallel(calc_eed,u))
        results.append(('eed',eeds))
    if dmax:
        print('Calculating dmax...')
        dmaxs = np.array(u.do_in_parallel(calc_dmax,u))
        results.append(('dmax',dmaxs))
    if cm_dist:
        for cm_dist_id in cm_dist:
            print(f'Calculating {cm_dist_id}...')
            cm_dists = np.array(u.do_in_parallel(calc_cm_dist,
                                                 u,
                                                 cm_dist[cm_dist_id][0],
                                                 cm_dist[cm_dist_id][1]))
            results.append((cm_dist_id,cm_dists))

    # Extract column names and values
    column_ids = []
    values = []
    for metric_id,metric_values in results:
        column_ids.append(metric_id)
        values.append(metric_values)

    # Create trajectory analysis DataFrame
    metrics_array = np.dstack(tuple(values))
    traj_analysis = pd.DataFrame(metrics_array.reshape(-1,metrics_array.shape[-1]),
                                 columns=column_ids)

    if output_path is not None:
        # Save structural metrics
        if os.path.isdir(output_path):
            traj_analysis.to_csv(os.path.join(output_path,'structural_metrics.csv'))
        elif output_path.endswith('.csv'):
            traj_analysis.to_csv(output_path)
        else:
            print(('Structural metrics DataFrame was not saved to disk, '
                   'output path must be a directory or .csv filepath!'))

    return traj_analysis


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
            metrics_fig.update_yaxes(showticklabels=False, # Remove trace names from boxplot y axis
                                     row=rownum,
                                     col=colnum)
            metrics_fig.update_xaxes(ticks='outside',
                                     row=rownum,
                                     col=colnum)

            colnum += 1

    # Store min and max values for each column to have all x_axis with the same dimensions
    min_max_values = {}

    # Add the Histogram and Scatter traces for last row
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

            metrics_fig.update_xaxes(ticks='outside',
                                     title_text=f'{metric_name} (&#197;)', # angstrom symbol
                                     row=nrows-1,
                                     col=colnum)
            metrics_fig.update_yaxes(ticks='outside',
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

    # Update Column titles position
    for annotation in metrics_fig.layout.annotations[:nmetrics]:
        annotation['yshift'] = 10

    # Update column titles size
    metrics_fig.update_annotations(font_size=34)

    # Update Figure Layout
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
                                    font_size=20,
                                    font_family='Helvetica',
                                    font_color='black',
                                    xanchor='left',
                                    x=1.02,
                                    yanchor='bottom',
                                    y=1.01)

    toggle_mean_lines_button =  dict(type='buttons',
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
                                     font_size=20,
                                     font_family='Helvetica',
                                     font_color='black',
                                     xanchor='left',
                                     x=1.02,
                                     yanchor='top',
                                     y=1)


    metrics_fig.update_layout(height=120 + (225 * nrows) + 0,
                              width=140 + (580 * nmetrics) + 220,
                              legend=dict(title='KDE Plots',
                                          yanchor='bottom',
                                          y=0),
                              font=dict(family='Helvetica',
                                        color='black',
                                        size=30),
                              plot_bgcolor='#FFFFFF',
                              paper_bgcolor='#FFFFFF',
                              updatemenus=[toggle_histograms_button, toggle_mean_lines_button],
                              modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                              margin=dict(t=120,
                                          l=140,
                                          r=220,
                                          b=0))

    metrics_fig.update_layout(legend_font_size=20)

    metrics_fig.update_xaxes(showline=True,
                             ticklen=10,
                             tickwidth=4,
                             linewidth=4,
                             linecolor='black',
                             color='black',
                             mirror=True,
                             title_standoff=30)

    metrics_fig.update_yaxes(showline=True,
                             ticks='',
                             linewidth=4,
                             linecolor='black',
                             color='black',
                             mirror=True,
                             title_standoff=15)

    # Fix yaxis title moving to the left when we try to
    # showticklabels=False by making its ticklabels invisible
    # see https://github.com/plotly/plotly.js/issues/6552
    metrics_fig.update_yaxes(tickfont=dict(color='rgba(0,0,0,0)',
                                           size=1))

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


def calculate_analysis_data(
    trajectories: list[str],
    topologies: list[str],
    trajectory_ids: list[str],
    output_directory: str | None = os.getcwd(),
    ramachandran_data: bool = True,
    distancematrices: bool = True,
    contactmatrices: bool = True,
    ssfrequencies: bool = True,
    rg: bool = True,
    dmax: bool = True,
    eed: bool = True,
    cm_dist: dict[str,tuple[str,str]] | None = None,
    ) -> dict[str,list[pd.DataFrame]]:
    """Calculate  structural data for each given pair of trajectory,topology files.

    Args:
        trajectories:
            list of paths to .xtc trajectory files.
        topologies:
            list of paths to .pdb topology files.
        trajectory_ids:
            prefix trajectory identifiers to distinguish between calculated data files.
        output_directory:
            path to directory where calculated data will be stored. Defaults to current
            working directory.
        ramachandran_data:
            whether to calculate a dihedral angles matrix for each trajectory,topology 
            file pair.
        distancematrices:
            whether to calculate an alpha carbon distance matrix for each trajectory,topology
            file pair.
        contactmatrices:
            whether to calculate a contact frequency matrix for each trajectory,topology
            file pair.
        ssfrequencies:
            whether to calculate a secondary structure assignment frequency matrix for each
            trajectory, topology file pair.
        rg:
            whether to calculate the radius of gyration for each trajectory,topology file pair.
        dmax:
            whether to calculate the maximum distance between any two alpha carbons for each
            trajectory,topology file pair.
        eed:
            whether to calculate the distance between the N and C terminal for each trajectory,
            topology file pair.
        cm_dist:
            mapping of identifiers to tuples with two selection strings for creating MDAnalysis
            AtomGroups, whose center mass distance will be calculated. If None, no center mass
            distances are calculated. See https://userguide.mdanalysis.org/stable/selections.html
            for more information about MDAnalysis selections. For example:
            {'inter_domain' : ('resid 1:30', 'resid 110:140')}

    Returns:
        data:
            mapping of data identifiers to lists of DataFrames with the calculated analysis data,
            one element for each given trajectory,topology,trajectory_id trio. For example:
            data = {'DistanceMatrices' : [DistanceMatrix1,DistanceMatrix2,DistanceMatrix3],
                    'ContactMatrices' : [ContactMatrix1,ContactMatrix2,ContactMatrix3],
                    'SecondaryStructureFrequencies' : [SSFrequency1,SSFrequency2,SSFrequency3],
                    'StructuralMetrics' : [StructuralMetrics1,StructuralMetrics2,
                                           StructuralMetrics3]}
    """
    # Calculate analysis data
    data = {'DistanceMatrices' : [],
            'ContactMatrices' : [],
            'SecondaryStructureFrequencies' : [],
            'StructuralMetrics' : [] }

    for trajectory_id,trajectory,topology in zip(trajectory_ids,trajectories,topologies):
        print(f'Analyzing {trajectory_id} trajectory...')

        # Analysis not meant for interactive figures
        if ramachandran_data:
            print(f'Calculating ramachandran data for {trajectory_id}...')
            rama_data_out = os.path.join(output_directory,
                                         f'{trajectory_id}_ramachandran_data.csv')

            calculate_ramachandran_data(trajectory=trajectory,
                                        topology=topology,
                                        output_path=rama_data_out)

        # Analysis meant for interactive figures
        if contactmatrices:
            print(f'Calculating contact matrix for {trajectory_id}...')
            cmatrix_out = os.path.join(output_directory,
                                       f'{trajectory_id}_contact_matrix.csv')

            cmatrix = calculate_contact_matrix(trajectory=trajectory,
                                               topology=topology,
                                               output_path=cmatrix_out)

            data['ContactMatrices'].append(cmatrix)

        if distancematrices:
            print(f'Calculating distance matrix for {trajectory_id}...')
            dmatrix_out = os.path.join(output_directory,
                                       f'{trajectory_id}_distance_matrix.csv')

            dmatrix = calculate_distance_matrix(trajectory=trajectory,
                                                topology=topology,
                                                output_path=dmatrix_out)

            data['DistanceMatrices'].append(dmatrix)

        if ssfrequencies:
            print('Calculating secondary structure assignment frequency matrix for '
                  f'{trajectory_id}...')
            ssfreq_out = os.path.join(output_directory,
                                      f'{trajectory_id}_ss_frequency.csv')

            ssfreq = calculate_ss_frequency(trajectory=trajectory,
                                            topology=topology,
                                            output_path=ssfreq_out)

            data['SecondaryStructureFrequencies'].append(ssfreq)

        if rg or dmax or eed or cm_dist:
            metrics_out = os.path.join(output_directory,
                                       f'{trajectory_id}_structural_metrics.csv')
            print(f'Calculating structural metrics data for {trajectory_id}...')
            metrics = calculate_metrics_data(trajectory=trajectory,
                                             topology=topology,
                                             output_path=metrics_out,
                                             rg=rg,
                                             dmax=dmax,
                                             eed=eed,
                                             cm_dist=cm_dist)

            data['StructuralMetrics'].append(metrics)

    return data


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
            Structural Metrics interactive dashboard.
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
        color_palette = ['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                         '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52']

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
