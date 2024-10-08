"""Auxiliary functions for the analysis module."""

# IMPORTS
## Standard Library Imports
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
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.utils import extract_pdb_info, kde
from ensemblify.analysis.third_party.mdreader_CHANGED import MDreader

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
    output_path: str = os.getcwd(),
    ) -> pd.DataFrame:
    """Calculate a dihedral angles matrix from trajectory and topology files.

    Trajectory and topology files are used to create a MDAnalysis Universe object, from
    which values for the Phi and Psi dihedral angles for each residue in each frame are
    calculated.
    Calculated matrix is saved to output_path, which defaults to the current working directory.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'ramachandran_data.csv'. Defaults to current working directory.

    Returns:
        rama_data:
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
    for row in rama_angles:
        for x,y in row:
            rama_xs.append(x)
            rama_ys.append(y)

    # Create DataFrame with data
    rama_data = pd.DataFrame({'Phi':rama_xs,'Psi':rama_ys})

    # Save dihedrals data
    if os.path.isdir(output_path):
        rama_data.to_csv(os.path.join(output_path,'ramachandran_data.csv'))
    elif output_path.endswith('.csv'):
        rama_data.to_csv(output_path)
    else:
        print(('Ramachandran data was not saved to disk, '
               'output path must be a directory or .csv filepath!'))

    return rama_data


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
    output_path: str | None = None,
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
    # Setup output directory
    if output_path is None:
        output_path = os.getcwd()

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
                                  columns=list(range(1,contact_matrix_array.shape[0]+1)))

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
    reweighted: bool = False,
    difference: bool = False,
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
    if isinstance(contact_matrix,str):
        assert contact_matrix.endswith('.csv'), 'Contact matrix file must be in .csv format!'
        contact_matrix = pd.read_csv(contact_matrix,index_col=0)

    # Extract info regarding chains and resnums
    top_info = extract_pdb_info(topology)
    resranges = {}
    chain_letters = []

    # Start from the last chain as the .pdb was also parsed from the last res
    for chain_number in range(len(top_info.keys()),0,-1):
        chain_info = top_info[chain_number] # (chain_letter, starting_res, chain_size)
        resrange = [ x for x in range(chain_info[1],chain_info[1] + chain_info[2])]
        resranges[chain_info[0]] = resrange
        chain_letters.append(chain_info[0])

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
                text = f'x: {xx}<br />y: {yy}<br />z: {round(contact_matrix.iat[yi,xi],3)}'
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
                                      colorscale=px.colors.diverging.RdBu,
                                      zmid=0,
                                      reversescale=True))

    # Setup chain dividers lines
    num_res = len(contact_matrix.columns)
    chain_ends = [] # to be used in tickvals
    chain_begins = [] # to be used in tickvals
    shapes = []
    cum_res = 0
    for chain_letter in chain_letters[:-1]:
        chain_size = len(resranges[chain_letter])
        chain_begins.append(cum_res)
        chain_end = cum_res + chain_size
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
        cum_res += chain_size
    chain_begins.append(num_res-len(resranges[chain_letters[-1]]))

    # Setup tick values
    tickvals = chain_begins
    curr_chain = 0
    curr_val = 0
    tick_step = num_res // len(chain_letters) // 4 # always 5 ticks per axis
    while curr_val <= num_res:
        try:
            chain_end = chain_ends[curr_chain]
        except IndexError:
            tickvals.append(curr_val)
        else:
            if chain_end - curr_val <= tick_step:
                curr_chain += 1
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

    # Add subtitle
    cmap_fig.add_annotation(text='Frequency of contacts between any atoms of a residue pair',
                            font=dict(family='Helvetica',
                                      color='#707070',
                                      size=24),
                            xref='paper',
                            yref='paper',
                            x=0.5,
                            y=1.05,
                            showarrow=False)

    cmap_fig.update_layout(width=900,
                           height=900,
                           plot_bgcolor='#FFFFFF',
                           font=dict(family='Helvetica',
                                     color='black',
                                     size=30),
                           modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                           title=dict(text=cmap_title,
                                      x=0.5),
                           xaxis=dict(title='Residue',
                                      tickvals=tickvals,
                                      ticks='outside',
                                      ticklen=10,
                                      tickwidth=4,
                                      showgrid=False,),
                           yaxis=dict(title='Residue',
                                      tickvals=tickvals,
                                      ticks='outside',
                                      ticklen=10,
                                      tickwidth=4,
                                      showgrid=False,
                                      title_standoff=5,
                                      scaleanchor='x'),
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
    output_path: str | None = None,
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
    # Setup output directory
    if output_path is None:
        output_path = os.getcwd()

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
                                   columns=list(range(1,distance_matrix_array.shape[0]+1)))

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
    reweighted: bool = False,
    difference: bool = False,
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
            maximum limit for the angstrom distance colorbar. Defaults to None, in which case it is
            derived from the data.
        reweighted:
            boolean stating whether we are creating a reweighted distance matrix figure or a default
            one. Defaults to False.
        difference:
            boolean stating whether we are creating a difference distance matrix figure or a default
            one. Defaults to False.

    Returns:
        dmatrix_fig:
            Ploty Figure object displaying a distance matrix.
    """
    if isinstance(distance_matrix,str):
        assert distance_matrix.endswith('.csv'), 'Distance matrix file must be in .csv format!'
        distance_matrix = pd.read_csv(distance_matrix,index_col=0)

    # Extract info regarding chains and resnums
    top_info = extract_pdb_info(topology)
    resranges = {}
    chain_letters = []

    # Start from the last chain as the .pdb was also parsed from the last res
    for chain_number in range(len(top_info.keys()),0,-1):
        chain_info = top_info[chain_number] # (chain_letter, starting_res, chain_size)
        resrange = [ x for x in range(chain_info[1],chain_info[1] + chain_info[2])]
        resranges[chain_info[0]] = resrange
        chain_letters.append(chain_info[0])

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
        dmatrix_fig.add_trace(go.Heatmap(z=distance_matrix,
                                         x=x_labels,
                                         y=y_labels,
                                         zmin=0,
                                         zmax=max_colorbar,
                                         colorbar_title='&#197;',
                                         transpose=True,
                                         hoverongaps=False,
                                         colorscale=px.colors.sequential.Reds))
    else:
        # Create hovertext
        hovertext = []
        for yi, yy in enumerate(x_labels):
            hovertext.append([])
            for xi, xx in enumerate(x_labels):
                text = f'x: {xx}<br />y: {yy}<br />z: {round(distance_matrix.iat[yi,xi],3)}'
                hovertext[-1].append(text)

        distance_matrix.replace(0,np.nan,inplace=True)
        dmatrix_fig.add_trace(go.Heatmap(z=distance_matrix,
                                         x=x_labels,
                                         y=y_labels,
                                         hoverongaps=False,
                                         colorbar_title='&#197;',
                                         hoverinfo='text',
                                         hovertext=hovertext,
                                         colorscale=px.colors.diverging.RdBu,
                                         zmid=0,
                                         reversescale=True))

    # Setup chain dividers lines
    num_res = len(distance_matrix.columns)
    chain_ends = [] # to be used in tickvals
    chain_begins = [] # to be used in tickvals
    shapes = []
    cum_res = 0
    for chain_letter in chain_letters[:-1]:
        chain_size = len(resranges[chain_letter])
        chain_begins.append(cum_res)
        chain_end = cum_res + chain_size
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
        cum_res += chain_size
    chain_begins.append(num_res-len(resranges[chain_letters[-1]]))

    # Setup tick values
    tickvals = chain_begins
    curr_chain = 0
    curr_val = 0
    tick_step = num_res // len(chain_letters) // 4 # always 5 ticks per axis
    while curr_val <= num_res:
        try:
            chain_end = chain_ends[curr_chain]
        except IndexError:
            tickvals.append(curr_val)
        else:
            if chain_end - curr_val <= tick_step:
                curr_chain += 1
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

    # Add subtitle
    dmatrix_fig.add_annotation(text='Average distance between alpha carbons of a residue pair',
                               font=dict(family='Helvetica',
                                         color='#707070',
                                         size=24),
                               xref='paper',
                               yref='paper',
                               x=0.5,
                               y=1.05,
                               showarrow=False)

    dmatrix_fig.update_layout(width=900,
                              height=900,
                              plot_bgcolor='#FFFFFF',
                              font=dict(family='Helvetica',
                                        color='black',
                                        size=30),
                              modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                              title=dict(text=dmatrix_title,
                                         x=0.5),
                              xaxis=dict(title='Residue',
                                         tickvals=tickvals,
                                         ticks='outside',
                                         ticklen=10,
                                         tickwidth=4,
                                         showgrid=False,),
                              yaxis=dict(title='Residue',
                                         tickvals=tickvals,
                                         ticks='outside',
                                         ticklen=10,
                                         tickwidth=4,
                                         showgrid=False,
                                         title_standoff=5,
                                         scaleanchor='x'),
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
    output_path: str = None,
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
    output_path: str = os.getcwd(),
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
            'ss_assignment.csv'. Defaults to current working directory.

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
        frequency.columns = [int(x) for x in frequency.columns]
    else:
        # Count the frequency of each secondary structure element and reweigh it
        frequency = pd.DataFrame([[0.0]*ss_assignment.shape[1]]*3,
                                 index=['C','E','H'],
                                 columns=[int(x) for x in ss_assignment.columns])
        for row_idx in tqdm(range(ss_assignment.shape[0]),
                            desc='Calculating reweighted secondary structure frequencies... ',
                            total=ss_assignment.shape[0]):
            row_series = ss_assignment.iloc[row_idx,:]
            weight = weights[row_idx]
            for col_idx in range(1,frequency.shape[1]+1):
                ssa_label = row_series.iloc[col_idx-1]
                frequency.loc[ssa_label,col_idx] += weight

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
    if isinstance(ss_frequency,str):
        assert ss_frequency.endswith('.csv'), ('Secondary structure assignment frequency '
                                                'matrix must be in .csv format!')
        ss_frequency = pd.read_csv(ss_frequency,index_col=0)

    # Extract info regarding chains and resnums
    top_info = extract_pdb_info(topology)
    resranges = {}
    chain_letters = []

    # Iterate through chains
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
            ss_freq_fig.add_trace(go.Scatter(x=list(range(ss_frequency.columns[0],
                                                          ss_frequency.columns[-1]+1)),
                                             y=ss_frequency.loc[structure],
                                             mode='lines',
                                             marker_color=color,
                                             line_width=4,
                                             name=structure,
                                             hoverinfo='text',
                                             hovertext=hovertext))
    else:
        for structure,color in zip(ss_frequency.index,colors):
            ss_freq_fig.add_trace(go.Scatter(x=ss_frequency.columns,
                                             y=ss_frequency.loc[structure],
                                             mode='lines',
                                             stackgroup='one', # remove for non stacked plot
                                             marker_color=color,
                                             line_width=0,
                                             name=structure))

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
                            x0=cumulative_residues+len(resranges[chain_letter])-1,
                            x1=cumulative_residues+len(resranges[chain_letter])-1,
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
            ss_freq_title = f'{trajectory_id} Difference Secondary Structure Frequencies'
        else:
            ss_freq_title = 'Difference Secondary Structure Frequencies'
    elif reweighted:
        if trajectory_id is not None:
            ss_freq_title = f'{trajectory_id} Reweighted Secondary Structure Frequencies'
        else:
            ss_freq_title = 'Reweighted Secondary Structure Frequencies'
    else:
        if trajectory_id is not None:
            ss_freq_title = f'{trajectory_id} Secondary Structure Frequencies'
        else:
            ss_freq_title = 'Secondary Structure Frequencies'

    # Add subtitle
    ss_freq_fig.add_annotation(text='Frequency of each secondary structure assignment code '
                                    'for each residue',
                               font=dict(family='Helvetica',
                                         color='#707070',
                                         size=24),
                               xref='paper',
                               yref='paper',
                               x=0.5,
                               y=1.07,
                               showarrow=False)

    if difference:
        range_y = [-1,1]
    else:
        range_y = [0,1]

    ss_freq_fig.update_layout(width=1000,
                              height=750,
                              font=dict(family='Helvetica',
                                        color='black',
                                        size=30),
                              plot_bgcolor = '#FFFFFF',
                              paper_bgcolor = '#FFFFFF',
                              modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                              title=dict(text=ss_freq_title,
                                         x=0.5),
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
    rg: bool = True,
    dmax: bool = True,
    eed: bool = True,
    cm_dist: dict[str,tuple[str,str]] | None = None,
    output_path: str | None = None,
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
            'structural_metrics.csv'.

    Returns:
        traj_analysis:
            DataFrame where columns are the desired structural metrics and rows are the frames
            of the trajectory.
    """
    # Initialize MDReader
    u = MDreader(f'-f {trajectory} -s {topology}'.split())

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
    traj_analysis = pd.DataFrame(metrics_array.reshape(-1,
                                                       metrics_array.shape[-1]),
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
    color: str = '#000000',
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
            hex code for the color the created traces will be. Defaults to black.

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
        hist_trace = go.Histogram(x=metrics[col_name],
                                  name=f'{trajectory_id}_{col_name}',
                                  xbins=dict(start=min(metrics[col_name]),
                                             # To include the last value in the last bin
                                             end=max(metrics[col_name]) + 1,
                                             size=1),
                                  histnorm='probability density',
                                  marker=dict(color=color),
                                  opacity=0.7,
                                  visible=False)

        hist_traces.append(hist_trace)

        kde_x, kde_y, avg, avg_stderr = kde(data=np.array(metrics[col_name]))

        scatter_trace = go.Scatter(x=kde_x,
                                   y=kde_y,
                                   mode='lines',
                                   name=f'{trajectory_id}_{col_name}',
                                   marker_color=color,
                                   line=dict(width=4),
                                   legend='legend',
                                   visible=True)

        scatter_traces.append(scatter_trace)
        avg_values.append(avg)
        avg_stderr_values.append(avg_stderr)

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
    nrows = len(trajectory_ids) + 1
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

    # Create metrics figure
    metrics_fig = make_subplots(rows=nrows,
                                cols=nmetrics,
                                row_heights=[0.75]*(nrows-1) + [1.75],
                                column_titles=col_titles,
                                row_titles=trajectory_ids+[''],
                                horizontal_spacing=0.085)


    # Add the Box traces for all rows except last
    for rownum in range(1,nrows):
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
                                  row=nrows,
                                  col=colnum)
            metrics_fig.add_trace(scatter_trace,
                                  row=nrows,
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
                                       row=nrows,
                                       col=colnum)

            # Allows for hovering the dashed line to get mean value
            metrics_fig.add_trace(go.Scatter(x=[mean_value],
                                             y=[x for x in np.arange(0,
                                                                     np.interp(mean_value,
                                                                               scatter_trace.x,
                                                                               scatter_trace.y)+0.001,
                                                                               0.001)],
                                             mode='markers',
                                             marker_color=hist_trace.marker.color,
                                             hovertext=f'Avg: {round(mean_value,2)} &plusmn; {round(mean_value_stderr,2)}',
                                             hoverinfo='text',
                                             hoverlabel_bgcolor=hist_trace.marker.color,
                                             fill='toself',
                                             name=f'Avg: {round(mean_value,2)} &plusmn; {round(mean_value_stderr,2)}',
                                             opacity=0,
                                             showlegend=False),
                                  row=nrows,
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
                                     row=nrows,
                                     col=colnum)
            metrics_fig.update_yaxes(ticks='outside',
                                     title_text=f'KDE ({metric_name})',
                                     rangemode='tozero',
                                     row=nrows,
                                     col=colnum)

    # Equalize x axis range
    for rownum in range(1,nrows+1):
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
                                    xanchor='center',
                                    x=1.06,
                                    y=1.13)

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
                                     xanchor='center',
                                     x=1.06,
                                     y=1.05)

    metrics_fig.update_layout(height=200 + 100 *nrows-1 + 350,
                              width=610*nmetrics+450,
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
                             ticklen=10,
                             tickwidth=4,
                             linewidth=4,
                             linecolor='black',
                             color='black',
                             mirror=True,
                             title_standoff=0)

    metrics_fig.update_yaxes(row=nrows,
                             col=1,
                             title_standoff=0)

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
    output_directory: str | None = None,
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
                    'StructuralMetrics' : [StructuralMetrics1,StructuralMetrics2,StructuralMetrics3]}
    """
    # Setup output directory
    if output_directory is None:
        output_directory = os.getcwd()

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
            ramachandran_data_output_path = os.path.join(output_directory,
                                                         f'{trajectory_id}_ramachandran_data.csv')

            calculate_ramachandran_data(trajectory=trajectory,
                                        topology=topology,
                                        output_path=ramachandran_data_output_path)

        # Analysis meant for interactive figures
        if distancematrices:
            print(f'Calculating distance matrix for {trajectory_id}...')
            distance_matrix_output_path = os.path.join(output_directory,
                                                      f'{trajectory_id}_distance_matrix.csv')

            distance_matrix = calculate_distance_matrix(trajectory=trajectory,
                                                        topology=topology,
                                                        output_path=distance_matrix_output_path)

            data['DistanceMatrices'].append(distance_matrix)

        if contactmatrices:
            print(f'Calculating contact matrix for {trajectory_id}...')
            contact_matrix_output_path = os.path.join(output_directory,
                                                      f'{trajectory_id}_contact_matrix.csv')

            contact_matrix = calculate_contact_matrix(trajectory=trajectory,
                                                      topology=topology,
                                                      output_path=contact_matrix_output_path)

            data['ContactMatrices'].append(contact_matrix)

        if ssfrequencies:
            print('Calculating secondary structure assignment frequency matrix for '
                  f'{trajectory_id}...')
            ss_frequency_output_path = os.path.join(output_directory,
                                                     f'{trajectory_id}_ss_frequency.csv')

            ss_frequency = calculate_ss_frequency(trajectory=trajectory,
                                                  topology=topology,
                                                  output_path=ss_frequency_output_path)

            data['SecondaryStructureFrequencies'].append(ss_frequency)

        if rg or dmax or eed or cm_dist:
            structural_metrics_output_path = os.path.join(output_directory,
                                                          f'{trajectory_id}_structural_metrics.csv')
            print(f'Calculating structural metrics data for {trajectory_id}...')
            structural_metrics = calculate_metrics_data(trajectory=trajectory,
                                                        topology=topology,
                                                        output_path=structural_metrics_output_path,
                                                        rg=rg,
                                                        dmax=dmax,
                                                        eed=eed,
                                                        cm_dist=cm_dist)

            data['StructuralMetrics'].append(structural_metrics)

    return data


def create_analysis_figures(
    analysis_data: dict[str,list[pd.DataFrame]] | None,
    topologies: list[str],
    trajectory_ids: list[str],
    output_directory: str | None = None,
    color_palette: list[str] | None = None,
    ) -> dict[str,list[go.Figure]]:
    """Create interactive figures given analysis data for one or more trajectory,topology
    pair of files.

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
            data = {'DistanceMatrices' : [DistanceMatrix1,DistanceMatrix2,DistanceMatrix3],
                    'ContactMaps' : [ContactMap1,ContactMap2,ContactMap3],
                    'SecondaryStructureFrequencies' : [SSFrequency1,SSFrequency2,SSFrequency3],
                    'StructuralMetrics' : [StructuralMetrics1,StructuralMetrics2,StructuralMetrics3] }
    """
    # Setup output directory
    if output_directory is None:
        output_directory = os.getcwd()
    
    # Setup color palette
    if color_palette is None:
        color_palette = ['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                         '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52']

    # If data is not available check output directory for it
    if analysis_data is None:
        analysis_data = {'DistanceMatrices' : [],
                         'ContactMatrices' : [],
                         'SecondaryStructureFrequencies' : [],
                         'StructuralMetrics' : [] }
        data_ids = ['distance_matrix.csv','contact_matrix.csv','ss_frequency.csv','structural_metrics.csv']
        data_ids_2_data = {'distance_matrix.csv': 'DistanceMatrices',
                           'contact_matrix.csv': 'ContactMatrices',
                           'ss_frequency.csv': 'SecondaryStructureFrequencies',
                           'structural_metrics.csv': 'StructuralMetrics'}
        for data_id in data_ids:
            for trajectory_id in trajectory_ids:
                try:
                    filename = os.path.join(output_directory,f'{trajectory_id}_{data_id}')
                    calculated_data = pd.read_csv(filename,index_col=0)
                    analysis_data[data_ids_2_data[data_id]].append(calculated_data)
                except FileNotFoundError:
                    continue
                else:
                    print(f'Found calculated {trajectory_id}_{data_id}')

    # Create figures
    figures = {'DistanceMatrices' : [],
               'ContactMaps' : [],
               'SecondaryStructureFrequencies' : [],
               'StructuralMetrics' : None }

    for i,(trajectory_id,topology,color) in enumerate(zip(trajectory_ids,topologies,color_palette)):
        print(f'Creating {trajectory_id} analysis figures...')

        try:
            distance_matrix = analysis_data['DistanceMatrices'][i]
        except (KeyError,IndexError):
            pass
        else:
            distance_matrix_fig_output_path = os.path.join(output_directory,
                                                       f'{trajectory_id}_distance_matrix.html')

            distance_matrix_fig = create_distance_matrix_fig(distance_matrix=distance_matrix,
                                                             trajectory_id=trajectory_id,
                                                             topology=topology,
                                                             output_path=distance_matrix_fig_output_path)

            figures['DistanceMatrices'].append(distance_matrix_fig)

        try:
            contact_matrix = analysis_data['ContactMatrices'][i]
        except (KeyError,IndexError):
            pass
        else:
            contact_map_fig_output_path = os.path.join(output_directory,
                                                       f'{trajectory_id}_contact_map.html')

            contact_map_fig = create_contact_map_fig(contact_matrix=contact_matrix,
                                                     trajectory_id=trajectory_id,
                                                     topology=topology,
                                                     output_path=contact_map_fig_output_path)

            figures['ContactMaps'].append(contact_map_fig)

        try:
            ss_frequency = analysis_data['SecondaryStructureFrequencies'][i]
        except (KeyError,IndexError):
            pass
        else:
            ss_frequency_fig_output_path = os.path.join(output_directory,
                                                        f'{trajectory_id}_ss_frequency.html')

            ss_frequency_fig = create_ss_frequency_figure(ss_frequency=ss_frequency,
                                                          topology=topology,
                                                          trajectory_id=trajectory_id,
                                                          output_path=ss_frequency_fig_output_path)

            figures['SecondaryStructureFrequencies'].append(ss_frequency_fig)

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
        metrics_fig_output_path = os.path.join(output_directory,
                                               'structural_metrics.html')

        metrics_fig = create_metrics_fig(trajectory_ids=trajectory_ids,
                                         total_box_traces=total_box_traces,
                                         total_hist_traces=total_hist_traces,
                                         total_scatter_traces=total_scatter_traces,
                                         total_avg_values=total_avg_values,
                                         total_avg_stderr_values=total_avg_stderr_values,
                                         output_path=metrics_fig_output_path)

        figures['StructuralMetrics'] = metrics_fig

    return figures
