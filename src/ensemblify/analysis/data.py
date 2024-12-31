"""Auxiliary functions for calculating structural properties data."""

# IMPORTS
## Standard Library Imports
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from functools import reduce

## Third Party Imports
import MDAnalysis as mda
import mdtraj
import numpy as np
import pandas as pd
import scipy
from MDAnalysis.analysis import dihedrals
from tqdm import tqdm

## Local Imports
from ensemblify.analysis.third_party.simple_mdreader import SimpleMDreader
from ensemblify.utils import extract_pdb_info

# FUNCTIONS
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

    # Setup multiprocessing variables
    frame_idxs = np.array(range(u.trajectory.n_frames))
    universes = [u] * u.trajectory.n_frames

    if weights is None:
        # Setup multiprocessing variables
        weights = np.array([1/u.trajectory.n_frames] * u.trajectory.n_frames)

        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            contact_matrix_array = reduce(lambda x,y: np.add(x,y),
                                          tqdm(ppe.map(calculate_contact_matrix_frame,
                                                       universes,
                                                       frame_idxs,
                                                       weights),
                                               desc='Calculating contact matrix...',
                                               total=u.trajectory.n_frames))
    else:
        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            contact_matrix_array = reduce(lambda x,y: np.add(x,y),
                                          tqdm(ppe.map(calculate_contact_matrix_frame,
                                                       universes,
                                                       frame_idxs,
                                                       weights),
                                               desc='Calculating reweighted contact matrix...',
                                               total=u.trajectory.n_frames))

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
    distance_matrix = scipy.spatial.distance.cdist(ca_coordinates,
                                                   ca_coordinates,
                                                   'euclidean')

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

    # Setup multiprocessing variables
    frame_idxs = np.array(range(u.trajectory.n_frames))
    universes = [u] * u.trajectory.n_frames

    if weights is None:
        weights = np.array([1/u.trajectory.n_frames] * u.trajectory.n_frames )

        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            distance_matrix_array = reduce(lambda x,y: np.add(x,y),
                                           tqdm(ppe.map(calculate_distance_matrix_frame,
                                                        universes,
                                                        frame_idxs,
                                                        weights),
                                                desc='Calculating distance matrix... ',
                                                total=u.trajectory.n_frames))
    else:
        # Calculate average distance matrix using multiprocessing
        with ProcessPoolExecutor() as ppe:
            distance_matrix_array = reduce(lambda x,y: np.add(x,y),
                                           tqdm(ppe.map(calculate_distance_matrix_frame,
                                                        universes,
                                                        frame_idxs,
                                                        weights),
                                                desc='Calculating reweighted distance matrix... ',
                                                total=u.trajectory.n_frames))

    # Convert calculated averaged matrix to DataFrame
    distance_matrix = pd.DataFrame(data=distance_matrix_array,
                                   index=list(range(1,distance_matrix_array.shape[0]+1)),
                                   columns=[str(x) for x in
                                            range(1,distance_matrix_array.shape[0]+1)])

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
