"""Reweight a conformational ensemble using experimental data."""

# IMPORTS
## Standard Library Imports
import glob
import os
import shutil
import subprocess
import sys

## Third Party Imports
import numpy as np
import pandas as pd
from plotly.offline import get_plotlyjs

## Local Imports
from ensemblify.analysis import (
    calculate_contact_matrix,
    calculate_distance_matrix,
    calculate_metrics_data,
    calculate_ss_frequency,
    create_contact_map_fig,
    create_distance_matrix_fig,
    create_ss_frequency_figure,
)
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.conversion import traj2saxs, calc_chi2_fit
from ensemblify.reweighting.data import (
    attempt_read_calculated_data,
    attempt_read_reweighting_data,
    average_saxs_profiles,
    bme_ensemble_reweighting,
    correct_exp_saxs_error,
    process_exp_saxs_data,
)
from ensemblify.reweighting.figures import (
    create_effective_frames_fit_fig,
    create_rw_saxs_fits_fig,
    create_reweighting_metrics_fig,
)
from ensemblify.utils import get_array_extremum, round_to_nearest_multiple

# FUNCTIONS
def reweight_ensemble(
    trajectory: str,
    topology: str,
    trajectory_id: str,
    exp_data: str | list[str],
    exp_type: str | list[str] | None = None,
    output_dir: str | None = None,
    thetas: list[int] | None = None,
    calculated_SAXS_data: np.ndarray | str | None = None,
    calculated_cmatrix: pd.DataFrame | str | None = None,
    calculated_dmatrix: pd.DataFrame | str | None = None,
    calculated_ss_frequency: pd.DataFrame | str | None = None,
    calculated_metrics_data: pd.DataFrame | str | None = None,
    compare_rg: bool | None = True,
    compare_dmax: bool | None = True,
    compare_eed: bool | None = True,
    compare_cmdist: bool | None = None,
    pepsi_saxs_opt: str | None = None,
    ):
    """Apply Bayesian Maximum Entropy (BME) reweighting to a conformational ensemble.

    If exp_type if not provided, it must be specified in the first line of the exp_data file, in
    the form of 'DATA=<type>', where <type> is one of the following: ['CS','JCOUPLINGS','RDC',
    'SAXS'].
    Inside the output directory, a directory named trajectory_id will be created where all output
    files will be stored.
    If calculated metrics data is provided, it will not be recalculated.
    If data for the center mass distance is to be taken from the given calculated metrics data,
    the compare_cmdist mapping must be provided. The identifiers of this mapping will be matched
    to column names of the given DataFrame, if present.
 
    Args:
        trajectory (str):
            Path to .xtc trajectory file where conformational ensemble is stored.
        topology (str):
            Path to .pdb topology file corresponding to any one frame of the ensemble.
        trajectory_id (str):
            Prefix trajectory identifier to be added to plotted traces and output files.
        exp_data (str | list[str]):
            Path(s) to .dat file(s) with experimental data with 3 columns {Label, Value, Error}.
        exp_type (str, optional):
            Type(s) of experimental data provided in the same order as exp_data. If None, the
            type(s) are taken from the first line of each exp_data file.
        output_dir (str, optional):
            Path to output directory. Is created if it does not exist. Defaults to current working
            directory. After output_dir is setup, a directory named trajectory_id is created inside
            it, where the interactive .html plots and reweighting output files will be stored.
        thetas (list[int], optional):
            List of values to try as the theta parameter in BME. The ensemble will be reweighted
            each time using a different theta value. The effect of different theta values can be
            analyzed in the created effective frames figure.
        calculated_SAXS_data (np.ndarray | str , optional): 
            Array with calculated SAXS profiles for each conformer in the ensemble or path to 
            this file in .dat format. Number of rows must match number of different conformers
            in the ensemble. Number of calculated data points (columns) must match the number
            of experimental data points provided in exp_data. If provided, this data is used
            instead of back-calculating it from the ensemble. Defaults to None.
        calculated_cmatrix (pd.DataFrame | str, optional):
            DataFrame with the calculated average contact matrix for the current trajectory or
            path to this file in .csv format. Defaults to None, and this data is calculated anew.
        calculated_dmatrix (pd.DataFrame | str, optional):
            DataFrame with the calculated average distance matrix for the current trajectory or
            path to this file in .csv format. Defaults to None, and this data is calculated anew.
        calculated_ss_frequency (pd.DataFrame | str, optional):
            DataFrame with the calculated secondary structure assignment frequency matrix for the
            current trajectory or path to this file in .csv format. Defaults to None, and this data
            is calculated anew.
        calculated_metrics_data (pd.DataFrame | str, optional):
            DataFrame with calculated structural metrics (columns) for each frame of the trajectory
            (rows) or path to this DataFrame in .csv format. Defaults to None, and this data is
            calculated anew.
        compare_rg (bool, optional):
            Whether to calculate/consider the radius of gyration when comparing structural metrics
            between uniform and reweighted conformational ensembles. Defaults to True.
        compare_dmax (bool, optional):
            Whether to calculate/consider the maximum distance between any two alpha carbons when
            comparing structural metrics between uniform and reweighted conformational ensembles.
            Defaults to True.
        compare_eed (bool, optional):
            Whether to calculate/consider the distance from the N to C terminal when comparing
            structural metrics between uniform and reweighted conformational ensembles. Defaults
            to True.
        compare_cmdist (bool, optional):
            Mapping of identifiers to tuples with two selection strings for creating MDAnalysis
            AtomGroups, whose center mass distance will be calculated. For example:

                {'inter_domain' : ('resid 1:30', 'resid 110:140')}

            If None, no center mass distances are calculated or compared.
            See https://userguide.mdanalysis.org/stable/selections.html for more information about
            MDAnalysis selections.
            Defaults to None.
        pepsi_saxs_opt (str, optional):
            If back-calculation of SAXS data from the ensemble will be attempted, this string
            will be passed onto Pepsi-SAXS as additional command line options. If None, default
            Pepsi-SAXS options are used instead.
    """
    # Setup experimental data inputs
    if isinstance(exp_data, str):
        exp_data_input = [exp_data]
    elif isinstance(exp_data, list):
        exp_data_input = [x for x in exp_data]
    else:
        raise TypeError('exp_data must be a string or a list of strings.')

    # Setup output directory
    if output_dir is None:
        output_dir = os.getcwd()
    elif not os.path.isdir(output_dir):
        os.mkdir(output_dir)
       
    # Setup experimental data type(s) if not provided
    if exp_type is None:
        exp_type = []
        for data_file in exp_data_input:
            # Read first line of experimental data file to get type
            with open(data_file, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
                if 'DATA=' in first_line:
                    for equality in first_line.split(' '):
                        if 'DATA=' in equality:
                            exp_type.append(equality.split('=')[1].upper())
                            break
                else:
                    raise ValueError('Experimental data type not specified in the first line of '
                                     'the experimental data file. Please add it or specify it '
                                     'using the exp_type argument.')

    # Copy experimental data file(s) into output dir before working on it
    exp_data = []
    for data_file, data_type in zip(exp_data_input,exp_type):
        exp_data_copy = os.path.join(output_dir,
                                     f'{trajectory_id}_{data_type}_exp_input.dat')
        shutil.copy(src=data_file,
                    dst=exp_data_copy)
        exp_data.append(exp_data_copy)

    # Setup directory for reweighting files
    reweighting_dir = os.path.join(output_dir,
                                   'bme_reweighting_results')
    if not os.path.isdir(reweighting_dir):
        os.mkdir(reweighting_dir)

    # Setup theta values
    if thetas is None:
        thetas = [1, 10, 20, 50, 75, 100, 200, 400, 750, 1000, 5000, 10000]

    # Check if we can skip some steps by reading previously computed data from output directory
    exp_data_BME, \
    calc_data_BME, \
    thetas_array, \
    stats, \
    weights = attempt_read_reweighting_data(reweighting_output_directory=output_dir,
                                            trajectory_id=trajectory_id,
                                            exp_type=exp_type)

    # Check read exp_data_BME
    for i, (data_file,data_type) in enumerate(zip(exp_data,exp_type)):
        if exp_data_BME[i] is None:

            # Process input experimental data
            print(f'Processing {trajectory_id} experimental {data_type} file...')

            # If SAXS, additional processing
            if data_type == 'SAXS':

                # Check units
                exp_saxs_file_processed = process_exp_saxs_data(experimental_data_path=data_file)

                # Correct errors
                try:
                    exp_file_BME = correct_exp_saxs_error(experimental_data_path=exp_saxs_file_processed)
                except subprocess.CalledProcessError: # if BIFT is not available
                    exp_file_BME = os.path.join(os.path.dirname(exp_saxs_file_processed),
                                                f'{trajectory_id}_SAXS_exp.dat')
                    np.savetxt(exp_file_BME,
                               np.loadtxt(exp_saxs_file_processed),
                               header=' DATA=SAXS')

                # Clean up intermediate files
                os.remove(data_file)
                os.remove(exp_saxs_file_processed)

            else:
                print(f'{data_type} data is currently not supported.')
                sys.exit(1)

            # Add processed experimental data file to list
            exp_data_BME[i] = exp_file_BME

    # Check read calc_data_BME
    for i, (data_file,data_type) in enumerate(zip(exp_data_BME,exp_type)):
        if calc_data_BME[i] is None:
            
            if data_type == 'SAXS':
                
                if calculated_SAXS_data:
                    all_calc_file = calculated_SAXS_data
                else:
                    # Calculate SAXS data from ensemble
                    all_calc_file, _, _ = traj2saxs(trajectory=trajectory,
                                                    topology=topology,
                                                    trajectory_id=trajectory_id,
                                                    exp_saxs_file=data_file,
                                                    pepsi_saxs_opt=pepsi_saxs_opt)

            else:
                print((f'No method to automatically back-calculate {data_type} data is currently '
                       'available. Please provide the file(s) with your calculated experimental '
                       f'data using the calculated_{data_type}_data parameter.'))
                sys.exit(1)

            # Add calculated data file to list
            calc_data_BME[i] = all_calc_file

    if thetas_array is None or stats is None or weights is None:
        # Reweigh ensemble using different theta values
        thetas_array = np.array(thetas)

        print((f'Applying BME reweighting to {trajectory_id} ensemble '
               f'with theta values {thetas} ...'))
        stats, weights = bme_ensemble_reweighting(exp_data=exp_data_BME,
                                                  exp_type=exp_type,
                                                  calc_data=calc_data_BME,
                                                  thetas=thetas_array,
                                                  output_dir=reweighting_dir)

    effective_frames_fit_fig = create_effective_frames_fit_fig(stats=stats,
                                                               thetas=thetas_array,
                                                               title_text=trajectory_id)

    effective_frames_fit_fig.write_html(os.path.join(output_dir,'effective_frames_fit.html'),
                                        config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

    print(('Please analyze the provided interactive figure (effective_frames_fit.html) and '
           'input the desired value(s) for the theta parameter.\nIf more than one '
           'value, please separate them using a comma.'),
           flush=True)

    # Capture chosen theta values
    input_choices = input('Choose theta(s): ').split(',')
    print(f'Chosen theta value(s): {", ".join(input_choices)}.')
    chosen_thetas = [int(x) for x in input_choices]

    # Plot L curve with chosen theta(s)
    chosen_thetas_fig = create_effective_frames_fit_fig(stats=stats,
                                                        thetas=thetas_array,
                                                        choices=chosen_thetas,
                                                        title_text=(f'{trajectory_id} BME '
                                                                    'Reweighting'))

    # Extract correct set(s) of weights using index(es) for chosen theta(s)
    choice_idxs = [np.where(thetas_array == x)[0][0] for x in chosen_thetas]
    chosen_weight_sets = [ weights[i] for i in choice_idxs]

    ##############################################################################################
    ############################# CALCULATE REWEIGHTING FIGURES DATA #############################
    ##############################################################################################

    # Calculate Experimental Fittings
    for exp_data_file, calc_data_file, data_type in zip(exp_data_BME, calc_data_BME, exp_type):
        print(f'Calculating reweighting figures data for {trajectory_id} experimental {data_type} '
              f'data...')
        
        if data_type == 'SAXS':

            # Calculate prior and posterior average SAXS intensities
            i_posts = []

            for chosen_theta,chosen_weight_set in zip(chosen_thetas,chosen_weight_sets):

                # Grab reweighted SAXS data using choice theta
                rw_calc_saxs_file = glob.glob(os.path.join(reweighting_dir,
                                                           f'ibme_t{chosen_theta}_*.calc.dat'))[0]

                # Calculate uniform (prior) and reweighted (posterior) average SAXS profiles
                common_i_prior, i_post = average_saxs_profiles(exp_saxs_file=exp_data_file,
                                                               calc_saxs_file=calc_data_file,
                                                               rw_calc_saxs_file=rw_calc_saxs_file,
                                                               weights=chosen_weight_set)

                i_posts.append(i_post)

                # Fit reweighted SAXS data to experimental SAXS data, save chisquare and residuals
                ## Get experimental intensities and errors
                q, i_exp, err = np.loadtxt(data_file,
                                           unpack=True)

                chi2, residuals = calc_chi2_fit(exp=np.hstack((i_exp.reshape(-1,1),
                                                               err.reshape(-1,1))),
                                                calc=np.loadtxt(rw_calc_saxs_file)[...,1:],
                                                sample_weights=chosen_weight_set)
                ## Save fitting data
                fitting_saxs_file = os.path.join(os.path.dirname(data_file),
                                     f'{trajectory_id}_SAXS_fitting_t{chosen_theta}.dat')
                chi2_header = f'Chi^2 = {chi2}\n'
                col_names = ['Momentum Vector', 'Experimental SAXS Intensity', 'Experimental Error',
                             'Calculated SAXS Intensity', 'Residuals of Fit']
                file_header = chi2_header + '\t'.join(col_names)
                
                np.savetxt(fitting_saxs_file,
                           np.hstack((q.reshape(-1,1),
                                      i_exp.reshape(-1,1),
                                      err.reshape(-1,1),
                                      i_post.reshape(-1,1),
                                      residuals.reshape(-1,1))),
                           header=file_header,
                           encoding='utf-8')

    # Calculate Contact Matrices
    ## Uniform
    cmatrix = attempt_read_calculated_data(data=calculated_cmatrix,
                                           data_msg_tag='cmatrix',
                                           calc_fn=calculate_contact_matrix,
                                           trajectory=trajectory,
                                           topology=topology,
                                           output_path=os.path.join(output_dir,
                                                                    (f'{trajectory_id}_contact_'
                                                                     'matrix.csv')))
    ## Reweighted
    rw_cmatrices = []
    diff_cmatrices = []
    for chosen_theta,chosen_weight_set in zip(chosen_thetas,chosen_weight_sets):
        rw_cmatrix = calculate_contact_matrix(trajectory=trajectory,
                                              topology=topology,
                                              weights=chosen_weight_set,
                                              output_path=os.path.join(output_dir,
                                                                       f'{trajectory_id}'
                                                                       '_contact_matrix_'
                                                                       f't{chosen_theta}_'
                                                                       'reweighted.csv'))
        rw_cmatrices.append(rw_cmatrix)
        diff_cmatrices.append(rw_cmatrix - cmatrix)

    # Calculate Distance Matrices
    ## Uniform
    dmatrix = attempt_read_calculated_data(data=calculated_dmatrix,
                                           data_msg_tag='dmatrix',
                                           calc_fn=calculate_distance_matrix,
                                           trajectory=trajectory,
                                           topology=topology,
                                           output_path=os.path.join(output_dir,
                                                                    (f'{trajectory_id}_distance_'
                                                                     'matrix.csv')))
    ## Reweighted
    rw_dmatrices = []
    diff_dmatrices = []
    for chosen_theta,chosen_weight_set in zip(chosen_thetas,chosen_weight_sets):
        rw_dmatrix = calculate_distance_matrix(trajectory=trajectory,
                                               topology=topology,
                                               weights=chosen_weight_set,
                                               output_path=os.path.join(output_dir,
                                                                        f'{trajectory_id}'
                                                                        '_distance_matrix_'
                                                                        f't{chosen_theta}_'
                                                                        'reweighted.csv'))
        rw_dmatrices.append(rw_dmatrix)
        diff_dmatrices.append(rw_dmatrix - dmatrix)

    # Calculate Secondary Structure Assignment Frequency Matrices
    ## Uniform
    ssfreq = attempt_read_calculated_data(data=calculated_ss_frequency,
                                          data_msg_tag='ss_freq',
                                          calc_fn=calculate_ss_frequency,
                                          trajectory=trajectory,
                                          topology=topology,
                                          output_path=os.path.join(output_dir,
                                                                   (f'{trajectory_id}_ss_frequency'
                                                                    '.csv')))

    ## Reweighted
    rw_ssfreqs = []
    diff_ssfreqs = []
    for chosen_theta,chosen_weight_set in zip(chosen_thetas,chosen_weight_sets):
        rw_ssfreq = calculate_ss_frequency(trajectory=trajectory,
                                           topology=topology,
                                           weights=chosen_weight_set,
                                           output_path=os.path.join(output_dir,
                                                                    f'{trajectory_id}'
                                                                    '_ss_frequency_'
                                                                    f't{chosen_theta}_'
                                                                    'reweighted.csv'))
        rw_ssfreqs.append(rw_ssfreq)
        diff_ssfreqs.append(rw_ssfreq - ssfreq)

    # Calculate Uniform Structural Metrics Distributions
    metrics = attempt_read_calculated_data(data=calculated_metrics_data,
                                           data_msg_tag='structural_metrics',
                                           calc_fn=calculate_metrics_data,
                                           trajectory=trajectory,
                                           topology=topology,
                                           rg=compare_rg,
                                           dmax=compare_dmax,
                                           eed=compare_eed,
                                           cm_dist=compare_cmdist,
                                           output_path=os.path.join(output_dir,
                                                                    f'{trajectory_id}'
                                                                    '_structural_metrics.csv'))

    ##############################################################################################
    ################################# CREATE REWEIGHTING FIGURES #################################
    ##############################################################################################
    # Create interactive figures
    print(f'Creating {trajectory_id} reweighted interactive figures...')

    # Experimental fitting figures
    rw_fits_figs = []
    for data_file, data_type in zip(exp_data_BME, exp_type):
        print(f'Creating {trajectory_id} reweighted interactive figures for experimental '
              f'{data_type} data...')
        
        if data_type == 'SAXS':
            # Read experimental data
            q, i_exp, err = np.loadtxt(data_file,
                                       unpack=True)

            # Create fittings Figure
            rw_saxs_fits_fig = create_rw_saxs_fits_fig(q=q,
                                                       i_exp=i_exp,
                                                       err=err,
                                                       i_prior=common_i_prior,
                                                       i_posts=i_posts,
                                                       title_text=(f'{trajectory_id} '
                                                                   'Reweighted SAXS '
                                                                   'Fittings'))
            rw_fits_figs.append(rw_saxs_fits_fig)

    # Create nested dictionary to split up reweighted figures according to theta values
    theta_2_reweighted_figures = {}
    for chosen_theta in chosen_thetas:
        theta_2_reweighted_figures[chosen_theta] = {'cmap': None,
                                                    'rw_cmap': None,
                                                    'diff_cmap': None,
                                                    'dmatrix': None,
                                                    'rw_dmatrix': None,
                                                    'diff_dmatrix': None,
                                                    'ssfreq': None,
                                                    'rw_ssfreq': None,
                                                    'diff_ssfreq': None
                                                    }
    # Contact Maps
    for chosen_theta, rw_cm, diff_cm in zip(chosen_thetas, rw_cmatrices,diff_cmatrices):
        ## Uniform
        cmap_fig = create_contact_map_fig(contact_matrix=cmatrix,
                                          topology=topology,
                                          trajectory_id=trajectory_id,
                                          output_path=output_dir)
        theta_2_reweighted_figures[chosen_theta]['cmap'] = cmap_fig

        ## Reweighted
        rw_cmap_fig = create_contact_map_fig(contact_matrix=rw_cm,
                                             topology=topology,
                                             trajectory_id=trajectory_id,
                                             output_path=os.path.join(output_dir,
                                                                      (f'{trajectory_id}_'
                                                                       'contact_map_'
                                                                       f't{chosen_theta}_'
                                                                       'reweighted.html')),
                                             reweighted=True)
        theta_2_reweighted_figures[chosen_theta]['rw_cmap'] = rw_cmap_fig

        diff_cmap_fig = create_contact_map_fig(contact_matrix=diff_cm,
                                               topology=topology,
                                               trajectory_id=trajectory_id,
                                               output_path=output_dir,
                                               difference=True)
        theta_2_reweighted_figures[chosen_theta]['diff_cmap'] = diff_cmap_fig

    # Distance Matrices
    ## Get maximum of colorbar
    max_data = get_array_extremum([dmatrix] + rw_dmatrices)
    max_colorbar = round_to_nearest_multiple(max_data,5)

    ## Get max/min of difference colorbar
    max_diff_data = get_array_extremum(diff_dmatrices)
    min_diff_data = get_array_extremum(diff_dmatrices,get_max=False)
    max_diff_colorbar = round_to_nearest_multiple(max_diff_data,2)
    min_diff_colorbar = round_to_nearest_multiple(min_diff_data,-2,up=False)

    if abs(max_diff_colorbar) > abs(min_diff_colorbar):
        min_diff_colorbar = - max_diff_colorbar
    else:
        max_diff_colorbar = - min_diff_colorbar

    ## Create figures
    for chosen_theta, rw_dm, diff_dm in zip(chosen_thetas,rw_dmatrices,diff_dmatrices):
        ## Uniform
        dmatrix_fig = create_distance_matrix_fig(distance_matrix=dmatrix,
                                                 topology=topology,
                                                 trajectory_id=trajectory_id,
                                                 output_path=output_dir,
                                                 max_colorbar=max_colorbar)
        theta_2_reweighted_figures[chosen_theta]['dmatrix'] = dmatrix_fig

        ## Reweighted
        rw_dmatrix_fig = create_distance_matrix_fig(distance_matrix=rw_dm,
                                                    topology=topology,
                                                    trajectory_id=trajectory_id,
                                                    output_path=os.path.join(output_dir,
                                                                             (f'{trajectory_id}_'
                                                                              'distance_matrix_'
                                                                              f't{chosen_theta}_'
                                                                              'reweighted.html')),
                                                    max_colorbar=max_colorbar,
                                                    reweighted=True)
        theta_2_reweighted_figures[chosen_theta]['rw_dmatrix'] = rw_dmatrix_fig

        diff_dmatrix_fig = create_distance_matrix_fig(distance_matrix=diff_dm,
                                                      topology=topology,
                                                      trajectory_id=trajectory_id,
                                                      output_path=output_dir,
                                                      max_colorbar=max_diff_colorbar,
                                                      min_colorbar=min_diff_colorbar,
                                                      difference=True)
        theta_2_reweighted_figures[chosen_theta]['diff_dmatrix'] = diff_dmatrix_fig

    # Secondary Structure Frequencies
    for chosen_theta, rw_ssf, diff_ssf in zip(chosen_thetas,rw_ssfreqs,diff_ssfreqs):
        ## Uniform
        ssfreq_fig = create_ss_frequency_figure(ss_frequency=ssfreq,
                                                topology=topology,
                                                trajectory_id=trajectory_id,
                                                output_path=output_dir)
        theta_2_reweighted_figures[chosen_theta]['ssfreq'] = ssfreq_fig

        ## Reweighted
        rw_ssfreq_fig = create_ss_frequency_figure(ss_frequency=rw_ssf,
                                                   topology=topology,
                                                   trajectory_id=trajectory_id,
                                                   output_path=os.path.join(output_dir,
                                                                            (f'{trajectory_id}_'
                                                                             'ss_frequency_'
                                                                             f't{chosen_theta}_'
                                                                             'reweighted.html')),
                                                   reweighted=True)
        theta_2_reweighted_figures[chosen_theta]['rw_ssfreq'] = rw_ssfreq_fig

        diff_ssfreq_fig = create_ss_frequency_figure(ss_frequency=diff_ssf,
                                                     topology=topology,
                                                     trajectory_id=trajectory_id,
                                                     output_path=output_dir,
                                                     difference=True)
        theta_2_reweighted_figures[chosen_theta]['diff_ssfreq'] = diff_ssfreq_fig

    # Structural Metrics Distributions (Uniform + Reweighted)
    rw_metrics_fig = create_reweighting_metrics_fig(metrics=metrics,
                                                    weight_sets=chosen_weight_sets,
                                                    title_text=f'{trajectory_id} Reweighted '
                                                                'Structural Metrics')

    ##############################################################################################
    ############################# BUILD REWEIGHTING FIGURES HTML DIVS ############################
    ##############################################################################################

    # Build HTML dashboard
    print(f'Building {trajectory_id} reweighting dashboard...')

    # Effective Frames/Chosen theta value
    chosen_theta_div = chosen_thetas_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                 full_html=False,
                                                 include_plotlyjs=False,
                                                 div_id='chosen_thetas')

    # Experimental Fittings
    rw_fits_divs = []
    for rw_fits_fig, data_type in zip(rw_fits_figs, exp_type):      
        if data_type == 'SAXS':
            rw_saxs_fits_div = rw_fits_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                   full_html=False,
                                                   include_plotlyjs=False,
                                                   div_id=f'rw_fits_{data_type}')
            rw_fits_divs.append(rw_saxs_fits_div)
    rw_fits_div_total = ''.join(rw_fits_divs)

    # Build Reweighted Figures Divs
    theta_2_reweighted_divs = {}
    for i,chosen_theta in enumerate(chosen_thetas):
        # Nested dict to split up reweighted divs according to theta values
        theta_2_reweighted_divs[chosen_theta] = {'cmap': None,
                                                 'rw_cmap': None,
                                                 'diff_cmap': None,
                                                 'dmatrix': None,
                                                 'rw_dmatrix': None,
                                                 'diff_dmatrix': None,
                                                 'ssfreq': None,
                                                 'rw_ssfreq': None,
                                                 'diff_ssfreq': None
                                                 }

        # Contact Maps
        ## Uniform
        cmap_fig = theta_2_reweighted_figures[chosen_theta]['cmap']
        cmap_div = cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                    full_html=False,
                                    include_plotlyjs=False,
                                    div_id=f'cmap_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['cmap'] = cmap_div

        ## Reweighted
        rw_cmap_fig = theta_2_reweighted_figures[chosen_theta]['rw_cmap']
        rw_cmap_div = rw_cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                          full_html=False,
                                          include_plotlyjs=False,
                                          div_id=f'rw_cmap_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['rw_cmap'] = rw_cmap_div

        diff_cmap_fig = theta_2_reweighted_figures[chosen_theta]['diff_cmap']
        diff_cmap_div = diff_cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                              full_html=False,
                                              include_plotlyjs=False,
                                              div_id=f'diff_cmap_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['diff_cmap'] = diff_cmap_div

        # Distance Matrices
        ## Uniform
        dmatrix_fig = theta_2_reweighted_figures[chosen_theta]['dmatrix']
        dmatrix_div = dmatrix_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                        full_html=False,
                                        include_plotlyjs=False,
                                        div_id=f'dmatrix_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['dmatrix'] = dmatrix_div

        ## Reweighted
        rw_dmatrix_fig = theta_2_reweighted_figures[chosen_theta]['rw_dmatrix']
        rw_dmatrix_div = rw_dmatrix_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                full_html=False,
                                                include_plotlyjs=False,
                                                div_id=f'rw_dmatrix_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['rw_dmatrix'] = rw_dmatrix_div

        diff_dmatrix_fig = theta_2_reweighted_figures[chosen_theta]['diff_dmatrix']
        diff_dmatrix_div = diff_dmatrix_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                    full_html=False,
                                                    include_plotlyjs=False,
                                                    div_id=f'diff_dmatrix_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['diff_dmatrix'] = diff_dmatrix_div

        # Secondary Structure Frequencies
        ## Uniform
        ssfreq_fig = theta_2_reweighted_figures[chosen_theta]['ssfreq']
        ssfreq_div = ssfreq_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                        full_html=False,
                                        include_plotlyjs=False,
                                        div_id=f'ssfreq_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['ssfreq'] = ssfreq_div

        ## Reweighted
        rw_ssfreq_fig = theta_2_reweighted_figures[chosen_theta]['rw_ssfreq']
        rw_ssfreq_div = rw_ssfreq_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                              full_html=False,
                                              include_plotlyjs=False,
                                              div_id=f'rw_ssfreq_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['rw_ssfreq'] = rw_ssfreq_div

        diff_ssfreq_fig = theta_2_reweighted_figures[chosen_theta]['diff_ssfreq']
        diff_ssfreq_div = diff_ssfreq_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                  full_html=False,
                                                  include_plotlyjs=False,
                                                  div_id=f'diff_ssfreq_{i+1}')
        theta_2_reweighted_divs[chosen_theta]['diff_ssfreq'] = diff_ssfreq_div

    # Structural Metrics Distributions
    rw_metrics_div = rw_metrics_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                            full_html=False,
                                            include_plotlyjs=False,
                                            div_id='rw_metrics')

    ## Build final dashboard divs
    theta_2_dashboard_divs = {}
    for chosen_theta in chosen_thetas:
        theta_2_dashboard_divs[chosen_theta] = {'title': '',
                                                'cmap_div': '',
                                                'dmatrix_div': '',
                                                'ssfreq_div': ''}

        theta_2_dashboard_divs[chosen_theta]['title'] += (f'{trajectory_id} '
                                                          'BME Reweighting '
                                                          '[\u03B8='
                                                          f'{chosen_theta}]')

        # Create Contact Map Div
        cmap_div = theta_2_reweighted_divs[chosen_theta]['cmap']
        rw_cmap_div = theta_2_reweighted_divs[chosen_theta]['rw_cmap']
        diff_cmap_div = theta_2_reweighted_divs[chosen_theta]['diff_cmap']
        cmap_div_str = cmap_div + ('&nbsp;' * 6) \
                       + rw_cmap_div + ('&nbsp;' * 6) \
                       + diff_cmap_div
        theta_2_dashboard_divs[chosen_theta]['cmap_div'] += cmap_div_str

        # Create Distance Matrix Div
        dmatrix_div = theta_2_reweighted_divs[chosen_theta]['dmatrix']
        rw_dmatrix_div = theta_2_reweighted_divs[chosen_theta]['rw_dmatrix']
        diff_dmatrix_div = theta_2_reweighted_divs[chosen_theta]['diff_dmatrix']
        dmatrix_div_str = dmatrix_div + ('&nbsp;' * 6) \
                          + rw_dmatrix_div + ('&nbsp;' * 6) \
                          + diff_dmatrix_div
        theta_2_dashboard_divs[chosen_theta]['dmatrix_div'] += dmatrix_div_str

        # Create Secondary Structure Frequency Div
        ssfreq_div = theta_2_reweighted_divs[chosen_theta]['ssfreq']
        rw_ssfreq_div = theta_2_reweighted_divs[chosen_theta]['rw_ssfreq']
        diff_ssfreq_div = theta_2_reweighted_divs[chosen_theta]['diff_ssfreq']
        ssfreq_div_str = ssfreq_div + ('&nbsp;' * 6) \
                         + rw_ssfreq_div + ('&nbsp;' * 6) \
                         + diff_ssfreq_div
        theta_2_dashboard_divs[chosen_theta]['ssfreq_div'] += ssfreq_div_str

    theta_divs = ''
    for i,chosen_theta in enumerate(chosen_thetas):
        div_str = f'''
            <div class="header-container"><br><br><br>
                <subheader>{theta_2_dashboard_divs[chosen_theta]['title']}</subheader>
            </div>
            <div class="flex-container" id="cmaps_container_{i+1}">
                {theta_2_dashboard_divs[chosen_theta]['cmap_div']}
            </div>
            <div class="flex-container" id="dmatrices_container_{i+1}">
                {theta_2_dashboard_divs[chosen_theta]['dmatrix_div']}
            </div>
            <div class="flex-container" id="ssfreqs_container_{i+1}">
                {theta_2_dashboard_divs[chosen_theta]['ssfreq_div']}
            </div>
            <br><br>
            '''
        theta_divs += div_str

    ##############################################################################################
    ######################### BUILD/SAVE REWEIGHTING FINAL HTML DASHBOARD ########################
    ##############################################################################################

    ## Build dashboard
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

            subheader {{
                padding-top:25px;
                padding-bottom: 25px;
                padding-right: 25px;
                padding-left: 25px;
                font-family: Helvetica;
                font-size: 60px;
                text-align: center;
                background: radial-gradient(circle, rgba(0,38,66,1) 25%, rgba(48,77,109,1) 100%);
                color: white;
                display: inline-block;
                width: auto;
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
            </style>
        </head>
        <body>
            <div class="header-container">
            <header>Ensemblify Reweighting</header>
            </div>
            <br>
            <div class="flex-container">
            {chosen_theta_div}
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            {rw_fits_div_total}
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            {rw_metrics_div}
            </div>
            <br>
            {theta_divs}
            <br>
        </body>
    </html>
    '''

    # Save dashboard
    with open(os.path.join(output_dir,'reweighting_dashboard.html'), 'w',encoding='utf-8') as f:
        f.write(dashboard_html)

    print('Ensemble reweighting has finished. Please refer to the interactive '
          'reweighting_dashboard.html figure for analysis.')
