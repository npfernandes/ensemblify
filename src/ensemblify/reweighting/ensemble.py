"""Reweight a conformational ensemble using experimental data."""

# IMPORTS
## Standard Library Imports
import glob
import os
import shutil
import subprocess

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
from ensemblify.conversion import traj2saxs
from ensemblify.reweighting.ensemble_utils import (
    attempt_read_data,
    attempt_read_reweighting_data,
    average_saxs_profiles,
    bme_ensemble_reweighting,
    create_effective_frames_fit_fig,
    create_reweighting_fits_fig,
    create_reweighting_metrics_fig,
    correct_exp_error,
    process_exp_data,
)
from ensemblify.utils import get_array_extremum, round_to_nearest_multiple

# FUNCTIONS
def reweight_ensemble(
    trajectory: str,
    topology: str,
    trajectory_id: str,
    exp_saxs_data: str,
    output_dir: str | None = os.getcwd(),
    thetas: list[int] | None = None,
    calculated_cmatrix: pd.DataFrame | str | None = None,
    calculated_dmatrix: pd.DataFrame | str | None = None,
    calculated_ss_frequency: pd.DataFrame | str | None = None,
    calculated_metrics_data: pd.DataFrame | str | None = None,
    compare_rg: bool | None = True,
    compare_dmax: bool | None = True,
    compare_eed: bool | None = True,
    compare_cmdist: bool | None = None,
    ):
    """Apply Bayesian Maximum Entropy (BME) reweighting to a conformational ensemble, given
    experimental SAXS data.

    If calculated metrics data is provided, it will not be recalculated.
    If data for the center mass distance is to be taken from the given calculated metrics data,
    the compare_cmdist mapping must be provided. The identifiers of this mapping will be matched
    to column names of the given DataFrame, if present.
 
    Args:
        trajectory:
            path to .xtc trajectory file where conformational ensemble is stored.
        topology:
            path to .pdb topology file corresponding to any one frame of the ensemble.
        trajectory_id:
            prefix trajectory identifier to be added to plotted traces and output files.
        exp_saxs_data:
            path to .dat file with experimental SAXS data.
        output_dir:
            path to directory where interactive .html plots and reweighting output files will be
            stored. Is created if it does not exist. Defaults to current working directory.
        thetas:
            list of values to try as the theta parameter in BME. The ensemble will be reweighted
            each time using a different theta value. The effect of different theta values can be
            analyzed in the created effective frames figure.
        calculated_cmatrix:
            DataFrame with the calculated average contact matrix for the current trajectory or
            path to this file in .csv format. Defaults to None, and this data is calculated anew.
        calculated_dmatrix:
            DataFrame with the calculated average distance matrix for the current trajectory or
            path to this file in .csv format. Defaults to None, and this data is calculated anew.
        calculated_ss_frequency:
            DataFrame with the calculated secondary structure assignment frequency matrix for the
            current trajectory or path to this file in .csv format. Defaults to None, and this data
            is calculated anew.
        calculated_metrics_data:
            DataFrame with calculated structural metrics (columns) for each frame of the trajectory
            (rows) or path to this DataFrame in .csv format. Defaults to None, and this data is
            calculated anew.
        compare_rg:
            whether to calculate/consider the radius of gyration when comparing structural metrics
            between uniform and reweighted conformational ensembles. Defaults to True.
        compare_dmax:
            whether to calculate/consider the maximum distance between any two alpha carbons when
            comparing structural metrics between uniform and reweighted conformational ensembles.
            Defaults to True.
        compare_eed:
            whether to calculate/consider the distance from the N to C terminal when comparing
            structural metrics between uniform and reweighted conformational ensembles. Defaults
            to True.
        compare_cmdist:
            mapping of identifiers to tuples with two selection strings for creating MDAnalysis
            AtomGroups, whose center mass distance will be calculated. For example:
                {'inter_domain' : ('resid 1:30', 'resid 110:140')}
            If None, no center mass distances are calculated or compared.
            See https://userguide.mdanalysis.org/stable/selections.html for more information about
            MDAnalysis selections.
            Defaults to None.
    """
    # Setup theta values
    if thetas is None:
        thetas = [1, 10, 20, 50, 75, 100, 200, 400, 750, 1000, 5000, 10000]

    # Setup output directory
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Setup directory for reweighting files
    reweighting_dir = os.path.join(output_dir,
                                   'bme_reweighting_results')
    if not os.path.isdir(reweighting_dir):
        os.mkdir(reweighting_dir)

    # Check if we can skip some steps by reading previously computed data from output directory
    exp_saxs_file,\
    calc_saxs_file,\
    thetas_array,\
    stats,\
    weights = attempt_read_reweighting_data(reweighting_output_directory=output_dir,
                                            trajectory_id=trajectory_id)

    if exp_saxs_file is None:
        # Copy experimental data file into output dir before working on it
        exp_saxs_data_copy = os.path.join(output_dir,
                                        f'{trajectory_id}_exp_saxs_input.dat')
        shutil.copy(src=exp_saxs_data,
                    dst=exp_saxs_data_copy)

        # Process input experimental data
        print(f'Processing {trajectory_id} experimental data file...')

        ## Check units
        processed_exp_data = process_exp_data(experimental_data_path=exp_saxs_data_copy)

        ## Correct error
        try:
            exp_saxs_file = correct_exp_error(experimental_data_path=processed_exp_data)
        except subprocess.CalledProcessError: # if BIFT is not available
            exp_saxs_file = os.path.join(os.path.split(processed_exp_data)[0],
                                        f'{trajectory_id}_exp_saxs.dat')
            np.savetxt(exp_saxs_file,
                       np.loadtxt(processed_exp_data),
                       header=' DATA=SAXS')

    if calc_saxs_file is None:
        # Calculate SAXS data from ensemble
        calc_saxs_file = traj2saxs(trajectory=trajectory,
                                   topology=topology,
                                   trajectory_id=trajectory_id,
                                   exp_saxs_file=exp_saxs_file)

    if thetas_array is None or stats is None or weights is None:
        # Reweigh ensemble using different theta values
        thetas_array = np.array(thetas)

        print(f'Applying BME reweighting to {trajectory_id} ensemble with different values for '
              'theta parameter...')
        stats, weights = bme_ensemble_reweighting(exp_saxs_file=exp_saxs_file,
                                                  calc_saxs_file=calc_saxs_file,
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
    input_choices = input('Choose theta(s):').split(',')
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

    ##################################################################### CALCULATE REWEIGHTING FIGURES DATA #####################################################################
    # Calculate prior and posterior average SAXS intensities
    common_i_prior = None
    i_posts = []

    for chosen_theta,chosen_weight_set in zip(chosen_thetas,chosen_weight_sets):
        # Grab reweighted SAXS data using choice theta
        rw_calc_saxs_file = glob.glob(os.path.join(reweighting_dir,
                                                   f'ibme_t{chosen_theta}_*.calc.dat'))[0]

        # Calculate uniform (prior) and  reweighted (posterior) average SAXS profiles
        i_prior, i_post = average_saxs_profiles(exp_saxs_file=exp_saxs_file,
                                                calc_saxs_file=calc_saxs_file,
                                                rw_calc_saxs_file=rw_calc_saxs_file,
                                                weights=chosen_weight_set)

        common_i_prior = i_prior
        i_posts.append(i_post)

    # Calculate Contact Matrices
    ## Uniform
    cmatrix = attempt_read_data(data=calculated_cmatrix,
                                data_msg_tag='cmatrix',
                                calc_fn=calculate_contact_matrix,
                                trajectory=trajectory,
                                topology=topology,
                                output_path=os.path.join(output_dir,
                                                         f'{trajectory_id}_contact_matrix.csv'))
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
    dmatrix = attempt_read_data(data=calculated_dmatrix,
                                data_msg_tag='dmatrix',
                                calc_fn=calculate_distance_matrix,
                                trajectory=trajectory,
                                topology=topology,
                                output_path=os.path.join(output_dir,
                                                         f'{trajectory_id}_distance_matrix.csv'))
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
    ssfreq = attempt_read_data(data=calculated_ss_frequency,
                               data_msg_tag='ss_freq',
                               calc_fn=calculate_ss_frequency,
                               trajectory=trajectory,
                               topology=topology,
                               output_path=os.path.join(output_dir,
                                                        f'{trajectory_id}_ss_frequency.csv'))

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
    metrics = attempt_read_data(data=calculated_metrics_data,
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

    ##################################################################### CREATE REWEIGHTING FIGURES #####################################################################
    # Create interactive figures
    print(f'Creating {trajectory_id} reweighted interactive figures...')

    # Experimental Fittings
    ## Read experimental data
    q, i_exp, err = np.loadtxt(exp_saxs_file,
                               unpack=True)
    ## Create Figure
    rw_fits_fig = create_reweighting_fits_fig(q=q,
                                              i_exp=i_exp,
                                              err=err,
                                              i_prior=common_i_prior,
                                              i_posts=i_posts,
                                              title_text=f'{trajectory_id} Reweighted Fittings')

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
                                             output_path=output_dir,
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
                                                    output_path=output_dir,
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
                                                   output_path=output_dir,
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

    ##################################################################### BUILD REWEIGHTING FIGURES HTML DIVS #####################################################################
    # Build HTML dashboard
    print(f'Building {trajectory_id} reweighting dashboard...')

    # Effective Frames/Chosen theta value
    chosen_theta_div = chosen_thetas_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                 full_html=False,
                                                 include_plotlyjs=False,
                                                 div_id='chosen_thetas')

    # Experimental Fittings
    rw_fits_div = rw_fits_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                      full_html=False,
                                      include_plotlyjs=False,
                                      div_id='rw_fits')

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
        cmap_div_str = cmap_div + ('&nbsp;' * 9) \
                       + rw_cmap_div + ('&nbsp;' * 9) \
                       + diff_cmap_div
        theta_2_dashboard_divs[chosen_theta]['cmap_div'] += cmap_div_str

        # Create Distance Matrix Div
        dmatrix_div = theta_2_reweighted_divs[chosen_theta]['dmatrix']
        rw_dmatrix_div = theta_2_reweighted_divs[chosen_theta]['rw_dmatrix']
        diff_dmatrix_div = theta_2_reweighted_divs[chosen_theta]['diff_dmatrix']
        dmatrix_div_str = dmatrix_div + ('&nbsp;' * 9) \
                          + rw_dmatrix_div + ('&nbsp;' * 9) \
                          + diff_dmatrix_div
        theta_2_dashboard_divs[chosen_theta]['dmatrix_div'] += dmatrix_div_str

        # Create Secondary Structure Frequency Div
        ssfreq_div = theta_2_reweighted_divs[chosen_theta]['ssfreq']
        rw_ssfreq_div = theta_2_reweighted_divs[chosen_theta]['rw_ssfreq']
        diff_ssfreq_div = theta_2_reweighted_divs[chosen_theta]['diff_ssfreq']
        ssfreq_div_str = ssfreq_div + ('&nbsp;' * 9) \
                         + rw_ssfreq_div + ('&nbsp;' * 9) \
                         + diff_ssfreq_div
        theta_2_dashboard_divs[chosen_theta]['ssfreq_div'] += ssfreq_div_str

    theta_divs = ''
    for chosen_theta in chosen_thetas:
        div_str = f'''
            <div>
            <p style="text-align:center; font-family:Helvetica; font-size:48px;">
            {theta_2_dashboard_divs[chosen_theta]['title']}
            </p>
            </div>
            <div class="flex-container">
                {theta_2_dashboard_divs[chosen_theta]['cmap_div']}
            </div>
            <br><br>
            <div class="flex-container">
                {theta_2_dashboard_divs[chosen_theta]['dmatrix_div']}
            </div>
            <br><br>
            <div class="flex-container">
                {theta_2_dashboard_divs[chosen_theta]['ssfreq_div']}
            </div>
            <br><br><br>
            '''
        theta_divs += div_str

    ##################################################################### BUILD/SAVE REWEIGHTING FINAL HTML DASHBOARD #####################################################################
    ## Build dashboard
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
            </style>
        </head>
        <body>
            <br>
            <div class="flex-container">
                {chosen_theta_div}
            </div>
            <br><br><br>
            <div class="flex-container">
                {rw_fits_div}
            </div>
            <br><br><br>
            <div class="flex-container">
                {rw_metrics_div}
            </div>
            <br><br><br>
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

if __name__ == '__main__':
    from ensemblify import update_config

    update_config({'PEPSI_SAXS_PATH': '/home/tiagogomes/software/Pepsi-SAXS',
                   'BIFT_PATH': '/home/tiagogomes/software/bift'})

    # reweight_ensemble(trajectory='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_trajectory.xtc',
    #                   topology='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_top.pdb',
    #                   trajectory_id='Hst5',
    #                   exp_saxs_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/proteins_plus_saxs/SAXS/bift_Hst5.dat',
    #                   output_dir='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/reweighting',
    #                   calculated_cmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/reweighting/Hst5_contact_matrix.csv',
    #                   calculated_dmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/reweighting/Hst5_distance_matrix.csv',
    #                   calculated_ss_frequency='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/reweighting/Hst5_ss_frequency.csv',
    #                   calculated_metrics_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/reweighting/Hst5_structural_metrics.csv'
    #                 )

    reweight_ensemble(trajectory='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/TRAJECTORIES/SMAD4/SMAD4_trajectory.xtc',
                      topology='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/TRAJECTORIES/SMAD4/SMAD4_top.pdb',
                      trajectory_id='SMAD4',
                      exp_saxs_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/proteins_plus_saxs/SAXS/smad4FL_clean_copy.dat',
                      output_dir='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/REWEIGHTING/SMAD4',
                      calculated_cmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/TRAJECTORY_ANALYSIS/SMAD4/SMAD4_contact_matrix.csv',
                      calculated_dmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/TRAJECTORY_ANALYSIS/SMAD4/SMAD4_distance_matrix.csv',
                      #calculated_ss_frequency='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/reweighting/Hst5_ss_frequency.csv',
                      calculated_metrics_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/TRAJECTORY_ANALYSIS/SMAD4/SMAD4_structural_metrics.csv'
                    )
    