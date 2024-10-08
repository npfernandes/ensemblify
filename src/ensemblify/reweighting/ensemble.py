"""Reweigh a conformational ensemble using experimental data."""

# IMPORTS
## Standard Library Imports
import glob
import math
import os
import shutil

## Third Party Imports
import subprocess
import numpy as np
import pandas as pd
from plotly.offline import get_plotlyjs

## Local Imports
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.analysis import (calculate_contact_matrix,
                                 calculate_distance_matrix,
                                 calculate_ss_frequency,
                                 calculate_metrics_data,
                                 create_contact_map_fig,
                                 create_distance_matrix_fig,
                                 create_ss_frequency_figure)
from ensemblify.conversion import traj2saxs
from ensemblify.reweighting.ensemble_utils import (attempt_read_reweigthing_data,
                                                   process_exp_data,
                                                   correct_exp_error,
                                                   bme_ensemble_reweighting,
                                                   create_effective_frames_fit_fig,
                                                   average_saxs_profiles,
                                                   attempt_read_data,
                                                   create_reweighting_fits_fig,
                                                   create_reweighting_metrics_fig)

# FUNCTIONS
def reweight_ensemble(
    trajectory: str,
    topology: str,
    trajectory_id: str,
    exp_saxs_data: str,
    output_dir: str = None,
    thetas: list[int] = None,
    calculated_cmatrix: pd.DataFrame | str = None,
    calculated_dmatrix: pd.DataFrame | str = None,
    calculated_ss_frequency: pd.DataFrame | str = None,
    calculated_metrics_data: pd.DataFrame | str = None,
    compare_rg: bool = True,
    compare_dmax: bool = True,
    compare_eed: bool = True,
    compare_cmdist: bool = None,
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
            stored. Defaults to current working directory.
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
    if output_dir is None:
        output_dir = os.getcwd()
    elif not os.path.isdir(output_dir):
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
    weights = attempt_read_reweigthing_data(reweighting_output_directory=output_dir,
                                            trajectory_id=trajectory_id)

    if exp_saxs_file is None or calc_saxs_file is None:
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

        # Calculate SAXS data from ensemble
        calc_saxs_file = traj2saxs(trajectory=trajectory,
                                   topology=topology,
                                   trajectory_id=trajectory_id,
                                   exp_saxs_file=exp_saxs_file)

    else:
        print('Found processed experimental SAXS data and calculated SAXS data.')

    if thetas_array is None or stats is None or weights is None:
        # Reweigh ensemble using different theta values
        thetas_array = np.array(thetas)

        print(f'Applying BME reweighting to {trajectory_id} ensemble with different values for '
              'theta parameter...')
        stats, weights = bme_ensemble_reweighting(exp_saxs_file=exp_saxs_file,
                                                  calc_saxs_file=calc_saxs_file,
                                                  thetas=thetas_array,
                                                  output_dir=reweighting_dir)

    else:
        print(f'Found calculated BME reweighting data with theta values: {list(thetas_array)}')

    print(('Please analyse the provided interactive figure (effective_frames_fit.html) and '
           'input the desired value for the theta parameter.'),
           flush=True)

    effective_frames_fit_fig = create_effective_frames_fit_fig(stats=stats,
                                                               thetas=thetas_array,
                                                               title_text=trajectory_id)

    effective_frames_fit_fig.write_html(os.path.join(output_dir,'effective_frames_fit.html'),
                                        config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

    # Capture chosen theta value
    chosen_theta = int(input('Choose theta:'))
    print(f'Chosen theta value: {chosen_theta}.')

    # Plot L curve with chosen theta(s)
    chosen_thetas_fig = create_effective_frames_fit_fig(stats=stats,
                                                        thetas=thetas_array,
                                                        choices=chosen_theta,
                                                        title_text=f'{trajectory_id} BME '
                                                                    'Reweighting')

    # Extract correct set of weights using index for chosen theta
    chosen_weight_set = weights[np.where(thetas_array == chosen_theta)[0][0]]

    # Grab reweighted SAXS data using choice theta
    rw_calc_saxs_file = glob.glob(os.path.join(reweighting_dir,
                                               f'ibme_t{chosen_theta}_*.calc.dat'))[0]

    # Calculate uniform (prior) and  reweighted (posterior) average SAXS profiles
    i_prior, i_post = average_saxs_profiles(exp_saxs_file=exp_saxs_file,
                                            calc_saxs_file=calc_saxs_file,
                                            rw_calc_saxs_file=rw_calc_saxs_file,
                                            weights=chosen_weight_set)

    # print(('Please analyse the provided interactive figure (effective_frames_fit.html) and '
    #        'input the desired value(s) for the theta parameter.\nIf more than one '
    #        'value, please separate them using a comma.'),
    #        flush=True)

    # effective_frames_fit_fig = create_effective_frames_fit_fig(stats=stats,
    #                                                            thetas=thetas_array,
    #                                                            choices=None,
    #                                                            title_text=trajectory_id)

    # effective_frames_fit_fig.write_html(os.path.join(output_dir,'effective_frames_fit.html'),
    #                                     config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

    # # Capture chosen theta values
    # input_choices = input('Choose theta(s):').split(',')
    # print(f'Chosen theta value(s): {",".join(input_choices)}.')
    # chosen_thetas = [int(x) for x in input_choices]

    # # Plot L curve with chosen theta(s)
    # chosen_thetas_fig = create_effective_frames_fit_fig(stats=stats,
    #                                                     thetas=thetas_array,
    #                                                     choices=chosen_thetas,
    #                                                     title_text=f'{trajectory_id} BME Reweighting')

    # # Extract correct set(s) of weights using index(es) for chosen theta(s)
    # choice_idxs = [np.where(thetas_array == x)[0][0] for x in chosen_thetas]
    # chosen_weight_sets = [ weights[i] for i in choice_idxs]

    # # Calculate reweighted SAXS curves
    # print(f'Calculating {trajectory_id} reweighted SAXS curve(s)...')

    # # Calculate reweighted curve for each theta + weights
    # i_posts = []
    # common_i_prior = None
    # for chosen_theta,chosen_weight_set in zip(chosen_thetas,chosen_weight_sets):

    #     # Take reweighted SAXS data using choice theta, and corresponding weights
    #     rw_calc_saxs_file = glob.glob(os.path.join(reweighting_dir,
    #                                                f'ibme_t{chosen_theta}_*.calc.dat'))[0]

    #     # Calculate prior and posterior SAXS profiles
    #     i_prior, i_post = average_saxs_profiles(exp_saxs_file=exp_saxs_file,
    #                                             calc_saxs_file=calc_saxs_file,
    #                                             rw_calc_saxs_file=rw_calc_saxs_file,
    #                                             weights=chosen_weight_set)

    #     common_i_prior = i_prior
    #     i_posts.append(i_post)


    # Calculate Contact Maps
    cmatrix = attempt_read_data(data=calculated_cmatrix,
                                data_msg_tag='cmatrix',
                                calc_fn=calculate_contact_matrix,
                                trajectory=trajectory,
                                topology=topology,
                                output_path=os.path.join(output_dir,
                                                         f'{trajectory_id}_contact_matrix.csv'))

    rw_cmatrix = calculate_contact_matrix(trajectory=trajectory,
                                          topology=topology,
                                          weights=chosen_weight_set,
                                          output_path=os.path.join(output_dir,
                                                                   f'{trajectory_id}_contact_'
                                                                    'matrix_reweighted.csv'))

    diff_cmatrix = rw_cmatrix - cmatrix

    # Calculate Distance Matrices
    dmatrix = attempt_read_data(data=calculated_dmatrix,
                                data_msg_tag='dmatrix',
                                calc_fn=calculate_distance_matrix,
                                trajectory=trajectory,
                                topology=topology,
                                output_path=os.path.join(output_dir,
                                                         f'{trajectory_id}_distance_matrix.csv'))

    rw_dmatrix = calculate_distance_matrix(trajectory=trajectory,
                                           topology=topology,
                                           weights=chosen_weight_set,
                                           output_path=os.path.join(output_dir,
                                                                    f'{trajectory_id}_distance_'
                                                                     'matrix_reweighted.csv'))

    diff_dmatrix = rw_dmatrix - dmatrix

    # Calculate Secondary Structure Assignment
    ssfreq = attempt_read_data(data=calculated_ss_frequency,
                               data_msg_tag='ss_freq',
                               calc_fn=calculate_ss_frequency,
                               trajectory=trajectory,
                               topology=topology,
                               output_path=os.path.join(output_dir,
                                                        f'{trajectory_id}_ss_frequency.csv'))

    rw_ssfreq = calculate_ss_frequency(trajectory=trajectory,
                                       topology=topology,
                                       weights=chosen_weight_set,
                                       output_path=os.path.join(output_dir,
                                                                f'{trajectory_id}_ss_frequency_'
                                                                 'reweighted.csv'))

    diff_ssfreq = rw_ssfreq - ssfreq

    # Calculate Structural Metrics Distributions
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
                                                         f'{trajectory_id}_structural_'
                                                          'metrics.csv'))

    # Create interactive figures
    print(f'Creating {trajectory_id} reweighted interactive figures...')

    # Experimental Fittings
    ## Read experimental data
    q, i_exp, err = np.loadtxt(exp_saxs_file,
                               unpack=True)
    rw_fits_fig = create_reweighting_fits_fig(q=q,
                                              i_exp=i_exp,
                                              err=err,
                                              i_prior=i_prior,
                                              i_posts=i_post,
                                              title_text=f'{trajectory_id} Reweighted Fittings')

    # Contact Maps
    cmap_fig = create_contact_map_fig(contact_matrix=cmatrix,
                                      topology=topology,
                                      trajectory_id=trajectory_id,
                                      output_path=output_dir)

    rw_cmap_fig = create_contact_map_fig(contact_matrix=rw_cmatrix,
                                         topology=topology,
                                         trajectory_id=trajectory_id,
                                         output_path=output_dir,
                                         reweighted=True)

    diff_cmaps_fig = create_contact_map_fig(contact_matrix=diff_cmatrix,
                                            topology=topology,
                                            trajectory_id=trajectory_id,
                                            output_path=output_dir,
                                            difference=True)

    # Distance Matrices
    ## Get maximum of colorbar
    max_data = max(np.max(dmatrix.values),
                   np.max(rw_dmatrix.values))
    max_colorbar = 5*(math.ceil(max_data/5))

    dmatrix_fig = create_distance_matrix_fig(distance_matrix=dmatrix,
                                             topology=topology,
                                             trajectory_id=trajectory_id,
                                             output_path=output_dir,
                                             max_colorbar=max_colorbar)

    rw_dmatrix_fig = create_distance_matrix_fig(distance_matrix=rw_dmatrix,
                                                topology=topology,
                                                trajectory_id=trajectory_id,
                                                output_path=output_dir,
                                                max_colorbar=max_colorbar,
                                                reweighted=True)

    diff_dmatrices_fig = create_distance_matrix_fig(distance_matrix=diff_dmatrix,
                                                    topology=topology,
                                                    trajectory_id=trajectory_id,
                                                    output_path=output_dir,
                                                    difference=True)

    # Secondary Structure Frequencies
    ssfreq_fig = create_ss_frequency_figure(ss_frequency=ssfreq,
                                            topology=topology,
                                            trajectory_id=trajectory_id,
                                            output_path=output_dir)

    rw_ssfreq_fig = create_ss_frequency_figure(ss_frequency=rw_ssfreq,
                                               topology=topology,
                                               trajectory_id=trajectory_id,
                                               output_path=output_dir,
                                               reweighted=True)

    diff_ssfreqs_fig = create_ss_frequency_figure(ss_frequency=diff_ssfreq,
                                                  topology=topology,
                                                  trajectory_id=trajectory_id,
                                                  output_path=output_dir,
                                                  difference=True)

    # Structural Metrics Distributions
    rw_metrics_fig = create_reweighting_metrics_fig(metrics=metrics,
                                                    weight_sets=chosen_weight_set,
                                                    title_text=f'{trajectory_id} Reweighted '
                                                                'Structural Metrics')

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

    # Contact Maps
    cmap_div = cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                full_html=False,
                                include_plotlyjs=False,
                                div_id='cmap')

    rw_cmap_div = rw_cmap_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                      full_html=False,
                                      include_plotlyjs=False,
                                      div_id='rw_cmap')

    diff_cmaps_div = diff_cmaps_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                            full_html=False,
                                            include_plotlyjs=False,
                                            div_id='diff_cmaps')

    # Distance Matrices
    dmatrix_div = dmatrix_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                      full_html=False,
                                      include_plotlyjs=False,
                                      div_id='dmatrix')

    rw_dmatrix_div = rw_dmatrix_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                            full_html=False,
                                            include_plotlyjs=False,
                                            div_id='rw_dmatrix')

    diff_dmatrices_div = diff_dmatrices_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                    full_html=False,
                                                    include_plotlyjs=False,
                                                    div_id='diff_dmatrices')

    # Secondary Structure Frequencies
    ssfreq_div = ssfreq_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                    full_html=False,
                                    include_plotlyjs=False,
                                    div_id='ssfreq')

    rw_ssfreq_div = rw_ssfreq_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                          full_html=False,
                                          include_plotlyjs=False,
                                          div_id='rw_ssfreq')

    diff_ssfreqs_div = diff_ssfreqs_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                full_html=False,
                                                include_plotlyjs=False,
                                                div_id='diff_ssfreqs')

    # Structural Metrics Distributions
    rw_metrics_div = rw_metrics_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                            full_html=False,
                                            include_plotlyjs=False,
                                            div_id='rw_metrics')

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
                {cmap_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {rw_cmap_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {diff_cmaps_div}
            </div>
            <br><br><br>
            <div class="flex-container">
                {dmatrix_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {rw_dmatrix_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {diff_dmatrices_div}
            </div>
            <br><br><br>
            <div class="flex-container">
                {ssfreq_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {rw_ssfreq_div}
                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                {diff_ssfreqs_div}
            </div>
            <br><br><br>
            <div class="flex-container">
                {rw_metrics_div}
            </div>
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

    reweight_ensemble(trajectory='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_trajectory.xtc',
                      topology='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_top.pdb',
                      trajectory_id='Hst5',
                      exp_saxs_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/proteins_plus_saxs/SAXS/bift_Hst5.dat',
                      output_dir='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/reweighting',
                    #   calculated_cmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/Hst5/Hst5_contact_matrix.csv',
                    #   calculated_dmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/Hst5/Hst5_distance_matrix.csv',
                    #   calculated_ss_frequency=None,#'/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/Hst5/Hst5_ss_assignment.csv',
                    #   calculated_metrics_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal/Hst5_structural_metrics.csv'
                    )
