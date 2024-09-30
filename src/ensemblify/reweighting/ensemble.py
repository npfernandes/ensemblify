"""Reweigh a conformational ensemble using experimental data."""

# IMPORTS
## Standard Library Imports
import glob
import os
import shutil

## Third Party Imports
import subprocess
import numpy as np
import pandas as pd
from plotly.offline import get_plotlyjs

## Local Imports
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.analysis import calculate_metrics_data
from ensemblify.conversion import traj2saxs
from ensemblify.reweighting.ensemble_utils import (process_exp_data,correct_exp_error,bme_ensemble_reweighting,
                                                   create_effective_frames_fit_fig,average_saxs_profiles,
                                                   create_reweighting_fits_fig,create_ss_frequency_difference_fig,
                                                   create_reweighting_metrics_fig)

# FUNCTIONS
def reweigh_ensemble(
    trajectory: str,
    topology: str,
    trajectory_id: str,
    exp_saxs_data: str,
    output_dir: str = None,
    thetas: list[int] = None,
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
        calculated_metrics_data:
            DataFrame with calculated structural metrics (columns) for each frame of the trajectory
            (rows) or path to this DataFrame in .csv format. Defaults to None, and this data is calculated anew.
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

    # Copy experimental data file into output dir before working on it
    exp_saxs_data_copy = os.path.join(output_dir,
                                      f'{trajectory_id}_exp_saxs_input.dat')
    shutil.copy(src=exp_saxs_data,
                dst=exp_saxs_data_copy)

    # Setup directory for reweighting files
    reweighting_dir = os.path.join(output_dir,
                                   'reweighting_results')
    if not os.path.isdir(reweighting_dir):
        os.mkdir(reweighting_dir)

    # Process input experimental data (check units)
    processed_exp_data = process_exp_data(experimental_data_path=exp_saxs_data_copy)

    # Process input experimental data (correct error)
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

    # Calculate metrics from trajectory
    if isinstance(calculated_metrics_data,str):
        assert calculated_metrics_data.endswith('.csv'), ('Calculated metrics data must be'
                                                          ' provided in .csv format!')
        metrics = pd.read_csv(calculated_metrics_data,index_col=0)
        print('Ensemble structural metrics data has been read from file.')
    elif isinstance(calculated_metrics_data,pd.DataFrame):
        print('Ensemble structural metrics data has been provided.')
        metrics = calculated_metrics_data
    else:
        print('No data was provided for ensemble structural metrics.')
        print('Calculating ensemble structural metrics...')
        metrics = calculate_metrics_data(trajectory=trajectory,
                                         topology=topology,
                                         output_path=os.path.join(output_dir,
                                                                  f'{trajectory_id}_structural_metrics.csv'),
                                         rg=compare_rg,
                                         dmax=compare_dmax,
                                         eed=compare_eed,
                                         cm_dist=compare_cmdist)

    # Reweigh ensemble using different theta values
    thetas_array = np.array(thetas)

    print('Reweighting ensemble with different values for theta parameter...')
    stats, weights = bme_ensemble_reweighting(exp_saxs_file=exp_saxs_file,
                                              calc_saxs_file=calc_saxs_file,
                                              thetas=thetas_array,
                                              output_dir=reweighting_dir)

    print(('Please analyse the provided interactive figure (effective_frames_fit.html) and '
           'input the desired value(s) for the theta parameter.\nIf more than one '
           'value, please separate them using a comma.'),
           flush=True)

    effective_frames_fit_fig = create_effective_frames_fit_fig(stats=stats,
                                                               thetas=thetas_array,
                                                               choices=None,
                                                               title_text=trajectory_id)

    effective_frames_fit_fig.write_html(os.path.join(output_dir,'effective_frames_fit.html'),
                                        config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'])

    # Capture chosen theta values
    input_choices = input('Chosen theta(s):').split(',')
    print(f'Chosen theta value(s): {",".join(input_choices)}.')
    choices = [int(x) for x in input_choices]

    # Index(es) for chosen theta(s), used to extract correct set(s) of weights
    choice_idxs = [np.where(thetas_array == x)[0][0] for x in choices]

    # Plot L curve with chosen theta(s)
    chosen_thetas_fig = create_effective_frames_fit_fig(stats=stats,
                                                        thetas=thetas_array,
                                                        choices=choices,
                                                        title_text=f'{trajectory_id} BME Reweighting')

    # Take experimental data
    q, i_exp, err = np.loadtxt(exp_saxs_file,
                               unpack=True)

    rw_weights = []
    i_posts = []
    i_priors = []

    print(f'Reweighting {trajectory_id} ensemble analysis...')

    for choice,choice_idx in zip(choices,choice_idxs):

        # Take reweighted SAXS data using choice theta, and corresponding weights
        rw_calc_saxs_file = glob.glob(os.path.join(reweighting_dir,
                                                   f'ibme_t{choice}_*.calc.dat'))
        choice_weights = weights[choice_idx]
        rw_weights.append(choice_weights)

        # Calculate prior and posterior SAXS profiles
        i_prior, i_post = average_saxs_profiles(exp_saxs_file=exp_saxs_file,
                                                calc_saxs_file=calc_saxs_file,
                                                rw_calc_saxs_file=rw_calc_saxs_file,
                                                weights=choice_weights)
        i_priors.append(i_prior)
        i_posts.append(i_post)

    # All the priors are the same
    common_i_prior = i_priors[0]

    # Plot uniform/reweighted fits
    rw_fits_fig = create_reweighting_fits_fig(q=q,
                                              i_exp=i_exp,
                                              err=err,
                                              i_prior=common_i_prior,
                                              i_posts=i_posts,
                                              title_text=f'{trajectory_id} Reweighted Fittings')

    # Compare contact maps between uniform/reweighted ensembles
    #rw_cmaps_fig = create_reweighting_contact_maps_fig()

    # Compare metrics between uniform/reweighted ensembles
    rw_metrics_fig = create_reweighting_metrics_fig(metrics=metrics,
                                                    rw_weights=rw_weights,
                                                    title_text=f'{trajectory_id} Reweighted Structural Metrics')

    # Offline Interactive Html
    chosen_theta_div = chosen_thetas_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                                 full_html=False,
                                                 include_plotlyjs=False,
                                                 div_id='chosen_thetas')

    rw_fits_div = rw_fits_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                      full_html=False,
                                      include_plotlyjs=False,
                                      div_id='rw_fits')

    rw_metrics_div = rw_metrics_fig.to_html(config=GLOBAL_CONFIG['PLOTLY_DISPLAY_CONFIG'],
                                            full_html=False,
                                            include_plotlyjs=False,
                                            div_id='rw_metrics')

    # Build html
    dashboard_html = f'''
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
        </body>
    </html>
    '''

    # Save dashboard
    with open(os.path.join(output_dir,'reweighting_dashboard.html'), 'w',encoding='utf-8') as f:
        f.write(dashboard_html)

    print('Ensemble reweighting has finished. Please refer to the interactive '
          'reweighting_dashboard.html figure for analysis.')
