"""Auxiliary functions for the reweighting module."""

# IMPORTS
## Standard Library Imports
import contextlib
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

## Third Party Imports
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import sklearn
import subprocess
from plotly.subplots import make_subplots
from tqdm import tqdm

## Local Imports
from ensemblify.config import GLOBAL_CONFIG
from ensemblify.utils import extract_pdb_info,kde
from ensemblify.reweighting.third_party.BME_main import BME
from ensemblify.analysis.trajectory_utils import calculate_ss_frequency

# FUNCTIONS
def process_exp_data(experimental_data_path: str) -> str:
    """Check formatting and units in input experimental data file.

    If values for q are in Ångstrom, they are converted to nanometer.
    Any q-values above 5nm^(-1) are removed, as SAXS calculations are not reliable in that
    range.
    Trajectory ID, the prefix identifier for the output processed experimental data file, is
    set as the text before the first underscore of experimental data path. If no underscores
    are present, it will just be the whole filename.

    Args:
        experimental_data_path:
            path to experimental SAXS data file.
    Returns:
        processed_exp_saxs:
            path to experimental SAXS data file with any applied changes.
    """
    # Read experimental data into array
    exp_saxs_input = np.loadtxt(experimental_data_path)
    assert np.shape(exp_saxs_input)[1] == 3, ('Unable to read file. Make sure the file only '
                                              'contains 3 columns (q,I,sigma) '
                                              'and #commented lines.')

    # Convert from A to nm if required
    if exp_saxs_input[...,0][-1] <  1:
        print('q is in Å units. Converting to nm.')
        exp_saxs_input[...,0] = exp_saxs_input[...,0]*10

    # Remove values from the SAXS profile related to noise or aggregation (> 5 nm^(-1))
    n_unreliable = (exp_saxs_input[...,0] >= 5).sum()
    if n_unreliable > 0:
        print(f'Found {n_unreliable} q-values above 5 nm^(-1). SAXS calculations are not reliable '
               'in that region of the spectrum. Those datapoints will be removed.')
        exp_saxs_input = exp_saxs_input[(exp_saxs_input[...,0] < 5)]

    # Get prefix for output filename
    exp_data_path_components = os.path.split(experimental_data_path)[1].split('_')
    no_exp = True
    for i,s in enumerate(exp_data_path_components):
        if s == 'exp':
            no_exp = False
            trajectory_id = '_'.join([x for x in exp_data_path_components[:i] ])
    if no_exp:
        exp_data_path_components[-1] = exp_data_path_components[-1][:-4]
        trajectory_id = '_'.join([x for x in exp_data_path_components])

    # Save processed data
    processed_exp_saxs = os.path.join(os.path.split(experimental_data_path)[0],
                                         f'{trajectory_id}_exp_saxs_input_processed.dat')
    np.savetxt(processed_exp_saxs, exp_saxs_input)

    return processed_exp_saxs


def correct_exp_error(experimental_data_path: str) -> str:
    """Correct experimental error of input experimental data file using BIFT.

    Bayesian Indirect Fourier Transformation (BIFT) can identify whether the
    experimental error in small-angle scattering data is over- or
    underestimated. The error values are then scaled accordingly.

    Reference:
        Larsen, A.H. and Pedersen, M.C. (2021), Experimental noise in small-angle scattering
        can be assessed using the Bayesian indirect Fourier transformation. J. Appl. Cryst.,
        54: 1281-1289. https://doi.org/10.1107/S1600576721006877

    Args:
        experimental_data_path:
            path to experimental SAXS data file.
    Returns:
        corrected_exp_saxs:
            path to experimental SAXS data file with corrected errors.
    """
    # Deduce working directory
    working_dir, experimental_data_file = os.path.split(experimental_data_path)

    # Prepare input file for BIFT
    input_file = os.path.join(working_dir,'inputfile.dat')
    with open(input_file,'w',encoding='utf-8') as f:
        f.write(f'{experimental_data_file}\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')

    # Run BIFT
    try:
        subprocess.run(f'{GLOBAL_CONFIG["BIFT_PATH"]} < {os.path.split(input_file)[1]}',
                       shell=True,
                       cwd=working_dir,
                       stdout=open(os.path.join(working_dir,'bift.log'),'w',encoding='utf-8'),
                       stderr=subprocess.STDOUT,
                       check=True)
    except subprocess.CalledProcessError as e:
        raise e

    # Get prefix for output filename
    exp_data_file_components = experimental_data_file.split('_')
    no_exp = True
    for i,s in enumerate(exp_data_file_components):
        if s == 'exp':
            no_exp = False
            trajectory_id = '_'.join([x for x in exp_data_file_components[:i] ])
    if no_exp:
        exp_data_file_components[-1] = exp_data_file_components[-1][:-4]
        trajectory_id = '_'.join([x for x in exp_data_file_components])

    # Save rescaled experimental data
    corrected_exp_saxs = os.path.join(working_dir,
                                      f'{trajectory_id}_exp_saxs.dat')
    np.savetxt(corrected_exp_saxs,
               np.loadtxt(os.path.join(working_dir,'rescale.dat')),
               header=' DATA=SAXS')

    # Save scale factor
    scale_factor = np.loadtxt(os.path.join(working_dir,'scale_factor.dat'))[0,1]

    # Save used parameters in log file
    with (open(os.path.join(working_dir,'parameters.dat'),'r',encoding='utf-8-sig') as f,
          open(os.path.join(working_dir,'bift.log'),'a',encoding='utf-8') as t):
        t.write('-------- Input Parameters --------\n')
        t.write(f.read())

    # Remove unnecessary bift output files
    for tmp_file in ['data.dat','dummy.dat','fit.dat','gs_pr.dat','gx_pr.dat','inputfile.dat',
                     'parameters.dat','pr.dat','rescale.dat','scale_factor.dat','st_pr.dat']:
        os.remove(os.path.join(working_dir,tmp_file))

    print( 'Experimental errors on SAXS intensities have been corrected with BIFT using scale '
          f'factor {scale_factor}.')

    return corrected_exp_saxs


def ibme(
    theta: int,
    exp_file: str,
    calc_file: str,
    output_dir: str,
    ) -> tuple[int,tuple[float,float,float],np.ndarray]:
    """Perform the Iterative Bayesian Maximum Entropy (BME) algorithm on calculated SAXS data,
    given a value for the theta parameter.

    The used algorithm is explained in:
        https://github.com/KULL-Centre/BME/blob/main/notebook/example_04.ipynb
    
    Reference:
        Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
        Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
        doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

    Args:
        theta:
            value for the theta parameter to be used in BME algorithm.
        exp_file:
            path to .dat file with experimental SAXS curve.
        calc_file:
            path to .dat file with SAXS curve calculated from an ensemble.
        output_dir:
            path to directory where all the files resulting from the reweighting procedure will be
            stored.

    Returns:
        A tuple (theta, stats, weights) where:
            theta:
                value for the theta parameter used in BME algorithm (same as input).
            stats:
                a tuple (chi2_before,chi2_after,phi) where:
                    chi2_before:
                        the value for the chisquare of fitting the ensemble with uniform
                        weights to the experimental data.
                    chi2_after:
                        the value for the chisquare of fitting the reweighted ensemble to
                        the experimental data.
                    phi:
                        the fraction of effective frames being used in the reweighted ensemble.
            weights:
                an array containing the new weights of the ensemble, one for each frame.

    """
    # Change current working directory
    old_cd = os.getcwd()
    os.chdir(output_dir)

    # Create reweight object
    rew = BME.Reweight(f'ibme_t{theta}')

    # Load files
    rew.load(exp_file=exp_file,
             calc_file=calc_file)

    # Do reweighting
    with contextlib.redirect_stdout(open(os.devnull, 'w',encoding='utf-8')):
        rew.ibme(theta=theta,
                 iterations=25,
                 ftol=0.001)

    # Restore working directory
    os.chdir(old_cd)

    weights = rew.get_ibme_weights()[-1] # get the final weights
    stats = rew.get_ibme_stats()[-1] # get the final stats

    return theta,stats,weights


def bme_ensemble_reweighting(
    exp_saxs_file: str,
    calc_saxs_file: str,
    thetas: list[int],
    output_dir: str,
    ) -> tuple[np.ndarray,np.ndarray]:
    """Perform Bayesian Maximum Entropy (BME) reweighting on calculated SAXS data based on
    experimental SAXS data of the same protein.

    Applies the iterative BME algorithm, explained in:
        https://github.com/KULL-Centre/BME/blob/main/notebook/example_04.ipynb
    The algorithm is applied using different theta values and the results for each value are stored.
    
    Reference:
        Bottaro S, Bengtsen T, Lindorff-Larsen K. Integrating Molecular Simulation and Experimental
        Data: A Bayesian/Maximum Entropy Reweighting Approach. Methods Mol Biol. 2020;2112:219-240.
        doi: 10.1007/978-1-0716-0270-6_15. PMID: 32006288.

    Args:
        exp_saxs_file:
            path to .dat file with experimental SAXS data.
        calc_saxs_file:
            path to .dat file with SAXS data calculated from a conformational ensemble.
        thetas:
            values of theta to try when applying iBME.
        output_dir:
            path to directory where output files from reweighting protocol will be stored.

    Returns:
        A tuple (stats,weights) where:
            stats:
                an array where each row corresponds to a different theta value with columns
                (chi2_before,chi2_after,phi) where:
                    chi2_before:
                        the value for the chisquare of fitting the ensemble with uniform
                        weights to the experimental data.
                    chi2_after:
                        the value for the chisquare of fitting the reweighted ensemble to
                        the experimental data.
                    phi:
                        the fraction of effective frames being used in the reweighted ensemble.
            weights:
                an array where each row corresponds to a different theta value with column
                contains the set of weights of
                the ensemble, one for each frame.
    """
    # Setup variables for parallel computing
    exp_file = os.path.abspath(exp_saxs_file)
    calc_file = os.path.abspath(calc_saxs_file)

    # Parallelize the computation across the theta values
    results = []
    with ProcessPoolExecutor() as ppe:
        futures = [ppe.submit(ibme, theta, exp_file, calc_file, output_dir) for theta in thetas]
        for future in tqdm(as_completed(futures),total=len(thetas),desc='Reweighting ensemble... '):
            results.append(future.result())

    # Extract the results
    results.sort()  # sort by theta (first element)
    thetas, stats, weights = zip(*results)

    # Convert to numpy arrays for easier indexing
    stats = np.array(stats)

    return stats,weights


def create_effective_frames_fit_fig(
    stats: np.ndarray,
    thetas: np.ndarray,
    choices: list[int] | None = None,
    title_text: str | None = None,
    colors: list[str] = None,
    ) -> go.Figure:
    """Create a Figure plotting the fraction of effective frames vs the chisquare value, resulting
    from applying BME using different theta values.
    
    The fraction of effective frames of an ensemble after reweighting is plotted agaisnt the
    chisquare value of the fitting of the data calculated from the reweighted ensemble to the
    experimental data.
    Each data point results from the application of the Bayesian Maximum Entropy (BME) algorithm to
    the calculated+experimental data using different values for the theta parameter.

    Args:
        stats:
            an array where each row corresponds to a different theta value with columns
            (chi2_before,chi2_after,phi) where:
                chi2_before:
                    the value for the chisquare of fitting the ensemble with uniform
                    weights to the experimental data.
                chi2_after:
                    the value for the chisquare of fitting the reweighted ensemble to
                    the experimental data.
                phi:
                    the fraction of effective frames being used in the reweighted ensemble.
        thetas:
            array of values for the theta parameter used when applying BME algorithm.
        choices:
            theta value(s) chosen for reweighting ensemble, corresponding data points will be
            highlighted in the created Figure. Defaults to None.
        title_text:
            title for the created Figure. Defaults to None.
        colors:
            Hexcodes for the colors to use for highlighting theta values. Defaults to ['#E69F00',
            '#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].

    Returns:
        fig:
            the created plot, optionally with data points corresponding to highlighted theta
            values in different colors.
    """
    # Setup color palette
    if colors is None:
        colors = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']

    # Create Figure
    fig = go.Figure()

    # Add data points
    fig.add_trace(go.Scatter(x=stats[...,2],
                             y=stats[...,1],
                             mode='markers',
                             marker=dict(color='black',
                                         size=25),
                             name='',
                             text=[f'Theta: {theta}' for theta in thetas],
                             showlegend=False))

    # Highlight data points related to chosen thetas (if any)
    if choices is not None:
        for choice,color in zip(choices,colors):
            ndx = np.where(thetas == choice)[0][0]
            x_val = stats[..., 2][ndx]
            y_val = stats[..., 1][ndx]
            fig.add_trace(go.Scatter(x=[x_val],
                                     y=[y_val],
                                     mode='markers',
                                     marker=dict(color=color,
                                                 size=25),
                                     name=('&#x1D719;<sub>eff</sub> = '
                                           f'{round(x_val,2)}'
                                           '<br>&#120594;<i><sup>2</sup></i> = '
                                           f'{round(y_val,2)}'
                                           '<br>\u03B8 = '
                                           f'{choice}')))

    # Update Figure layout
    fig.update_layout(plot_bgcolor='#FFFFFF',
                      font=dict(family='Arial',
                                color='black',
                                size=34),
                      modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                      xaxis=dict(title=dict(text='<i>&#x1D719;<sub>eff</sub><i>',
                                            standoff=30),
                                 showgrid=False,
                                 showline=True,
                                 linewidth=4,
                                 linecolor='black',
                                 color='black',
                                 mirror=True,
                                 ticks='outside',
                                 ticklen=20,
                                 tickwidth=5,
                                 range=[0,1.02]),
                      yaxis=dict(title=dict(text='Reduced <i>&#967;<sup>2</sup><sub></sub></i>',
                                            standoff=30),
                                 showgrid=False,
                                 showline=True,
                                 linewidth=5,
                                 linecolor='black',
                                 color='black',
                                 mirror=True,
                                 ticks='outside',
                                 ticklen=20,
                                 tickwidth=5),
                      width=1200,
                      height=800,
                      margin=dict(t=75,l=80,r=0,b=0),
                      legend=dict(x=0.02,
                                  y=0.98,
                                  borderwidth=5,
                                  bordercolor='black',
                                  itemsizing='constant'),)

    # Add title
    if title_text:
        fig.update_layout(title=dict(text=title_text,
                                     x=0.5,
                                     pad_b=0))
    else:
        fig.update_layout(margin_t=40)

    return fig


def average_saxs_profiles(
    exp_saxs_file: str,
    calc_saxs_file: str,
    rw_calc_saxs_file: str,
    weights: np.ndarray,
    ) -> tuple[float,float]:
    """Average the SAXS intensities for uniform and reweighted calculated SAXS data.
    The uniform data is then scaled and offset by linear regression fitting to experimental data.

    Args:
        exp_saxs_file:
            path to .dat file with experimental SAXS data.
        calc_saxs_file:
            path to .dat file with SAXS data calculated from a conformational ensemble.
        rw_calc_saxs_file:
            path to .dat file with SAXS data calculated from a conformational ensemble considering
            the weights (from iBME) for each frame.
        weights:
            array resulting from iBME with weights for each data point. Defaults to uniform weights.

    Returns:
        A tuple i_prior,i_post where:
            i_prior:
                an array of SAXS intensities averaged over all the frames of a SAXS data file
                calculated from a conformational ensemble with uniform weights.
            i_post:
                an array of SAXS intensities averaged over all the frames of a SAXS data file
                calculated from a conformational ensemble with the provided set of weights.
    """
    # Average the uniform and reweighted calculated saxs intensities
    i_prior = np.average(np.loadtxt(calc_saxs_file)[...,1:],
                         axis=0)
    i_post = np.average(np.loadtxt(rw_calc_saxs_file)[...,1:],
                        axis=0,
                        weights=weights)

    # Get experimental intensities and errors
    _, i_exp, err = np.loadtxt(exp_saxs_file,
                               unpack=True)

    # Perform ordinary least squares Linear Regression, fitting the calculated SAXS data to the
    # experimental SAXS data, taking into account the experimental error of each data point
    model = sklearn.linear_model.LinearRegression()
    model.fit(X=i_prior.reshape(-1,1),
              y=i_exp,
              sample_weight=1/(err**2))

    # Adjust uniform saxs profile by linear fitting of scale and offset from experimental data
    a = model.coef_[0]
    b = model.intercept_
    i_prior = a * i_prior + b

    return i_prior,i_post


def create_reweighting_fits_fig(
    q: np.ndarray,
    i_exp: np.ndarray,
    err: np.ndarray,
    i_prior: np.ndarray,
    i_posts: np.ndarray,
    title_text: str | None = None,
    colors: list[str] = None,
    ) -> go.Figure:
    """Create a multiplot Figure showcasing the differences between uniform and reweighted
    calculated SAXS data, when fit to experimental data.

    Args:
        q:
            an array with momentum transfer values, common to all SAXS curves being deal with.
        i_exp:
            an array with experimentally measured SAXS intensities.
        err:
            an array with the experimental error of the provided experimental SAXS intensities.
        i_prior:
            an array of SAXS intensities averaged over all the frames of a SAXS data file
            calculated from a conformational ensemble with uniform weights.
        i_posts:
            a list of arrays of SAXS intensities averaged over all the frames of a SAXS data file
            calculated from a conformational ensemble with the provided set of weights.
        title_text:
            a title for the created multiplot Figure. Defaults to None.
        colors:
            color to attribute to the plotted prior and posterior traces, in order of input.
            Defaults to ['#E69F00', '#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].

    Returns:
        fig:
            a multiplot Figure containing four plots:
                - the fitting of i_prior and i_post(s) to the experimental SAXS data i_exp.
                - the previous plot in log scale.
                - Kraty plot for i_prior and i_post fitted to experimental data.
                - residuals between i_prior/i_post(s) and i_exp.
    """
    # Setup color palette
    if colors is None:
        colors = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']
    
    # Create Figure
    fig = make_subplots(rows=2,
                        cols=2,
                        horizontal_spacing=0.25/2,
                        vertical_spacing=0.35/2)

    # Add exp and prior traces
    ## Exp SAXS vs SAXS reweighted (i_post) and unreweighted (i_prior). No log scale
    fig.add_traces([go.Scatter(x=q,
                               y=i_exp,
                               error_y=go.scatter.ErrorY(array=err,
                                                         color='#A9A9A9',
                                                         width=2),
                               line=dict(width=4,
                                         color='#808080'),
                               opacity=0.5,
                               name='Experimental Data',
                               showlegend=False),
                    go.Scatter(x=q,
                               y=i_prior,
                               line=dict(width=3,
                                         color='black'),
                               name='Uniform',
                               showlegend=False)],
                    rows=1,
                    cols=1)

    ## Exp SAXS vs SAXS reweighted and unreweighted. yy axis log scale
    fig.add_traces([go.Scatter(x=q,
                               y=i_exp,
                               error_y=go.scatter.ErrorY(array=err,
                                                         color='#A9A9A9',
                                                         width=2),
                               line=dict(width=4,
                                         color='#808080'),
                               opacity=0.5,
                               name='Experimental Data',
                               showlegend=False),
                    go.Scatter(x=q,
                               y=i_prior,
                               line=dict(width=3,
                                         color='black'),
                               name='Uniform',
                               showlegend=False)],
                    rows=1,
                    cols=2)

    ## Kratky plots
    fig.add_traces([go.Scatter(x=q,
                               y=(q**2)*i_exp,
                               error_y=go.scatter.ErrorY(array=(q**2)*err,
                                                         color='#A9A9A9',
                                                         width=2),
                               line=dict(width=4,
                                         color='#808080'),
                               opacity=0.5,
                               name='Experimental Data',
                               showlegend=False),
                    go.Scatter(x=q,
                               y=(q**2)*i_prior,
                               line=dict(width=3,
                                         color='black'),
                               name='Uniform',
                               showlegend=False)],
                    rows=2,
                    cols=1)

    ## Residuals
    fig.add_trace(go.Scatter(x=q,
                             y=(i_exp-i_prior)/err,
                             line=dict(width=3,
                                       color='black'),
                             name='Uniform',
                             opacity=0.8,
                             showlegend=False),
                  row=2,
                  col=2)

    # Add post traces
    for i_post,color in zip(i_posts,colors):

        # Exp SAXS vs SAXS reweighted (i_post) and unreweighted (i_prior). No log scale
        fig.add_trace(go.Scatter(x=q,
                                 y=i_post,
                                 line=dict(width=3,
                                           dash='dash',
                                           color=color),
                                 name='Reweighted',
                                 showlegend=False),
                      row=1,
                      col=1)

        # Exp SAXS vs SAXS reweighted and unreweighted. yy axis log scale
        fig.add_trace(go.Scatter(x=q,
                                 y=i_post,
                                 line=dict(width=3,
                                           dash='dash',
                                           color=color),
                                 name='Reweighted',
                                 showlegend=False),
                      row=1,
                      col=2)

        # Kratky plots
        fig.add_trace(go.Scatter(x=q,
                                 y=(q**2)*i_post,
                                 line=dict(width=3,
                                           dash='dash',
                                           color=color),
                                 name='Reweighted',
                                 showlegend=False),
                      row=2,
                      col=1)

        # Residuals
        fig.add_trace(go.Scatter(x=q,
                                 y=(i_exp-i_post)/err,
                                 line=dict(width=3,
                                           color=color),
                                 opacity=0.8,
                                 name='Reweighted',
                                 showlegend=False),
                       row=2,
                       col=2)

    # Update axes and layout
    ## Exp SAXS vs SAXS reweighted (I_post) and unreweighted (I_prior). No log scale
    fig.update_yaxes(row=1,
                     col=1,
                     title_text='Intensity',
                     ticks='')

    ## Exp SAXS vs SAXS reweighted and unreweighted. yy axis log scale
    fig.update_yaxes(row=1,
                     col=2,
                     title_text='log (Intensity)',
                     type='log',
                     ticks='')

    ## Kratky plots
    fig.update_yaxes(row=2,
                     col=1,
                     title_text='q<sup>2</sup>Intensity',
                     range=[0,
                            np.max((q**2)*i_prior) + 0.1*np.max((q**2)*i_prior)],
                     ticks='outside',
                     ticklen=20,
                     tickwidth=4,)

    ## Residuals
    fig.update_yaxes(row=2,
                     col=2,
                     title_text='(I<sup>EXP</sup>-I<sup>CALC</sup>)/&#963;',
                     ticks='outside',
                     ticklen=20,
                     tickwidth=4,)

    fig.add_shape(row=2,
                  col=2,
                  name='zeroline',
                  type='line',
                  x0=np.min(q),
                  x1=np.max(q),
                  y0=0,
                  y1=0,
                  line=dict(dash='dash',
                            color='black',
                            width=3),
                  opacity=0.6)

    ## All
    fig.update_yaxes(showline=True,
                     linewidth=4,
                     linecolor='black',
                     color='black',
                     mirror=True,
                     title_standoff=0)

    fig.update_xaxes(title_text='q (nm<sup>-1</sup>)',
                     ticks='outside',
                     ticklen=20,
                     tickwidth=4,
                     showline=True,
                     linewidth=4,
                     linecolor='black',
                     color='black',
                     mirror=True,
                     title_standoff=20)

    # Fix yaxis title moving to the left
    # see https://github.com/plotly/plotly.js/issues/6552
    fig.update_yaxes(row=1,
                     tickfont=dict(color='rgba(0,0,0,0)',
                                   size=1))

    ## Layout
    fig.update_layout(plot_bgcolor='#FFFFFF',
                      font=dict(family='Arial',
                                color='black',
                                size=34),
                      modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                      width=1600,
                      height=1200,
                      margin=dict(t=80,
                                  b=0,
                                  r=0,
                                  l=160))
    if title_text:
        fig.update_layout(title=dict(text=title_text,
                                     x=0.5))
    else:
        fig.update_layout(margin_t=40)

    return fig


def create_ss_frequency_difference_fig(
    ss_assignment: pd.DataFrame | str,
    topology: str,
    weights: np.ndarray | None = None,
    trajectory_id: str | None = None,
    output_path: str | None = None,
    ) -> go.Figure:
    """Create a secondary structure frequency Figure from a secondary structure assignment matrix.

    Args:
        ss_assignment:
            calculated secondary structure assignment matrix DataFrame or path to calculated matrix
            in .csv format.
        topology:
            path to topology .pdb file.
        trajectory_id:
            used on Figure title and prefix for saved ss_frequency filename. Defaults to None.
        output_path:
            path to output .html file or output directory where created Figure will be stored.
            If directory, written file is named 'ss_frequency.html', optionally with
            trajectory_id prefix. Defaults to None.

    Returns:
        ss_freq_fig:
            a stacked line plot with the secondary structure frequencies of each secondary structure
            type for each residue in the structure.
    """
    if isinstance(ss_assignment,str):
        assert ss_assignment.endswith('.csv'), ('Secondary structure assignment '
                                                'matrix must be in .csv format!')
        ss_assignment = pd.read_csv(ss_assignment,index_col=0)

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

    # Count the frequency of each secondary structure element
    uniform_frequency = calculate_ss_frequency(ss_assignment=ss_assignment,
                                               weights=None)

    # Count the frequency of each secondary structure element with reweighting
    reweighted_frequency = calculate_ss_frequency(ss_assignment=ss_assignment,
                                                  weights=weights)

    # Calculate the frequency difference
    difference_frequency = uniform_frequency - reweighted_frequency

    # Create Figure
    ss_freq_fig = go.Figure()

    # Adding traces for each secondary structure type
    colors = ['#1f77b4',  # Blue
              '#2ca02c',  # Green
              '#d62728'  # Red
              ]

    for structure,color in zip(difference_frequency.index,colors):
        ss_freq_fig.add_trace(go.Scatter(x=list(range(difference_frequency.columns[0],
                                                      difference_frequency.columns[-1]+1)),
                                        y=difference_frequency.loc[structure],
                                        mode='lines',
                                        marker_color=color,
                                        marker_size=50,
                                        line_width=2,
                                        name=structure,
                                        hoverinfo='text',
                                        hovertext=[ f'{x_label}, \
                                                   {difference_frequency.loc[structure].iloc[i]}' \
                                                   for i,x_label in enumerate(x_labels) ]))

    # Setup chain dividers lines
    num_res = len(difference_frequency.columns)
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
                                      width=0.5)))

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
    if trajectory_id is not None:
        ss_freq_title = f'{trajectory_id} Secondary Structure Frequencies'
    else:
        ss_freq_title = 'Secondary Structure Frequencies'
    ss_freq_fig.update_layout(width=1000,
                              height=750,
                              font=dict(family='Arial',
                                        color='black',
                                        size=18),
                              plot_bgcolor = '#FFFFFF',
                              paper_bgcolor = '#FFFFFF',
                              modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                              title=dict(text=ss_freq_title,
                                         x=0.5),
                              xaxis=dict(title='Residue',
                                         ticks='outside',
                                         tickvals=tickvals,
                                         ticktext=x_text,
                                         showgrid=True),
                              yaxis=dict(title='Frequency',
                                         ticks='outside',
                                         range=[0,1],
                                         showgrid=False,
                                         title_standoff=5),
                              shapes=shapes)

    if output_path is not None:
        # Save Secondary Structure frequency
        if os.path.isdir(output_path):
            if trajectory_id is not None:
                output_filename = f'{trajectory_id}_ss_frequency.html'
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


def create_reweighting_metrics_fig(
    metrics: pd.DataFrame,
    rw_weights: list[np.ndarray],
    title_text: str | None = None,
    colors: list[str] = None,
    ) -> go.Figure:
    """Create a Figure with probability distribution plots for calculated structural metrics, using
    uniform or unequal weights.

    Args:
        metrics:
            a DataFrame with the calculated structural metrics, one row per frame in the
            conformational ensemble.
        rw_weights:
            a list of arrays containing the weights for calculating the probability distributions of
            each structural metric, for each set of weights.
        title_text:
            title for the created Figure. Defaults to None.
        colors:
            hexcodes for colors to use for the traces relative to each i_post, in order of input.
            Defaults to ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].

    Returns:
        fig:
            a Figure plotting the structural metrics distributions for uniformly and unequally
            weighted conformational ensembles.
    """
    # Setup color palette
    if colors is None:
        colors = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']
    
    # Setup axis titles
    axis_titles = {'rg': 'R<sub>g</sub>',
                   'dmax': 'D<sub>max</sub>',
                   'eed': 'D<sub>ee</sub>'}

    # Create Figure
    nrows = math.ceil(len(metrics.columns)/2)
    fig = make_subplots(rows=nrows,
                        cols=2,
                        horizontal_spacing=0.35/2,
                        vertical_spacing=0.45/2)

    # Iterate over the different metrics
    row_num = 1
    col_num = 1
    for metric in metrics.columns:
        # Calculate KDE for each metric and add it to Figure, along with average
        x, p_x, av = kde(data=metrics[metric])
        uniform_trace = go.Scatter(x=x,
                                   y=p_x,
                                   line=dict(width=3,
                                             color='black'),
                                   name=f'{metric}_uniform')
        fig.add_trace(uniform_trace,
                      row=row_num,
                      col=col_num)

        fig.add_shape(dict(name=uniform_trace.name,
                           type='line',
                           xref='x',
                           x0=av,
                           x1=av,
                           y0=0,
                           y1=np.interp(av,x,p_x),
                           line=dict(dash='dash',
                                     color=uniform_trace.line.color,
                                     width=3)),
                           legend='legend',
                           row=row_num,
                           col=col_num)

        # Allows for hovering the dashed line to get average value
        fig.add_trace(go.Scatter(x=[av],
                                 y=[x for x in np.arange(0,
                                                         np.interp(av,x,p_x)+0.001,
                                                         0.001)],
                                 mode='markers',
                                 marker_color='black',
                                 hovertext=f'Avg: {round(av,2)}',
                                 hoverinfo='text',
                                 fill='toself',
                                 name=f'Avg: {round(av,2)}',
                                 opacity=0),
                      row=row_num,
                      col=col_num)

        # Iterate over each set of weights
        x_rews = []
        p_x_rews = []

        for weights,color in zip(rw_weights,colors):
            # Recalculate KDE for each metric and add it to Figure, along with average
            x_rew, p_x_rew, av_rew = kde(data=metrics[metric],
                                         weights=weights)
            x_rews.append(x_rew)
            p_x_rews.append(p_x_rew)
            reweighted_trace = go.Scatter(x=x_rew,
                                          y=p_x_rew,
                                          line=dict(width=3,
                                                    color=color),
                                          name=f'{metric}_reweighted')
            fig.add_trace(reweighted_trace,
                          row=row_num,
                          col=col_num)

            fig.add_shape(dict(name=reweighted_trace.name,
                               type='line',
                               xref='x',
                               x0=av_rew,
                               x1=av_rew,
                               y0=0,
                               y1=np.interp(av_rew,x_rew,p_x_rew),
                               line=dict(dash='dash',
                                         color=reweighted_trace.line.color,
                                         width=3)),
                               legend='legend',
                               row=row_num,
                               col=col_num)

            # Allows for hovering the dashed line to get average value
            fig.add_trace(go.Scatter(x=[av_rew],
                                     y=[x for x in np.arange(0,
                                                             np.interp(av_rew,x_rew,p_x_rew)+0.001,
                                                             0.001)],
                                     mode='markers',
                                     marker_color=reweighted_trace.line.color,
                                     hovertext=f'Avg: {round(av_rew,2)}',
                                     hoverinfo='text',
                                     fill='toself',
                                     name=f'Avg: {round(av_rew,2)}',
                                     opacity=0),
                          row=row_num,
                          col=col_num)

        # Set axis limits
        x_min = min(np.min(x),
                    *[np.min(x_rew) for x_rew in x_rews])
        x_max = max(np.max(x),
                    *[np.max(x_rew) for x_rew in x_rews])
        y_min = 0
        max_p_x_rews = max([np.max(p_x_rew) for p_x_rew in p_x_rews])
        y_max = max(np.max(p_x)+np.max(p_x)*0.1,
                    max_p_x_rews+max_p_x_rews*0.1)

        # Set axis titles
        try:
            x_title = f'{axis_titles[metric]} (&#197;)'
        except KeyError:
            x_title = f'<i>{metric}</i> (&#197;)'
        try:
            y_title = f'KDE({axis_titles[metric]})'
        except KeyError:
            y_title = f'KDE(<i>{metric}</i>)'

        fig.update_xaxes(title_text=x_title,
                         range=[x_min,x_max],
                         row=row_num,
                         col=col_num)
        fig.update_yaxes(title_text=y_title,
                         range=[y_min,y_max],
                         row=row_num,
                         col=col_num)

        # Go back to first col if changing rows
        if col_num == 2:
            col_num = 1
            row_num += 1
        else:
            col_num += 1

    # Update Figure layout
    fig.update_yaxes(ticks='outside',
                     ticklen=20,
                     tickwidth=4,
                     showline=True,
                     linewidth=4,
                     linecolor='black',
                     color='black',
                     mirror=True,
                     title_standoff=0)

    fig.update_xaxes(ticks='outside',
                     ticklen=20,
                     tickwidth=4,
                     showline=True,
                     linewidth=4,
                     linecolor='black',
                     color='black',
                     mirror=True,
                     title_standoff=20)

    fig.update_layout(plot_bgcolor='#FFFFFF',
                      font=dict(family='Arial',
                                color='black',
                                size=34),
                      modebar_remove=['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
                      width=1400,
                      height=1000,
                      margin=dict(t=80,
                                  b=0,
                                  r=0,
                                  l=180),
                      showlegend=False)

    # Add Figure title
    if title_text:
        fig.update_layout(title=dict(text=title_text,
                                     x=0.5,
                                     xanchor='center',
                                     pad_b=0))
    else:
        fig.update_layout(margin_t=40)

    return fig
