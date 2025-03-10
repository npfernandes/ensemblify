"""Auxiliary functions for creating reweighting figures."""

# IMPORTS
## Third Party Imports
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

## Local Imports
from ensemblify.utils import kde

# CONSTANTS
DEFAULT_REWEIGHTING_TRACE_COLOR_PALETTE = ['#E69F00','#56B4E9','#009E73','#F0E442',
                                           '#0072B2','#D55E00','#CC79A7']
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
        'family': 'Arial',
        'color': 'black',
        'size': 34,
    },
    'modebar_remove': ['zoom','pan','select','lasso2d','zoomIn','zoomOut'],
    'xaxis': DEFAULT_AXIS,
    'yaxis': DEFAULT_AXIS,
}

# FUNCTIONS
def create_effective_frames_fit_fig(
    stats: np.ndarray,
    thetas: np.ndarray,
    choices: int | list[int] | None = None,
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
        stats (np.ndarray):
            An array where each row corresponds to a different theta value with columns
            (chi2_before,chi2_after,phi) where:
                chi2_before:
                    The value for the chisquare of fitting the ensemble with uniform
                    weights to the experimental data.
                chi2_after:
                    The value for the chisquare of fitting the reweighted ensemble to
                    the experimental data.
                phi:
                    The fraction of effective frames being used in the reweighted ensemble.
        thetas (np.ndarray):
            Array of values for the theta parameter used when applying BME algorithm.
        choices (int | list[int] | None):
            Theta value(s) chosen for reweighting ensemble, corresponding data points will be
            highlighted in the created Figure. Defaults to None.
        title_text (str | None):
            Title for the created Figure. Defaults to None.
        colors (list[str]):
            Hexcodes for the colors to use for highlighting theta values. Defaults to ['#E69F00',
            '#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].

    Returns:
        go.Figure:
            The created plot, optionally with data points corresponding to highlighted theta
            values in different colors.
    """
    # Setup choices
    if isinstance(choices,int):
        choices = [choices]

    # Setup color palette
    if colors is None:
        colors = DEFAULT_REWEIGHTING_TRACE_COLOR_PALETTE

    # Create Figure
    fig = go.Figure(layout=DEFAULT_LAYOUT)

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
    fig.update_layout(width=1200,
                      height=800,
                      xaxis=dict(title=dict(text='<i>&#x1D719;<sub>eff</sub><i>',
                                            standoff=30),
                                 range=[0,1.02]),
                      yaxis=dict(title=dict(text='Reduced <i>&#967;<sup>2</sup><sub></sub></i>',
                                            standoff=45)),
                      margin=dict(t=75,
                                  l=80,
                                  r=0,
                                  b=0),
                      legend=dict(x=0.02,
                                  y=0.98,
                                  borderwidth=3,
                                  bordercolor='black',
                                  itemsizing='constant'))

    # Add title
    if title_text:
        fig.update_layout(title=dict(text=title_text,
                                     x=0.5,
                                     pad_b=0))
    else:
        fig.update_layout(margin_t=40)

    return fig


def create_reweighting_fits_fig(
    q: np.ndarray,
    i_exp: np.ndarray,
    err: np.ndarray,
    i_prior: np.ndarray,
    i_posts: np.ndarray | list[np.ndarray],
    title_text: str | None = None,
    colors: list[str] = None,
    ) -> go.Figure:
    """Create a multiplot Figure showcasing the differences between uniform and reweighted
    calculated SAXS data, when fit to experimental data.

    Args:
        q (np.ndarray):
            An array with momentum transfer values, common to all SAXS curves being deal with.
        i_exp (np.ndarray):
            An array with experimentally measured SAXS intensities.
        err (np.ndarray):
            An array with the experimental error of the provided experimental SAXS intensities.
        i_prior (np.ndarray):
            An array of SAXS intensities averaged over all the frames of a SAXS data file
            calculated from a conformational ensemble with uniform weights.
        i_posts (np.ndarray | list[np.ndarray]):
            An array or list of arrays of SAXS intensities averaged over all the frames of a SAXS
            data file calculated from a conformational ensemble with the provided set of weights.
        title_text (str, optional):
            A title for the created multiplot Figure. Defaults to None.
        colors (list[str], optional):
            Color to attribute to the plotted prior and posterior traces, in order of input.
            Defaults to ['#E69F00', '#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].

    Returns:
        go.Figure:
            A multiplot Figure containing four plots:
                - the fitting of i_prior and i_post(s) to the experimental SAXS data i_exp.
                - the previous plot in log scale.
                - Kraty plot for i_prior and i_post fitted to experimental data.
                - residuals between i_prior/i_post(s) and i_exp.
    """
    # Setup i_posts
    if isinstance(i_posts,np.ndarray):
        i_posts = [i_posts]

    # Setup color palette
    if colors is None:
        colors = DEFAULT_REWEIGHTING_TRACE_COLOR_PALETTE

    # Create Figure
    fig = make_subplots(rows=2,
                        cols=2,
                        horizontal_spacing=0.25/2,
                        vertical_spacing=0.35/2)

    # Add exp and prior traces
    ## Exp SAXS vs SAXS reweighted (i_post) and uniform (i_prior). No log scale
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

    ## Exp SAXS vs SAXS reweighted and uniform. yy axis log scale
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

        # Exp SAXS vs SAXS reweighted (i_post) and uniform (i_prior). No log scale
        fig.add_trace(go.Scatter(x=q,
                                 y=i_post,
                                 line=dict(width=3,
                                           dash='dash',
                                           color=color),
                                 name='Reweighted',
                                 showlegend=False),
                      row=1,
                      col=1)

        # Exp SAXS vs SAXS reweighted and uniform. yy axis log scale
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

    # Update Figure Layout
    fig.update_layout(**DEFAULT_LAYOUT)
    fig.update_xaxes(**DEFAULT_AXIS)
    fig.update_yaxes(**DEFAULT_AXIS)

    ## Exp SAXS vs SAXS reweighted (I_post) and uniform (I_prior). No log scale
    fig.update_yaxes(row=1,
                     col=1,
                     title_text='Intensity',
                     ticks='')

    ## Exp SAXS vs SAXS reweighted and uniform. yy axis log scale
    fig.update_yaxes(row=1,
                     col=2,
                     title_text='log (Intensity)',
                     type='log',
                     ticks='')

    ## Kratky plots
    fig.update_yaxes(row=2,
                     col=1,
                     title_text='q<sup>2</sup>Intensity',
                     range=[0, np.max((q**2)*i_prior) + 0.1*np.max((q**2)*i_prior)],
                     ticks='',)

    ## Residuals
    fig.update_yaxes(row=2,
                     col=2,
                     title_text='(I<sup>EXP</sup>-I<sup>CALC</sup>)/&#963;')

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
    fig.update_yaxes(title_standoff=0)

    fig.update_xaxes(title_text='q (nm<sup>-1</sup>)',
                     title_standoff=20)

    # Fix yaxis title moving to the left when we try to
    # showticklabels=False by making its ticklabels invisible
    # see https://github.com/plotly/plotly.js/issues/6552
    fig.update_yaxes(row=1,
                     tickfont=dict(color='rgba(0,0,0,0)',
                                   size=1))
    fig.update_yaxes(row=2,
                     col=1,
                     tickfont=dict(color='rgba(0,0,0,0)',
                                   size=1))

    ## Layout
    fig.update_layout(width=1600,
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


def create_reweighting_metrics_fig(
    metrics: pd.DataFrame,
    weight_sets: np.ndarray | list[np.ndarray],
    title_text: str | None = None,
    colors: list[str] = None,
    ) -> go.Figure:
    """Create a Figure with probability distribution plots for calculated structural metrics, using
    uniform or unequal weights.

    Args:
        metrics (pd.DataFrame):
            A DataFrame with the calculated structural metrics, one row per frame in the
            conformational ensemble.
        weight_sets (np.ndarray | list[np.ndarray]):
            An array or list of arrays containing the weights for calculating the probability
            distributions of each structural metric, for each set of weights.
        title_text (str | None):
            Title for the created Figure. Defaults to None.
        colors (list[str]):
            Hexcodes for colors to use for the traces relative to each i_post, in order of input.
            Defaults to ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'].

    Returns:
        go.Figure:
            A Figure plotting the structural metrics distributions for uniformly and unequally
            weighted conformational ensembles.
    """
    # Setup weights
    if isinstance(weight_sets,np.ndarray):
        weight_sets = [weight_sets]

    # Setup color palette
    if colors is None:
        colors = DEFAULT_REWEIGHTING_TRACE_COLOR_PALETTE

    # Setup axis titles
    axis_titles = {'rg': 'R<sub>g</sub>',
                   'dmax': 'D<sub>max</sub>',
                   'eed': 'D<sub>ee</sub>'}

    # Create Figure
    nrows = len(metrics.columns) // 2 + (len(metrics.columns) % 2 > 0)
    fig = make_subplots(rows=nrows,
                        cols=2,
                        horizontal_spacing=0.3/2,
                        vertical_spacing=0.45/2)

    # Iterate over the different metrics
    row_num = 1
    col_num = 1
    for metric in metrics.columns:
        # Calculate KDE for each metric and add it to Figure, along with average
        x, p_x, av, av_stderr = kde(data=metrics[metric])
        uniform_trace = go.Scatter(x=x,
                                   y=p_x,
                                   line=dict(width=4,
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
                                     width=4)),
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
                                 hovertext=f'Avg: {round(av,2)} &plusmn; {round(av_stderr,2)}',
                                 hoverinfo='text',
                                 fill='toself',
                                 name=f'Avg: {round(av,2)} &plusmn; {round(av_stderr,2)}',
                                 opacity=0),
                      row=row_num,
                      col=col_num)

        # Iterate over each set of weights
        x_rews = []
        p_x_rews = []

        for weights,color in zip(weight_sets,colors):
            # Recalculate KDE for each metric and add it to Figure, along with average
            x_rew, p_x_rew, av_rew, av_rew_stderr = kde(data=metrics[metric],
                                                        weights=weights)
            x_rews.append(x_rew)
            p_x_rews.append(p_x_rew)
            reweighted_trace = go.Scatter(x=x_rew,
                                          y=p_x_rew,
                                          line=dict(width=4,
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
                                         width=4)),
                               legend='legend',
                               row=row_num,
                               col=col_num)

            # Allows for hovering the dashed line to get average value
            fig.add_trace(go.Scatter(x=[av_rew],
                                     y=[x for x in
                                        np.arange(0,
                                                  np.interp(av_rew,x_rew,p_x_rew)+0.001,
                                                  0.001)],
                                     mode='markers',
                                     marker_color=reweighted_trace.line.color,
                                     hovertext=(f'Avg: {round(av_rew,2)} &plusmn; '
                                                f'{round(av_rew_stderr,2)}'),
                                     hoverinfo='text',
                                     fill='toself',
                                     name=(f'Avg: {round(av_rew,2)} &plusmn; '
                                           f'{round(av_rew_stderr,2)}'),
                                     opacity=0),
                          row=row_num,
                          col=col_num)

        # Set axis limits
        x_min = min(np.min(x), *[np.min(x_rew) for x_rew in x_rews])
        x_max = max(np.max(x), *[np.max(x_rew) for x_rew in x_rews])
        y_min = 0
        max_p_x_rews = max([np.max(p_x_rew) for p_x_rew in p_x_rews])
        y_max = max(np.max(p_x) + np.max(p_x)*0.1, max_p_x_rews + max_p_x_rews*0.1)

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
    fig.update_layout(**DEFAULT_LAYOUT)
    fig.update_xaxes(**DEFAULT_AXIS)
    fig.update_yaxes(**DEFAULT_AXIS)

    fig.update_xaxes(title_standoff=30)
    fig.update_yaxes(ticks='',
                     title_standoff=0)

    # Fix yaxis title moving to the left when we try to
    # showticklabels=False by making its ticklabels invisible
    # see https://github.com/plotly/plotly.js/issues/6552
    fig.update_yaxes(tickfont=dict(color='rgba(0,0,0,0)',
                                   size=1))

    fig.update_layout(width=1200,
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
