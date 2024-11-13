"""Auxiliary functions and classes with miscellaneous use."""

# IMPORTS
## Third Party Imports
import math
import numpy as np
import scipy

# FUNCTIONS
def kde(
    data: np.ndarray,
    weights: list | None = None,
    ) -> tuple[np.ndarray,np.ndarray,float]:
    """Calculate a Kernel Density Estimate (KDE) distribution for a given dataset.

    Weights for the given dataset can be provided to alter the contribution of each
    data point to the KDE distribution.

    Args:
        data:
            dataset to calculate KDE values for.
        weights:
            array with weights for each data point. Defaults to uniform weights.

    Returns:
        A tuple (x_coords,norm_kde,avg) where:
            x_coords:
                x axis coordinates corresponding to the calculated kde distribution.
            norm_kde:
                normalized kde distribution.
            avg:
                average of the given dataset.
            avg_stderr:
                the standard error of the calculated average.

    Adapted from:
        https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb
    """
    # By default weights are uniform
    if weights is None:
        weights = np.full(len(data),1/len(data))

    # Calculate the probability density function of our data
    # through Kernel Density Estimation using the provided weights
    pdf = scipy.stats.gaussian_kde(data,
                                   bw_method='silverman',
                                   weights=weights)

    # Create 50 point from min to max
    x_coords = np.linspace(start=np.min(data),
                           stop=np.max(data),
                           num=50)

    # Get the KDE distribution of our data by evaluating
    # the calculated pdf over the 50 points
    kde_dist = pdf.evaluate(x_coords)

    # Make sure the kde distribution is normalized
    norm_kde = kde_dist/np.sum(kde_dist)

    # Get the average of our data according to given weights
    avg = np.average(data,
                     weights=weights)

    # Get the standard error of the calculated mean
    avg_stderr = np.std(data, ddof=1) / np.sqrt(np.size(data)) # SEM

    return x_coords,norm_kde,avg,avg_stderr


def get_array_extremum(arrays: list[np.ndarray], get_max: bool | None =True) -> float:
    """Get maximum or minimum value of all elements of all provided arrays.

    Args:
        arrays:
            list of arrays to analyze.
        get_max:
            Whether to get the maximum or minimum (if False) value. Defaults to True.

    Returns:
        ext:
            maximum or minimum value.
    """
    if get_max:
        ext = max(list(map(np.max,arrays)))
    else:
        ext = min(list(map(np.min,arrays)))
    return ext


def round_to_nearest_multiple(n: int, factor: int, up: bool | None = True) -> int:
    """Round a number to the nearest (up or down) multiple of a given factor.

    Args:
        n:
            number to round.
        factor:
            n will be rounded to multiple of this number.
        up: 
            whether to round up or down (if False). Defaults to True.

    Returns:
        rounded:
            rounded number.
    """
    if up:
        rounded = factor*(math.ceil(n/factor))
    else:
        rounded = factor*(math.floor(n/factor))
    return rounded
