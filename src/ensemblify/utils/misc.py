"""Auxiliary functions and classes with miscellaneous use."""

# IMPORTS
## Third Party Imports
import numpy as np
import scipy

# CLASSES
class HashableDict(dict):
    """Take a Python dictionary and make it hashable.
    
    Appropriate for when we will NOT ever modify the dictionary after hashing.
    
    Reference:
        https://stackoverflow.com/a/1151686
    """
    def __key(self):
        return tuple((k,self[k]) for k in sorted(self))
    def __hash__(self):
        return hash(self.__key())
    def __eq__(self, other):
        return self.__key() == other.__key()


# FUNCTIONS
def kde(
    data: np.ndarray,
    weights: list | None = None,
    ) -> tuple[np.ndarray,np.ndarray,float]:
    """Create an array with the Kernel Density Estimate (KDE) distribution for a
    given dataset.

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

    """
    # By default weights are uniform
    if weights is None:
        weights = np.full(len(data),
                          1)

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
