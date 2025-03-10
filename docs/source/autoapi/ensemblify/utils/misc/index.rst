ensemblify.utils.misc
=====================

.. py:module:: ensemblify.utils.misc

.. autoapi-nested-parse::

   Auxiliary functions and classes with miscellaneous use.



Functions
---------

.. autoapisummary::

   ensemblify.utils.misc.kde
   ensemblify.utils.misc.get_array_extremum
   ensemblify.utils.misc.round_to_nearest_multiple


Module Contents
---------------

.. py:function:: kde(data: numpy.ndarray, weights: list | None = None) -> tuple[numpy.ndarray, numpy.ndarray, float]

   Calculate a Kernel Density Estimate (KDE) distribution for a given dataset.

   Weights for the given dataset can be provided to alter the contribution of each
   data point to the KDE distribution.

   :param data: Dataset to calculate KDE values for.
   :type data: :py:class:`np.ndarray`
   :param weights: Array with weights for each data point. Defaults to uniform weights.
   :type weights: :py:class:`list, optional`

   :returns:

                 x_coords (np.ndarray):
                     X axis coordinates corresponding to the calculated kde distribution.
                 norm_kde (np.ndarray):
                     Normalized kde distribution.
                 weighted_average (float):
                     Weighted average of the given dataset.
                 weighted_standard_error (float):
                     Weighted standard error of the calculated average.
   :rtype: tuple[np.ndarray,np.ndarray,float,float]

   Adapted from:
       https://github.com/FrPsc/EnsembleLab/blob/main/EnsembleLab.ipynb

   Reference for standard error of weighted average:
       https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf
       Copyright © 2006 Prof. James Kirchner


.. py:function:: get_array_extremum(arrays: list[numpy.ndarray], get_max: bool | None = True) -> float

   Get maximum or minimum value of the set with all values of all provided arrays.

   :param arrays: List of arrays to analyze.
   :type arrays: :py:class:`list[np.ndarray]`
   :param get_max: Whether to get the maximum or minimum (if False) value. Defaults to True.
   :type get_max: :py:class:`bool, optional`

   :returns:     Maximum or minimum value.
   :rtype: float


.. py:function:: round_to_nearest_multiple(n: int, factor: int, up: bool | None = True) -> int

   Round a number to the nearest (up or down) multiple of a given factor.

   :param n: Number to round.
   :type n: :py:class:`int`
   :param factor: The number n will be rounded to a multiple of factor.
   :type factor: :py:class:`int`
   :param up: Whether to round up or down (if False). Defaults to True.
   :type up: :py:class:`bool, optional`

   :returns:     rounded number.
   :rtype: rounded


