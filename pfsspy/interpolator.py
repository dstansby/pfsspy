"""
Interpolator functions

This module contains a version of scipy.interpolate.RegularGridInterpolator,
which has been edited for performance.

THE CODE HERE IS LIABLE TO CHANGE/BREAK AT ANY TIME. Do not use this code
outside of pfsspy.
"""
import itertools
import numpy as np


class RegularGridInterpolator(object):
    """
    Interpolation on a regular grid in arbitrary dimensions
    The data must be defined on a regular grid; the grid spacing however may be
    uneven.  Linear interpolation is performed.
    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, )
        The points defining the regular grid in n dimensions.
    values : array_like, shape (m1, ..., mn, ...)
        The data on the regular grid in n dimensions.
    fill_value : number, optional
        If provided, the value to use for points outside of the
        interpolation domain. If None, values outside
        the domain are extrapolated.
    Methods
    -------
    __call__
    Notes
    -----
    Contrary to LinearNDInterpolator and NearestNDInterpolator, this class
    avoids expensive triangulation of the input data by taking advantage of the
    regular grid structure.
    If any of `points` have a dimension of size 1, linear interpolation will
    return an array of `nan` values. Nearest-neighbor interpolation will work
    as usual in this case.
    .. versionadded:: 0.14
    Examples
    --------
    Evaluate a simple example function on the points of a 3D grid:
    >>> from scipy.interpolate import RegularGridInterpolator
    >>> def f(x, y, z):
    ...     return 2 * x**3 + 3 * y**2 - z
    >>> x = np.linspace(1, 4, 11)
    >>> y = np.linspace(4, 7, 22)
    >>> z = np.linspace(7, 9, 33)
    >>> data = f(*np.meshgrid(x, y, z, indexing='ij', sparse=True))
    ``data`` is now a 3D array with ``data[i,j,k] = f(x[i], y[j], z[k])``.
    Next, define an interpolating function from this data:
    >>> my_interpolating_function = RegularGridInterpolator((x, y, z), data)
    Evaluate the interpolating function at the two points
    ``(x,y,z) = (2.1, 6.2, 8.3)`` and ``(3.3, 5.2, 7.1)``:
    >>> pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])
    >>> my_interpolating_function(pts)
    array([ 125.80469388,  146.30069388])
    which is indeed a close approximation to
    ``[f(2.1, 6.2, 8.3), f(3.3, 5.2, 7.1)]``.
    See also
    --------
    NearestNDInterpolator : Nearest neighbour interpolation on unstructured
                            data in N dimensions
    LinearNDInterpolator : Piecewise linear interpolant on unstructured data
                           in N dimensions
    References
    ----------
    .. [1] Python package *regulargrid* by Johannes Buchner, see
           https://pypi.python.org/pypi/regulargrid/
    .. [2] Wikipedia, "Trilinear interpolation",
           https://en.wikipedia.org/wiki/Trilinear_interpolation
    .. [3] Weiser, Alan, and Sergio E. Zarantonello. "A note on piecewise linear
           and multilinear table interpolation in many dimensions." MATH.
           COMPUT. 50.181 (1988): 189-196.
           https://www.ams.org/journals/mcom/1988-50-181/S0025-5718-1988-0917826-0/S0025-5718-1988-0917826-0.pdf
    """
    # this class is based on code originally programmed by Johannes Buchner,
    # see https://github.com/JohannesBuchner/regulargrid

    def __init__(self, points, values,
                 fill_value=np.nan):

        if not hasattr(values, 'ndim'):
            # allow reasonable duck-typed values
            values = np.asarray(values)

        if len(points) > values.ndim:
            raise ValueError("There are %d point arrays, but values has %d "
                             "dimensions" % (len(points), values.ndim))

        if hasattr(values, 'dtype') and hasattr(values, 'astype'):
            if not np.issubdtype(values.dtype, np.inexact):
                values = values.astype(float)

        self.fill_value = fill_value
        if fill_value is not None:
            fill_value_dtype = np.asarray(fill_value).dtype
            if (hasattr(values, 'dtype') and not
                    np.can_cast(fill_value_dtype, values.dtype,
                                casting='same_kind')):
                raise ValueError("fill_value must be either 'None' or "
                                 "of a type compatible with values")

        for i, p in enumerate(points):
            if not np.all(np.diff(p) > 0.):
                raise ValueError("The points in dimension %d must be strictly "
                                 "ascending" % i)
            if not np.asarray(p).ndim == 1:
                raise ValueError("The points in dimension %d must be "
                                 "1-dimensional" % i)
            if not values.shape[i] == len(p):
                raise ValueError("There are %d points and %d values in "
                                 "dimension %d" % (len(p), values.shape[i], i))
        self.grid = tuple([np.asarray(p) for p in points])
        self.values = values

    def __call__(self, xi):
        """
        Interpolation at coordinates
        Parameters
        ----------
        xi : ndarray of shape (..., ndim)
            The coordinates to sample the gridded data at
        """
        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError("The requested sample points xi have dimension "
                             "%d, but this RegularGridInterpolator has "
                             "dimension %d" % (xi.shape[1], ndim))

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        indices, norm_distances, out_of_bounds = _find_indices(xi.T, self.grid)
        edges = np.array(list(itertools.product(*[[i, i + 1] for i in indices])))[:, :, 0]
        result = _evaluate_linear(self.values, indices, norm_distances, edges)
        if self.fill_value is not None:
            result[out_of_bounds] = self.fill_value

        return result.reshape(xi_shape[:-1] + self.values.shape[ndim:])


def _find_indices(xi, grid):
    # find relevant edges between which xi are situated
    indices = []
    # compute distance to lower edge in unity units
    norm_distances = []
    # check for out of bounds xi
    out_of_bounds = np.zeros((xi.shape[1]), dtype=np.bool)
    # iterate through dimensions
    for x, grid in zip(xi, grid):
        i = np.searchsorted(grid, x) - 1
        i[i < 0] = 0
        i[i > grid.size - 2] = grid.size - 2
        indices.append(i)
        norm_distances.append((x - grid[i]) /
                              (grid[i + 1] - grid[i]))
        out_of_bounds += x < grid[0]
        out_of_bounds += x > grid[-1]

    return indices, norm_distances, out_of_bounds


def _evaluate_linear(values_in, indices, norm_distances, edges):
    values = 0
    for j in range(edges.shape[0]):
        edge_indices = edges[j, :]
        weight = 1
        for ei, i, yi in zip(edge_indices, indices, norm_distances):
            if ei == i:
                weight *= 1 - yi
            else:
                weight *= yi
        values += values_in[edge_indices[0], edge_indices[1], edge_indices[2], :] * weight
    return np.atleast_2d(values)


def _ndim_coords_from_arrays(points, ndim=None):
    """
    Convert a tuple of coordinate arrays to a (..., ndim)-shaped array.
    """
    if isinstance(points, tuple) and len(points) == 1:
        # handle argument tuple
        points = points[0]
    if isinstance(points, tuple):
        p = np.broadcast_arrays(*points)
        n = len(p)
        for j in range(1, n):
            if p[j].shape != p[0].shape:
                raise ValueError("coordinate arrays do not have the same shape")
        points = np.empty(p[0].shape + (len(points),), dtype=float)
        for j, item in enumerate(p):
            points[...,j] = item
    else:
        points = np.asanyarray(points)
        if points.ndim == 1:
            if ndim is None:
                points = points.reshape(-1, 1)
            else:
                points = points.reshape(-1, ndim)
    return points
