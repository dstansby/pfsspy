import warnings

import numpy as np
import sunpy.map

from pfsspy.grid import Grid
import pfsspy.utils


class Input:
    r"""
    Input to PFSS modelling.

    .. warning::
        The input must be on a regularly spaced grid in :math:`\phi` and
        :math:`s = \cos (\theta)`. See `pfsspy.grid` for more
        information on the coordinate system.

    Parameters
    ----------
    br : sunpy.map.GenericMap
        Boundary condition of radial magnetic field at the inner surface.
        Note that the data *must* have a cylindrical equal area projection.

    nr : int
        Number of cells in the radial direction to calculate the PFSS solution
        on.

    rss : float
        Radius of the source surface, as a fraction of the solar radius.
    """
    def __init__(self, br, nr, rss):
        if not isinstance(br, sunpy.map.GenericMap):
            raise ValueError('br must be a SunPy Map')
        if np.any(~np.isfinite(br.data)):
            raise ValueError('At least one value in the input is NaN or '
                             'non-finite. The input must consist solely of '
                             'finite values.')
        if np.mean(br.data) > 1e-10:
            warnings.warn('Input data has a non-zero mean. '
                          'pfsspy will ignore this non-zero monopole term '
                          'when calculating the PFSS solution.')

        pfsspy.utils.is_cea_map(br, error=True)
        pfsspy.utils.is_full_sun_synoptic_map(br, error=True)

        self._map_in = br
        self.dtime = br.date
        self.br = br.data

        # Force some nice defaults
        self._map_in.plot_settings['cmap'] = 'RdBu'
        lim = np.nanmax(np.abs(self._map_in.data))
        self._map_in.plot_settings['vmin'] = -lim
        self._map_in.plot_settings['vmax'] = lim

        ns = self.br.shape[0]
        nphi = self.br.shape[1]
        self.grid = Grid(ns, nphi, nr, rss)

    @property
    def map(self):
        """
        :class:`sunpy.map.GenericMap` representation of the input.
        """
        return self._map_in
