import astropy.coordinates as coord
import astropy.constants as const
import sunpy.coordinates.frames as frames

import numpy as np
import scipy.interpolate


class FieldLine:
    """
    A single magnetic field line.

    Parameters
    ----------
    x, y, z : array
        Field line coordinates in a Carrington frame of reference. Must be in
        units of solar radii.
    dtime : astropy.time.Time
        Time at which the field line was traced. Needed for transforming the
        field line coordinates to other coordinate frames.
    output : Output
        The PFSS output through which this field line was traced.

    Attributes
    ----------
    coords : astropy.coordinates.SkyCoord
        Field line coordinates.
    """
    def __init__(self, x, y, z, dtime, output):
        self.coords = coord.SkyCoord(x=x * const.R_sun,
                                     y=y * const.R_sun,
                                     z=z * const.R_sun,
                                     frame=frames.HeliographicCarrington,
                                     obstime=dtime,
                                     representation_type='cartesian')
        self._is_open = None
        self._polarity = None
        self._expansion_factor = None
        self._output = output

    @property
    def is_open(self):
        """
        Returns ``True`` if one of the field line is connected to the solar
        surface and one to the outer boundary, ``False`` otherwise.
        """
        if self._is_open is None:
            r = coord.SkyCoord(
                self.coords, representation_type='spherical').radius
            rtol = 0.1
            if np.abs(r[0] - r[-1]) < r[0] * rtol:
                self._is_open = False
            else:
                self._is_open = True
        return self._is_open

    @property
    def polarity(self):
        """
        Magnetic field line polarity.

        Returns
        -------
        pol : int
            0 if the field line is closed, otherwise sign(Br) of the magnetic
            field on the solar surface.
        """
        if self._polarity is None:
            if not self.is_open:
                self._polarity = 0
            # Because the field lines are integrated forwards, if the end point
            # is on the outer boundary the field is outwards
            foot1 = coord.SkyCoord(
                self.coords[0], representation_type='spherical')
            foot2 = coord.SkyCoord(
                self.coords[-1], representation_type='spherical')

            if foot2.radius - foot1.radius > 0:
                self._polarity = 1
            else:
                self._polarity = -1
        return self._polarity

    @property
    def solar_footpoint(self):
        """
        Solar surface magnetic field footpoint.

        This is the ends of the magnetic field line that lies on the solar
        surface.

        Returns
        -------
        footpoint : :class:`~astropy.coordinates.SkyCoord`

        Notes
        -----
        For a closed field line, both ends lie on the solar surface. This
        method returns the field line pointing out from the solar surface in
        this case.
        """
        if self.polarity == 1 or not self.is_open:
            return coord.SkyCoord(self.coords[0])
        else:
            return coord.SkyCoord(self.coords[-1])

    @property
    def source_surface_footpoint(self):
        """
        Solar surface magnetic field footpoint.

        This is the ends of the magnetic field line that lies on the solar
        surface.

        Returns
        -------
        footpoint : :class:`~astropy.coordinates.SkyCoord`

        Notes
        -----
        For a closed field line, both ends lie on the solar surface. This
        method returns the field line pointing out from the solar surface in
        this case.
        """
        if self.polarity == 1 or not self.is_open:
            return coord.SkyCoord(self.coords[-1])
        else:
            return coord.SkyCoord(self.coords[0])

    @property
    def expansion_factor(self):
        r"""
        Magnetic field expansion factor.

        The expansion factor is defnied as
        :math:`(r_{\odot}^{2} B_{\odot}) / (r_{ss}^{2} B_{ss}))`

        Returns
        -------
        exp_fact : float
            Field line expansion factor.
            If field line is closed, returns np.nan.
        """
        if self._expansion_factor is not None:
            return self._expansion_factor
        import scipy.interpolate

        if not self.is_open:
            return np.nan
        # Extract ends of magnetic field line, and get them in spherical coords
        foot1 = coord.SkyCoord(self.coords[0], representation_type='spherical')
        foot2 = coord.SkyCoord(self.coords[-1], representation_type='spherical')
        if foot1.radius > foot2.radius:
            solar_foot = foot2
            source_foot = foot1
        else:
            solar_foot = foot1
            source_foot = foot2

        def interp(map, coord):
            phi = coord.lon
            s = np.sin(coord.lat)
            interpolator = scipy.interpolate.RectBivariateSpline(
                self._output.grid.pg, self._output.grid.sg, map)
            return interpolator(phi, s)

        # Get output magnetic field, and calculate |B|
        br, bs, bphi = self._output.bg
        modb = np.sqrt(br**2 + bs**2 + bphi**2)
        # Interpolate at each end of field line
        b_solar = interp(modb[:, :, 0], solar_foot)[0, 0]
        b_source = interp(modb[:, :, -1], source_foot)[0, 0]
        self._expansion_factor = ((1**2 * b_solar) /
                                  (self._output.grid.rss**2 * b_source))
        return self._expansion_factor
