import astropy.coordinates as coord
import astropy.constants as const
import sunpy.coordinates.frames as frames

import numpy as np
import scipy.interpolate

import functools


class FieldLines:
    """
    A collection of :class:`FieldLine`.

    Parameters
    ----------
    field_lines : list of FieldLine
    """
    def __init__(self, field_lines):
        self.field_lines = field_lines

    def __getitem__(self, idx):
        return self.field_lines[idx]

    @property
    @functools.lru_cache()
    def solar_feet(self):
        """
        Coordinates of the solar footpoints.

        Notes
        -----
        For closed field lines, the footpoint pointing out from the solar
        surface is returned.
        """
        solar_feet = [fline.solar_footpoint for fline in self.field_lines]

        if len(solar_feet) == 1:
            solar_feet = solar_feet[0]
        else:
            solar_feet = coord.concatenate(solar_feet)

        return solar_feet

    @property
    @functools.lru_cache()
    def polarities(self):
        """
        Magnetic field line polarities.
        """
        polarities = [fline.polarity for fline in self.field_lines]
        return np.array(polarities)

    @property
    @functools.lru_cache()
    def source_surface_feet(self):
        """
        Coordinates of the source suface footpoints.

        Notes
        -----
        For closed field lines, there is no source surface footpoint, but
        instead the solar footpoint pointing in towards the solar surface is
        returned.
        """
        source_surface_feet = [fline.source_surface_footpoint for
                               fline in self.field_lines]

        if len(source_surface_feet) == 1:
            source_surface_feet = source_surface_feet[0]
        else:
            source_surface_feet = coord.concatenate(source_surface_feet)

        return source_surface_feet


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
        self._output = output
        # Set _is_open
        r = np.sqrt(np.array(x)**2 + np.array(y)**2 + np.array(z)**2)
        rtol = 0.1
        self._is_open = np.abs(r[0] - r[-1]) > 1 * rtol
        # Set _polarity
        self._polarity = -np.sign(r[0] - r[-1]) * self._is_open

    @property
    def is_open(self):
        """
        Returns ``True`` if one of the field line is connected to the solar
        surface and one to the outer boundary, ``False`` otherwise.
        """
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
    @functools.lru_cache()
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
        expansion_factor = ((1**2 * b_solar) /
                            (self._output.grid.rss**2 * b_source))
        return expansion_factor
