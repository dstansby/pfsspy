import astropy.coordinates as coord
import astropy.constants as const
import astropy.units as u

import numpy as np

from pfsspy import coords

import functools


class FieldLines:
    """
    A collection of :class:`FieldLine`.

    Parameters
    ----------
    field_lines : list of `FieldLine`.
    """
    def __init__(self, field_lines):
        self.field_lines = np.array(field_lines)

    def __getitem__(self, idx):
        return self.field_lines[idx]

    def __len__(self):
        return len(self.field_lines)

    @property
    @functools.lru_cache()
    def polarities(self):
        """
        Magnetic field line polarities. ``0`` for closed, otherwise sign(Br) on
        the solar surface.
        """
        polarities = [fline.polarity for fline in self.field_lines]
        return np.array(polarities, dtype=int)

    @property
    def connectivities(self):
        """
        Field line connectivities. ``1`` for open, ``0`` for closed.
        """
        return np.abs(self.polarities)

    @property
    def expansion_factors(self):
        """
        Expansion factors. Set to NaN for closed field lines.
        """
        return np.array([fline.expansion_factor for fline in self.field_lines])

    @property
    def open_field_lines(self):
        """
        An `OpenFieldLines` object containing open field lines.
        """
        open_idxs = np.where(self.connectivities == 1)[0]
        return OpenFieldLines(np.array(self.field_lines)[open_idxs])

    @property
    def closed_field_lines(self):
        """
        An `ClosedFieldLines` object containing open field lines.
        """
        closed_idxs = np.where(self.connectivities == 0)[0]
        return ClosedFieldLines(self.field_lines[closed_idxs])


class OpenFieldLines(FieldLines):
    """
    A set of open field lines.
    """
    def __init__(self, field_lines):
        super().__init__(field_lines)
        if not np.all(self.connectivities):
            raise ValueError('Not all field lines are open')

    @property
    @functools.lru_cache()
    def source_surface_feet(self):
        """
        Coordinates of the source suface footpoints.
        """
        x = np.array([fline._x[fline._ss_coord_index] for fline in self.field_lines])
        y = np.array([fline._y[fline._ss_coord_index] for fline in self.field_lines])
        z = np.array([fline._z[fline._ss_coord_index] for fline in self.field_lines])

        return FieldLine._coords(x, y, z, self.field_lines[0]._output)

    @property
    @functools.lru_cache()
    def solar_feet(self):
        """
        Coordinates of the solar footpoints.
        """
        x = np.array([fline._x[fline._solar_coord_index] for fline in self.field_lines])
        y = np.array([fline._y[fline._solar_coord_index] for fline in self.field_lines])
        z = np.array([fline._z[fline._solar_coord_index] for fline in self.field_lines])

        return FieldLine._coords(x, y, z, self.field_lines[0]._output)


class ClosedFieldLines(FieldLines):
    """
    A set of closed field lines.
    """
    def __init__(self, field_lines):
        super().__init__(field_lines)
        if np.any(self.connectivities):
            raise ValueError('Not all field lines are closed')


class FieldLine:
    """
    A single magnetic field line.

    Parameters
    ----------
    x, y, z :
        Field line coordinates in cartesian coordinates.
    output : Output
        The PFSS output through which this field line was traced.
    """
    def __init__(self, x, y, z, output):
        self._x = np.array(x)
        self._y = np.array(y)
        self._z = np.array(z)
        self._r = np.sqrt(self._x**2 + self._y**2 + self._z**2)
        self._output = output
        # Set _is_open
        atol = 0.1
        self._is_open = np.abs(self._r[0] - self._r[-1]) > atol
        # Set _polarity
        self._polarity = -np.sign(self._r[0] - self._r[-1]) * self._is_open

    @property
    def coords(self):
        """
        Field line `~astropy.coordinates.SkyCoord`.
        """
        return self._coords(self._x, self._y, self._z, self._output)

    @staticmethod
    def _coords(x, y, z, output):
        r, lat, lon = coord.cartesian_to_spherical(x, y, z)
        r *= const.R_sun
        lon += output._lon0 + 180 * u.deg
        coords = coord.SkyCoord(
            lon, lat, r, frame=output.coordinate_frame)
        return coords

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
        return self.coords[self._solar_coord_index]

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
        return self.coords[self._ss_coord_index]

    @property
    def _ss_coord_index(self):
        """
        Return 0 or -1 depending on which end of the coordinate array is the
        source surface footpoint.
        """
        if self.polarity == 1 or not self.is_open:
            return -1
        else:
            return 0

    @property
    def _solar_coord_index(self):
        return -1 - self._ss_coord_index

    @property
    def b_along_fline(self):
        """
        The magnetic field vectors along the field line.
        """
        coords = self.coords
        return self._output.get_bvec(coords)

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
        """
        import scipy.interpolate

        if not self.is_open:
            return np.nan

        solar_foot = self._solar_coord_index
        source_foot = self._ss_coord_index

        def interp(map, idx):
            x, y, z = self._x[idx], self._y[idx], self._z[idx]
            rho, s, phi = coords.cart2strum(x, y, z)
            interpolator = scipy.interpolate.RegularGridInterpolator(
                (self._output.grid.pg, self._output.grid.sg), map, method='linear')
            return interpolator((phi, s))

        # Get output magnetic field, and calculate |B|
        modb = self._output._modbg
        # Interpolate at each end of field line
        b_solar = interp(modb[:, :, 0], solar_foot)
        b_source = interp(modb[:, :, -1], source_foot)
        expansion_factor = ((1**2 * b_solar) /
                            (self._output.grid.rss**2 * b_source))
        return expansion_factor
