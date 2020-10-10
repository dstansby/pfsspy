import functools
import warnings

from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import sunpy.map

import pfsspy.coords
from pfsspy.grid import Grid

# Default colourmap for magnetic field maps
_MAG_CMAP = 'RdBu'


def load_output(file):
    """
    Load a saved output file.

    Loads a file saved using :meth:`Output.save`.

    Parameters
    ----------
    file : str, file, :class:`~pathlib.Path`
        File to load.

    Returns
    -------
    :class:`Output`
    """
    with np.load(file) as f:
        alr = f['alr']
        als = f['als']
        alp = f['alp']
        rss = f['rss']
    ns = als.shape[1]
    nphi = alp.shape[0]
    nr = alr.shape[2]
    rss = rss[0]
    grid = Grid(ns, nphi, nr, rss)
    return Output(alr, als, alp, grid)


class Output:
    '''
    Output of PFSS modelling.

    Parameters
    ----------
    alr :
        Vector potential * grid spacing in radial direction.
    als :
        Vector potential * grid spacing in elevation direction.
    alp :
        Vector potential * grid spacing in azimuth direction.
    grid : Grid
        Grid that the output was caclulated on.
    input_map : sunpy.map.GenericMap
        The input map.

    Notes
    -----
    Instances of this class are intended to be created by `pfsspy.pfss`, and
    not by users.
    '''
    def __init__(self, alr, als, alp, grid, input_map=None):
        self._alr = alr
        self._als = als
        self._alp = alp
        self.grid = grid
        self.input_map = input_map

        # Cache attributes
        self._common_b_cache = None
        self._rgi = None

    def save(self, file):
        """
        Save the output to file.

        This saves the required information to reconstruct an Output object
        in a compressed binary numpy file (see :func:`numpy.savez_compressed`
        for more information). The file extension is ``.npz``, and is
        automatically added if not present.

        Parameters
        ----------
        file : str, file, :class:`~pathlib.Path`
            File to save to. If ``.npz`` extension isn't present it is added
            when saving the file.
        """
        np.savez_compressed(
            file, alr=self._alr, als=self._als, alp=self._alp,
            rss=np.array([self.grid.rss]))

    def _wcs_header(self):
        """
        Construct a world coordinate system describing the pfsspy solution.
        """
        return self.input_map.wcs

    @property
    def _lon0(self):
        """Longitude offset of the map."""
        return (self.input_map.meta['crval1'] *
                u.Unit(self.input_map.meta['cunit1']))

    @property
    def coordinate_frame(self):
        """
        The coordinate frame that the PFSS solution is in.

        Notes
        -----
        This is either a `~sunpy.coordinates.frames.HeliographicCarrington` or
        `~sunpy.coordinates.frames.HeliographicStonyhurst` frame, depending on
        the input map.
        """
        return self.input_map.coordinate_frame

    @property
    def dtime(self):
        return self.input_map.date

    @property
    def bunit(self):
        """
        `~astropy.units.Unit` of the input map data.
        """
        # Note that this can be removed once pfsspy depends on sunpy>=2.1, see
        # https://github.com/sunpy/sunpy/pull/4451
        unit_str = self.input_map.meta.get('bunit', None)
        if unit_str is None:
            return

        unit = u.Unit(unit_str, format='fits', parse_strict='silent')
        if isinstance(unit, u.UnrecognizedUnit):
            warnings.warn(f'Could not parse unit string "{unit_str}" as a valid FITS unit.\n'
                          'See https://fits.gsfc.nasa.gov/fits_standard.html '
                          'for the FITS unit standards.')
            unit = None
        return unit

    @property
    def source_surface_br(self):
        """
        Br on the source surface.

        Returns
        -------
        :class:`sunpy.map.GenericMap`
        """
        # Get radial component at the top
        br = self.bc[0][:, :, -1]
        # Remove extra ghost cells off the edge of the grid
        m = sunpy.map.Map((br.T, self._wcs_header()))
        vlim = np.max(np.abs(br))
        m.plot_settings['cmap'] = _MAG_CMAP
        m.plot_settings['vmin'] = -vlim
        m.plot_settings['vmax'] = vlim
        return m

    @property
    def source_surface_pils(self):
        """
        Coordinates of the polarity inversion lines on the source surface.

        Notes
        -----
        This is always returned as a list of coordinates, as in general there
        may be more than one polarity inversion line.
        """
        from skimage import measure
        m = self.source_surface_br
        contours = measure.find_contours(m.data, 0)
        contours = [m.wcs.pixel_to_world(c[:, 1], c[:, 0]) for c in contours]
        return contours

    @property
    def _brgi(self):
        """
        Regular grid interpolator for B.
        """
        from pfsspy.interpolator import RegularGridInterpolator as rgi
        if self._rgi is not None:
            return self._rgi

        f32 = np.float32
        # - (rho,s,phi) coordinates:
        rho = self.grid.rg.astype(f32)
        s = self.grid.sg.astype(f32)
        phi = self.grid.pg.astype(f32)
        br, bth, bph = self.bg[..., 2], self.bg[..., 1], self.bg[..., 0]

        # Because we need the cartesian grid to stretch just beyond r=rss,
        # add an extra dummy layer of magnetic field pointing radially outwards
        rho = np.append(rho, rho[-1] + 0.01)
        extras = np.ones(br.shape[0:2] + (1, ))
        br = np.concatenate((br, extras), axis=2).astype(f32)
        bth = np.concatenate((bth, 0 * extras), axis=2).astype(f32)
        bph = np.concatenate((bph, 0 * extras), axis=2).astype(f32)

        # - convert to Cartesian components and make interpolator on
        # (rho,s,phi) grid:
        ph3, s3, rh3 = np.meshgrid(phi, s, rho, indexing='ij')
        sin_th = np.sqrt(1 - s3**2)
        cos_th = s3
        sin_ph = np.sin(ph3)
        cos_ph = np.cos(ph3)

        # Directly stack the expressions below, to save a bit of memory
        # bx = (sin_th * cos_ph * br) + (cos_th * cos_ph * bth) - (sin_ph * bph)
        # by = (sin_th * sin_ph * br) + (cos_th * sin_ph * bth) + (cos_ph * bph)
        # bz = (cos_th * br) - (sin_th * bth)
        bstack = np.stack(((sin_th * cos_ph * br) + (cos_th * cos_ph * bth) - (sin_ph * bph),
                           (sin_th * sin_ph * br) + (cos_th * sin_ph * bth) + (cos_ph * bph),
                           (cos_th * br) - (sin_th * bth)),
                          axis=-1)

        self._rgi = rgi((phi, s, rho), bstack)
        return self._rgi

    def _bTrace(self, t, coord, direction):
        """
        Return B/|B| for use by the field line tracer.
        """
        x, y, z = coord
        # (ph, s, rh) coordinates of current point:
        rho, s, phi = pfsspy.coords.cart2strum(x, y, z)

        # Check if position vector is outside the data limits
        if rho < 0 or rho > np.log(self.grid.rss):
            return np.array([0, 0, 0])
            # raise _OutOfBoundsError

        b1 = self._brgi(np.stack((phi, s, rho)))[0]
        return direction * b1 / np.linalg.norm(b1)

    def trace(self, tracer, seeds):
        """
        Parameters
        ----------
        tracer : tracing.Tracer
            Field line tracer.
        seeds : astropy.coordinates.SkyCoord
            Starting coordinates.
        """
        return tracer.trace(seeds, self)

    def _integrate_one_way(self, dt, start_point, rtol, atol):
        import scipy.integrate

        direction = np.sign(dt)

        def finish_integration(t, coord):
            r = np.linalg.norm(coord)
            ret = (r - 1) * (r - self.grid.rss)
            return ret

        finish_integration.terminal = True
        # The integration domain is deliberately huge, because the
        # the interation automatically stops when an out of bounds error
        # is thrown
        t_span = (0, 1e4)

        def fun(t, y):
            return self._bTrace(t, y, direction)

        res = scipy.integrate.solve_ivp(
            fun, t_span, start_point, method='RK23',
            rtol=rtol, atol=atol, events=finish_integration)

        xout = res.y
        return xout

    @property
    def _al(self):
        """
        Vector potential times cell edge lenghts.

        Returns ar*Lr, as*Ls, ap*Lp on cell edges.
        """
        return self._alr, self._als, self._alp

    @property
    def bc(self):
        """
        B on the centres of the cell faces.

        Returns
        -------
        br
        btheta
        bphi
        """
        br, bs, bp, Sbr, Sbs, Sbp = self._common_b()
        # Remove area factors:
        br = br.copy()
        bs = bs.copy()
        bp = bp.copy()
        for i in range(self.grid.nphi + 2):
            br[i, :, :] = br[i, :, :] / Sbr
            bs[i, :, :] = bs[i, :, :] / Sbs
        for i in range(self.grid.nphi + 1):
            bp[i, :, :] = bp[i, :, :] / Sbp

        # Slice to remove ghost cells
        return br[1:-1, 1:-1, :], -bs[1:-1, :, 1:-1], bp[:, 1:-1, 1:-1]

    @property
    @functools.lru_cache(maxsize=1)
    def bg(self):
        """
        B as a (weighted) averaged on grid points.

        Returns
        -------
        numpy.ndarray
            A ``(nphi + 1, ns + 1, nrho + 1, 3)`` shaped array.
            The last index gives the corodinate axis, 0 for Bphi, 1 for Bs, 2
            for Brho.
        """
        br, bs, bp, Sbr, Sbs, Sbp = self._common_b()
        # Weighted average to grid points:
        brg = br[:-1, :-1, :] + br[1:, :-1, :] + br[1:, 1:, :] + br[:-1, 1:, :]
        bsg = bs[:-1, :, :-1] + bs[1:, :, :-1] + bs[1:, :, 1:] + bs[:-1, :, 1:]
        bpg = bp[:, :-1, :-1] + bp[:, 1:, :-1] + bp[:, 1:, 1:] + bp[:, :-1, 1:]
        for i in range(self.grid.nphi + 1):
            brg[i, :, :] /= 2 * (Sbr[:-1, :] + Sbr[1:, :])
            bsg[i, :, :] /= 2 * (Sbs[:, :-1] + Sbs[:, 1:])
        for i in range(self.grid.nphi + 1):
            bpg[i, :, :] /= (Sbp[:-1, :-1] + Sbp[1:, :-1] +
                             Sbp[1:, 1:] + Sbp[:-1, 1:])
        bsg *= -1
        out = np.stack((bpg, bsg, brg), axis=-1)
        out.flags.writeable = False
        return out

    @property
    @functools.lru_cache(maxsize=1)
    def _modbg(self):
        return np.linalg.norm(self.bg, axis=-1)

    def _common_b(self):
        """
        Common code needed to calculate magnetic field from vector potential.
        """
        if self._common_b_cache is not None:
            return self._common_b_cache

        dr = self.grid.dr
        ds = self.grid.ds
        dp = self.grid.dp

        nr = self.grid.nr
        ns = self.grid.ns
        nphi = self.grid.nphi

        rss = self.grid.rss

        rc = self.grid.rc
        sc = self.grid.sc

        rg = self.grid.rg
        sg = self.grid.sg

        alr, als, alp = self._al

        # Centre of cells in rho (including ghost cells)
        rc = np.linspace(-0.5 * dr, np.log(rss) + 0.5 * dr, nr + 2)
        rrc = np.exp(rc)
        thc = np.zeros(ns + 2) - 1
        thc[1:-1] = np.arccos(sc)
        # Centre of cells in phi (including ghost cells)
        pc = np.linspace(-0.5 * dp, 2 * np.pi + 0.5 * dp, nphi + 2)

        # Required face normals:
        dnp = np.zeros((ns + 2, 2))
        dns = np.zeros((ns + 1, 2))
        dnr = np.zeros(ns + 2)
        for k in range(2):
            for j in range(1, ns + 1):
                dnp[j, k] = rrc[k] * np.sqrt(1 - sc[j - 1]**2) * dp
            dnp[0, k] = dnp[1, k]
            dnp[-1, k] = dnp[-2, k]
            for j in range(1, ns):
                dns[j, k] = rrc[k] * (np.arcsin(sc[j]) - np.arcsin(sc[j - 1]))
            dns[0, k] = dns[1, k]
            dns[-1, k] = dns[-2, k]
        for j in range(ns + 2):
            dnr[j] = rrc[0] * (np.exp(dr) - 1)
        dnr[0] = -dnr[0]
        dnr[-1] = -dnr[-1]

        # Required area factors:
        Sbr = np.zeros((ns + 2, nr + 1))
        for k in range(nr + 1):
            Sbr[1:-1, k] = np.exp(2 * rg[k]) * ds * dp
            Sbr[0, k] = Sbr[1, k]
            Sbr[-1, k] = Sbr[-2, k]
        Sbs = np.zeros((ns + 1, nr + 2))
        for k in range(nr + 2):
            for j in range(1, ns):
                Sbs[j, k] = 0.5 * np.exp(2 * rc[k] - dr) * dp * (np.exp(2 * dr) - 1) * np.sqrt(1 - sg[j]**2)
            Sbs[0, k] = Sbs[1, k]
            Sbs[-1, k] = Sbs[-2, k]
        Sbp = np.zeros((ns + 2, nr + 2))
        for k in range(nr + 2):
            for j in range(1, ns + 1):
                Sbp[j, k] = 0.5 * np.exp(2 * rc[k] - dr) * (np.exp(2 * dr) - 1) * (np.arcsin(sg[j]) - np.arcsin(sg[j - 1]))
            Sbp[0, k] = Sbp[1, k]
            Sbp[-1, k] = Sbp[-2, k]

        # Compute br*Sbr, bs*Sbs, bp*Sbp at cell centres by Stokes theorem:
        br = np.zeros((nphi + 2, ns + 2, nr + 1))
        bs = np.zeros((nphi + 2, ns + 1, nr + 2))
        bp = np.zeros((nphi + 1, ns + 2, nr + 2))
        br[1:-1, 1:-1, :] = als[1:, :, :] - als[:-1, :, :] + alp[:, :-1, :] - alp[:, 1:, :]
        bs[1:-1, :, 1:-1] = alp[:, :, 1:] - alp[:, :, :-1]
        bp[:, 1:-1, 1:-1] = als[:, :, :-1] - als[:, :, 1:]

        # Fill ghost values with boundary conditions:
        # - zero-gradient at outer boundary:
        bs[1:-1, :, -1] = 2 * bs[1:-1, :, -2] - bs[1:-1, :, -3]
        bp[:, 1:-1, -1] = 2 * bp[:, 1:-1, -2] - bp[:, 1:-1, -3]
        # - periodic in phi:
        bs[0, :, :] = bs[-2, :, :]
        bs[-1, :, :] = bs[1, :, :]
        br[0, :, :] = br[-2, :, :]
        br[-1, :, :] = br[1, :, :]
        # js = jp = 0 at photosphere:
        for i in range(nphi + 1):
            bp[i,:,0] = Sbp[:,0]/dnp[:,0]*(bp[i,:,1]*dnp[:,1]/Sbp[:,1] + br[i,:,0]*dnr[:]/Sbr[:,0] - br[i+1,:,0]*dnr[:]/Sbr[:,0])
        for i in range(nphi + 2):
            bs[i,:,0] = Sbs[:,0]/dns[:,0]*(bs[i,:,1]*dns[:,1]/Sbs[:,1] + br[i,:-1,0]*dnr[:-1]/Sbr[:-1,0] - br[i,1:,0]*dnr[1:]/Sbr[1:,0])
        # - polar boundaries as in dumfric:
        for i in range(nphi + 2):
            i1 = (i + nphi//2) % nphi
            br[i,-1,:] = br[i1,-2,:]
            br[i,0,:] = br[i1,1,:]
            bs[i, -1, :] = 0.5 * (bs[i, -2, :] - bs[i1, -2, :])
            bs[i, 0, :] = 0.5 * (bs[i, 1, :] - bs[i1, 1, :])
        for i in range(nphi + 1):
            i1 = (i + nphi // 2) % nphi
            bp[i, -1, :] = -bp[i1, -2, :]
            bp[i, 0, :] = -bp[i1, 1, :]

        self._common_b_cache = br, bs, bp, Sbr, Sbs, Sbp
        return self._common_b_cache

    def get_bvec(self, coords, out_type="spherical"):
        """
        Evaluate magnetic vectors in pfss model.

        Method which takes an arbitrary astropy SkyCoord and
        returns a numpy array containing magnetic field vectors
        evaluated from the parent pfsspy.Output pfss model at
        the locations specified by the SkyCoords

        Parameters
        ----------
        coords : `astropy.SkyCoord`
            An arbitary point or set of points (length N >= 1)
            in the PFSS model domain (1Rs < r < Rss)

        out_type : str, optional
            Takes values 'spherical' (default) or 'cartesian'
            and specifies whether the output vector is in
            spherical coordinates (B_r,B_theta,B_phi) or
            cartesian (B_x,B_y,B_z)

        Returns
        -------
        bvec : ndarray
            Magnetic field vectors at the requested locations
            ndarray.shape = (N,3), units nT)

        Notes
        -----
        The output coordinate system is defined by the input
        magnetogram with x-z plane equivalent to the plane
        containing the Carrington meridian (0 deg longitude)

        The spherical coordinates follow the physics convention:
        https://upload.wikimedia.org/wikipedia/commons/thumb/4/4f/3D_Spherical.svg/240px-3D_Spherical.svg.png)
        Therefore the polar angle (theta) is the co-latitude, rather
        than the latitude, with range 0 (north pole) to 180 degrees
        (south pole)

        The conversion which relates the spherical and cartesian
        coordinates is as follows:

        .. math:: B_R = sin\\theta cos\\phi B_x + sin\\theta sin\\phi B_y + cos\\theta B_z
        .. math:: B_\\theta = cos\\theta cos\\phi B_x + cos\\theta sin\\phi B_y - sin\\theta B_z
        .. math:: B_\\phi = -sin\\phi B_x + cos\\phi B_y

        The above equations may be written as a (3x3) matrix and
        inverted to retrieve the inverse transformation (cartesian from spherical)
        """

        # Assert skycoord is type astropy.coordinates.SkyCoord
        if not isinstance(coords, SkyCoord):
            raise ValueError("coords must be of type "
                             "astropy.coordinates.SkyCoord")

        # Ensure representation type is spherical for input to interpolator
        coords.representation_type = "spherical"

        # Check coord_type is cartesian or spherical
        if out_type not in ["cartesian", "spherical"]:
            raise ValueError("out_type must be 'cartesian' or 'spherical' "
                             f"(got {out_type})")

        # Raise warning if input skycoord obstime does not match
        # self.coordinate_frame.obstime
        if np.any(self.coordinate_frame.obstime.to_datetime() !=
                  coords.obstime.to_datetime()):
            warnings.warn("The obstime of one of more input coordinates "
                          "do not match the pfss model obstime.")

        # Convert SkyCoord to pfsspy.Output coordinate frame
        coords.transform_to(self.coordinate_frame)

        # Do interpolation (returns cartesian vector)
        bvecs = self._brgi(np.array([coords.lon.to("rad").value,
                                     np.sin(coords.lat).value,
                                     np.log(coords.radius.to("R_sun").value)]).T
                           )

        # Convert to spherical if requested
        if out_type == "spherical":
            # Generate vector of 3x3 rotation matrices
            M = np.array([
                [np.cos(coords.lat).value * np.cos(coords.lon).value,
                 np.cos(coords.lat).value * np.sin(coords.lon).value,
                 np.sin(coords.lat).value],
                [np.sin(coords.lat).value * np.cos(coords.lon).value,
                 np.sin(coords.lat).value * np.sin(coords.lon).value,
                 -np.cos(coords.lat).value],
                [-np.sin(coords.lon).value,
                 np.cos(coords.lon).value,
                 np.zeros(len(coords))],
            ])
            bvecs = np.array([np.dot(M_.T,v) for M_,v in zip(M.T,bvecs)])
        if self.bunit is not None:
            bvecs *= self.bunit
        return bvecs
