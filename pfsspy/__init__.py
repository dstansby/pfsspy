import distutils.version
import functools

import astropy
import astropy.coordinates as coord
import astropy.constants as const
import astropy.time
import astropy.units as u
import astropy.wcs

from sunpy.coordinates import frames
import sunpy.map

import numpy as np
import scipy.linalg as la
import pfsspy.coords
import pfsspy.tracing
import pfsspy.fieldline
# Import here to register map sources
import pfsspy.map
# Utility Functions
import pfsspy.utils

from .grid import Grid

HAS_NUMBA = False
try:
    import numba
    HAS_NUMBA = True
except Exception:
    pass

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# Do a version check for astropy
if (distutils.version.LooseVersion(astropy.__version__) <
        distutils.version.LooseVersion("3")):
    raise RuntimeError('pfsspy requires astropy v3 to run ' +
                       f'(found version {astropy.__version__} installed)')

# Default colourmap for magnetic field maps
_MAG_CMAP = 'RdBu'


class Input:
    r"""
    Input to PFSS modelling.

    .. warning::
        The input must be on a regularly spaced grid in :math:`\phi` and
        :math:`s = \cos (\theta)`. See :mod:`pfsspy.coords` for more
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

        self._validate_cea(br.meta)

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

    @staticmethod
    def _validate_cea(meta):
        for i in ('1', '2'):
            proj = meta[f'ctype{i}'][5:8]
            if proj != 'CEA':
                raise ValueError(f'Projection type in CTYPE{i} keyword '
                                 f'must be CEA (got "{proj}")')

    @property
    def map(self):
        """
        :class:`sunpy.map.GenericMap` representation of the input.
        """
        return self._map_in


def carr_cea_wcs_header(dtime, shape):
    """
    Create a Carrington WCS header for a Cylindrical Equal Area (CEA)
    projection. See [1]_ for information on how this is constructed.

    dtime : datetime, None
        Datetime to associate with the map.
    data : [ntheta, nphi]
        Map data. First entry is latitude, second entry is longitude.

    References
    ----------
    .. [1] W. T. Thompson, "Coordinate systems for solar image data",
       https://doi.org/10.1051/0004-6361:20054262
    """
    # If datetime is None, put in a dummy value here to make
    # make_fitswcs_header happy, then strip it out at the end
    obstime = dtime or astropy.time.Time('2000-1-1')

    frame_out = coord.SkyCoord(
        0 * u.deg, 0 * u.deg, obstime=obstime,
        frame="heliographic_carrington", observer='sun')
    # Construct header
    header = sunpy.map.make_fitswcs_header(
        shape, frame_out,
        scale=[180 / shape[0],
               360 / shape[1]] * u.deg / u.pix,
        reference_pixel=[(shape[1] / 2) + 0.5, (shape[0] / 2) + 0.5] * u.pix,
        projection_code="CEA")

    # Fill in these missing values
    header['PV1_1'] = 1
    header['PV2_1'] = 1
    # Fix CELT for lat axis
    header['CDELT2'] = (180 / np.pi) * (2 / shape[0])
    # pop out the time if it isn't supplied
    if dtime is None:
        header.pop('date-obs')
    return header


class _OutOfBoundsError(RuntimeError):
    pass


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
        br = br[1:-1, 1:-1].T
        m = sunpy.map.Map((br, self._wcs_header()))
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
        dt = np.abs(dt)
        t = 0.0
        xout = np.atleast_2d(start_point.copy())

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
    def al(self):
        """
        Vector potential times cell edge lenghts.

        Returns ar*Lr, as*Ls, ap*Lp on cell edges.
        """
        return self._alr, self._als, self._alp

    @property
    def bc(self):
        """
        B on the centres of the cell faces.
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

        return br, -bs, bp

    @property
    @functools.lru_cache(maxsize=1)
    def bg(self):
        """
        B as a (weighted) averaged on grid points.

        Returns
        -------
        numpy.ndarray
            A (nphi, ns, nrho, 3) shaped array. The last index gives the
            corodinate axis, 0 for Bphi, 1 for Bs, 2 for Brho.
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

        alr, als, alp = self.al

        rc = np.linspace(-0.5 * dr, np.log(rss) + 0.5 * dr, nr + 2)
        rrc = np.exp(rc)
        thc = np.zeros(ns + 2) - 1
        thc[1:-1] = np.arccos(sc)
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
                Sbs[j,k] = 0.5*np.exp(2*rc[k] - dr)*dp*(np.exp(2*dr)-1)*np.sqrt(1 - sg[j]**2)
            Sbs[0,k] = Sbs[1,k]
            Sbs[-1,k] = Sbs[-2,k]
        Sbp = np.zeros((ns+2,nr+2))
        for k in range(nr+2):
            for j in range(1,ns+1):
                Sbp[j,k] = 0.5*np.exp(2*rc[k] - dr)*(np.exp(2*dr) - 1)*(np.arcsin(sg[j]) - np.arcsin(sg[j-1]))
            Sbp[0,k] = Sbp[1,k]
            Sbp[-1,k] = Sbp[-2,k]

        # Compute br*Sbr, bs*Sbs, bp*Sbp at cell centres by Stokes theorem:
        br = np.zeros((nphi+2,ns+2,nr+1))
        bs = np.zeros((nphi+2,ns+1,nr+2))
        bp = np.zeros((nphi+1,ns+2,nr+2))
        br[1:-1,1:-1,:] = als[1:,:,:] - als[:-1,:,:] + alp[:,:-1,:] - alp[:,1:,:]
        bs[1:-1,:,1:-1] = alp[:,:,1:] - alp[:,:,:-1]
        bp[:,1:-1,1:-1] = als[:,:,:-1] - als[:,:,1:]

        # Fill ghost values with boundary conditions:
        # - zero-gradient at outer boundary:
        bs[1:-1,:,-1] = 2*bs[1:-1,:,-2] - bs[1:-1,:,-3]
        bp[:,1:-1,-1] = 2*bp[:,1:-1,-2] - bp[:,1:-1,-3]
        # - periodic in phi:
        bs[0,:,:] = bs[-2,:,:]
        bs[-1,:,:] = bs[1,:,:]
        br[0,:,:] = br[-2,:,:]
        br[-1,:,:] = br[1,:,:]
        # js = jp = 0 at photosphere:
        for i in range(nphi+1):
            bp[i,:,0] = Sbp[:,0]/dnp[:,0]*(bp[i,:,1]*dnp[:,1]/Sbp[:,1] + br[i,:,0]*dnr[:]/Sbr[:,0] - br[i+1,:,0]*dnr[:]/Sbr[:,0])
        for i in range(nphi+2):
            bs[i,:,0] = Sbs[:,0]/dns[:,0]*(bs[i,:,1]*dns[:,1]/Sbs[:,1] + br[i,:-1,0]*dnr[:-1]/Sbr[:-1,0] - br[i,1:,0]*dnr[1:]/Sbr[1:,0])
        # - polar boundaries as in dumfric:
        for i in range(nphi+2):
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


def _eigh(A):
    return np.linalg.eigh(A)


def _compute_r_term(m, k, ns, Q, brt, lam, ffm, nr, ffp, psi, psir):
    for l in range(ns):
        # - sum c_{lm} + d_{lm}
        cdlm = np.dot(Q[:, l], brt[:, m]) / lam[l]
        # - ratio c_{lm}/d_{lm} [numerically safer this way up]
        ratio = (ffm[l]**(nr - 1) - ffm[l]**nr) / (ffp[l]**nr - ffp[l]**(nr - 1))
        dlm = cdlm / (1.0 + ratio)
        clm = ratio * dlm
        psir[:, l] = clm * ffp[l]**k + dlm * ffm[l]**k

    # - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}**j
    psi[:, :, m] = np.dot(psir, Q.T)
    return psi, psir


def _als_alp(nr, nphi, Fs, psi, Fp, als, alp):
    for j in range(nr + 1):
        for i in range(nphi + 1):
            als[i, :, j] = Fs * (psi[j, :, ((i - 1) % nphi)] - psi[j, :, ((i) % nphi)])
        for i in range(nphi):
            alp[i, 1:-1, j] = Fp[1:-1] * (psi[j, 1:, i] - psi[j, :-1, i])
    return als, alp


def _A_diag(A, ns, Vg, Uc, mu, m):
    for j in range(ns):
        A[j, j] = Vg[j] + Vg[j + 1] + Uc[j] * mu[m]
    return A


if HAS_NUMBA:
    _eigh = numba.jit(nopython=True)(_eigh)
    _compute_r_term = numba.jit(nopython=True)(_compute_r_term)
    _als_alp = numba.jit(nopython=True)(_als_alp)
    _A_diag = numba.jit(nopython=True)(_A_diag)


def pfss(input):
    r"""
    Compute PFSS model.

    Extrapolates a 3D PFSS using an eigenfunction method in :math:`r,s,p`
    coordinates, on the dumfric grid
    (equally spaced in :math:`\rho = \ln(r/r_{sun})`,
    :math:`s= \cos(\theta)`, and :math:`p=\phi`).

    The output should have zero current to machine precision,
    when computed with the DuMFriC staggered discretization.


    Parameters
    ----------
    input : :class:`Input`
        Input parameters.

    Returns
    -------
    out : :class:`Output`
    """
    br0 = input.br
    nr = input.grid.nr
    ns = input.grid.ns
    nphi = input.grid.nphi
    rss = input.grid.rss

    # Coordinates:
    ds = input.grid.ds
    dp = input.grid.dp
    dr = input.grid.dr

    rg = input.grid.rg
    rc = input.grid.rc

    sg = input.grid.sg
    sc = input.grid.sc

    k = np.linspace(0, nr, nr + 1)

    Fp = sg * 0  # Lp/Ls on p-ribs
    Fp[1:-1] = np.sqrt(1 - sg[1:-1]**2) / (np.arcsin(sc[1:]) - np.arcsin(sc[:-1])) * dp
    Vg = Fp / ds / dp
    Fs = (np.arcsin(sg[1:]) - np.arcsin(sg[:-1])) / np.sqrt(1 - sc**2) / dp  # Ls/Lp on s-ribs
    Uc = Fs / ds / dp

    # FFT in phi of photospheric distribution at each latitude:
    brt = np.fft.rfft(br0, axis=1)
    brt = brt.astype(np.complex128)

    # Prepare tridiagonal matrix:
    # - create off-diagonal part of the matrix:
    A = np.zeros((ns, ns))
    for j in range(ns - 1):
        A[j, j+1] = -Vg[j+1]
        A[j+1, j] = A[j, j+1]
    # - term required for m-dependent part of matrix:
    mu = np.fft.fftfreq(nphi)
    mu = 4 * np.sin(np.pi * mu)**2
    # - initialise:
    psir = np.zeros((nr + 1, ns), dtype='complex')
    psi = np.zeros((nr + 1, ns, nphi), dtype='complex')
    e1 = np.exp(dr)
    fact = np.sinh(dr) * (e1 - 1)

    # Loop over azimuthal modes (positive m):
    for m in range(nphi // 2 + 1):
        # - set diagonal terms of matrix:
        A = _A_diag(A, ns, Vg, Uc, mu, m)

        # - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
        #   (note that A is symmetric so use special solver)
        lam, Q = _eigh(A)
        Q = Q.astype(np.complex128)
        # - solve quadratic:
        Flm = 0.5 * (1 + e1 + lam * fact)
        ffp = Flm + np.sqrt(Flm**2 - e1)
        ffm = e1 / ffp

        # - compute radial term for each l (for this m):
        psi, psir = _compute_r_term(m, k, ns, Q, brt, lam, ffm, nr, ffp, psi, psir)

        if (m > 0):
            psi[:, :, nphi - m] = np.conj(psi[:, :, m])

    # Past this point only psi, Fs, Fp are needed
    # Compute psi by inverse fft:
    psi = np.real(np.fft.ifft(psi, axis=2))

    # Hence compute vector potential [note index order, for netcdf]:
    alr = np.zeros((nphi + 1, ns + 1, nr))
    als = np.zeros((nphi + 1, ns, nr + 1))
    alp = np.zeros((nphi, ns + 1, nr + 1))

    als, alp = _als_alp(nr, nphi, Fs, psi, Fp, als, alp)

    r = np.exp(rg)
    th = np.arccos(sg)
    ph = np.linspace(0, 2 * np.pi, nphi + 1)

    return Output(alr, als, alp, input.grid, input.map)
