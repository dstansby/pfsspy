import astropy.coordinates as coord
import astropy.constants as const
import astropy.units as u
from sunpy.coordinates import frames

import numpy as np
import scipy.linalg as la
import pfsspy.plot


class Grid:
    """
    Grid on which the solution is calculated.

    The grid is evenly spaced in (cos(theta), phi, log(r)).
    """
    def __init__(self, ns, nphi, nr, rss):
        self.ns = ns
        self.nphi = nphi
        self.nr = nr
        self.rss = rss

    @property
    def ds(self):
        """
        Cell size in cos(theta).
        """
        return 2.0 / self.ns

    @property
    def dr(self):
        """
        Cell size in log(r).
        """
        return np.log(self.rss) / self.nr

    @property
    def dp(self):
        """
        Cell size in phi.
        """
        return 2 * np.pi / self.nphi

    @property
    def rc(self):
        """
        Location of the centre of cells in log(r).
        """
        return np.linspace(0.5 * self.dr, np.log(self.rss) - 0.5 * self.dr, self.nr)

    @property
    def sc(self):
        """
        Location of the centre of cells in cos(theta).
        """
        return np.linspace(-1 + 0.5 * self.ds, 1 - 0.5 * self.ds, self.ns)

    @property
    def sp(self):
        """
        Location of the centre of cells in phi.
        """
        return np.linspace(0.5 * self.dp, 2 * np.pi - 0.5 * self.dp, self.nphi)

    @property
    def rg(self):
        """
        Location of the edges of grid cells in log(r).
        """
        return np.linspace(0, np.log(self.rss), self.nr + 1)

    @property
    def sg(self):
        """
        Location of the edges of grid cells in cos(theta).
        """
        return np.linspace(-1, 1, self.ns + 1)

    @property
    def pg(self):
        """
        Location of the edges of grid cells in phi.
        """
        return np.linspace(0, 2 * np.pi, self.nphi + 1)


class Input:
    """
    Input to PFSS modelling.

    Parameters
    ----------
    br0 : array
        Boundary condition of radial magnetic field at the inner surface.

    nr : int
        Number of cells in the radial direction.

    ns : int
        Number of cells in the polar direction.

    np : int
        Number of cells in the azimuthal direction.

    rss : float
        Radius of the source surface, as a fraction of the solar radius.
    """
    def __init__(self, br, nr, rss):
        self.br = br
        ns = self.br.shape[0]
        nphi = self.br.shape[1]
        self.grid = Grid(ns, nphi, nr, rss)

    def plot_input(self, ax=None):
        """
        Plot a 2D image of the magnetic field boundary condition.

        Parameters
        ----------
        ax : Axes
            Axes to plot to. If ``None``, creates a new figure.
        """
        mesh = pfsspy.plot.radial_cut(self.grid.sp, self.grid.sc, self.br, ax)
        return mesh


class Output:
    '''
    Output of PFSS modelling.

    Parameters
    ----------
    r :
        Radial cell edges

    th :
        Elevation cell edges

    ph :
        Azimuthal cell edges

    alr :
        Vector potential * grid spacing in radial direction.

    als :
        Vector potential * grid spacing in elevation direction.

    alp :
        Vector potential * grid spacing in azimuth direction.

    '''
    def __init__(self, r, th, ph, alr, als, alp, input):
        self.r = r
        self.th = th
        self.ph = ph
        self._alr = alr
        self._als = als
        self._alp = alp
        self.grid = input.grid
        self.input = input

        # Cache attributes
        self._common_b_cache = None
        self._rgi = None
        self._bg = None

    def plot_source_surface(self, ax=None):
        """
        Plot a 2D image of the magnetic field at the source surface.

        Parameters
        ----------
        ax : Axes
            Axes to plot to. If ``None``, creates a new figure.
        """
        br, _, _ = self.bg
        mesh = pfsspy.plot.radial_cut(
            self.ph, np.cos(self.th), br[:, :, -1].T, ax)
        return mesh

    def plot_pil(self, ax=None):
        """
        Plot the polarity inversion line on the source surface.

        The PIL is where Br = 0.

        Parameters
        ----------
        ax : Axes
            Axes to plot to. If ``None``, creates a new figure.
        """
        br, _, _ = self.bg
        phi, theta = np.meshgrid(self.ph, np.cos(self.th))
        ax.contour(np.rad2deg(phi), theta, br[:, :, -1].T, levels=[0])

    @property
    def _brgi(self):
        """
        Regular grid interpolator for B.
        """
        from scipy.interpolate import RegularGridInterpolator as rgi
        if self._rgi is not None:
            return self._rgi

        # - (rho,s,phi) coordinates:
        rho = np.log(self.r)
        s = np.cos(self.th)
        phi = self.ph
        br, bth, bph = self.bg

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

    def _bTrace(self, t, x):
        """
        Return B/|B| for use by the field line tracer.
        """
        Bx, By, Bz = x
        # (ph, s, rh) coordinates of current point:
        phi = (np.arctan2(By, Bx) + 2 * np.pi) % (2 * np.pi)
        r = np.linalg.norm(x)
        s = Bz / r  # = cos(theta)
        rh = np.log(r)
        b1 = self._brgi(np.stack((phi, s, rh)))
        return b1 / np.linalg.norm(b1)

    def trace(self, x0, dtf=1e-2, tol=1e-2, nrefine=1):
        """
        Traces a field-line from x0.

        Uses sing scipy.integrate.ode, with an implicit Adams method
        (up to order 12).

        Parameters
        ----------
        x0 : array
            Starting coordinate, in cartesian coordinates.
        dtf : float, optional
            The maximum step-size, which will be the output resolution of the
            field line in most of the domain.
        tol : float, optional
            Relative tolerance of the tracing.
            Absolute tolerance is calculated as ``tol * dtf``.
        nrefine : int, optional
            Number of times to refine the step size near the boundary. Each
            refinement adds extra points to the ends field line to try and
            get close to the inner and outer boundaries.

        Returns
        -------
        fl : :class:`FieldLine`
        """
        from scipy.integrate import ode

        def integrate(dt, start_point):
            t = 0.0
            xout = np.atleast_2d(start_point.copy())
            solver = ode(self._bTrace).set_integrator(
                'vode', method='adams', atol=tol * np.abs(dt))
            solver.set_initial_value(start_point, t)
            while True:
                try:
                    solver.integrate(solver.t + dt)
                    if dt < 0:
                        xout = np.row_stack((solver.y, xout))
                    else:
                        xout = np.row_stack((xout, solver.y))
                except ValueError as e:  # reached boundary
                    if 'One of the requested xi is out of bounds' in str(e):
                        break
                    raise e
            return xout

        xback = integrate(-dtf, x0)
        xforw = integrate(dtf, x0)
        for i in range(nrefine):
            dtf /= 10
            xback = np.row_stack((integrate(-dtf, xback[0, :]), xback))
            xforw = np.row_stack((xforw, integrate(dtf, xforw[-1, :])))

        xout = np.row_stack((xback, xforw))
        return FieldLine(xout[:, 0], xout[:, 1], xout[:, 2], self)

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
        for i in range(self.grid.nphi + 2):
            br[i, :, :] = br[i, :, :] / Sbr
            bs[i, :, :] = bs[i, :, :] / Sbs
        for i in range(self.grid.nphi + 1):
            bp[i, :, :] = bp[i, :, :] / Sbp

        return br, -bs, bp

    @property
    def bg(self):
        """
        B as a (weighted) averaged on grid points.
        """
        if self._bg is not None:
            return self._bg

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

        self._bg = brg, bsg, bpg
        return self._bg

    def save_a(self, fname):
        """
        Save vector potential * edge lengths to a file.
        """
        from scipy.io import netcdf
        r = self.r
        th = self.th
        ph = self.ph
        apr, aps, app = self.al

        nr = np.size(r) - 1
        ns = np.size(th) - 1
        nphi = np.size(ph) - 1

        fid = netcdf.netcdf_file(fname, 'w')
        fid.createDimension('rc', nr)
        fid.createDimension('r', nr + 1)
        fid.createDimension('thc', ns)
        fid.createDimension('th', ns + 1)
        fid.createDimension('phc', nphi)
        fid.createDimension('ph', nphi + 1)
        vid = fid.createVariable('r', 'd', ('r',))
        vid[:] = r
        vid = fid.createVariable('th', 'd', ('th',))
        vid[:] = th
        vid = fid.createVariable('ph', 'd', ('ph',))
        vid[:] = ph
        vid = fid.createVariable('ar', 'd', ('ph', 'th', 'rc'))
        vid[:] = apr
        vid = fid.createVariable('as', 'd', ('ph', 'thc', 'r'))
        vid[:] = aps
        vid = fid.createVariable('ap', 'd', ('phc', 'th', 'r'))
        vid[:] = app
        fid.close()
        print('Wrote A*L to file ' + fname)

    def save_bg(self, fname):
        """
        Save magnetic field components co-located at grid points.
        """
        from scipy.io import netcdf
        r = self.r
        th = self.th
        ph = self.ph
        brg, bsg, bpg = self.bg

        nr = np.size(r) - 1
        ns = np.size(th) - 1
        nphi = np.size(ph) - 1

        fid = netcdf.netcdf_file(fname, 'w')
        fid.createDimension('r', nr + 1)
        fid.createDimension('th', ns + 1)
        fid.createDimension('ph', nphi + 1)
        vid = fid.createVariable('r', 'd', ('r',))
        vid[:] = r
        vid = fid.createVariable('th', 'd', ('th',))
        vid[:] = th
        vid = fid.createVariable('ph', 'd', ('ph',))
        vid[:] = ph
        vid = fid.createVariable('br', 'd', ('ph', 'th', 'r'))
        vid[:] = brg
        vid = fid.createVariable('bth', 'd', ('ph', 'th', 'r'))
        vid[:] = -bsg
        vid = fid.createVariable('bph', 'd', ('ph', 'th', 'r'))
        vid[:] = bpg
        fid.close()
        print('Wrote B at grid points to file ' + fname)

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
            bs[i,-1,:] = 0.5*(bs[i,-2,:] - bs[i1,-2,:])
            bs[i,0,:] = 0.5*(bs[i,1,:] - bs[i1,1,:])
        for i in range(nphi + 1):
            i1 = (i + nphi // 2) % nphi
            bp[i, -1, :] = -bp[i1, -2, :]
            bp[i, 0, :] = -bp[i1, 1, :]

        self._common_b_cache = br, bs, bp, Sbr, Sbs, Sbp
        return self._common_b_cache


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
        for j in range(ns):
            A[j, j] = Vg[j] + Vg[j + 1] + Uc[j] * mu[m]
        # - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
        #   (note that A is symmetric so use special solver)
        lam, Q = la.eigh(A)
        # - solve quadratic:
        Flm = 0.5 * (1 + e1 + lam * fact)
        ffp = Flm + np.sqrt(Flm**2 - e1)
        ffm = e1 / ffp
        # - compute radial term for each l (for this m):
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
        if (m > 0):
            psi[:, :, nphi - m] = np.conj(psi[:, :, m])

    # Compute psi by inverse fft:
    psi = np.real(np.fft.ifft(psi, axis=2))

    # Hence compute vector potential [note index order, for netcdf]:
    alr = np.zeros((nphi + 1, ns + 1, nr))
    als = np.zeros((nphi + 1, ns, nr + 1))
    alp = np.zeros((nphi, ns + 1, nr + 1))

    for j in range(nr + 1):
        for i in range(nphi + 1):
            als[i, :, j] = Fs * (psi[j, :, ((i - 1) % nphi)] - psi[j, :, ((i) % nphi)])
        for i in range(nphi):
            alp[i, 1:-1, j] = Fp[1:-1] * (psi[j, 1:, i] - psi[j, :-1, i])

    # Output to netcdf file:
    r = np.exp(rg)
    th = np.arccos(sg)
    ph = np.linspace(0, 2 * np.pi, nphi + 1)

    return Output(r, th, ph, alr, als, alp, input)


class FieldLine(coord.SkyCoord):
    """
    A single magnetic field line.

    Parameters
    ----------
    x :
    y :
    z :
    output : :class:`Output`
    """
    def __init__(self, x, y, z, output):
        self._output = output
        super().__init__(x=x * const.R_sun,
                         y=y * const.R_sun,
                         z=z * const.R_sun,
                         frame=frames.HeliographicCarrington,
                         representation_type='cartesian')

    @property
    def is_open(self):
        """
        Returns ``True`` if one of the field line is connected to the solar
        surface and one to the outer boundary, ``False`` otherwise.
        """
        r = coord.SkyCoord(self, representation_type='spherical').radius
        rtol = 0.1
        if np.abs(r[0] - r[-1]) < r[0] * rtol:
            return False
        return True

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
        if not self.is_open:
            return 0
        # Because the field lines are integrated forwards, if the end point
        # is on the outer boundary the field is outwards
        foot1 = coord.SkyCoord(self[0], representation_type='spherical')
        foot2 = coord.SkyCoord(self[-1], representation_type='spherical')
        if foot2.radius - foot1.radius > 0:
            return 1
        else:
            return -1

    @property
    def expansion_factor(self):
        r"""
        Magnetic field expansion factor.

        The expansion factor is defnied as
        :math:`(r_{\odot} B_{\odot} / (r_{ss} B_{ss}))^{2}`

        Returns
        -------
        exp_fact : float
            Field line expansion factor.
            If field line is closed, returns ``None``.
        """
        if not self.is_open:
            return
        # Extract ends of magnetic field line, and get them in spherical coords
        foot1 = coord.SkyCoord(self[0], representation_type='spherical')
        foot2 = coord.SkyCoord(self[-1], representation_type='spherical')
        if foot1.radius > foot2.radius:
            solar_foot = foot2
            source_foot = foot1
        else:
            solar_foot = foot1
            source_foot = foot2

        def b_at_coord(coord):
            interp_coords = [coord.lon / u.rad,
                             np.sin(coord.lat),
                             np.log(coord.radius / const.R_sun)]
            b = self._output._brgi(interp_coords)
            return np.linalg.norm(b)

        b_solar = b_at_coord(solar_foot)
        b_source = b_at_coord(source_foot)
        return ((solar_foot.radius * b_solar) / (source_foot.radius * b_source))**2
