import numpy as n
from scipy.io import netcdf
from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.integrate import ode
import matplotlib.pyplot as plt


def plot_2d():
    """
    Script for reading magnetic field from netcdf file, tracing some magnetic field
    lines, and plotting the result in the plane-of-sky.

    Note that the netcdf file must be of grid point type (i.e. br, bth, bph all
    co-located at the mesh points).
    """
    # SPECIFY B FILE:
    bfile = 'bPF.nc'

    # SPECIFY DIRECTION OF VIEWER:
    lon0 = n.deg2rad(290.7)
    lat0 = n.deg2rad(6.93)

    # SPECIFY MAXIMUM OF COLOUR SCALE FOR Br ON SOLAR SURFACE:
    bmax = 20.0

    # SPECIFY FIELD LINE START POINTS:
    # - these should be 1d arrays in theta, phi and r (usual spherical coordinates):
    nth = 32
    th00 = n.linspace(0.02*n.pi, 0.98*n.pi, nth)
    th0 = n.concatenate((th00, th00), axis=0)
    ph0 = n.concatenate((th00*0 + (lon0 + 0.5*n.pi), th00*0 + (lon0 - 0.5*n.pi)), axis=0)
    r0 = th0*0 + 2.5

    #-----------------------------------------------------------------------
    def bTrace(t, x):
        """
        Return B/|B| for use by the field line tracer.
        """
        # (ph,s,rh) coordinates of current point:
        ph = (n.arctan2(x[1], x[0]) + 2*n.pi) % (2*n.pi)
        r = n.sqrt(n.sum(x**2))
        s = x[2]/r
        rh = n.log(r)
        b1 = brgi( n.stack((ph, s, rh)) )
        return b1/n.linalg.norm(b1)

    def trace(x0, dtf=1e-2, tol=1e-2, nrefine=3):
        """
        Trace the fieldline starting from x0, using scipy.integrate.ode.
        - uses implicit Adams method (up to order 12).
        - the parameter dtf is the maximum step-size, which
        will be the output resolution of the field line in most of the domain.
        - the tolerance for tracing is tol*dt.
        - nrefine is the number of times to refine the step-size to get close
        to the boundary.
        """
        xl = x0.copy()

        # Backwards:
        t = 0.0
        dt = dtf
        for j in range(nrefine):
            solver = ode(bTrace).set_integrator('vode', method='adams', atol=tol*dt)
            solver.set_initial_value(xl[:,0:1], t)
            while True:
                try:
                    solver.integrate(solver.t - dt)
                    xl = n.insert(xl, [0], solver.y, axis=1)
                except ValueError: # reached boundary
                    break
            t = solver.t
            dt /= 10.0

        # Forwards:
        t = 0.0
        dt = dtf
        for j in range(nrefine):
            solver = ode(bTrace).set_integrator('vode', method='adams', atol=tol*dt)
            solver.set_initial_value(xl[:,-1:], t)
            while True:
                try:
                    solver.integrate(solver.t + dt)
                    xl = n.append(xl, solver.y, axis=1)
                except ValueError: # reached boundary
                    break
            t = solver.t
            dt /= 10.0
        return xl[0,:], xl[1,:], xl[2,:]

    #-----------------------------------------------------------------------
    # COMPUTE FIELD LINES:
    # - read in magnetic field:
    fh = netcdf.netcdf_file(bfile, 'r', mmap=False)
    r = fh.variables['r'][:]
    th = fh.variables['th'][:]
    ph = fh.variables['ph'][:]
    br = fh.variables['br'][:]
    bth = fh.variables['bth'][:]
    bph = fh.variables['bph'][:]
    fh.close()
    # - (rho,s,phi) coordinates:
    rh = n.log(r)
    s = n.cos(th)
    # - convert to Cartesian components and make interpolator on (rho,s,phi) grid:
    ph3, s3, rh3 = n.meshgrid(ph, s, rh, indexing='ij')
    bx = n.sqrt(1-s3**2)*n.cos(ph3)*br + s3*n.cos(ph3)*bth - n.sin(ph3)*bph
    by = n.sqrt(1-s3**2)*n.sin(ph3)*br + s3*n.sin(ph3)*bth + n.cos(ph3)*bph
    bz = s3*br - n.sqrt(1-s3**2)*bth
    del(br, bth, bph)
    bstack = n.stack((bx,by,bz),axis=3)
    del(bx, by, bz)
    brgi = rgi((ph, s, rh), bstack)
    del(bstack)
    # - convert starting points to Cartesian coordinates:
    x0 = n.stack((r0*n.cos(ph0)*n.sin(th0), r0*n.sin(ph0)*n.sin(th0), r0*n.cos(th0)), axis=0)

    # MAKE PROJECTION OF Br ON SOLAR SURFACE:
    fh = netcdf.netcdf_file(bfile, 'r', mmap=False)
    th = fh.variables['th'][:]
    ph = fh.variables['ph'][:]
    br = fh.variables['br'][:][:,:,0]
    fh.close()
    # - make regular grid interpolator:
    s = n.cos(th)
    br2 = rgi((ph,s), br, bounds_error=False, fill_value=-bmax)
    # - project on to sphere:
    xx = n.linspace(-1, 1, 256)
    y2, z2 = n.meshgrid(xx, xx, indexing='ij')
    r2 = n.sqrt(y2**2 + z2**2)
    x2 = y2*0
    x2[r2 <= 1] = n.sqrt(1.0 - r2[r2 <= 1]**2)
    x2, z2 = n.cos(lat0)*x2 - n.sin(lat0)*z2, n.sin(lat0)*x2 + n.cos(lat0)*z2
    x2, y2 = n.cos(lon0)*x2 - n.sin(lon0)*y2, n.sin(lon0)*x2 + n.cos(lon0)*y2
    x2[r2 > 1] = 0
    y2[r2 > 1] = 0
    sp = z2
    php = (n.arctan2(y2, x2) + 2*n.pi) % (2*n.pi)
    brp = br2(n.stack((php, sp), axis=2))
    brp[r2 > 1] = -bmax

    # PLOT:
    # - set up figure:
    plt.figure(figsize=(6,6))
    ax = plt.subplot(111)
    # ax.set_axis_bgcolor('k')
    cmap0 = plt.cm.get_cmap('gray')
    cmap = plt.cm.get_cmap('bwr')
    # - plot br on the solar surface:
    lev = n.linspace(-bmax, bmax, 128)
    y2, z2 = n.meshgrid(xx, xx, indexing='ij')
    plt.contourf(y2, z2, brp, lev, cmap=cmap0, extend='both')
    # - trace and plot field lines:
    nl = n.size(r0)
    for j in range(nl):
        xl, yl, zl = trace(x0[:,j:j+1])
        #   - rotate to correct viewing direction:
        xl, yl = xl*n.cos(lon0) + yl*n.sin(lon0), -xl*n.sin(lon0) + yl*n.cos(lon0)
        xl, zl = xl*n.cos(lat0) + zl*n.sin(lat0), -xl*n.sin(lat0) + zl*n.cos(lat0)
        #   - remove points that are behind the Sun:
        rl = n.sqrt(yl**2 + zl**2)
        ind = (rl >= 1) | (xl > 0)
        #   - plot this line:
        ax.plot(yl[ind], zl[ind], linewidth=1)
    # - tidy plot window:
    rmax = 2.5
    ax.set_xlim(-rmax, rmax)
    ax.set_ylim(-rmax, rmax)
    ax.text(1.5, -2.25, '$\phi = $%4.2f$^\circ$' % (lon0*180/n.pi), family='serif', color='white', fontsize=9)
    ax.text(1.5, -2.42, '$\lambda = $%4.2f$^\circ$' % (lat0*180/n.pi), family='serif', color='white', fontsize=9)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', \
                    left='off', right='off', labelbottom='off', \
                    labelleft='off')
    # - save figure to file:
    # plt.savefig('fl.png', bbox_inches='tight')
    # - display on screen:
    plt.show()
