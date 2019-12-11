"""
Open/closed field map
=====================

Creating an open/closed field map on the solar surface.
"""

###############################################################################
# First, import required modules
import os
import astropy.constants as const
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor

import numpy as np
import sunpy.map

import pfsspy
from pfsspy import coords
from pfsspy import tracing


###############################################################################
# Load a GONG magnetic field map. If 'gong.fits' is present in the current
# directory, just use that, otherwise download a sample GONG map.
if not os.path.exists('190310t0014gong.fits') and not os.path.exists('190310t0014gong.fits.gz'):
    import urllib.request
    urllib.request.urlretrieve(
        'https://gong2.nso.edu/oQR/zqs/201903/mrzqs190310/mrzqs190310t0014c2215_333.fits.gz',
        '190310t0014gong.fits.gz')

if not os.path.exists('190310t0014gong.fits'):
    import gzip
    with gzip.open('190310t0014gong.fits.gz', 'rb') as f:
        with open('190310t0014gong.fits', 'wb') as g:
            g.write(f.read())

###############################################################################
# We can now use SunPy to load the GONG fits file, and extract the magnetic
# field data.
#
# The mean is subtracted to enforce div(B) = 0 on the solar surface: n.b. it is
# not obvious this is the correct way to do this, so use the following lines
# at your own risk!
[[br, header]] = sunpy.io.fits.read('190310t0014gong.fits')
br = br - np.mean(br)
###############################################################################
# GONG maps have their LH edge at -180deg in Carrington Longitude,
# so roll to get it at 0deg. This way the input magnetic field is in a
# Carrington frame of reference, which matters later when lining the field
# lines up with the AIA image.
br = np.roll(br, header['CRVAL1'] + 180, axis=1)


###############################################################################
# The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
# rho = ln(r), and r is the standard spherical radial coordinate. We need to
# define the number of rho grid points, and the source surface radius.
nrho = 60
rss = 2.5

###############################################################################
# From the boundary condition, number of radial grid points, and source
# surface, we now construct an Input object that stores this information
input = pfsspy.Input(br, nrho, rss)
output = pfsspy.pfss(input)


###############################################################################
# Finally, using the 3D magnetic field solution we can trace some field lines.
# In this case a grid of 90 x 180 points equally gridded in theta and phi are
# chosen and traced from the source surface outwards.
#
# First, set up the tracing seeds

r = 1
# Number of steps in cos(latitude)
nsteps = 90
phi_1d = np.linspace(0, 2 * np.pi, nsteps * 2 + 1)
theta_1d = np.arccos(np.linspace(-1, 1, nsteps + 1))
phi, theta = np.meshgrid(phi_1d, theta_1d, indexing='ij')
seeds = np.array(coords.sph2cart(r, theta.ravel(), phi.ravel())).T

###############################################################################
# Trace the field lines
print('Tracing field lines...')
tracer = tracing.FortranTracer(1000, 0.01)
field_lines = tracer.trace(seeds, output)
print('Finished tracing field lines')

###############################################################################
# Plot the result. The to plot is the input magnetogram, and the bottom plot
# shows a contour map of the the footpoint polarities, which are +/- 1 for open
# field regions and 0 for closed field regions.
fig, axs = plt.subplots(nrows=2, sharex=True, sharey=True)
ax = axs[0]
input.plot_input(ax)
ax.set_title('Input GONG magnetogram')

ax = axs[1]
cmap = mcolor.ListedColormap(['tab:red', 'black', 'tab:blue'])
norm = mcolor.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ncolors=3)
pols = field_lines.polarities.reshape(2 * nsteps + 1, nsteps + 1).T
ax.contourf(np.rad2deg(phi_1d), np.cos(theta_1d), pols, norm=norm, cmap=cmap)

ax.set_title('Open (blue/red) and closed (black) field')
ax.set_aspect(0.5 * 360 / 2)
plt.show()
