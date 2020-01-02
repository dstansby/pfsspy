"""
Open/closed field map
=====================

Creating an open/closed field map on the solar surface.
"""

###############################################################################
# First, import required modules
import os

import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
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
# GONG maps have their LH edge at varying Carrington longitudes,
# so roll to get it at -180deg, which is what the carr_cea_wcs_header function
# expects.
br = np.roll(br, header['CRVAL1'], axis=1)

header = pfsspy.carr_cea_wcs_header(header['DATE'], br.shape)
gong_map = sunpy.map.Map((br, header))


###############################################################################
# Set the model parameters
nrho = 60
rss = 2.5

###############################################################################
# Construct the input, and calculate the output solution
input = pfsspy.Input(gong_map, nrho, rss)
output = pfsspy.pfss(input)


###############################################################################
# Finally, using the 3D magnetic field solution we can trace some field lines.
# In this case a grid of 90 x 180 points equally gridded in theta and phi are
# chosen and traced from the source surface outwards.
#
# First, set up the tracing seeds

r = const.R_sun
# Number of steps in cos(latitude)
nsteps = 90
lon_1d = np.linspace(0, 2 * np.pi, nsteps * 2 + 1)
lat_1d = np.arcsin(np.linspace(-1, 1, nsteps + 1))
lon, lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')
lon, lat = lon*u.rad, lat*u.rad
seeds = SkyCoord(lon.ravel(), lat.ravel(), r, frame=output.coordinate_frame)

###############################################################################
# Trace the field lines
print('Tracing field lines...')
tracer = tracing.FortranTracer(max_steps=2000)
field_lines = tracer.trace(seeds, output)
print('Finished tracing field lines')

###############################################################################
# Plot the result. The to plot is the input magnetogram, and the bottom plot
# shows a contour map of the the footpoint polarities, which are +/- 1 for open
# field regions and 0 for closed field regions.
fig = plt.figure()
m = input.map
ax = fig.add_subplot(2, 1, 1, projection=m)
m.plot()
ax.set_title('Input GONG magnetogram')

ax = fig.add_subplot(2, 1, 2)
cmap = mcolor.ListedColormap(['tab:red', 'black', 'tab:blue'])
norm = mcolor.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ncolors=3)
pols = field_lines.polarities.reshape(2 * nsteps + 1, nsteps + 1).T
ax.contourf(np.rad2deg(lon_1d), np.sin(lat_1d), pols, norm=norm, cmap=cmap)
ax.set_ylabel('sin(latitude)')

ax.set_title('Open (blue/red) and closed (black) field')
ax.set_aspect(0.5 * 360 / 2)

plt.show()
