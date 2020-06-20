"""
Overplotting field lines on AIA maps
====================================

This example shows how to take a PFSS solution, trace some field lines, and
overplot the traced field lines on an AIA 193 map.
"""

###############################################################################
# First, we import the required modules
from datetime import datetime
import os

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.io.fits

import pfsspy
import pfsspy.coords as coords
import pfsspy.tracing as tracing
from pfsspy.sample_data import get_gong_map


###############################################################################
# Load a GONG magnetic field map. If 'gong.fits' is present in the current
# directory, just use that, otherwise download a sample GONG map.
gong_fname = get_gong_map()

###############################################################################
# We can now use SunPy to load the GONG fits file, and extract the magnetic
# field data.
#
# The mean is subtracted to enforce div(B) = 0 on the solar surface: n.b. it is
# not obvious this is the correct way to do this, so use the following lines
# at your own risk!
gong_map = sunpy.map.Map(gong_fname)
# Remove the mean
gong_map = sunpy.map.Map(gong_map.data - np.mean(gong_map.data), gong_map.meta)

###############################################################################
# Load the corresponding AIA 193 map
if not os.path.exists('AIA20190310.fits'):
    import urllib.request
    urllib.request.urlretrieve(
        'http://jsoc2.stanford.edu/data/aia/synoptic/2019/03/10/H0000/AIA20190310_0000_0193.fits',
        'AIA20190310.fits')

aia = sunpy.map.Map('AIA20190310.fits')
dtime = aia.date

###############################################################################
# The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
# rho = ln(r), and r is the standard spherical radial coordinate. We need to
# define the number of grid points in rho, and the source surface radius.
nrho = 25
rss = 2.5

###############################################################################
# From the boundary condition, number of radial grid points, and source
# surface, we now construct an `Input` object that stores this information
input = pfsspy.Input(gong_map, nrho, rss)

###############################################################################
# Using the `Input` object, plot the input photospheric magnetic field
m = input.map
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot()
plt.colorbar()
ax.set_title('Input field')

###############################################################################
# We can also plot the AIA map to give an idea of the global picture. There
# is a nice active region in the top right of the AIA plot, that can also
# be seen in the top left of the photospheric field plot above.
ax = plt.subplot(1, 1, 1, projection=aia)
aia.plot(ax)


###############################################################################
# Now we construct a 5 x 5 grid of footpoitns to trace some magnetic field
# lines from. These coordinates are defined in the native Carrington
# coordinates of the input magnetogram.

# Create 5 points spaced between sin(lat)={0.1, 0.2}
s = np.linspace(0.1, 0.2, 5)
# Create 5 points spaced between long={55, 65} degrees
phi = np.linspace(55, 65, 5)
print(f's = {s}')
print(f'phi = {phi}')
# Make a 2D grid from these 1D points
s, phi = np.meshgrid(s, phi)

# Now convert the points to a coordinate object
lat = np.arcsin(s) * u.rad
lon = phi * u.deg
seeds = SkyCoord(lon.ravel(), lat.ravel(), 1.01 * const.R_sun,
                 frame=gong_map.coordinate_frame)

###############################################################################
# Plot the magnetogram and the seed footpoints The footpoints are centered
# around the active region metnioned above.
m = input.map
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot()
plt.colorbar()

ax.plot_coord(seeds, color='black', marker='o', linewidth=0, markersize=2)

# Set the axes limits. These limits have to be in pixel values
ax.set_xlim(0, 180)
ax.set_ylim(45, 135)
ax.set_title('Field line footpoints')
ax.set_ylim(bottom=0)

#######################################################################
# Compute the PFSS solution from the GONG magnetic field input
output = pfsspy.pfss(input)

###############################################################################
# Trace field lines from the footpoints defined above.
tracer = tracing.PythonTracer()
flines = tracer.trace(seeds, output)

###############################################################################
# Plot the input GONG magnetic field map, along with the traced mangetic field
# lines.
m = input.map
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot()
plt.colorbar()

for fline in flines:
    ax.plot_coord(fline.coords, color='black', linewidth=1)

# Set the axes limits. These limits have to be in pixel values
ax.set_xlim(0, 180)
ax.set_ylim(45, 135)
ax.set_title('Photospheric field and traced field lines')
###############################################################################
# Plot the AIA map, along with the traced magnetic field lines. Inside the
# loop the field lines are converted to the AIA observer coordinate frame,
# and then plotted on top of the map.
fig = plt.figure()
ax = plt.subplot(1, 1, 1, projection=aia)
aia.plot(ax)
for fline in flines:
    ax.plot_coord(fline.coords, alpha=0.8, linewidth=1, color='black')

ax.set_xlim(500, 900)
ax.set_ylim(400, 800)
plt.show()

# sphinx_gallery_thumbnail_number = 5
