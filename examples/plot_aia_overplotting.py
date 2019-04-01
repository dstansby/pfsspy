"""
GONG PFSS extrapolation
=======================

Calculating PFSS solution for a GONG synoptic magnetic field map.
"""

###############################################################################
# First, import required modules
import os
from datetime import datetime
import astropy.units as u
import astropy.constants as const
import astropy.coordinates
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sunpy.map
import sunpy.io.fits
from sunpy.coordinates import frames

import pfsspy
import pfsspy.coords as coords


###############################################################################
# Load a GONG magnetic field map. If 'gong.fits' is present in the current
# directory, just use that, otherwise download a sample GONG map.
#
# The map date is 06/03/2019
if not os.path.exists('gong.fits') and not os.path.exists('gong.fits.gz'):
    import urllib.request
    urllib.request.urlretrieve(
        'https://gong2.nso.edu/oQR/zqs/201903/mrzqs190310/mrzqs190310t0014c2215_333.fits.gz',
        'gong.fits.gz')

if not os.path.exists('gong.fits'):
    import gzip
    with gzip.open('gong.fits.gz', 'rb') as f:
        with open('gong.fits', 'wb') as g:
            g.write(f.read())

###############################################################################
# Load an AIA map
if not os.path.exists('aia.fits'):
    import urllib.request
    urllib.request.urlretrieve(
        'http://jsoc2.stanford.edu/data/aia/synoptic/2019/03/10/H0000/AIA20190310_0000_0193.fits',
        'aia.fits')

aia = sunpy.map.Map('aia.fits')
dtime = datetime(2019, 3, 10)

###############################################################################
# We can now use SunPy to load the .fits file, and extract the magnetic field
# data.
#
# The mean is subtracted to enforce div(B) = 0 on the solar surface: n.b. it is
# not obvious this is the correct way to do this, so use the following lines
# at your own risk!
[[br, header]] = sunpy.io.fits.read('gong.fits')
br = br - np.mean(br)
# GONG maps have their LH edge at -180deg, so roll to get it at 0deg
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

###############################################################################
# Using the Input object, plot the input field
fig, ax = plt.subplots()
mesh = input.plot_input(ax)
fig.colorbar(mesh)
ax.set_title('Input field')

ax = plt.subplot(1, 1, 1, projection=aia)
aia.plot(ax)


###############################################################################
# Calculate a 10 x 10 grid of footpoitns to trace magnetic field lines from
# The figure shows a zoom in of the magnetic field map, with the footpoints
# overplotted.
s, phi = np.meshgrid(np.linspace(0.1, 0.2, 10),
                     np.deg2rad(np.linspace(55, 65, 10)))

fig, ax = plt.subplots()
mesh = input.plot_input(ax)
fig.colorbar(mesh)
ax.set_xlim(40, 70)
ax.set_ylim(0, 0.35)
ax.scatter(np.rad2deg(phi), s, color='k', s=1)
ax.set_title('Field line footpoints')


###############################################################################
#
output = pfsspy.pfss(input)
flines = []
for s, phi in zip(s.ravel(), phi.ravel()):
    x0 = np.array(pfsspy.coords.strum2cart(0.01, s, phi))
    flines.append(output.trace(x0, atol=1e-6))


fig, ax = plt.subplots()
mesh = input.plot_input(ax)
ax.set_xlim(60, 90)
ax.set_ylim(0, 0.35)

for fline in flines:
    fline.representation_type = 'spherical'
    ax.plot(fline.lon / u.deg, np.sin(fline.lat), color='black', linewidth=1)


def transform_to_proj(fline, dtime):
    """
    Transform a field line into a projective frame. Assumes the observer is
    Earth.
    """
    proj_frame = frames.Helioprojective(obstime=dtime, observer='earth')
    fline.representation_type = 'spherical'
    fline = astropy.coordinates.SkyCoord(lon=fline.lon, lat=fline.lat, radius=fline.radius,
                                         frame=frames.HeliographicCarrington,
                                         obstime=dtime)
    fline = fline.transform_to(proj_frame)
    return fline


fig = plt.figure()
ax = plt.subplot(1, 1, 1, projection=aia)
aia.plot(ax)
for fline in flines:
    fline = transform_to_proj(fline, dtime)
    Tx = fline.Tx.to(u.deg)
    Ty = fline.Ty.to(u.deg)
    # trans = ax.get_transform('world')
    ax.plot(Tx, Ty, color='w', alpha=0.8,
            transform=ax.get_transform('world'), linewidth=1)

plt.show()
# sphinx_gallery_thumbnail_number = 3
