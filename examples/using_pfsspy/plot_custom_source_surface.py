"""
Dipole source solution: Custom Source Surface
=============================================

A simple example comparing the solutions to a dipole source field when different
source surface outer boundary conditions are provided: radial (default) and closed.
"""
import astropy.constants as const
import astropy.units as u
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.time import Time

import pfsspy

###############################################################################
# To start with we need to construct an input for the PFSS model. To do this,
# first set up a regular 2D grid in (phi, s), where s = cos(theta) and
# (phi, theta) are the standard spherical coordinate system angular
# coordinates. In this case the resolution is (360 x 180).
nphi = 360
ns = 180

phi = np.linspace(0, 2 * np.pi, nphi)
s = np.linspace(-1, 1, ns)
s, phi = np.meshgrid(s, phi)


###############################################################################
# Now we can take the grid and calculate the boundary condition magnetic fields.
# The radial outer boundary condition requiring the magnetic field to be fully
# radial at the source surface will be applied by default, but the closed outer
# boundary condition requiring the radial component to vanish will need to be
# provided as a custom outer boundary condition. Therefore, we must construct
# an input map where the radial component of the magnetic field is everywhere
# zero along with the dipole source field.
def dipole_Br(r, s):
    return 2 * s / r**3


# Dipole inner boundary:
br_dipole = dipole_Br(1, s)
header_dipole = pfsspy.utils.carr_cea_wcs_header(Time('2020-1-1'), br_dipole.shape)
map_dipole = sunpy.map.Map((br_dipole.T, header_dipole))

# Closed outer boundary:
br_zeros = np.zeros((nphi, ns))
header_zeros = pfsspy.utils.carr_cea_wcs_header(Time('2020-1-1'), br_zeros.shape)
map_zeros = sunpy.map.Map((br_zeros.T, header_zeros))

###############################################################################
# The PFSS solutions are calculated on a regular 3D grid in (phi, s, rho), where
# rho = ln(r), and r is the standard spherical radial coordinate. We need to
# define the number of rho grid points, and the source surface radius.
nrho = 30
rss = 2.5

###############################################################################
# From the boundary conditions, number of radial grid points, and source
# surfaces, we now construct Input objects that store this information.
pfss_in_radial = pfsspy.Input(map_dipole, nrho, rss)  # 'radial'
pfss_in_closed = pfsspy.Input(map_dipole, nrho, rss, map_zeros)

###############################################################################
# Now calculate the PFSS solutions.
pfss_out_radial = pfsspy.pfss(pfss_in_radial)
pfss_out_closed = pfsspy.pfss(pfss_in_closed)

###############################################################################
# Finally, using the 3D magnetic field solutions we can trace some field lines.
# In this case 32 points equally spaced in theta are chosen and traced from
# the source surfaces outwards.
fig = plt.figure()

# Add the radial solution subplot
ax = fig.add_subplot(1, 2, 1)
ax.set_aspect('equal')

# Take 32 start points spaced equally in theta
r = 1.01 * const.R_sun
lon = np.pi / 2 * u.rad
lat = np.linspace(-np.pi / 2, np.pi / 2, 33) * u.rad
seeds = SkyCoord(lon, lat, r, frame=pfss_out_radial.coordinate_frame)

tracer = pfsspy.tracing.FortranTracer()
field_lines_radial = tracer.trace(seeds, pfss_out_radial)

for field_line in field_lines_radial:
    coords = field_line.coords
    coords.representation_type = 'cartesian'
    color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
    ax.plot(coords.y / const.R_sun,
            coords.z / const.R_sun, color=color)

# Add inner and outer boundary circles
ax.add_patch(mpatch.Circle((0, 0), 1, color='k', fill=False))
ax.add_patch(mpatch.Circle((0, 0), pfss_in_radial.grid.rss, color='k', linestyle='--',
                           fill=False))
ax.set_title('Radial source surface')

# Add the closed solution subplot
ax = fig.add_subplot(1, 2, 2)
ax.set_aspect('equal')

field_lines_closed = tracer.trace(seeds, pfss_out_closed)

for field_line in field_lines_closed:
    coords = field_line.coords
    coords.representation_type = 'cartesian'
    color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
    ax.plot(coords.y / const.R_sun,
            coords.z / const.R_sun, color=color)

# Add inner and outer boundary circles
ax.add_patch(mpatch.Circle((0, 0), 1, color='k', fill=False))
ax.add_patch(mpatch.Circle((0, 0), pfss_in_closed.grid.rss, color='k', linestyle='--',
                           fill=False))
ax.set_title('Closed source surface')
plt.suptitle('PFSS solution for a dipole source field', fontsize=16)

plt.show()
