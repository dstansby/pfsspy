"""
Dipole source solution
======================

A simple example showing how to use pfsspy to compute the solution to a dipole
source field.
"""

###############################################################################
# First, import required modules
import astropy.constants as const
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import numpy as np
import pfsspy


###############################################################################
# Set up a 1degree by 1degree grid in theta and phi
nphi = 360
ntheta = 180

phi = np.linspace(0, 2 * np.pi, nphi)
theta = np.linspace(-np.pi / 2, np.pi / 2, ntheta)
theta, phi = np.meshgrid(theta, phi)

###############################################################################
# Define the number of radial grid points and the source surface radius
nr = 50
rss = 2.5


###############################################################################
# Compute radial component ofa dipole field
def dipole_Br(r, theta):
    return 2 * np.sin(theta) / r**3


br = dipole_Br(1, theta).T

###############################################################################
# Create PFSS input object
input = pfsspy.Input(br, nr, rss)

###############################################################################
# Plot input magnetic field
fig, ax = plt.subplots()
mesh = input.plot_input(ax)
fig.colorbar(mesh)
ax.set_title('Input dipole field')

###############################################################################
# Calculate PFSS solution
output = pfsspy.pfss(input)

###############################################################################
# Plot output field
fig, ax = plt.subplots()
mesh = output.plot_source_surface(ax)
fig.colorbar(mesh)
output.plot_pil(ax)
ax.set_title('Source surface magnetic field')

###############################################################################
# Trace some field lines
br, btheta, bphi = output.bg

fig, ax = plt.subplots()
ax.set_aspect('equal')

# Take 32 start points spaced equally in theta
r = 1.01
for theta in np.linspace(0, np.pi, 33):
    x0 = np.array([0, r * np.sin(theta), r * np.cos(theta)])
    field_line = output.trace(x0)
    color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
    ax.plot(field_line.y / const.R_sun,
            field_line.z / const.R_sun, color=color)

# Add inner and outer boundary circles
ax.add_patch(mpatch.Circle((0, 0), 1, color='k', fill=False))
ax.add_patch(mpatch.Circle((0, 0), input.rss, color='k', linestyle='--',
                           fill=False))
ax.set_title('PFSS solution for a dipole source field')
plt.show()
