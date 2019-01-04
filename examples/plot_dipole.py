"""
Dipole source solution
======================

A simple example showing how to use pfsspy to compute the solution to a dipole
source field.
"""

###############################################################################
# First, import required modules
import matplotlib.pyplot as plt
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
nr = 20
rss = 2.5


###############################################################################
# Compute radial component ofa dipole field
def dipole_Br(r, theta):
    return 2 * np.sin(theta) / r**3


br = dipole_Br(1, theta).T

###############################################################################
# Create PFSS input object
input = pfsspy.Input(br, nr, ntheta, nphi, rss)

###############################################################################
# Plot input magnetic field
fig, ax = plt.subplots()
input.plot_input(ax)
ax.set_title('Input dipole field')
plt.show()

###############################################################################
# Calculate PFSS solution
# output = pfsspy.pfss(input)
