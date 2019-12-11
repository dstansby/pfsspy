"""
Different field line tracers
============================

pfsspy comes with two built in field line tracers; a python one that works
natively in python, and a fortran one that uses the streamtracer package.

Unsurprisingly the fortran tracer is much faster, but requires a fortran
compiler to install.

In this example we show how easy it is to switch between the two different
tracers, and how much faster the fortran tracer is.
"""
###############################################################################
# First, import required modules
import time

import astropy.constants as const
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import numpy as np
import pfsspy
import pfsspy.coords as coords


###############################################################################
# Setup a grid for the input.
nphi = 360
ns = 180

phi = np.linspace(0, 2 * np.pi, nphi)
s = np.linspace(-1, 1, ns)
s, phi = np.meshgrid(s, phi)


###############################################################################
# Now we can take the grid and calculate the boundary condition magnetic field.
def dipole_Br(r, s):
    return 2 * s / r**3


br = dipole_Br(1, s).T


###############################################################################
# Define the number of radial grid points, and source surface
nrho = 50
rss = 2.5

###############################################################################
# Construct an input object and calculate the PFSS solution
input = pfsspy.Input(br, nrho, rss)
output = pfsspy.pfss(input)

python_tracer = pfsspy.tracing.PythonTracer()
fortran_tracer = pfsspy.tracing.FortranTracer(1000, 0.01)

###############################################################################
# Finally, using the 3D magnetic field solution we can trace some field lines.
# In this case 32 points equally spaced in theta are chosen and traced from
# the source surface outwards.

# Take 32 start points spaced equally in theta
r = 1.01
phi = np.pi / 2
theta = np.linspace(0, np.pi, 33)
seeds = np.array(coords.sph2cart(r, theta, phi)).T

flines = {}
for tracer, name in zip([python_tracer, fortran_tracer],
                        ['python', 'fortran']):
    print(f'Tracing with {name} tracer...')
    t1 = time.time()
    flines[name] = tracer.trace(seeds, output)
    t2 = time.time()
    dt = t2 - t1
    print(f'Tracing took {dt} seconds')

###############################################################################
# Plot the field lines
fig, ax = plt.subplots()
ax.set_aspect('equal')
linestyle = {'python': '-', 'fortran': '--'}

for key in flines:
    for field_line in flines[key]:
        color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
        ax.plot(field_line.coords.y / const.R_sun,
                field_line.coords.z / const.R_sun,
                color=color, linestyle=linestyle[key])

# Add inner and outer boundary circles
ax.add_patch(mpatch.Circle((0, 0), 1, color='k', fill=False))
ax.add_patch(mpatch.Circle((0, 0), input.grid.rss, color='k', linestyle='--',
                           fill=False))
ax.set_title('PFSS solution for a dipole source field')
plt.show()
