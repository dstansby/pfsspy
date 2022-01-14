"""
Analytic dipole field lines
===========================
"""

###############################################################################
# First, import required modules
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import pandas as pd
import numpy as np
import sunpy.map
import pfsspy
import pfsspy.coords as coords
from pfsspy import analytic, tracing

from helpers import pffspy_output, fr, theta_fline_coords, phi_fline_coords

###############################################################################
# Compare the the pfsspy solution to the analytic solutions. Cuts are taken
# on the source surface at a constant phi value to do a 1D comparison.
l = 3
m = 3
nphi = 360
ns = 180
nr = 40
rss = 2


###############################################################################
# Calculate PFSS solution
pfsspy_out = pffspy_output(nphi, ns, nr, rss, l, m)

rss = rss * const.R_sun
###############################################################################
# Trace some field lines
n = 90
# Create 1D theta, phi arrays
phi = np.linspace(0, 360, n * 2)
phi = phi[:-1] + np.diff(phi) / 2
theta = np.arcsin(np.linspace(-0.98, 0.98, n, endpoint=False) + 1/n)
# Mesh into 2D arrays
theta, phi = np.meshgrid(theta, phi, indexing='ij')
theta, phi = theta * u.rad, phi * u.deg
seeds = SkyCoord(radius=rss, lat=theta.ravel(), lon=phi.ravel(),
                 frame=pfsspy_out.coordinate_frame)

step_sizes = [32, 16, 8, 4, 2, 1, 0.5]
dthetas = []
dphis = []
for step_size in step_sizes:
    print(f'Tracing {step_size}...')
    # Trace
    tracer = tracing.FortranTracer(step_size=step_size)
    flines = tracer.trace(seeds, pfsspy_out)
    # Set a mask of open field lines
    mask = flines.connectivities.astype(bool).reshape(theta.shape)

    # Get solar surface latitude
    phi_solar = np.ones_like(phi) * np.nan
    phi_solar[mask] = flines.open_field_lines.solar_feet.lon
    theta_solar = np.ones_like(theta) * np.nan
    theta_solar[mask] = flines.open_field_lines.solar_feet.lat
    r_out = np.ones_like(theta.value) * const.R_sun * np.nan
    r_out[mask] = flines.open_field_lines.solar_feet.radius

    ###########################################################################
    # Calculate analytical solution
    theta_analytic = theta_fline_coords(r_out, rss, l, m, theta)
    dtheta = (theta_solar - theta_analytic).to_value(u.deg)
    phi_analytic = phi_fline_coords(r_out, rss, l, m, theta, phi)
    dphi = (phi_solar - phi_analytic).to_value(u.deg)
    # Wrap phi values
    dphi[dphi > 180] -= 360
    dphi[dphi < -180] += 360

    dthetas.append(dtheta.ravel())
    dphis.append(dphi.ravel())

dthetas = pd.DataFrame(data=np.array(dthetas), index=step_sizes)
dthetas.to_hdf(f'results/dthetas_{l}{m}.hdf', 'table', mode='w')
dphis = pd.DataFrame(data=np.array(dphis), index=step_sizes)
dphis.to_hdf(f'results/dphis_{l}{m}.hdf', 'table', mode='w')
