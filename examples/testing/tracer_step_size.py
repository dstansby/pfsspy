"""
Tracer step size (calculations)
===============================
"""
import astropy.constants as const
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from helpers import pffspy_output, phi_fline_coords, theta_fline_coords

from pfsspy import tracing

l = 3
m = 3
nphi = 360
ns = 180
nr = 40
rss = 2


###############################################################################
# Calculate PFSS solution
pfsspy_out = pffspy_output(nphi, ns, nr, rss, l, m)

###############################################################################
# Trace an array of field lines from the source surface down to the solar
# surface
n = 90
# Create 1D theta, phi arrays
phi = np.linspace(0, 360, n * 2)
phi = phi[:-1] + np.diff(phi) / 2
theta = np.arcsin(np.linspace(-0.98, 0.98, n, endpoint=False) + 1/n)
# Mesh into 2D arrays
theta, phi = np.meshgrid(theta, phi, indexing='ij')
theta, phi = theta * u.rad, phi * u.deg
rss = rss * const.R_sun
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


###############################################################################
# Save results. This saves the maximum error in both phi and theta as a
# function of thet tracer step size.
dthetas = pd.DataFrame(data=np.array(dthetas), index=step_sizes)
dthetas = dthetas.mask(np.abs(dthetas) > 30).max(axis=1)

dphis = pd.DataFrame(data=np.array(dphis), index=step_sizes)
dphis = dphis.mask(np.abs(dphis) > 30).max(axis=1)

dthetas.to_csv(f'results/flines/dthetas_{l}{m}.csv', header=False)
dphis.to_csv(f'results/flines/dphis_{l}{m}.csv', header=False)
