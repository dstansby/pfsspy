"""
Script for reading in map of Br(theta, phi) on the solar surface, computing a
PFSS extrapolation, and outputting to a netcdf file.
"""
import sys
import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
from pfss import pfss

# DEFINE GRID FOR PFSS COMPUTATION
# - equally spaced in rho=log(r/Rsun), s=cos(theta) and phi.
# - specify number of grid points:
nr = 60     # Radial
ns = 180    # Cos(theta)
nphi = 360    # Phi
if ((nphi % 2) != 0):
    raise RuntimeError(
        'Number of points in phi must be even (for polar boundary conditions)')
# - choose source surface radius r (in Rsun):
rss = 2.5

# INPUT MAP OF br(theta,phi)
d = 180.0 / np.pi
# - read in map from SFT simulation:
br0 = np.loadtxt('test.dat')    # DEMO -- replace this with your own map
ns0 = np.size(br0, axis=0)
np0 = np.size(br0, axis=1)

# - original grid:
th0 = np.linspace(0, np.pi, ns0)
ph0 = np.linspace(0, 2 * np.pi, np0)

# - pfss grid:
dphi = 2 * np.pi / nphi
ds = 2 / ns

pc = np.linspace(0.5 * dp, 2 * np.pi - 0.5 * dp, np)
sc = np.linspace(-1 + 0.5 * ds, 1 - 0.5 * ds, ns)
thc = np.flip(np.arccos(sc), 0)

# - rotate to Carrington frame:
clon = 290.7 / d
pc1 = (pc - clon + 2 * np.pi) % (2 * np.pi)
# - interpolate input map to pfss grid:
bri = interp2d(ph0, th0, br0, kind='cubic', copy=True,
               bounds_error=False, fill_value=0)
br = np.zeros((ns, nphi))
for i in range(np):
    br[:, i] = bri(pc1[i], thc).flatten()

# - correct flux balance (otherwise PFSS solver won't work):
ipos = np.where(br > 0)
ineg = np.where(br < 0)
fluxp = np.abs(np.sum(br[ipos]))
fluxn = np.abs(np.sum(br[ineg]))
fluxmn = 0.5 * (fluxp + fluxn)
br[ineg] = br[ineg] * fluxmn / fluxn
br[ipos] = br[ipos] * fluxmn / fluxp


# PLOT INPUT MAP BEFORE AND AFTER:
# - threshold for colour scales:
bmax = 10

plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax = plt.subplot(211)
lat0 = 0.5 * np.pi - th0
pm = ax.pcolormesh(ph0*d, lat0 * d, br0, cmap='bwr')
pm.set_clim(vmin=-bmax, vmax=bmax)
cb1 = plt.colorbar(pm)
# cb1.set_clim(vmin=-bmax, vmax=bmax)
ax.set_xlabel('Carrington longitude')
ax.set_ylabel('Latitude')
ax.set_title('Original input data')

ax = plt.subplot(212)
lat = 0.5 * np.pi - thc
pm = ax.pcolormesh(pc * d, lat * d, br, cmap='bwr')
pm.set_clim(vmin=-bmax, vmax=bmax)
cb1 = plt.colorbar(pm)
# cb1.set_clim(vmin=-bmax, vmax=bmax)
ax.set_xlabel('Carrington longitude')
ax.set_ylabel('Latitude')
ax.set_title('Mapped to PFSS grid')

plt.subplots_adjust(hspace=.5)

plt.show()

# COMPUTE POTENTIAL FIELD:
# (note that output='bg' gives the output B averaged to grid points)
print('Computing PFSS...')
br = br[::-1, :]
pfss(br, nr, ns, np, rss, filename='./bPF.nc', output='bg', testQ=False)
