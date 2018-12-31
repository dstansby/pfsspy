"""
    Script for reading in map of Br(theta, phi) on the solar surface, computing a
    PFSS extrapolation, and outputting to a netcdf file.
    
    Copyright (C) Anthony R. Yeates, Durham University 29/8/17

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import sys
import numpy as n
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
from pfss import pfss

# DEFINE GRID FOR PFSS COMPUTATION
# - equally spaced in rho=log(r/Rsun), s=cos(theta) and phi.
# - specify number of grid points:
nr = 60
ns = 180
np = 360
if ((np%2)!=0):
    sys.exit('ERROR: np must be even (for polar boundary conditions)')
# - choose source surface radius r (in Rsun):
rss = 2.5

# INPUT MAP OF br(theta,phi)
d = 180.0/n.pi
# - read in map from SFT simulation:
br0 = n.loadtxt('test.dat')    # DEMO -- replace this with your own map
ns0 = n.size(br0, axis=0)
np0 = n.size(br0, axis=1)
# - original grid:
th0 = n.linspace(0, n.pi, ns0)
ph0 = n.linspace(0, 2*n.pi, np0)
# - pfss grid:
dp = 2*n.pi/np
ds = 2/ns
pc = n.linspace(0.5*dp, 2*n.pi - 0.5*dp, np)
sc = n.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns)
thc = n.flip(n.arccos(sc), 0)
# - rotate to Carrington frame:
clon = 290.7/d
pc1 = (pc - clon + 2*n.pi) % (2*n.pi)
# - interpolate input map to pfss grid:
bri = interp2d(ph0, th0, br0, kind='cubic', copy=True, bounds_error=False, fill_value=0)
br = n.zeros((ns,np))
for i in range(np):
    br[:,i] = bri(pc1[i], thc).flatten()
# - correct flux balance (otherwise PFSS solver won't work):
ipos = n.where(br > 0)
ineg = n.where(br < 0)
fluxp = n.abs(n.sum(br[ipos]))
fluxn = n.abs(n.sum(br[ineg]))
fluxmn = 0.5*(fluxp + fluxn)
br[ineg] = br[ineg]*fluxmn/fluxn
br[ipos] = br[ipos]*fluxmn/fluxp


# PLOT INPUT MAP BEFORE AND AFTER:
# - threshold for colour scales:
bmax = 10

plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax = plt.subplot(211)
lat0 = 0.5*n.pi - th0
pm = ax.pcolormesh(ph0*d, lat0*d, br0, cmap='bwr')
pm.set_clim(vmin=-bmax, vmax=bmax)
cb1 = plt.colorbar(pm)
#cb1.set_clim(vmin=-bmax, vmax=bmax)
ax.set_xlabel('Carrington longitude')
ax.set_ylabel('Latitude')
ax.set_title('Original input data')

ax = plt.subplot(212)
lat = 0.5*n.pi - thc
pm = ax.pcolormesh(pc*d, lat*d, br, cmap='bwr')
pm.set_clim(vmin=-bmax, vmax=bmax)
cb1 = plt.colorbar(pm)
#cb1.set_clim(vmin=-bmax, vmax=bmax)
ax.set_xlabel('Carrington longitude')
ax.set_ylabel('Latitude')
ax.set_title('Mapped to PFSS grid')

plt.subplots_adjust(hspace=.5)

plt.show()

# COMPUTE POTENTIAL FIELD:
# (note that output='bg' gives the output B averaged to grid points)
print('Computing PFSS...')
br = br[::-1,:]
pfss(br, nr, ns, np, rss, filename='./bPF.nc', output='bg', testQ=False)
