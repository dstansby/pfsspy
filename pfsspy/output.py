"""
Routines for writing netcdf files.
"""

from scipy.io import netcdf
import numpy as n


def bc(filename, r, th, ph, rc, thc, phc, br, bs, bp):
    """
    Magnetic field components on cell faces, including ghost cells.
    """

    nr = n.size(r) - 1
    ns = n.size(th) - 1
    np = n.size(ph) - 1

    fid = netcdf.netcdf_file(filename, 'w')
    fid.createDimension('rc', nr + 2)
    fid.createDimension('r', nr + 1)
    fid.createDimension('thc', ns + 2)
    fid.createDimension('th', ns + 1)
    fid.createDimension('phc', np + 2)
    fid.createDimension('ph', np + 1)
    vid = fid.createVariable('r', 'd', ('r',))
    vid[:] = r
    vid = fid.createVariable('th', 'd', ('th',))
    vid[:] = th
    vid = fid.createVariable('ph', 'd', ('ph',))
    vid[:] = ph
    vid = fid.createVariable('rc', 'd', ('rc',))
    vid[:] = rc
    vid = fid.createVariable('thc', 'd', ('thc',))
    vid[:] = thc
    vid = fid.createVariable('phc', 'd', ('phc',))
    vid[:] = phc
    vid = fid.createVariable('br', 'd', ('phc', 'thc', 'r'))
    vid[:] = br
    vid = fid.createVariable('bth', 'd', ('phc', 'th', 'rc'))
    vid[:] = -bs
    vid = fid.createVariable('bph', 'd', ('ph', 'thc', 'rc'))
    vid[:] = bp
    fid.close()
    print('Wrote B on faces to file ' + filename)
