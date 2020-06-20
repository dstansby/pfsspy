from astropy.time import Time
import numpy as np
from sunpy.map import Map
from sunpy.data import manager
import pytest

import pfsspy


@pytest.fixture
def zero_map():
    # Test a completely zero input
    ns = 30
    nphi = 20
    nr = 10
    rss = 2.5
    br = np.zeros((nphi, ns))
    header = pfsspy.utils.carr_cea_wcs_header(Time('1992-12-21'), br.shape)
    input_map = Map((br.T, header))

    input = pfsspy.Input(input_map, nr, rss)
    output = pfsspy.pfss(input)
    return input, output


@pytest.fixture
def dipole_map():
    ntheta = 30
    nphi = 20

    phi = np.linspace(0, 2 * np.pi, nphi)
    theta = np.linspace(-np.pi / 2, np.pi / 2, ntheta)
    theta, phi = np.meshgrid(theta, phi)

    def dipole_Br(r, theta):
        return 2 * np.sin(theta) / r**3

    br = dipole_Br(1, theta)
    header = pfsspy.utils.carr_cea_wcs_header(Time('1992-12-21'), br.shape)
    return Map((br.T, header))


@pytest.fixture
def dipole_result(dipole_map):
    nr = 10
    rss = 2.5

    input = pfsspy.Input(dipole_map, nr, rss)
    output = pfsspy.pfss(input)
    return input, output


@pytest.fixture
@manager.require('gong_map',
                 'https://gong2.nso.edu/oQR/zqs/201903/mrzqs190310/'
                 'mrzqs190310t0014c2215_333.fits.gz',
                 '712c7543fa1d964d03e73523ec225676'
                 '6348db3a2dd5a8406e2f0d711b666cb6')
def gong_map():
    """
    Automatically download and unzip a sample GONG synoptic map.
    """
    return manager.get('gong_map')


@pytest.fixture
@manager.require('adapt_map',
                 'https://gong.nso.edu/adapt/maps/gong/2020/'
                 'adapt40311_03k012_202001010000_i00005600n1.fts.gz',
                 'fd8f3a23059b2118d0097c448df3ccbd'
                 'b5388a0031536dd3a6f61fa0e08a9bb5')
def adapt_map():
    """
    Automatically download and unzip a sample GONG synoptic map.
    """
    return manager.get('adapt_map')
