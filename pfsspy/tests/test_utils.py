import os
import numpy as np

import sunpy.map
import pytest

import pfsspy
from pfsspy import utils
from .example_maps import adapt_map


def test_load_adapt(adapt_map):
    adaptMapSequence = utils.load_adapt(adapt_map)
    assert isinstance(adaptMapSequence, sunpy.map.MapSequence)
    for map_ in adaptMapSequence:
        assert map_.meta['model'] == "ADAPT"


def test_header_generation():
    dtime = None
    ntheta = 180
    nphi = 360
    shape = [nphi, ntheta]
    header = pfsspy.utils.carr_cea_wcs_header(dtime, shape)
    assert header['LONPOLE'] == 0
    assert header['CTYPE1'] == 'CRLN-CEA'
    assert header['CTYPE2'] == 'CRLT-CEA'
    assert header['PV1_1'] == 1
    assert header['PV2_1'] == 1
    assert header['CDELT1'] == 360 / nphi
    np.testing.assert_almost_equal(
        header['CDELT2'], (180 / np.pi) * (2 / ntheta))

    # Extra + 1s are for FITS counting from 1 indexing
    assert header['CRPIX1'] == (nphi / 2) + 0.5 + 1
    assert header['CRPIX2'] == (ntheta / 2) + 0.5 + 1
    assert header['CRVAL1'] == 0
    assert header['CRVAL2'] == 0
    assert header['CUNIT1'] == 'deg'
    assert header['CUNIT2'] == 'deg'


@pytest.mark.parametrize('error', [True, False])
def test_validation(dipole_map, error):
    assert utils.is_cea_map(dipole_map, error)
    assert utils.is_full_sun_synoptic_map(dipole_map, error)
