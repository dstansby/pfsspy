import os

import astropy.units as u
import sunpy.map

import pfsspy.map
from .example_maps import gong_map, adapt_map

import pytest


def test_gong_source(gong_map):
    m = sunpy.map.Map(gong_map)
    assert isinstance(m, pfsspy.map.GongSynopticMap)
    # Check round-trip is robust against sunpy changes to the meta
    m = sunpy.map.Map(m.data, m.meta)
    assert m.date.isot == '2020-09-01T13:04:00.000'
    # Construct a WCS to check no warnings are thrown
    m.wcs
    # Check observer coordinate is populated
    observer = m.coordinate_frame.observer
    assert observer.obstime.isot == m.date.isot
    assert observer.lon == 0 * u.deg
    assert u.allclose(observer.lat, 7.20584924 * u.deg)
    assert u.allclose(observer.radius, 1.50953137e+11 * u.m)


def test_adapt_map(adapt_map):
    import sunpy.io
    adapt_fits = sunpy.io.fits.read(adapt_map)
    for map_slice in adapt_fits[0].data:
        m = sunpy.map.Map((map_slice, adapt_fits[0].header))
        assert isinstance(m, pfsspy.map.ADAPTMap)
