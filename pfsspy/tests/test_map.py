import os

import sunpy.map
import pytest

import pfsspy.map
from .example_maps import gong_map, adapt_map


def test_gong_source(gong_map):
    m = sunpy.map.Map(gong_map)
    assert isinstance(m, pfsspy.map.GongSynopticMap)
    # Check round-trip is robust against sunpy changes to the meta
    m = sunpy.map.Map(m.data, m.meta)
    assert m.date.isot == '2019-03-10T00:14:00.000'
    # Construct a WCS to check no warnings are thrown
    m.wcs


def test_adapt_map(adapt_map):
    import sunpy.io
    adapt_fits = sunpy.io.fits.read(adapt_map)
    for map_slice in adapt_fits[0].data:
        m = sunpy.map.Map((map_slice, adapt_fits[0].header))
        assert isinstance(m, pfsspy.map.ADAPTMap)
