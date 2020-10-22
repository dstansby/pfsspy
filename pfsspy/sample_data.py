"""
Functions to get sample data. The data is fetched and cached automatically
using `sunpy.data.manager`.
"""

from sunpy.data import manager


@manager.require('gong_map',
                 'https://gong2.nso.edu/oQR/zqs/202009/mrzqs200901/'
                 'mrzqs200901t1304c2234_022.fits.gz',
                 'aad927d8f617f32b72255b862c4910f13640fc7ca13edf982'
                 '88cd0735a2db6a0')
def get_gong_map():
    """
    Automatically download and unzip a sample GONG synoptic map.
    """
    return manager.get('gong_map')


@manager.require('adapt_map',
                 'https://gong.nso.edu/adapt/maps/gong/2020/'
                 'adapt40311_03k012_202001010000_i00005600n1.fts.gz',
                 'fd8f3a23059b2118d0097c448df3cc'
                 'bdb5388a0031536dd3a6f61fa0e08a9bb5')
def get_adapt_map():
    return manager.get('adapt_map')
