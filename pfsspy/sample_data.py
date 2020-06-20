"""
Functions to get sample data. The data is fetched and cached automatically
using `sunpy.data.manager`.
"""

from sunpy.data import manager


@manager.require('gong_map',
                 'https://gong2.nso.edu/oQR/zqs/201903/mrzqs190310/'
                 'mrzqs190310t0014c2215_333.fits.gz',
                 '712c7543fa1d964d03e73523ec2256766348d'
                 'b3a2dd5a8406e2f0d711b666cb6')
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
