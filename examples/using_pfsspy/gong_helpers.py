"""
GONG helper functions
=====================
"""
from sunpy.data import manager


@manager.require('gong_map',
                 'https://gong2.nso.edu/oQR/zqs/201903/mrzqs190310/mrzqs190310t0014c2215_333.fits.gz',
                 '712c7543fa1d964d03e73523ec2256766348db3a2dd5a8406e2f0d711b666cb6')
def get_gong_map():
    """
    Automatically download and unzip a sample GONG synoptic map.
    """
    return manager.get('gong_map')
