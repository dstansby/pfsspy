"""
ADAPT helper functions
======================
"""
from sunpy.data import manager


@manager.require('adapt_map',
                 'https://gong.nso.edu/adapt/maps/gong/2020/adapt40311_03k012_202001010000_i00005600n1.fts.gz',
                 'fd8f3a23059b2118d0097c448df3ccbdb5388a0031536dd3a6f61fa0e08a9bb5')
def example_adapt_map():
    """
    Automatically download and unzip a sample GONG synoptic map.
    """
    return manager.get('adapt_map')
