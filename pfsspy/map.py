"""
Custom `sunpy.map.GenericMap` sub-classes for different magnetogram sources.
"""
import numpy as np
import sunpy.map
from astropy.time import Time


class GongSynopticMap(sunpy.map.GenericMap):
    def __init__(self, data, header, **kwargs):
        if 'KEYCOMMENTS' in header:
            if 'deg' in header['KEYCOMMENTS']['CDELT1']:
                header['CUNIT1'] = 'deg'
            if header['KEYCOMMENTS']['CDELT2'] == 'Sine-lat step':
                header['CUNIT2'] = 'deg'
                # Instead of the spacing in sin(lat), this should be 180/pi times
                # that value (see Thompson 2005)
                header['CDELT2'] = 180 / np.pi * header['CDELT2']
                header['KEYCOMMENTS']['CDELT2'] = '180 / pi * sine-lat step'
            if 'time-obs' in header:
                header['date-obs'] = (header['date-obs'] + ' ' +
                                      header.pop('time-obs'))
                header['date-obs'] = Time(header['date-obs']).isot
        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an GONG map."""
        return (str(header.get('TELESCOP', '')).endswith('GONG') and
                str(header.get('CTYPE1', '').startswith('CRLN')))


class ADAPTMap(sunpy.map.GenericMap):
    def __init__(self, data, header, **kwargs):
        if 'date-obs' not in header:
            header['date-obs'] = header['maptime']
        # Fix CTYPE
        if header['ctype1'] == 'Long':
            header['ctype1'] = 'CRLN-CAR'
        if header['ctype2'] == 'Lat':
            header['ctype2'] = 'CRLT-CAR'

        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an ADAPT map."""
        return header.get('model') == 'ADAPT'
