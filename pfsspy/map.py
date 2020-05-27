"""
Custom `sunpy.map.GenericMap` sub-classes for different magnetogram sources.
"""
import numpy as np
import sunpy.map


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
        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an GONG map."""
        return (str(header.get('TELESCOP', '')).endswith('GONG') and
                str(header.get('CTYPE1', '').startswith('CRLN')))


class ADAPTMap(sunpy.map.GenericMap):
    def __init__(self, data, header, **kwargs):
        header['date-obs'] = header['maptime']
        if not ((header['cunit2'] == 'deg') &
                (header['naxis2']*header['cdelt2'] == 180)) :
            raise AssertionError("Latitude metadata doesn't add to 180deg")
        header['ctype1'] = 'CRLN-CAR'
        header['ctype2'] = 'CRLT-CAR'
        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an ADAPT map."""
        return header.get('model') == 'ADAPT'
