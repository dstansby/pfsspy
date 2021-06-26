import numpy as np 
import sunpy 

class CoordObj:
    def __init__(self,coord):
        self.lat = coord[0]
        self.lon = coord[1]
class ExtentObj:
    def __init__(self,coordinates,shape):
        self.bottom_left_coord = CoordObj(coordinates[0])
        self.top_right_coord = CoordObj(coordinates[1])
        self.size_lat = shape[0]
        self.size_lon = shape[1]

def get_units(mapobj):
    return mapobj.meta['cunit1'],mapobj.meta['cunit1']
def Mm_to_unit(unit):
    """Convert Megameters to degree or arcsec.
        1arcsec = 725 km on the Sun. 1 degree= 3600 arcsec

    Args:
        unit (string): "degree" or "arcsec"

    Returns:
        (float): Conversion factor FROM Mm TO unit
    """
    if unit=='deg':
        return 1000.0/(725*3600)
    elif unit=='arcsec':
        return 1000.0/(725)
    else:
        raise ValueError("Units should be degree or arcsec.")

def calc_mu(x,y):
    """Calculates the heliocentric cosine angle given a map object.
      The cosine value, called "mu" is defined as :
        mu = cos(theta), where
        
        sin(theta) = sqrt(x**2+y**2)/R, 
        
        where x,y are longitude and latitude resp in arcsec, and R = 1919/2 = 959.6 arcsec.
        
        By simple trignometry, we get:
        cos(theta) = sqrt(R**2-x**2-y**2)/R
        
        Reference: arXiv:1307.2410v1 
    Args:
        x (float): latitude value in arcsec
        y (float): longitude value in arcsec
    Returns:
        mu = sqrt(R**2-x**2-y**2)/R
    """
    R=1919.0/2 #In arcsec.
    mu = np.sqrt(R**2-np.square(x)-np.square(y))/R
    return mu 

def calc_cosine_heliocentric_angle(map_ob):
    """Calculates the heliocentric cosine angle given a map object or an ExtentObj.

        Only use this function if the map object is known. If you have an array of 
        latitudes or longitudes, use the function calc_mu.
    Args:
        map_ob (sunpy.map.GenericMap,ExtentObj ): a sunpy map object or ExtentObj
    """
    if isinstance(map_ob,sunpy.map.GenericMap):
        lat_r = map_ob.fits_header['CRVAL1']+map_ob.fits_header['CDELT1']*map_ob.fits_header['NAXIS2']/2
        lat_l = map_ob.fits_header['CRVAL1']-map_ob.fits_header['CDELT1']*map_ob.fits_header['NAXIS2']/2
        lon_r = map_ob.fits_header['CRVAL2']+map_ob.fits_header['CDELT2']*map_ob.fits_header['NAXIS1']/2
        lon_l = map_ob.fits_header['CRVAL2']-map_ob.fits_header['CDELT2']*map_ob.fits_header['NAXIS1']/2
        latitude = np.linspace(lat_l,lat_r,map_ob.fits_header['NAXIS2']).reshape([-1,1])/0.000277778
        longitude = np.linspace(lon_l,lon_r,map_ob.fits_header['NAXIS1']).reshape([1,-1])/0.000277778
        mu = calc_mu(latitude,longitude)
    elif isinstance(map_ob,ExtentObj):
        latitude = np.linspace(map_ob.bottom_left_coord.lat,map_ob.top_right_coord.lat,map_ob.size_lat).reshape([-1,1])
        longitude = np.linspace(map_ob.bottom_left_coord.lon,map_ob.top_right_coord.lon,map_ob.size_lon).reshape([1,-1])
        mu = calc_mu(latitude,longitude)
    else:
        raise TypeError("Please use either a sunpy map or an ExtenObj instance.")
    return mu 

