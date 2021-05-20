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
def calc_cosine_heliocentric_angle(map_ob):
    """Calculates the heliocentric cosine angle given a map object.
      The cosine value, called "mu" is defined as :
        mu = cos(theta), where
        
        sin(theta) = sqrt(x**2+y**2)/R, 
        
        where x,y are longitude and latitude resp in arcsec, and R = 1919/2 = 959.6 arcsec.
        
        By simple trignometry, we get:
        cos(theta) = sqrt(R**2-x**2-y**2)/R
        
        Reference: arXiv:1307.2410v1 

    Args:
        map_ob (sunpy.map.GenericMap): a sunpy map object
    """
    R=1919.0/2
    if isinstance(map_ob,sunpy.map):
        latitude = np.linspace(map_ob.bottom_left_coord.data.lat,map_ob.top_right_coord.data.lat,map_ob.data.shape[0]).reshape([-1,1]).value
        longitude = np.linspace(map_ob.bottom_left_coord.data.lon,map_ob.top_right_coord.data.lon,map_ob.data.shape[1]).reshape([1,-1]).value
    elif isinstance(map_ob,ExtentObj):
        latitude = np.linspace(map_ob.bottom_left_coord.lat,map_ob.top_right_coord.lat,map_ob.size_lat).reshape([-1,1])
        longitude = np.linspace(map_ob.bottom_left_coord.lon,map_ob.top_right_coord.lon,map_ob.size_lon).reshape([1,-1])
    else:
        raise TypeError("Please use either a sunpy map or an ExtenObj instance.")
    mu = np.sqrt(R**2-np.square(latitude)-np.square(longitude))/R
    return mu 
