import sunpy.map,sunpy.io
import numpy as np
import os, warnings,copy
from scipy import interpolate
from astropy import units as u

def roll_to_CRLN(synmap) :
    """ Roll sunpy.map.Map synoptic map so LH edge is 0 Carr Lon.

    Parameters
    ----------
    synmap : `sunpy.map.GenericMap` 
        a `sunpy.map.GenericMap` instance containing a synoptic map
    
    Returns
    -------
    synmap_rolled : `sunpy.map.GenericMap` 
        a `sunpy.map.GenericMap` instance with it's data appropriately 
        rolled to align Carrington Longitude 0 with LH edge

    """
    data,header = synmap.data,synmap.meta
    rolled_header = copy.copy(header)

    # Identify reference pixel and calculate what it's longitude
    # should be after the roll
    rolled_crval1 = (header['CRPIX1']-0.5)*header['CDELT1']

    # Calculate longitude roll needed
    roll_lng = rolled_crval1 - header['CRVAL1']

    # Convert to integer pixels
    roll_px = int(roll_lng*header['CDELT1'])
    
    rolled_data = np.roll(data,int(roll_px),axis=1)
    rolled_header['CRVAL1'] = (header['CRPIX1']-0.5)*header['CDELT1']

    return sunpy.map.Map(rolled_data,rolled_header)


def interp_CAR2CEA(synmap) :
    """ Interpolate Plate-Carree synoptic maps to CEA

    Description
    -----------
    Plate-Carree Projected maps [CAR] are binned uniformly in lat
    and lon. Cylindrical equal area [CEA] maps have their latitude
    binned uniformly in sin(latitude). This is the required 
    projection for `sunpy.map.Map` arguments passed to 
    `pfsspy.Input`. 

    This function performs this interpolation on a given 
    `sunpy.map.Map` object with initial `.meta['ctype2')`
    

    Parameters
    ----------
    synmap : `sunpy.map.Map` 
        a `sunpy.map.Map` instance containing a synoptic map
    
    Returns
    -------
    synmap_CEA : `sunpy.map.Map` 
        a `sunpy.map.Map` instance with it's data reprojected
        to CEA coordinates.

    """

    ## Check Input Map Header implies Plate-Carree projected synoptic map
    lat_type = synmap.meta['ctype2']
    if lat_type == "CRLT-CAR" : ""
    elif lat_type == "CRLT-CEA" :
        warnings.warn("Input map already in CEA projection, returning unchanged")
        return synmap
    else : 
        raise TypeError(f"Header 'ctype2' {lat_type} not recognized, "\
                        +"map may not be synoptic map")

    ## Original Data Axes # Note we ignore any roll of the map from CRLN=0
    lons = np.linspace(0,360,synmap.meta['naxis1'])
    lats = np.linspace(-90,90,synmap.meta['naxis2'])

    ## Get data and header
    data,header = synmap.data,synmap.meta

    ## Array of lats to interpolate to for CEA Binning
    sinlats = np.linspace(-1,1,synmap.meta['naxis2'])
    lats_interp = np.degrees(np.arcsin(sinlats))

    ## Create interpolator function 
    synmap_interpolator = interpolate.RectBivariateSpline(lons,lats,data.T)

    ## Perform interpolation
    data_interpolated = synmap_interpolator(lons,lats_interp).T

    ## Adjust Header
    header_interpolated = copy.copy(header)
    header_interpolated['ctype1'] = "CRLN-CEA"
    header_interpolated['ctype2'] = "CRLT-CEA"
    header_interpolated['cdelt2'] = np.diff(sinlats)[0]*180/np.pi 
    # Hack to stop ADAPT GenericMap subclass overriding ctypes - this
    # makes the output a raw GenericMap (no subclass). 
    header_interpolated['model'] = header.get('model',"")+"-CEA"
    
    return sunpy.map.Map(data_interpolated,header_interpolated)

def load_adapt(adapt_path) :
    """ Parse adapt .fts file as a `sunpy.map.MapSequence`

    Description
    -----------
    ADAPT magnetograms contain 12 realizations and their data
    attribute consists of a 3D data cube where each slice is
    the data corresponding to a separate realization of the
    magnetogram. This function loads the raw fits file and 
    parses it to a `sunpy.map.MapSequence` object containing
    a `sunpy.map.Map` instance for each realization.
    
    Parameters
    ----------
    adapt_path : `str` 
        filepath corresponding to an ADAPT .fts file
    
    Returns
    -------
    adaptMapSequence : `sunpy.map.MapSequence` 
        a `sunpy.map.MapSequence` instance containing a
        `sunpy.map.Map` instance for each ADAPT realization.

    """
    adapt_fits = sunpy.io.fits.read(adapt_path)
    assert adapt_fits[0].header.get('MODEL') == 'ADAPT', \
        f"{os.path.basename(adapt_path)} header['MODEL'] is not 'ADAPT' "
    data_header_pairs = [
                        (map_slice,adapt_fits[0].header) 
                         for map_slice in adapt_fits[0].data
                        ]
    adaptMapSequence = sunpy.map.Map(data_header_pairs,sequence=True) 
    return adaptMapSequence



