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
        rolled to align Carrington Longitude 0 with LH edge of pixel 1

    """

    ## Check Map is Synoptic Type
    assert synmap.meta['CTYPE1'].startswith("CRLN"), \
           "Not synoptic map, FITS header CTYPE1 is not CRLN type"

    data,header = synmap.data,synmap.meta
    rolled_header = copy.copy(header)


    # Determine if CR0 is initially aligned with a pixel edge, if not 
    # interpolate to make this true. Note pixel edges are located at 
    # half integers in the FITS standard 
    
    # Pull out header values which need to be mutable :
    crval1 = header['CRVAL1']

    # First calculate initial position of CR0 in pixel units
    CR0PIX_init = header['CRPIX1'] - crval1/header['CDELT1']
    # ensure  pixel value of CR0 is in the range [0.5, NAXIS1 + 0.5]
    CR0PIX_init = (CR0PIX_init - 0.5) % header['NAXIS1'] + 0.5

    # If CR0 not on pixel edge, perform interpolation of data to a grid where
    # CR0 is enforced to be on a pixel edge
    if not np.isclose(CR0PIX_init % 1.0, 0.5) :
        ## Coordinates are position of pixel **centers**
        lonpx = np.linspace(0,header['NAXIS1']+1,header['NAXIS1']+2)
        latpx = np.linspace(1,header['NAXIS2'],header['NAXIS2'])
        ## To account for the periodic boundary, we agument the original
        ## data with the left hand most column appended to the right, 
        ## and vice versa
        data_aug = np.append(data[:,-1][:,None],
                            np.append(data,data[:,0][:,None],axis=1),
                            axis=1)
        interpolator = interpolate.RectBivariateSpline(lonpx,latpx,data_aug.T)

        ## Interpolate to new lonpx such that CR0 is located at 
        ## the **next edge ABOVE it's initial value**
        shift = (CR0PIX_init + 0.5) % 1.0
        lonpx_interp = np.linspace(1-shift,header['NAXIS1']-shift,header['NAXIS1'])
        data_ = interpolator(lonpx_interp,latpx).T

        # update crval1 to reflect interpolation
        crval1 -= (1-shift)

        # Assert CR0 is now on pixel edge :
        CR0PIX = header['CRPIX1'] - crval1/header['CDELT1']
        CR0PIX = (CR0PIX - 0.5) % header['NAXIS1'] + 0.5
        assert np.isclose(CR0PIX % 1.0, 0.5), f"CR0 not on pixel edge (pixel coord = {CR0PIX})"

    else : data_ = data.copy()

    # The roll we need is such that 0 longitude is at pixel position 0.5
    # which is the LH edge of pixel 1.0. 
    # 
    # Our transformation should maintain naxis1 and cdelt1
    #
    # Calculate the new value of crval1 such that carrington 0 is located
    # at the LH edge of pixel 1. We need to add 0.5*cdelt1 to crpix1*cdelt1
    # to achieve this. crval1 = crpix1*cdelt1 would place 0 at the center of
    # pixel 1.
    #
    crval1_rolled = (header["CRPIX1"]-0.5)*header['CDELT1']
    
    # Calculate longitude roll needed
    roll_lng = crval1_rolled - crval1

    # Convert to pixel units
    roll_px = roll_lng/header['CDELT1']

    # Roll the data array by int(roll_px) to get 
    rolled_data = np.roll(data_,int(roll_px),axis=1)

    ## Our transformation changes only "CRVAL1"
    #  "CRPIX1" and "NAXIS1" are unchanged    
    rolled_header['CRVAL1'] = crval1_rolled

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

    # Get spacing of cell centers of original map
    dlo,dla = synmap.meta['cdelt1'],synmap.meta['cdelt2'] #delta_lon,delta_lat

    ## Original Coordinates of *Cell Centers*  # Note we ignore any roll of the map from CRLN=0
    lons = np.linspace(0+dlo/2,360-dlo/2,synmap.meta['naxis1'])
    lats = np.linspace(-90+dla/2,90-dla/2,synmap.meta['naxis2'])

    ## Get data and header
    data,header = synmap.data,synmap.meta

    # Calculate cell spacing of y-axis of desired CEA map
    dsl = 2.0/synmap.meta['naxis2'] # delta_sinlat

    ## Array of lats to interpolate to for CEA Binning
    sinlats = np.linspace(-1+dsl/2,1-dsl/2,synmap.meta['naxis2'])
    lats_interp = np.degrees(np.arcsin(sinlats))

    ## Create interpolator function 
    synmap_interpolator = interpolate.RectBivariateSpline(lons,lats,data.T)

    ## Perform interpolation
    data_interp = synmap_interpolator(lons,lats_interp).T

    ## Adjust Header
    header_interp = copy.copy(header)
    header_interp['ctype1'] = "CRLN-CEA"
    header_interp['ctype2'] = "CRLT-CEA"
    header_interp['cdelt2'] = np.diff(sinlats)[0]*180/np.pi 

    # CRVAL2 only needs to be adjusted if it is not initially 0.0
    if header_interp['crval2'] != 0.0 : 
        header_interp['crval2'] = header_interp['crpix2']*header_interp['cdelt2'] - sinlats[0]

    # Hack to stop ADAPT GenericMap subclass overriding ctypes - this
    # makes the output a raw GenericMap (no subclass). 
    header_interp['model'] = header.get('model',"")+"-CEA"
    
    return sunpy.map.Map(data_interp,header_interp)

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



