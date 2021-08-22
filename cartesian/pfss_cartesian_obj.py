import warnings

import numpy as np
import sunpy.map
from astropy.io import fits
import utils

class Input:
    r"""
    Input to PFSS modelling.

    This code uses a Cartesian geometry (x,y,z) for performing the extrapolation. 
    Thus, make sure the grid defined is cartesian, and that solutions obtained are 
    valid at all. 
    Parameters
    ----------
    bz : sunpy.map.GenericMap
        Boundary condition of LOS magnetic field at the inner surface.

    Z_max : float
        Max height in vertical direction to which extrapolation is performed, in Mm (Mega meters)

    dz : float
        Grid spacing. 
    """
    def __init__(self, bz, Z_max, dz):
        if not isinstance(bz, sunpy.map.GenericMap):
            raise ValueError('br must be a SunPy Map')
        if np.any(~np.isfinite(bz.data)):
            raise ValueError('At least one value in the input is NaN or '
                             'non-finite. The input must consist solely of '
                             'finite values.')
        if np.mean(bz.data) > 1e-10:
            warnings.warn('Input data has a non-zero mean. '
                          'pfsspy will ignore this non-zero monopole term '
                          'when calculating the PFSS solution.')
        
        #Check if the x and y coordinates are in degrees. 
        if len(bz.data.shape)!=2:
            raise ValueError('The data should be a 2-D map of LOS B field.')
        
        unit1,unit2 = utils.get_units(bz)
        assert unit1==unit2, "Units of x and y coordinates must be the same"
        
        self.unit = unit1
        self.factor = utils.Mm_to_unit(unit1)

        self._map_in = bz
        self.dtime = bz.date
        self.bz = bz.data
        self.gridz = np.arange(0,Z_max,dz)*self.factor 
        self.dz = dz 
        self.Z_max = Z_max

        # Force some nice defaults
        self._map_in.plot_settings['cmap'] = 'PuOr'
        lim = np.nanmax(np.abs(self._map_in.data))
        self._map_in.plot_settings['vmin'] = -lim
        self._map_in.plot_settings['vmax'] = lim

    @property
    def map(self):
        """
        :class:`sunpy.map.GenericMap` representation of the input.
        """
        return self._map_in

class Output:
    '''
    Output of PFSS modelling.

    Parameters
    ----------
    4-D array of floats: array([Bx,By,Bz]) over the computed volume. 

    input_obj : Input object instance. 

    Notes
    -----
    Instances of this class are intended to be created by `pfsspy.pfss`, and
    not by users.
    '''
    def __init__(self,B_field,input_obj):
        """Bx, By, Bz arranged as (z,y,x), to ease fits saving.
           FITS reverses the axis of array from (x,y,z) --> (z,y,x) while reading. 
        """
        self.B_field = B_field
        self.input_map = input_obj.map
        self.grid_z = input_obj.gridz 
        self.dz = input_obj.dz
        self.unit = input_obj.unit 
        
    @property
    def Bx(self):
        return self.B_field[0]
    @property
    def By(self):
        return self.B_field[1]
    @property
    def Bz(self):
        return self.B_field[2]

    def Add_keys_extrapolation(self,b_hdu):
        b_hdu.header['CUNIT3'] = b_hdu.header['CUNIT1']
        b_hdu.header['CTYPE3'] = 'SOLZ'
        b_hdu.header['CDELT3'] =  self.dz #dz
        b_hdu.header['CRVAL3'] = np.mean(self.grid_z) #mean height
        b_hdu.header['CRPIX3'] = np.argmin(self.grid_z-np.mean(self.grid_z))#mean height location in array
        b_hdu.header['COMMENT'] = "Use the metadata from HMI BLOS corresponding to this extrapolation for the spatial keywords. "
        return b_hdu 
    def Header_extrapolation(self): 
        ref_head = self.input_map.wcs.to_header()
        bx_hdu = fits.PrimaryHDU(self.B_field[0,:,:,:],ref_head)
        bx_hdu.header['EXTNAME']="Bx"
        by_hdu = fits.ImageHDU(self.B_field[1,:,:,:],ref_head,name="By")
        bz_hdu = fits.ImageHDU(self.B_field[2,:,:,:],ref_head,name="Bz")
        
        bx_hdu=self.Add_keys_extrapolation(bx_hdu)
        by_hdu=self.Add_keys_extrapolation(by_hdu)
        bz_hdu=self.Add_keys_extrapolation(bz_hdu)
        self.pfss_hdulist=fits.HDUList([bx_hdu,by_hdu,bz_hdu])

    def save(self,path,**kwargs):
        """Save the PFSS extrapolation.
           
           Saves the extrapolation as a FITS file which can be loaded using atropy fits opener. 

        Args:
            path (str): path to save the extrapolation
            **kwargs: Arguments to fits saver. 
        """
        self.Header_extrapolation()
        self.pfss_hdulist.writeto(path,**kwargs)
