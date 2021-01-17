import pfss_cartesian_obj as pieck
"""
    Attack on Titan trivia: Pieck-san is the "Cart" titan.
    Get it? Cartesian --> Cart?

"""
import numpy as np 

def fft_pot_decompose_open(Bz):
    """
        Just computed the Fourier transform of Bz at z=0.
    """
    coeffs = np.fft.fft2(Bz)
    coeffs[0,0] = 0.0 #Can we even do this? Shouldn't our FOV be actually flux balanced?
    return coeffs

def fft_pot_field_plane_open(coeff,kx,ky,k,z):
    bx=np.real(np.fft.ifft2((-1j*kx/k)*coeff*np.exp(-k*z)))
    by=np.real(np.fft.ifft2((-1j*ky/k)*coeff*np.exp(-k*z)))
    bz=np.real(np.fft.ifft2(coeff*np.exp(-k*z)))
    return np.array([bx,by,bz])

def fft_pot_open(input):
    """
        Bz: sunpy.map object of Bz at z=0.
        Z_max: max Z value in Mm. 1arcsec = 725 km on the Sun. 1 degree= 3600 arcsec
        dz: Grid spacing in Z, in Mm.
    """
    Bz = input.map
    latitude = np.cos(np.linspace(Bz.bottom_left_coord.data.lat,Bz.top_right_coord.data.lat,Bz.data.shape[0]))
    print("Removing monopole term...")
    coeff = fft_pot_decompose_open((Bz.data-np.mean(Bz.data))/latitude.reshape([-1,1]))
    fits_keys = Bz.fits_header
    Z_range = input.gridz
    
    kx=2*np.pi*np.fft.fftfreq(coeff.shape[0],d=fits_keys['CDELT2']).reshape([-1,1])
    ky=2*np.pi*np.fft.fftfreq(coeff.shape[1],d=fits_keys['CDELT1']).reshape([1,-1])
    k=np.sqrt(kx**2+ky**2)
    k[0,0]=1.0
    B_field=np.asarray([fft_pot_field_plane_open(coeff,kx,ky,k,z_i) for z_i in Z_range]).transpose([1,0,3,2])
    return pieck.Output(B_field,input)