import pfss_cartesian_obj as pieck
from utils import ExtentObj,calc_cosine_heliocentric_angle
import warnings
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

def solve_cartesian(Blos,Zmax,dx,dy,dz,coordinates = None):
    """Generate a 3D PFSS extrapolated cube given line of sight B field.

    Args:
        Blos (2-D array): A 2-D array of B LoS measurement (signed)
        Zmax (float): max Z value in Mm. 1arcsec = 725 km on the Sun. 1 degree= 3600 arcsec
        dx (float): Grid spacing in x, in Mm
        dy (float): Grid spacing in y, in Mm
        dz (float): Grid spacing in Z, in Mm
        coordinates (ndarray): A 2x2 matrix of the form [bottom_left(lat,lon),top_right(lat,lon)], indicating the bottom left 
                                and top right coordinates of the map in two rows.

    """
    if coordinates is None:
        mu = 1.0 
    else:
        coordinates = ExtentObj(coordinates,Blos.shape)
        mu = calc_cosine_heliocentric_angle(coordinates)
    monopole = np.mean(Blos)
    #Ad hoc check for monopole term.
    if np.abs(monopole)>3:
        warnings.warn("Monopole term > 3 Gauss. The code removed monopole term before compuation of PFSS. You have been warned....")
    
    coeff = fft_pot_decompose_open((Blos-monopole)/mu)
    Z_range = np.arange(0,Zmax,dz)
    kx=2*np.pi*np.fft.fftfreq(coeff.shape[0],d=dx).reshape([-1,1])
    ky=2*np.pi*np.fft.fftfreq(coeff.shape[1],d=dy).reshape([1,-1])
    k=np.sqrt(kx**2+ky**2)
    k[0,0]=1.0
    # TODO: Possibly parallelize this loop. 
    B_field=np.asarray([fft_pot_field_plane_open(coeff,kx,ky,k,z_i) for z_i in Z_range]).transpose([1,0,3,2])
    return B_field

def fft_pot_open(input):
    """
        Bz: sunpy.map object of Bz at z=0.
        Z_max: max Z value in Mm. 1arcsec = 725 km on the Sun. 1 degree= 3600 arcsec
        dz: Grid spacing in Z, in Mm.
    """
    Blos = input.map
    mu = calc_cosine_heliocentric_angle(Blos)
    #latitude = np.cos(np.linspace(Bz.bottom_left_coord.data.lat,Bz.top_right_coord.data.lat,Bz.data.shape[0]))
    print("Removing monopole term...")
    coeff = fft_pot_decompose_open((Blos.data-np.mean(Blos.data))/mu)
    fits_keys = Blos.fits_header
    Z_range = input.gridz
    
    kx=2*np.pi*np.fft.fftfreq(coeff.shape[0],d=fits_keys['CDELT2']).reshape([-1,1])
    ky=2*np.pi*np.fft.fftfreq(coeff.shape[1],d=fits_keys['CDELT1']).reshape([1,-1])
    k=np.sqrt(kx**2+ky**2)
    k[0,0]=1.0
    B_field=np.asarray([fft_pot_field_plane_open(coeff,kx,ky,k,z_i) for z_i in Z_range]).transpose([1,0,3,2])
    return pieck.Output(B_field,input)