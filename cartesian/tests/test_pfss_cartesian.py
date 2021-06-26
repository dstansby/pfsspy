# import sys
# sys.path.append('../')
from cartesian import solve_cartesian
import sunpy.map as mp
import numpy as np 

def test_arcade():
    grid_lims=[-20,20]
    n_points=10
    gridx = np.linspace(grid_lims[0],grid_lims[1],n_points,endpoint=False)
    gridy = np.copy(gridx)
    grid_lims_z=[0,40]
    gridz = np.linspace(grid_lims_z[0],grid_lims_z[1],n_points,endpoint=False)
    dx = np.gradient(gridx)[3]
    dy = np.gradient(gridy)[3]
    dz = np.gradient(gridz)[3]
    xx,yy,zz = np.meshgrid(gridx,gridy,gridz,indexing='ij')
    k=np.pi/10
    l=k
    Bx=-(l/k)*np.cos(k*xx)*np.exp(-l*zz)
    By=-np.sqrt(1-(l/k)**2)*np.cos(k*xx)*np.exp(-l*zz)
    Bz=np.sin(k*xx)*np.exp(-l*zz)

    B_analytical=np.asarray([Bx,By,Bz]).transpose([0,3,2,1])
    #Blos = Bz, so we check our PFSS solution on this. But use the array at z=0.
    B_field = solve_cartesian(Bz[:,:,0],grid_lims_z[-1],dx,dy,dz,coordinates = None)

    abs_diff = np.mean(np.abs(B_analytical-B_field))
    print(abs_diff)

    assert abs_diff<1e-3

# #data url: https://drive.google.com/drive/folders/150sI8Bu3JFmoViN3cpZiocbtWBweeDhj?usp=sharing
# BASE = "/home/uvishal/Desktop/IUCAA/pfsspy/test_data/"
# path_data = f"{BASE}raster_0_hmi.fits"
# hmi_blos=mp.Map(path_data)
# print(hmi_blos.data.shape)
# dz = 0.1
# Z_max = 10
# input_obj = pieck.Input(hmi_blos,Z_max,dz)
# extrapolated_obj = panzer.fft_pot_open(input_obj)
# extrapolated_obj.save(f"{BASE}extrapolation.fits")
