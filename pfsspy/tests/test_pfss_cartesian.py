import sys 
sys.path.append('../')
import pfss_cartesian_obj as pieck
import pfss_cartesian as panzer 
""" 
    Attack on Titan trivia: The panzer unit is the artillery unit 
    which rides on the Cart titan.
"""
import sunpy.map as mp
#data url: https://drive.google.com/drive/folders/150sI8Bu3JFmoViN3cpZiocbtWBweeDhj?usp=sharing
BASE = "/home/uvishal/Desktop/IUCAA/pfsspy/test_data/"
path_data = f"{BASE}raster_0_hmi.fits"
hmi_blos=mp.Map(path_data)
print(hmi_blos.data.shape)
dz = 0.1
Z_max = 10
input_obj = pieck.Input(hmi_blos,Z_max,dz)
extrapolated_obj = panzer.fft_pot_open(input_obj)
extrapolated_obj.save(f"{BASE}extrapolation.fits")
