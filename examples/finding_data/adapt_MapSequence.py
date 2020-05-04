"""
Parsing ADAPT Ensemble .fits files
==================================

Parse an ADAPT FITS file into a `sunpy.map.MapSequence`. 
"""

###############################################################################
# Necessary imports
from adapt_helpers import example_adapt_map
import sunpy.map,sunpy.io
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec

###############################################################################
# Load an example ADAPT fits file, utility stored in adapt_helpers.py
adapt_fname = example_adapt_map()

###############################################################################
# ADAPT synoptic magnetograms contain 12 realizations of synoptic magnetograms 
# output as a result of varying model assumptions. See [here](https://www.swpc.noaa.gov/sites/default/files/images/u33/SWW_2012_Talk_04_27_2012_Arge.pdf))
# 
# Because the fits data is 3D, it cannot be passed directly to `sunpy.map.Map`, 
# because this will take the first slice only and the other realizations are 
# lost. We want to end up with a `sunpy.map.MapSequence` containing all these
# realiations as individual maps. These maps can then be individually accessed
# and PFSS solutions generated from them.
#
# We first read in the fits file using `sunpy.io` :
adapt_fits = sunpy.io.fits.read(adapt_fname)

###############################################################################
# `adapt_fits` is a list of `HDPair` objects. The first of these contains the 
# 12 realizations data and a header with sufficient information to build the
# `MapSequence`. We unpack this `HDPair` into a list of `(data,header)` tuples
# where `data` are the different adapt realizations.
data_header_pairs = [(map_slice,adapt_fits[0].header) 
					 for map_slice in adapt_fits[0].data
					 ]

###############################################################################
# Next, pass this list of tuples as the argument to sunpy.map.Map to create the
# map sequence :
adaptMapSequence = sunpy.map.Map(data_header_pairs,sequence=True) 

###############################################################################
# `adapt_map_sequence` is now a list of our individual adapt realizations. Note
# the `.peek()` and `.plot()` methods of `MapSequence` returns instances of
# `sunpy.visualization.MapSequenceAnimator` and 
# `matplotlib.animation.FuncAnimation1` instances. Here, we generate a static
# plot accessing the individual maps in turn :
fig = plt.figure(figsize=(7,8)) 
gs = gridspec.GridSpec(4,3,figure=fig)
for ii,aMap in enumerate(adaptMapSequence) : 
    ax=fig.add_subplot(gs[ii],projection=aMap) 
    aMap.plot(axes=ax,cmap='bwr',vmin=-2,vmax=2,title=f"Realization {1+ii:02d}") 
plt.tight_layout(pad=5,h_pad=2)
plt.show()   
