"""
HMI data
--------
How to search for HMI data.

This example shows how to search for, download, and load HMI data, using the
`sunpy.net.Fido` interface. HMI data is available via. the Joint Stanford
Operations Center (`JSOC`_), and the radial magnetic field synoptic maps come
in two sizes:

- 'hmi.Synoptic_Mr_720s': 3600 x 1440 in (lon, lat)
- 'hmi.mrsynop_small_720s': 720 x 360 in (lon, lat)

For more information on the maps, see the `synoptic maps page`_ on the JSOC site.

.. _JSOC: http://jsoc.stanford.edu/
.. _synoptic maps page: http://jsoc.stanford.edu/HMI/LOS_Synoptic_charts.html
"""

###############################################################################
# First import the required modules
import pfsspy
from sunpy.net import Fido, attrs as a
import sunpy.map

###############################################################################
# Set up the search.
#
# Note that for SunPy versions earlier than 2.0, a time attribute is needed to
# do the search, even if (in this case) it isn't used, as the synoptic maps are
# labelled by Carrington rotation number instead of time
time = a.Time('2010/01/01', '2010/01/01')
series = a.jsoc.Series('hmi.mrsynop_small_720s')

###############################################################################
# Do the search. This will return all the maps in the 'hmi_mrsynop_small_720s
# series.'
result = Fido.search(time, series)
print(result)

###############################################################################
# If we just want to download a specific map, we can specify a Carrington
# rotation number
crot = a.jsoc.PrimeKey('CAR_ROT', 2210)
result = Fido.search(time, series, crot, a.jsoc.Notify("jsoc@cadair.com"))
print(result)

###############################################################################
# Download the files
#
# This will download the files to the default sunpy download directory.
files = Fido.fetch(result)

###############################################################################
# Read in a file. This will read in the first file downloaded to a sunpy Map
# object.
hmi_map = sunpy.map.Map(files[0])
print(hmi_map)
