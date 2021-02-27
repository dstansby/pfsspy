"""
MDI data
--------
How to search for MDI data.

This example shows how to search for, download, and load HMI data, using the
`sunpy.net.Fido` interface. HMI data is available via. the Joint Stanford
Operations Center (JSOC_), and the radial magnetic field synoptic maps come
in two sizes:

- 'hmi.Synoptic_Mr_720s': 3600 x 1440 in (lon, lat)
- 'hmi.mrsynop_small_720s': 720 x 360 in (lon, lat)

For more information on the maps, see the synoptic maps page_ on the JSOC site.

_JSOC
_synoptic maps page http://jsoc.stanford.edu/HMI/LOS_Synoptic_charts.html
"""
import pfsspy
from sunpy.net import Fido, attrs as a
import sunpy.map
import matplotlib.pyplot as plt


# Note that for SunPy versions earlier than 2.0, a time attribute is needed to
# do the search, even if it isn't used
time = a.Time('2010/01/01', '2010/01/01')
series = a.jsoc.Series('mdi.synoptic_mr_small_96m')

result = Fido.search(time, series)

crot = a.jsoc.PrimeKey('CAR_ROT', 2099)
result = Fido.search(time, series, crot, a.jsoc.Notify("jsoc@cadair.com"))
print(result)

files = Fido.fetch(result)
print(files)

mdi_map = sunpy.map.Map(files[0])
mdi_map.plot()
plt.show()
