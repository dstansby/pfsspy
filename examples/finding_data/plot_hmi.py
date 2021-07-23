"""
HMI data
--------
How to search for HMI data.

This example shows how to search for, download, and load HMI data, using the
`sunpy.net.Fido` interface. HMI data is available via. the Joint Stanford
Operations Center (`JSOC`_).

The polar filled radial magnetic field synoptic maps
are obtained using the 'hmi.synoptic_mr_polfil_720s' series keyword. Note that
they are large (1440 x 720), so you will probably want to downsample them to
a smaller resolution to use them to calculate PFSS solutions.

For more information on the maps, see the `synoptic maps page`_ on the JSOC
site.

.. _JSOC: http://jsoc.stanford.edu/
.. _synoptic maps page: http://jsoc.stanford.edu/HMI/LOS_Synoptic_charts.html
"""

###############################################################################
# First import the required modules
import os

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import pfsspy.utils

###############################################################################
# Set up the search.
#
# Note that for SunPy versions earlier than 2.0, a time attribute is needed to
# do the search, even if (in this case) it isn't used, as the synoptic maps are
# labelled by Carrington rotation number instead of time
time = a.Time('2010/01/01', '2010/01/01')
series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')

###############################################################################
# Do the search. This will return all the maps in the 'hmi_mrsynop_small_720s
# series.'
result = Fido.search(time, series)
print(result)

###############################################################################
# If we just want to download a specific map, we can specify a Carrington
# rotation number. In addition, downloading files from JSOC requires a
# notification email. If you use this code, please replace this email address
# with your own one, registered here:
# http://jsoc.stanford.edu/ajax/register_email.html
crot = a.jsoc.PrimeKey('CAR_ROT', 2210)
result = Fido.search(time, series, crot,
                     a.jsoc.Notify(os.environ['JSOC_EMAIL']))
print(result)

###############################################################################
# Download the files. This downloads files to the default sunpy download
# directory.
files = Fido.fetch(result)
print(files)

###############################################################################
# Read in a file. This will read in the first file downloaded to a sunpy Map
# object. Note that HMI maps have several bits of metadata that do not comply
# to the FITS standard, so we need to fix them first.
hmi_map = sunpy.map.Map(files[0])
pfsspy.utils.fix_hmi_meta(hmi_map)
hmi_map.peek()
