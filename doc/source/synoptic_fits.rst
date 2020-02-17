Synoptic map FITS conventions
=============================

FITS is the most common filetype used for the storing of solar images. On this
page the FITS metadata conventions for synoptic maps are collected. All of this
information can be found in, and is taken from,
"Coordinate systems for solar image data (Thompson, 2005)".

=================   ======
   Keyword          Output
-----------------   ------
CRPIX\ *n*          Reference pixel to subtract along axis *n*. Counts from 1 to N.
                    Integer values refer to the centre of the pixel.
CRVAL\ *n*          Coordinate value of the reference pixel along axis *n*.
CDELT\ *n*          Pixel spacing along axis *n*.
CTYPE\ *n*          Coordinate axis label for axis *n*.
PV\ *i*\_\ *m*      Additional parameters needed for some coordinate systems.
=================   ======

Note that *CROTAn* is ignored in this short guide.

Cylindrical equal area projection
---------------------------------
In this projection, the latitude pixels are equally spaced in sin(latitude).
The reference pixel has to be on the equator, to facilitate alignment with the
solar rotation axis.

- CDELT2 is set to 180/pi times the pixel spacing in sin(latitude).
- CTYPE1 is either 'HGLN-CEA' or 'CRLN-CEA'.
- CTYPE2 is either 'HGLT-CEA' or 'CRLT-CEA'.
- PV\ *i*\ _1 is set to 1.
- LONPOLE is 0.

The abbreviations are "Heliographic Longitude - Cylindrical Equal Area" etc.
If the system is heliographic the observer must also be defined in the
metadata.
