Changelog
=========

0.5.1
-----
- Fixed calculations in `pfsspy.carr_cea_wcs_header`, and clarified in the
  docstring that the input shape must be in ``[nlon, nlat]`` order.
- Added validation to `pfsspy.Input` to check that the inputted map covers the
  whole solar surface.

0.5.0
-----

Changes to outputted maps
~~~~~~~~~~~~~~~~~~~~~~~~~
This release largely sees a transition to leveraging Sunpy Map objects. As such,
the following changes have been made:

`pfsspy.Input` now *must* take a `sunpy.map.GenericMap` as an
input boundary condition (as opposed to a numpy array). To convert a numpy array
to a `~sunpy.map.GenericMap`, the helper function
:func:`pfsspy.carr_cea_wcs_header` can be used::

  map_date = datetime(...)
  br = np.array(...)
  header = pfsspy.carr_cea_wcs_header(map_date, br.shape)

  m = sunpy.map.Map((br, header))
  pfss_input = pfsspy.Input(m, ...)


`pfsspy.Output.source_surface_br` now returns a `~sunpy.map.GenericMap`
instead of an array. To get the data array use ``source_surface_br.data``.

The new `pfsspy.Output.source_surface_pils` returns the coordinates of
the polarity inversion lines on the source surface.

In favour of directly using the plotting functionality built into SunPy,
the following plotting functionality has been removed:

- ``pfsspy.Input.plot_input``. Instead `~pfsspy.Input` has a new
  `~pfsspy.Input.map`  property, which returns a SunPy map, which can easily
  be plotted using `sunpy.map.GenericMap.plot`.
- ``pfsspy.Output.plot_source_surface``. A map of :math:`B_{r}` on the source
  surface can now be obtained using `pfsspy.Output.source_surface_br`, which
  again returns a SunPy map.
- ``pfsspy.Output.plot_pil``. The coordinates of the polarity inversion lines
  on the source surface can now be obtained using
  `pfsspy.Output.source_surface_pils`, which can then be plotted using
  ``ax.plot_coord(pil[0])`` etc. See the examples section for an example.

Specifying tracing seeds
~~~~~~~~~~~~~~~~~~~~~~~~
In order to make specifying seeds easier, they must now be a
`~astropy.coordinates.SkyCoord` object. The coordinates are internally
transformed to the Carrington frame of the PFSS solution, and then traced.

This should make specifying coordinates easier, as lon/lat/r coordinates can
be created using::

  seeds = astropy.coordinates.SkyCoord(lon, lat, r, frame=output.coordinate_frame)

To convert from the old x, y, z array used for seeds, do::

  r, lat, lon = pfsspy.coords.cart2sph
  r = r * astropy.constants.R_sun
  lat = (lat - np.pi / 2) * u.rad
  lon = lon * u.rad

  seeds = astropy.coordinates.SkyCoord(lon, lat, r, frame=output.coordinate_frame)

Note that the latitude must be in the range :math:`[-\pi/2, \pi/2]`.

GONG and ADAPT map sources
~~~~~~~~~~~~~~~~~~~~~~~~~~
pfsspy now comes with built in `sunpy` map sources for GONG and ADAPT synoptic
maps, which automatically fix some non-compliant FITS header values. To use
these, just import ``pfsspy`` and load the .FITS files as normal with sunpy.

Tracing seeds
~~~~~~~~~~~~~
`pfsspy.tracing.Tracer` no longer has a ``transform_seeds`` helper method, which
has been replaced by `~pfsspy.tracing.Tracer.coords_to_xyz` and
`~pfsspy.tracing.Tracer.xyz_to_coords`. These new methods convert
between `~astropy.coordinates.SkyCoord` objects, and Cartesian xyz coordinates
of the internal magnetic field grid.

0.4.3
-----

- Improved the error thrown when trying to use
  :class`pfsspy.tracing.FotranTracer` without the ``streamtracer`` module
  installed.
- Fixed some layout issues in the documentation.

0.4.2
-----

- Fix a bug where :class`pfsspy.tracing.FotranTracer` would overwrite the
  magnetic field values in an `~pfsspy.Output` each time it was used.

0.4.1
-----

- Reduced the default step size for the `~pfsspy.tracing.FortranTracer`
  from 0.1 to 0.01 to give more resolved field lines by default.

0.4.0
-----

New fortran field line tracer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:mod:`pfsspy.tracing` contains a new tracer,
`~pfsspy.tracing.FortranTracer`. This requires and uses the
`streamtracer <https://streamtracer.readthedocs.io/en/stable/>`_ package
which does streamline tracing rapidly in python-wrapped
fortran code. For large numbers of field lines this results in an ~50x
speedup compared to the `~pfsspy.tracing.PythonTracer`.

Changing existing code to use the new tracer is as easy as swapping out
``tracer = pfsspy.tracer.PythonTracer()`` for
``tracer = pfsspy.tracer.FortranTracer()``. If you notice any issues with the
new tracer, please report them at https://github.com/dstansby/pfsspy/issues.

Changes to field line objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``pfsspy.FieldLines`` and ``pfsspy.FieldLine`` have moved to
  `pfsspy.fieldline.FieldLines` and
  `pfsspy.fieldline.FieldLine`.
- `~pfsspy.fieldline.FieldLines` no longer has ``source_surface_feet``
  and ``solar_feet`` properties. Instead these have moved to the new
  `pfsspy.fieldline.OpenFieldLines` class. All the open field lines
  can be accessed from a `~pfsspy.fieldline.FieldLines` instance using
  the new `~pfsspy.fieldline.FieldLines.open_field_lines`
  property.

Changes to `~pfsspy.Output`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- `pfsspy.Output.bg` is now returned as a 4D array instead of three 3D
  arrays. The final index now indexes the vector components; see the docstring
  for more information.

0.3.2
-----
- Fixed a bug in ``pfsspy.FieldLine.is_open``, where some open field lines
  were incorrectly calculated to be closed.

0.3.1
-----
- Fixed a bug that incorrectly set closed line field polarities to -1 or 1
  (instead of the correct value of zero).
- ``FieldLine.footpoints`` has been removed in favour of the new
  ``pfsspy.FieldLine.solar_footpoint`` and
  ``pfsspy.FieldLine.source_surface_footpoint``. These each return a single
  footpoint. For a closed field line, see the API docs for further details
  on this.
- ``pfsspy.FieldLines`` has been added, as a convenience class to store a
  collection of field lines. This means convenience attributes such as
  ``pfsspy.FieldLines.source_surface_feet`` can be used, and their values are
  cached greatly speeding up repeated use.

0.3.0
-----

- The API for doing magnetic field tracing has changed.
  The new :mod:`pfsspy.tracing` module contains `~pfsspy.tracing.Tracer`
  classes that are used to perform the tracing. Code needs to be changed from::

    fline = output.trace(x0)

  to::

    tracer = pfsspy.tracing.PythonTracer()
    tracer.trace(x0, output)
    flines = tracer.xs

  Additionally ``x0`` can be a 2D array that contains multiple seed
  points to trace, taking advantage of the parallelism of some solvers.
- The ``pfsspy.FieldLine`` class no longer inherits from
  `~astropy.coordinates.SkyCoord`, but the
  `~astropy.coordinates.SkyCoord` coordinates are now stored in
  ``pfsspy.FieldLine.coords`` attribute.
- ``pfsspy.FieldLine.expansion_factor`` now returns ``np.nan`` instead of
  ``None`` if the field line is closed.
- ``pfsspy.FieldLine`` now has a ``~pfsspy.FieldLine.footpoints``
  attribute that returns the footpoint(s) of the field line.

0.2.0
-----

- `pfsspy.Input` and `pfsspy.Output` now take the optional keyword
  argument *dtime*, which stores the datetime on which the magnetic field
  measurements were made. This is then propagated to the *obstime* attribute
  of computed field lines, allowing them to be transformed in to coordinate
  systems other than Carrington frames.
- ``pfsspy.FieldLine`` no longer overrrides the SkyCoord ``__init__``;
  this should not matter to users, as FieldLine objects are constructed
  internally by calling `pfsspy.Output.trace`

0.1.5
-----

- ``Output.plot_source_surface`` now accepts keyword arguments that are given to
  Matplotlib to control the plotting of the source surface.

0.1.4
-----

- Added more explanatory comments to the examples
- Corrected the dipole solution calculation
- Added :func:`pfsspy.coords.sph2cart` to transform from spherical to cartesian
  coordinates.

0.1.3
-----

- ``pfsspy.Output.plot_pil`` now accepts keyword arguments that are given
  to Matplotlib to control the style of the contour.
- ``pfsspy.FieldLine.expansion_factor`` is now cached, and is only
  calculated once if accessed multiple times.
