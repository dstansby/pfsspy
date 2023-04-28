.. _changelog:

Changelog
=========

1.1.2
-----
- Added project status documentation.
- Bumped the minimum version of astropy to 5.0.
- Fixed ADAPT map reading with sunpy >= 4.0.

1.1.1
-----
Fixed imports so pfsspy does not depend on ``sympy`` as a runtime dependency.
(``sympy`` is still needed for the ``analytic`` module however).

1.1.0
-----
New requirements
~~~~~~~~~~~~~~~~
pfsspy now depends on Python >= 3.8, and is officially supported with
Python 3.10

New examples
~~~~~~~~~~~~
A host of new examples comparing pfsspy results to analytic solutions have
been added to the example gallery.

Bug fixes
~~~~~~~~~
- Updated the sunpy package requirement to include all packages needed to use
  sunpy maps.
- Any traced field line points that are out of bounds in latitude (ie. have a
  latitude > 90 deg) are now filtered out. This was previously only an issue
  for very low tracing step sizes.

1.0.1
-----
Bug fixes
~~~~~~~~~
- Fixed compatibility of map validity checks with sunpy 3.1.
- Updated this changelog to make it clear that pfsspy 1.0.0 depends on
  sunpy >= 3.0.

1.0.0
-----

New requirements
~~~~~~~~~~~~~~~~
pfsspy now depends on python >= 3.7, sunpy >=3,
and now does *not* depend on Matplotlib.

New features
~~~~~~~~~~~~
- The ``max_steps`` argument to `pfsspy.tracers.FortranTracer` now defaults to
  ``'auto'`` and automatically sets the maximum number of steps to four times the
  number of steps that are needed to span radially from the solar to source
  surface. ``max_steps`` can still be manually specified as a number if more
  or less steps are desired.
- `~pfsspy.fieldline.FieldLines` now has a ``__len__`` method, meaning one
  can now do ``n_field_lines = len(my_field_lines)``.
- Added :func:`pfsspy.utils.roll_map` to roll a map in the longitude direction.
  This is particularly helpful to modify GONG maps so they have a common
  longitude axis.
- Added the `pfsspy.analytic` sub-module that provides functions to sample
  analytic solutions to the PFSS equations.

Bug fixes
~~~~~~~~~
- :func:`pfsspy.utils.carr_cea_wcs_header` now works with versions of sunpy
  >= 2.0.
- GONG synoptic maps now automatically have their observer information corrected
  (by assuming an Earth observer) when loaded by `sunpy.map.Map`.
- The plot settings of input maps are no longer modified in place.

Breaking changes
~~~~~~~~~~~~~~~~
- The interpretation of the ``step_size`` to `pfsspy.tracers.FortranTracer` has
  been corrected so that it is the step size relative to the radial cell size.
  A step size of 0.01 specified in pfsspy<1.0 is approximately equivalent to a
  step size of 1 in pfsspy 1.0, so you will need to adjust any custom step
  sizes accordingly.
- Any points on field lines that are out of bounds (ie. below the solar surface
  or above the source surface) are now removed by the
  `~pfsspy.tracing.FortranTracer`.
- :func:`pfsspy.pfss` no longer warns if the mean of the input data is non-zero,
  and silently ignores the monopole component.

Removals
~~~~~~~~
- Saving and load PFSS solutions is no longer possible. This was poorly tested,
  and possibly broken. If you have interest in saving and loading being added
  as a new feature to pfsspy, please open a new issue at
  https://github.com/dstansby/pfsspy/issues.

0.6.6
-----
Two bugs have been fixed in `pfsspy.utils.carr_cea_wcs_header`:

  - The reference pixel was previously one pixel too large in both longitude and latitude.
  - The longitude coordinate was previously erroneously translated by one degree.

Both of these are now fixed.

0.6.5
-----
This release improves documentation and handling of HMI maps. In particular:

  - The HMI map downloading example has been updated to use the polar filled
    data product, which does not have any data missing at the poles.
  - :func:`pfsspy.utils.fix_hmi_meta` has been added to fix metadata issues in
    HMI maps. This modifies the metadata of a HMI map to make it FITS compliant,
    allowing it to be used with pfsspy.

0.6.4
-----
This release adds citation information to the documentation.

0.6.3
-----
This release contains the source for the accepted JOSS paper describing pfsspy.

0.6.2
-----
This release includes several small fixes in response to a review of pfsspy
for the Journal of Open Source Software. Thanks to Matthieu Ancellin and
Simon Birrer for their helpful feedback!

- A permanent code of conduct file has been added to the repository.
- Information on how to contribute to pfsspy has been added to the docs.
- The example showing the performance of different magnetic field tracers has
  been fixed.
- The docs are now clearer about optional dependencies that can increase
  performance.
- The GONG example data has been updated due to updated data on the remote
  GONG server.

0.6.1
-----

Bug fixes
~~~~~~~~~

- Fixed some messages in errors raised by functions in `pfsspy.utils`.

0.6.0
-----

New features
~~~~~~~~~~~~
- The `pfsspy.utils` module has been added, and contains various tools for
  loading and working with synoptic maps.
- `pfsspy.Output` has a new `~pfsspy.Output.bunit` property, which returns the
  `~astropy.units.Unit` of the input map.
- Added :meth:`pfsspy.Output.get_bvec`, to sample the magnetic field solution
  at arbitrary coordinates.
- Added the `pfsspy.fieldline.FieldLine.b_along_fline` property, to sample the
  magnetic field along a traced field line.
- Added a guide to the numerical methods used by pfsspy.

Breaking changes
~~~~~~~~~~~~~~~~
- The ``.al`` property of `pfsspy.Output` is now private, as it is not intended
  for user access. If you *really* want to access it, use ``._al`` (but this is
  now private API and there is no guarantee it will stay or return the same thing
  in the future).
- A `ValueError` is now raised if any of the input data to `pfsspy.Input` is
  non-finite or NaN. Previously the PFSS computation would run fine, but the
  output would consist entirely of NaNs.

Behaviour changes
~~~~~~~~~~~~~~~~~
- The monopole term is now ignored in the PFSS calculation. Previously a
  non-zero (but small) monopole term would cause floating point precision issues,
  leading to a very noisy result. Now the monopole term is explicitly removed
  from the calculation. If your input has a non-zero mean value, pfsspy will
  issue a warning about this.
- The data downloaded by the examples is now automatically downloaded and
  cached with `sunpy.data.manager`. This means the files used for running the
  examples will be downloaded and stored in your `sunpy` data directory if
  they are required.
- The observer coordinate information in GONG maps is now automatically set
  to the location of Earth at the time in the map header.

Bug fixes
~~~~~~~~~
- The ``date-obs`` FITS keyword in GONG maps is now correctly populated.

0.5.3
-----
- Improved descriptions in the AIA overplotting example.
- Fixed the 'date-obs' keyword in GONG metadata. Previously this just stored
  the date and not the time; now both the date and time are properly stored.
- Drastically sped up the calculation of source surface and solar surface
  magnetic field footpoints.

0.5.2
-----
- Fixed a bug in the GONG synoptic map source where a map failed to load once
  it had already been loaded once.

0.5.1
-----
- Fixed some calculations in ``pfsspy.carr_cea_wcs_header``, and clarified in the
  docstring that the input shape must be in ``[nlon, nlat]`` order.
- Added validation to `pfsspy.Input` to check that the inputted map covers the
  whole solar surface.
- Removed ghost cells from `pfsspy.Output.bc`. This changes the shape of the
  returned arrays by one along some axes.
- Corrected the shape of `pfsspy.Output.bg` in the docstring.
- Added an example showing how to load ADAPT ensemble maps into a
  `~sunpy.map.CompositeMap`
- Sped up field line expansion factor calculations.

0.5.0
-----

Changes to outputted maps
~~~~~~~~~~~~~~~~~~~~~~~~~
This release largely sees a transition to leveraging Sunpy Map objects. As such,
the following changes have been made:

`pfsspy.Input` now *must* take a `sunpy.map.GenericMap` as an
input boundary condition (as opposed to a numpy array). To convert a numpy array
to a `~sunpy.map.GenericMap`, the helper function
``pfsspy.carr_cea_wcs_header`` can be used::

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
``pfsspy.tracing.Tracer.xyz_to_coords``. These new methods convert
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
- Added ``pfsspy.coords.sph2cart`` to transform from spherical to cartesian
  coordinates.

0.1.3
-----

- ``pfsspy.Output.plot_pil`` now accepts keyword arguments that are given
  to Matplotlib to control the style of the contour.
- ``pfsspy.FieldLine.expansion_factor`` is now cached, and is only
  calculated once if accessed multiple times.
