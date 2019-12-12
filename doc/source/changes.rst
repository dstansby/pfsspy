Changelog
=========

0.4.0
-----
- :class:`pfsspy.FieldLines` no longer has ``source_surface_feet`` and
  ``solar_feet`` properties. Instad these have moved to the new
  :class:`pfsspy.OpenFieldLines` class. All the open field lines can be accessed
  using the new :attr:`pfsspy.FieldLines.open_field_lines`.

0.3.2
-----
- Fixed a bug in :attr:`pfsspy.FieldLine.is_open`, where some open field lines
  were incorrectly calculated to be closed.

0.3.1
-----
- Fixed a bug that incorrectly set closed line field polarities to -1 or 1
  (instead of the correct value of zero).
- :attr:`FieldLine.footpoints` has been removed in favour of the new
  :attr:`pfsspy.FieldLine.solar_footpoint` and
  :attr:`pfsspy.FieldLine.source_surface_footpoint`. These each return a single
  footpoint. For a closed field line, see the API docs for further details
  on this.
- :class:`pfsspy.FieldLines` has been added, as a convenience class to store a
  collection of field lines. This means convenience attributes such as
  :attr:`pfsspy.FieldLines.source_surface_feet` can be used, and their values are
  cached greatly speeding up repeated use.

0.3.0
-----

- The API for doing magnetic field tracing has changed.
  The new :mod:`pfsspy.tracing` module contains :class:`~pfsspy.tracing.Tracer`
  classes that are used to perform the tracing. Code needs to be changed from::

    fline = output.trace(x0)

  to::

    tracer = pfsspy.tracing.PythonTracer()
    tracer.trace(x0, output)
    flines = tracer.xs

  Additionally ``x0`` can be a 2D array that contains multiple seed
  points to trace, taking advantage of the parallelism of some solvers.
- The :class:`pfsspy.FieldLine` class no longer inherits from
  :class:`~astropy.coordinates.SkyCoord`, but the
  :class:`~astropy.coordinates.SkyCoord` coordinates are now stored in
  :attr:`pfsspy.FieldLine.coords` attribute.
- :attr:`pfsspy.FieldLine.expansion_factor` now returns ``np.nan`` instead of
  ``None`` if the field line is closed.
- :class:`pfsspy.FieldLine` now has a :attr:`~pfsspy.FieldLine.footpoints`
  attribute that returns the footpoint(s) of the field line.

0.2.0
-----

- :class:`pfsspy.Input` and :class:`pfsspy.Output` now take the optional keyword
  argument *dtime*, which stores the datetime on which the magnetic field
  measurements were made. This is then propagated to the *obstime* attribute
  of computed field lines, allowing them to be transformed in to coordinate
  systems other than Carrignton frames.
- :class:`pfsspy.FieldLine` no longer overrrides the SkyCoord ``__init__``;
  this should not matter to users, as FieldLine objects are constructed
  internally by calling :meth:`pfsspy.Output.trace`

0.1.5
-----

- `Output.plot_source_surface` now accepts keyword arguments that are given to
  Matplotlib to control the plotting of the source surface.

0.1.4
-----

- Added more explanatory comments to the examples
- Corrected the dipole solution calculation
- Added :func:`pfsspy.coords.sph2cart` to transform from spherical to cartesian
  coordinates.

0.1.3
-----

- :meth:`pfsspy.Output.plot_pil` now accepts keyword arguments that are given
  to Matplotlib to control the style of the contour.
- :attr:`pfsspy.FieldLine.expansion_factor` is now cached, and is only
  calculated once if accessed multiple times.
