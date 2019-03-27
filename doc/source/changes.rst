Changelog
=========

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
