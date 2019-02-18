Changelog
=========

0.1.3
-----

- :meth:`pfsspy.Output.plot_pil` now accepts keyword arguments that are given
  to Matplotlib to control the style of the contour.
- :attr:`pfsspy.FieldLine.expansion_factor` is now cached, and is only
  calculated once if accessed multiple times.
