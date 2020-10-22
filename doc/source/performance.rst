Improving performance
---------------------

numba
~~~~~
pfsspy automatically detects an installation of `numba`_, which compiles
some of the numerical code to speed up pfss calculations. To enable this
simply `install numba`_  and use pfsspy as normal.

Streamline tracing
~~~~~~~~~~~~~~~~~~
pfsspy has two streamline tracers: a pure python `pfsspy.tracing.PythonTracer`
and a FORTRAN `pfsspy.tracing.FortranTracer`. The FORTRAN version is
significantly faster, using the `streamtracer`_ package.


.. _numba: https://numba.pydata.org
.. _install numba: http://numba.pydata.org/numba-doc/latest/user/installing.html
.. _streamtracer: https://github.com/dstansby/streamtracer
