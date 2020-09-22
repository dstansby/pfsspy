Improving performance
---------------------

numba
~~~~~
pfsspy automatically detects an installation of
`numba <https://numba.pydata.org/>`_, which compiles
some of the numerical code to speed up pfss calculations. To enable this
simply `install numba <http://numba.pydata.org/numba-doc/latest/user/installing.html>`_
and use pfsspy as normal.

Streamline tracing
~~~~~~~~~~~~~~~~~~
pfsspy has two streamline tracers: a pure python `pfsspy.tracing.PythonTracer`
and a FORTRAN `pfsspy.tracing.FortranTracer`. The FORTRAN version is
significantly faster.
