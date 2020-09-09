.. image:: ../../logo/logo_rectangle.svg

pfsspy is a python package for carrying out Potential Field Source Surface
modelling, a commonly used magnetic field model of the Sun and other stars.
For more information on the actually PFSS calculation see
`this document <https://github.com/antyeates1983/pfss/blob/master/pfss-manual.pdf>`_.

.. note::

  pfsspy is a very new package, so elements of the API are liable to change with
  the first few releases.
  If you find any bugs or have any suggestions for improvement, please raise
  an issue here: https://github.com/dstansby/pfsspy/issues

Installing
----------

pfsspy can be installed from PyPi using

.. code::

    pip install pfsspy

Improving performance
---------------------

numba
~~~~~
pfsspy automatically detects an installation of
`numba <https://numba.pydata.org/>`_, which compiles
some of the numerical code to speed up pfss calculations. To enable this
simply `install numba <http://numba.pydata.org/numba-doc/latest/user/installing.html>`_
and use pfsspy as normal.

Citing
------

If you use pfsspy in work that results in publication, please cite the archived
code at *both*

  - https://zenodo.org/record/2566462
  - https://zenodo.org/record/1472183

Citation details can be found at the lower
right hand of each web page.

Code reference
--------------

For the main user-facing code and a changelog see

.. toctree::
   :maxdepth: 1

   pfsspy
   changes

for usage examples see

.. toctree::
   :maxdepth: 2

   auto_examples/index

for various helper functions for working with synoptic maps:

.. toctree::
   :maxdepth: 2

   utils

and for a quick reference guide to synoptic map FITS conventions see

.. toctree::
   :maxdepth: 1

   synoptic_fits

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
