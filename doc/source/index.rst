.. image:: ../../logo/logo_rectangle.svg

pfsspy is a python package for carrying out Potential Field Source Surface
modelling, a commonly used magnetic field model of the Sun and other stars.

.. note::

  Until pfsspy 1.0 is released, elements of the API are liable to change between
  versions. A full changelog that lists breaking changes, and how to adapt
  your code for them can be found at :ref:`changelog`.

.. note::
  If you find any bugs or have any suggestions for improvement, please raise
  an issue here: https://github.com/dstansby/pfsspy/issues

Installing
----------

pfsspy requires python >= 3.7, and can be installed from PyPi using

.. code::

    pip install pfsspy

This will install pfsspy and all of its dependencies. In addition to the core
dependencies, there are two optional dependencies (numba, streamtracer) that
improve code performance. These can be installed with

.. code::

    pip install pfsspy[performance]

.. include:: ../../CITING.rst

Contents
--------

.. toctree::
   :maxdepth: 1

   auto_examples/index
   pfsspy
   performance
   changes
   contributing
   numerical_method
   synoptic_fits

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
