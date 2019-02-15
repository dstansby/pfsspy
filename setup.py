from setuptools import setup
import os
import sys
import versioneer

if sys.version_info < (3, 5):
    sys.exit('Python versions older than 3.5 are not supported.')

setup(name='pfsspy',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Potential Field Source Surface model package',
      author='David Stansby',
      license='GPL3',
      author_email='dstansby@gmail.com',
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Physics'],
      url='https://github.com/dstansby/pfsspy',
      install_requires=['numpy', 'scipy', 'matplotlib', 'sunpy', 'astropy>=3'],
      python_requires='>=3.5',
      packages=['pfsspy'],
      )
