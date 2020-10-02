from setuptools import setup
import sys
import versioneer

if sys.version_info < (3, 5):
    sys.exit('Python versions older than 3.5 are not supported.')

setup(version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass())
