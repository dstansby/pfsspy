import astropy
import distutils.version

# Import pfsspy sub-modules to have them available through pfsspy.{name}
import pfsspy.coords
import pfsspy.tracing
import pfsspy.fieldline
# Import this to register map sources
import pfsspy.map
import pfsspy.utils
import pfsspy.sample_data

from .pfss import pfss
from .output import Output, load_output
from .input import Input

__all__ = ['Input', 'Output', 'load_output', 'pfss']


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# Do a version check for astropy
if (distutils.version.LooseVersion(astropy.__version__) <
        distutils.version.LooseVersion("3")):
    raise RuntimeError('pfsspy requires astropy v3 to run ' +
                       f'(found version {astropy.__version__} installed)')
