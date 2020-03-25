"""Central place for package metadata."""

from pkg_resources import DistributionNotFound, get_distribution  # pylint: disable=unused-import

# NOTE: We use __title__ instead of simply __name__ since the latter would
#       interfere with a global variable __name__ denoting object's name.
__title__ = 'imaps'
__summary__ = 'Bioinformatic pipelines for analysis of CLIP data'
__url__ = 'https://github.com/jernejule/imaps'

try:
    # __version__ = get_distribution(__title__).version
    __version__ = '0.2.3'
except DistributionNotFound:
    # Package is not (yet) installed.
    pass

__author__ = 'Ule lab and Genialis, Inc.'
__email__ = ' jernej.ule@crick.ac.uk'

__license__ = 'Apache 2.0'
__copyright__ = '2019, ' + __author__

__all__ = (
    '__title__', '__summary__', '__url__', '__version__', '__author__',
    '__email__', '__license__', '__copyright__',
)
