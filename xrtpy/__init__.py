"""
A Python data analysis package for the `X-Ray Telescope`_ (XRT) on
Hinode_.
"""

import warnings

from xrtpy import response

try:
    from xrtpy.version import __version__
except ImportError:
    warnings.warn("version not found.", stacklevel=3)
    __version__ = "0.0.0"

# Then you can be explicit to control what ends up in the namespace,
__all__ = ["response", "__version__"]
