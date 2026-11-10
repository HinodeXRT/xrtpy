"""
A Python data analysis package for the `X-Ray Telescope`_ (XRT) on
Hinode_.
"""

import warnings
from importlib.metadata import PackageNotFoundError, version

from xrtpy import response

try:
    __version__ = version("xrtpy")
except PackageNotFoundError:
    warnings.warn("version not found.", stacklevel=3)
    __version__ = "0.0.0"

# Then you can be explicit to control what ends up in the namespace,
__all__ = ["__version__", "response"]
