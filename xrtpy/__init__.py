"""
A Python data analysis package for the `X-Ray Telescope`_ (XRT) on
Hinode_.
"""


import warnings

from xrtpy import response

try:
    from xrtpy.version import __version__  # noqa
except ImportError:
    warnings.warn("version not found.")


# Then you can be explicit to control what ends up in the namespace,
__all__ = ["response"]
