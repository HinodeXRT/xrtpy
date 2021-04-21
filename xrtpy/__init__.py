# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
from . import response

try:
    from .version import __version__
except (ImportError, ModuleNotFoundError):
    warnings.warn("version not found.")


# Then you can be explicit to control what ends up in the namespace,
__all__ = ['response']

