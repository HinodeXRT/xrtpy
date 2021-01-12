# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

try:
    from .version import __version__
except (ImportError, ModuleNotFoundError):
    warnings.warn("version not found.")

from .example_mod import do_primes
# Then you can be explicit to control what ends up in the namespace,
__all__ = ['do_primes']

