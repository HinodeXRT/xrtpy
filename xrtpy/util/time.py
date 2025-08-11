"""
Provides the Hinode/XRT mission epoch as an astropy Time object.

The XRT epoch corresponds to mission elapsed time zero:
September 22, 2006 at 21:36:00 UTC.
"""

from astropy.time import Time

__all__ = ["epoch"]

# XRT mission epoch
epoch = Time("2006-09-22 21:36:00", format="iso", scale="utc")
