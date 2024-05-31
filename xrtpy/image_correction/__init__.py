"""Subpackage for image correction"""

__all__ = [
    "deconvolve",
    "remove_lightleak",
]

from xrtpy.image_correction import deconvolve, remove_lightleak

_SSW_MIRRORS = [
    "https://sohoftp.nascom.nasa.gov/solarsoft/",
    "https://hesperia.gsfc.nasa.gov/ssw/",
]
