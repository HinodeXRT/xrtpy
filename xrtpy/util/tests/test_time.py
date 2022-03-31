__all__ = []

import pkg_resources
import sunpy.io.special
import numpy as np

from astropy import units as u


_channel_name_to_index_mapping = {
    "Al-mesh": 0,
    "Al-poly": 1,
    "C-poly": 2,
    "Ti-poly": 3,
    "Be-thin": 4,
    "Be-med": 5,
    "Al-med": 6,
    "Al-thick": 7,
    "Be-thick": 8,
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10,
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12,
    "C-poly/Ti-poly": 13,
}