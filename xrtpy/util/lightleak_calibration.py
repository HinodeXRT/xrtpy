"""
Functionality for subtracting light leak (visible stray light) image from XRT synoptic composite images.
"""
# __all__ = ["lightleak_correction"]

import numpy as np

# Import necessary libraries
from pathlib import Path
from scipy.signal import convolve2d
