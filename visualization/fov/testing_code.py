import time as timer
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.time
import astropy.units as u

import sunpy
import ipywidgets as widgets
from ipywidgets import Layout, interact, IntSlider,IntProgress, RadioButtons, FloatSlider,FloatRangeSlider
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact
import numpy as np
from sunpy.net import Fido
from sunpy.net import attrs as a
import matplotlib.dates as mdates

import sys
sys.path.append('/Users/ntrueba/SOLAR/code/GIT/xrtpy/xrtpy/visualization/fov/') #
import xrt_metadata as xmet
import xrt_metadata_plot as xplt


# Define the time range of interest for solar observations
#time_range = a.Time("2011-06-07 06:00:00", "2011-06-07 07:30:54")
time_range = a.Time("2011-06-07 06:00:00", "2011-06-08 07:30:54")
# Specify the instrument as 'xrt' to search for Hinode X-Ray Telescope data
instrument = a.Instrument("xrt")


# This will return a catalog of available XRT data during t≠he specified period
start = timer.time()
xrt_downloaded_files = Fido.search(time_range, instrument)
end = timer.time()
length = end - start
print("It took", length, "seconds!")



###

xmeta = xmet.xrt_meta(xrt_downloaded_files)

####

#%matplotlib inline

xmeta.plot_preview()