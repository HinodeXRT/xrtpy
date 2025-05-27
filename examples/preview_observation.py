"""
=============================
FOV tool test
=============================

We are just testing stuff here
"""

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
import metadata_manager as ObsMeta

##############################################################################
# FIDO SEARCH
# The tool is meant to preview metadata within the fido ecosystem - this means that fido is the slowest part of the code unless you are looking at a very small time frame

# Define the time range of interest for solar observations
time_range_1 = a.Time("2011-06-07 06:00:00", "2011-06-07 06:45:54")
time_range_2 = a.Time("2007-12-17 10:40:00", "2007-12-17 13:00:54")

# Specify the instrument as 'xrt' to search for Hinode X-Ray Telescope data
instrument = a.Instrument("xrt")


# This will return a catalog of available XRT data during t≠he specified period
xrt_downloaded_files_1 = Fido.search(time_range_1, instrument)
xrt_downloaded_files_2 = Fido.search(time_range_2, instrument)



##############################################################################
#Metadata extraction and the XRT_meta structure
#The first function accepts the output of the fido search as an input and retrieves the corresponding metadata without having to download the data.
#By default, it downloads the level0 metadata (same as SSWIDL) which is very quick when dealing with many files. If you want the level 1 metadata, the syntax is (.., fast_bool = False) - this is only recommended for short observations (< 1hr), as it can retrieve metadata at a rate 5 obs/second
#The second function creates a handy metadata object for all the observations
#it contains a list of headers (xmeta.head_lis)
#and a dictonary (xmeta.metadata) containing important filter-separated quantities
#This could be super handy, as you get much more information than what fido gives you, so you can download observations that meet very specific conditions - let me know if you want something specific included here



xrt_dset1 = ObsMeta.DatasetMetaManager(xrt_downloaded_files_1)
xrt_dset2 = ObsMeta.DatasetMetaManager(xrt_downloaded_files_2)##

##############################################################################
#Plotting backend
#The %matplotlib inline works best when using notebooks to avoid flickering. Ipywidgets is not available with HTML, but you can uncomment this line when you download the notebook.
#%matplotlib inline
###
ani = xrt_dset2.plot_preview(ani_bool = False, d_mode=False)
plt.show()

##############################################################################
###
#We can also do nightmode 

ani = xrt_dset1.plot_preview(ani_bool = False, d_mode=True, vertical_plot=True)##
plt.show()

##############################################################################
###
#We can also do a horizontal version

ani = xrt_dset2.plot_preview(ani_bool = False, d_mode=True, vertical_plot=False)
plt.show()