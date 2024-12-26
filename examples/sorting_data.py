"""
=========================
Filtering and Visualizing
=========================

This example provides a simple overview of filtering and visualizing XRT data.
"""

import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from astropy.visualization import ImageNormalize, SqrtStretch
from sunpy.net import Fido
from sunpy.net import attrs as a

##############################################################################
# To start we will download a range of XRT data from the Virtual Solar Observatory (VSO).
# The goal is to acquire a large set of files we can sort through and visualize.

query = Fido.search(
    a.Time("2021-05-21 18:51:00", "2021-05-22 00:00:00"), a.Instrument("xrt")
)
print(query)

##############################################################################
# This query will return a large number of files, this is due to the fact we do not
# specify any additional filters. We can filter the data by specifying additional
# attributes in the query.
#
# For wavelength, we use a range that focuses the data to return only Al-Poly filter images.
# This will cut the results down in half.

query = Fido.search(
    a.Time("2021-05-21 20:51:00", "2021-05-22 00:00:00"),
    a.Instrument("xrt"),
    a.Wavelength(4 * u.nm, 5 * u.nm),
)
print(query)

##############################################################################
# Now we will download the data.
# As this is still over 60 files, this process can take some time.

xrt_files = Fido.fetch(query)

##############################################################################
# We can now load the data into a `~sunpy.map.MapSequence` and create a animation.

xrt_seq = sunpy.map.Map(xrt_files, sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=xrt_seq.maps[0])
ani = xrt_seq.plot(
    axes=ax, norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch())
)

##############################################################################
# You might notice that there is a jump in the sequence.
# The size of the data and the pointing changes.
# We can exclude these images by filtering the data further.

xrt_seq_filtered_shape = sunpy.map.Map(
    [m for m in xrt_seq if m.data.shape == (384, 384)], sequence=True
)

fig = plt.figure()
ax = fig.add_subplot(projection=xrt_seq.maps[0])
ani = xrt_seq_filtered_shape.plot(
    axes=ax, norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch())
)

##############################################################################
# In fact, `sunpy.map.Map` provides many attributes that can be used to filter the data.
# This provides a lot of flexibility in how you can filter the data for your science objective.
#
# For example, we can filter the data by the exposure time or the detector.

xrt_seq_filtered_exp_time = sunpy.map.Map(
    [m for m in xrt_seq_filtered_shape if m.exposure_time < 0.1 * u.s], sequence=True
)

fig = plt.figure()
ax = fig.add_subplot(projection=xrt_seq.maps[0])
ani = xrt_seq_filtered_exp_time.plot(
    axes=ax, norm=ImageNormalize(vmin=0, vmax=5e3, stretch=SqrtStretch())
)

##############################################################################
# If you want to save this animation to a file, you can use the ``save`` method.
# For more information on how to use this method, `see the matplotlib documentation <https://matplotlib.org/stable/api/animation_api.html#animation>`__.

plt.show()
