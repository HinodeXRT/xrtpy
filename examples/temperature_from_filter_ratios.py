"""
================================================
Calculating the temperature and emission measure
================================================

In this example, we will showcase how to use the filter method to calculate
the temperature and emission measure of the X-ray Telescope (XRT) on Hinode.
"""

import matplotlib.pyplot as plt
import sunpy.map
from astropy.visualization import ImageNormalize, LogStretch
from sunpy.net import Fido
from sunpy.net import attrs as a

from xrtpy.response import temperature_from_filter_ratio

##############################################################################
# To start, we will get XRT data via ``sunpy``.
#
# It is important to use images that same size and with the smallest time separation.
# Note that not all filter ratios produce good results.

query = Fido.search(
    a.Time("2011-01-28 01:31:55", "2011-01-28 01:32:05"), a.Instrument("xrt")
)
data_files = Fido.fetch(query)
xrt_map_1 = sunpy.map.Map(data_files[0])
xrt_map_2 = sunpy.map.Map(data_files[1])

##############################################################################
# The `xrtpy.response.temperature_from_filter_ratio` function has several options, mirroring
# the IDL routine xrt_teem.pro in SolarSoft in most respects.A simple call with
# no extra parameters calculates the temperature and (volume) emission measure
# for the two images without any binning or masking of the data.

T_EM = temperature_from_filter_ratio(xrt_map_1, xrt_map_2)

##############################################################################
# The output is a namedtuple with attributes ``Tmap``, ``EMmap``, ``Terrmap``, and ``EMerrmap``.
# As with the SolarSoft IDL routine xrt_teem.pro, the output images are logs of the quantities.
#
# ``Tmap`` is the electron temperature, ``EMmap`` is the volume emission measure, ``Terrmap``
# is a measure of the uncertainties in the temperature determined for each pixel and ``EMerrmap``
# is the same for the emission measure.

T_e = T_EM.Tmap

fig = plt.figure()

ax = plt.subplot(projection=T_e)
T_e.plot(
    title="Derived Temperature",
    cmap="inferno",
    norm=ImageNormalize(vmin=6, vmax=7, stretch=LogStretch(10)),
)
T_e.draw_limb()
T_e.draw_grid()
plt.colorbar(label="T (K)")
plt.tight_layout()

plt.show()

##############################################################################
# If you want to do the same for Level 2 synoptic composite images, you have to use
# `~.make_exposure_map` to generate the exposure maps for the composite images.
# This is then passed to `xrtpy.response.temperature_from_filter_ratio` as the
# ``expmap1`` and ``expmap2`` arguments.
# Otherwise without accounting for the different exposure time per pixel,
# the temperature and emission measure will be incorrect.
