"""
===================
Removing Light Leak
===================

In this example, we show how to remove the light leak (visible stray light)
from XRT synoptic composite images.
"""

from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from astropy.utils.data import get_pkg_data_path

from xrtpy.image_correction import remove_lightleak

##############################################################################
# This example will be using XRT synoptic data from the first day of summer of 2015.
# This is stored in the ``example_data`` directory of the `xrtpy` package.

directory = get_pkg_data_path("data/example_data", package="xrtpy.image_correction")
data_file = Path(directory) / "comp_XRT20150621_055911.7.fits"
xrt_map = sunpy.map.Map(data_file)

##############################################################################
# Removing the light leak from the composite image is done using the `xrtpy.image_correction.remove_lightleak` function.

lightleak_map = remove_lightleak(xrt_map)

##############################################################################
# Finally, we plot the original and light leak subtracted images side by side.

fig = plt.figure(figsize=(12, 6))

ax = fig.add_subplot(121, projection=xrt_map)
xrt_map.plot(axes=ax, title="Original", clip_interval=(1, 99.9) * u.percent)
ax1 = fig.add_subplot(122, projection=lightleak_map)
lightleak_map.plot(
    axes=ax1, title="Light Leak Subtracted", clip_interval=(1, 99.9) * u.percent
)

ax1.coords[1].set_ticks_visible(False)
ax1.coords[1].set_ticklabel_visible(False)
fig.tight_layout()

plt.show()
