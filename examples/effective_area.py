"""
=======================
Effective Area Analysis
=======================

In this example, we will explore the effective areas for different XRT filter channels.
Understanding the effective areas is important for accurately interpreting and quantifying the data.
"""

import matplotlib.pyplot as plt

import xrtpy

##############################################################################
# Let us begin by defining a filter channel using its abbreviation.
# For example, if we want to explore the effective area for an aluminum-on-polyimide filter
# channel, we need to specify the relevant abbreviation.

xrt_filter = "Al-poly"

##############################################################################
# `~.EffectiveAreaFundamental` allows us to accurately determine the effective area
# based on the specified filter channel, date, and time.

date_time = "2023-09-22T22:59:59"
eaf = xrtpy.response.EffectiveAreaFundamental(xrt_filter, date_time)

##############################################################################
# To actually calculate the effective area function we can call :meth:`~xrtpy.response.EffectiveAreaFundamental.effective_area`.

effective_area = eaf.effective_area()
print("Effective Area:\n", effective_area)

##############################################################################
# Differences overtime arise from an increase of the contamination layer on the
# CCD which blocks some of the X-rays thus reducing the effective area.
# For detailed information about the calculation of the XRT CCD contaminant layer thickness,
# you can refer to
# `Montana State University <http://solar.physics.montana.edu/HINODE/XRT/xrt_contam_db.html>`__.
#
# Additional information is provided by
# `Narukage et. al. (2011) <https://doi.org/10.1007/s11207-010-9685-2>`__.

relative_launch_date_time = "2006-09-22T22:59:59"
eaf_launch = xrtpy.response.EffectiveAreaFundamental(
    xrt_filter, relative_launch_date_time
)
launch_effective_area = eaf_launch.effective_area()

##############################################################################
# Finally, we can plot how the effective area has changed over time.

plt.figure()

plt.plot(
    eaf.channel_wavelength,
    effective_area,
    label=f"{date_time}",
)
plt.plot(
    eaf.channel_wavelength,
    launch_effective_area,
    label=f"{relative_launch_date_time}",
)

plt.title("XRT Effective Area - Al-Poly")
plt.xlabel("Wavelength (Å)")
plt.ylabel("Effective Area ($cm^{2}$)")
plt.legend()
plt.xlim(0, 60)

plt.grid(color="lightgrey")
plt.tight_layout()

plt.show()
