"""
====================
Temperature Response
====================

In this example, we will explore the temperature response of the filters on XRT.
The temperature response provides important information on how XRT responds to
the different temperatures of X-ray emissions.
"""

import matplotlib.pyplot as plt
import numpy as np

import xrtpy

##############################################################################
# A filter channel is defined by its common abbreviation, which represents
# a specific type of filter used to modify the X-ray radiation observed.
# In this example, we will explore the carbon-on-polyimide filter (abbreviated as "C-poly").

xrt_filter = "C-poly"

##############################################################################
# `~.TemperatureResponseFundamental` provides the functions and properties for
# calculating the temperature response.

date_time = "2023-09-22T21:59:59"
tpf = xrtpy.response.TemperatureResponseFundamental(
    xrt_filter, date_time, abundance_model="Photospheric"
)

##############################################################################
# To calculate the temperature response,we  can do the following:

temperature_response = tpf.temperature_response()
print("Temperature Response:\n", temperature_response)

##############################################################################
# We will now visualize the temperature response function using CHIANTI.
# These temperatures are of the plasma and are independent of the channel filter.
#
# We use the log of the these temperatures, to enhance the visibility of the
# variations at lower temperatures.

chianti_temperature = np.log10(tpf.CHIANTI_temperature.to_value())

##############################################################################
# Differences overtime arise from an increase of the contamination layer on the
# CCD which blocks some of the X-rays thus reducing the effective area.
# For detailed information about the calculation of the XRT CCD contaminant layer thickness,
# you can refer to
# `Montana State University <http://solar.physics.montana.edu/HINODE/XRT/xrt_contam_db.html>`__.
#
# Additional information is provided by
# `Narukage et. al. (2011) <https://doi.org/10.1007/s11207-010-9685-2>`__.


launch_datetime = "2006-09-22T23:59:59"

launch_temperature_response = xrtpy.response.TemperatureResponseFundamental(
    xrt_filter, launch_datetime, abundance_model="Photospheric"
).temperature_response()

##############################################################################
# Now we can plot the temperature response versus the log of the CHIANTI temperature
# and compare the results for the launch date and the chosen date.

plt.figure()

plt.plot(
    chianti_temperature,
    np.log10(temperature_response.value),
    label=f"{date_time}",
)
plt.plot(
    chianti_temperature,
    np.log10(launch_temperature_response.value),
    label=f"{launch_datetime}",
    color="red",
)

plt.title("XRT Temperature Response")
plt.xlabel("Log(T) ($K$)")
plt.ylabel("$DN$ $cm^5$ $ s^-1$ $pix^-1$")
plt.legend()
plt.grid()

plt.show()
