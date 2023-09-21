__all__ = [
    "EffectiveAreaFundamental",
    "effective_area",
]


import astropy.time
import datetime
import math
import numpy as np
import os
import scipy.io
import sunpy.io.special
import sunpy.time

from astropy import units as u
from astropy.utils.data import get_pkg_data_filename
from functools import cached_property
from pathlib import Path
from scipy import interpolate

from xrtpy.response.channel import Channel, resolve_filter_name
from xrtpy.util.time import epoch

index_mapping_to_fw1_name = {
    "Open": 0,
    "Al-poly": 1,
    "C-poly": 2,
    "Be-thin": 3,
    "Be-med": 4,
    "Al-med": 5,
}

index_mapping_to_fw2_name = {
    "Open": 0,
    "Al-mesh": 1,
    "Ti-poly": 2,
    "G-band": 3,
    "Al-thick": 4,
    "Be-thick": 5,
}

_ccd_contam_filename = (
    Path(__file__).parent.absolute() / "data" / "xrt_contam_on_ccd.geny"
)

_filter_contam_filename = (
    Path(__file__).parent.absolute() / "data" / "xrt_contam_on_filter.geny"
)

_ccd_contam_file = scipy.io.readsav(_ccd_contam_filename)
_filter_contam_file = scipy.io.readsav(_filter_contam_filename)

# CCD contam geny files keys for time and date.
_ccd_contamination_file_time = astropy.time.Time(
    _ccd_contam_file["p1"], format="utime", scale="utc"
)
_ccd_contamination = _ccd_contam_file["p2"]

# Filter contam geny files keys for time and date.
_filter_contamination_file_time = astropy.time.Time(
    _filter_contam_file["p1"], format="utime", scale="utc"
)
_filter_contamination = _filter_contam_file["p2"]


class EffectiveAreaFundamental:
    """
    Class for calculating the effective area.

    Parameters
    -----------
    filter_name : str
        The name of the filter.

    observation_date: str
        The date of the observation.  For valid date formats, look at the documentation for
        `sunpy.time.parse_time`.
    """

    def __init__(self, filter_name, observation_date):
        self._name = resolve_filter_name(filter_name)
        self.observation_date = observation_date
        self._channel = Channel(self.name)

    @property
    def name(self) -> str:
        """Name of XRT X-Ray channel filter."""
        return self._name

    @property
    def observation_date(self) -> str:
        """Date of observation."""
        return self._observation_date

    @observation_date.setter
    def observation_date(self, date):
        """Validating users requested observation date."""

        observation_date = sunpy.time.parse_time(date)

        if observation_date <= epoch:
            raise ValueError(
                f"\nInvalid date: {observation_date.datetime}.\n"
                f"Date must be after {epoch}."
            )

        modified_time_path = os.path.getmtime(_ccd_contam_filename)
        modified_time = astropy.time.Time(modified_time_path, format="unix")
        latest_available_ccd_data = _ccd_contamination_file_time[-1].datetime.strftime(
            "%Y/%m/%d"
        )
        modified_time_datetime = datetime.datetime.fromtimestamp(
            modified_time_path
        ).strftime("%Y/%m/%d")

        if observation_date > modified_time:
            raise ValueError(
                "\nNo contamination data is presently available for "
                f"{observation_date.datetime}.\n The latest available data is on "
                f"{latest_available_ccd_data}.\n Contamination data is "
                "updated periodically. The last update was on "
                f"{modified_time_datetime}. If this is more "
                "than one month ago, please raise an issue at: "
                "https://github.com/HinodeXRT/xrtpy/issues/new"
            )

        self._observation_date = observation_date

    @property
    def contamination_on_CCD(self):
        """Calculation of contamination layer on the CCD, thickness given in Angstrom (Å)."""
        interpolater = scipy.interpolate.interp1d(
            _ccd_contamination_file_time.utime, _ccd_contamination, kind="linear"
        )
        return interpolater(self.observation_date.utime)

    @property
    def filter_index_mapping_to_name(self):
        """Returns filter's corresponding number value."""
        if self.name in index_mapping_to_fw1_name:
            return index_mapping_to_fw1_name.get(self.name)
        elif self.name in index_mapping_to_fw2_name:
            return index_mapping_to_fw2_name.get(self.name)

    @property
    def filter_wheel_number(self):
        """Defining chosen filter to its corresponding filter wheel."""
        return 0 if self.name in index_mapping_to_fw1_name else 1

    @property
    def combo_filter_name_split(self):
        """Defining chosen filters to its corresponding filter wheel."""
        name = (self.name).split("/")
        filter1, filter2 = name[0], name[1]
        return filter1, filter2

    @property
    def combo_filter1_wheel_number(self):
        """Defining chosen filter to its corresponding filter wheel."""
        filter1, _ = self.combo_filter_name_split
        return 0 if filter1 in index_mapping_to_fw1_name else 1

    @property
    def combo_filter2_wheel_number(self):
        """Defining chosen filter to its corresponding filter wheel."""
        _, filter2 = self.combo_filter_name_split
        return 0 if filter2 in index_mapping_to_fw1_name else 1

    @property
    def combo_filter_index_mapping_to_name_filter1(self):
        """Returns filter's corresponding number value."""
        filter1, _ = self.combo_filter_name_split

        if filter1 in index_mapping_to_fw1_name:
            return index_mapping_to_fw1_name.get(filter1)
        elif filter1 in index_mapping_to_fw2_name:
            return index_mapping_to_fw2_name.get(filter1)

    @property
    def combo_filter_index_mapping_to_name_filter2(self):
        """Returns filter's corresponding number value."""
        filter1, filter2 = self.combo_filter_name_split

        if filter2 in index_mapping_to_fw1_name:
            return index_mapping_to_fw1_name.get(filter2)
        elif filter2 in index_mapping_to_fw2_name:
            return index_mapping_to_fw2_name.get(filter2)

    @property
    def combo_filter1_data(self):
        """Collecting filter data."""
        return _filter_contamination[self.combo_filter_index_mapping_to_name_filter1][
            self.combo_filter1_wheel_number
        ]

    @property
    def combo_filter2_data(self):
        """Collecting filter data."""
        return _filter_contamination[self.combo_filter_index_mapping_to_name_filter2][
            self.combo_filter2_wheel_number
        ]

    @property
    def contamination_on_filter1_combo(self) -> u.angstrom:
        """
        Thickness of the contamination layer on a filter."""

        interpolater = scipy.interpolate.interp1d(
            self.filter_data_dates_to_seconds, self.combo_filter1_data, kind="linear"
        )
        return interpolater(self.filter_observation_date_to_seconds)

    @property
    def contamination_on_filter2_combo(self) -> u.angstrom:
        """
        Thickness of the contamination layer on a filter."""

        interpolater = scipy.interpolate.interp1d(
            self.filter_data_dates_to_seconds, self.combo_filter2_data, kind="linear"
        )
        return interpolater(self.filter_observation_date_to_seconds)

    @property
    def filter_data(self):
        """Collecting filter contamination data."""
        return _filter_contamination[self.filter_index_mapping_to_name][
            self.filter_wheel_number
        ]

    @property
    def contamination_on_filter(self) -> u.angstrom:
        """Thickness of the contamination layer on a filter."""
        interpolater = scipy.interpolate.interp1d(
            _filter_contamination_file_time.utime, self.filter_data, kind="linear"
        )
        return interpolater(self.observation_date.utime)

    @cached_property
    def n_DEHP_attributes(self):
        """(Diethylhexylphthalate) Wavelength (nm), Delta, Beta."""
        _n_DEHP_filename = get_pkg_data_filename(
            "data/n_DEHP.txt", package="xrtpy.response"
        )

        with open(_n_DEHP_filename) as n_DEHP:
            list_of_DEHP_attributes = []
            for line in n_DEHP:
                stripped_line = line.strip()
                line_list = stripped_line.split()
                list_of_DEHP_attributes.append(line_list)

        return list_of_DEHP_attributes

    @cached_property
    def n_DEHP_wavelength(self):
        """(Diethylhexylphthalate) Wavelength given in Angstrom (Å)."""

        # Convert wavelength values from nanometers to Angstroms
        wavelength_str = [
            self.n_DEHP_attributes[i][0] for i in range(2, len(self.n_DEHP_attributes))
        ]

        return np.array([float(i) * 10 for i in wavelength_str])

    @cached_property
    def n_DEHP_delta(self):
        """(Diethylhexylphthalate) Delta."""

        delta_str = [
            self.n_DEHP_attributes[i][1] for i in range(2, len(self.n_DEHP_attributes))
        ]

        # Converting from str to float
        delta_float = np.array(
            [float(delta_str[i]) for i in range(len(self.n_DEHP_wavelength))]
        )

        # Interpolate so ranges are the same
        return interpolate.interp1d(self.n_DEHP_wavelength, delta_float)(
            self.n_DEHP_wavelength
        )

    @cached_property
    def n_DEHP_beta(self):
        """(Diethylhexylphthalate) Beta."""

        beta_str = [
            self.n_DEHP_attributes[i][2] for i in range(2, len(self.n_DEHP_attributes))
        ]

        # Converting from str to float
        beta_float = np.array(
            [float(beta_str[i]) for i in range(len(self.n_DEHP_wavelength))]
        )

        # Interpolate so ranges are the same
        return interpolate.interp1d(self.n_DEHP_wavelength, beta_float)(
            self.n_DEHP_wavelength
        )

    @cached_property
    def transmission_equation(self):
        """Defining equations that will be used to calculate the effective area.
        REFERENCES: G.R Fowles, Intro to Modern Optics 2nd Edition, pp 96-101."""

        n_o = 1.0  # index of medium at entrance of filter (assumed vacuum)
        n_t = 1.0  # index of medium at exit of filter (assumed vacuum)

        incidence_angle = 0  # Angle of incidence on Filter in radians

        wavelength_max = 4000  # Max wavelength in Angstroms

        index = [
            (complex((1 - self.n_DEHP_delta[i]), (1.0 * self.n_DEHP_beta[i])))
            for i in range(4000)
        ]

        # Snell's law
        sin_a = (n_o * np.sin(incidence_angle)) / index

        cos_a = 1

        return (index, sin_a, cos_a, wavelength_max, n_o, n_t, incidence_angle)

    @cached_property
    def angular_wavenumber_CCD(self):
        """Define angular wavenumber on CCD."""

        index, _, cos_a, wavelength_max, _, _, _ = self.transmission_equation

        # Define wavevector
        angular_wavenumber = np.array(
            [
                (2.0 * math.pi * index[i] * cos_a) / self.n_DEHP_wavelength[i]
                for i in range(4000)
            ]
        )

        # Multiply by thickness
        angular_wavenumber_thickness = angular_wavenumber * self.contamination_on_CCD

        real_angular_wavenumber = angular_wavenumber_thickness.real
        imaginary_angular_wavenumber = angular_wavenumber_thickness.imag

        return [
            (complex(real_angular_wavenumber[i], imaginary_angular_wavenumber[i]))
            for i in range(4000)
        ]

    @cached_property
    def filterwheel_angular_wavenumber(self):
        """Define angular wavenumber for a filter."""
        index, _, cos_a, _, _, _, _ = self.transmission_equation

        # Define wavevector
        angular_wavenumber = np.array(
            [
                (2.0 * math.pi * index[i] * cos_a) / self.n_DEHP_wavelength[i]
                for i in range(4000)
            ]
        )

        # Multiply by thickness
        angular_wavenumber_thickness = angular_wavenumber * self.contamination_on_filter

        real_angular_wavenumber = angular_wavenumber_thickness.real
        imaginary_angular_wavenumber = angular_wavenumber_thickness.imag

        return [
            (complex(real_angular_wavenumber[i], imaginary_angular_wavenumber[i]))
            for i in range(4000)
        ]

    @cached_property
    def CCD_contamination_transmission(self):
        """Calculate transmission matrix coefficient and transmittance on the CCD."""

        index, _, _, _, n_o, n_t, _ = self.transmission_equation

        i_i = complex(0, 1)  # Define complex number

        # Define transfer matrix
        M = [
            [
                [
                    np.cos(self.angular_wavenumber_CCD[i]),
                    (-i_i * np.sin(self.angular_wavenumber_CCD[i])) / index[i],
                ],
                [
                    -i_i * np.sin(self.angular_wavenumber_CCD[i]) * index[i],
                    np.cos(self.angular_wavenumber_CCD[i]),
                ],
            ]
            for i in range(4000)
        ]

        transmittance = [
            2
            * n_o
            / (
                (M[i][0][0] * n_o)
                + (M[i][0][1] * n_o * n_t)
                + (M[i][1][0])
                + (M[i][1][1] * n_t)
            )
            for i in range(4000)
        ]

        return np.array([abs(transmittance[i] ** 2) for i in range(4000)])

    @property
    def channel_wavelength(self):
        """Array of wavelengths for every X-ray channel in Angstroms (Å)."""
        return Channel(self.name).wavelength

    @property
    def channel_geometry_aperture_area(self):
        """XRT flight model geometry aperture area."""
        return Channel(self.name).geometry.geometry_aperture_area

    @property
    def channel_transmission(self):
        """XRT channel transmission."""
        return Channel(self.name).transmission

    @property
    def interpolated_CCD_contamination_transmission(self):
        """Interpolate filter contam transmission to the wavelength."""
        CCD_contam_transmission = interpolate.interp1d(
            self.n_DEHP_wavelength, self.CCD_contamination_transmission
        )
        return CCD_contam_transmission(self.channel_wavelength)

    @cached_property
    def filter_contamination_transmission(self):
        """Calculate transmission matrix coefficient and transmittance on a filter."""

        index, _, _, _, n_o, n_t, _ = self.transmission_equation

        i_i = complex(0, 1)  # Define complex number

        # Define transfer matrix
        M = [
            [
                [
                    np.cos(self.filterwheel_angular_wavenumber[i]),
                    (-i_i * np.sin(self.filterwheel_angular_wavenumber[i])) / index[i],
                ],
                [
                    -i_i * np.sin(self.filterwheel_angular_wavenumber[i]) * index[i],
                    np.cos(self.filterwheel_angular_wavenumber[i]),
                ],
            ]
            for i in range(4000)
        ]

        transmittance = [
            2
            * n_o
            / (
                (M[i][0][0] * n_o)
                + (M[i][0][1] * n_o * n_t)
                + (M[i][1][0])
                + (M[i][1][1] * n_t)
            )
            for i in range(4000)
        ]

        return [abs(transmittance[i] ** 2) for i in range(4000)]

    @property
    def interpolated_filter_contamination_transmission(self):
        """Interpolate filter contam transmission to the wavelength."""
        Filter_contam_transmission = interpolate.interp1d(
            self.n_DEHP_wavelength, self.filter_contamination_transmission
        )
        return Filter_contam_transmission(self.channel_wavelength)

    @u.quantity_input
    def effective_area(self) -> u.cm**2:
        """Calculation of the Effective Area."""
        return (
            self.channel_geometry_aperture_area
            * self.channel_transmission
            * self.interpolated_CCD_contamination_transmission
            * self.interpolated_filter_contamination_transmission
        )


def effective_area(filter_name, observation_date):
    EAP = EffectiveAreaFundamental(filter_name, observation_date)
    return EAP.effective_area()
