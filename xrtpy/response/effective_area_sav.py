__all__ = [
    "EffectiveAreaFundamental",
]

import datetime
import math
from functools import cached_property
from pathlib import Path
from typing import NamedTuple
from functools import cached_property

import astropy.time
import numpy as np
import scipy.interpolate
import scipy.io
import sunpy.io.special
import sunpy.time
from astropy import units as u
from astropy.utils.data import get_pkg_data_filename
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

index_mapping_to_multi_filter = {
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10,
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12,
    "C-poly/Ti-poly": 13,
}

####
#from scipy.io import readsav
#data = readsav(Path(__file__).parent.absolute() / "data" / "xrt_contam_on_ccd.sav")
# _ccd_contam_filename = (
#     Path(__file__).parent.absolute() / "data" / "xrt_contam_on_ccd.sav"
# )

# print(data.keys())

_ccd_contam_filename = (
    Path(__file__).parent.absolute() / "data" / "xrt_contam_on_ccd.sav"
)
#import pdb; pdb.set_trace()
######

# _ccd_contam_filename = (
#     Path(__file__).parent.absolute() / "data" / "xrt_contam_on_ccd.geny"
# )

_filter_contam_filename = (
    Path(__file__).parent.absolute() / "data" / "xrt_contam_on_filter.geny"
)

#_ccd_contam_file = scipy.io.readsav(_ccd_contam_filename)

# CCD contam geny files keys for time and date.
# _ccd_contamination_file_time = astropy.time.Time(
#     _ccd_contam_file["p1"], format="utime", scale="utc"
# )
# _ccd_contamination = _ccd_contam_file["p2"]

#_filter_contam_file = scipy.io.readsav(_filter_contam_filename)
@cached_property
def _filter_contam_file(self):
    import scipy.io
    return scipy.io.readsav(_filter_contam_filename)


# Filter contam geny files keys for time and date.
_filter_contamination_file_time = astropy.time.Time(
    _filter_contam_file["p1"], format="utime", scale="utc"
)
_filter_contamination = _filter_contam_file["p2"]


class ParsedFilter(NamedTuple):
    filter1: str
    filter2: str | None
    is_combo: bool


def parse_filter_input(filter_string):
    """
    Parses and validates a filter input string using official filter mappings.

    Parameters
    ----------
    filter_string : str
        A string representing either a single filter (e.g., "Al-poly")
        or a filter combo (e.g., "Al-poly/Ti-poly").

    Returns
    -------
    ParsedFilter
        Named tuple with 'filter1', 'filter2', and 'is_combo' boolean.

    Raises
    ------
    ValueError
        If the filter or combo is invalid.
    """
    if not isinstance(filter_string, str):
        raise TypeError("Filter name must be a string.")

    standardized = resolve_filter_name(filter_string.strip())

    if "/" in standardized:
        if standardized not in index_mapping_to_multi_filter:
            raise ValueError(
                f"'{standardized}' is not a valid filter combination.\n"
                f"Valid combinations are: {sorted(index_mapping_to_multi_filter)}"
            )
        f1, f2 = standardized.split("/")
        return ParsedFilter(filter1=f1, filter2=f2, is_combo=True)

    elif (
        standardized in index_mapping_to_fw1_name
        or standardized in index_mapping_to_fw2_name
    ):
        return ParsedFilter(filter1=standardized, filter2=None, is_combo=False)

    else:
        raise ValueError(
            f"'{standardized}' is not a recognized filter in either filter wheel.\n"
            f"Valid FW1 filters: {sorted(index_mapping_to_fw1_name)}\n"
            f"Valid FW2 filters: {sorted(index_mapping_to_fw2_name)}"
        )


class EffectiveAreaFundamental:
    """
    Class for calculating the effective area for an XRT filter at a specific observation date.

    This class handles the calculations required to determine the effective area of a filter
    used in the X-Ray Telescope (XRT) on the Hinode satellite. It considers various factors such as
    contamination on the CCD and filters, as well as the geometry and transmission of the XRT channel.

    Parameters
    ----------
    filter_name : str
        The name of the filter.
    observation_date : str or datetime.datetime
        The date of the observation. Acceptable formats include any string or datetime object
        that can be parsed by `sunpy.time.parse_time`.
    """

    def __init__(self, filter_name, observation_date):
        self._raw_input_name = filter_name
        self._parsed_filter = parse_filter_input(filter_name)

        self._filter1_name = self._parsed_filter.filter1
        self._filter2_name = self._parsed_filter.filter2
        self._is_combo = self._parsed_filter.is_combo

        # Store the standardized/resolved name
        if self._is_combo:
            self._name = f"{self._filter1_name}/{self._filter2_name}"
        else:
            self._name = self._filter1_name

        self._observation_date = sunpy.time.parse_time(observation_date)
        self._channel = Channel(self.name)

        self.observation_date = observation_date
        self._channel = Channel(self.name)

    @property
    def name(self) -> str:
        """
        The resolved name of the filter or filter combination.
        """
        return self._name

    @property
    def filter1_name(self) -> str:
        """
        Name of the first filter given.
        """
        return self._filter1_name

    @property
    def filter2_name(self) -> str:
        """
        Name of the second filter given.
        """
        return self._filter2_name

    @property
    def is_combo(self) -> bool:
        """
        True if the filter is a combination (e.g., 'Al-poly/Ti-poly'), False otherwise.
        """
        return self._is_combo

    @property
    def observation_date(self) -> str:
        """
        Date of observation.
        """
        return self._observation_date

    @observation_date.setter
    def observation_date(self, date):
        """
        Validates the user's requested observation date.
        """

        observation_date = sunpy.time.parse_time(date)

        if observation_date <= epoch:
            raise ValueError(
                f"\nInvalid date: {observation_date.datetime}.\n"
                f"Date must be after {epoch}."
            )

        modified_time_path = Path(_ccd_contam_filename).stat().st_mtime
        modified_time = astropy.time.Time(modified_time_path, format="unix")
        #REMOVE
        # latest_available_ccd_data = _ccd_contamination_file_time[-1].datetime.strftime(
        #     "%Y/%m/%d"
        # )
        latest_available_ccd_data = self._ccd_contamination_file_time[-1].datetime.strftime(
            "%Y/%m/%d"
        )

        modified_time_datetime = datetime.datetime.fromtimestamp(
            modified_time_path
        ).strftime("%Y/%m/%d")

        if observation_date > modified_time:
            raise ValueError(
                "\nNo contamination data is presently available for "
                f"{observation_date.datetime}.\nThe latest available data is on "
                f"{latest_available_ccd_data}.\nContamination data is "
                "updated periodically. The last update was on "
                f"{modified_time_datetime}. If this is more "
                "than one month ago, please raise an issue at: "
                "https://github.com/HinodeXRT/xrtpy/issues/new"
            )

        self._observation_date = observation_date

    @cached_property
    def ccd_contam_data(self):
        import scipy.io
        return scipy.io.readsav(_ccd_contam_filename)
    
    @cached_property
    def _ccd_contamination_file_time(self):
        """
        CCD contamination time axis from the .sav file.
        """
        return astropy.time.Time(
            self.ccd_contam_data["p1"], format="utime", scale="utc"
        )

    @cached_property
    def _ccd_contamination(self):
        """
        CCD contamination thickness values from the .sav file.
        """
        return self.ccd_contam_data["p2"]


    
    @property
    def contamination_on_CCD(self):
        """
        Calculate the thickness of the contamination layer on the CCD.

        This property interpolates the contamination data over time to determine the thickness
        of the contamination layer on the CCD at the observation date. The contamination layer
        is measured in Angstroms (Å).

        Returns
        -------
        astropy.units.Quantity
            The thickness of the contamination layer on the CCD, in Angstroms.

        Notes
        -----
        The interpolation is performed using a linear interpolation method over the
        available contamination data points. The `observation_date` attribute is used to
        provide the point at which to evaluate the interpolation.

        Raises
        ------
        ValueError
            If the observation date is outside the range of the available contamination data.
        """
        #REMOVE
        # interpolater = scipy.interpolate.interp1d(
        #     _ccd_contamination_file_time.utime, _ccd_contamination, kind="linear"
        # )
        interpolater = scipy.interpolate.interp1d(
            self._ccd_contamination_file_time.utime, self.filter2_name_ccd_contamination, kind="linear"
        )
        return interpolater(self.observation_date.utime)

    @staticmethod
    def _get_filter_wheel_number(filter_name: str) -> int:  # *********
        """Returns 0 if the filter is in FW1, 1 if in FW2. Raises error if unknown."""
        if filter_name in index_mapping_to_fw1_name:
            return 0
        elif filter_name in index_mapping_to_fw2_name:
            return 1
        else:
            raise ValueError(f"Filter '{filter_name}' not found in FW1 or FW2.")

    @property
    def _filter_wheel_number(self):
        """Defining chosen filter to its corresponding filter wheel."""
        return 0 if self.name in index_mapping_to_fw1_name else 1

    @property
    def _combo_filter_name_split(self):
        """Returns (filter1, filter2) even for a single filter (filter2 is None)."""
        return self.filter1_name, self.filter2_name

    @property
    def _filter1_wheel(self):
        """Returns filter wheel number for filter1 (0 or 1)."""
        return self._get_filter_wheel_number(self.filter1_name)

    @property
    def _filter2_wheel(self):
        """Returns filter wheel number for filter2 if combo; otherwise None."""
        if not self.is_combo or self.filter2_name is None:
            return None
        return self._get_filter_wheel_number(self.filter2_name)

    @staticmethod
    def _get_filter_index(filter_name: str) -> int:
        """Returns the index value for a filter name from either wheel mapping."""
        if filter_name in index_mapping_to_fw1_name:
            return index_mapping_to_fw1_name[filter_name]
        elif filter_name in index_mapping_to_fw2_name:
            return index_mapping_to_fw2_name[filter_name]
        else:
            raise ValueError(
                f"Filter '{filter_name}' not found in filter index mappings."
            )

    @property
    def filter_data_dates_to_seconds(self):
        """Returns the contamination file time axis in seconds (utime)."""
        return _filter_contamination_file_time.utime

    @property
    def filter_observation_date_to_seconds(self):
        """Returns the observation date in seconds (utime)."""
        return self.observation_date.utime

    @property
    def _filter1_data(self):
        """Returns filter contamination data for filter 1."""
        filter1_index = (
            index_mapping_to_fw1_name.get(self.filter1_name)
            if self._filter1_wheel == 0
            else index_mapping_to_fw2_name.get(self.filter1_name)
        )
        return _filter_contamination[filter1_index][self._filter1_wheel]

    @property
    def _filter2_data(self):
        """Returns filter contamination data for filter 2, if combo."""
        if not self.is_combo or self.filter2_name is None:
            return None
        filter2_index = (
            index_mapping_to_fw1_name.get(self.filter2_name)
            if self._filter2_wheel == 0
            else index_mapping_to_fw2_name.get(self.filter2_name)
        )
        return _filter_contamination[filter2_index][self._filter2_wheel]

    @property
    def contamination_on_filter1(self) -> u.angstrom:
        """
        Interpolates contamination layer thickness on filter1.
        """
        interpolater = scipy.interpolate.interp1d(
            _filter_contamination_file_time.utime, self._filter1_data, kind="linear"
        )
        return interpolater(self.observation_date.utime)

    @property
    def contamination_on_filter2(self) -> u.Quantity | None:
        """
        Interpolates the contamination thickness on filter2, if using a combo filter.

        Returns
        -------
        astropy.units.Quantity or None
            Contamination thickness in Angstroms for filter2 if present; otherwise None.
        """
        if not self.is_combo or self._filter2_data is None:
            return None

        interpolater = scipy.interpolate.interp1d(
            _filter_contamination_file_time.utime, self._filter2_data, kind="linear"
        )
        return interpolater(self.observation_date.utime)

    @cached_property
    def n_DEHP_attributes(self):
        """(Diethylhexylphthalate) Wavelength (nm), Delta, Beta."""
        _n_DEHP_filename = get_pkg_data_filename(
            "data/n_DEHP.txt", package="xrtpy.response"
        )

        with Path(_n_DEHP_filename).open() as n_DEHP:
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
    def _transmission_equation(self):
        """
        Define equations used to calculate the effective area.

        This method sets up the necessary equations to calculate the effective area of a filter
        in the X-Ray Telescope (XRT) on the Hinode satellite. The calculations are based on the
        principles of modern optics as described in G.R Fowles, "Introduction to Modern Optics,"
        2nd Edition, pp 96-101.

        Returns
        -------
        tuple
            A tuple containing the index of refraction, sine of the angle of incidence, cosine of the angle of incidence,
            maximum wavelength, index of refraction of the medium at the entrance, index of refraction of the medium at the exit,
            and the angle of incidence.
        """
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
    def _angular_wavenumber_CCD(self):
        """Define angular wavenumber on CCD."""

        index, _, cos_a, wavelength_max, _, _, _ = self._transmission_equation

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
    def _filter_contamination_angular_wavenumber(self):
        """Define angular wavenumber for the filter(s), considering contamination thickness."""
        index, _, cos_a, _, _, _, _ = self._transmission_equation

        angular_wavenumber = np.array(
            [
                (2.0 * math.pi * index[i] * cos_a) / self.n_DEHP_wavelength[i]
                for i in range(4000)
            ]
        )

        if self.is_combo:
            angular_wavenumber_thickness_f1 = (
                angular_wavenumber * self.contamination_on_filter1
            )
            angular_wavenumber_thickness_f2 = (
                angular_wavenumber * self.contamination_on_filter2
            )

            real_angular_f1 = angular_wavenumber_thickness_f1.real
            imag_angular_f1 = angular_wavenumber_thickness_f1.imag
            real_angular_f2 = angular_wavenumber_thickness_f2.real
            imag_angular_f2 = angular_wavenumber_thickness_f2.imag

            angular_wavenumber_combo = [
                complex(
                    real_angular_f1[i] + real_angular_f2[i],
                    imag_angular_f1[i] + imag_angular_f2[i],
                )
                for i in range(4000)
            ]
            return angular_wavenumber_combo
        else:
            angular_wavenumber_thickness = (
                angular_wavenumber * self.contamination_on_filter1
            )
            real_angular = angular_wavenumber_thickness.real
            imag_angular = angular_wavenumber_thickness.imag
            return [complex(real_angular[i], imag_angular[i]) for i in range(4000)]

    @cached_property
    def _CCD_contamination_transmission(self):
        """Calculate transmission matrix coefficient and transmittance on the CCD."""

        index, _, _, _, n_o, n_t, _ = self._transmission_equation

        i_i = complex(0, 1)  # Define complex number

        # Define transfer matrix
        M = [
            [
                [
                    np.cos(self._angular_wavenumber_CCD[i]),
                    (-i_i * np.sin(self._angular_wavenumber_CCD[i])) / index[i],
                ],
                [
                    -i_i * np.sin(self._angular_wavenumber_CCD[i]) * index[i],
                    np.cos(self._angular_wavenumber_CCD[i]),
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
    def wavelength(self):
        """Array of wavelengths for every X-ray channel in Angstroms (Å)."""
        _wave = self._channel.wavelength.to_value("AA")
        return _wave * u.Angstrom

    @property
    def channel_geometry_aperture_area(self):
        """XRT flight model geometry aperture area."""
        return self._channel.geometry.geometry_aperture_area

    @property
    def channel_transmission(self):
        """XRT channel transmission."""
        return np.interp(
            self.wavelength, self._channel.wavelength, self._channel.transmission
        )

    def _contamination_interpolator(self, x, y):
        return np.interp(self.wavelength.to_value("Angstrom"), x, y)

    @property
    def _interpolated_CCD_contamination_transmission(self):
        """Interpolate filter contam transmission to the wavelength."""
        return self._contamination_interpolator(
            self.n_DEHP_wavelength,
            self._CCD_contamination_transmission,
        )

    @cached_property
    def _filter_contamination_transmission(self):
        """Calculate transmission matrix coefficient and transmittance on the filter(s)."""

        index, _, _, _, n_o, n_t, _ = self._transmission_equation
        i_i = complex(0, 1)

        # angular_wave = self._filterwheel_angular_wavenumber  #handles both cases - filter 1 or 2 or
        angular_wave = self._filter_contamination_angular_wavenumber

        M = [
            [
                [
                    np.cos(angular_wave[i]),
                    (-i_i * np.sin(angular_wave[i])) / index[i],
                ],
                [
                    -i_i * np.sin(angular_wave[i]) * index[i],
                    np.cos(angular_wave[i]),
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
    def _interpolated_filter_contamination_transmission(self):
        """Interpolate filter contam transmission to the wavelength."""
        return self._contamination_interpolator(
            self.n_DEHP_wavelength,
            self._filter_contamination_transmission,
        )

    @u.quantity_input
    def effective_area(self) -> u.cm**2:
        r"""
        Calculate the Effective Area.

        The effective area is calculated by considering the geometry of the XRT flight model,
        the channel transmission, and the contamination layers on both the CCD and the filter.

        Returns
        -------
        astropy.units.Quantity
            Effective area in cm\ :math:`^2`.

        Notes
        -----
        The effective area is a crucial parameter for determining the sensitivity of the XRT
        to X-ray emissions. This method combines various factors, including the physical
        properties of the filter and CCD, to compute the total effective area.
        """
        return (
            self.channel_geometry_aperture_area
            * self.channel_transmission
            * self._interpolated_CCD_contamination_transmission
            * self._interpolated_filter_contamination_transmission
        )
