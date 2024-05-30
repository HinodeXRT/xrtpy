"""Classes for describing channels on Hinode/XRT."""

__all__ = [
    "Geometry",
    "EntranceFilter",
    "Mirror",
    "Filter",
    "CCD",
    "Channel",
    "resolve_filter_name",
]

from pathlib import Path

import numpy as np
import sunpy.io.special
import sunpy.time
from astropy import units as u

filename = Path(__file__).parent.absolute() / "data" / "xrt_channels_v0016.genx"

_channel_name_to_index_mapping = {
    "Al-mesh": 0,
    "Al-poly": 1,
    "C-poly": 2,
    "Ti-poly": 3,
    "Be-thin": 4,
    "Be-med": 5,
    "Al-med": 6,
    "Al-thick": 7,
    "Be-thick": 8,
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10,
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12,
    "C-poly/Ti-poly": 13,
}


_genx_file = sunpy.io.special.genx.read_genx(filename)["SAVEGEN0"]


def resolve_filter_name(name):
    """
    Formats the user's filter name to match the expected format.

    Parameters
    ----------
    name : str
        The filter name provided by the user.

    Returns
    -------
    str
        The formatted filter name.

    Raises
    ------
    TypeError
        If the provided name is not a string.
    """
    if not isinstance(name, str):
        raise TypeError("name must be a string")
    name = name.replace("_", "-")
    parts: list = name.split("/")
    new_parts: list = [part.capitalize() for part in parts]
    name: str = "/".join(new_parts)
    return name


class Geometry:
    """
    The physical geometric parameters of the X-Ray Telescope (XRT) on board the Hinode spacecraft.

    Parameters
    ----------
    index : int
        The index of the channel in the GENX file.

    Attributes
    ----------
        Channel.Geometry.geometry_name : str
            Hinode/XRT flight model geometry.
        Channel.Geometry.geometry_focal_len : astropy.units.Quantity
            XRT flight model geometry focal length in cm.
        Channel.Geometry.geometry_aperture_area : astropy.units.Quantity
            XRT flight model geometry aperture area in cm^2.
    """

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._geom_data = self._genx_file[self._channel_index]["GEOM"]

    @property
    def geometry_name(
        self,
    ) -> str:
        """Hinode/XRT flight model geometry."""
        return self._geom_data["LONG_NAME"]

    @property
    @u.quantity_input
    def geometry_focal_len(self) -> u.cm:
        """XRT flight model geometry focal length."""
        return u.Quantity(self._geom_data["FOC_LEN"], u.cm)

    @property
    @u.quantity_input
    def geometry_aperture_area(self) -> u.cm**2:
        """XRT flight model geometry aperture area."""
        return u.Quantity(self._geom_data["APERTURE_AREA"], u.cm**2)


class EntranceFilter:
    """
    Represents the entrance filter of the X-Ray Telescope (XRT) on the Hinode spacecraft.

    The entrance filter covers the annular entrance aperture of the XRT, serving to reduce
    both visible light and heat load entering the instrument. These filters are critical
    for maintaining the necessary optical conditions inside the telescope.

    Parameters
    ----------
    index : int
        The index of the channel within the GENX file that contains data specific to this filter.

    Attributes
    ----------
        Channel.EntranceFilter.entrancefilter_density : astropy.units.Quantity
            The density of the entrance filter material in g/cm³.
        Channel.EntranceFilter.entrancefilter_material : str
            The material composition of the entrance filter.
        Channel.EntranceFilter.entrancefilter_mesh_transmission : float
            The percentage transmission of the mesh part of the filter.
        Channel.EntranceFilter.entrancefilter_name : str
            The descriptive name of the entrance filter.
        Channel.EntranceFilter.entrancefilter_substrate : str
            The substrate material of the entrance filter.
        Channel.EntranceFilter.entrancefilter_thickness : astropy.units.Quantity
            The thickness of the entrance filter material measured in Angstroms.
        Channel.EntranceFilter.entrancefilter_transmission : numpy.ndarray
            The transmission efficiency of the entrance filter across different wavelengths.
        Channel.EntranceFilter.entrancefilter_wavelength : astropy.units.Quantity
            The wavelengths at which the transmission data of the filter are measured, in Angstroms.
    """

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._en_filter_data = self._genx_file[self._channel_index]["EN_FILTER"]

    @property
    def entrancefilter_density(self) -> u.g * u.cm**-3:
        r"""XRT entrance filter material density in g/cm\ :sup:`3`\ ."""
        return u.Quantity(self._en_filter_data["DENS"], u.g * u.cm**-3)

    @property
    def entrancefilter_material(self) -> str:
        """XRT entrance filter material."""
        return self._en_filter_data["MATERIAL"]

    @property
    def entrancefilter_mesh_transmission(self):
        """Transmission of mesh filter substrate."""
        return self._en_filter_data["MESH_TRANS"]

    @property
    def entrancefilter_name(self) -> str:
        """Entrance filter name."""
        return self._en_filter_data["LONG_NAME"]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._en_filter_data["LENGTH"]

    @property
    def entrancefilter_substrate(self) -> str:
        """XRT entrance filter substrate."""
        return self._en_filter_data["SUBSTRATE"]

    @property
    @u.quantity_input
    def entrancefilter_wavelength(self) -> u.angstrom:
        """Array of wavelengths for entrance filter transmission in angstroms."""
        return u.Quantity(self._en_filter_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    @u.quantity_input
    def entrancefilter_thickness(self) -> u.angstrom:
        """XRT entrance filter material thickness in angstroms."""
        return u.Quantity(self._en_filter_data["THICK"], u.angstrom)

    @property
    def entrancefilter_transmission(self):
        """Entrance filter transmission."""
        return self._en_filter_data["TRANS"][: self.number_of_wavelengths]


class Mirror:
    """
    Defines a grazing incidence mirror used in the X-Ray Telescope (XRT) for imaging in soft X-rays.

    XRT utilizes a two-bounce telescope design requiring at least two mirrors. This class provides
    properties of one of these mirrors.

    Parameters
    ----------
    index : int
        Index of the channel within the GENX data file.
    mirror_number : int
        Specifies whether this is the first or second mirror (1 or 2).

    Attributes
    ----------
        Channel.Mirror.mirror_density : astropy.units.Quantity
            The mass density of the mirror in g/cm³.
        Channel.Mirror.mirror_graze_angle : astropy.units.Quantity
            The graze angle of the mirror during operation, measured in degrees.
        Channel.Mirror.mirror_material : str
            The material composition of the mirror.
        Channel.Mirror.mirror_name : str
            The name or identifier of the mirror.
        Channel.Mirror.mirror_reflection : numpy.ndarray
            The reflectance values of the mirror at different wavelengths.
        Channel.Mirror.mirror_wavelength : astropy.units.Quantity
            The wavelengths at which the mirror's reflectance is measured, in Angstroms.
    """

    _genx_file = _genx_file

    def __init__(self, index, mirror_number):
        self._channel_index = index
        self._mirror_data = self._genx_file[self._channel_index][
            f"MIRROR{mirror_number}"
        ]

    @property
    @u.quantity_input
    def mirror_density(self) -> u.g * u.cm**-3:
        """Mirror mass density."""
        return u.Quantity(self._mirror_data["DENS"], u.g * u.cm**-3)

    @property
    @u.quantity_input
    def mirror_graze_angle(self) -> u.deg:
        """Mirror graze angle in units of degrees."""
        return u.Quantity(self._mirror_data["GRAZE_ANGLE"], u.deg)

    @property
    def mirror_name(self) -> str:
        """Hinode/XRT flight model mirror."""
        return self._mirror_data["LONG_NAME"]

    @property
    def mirror_material(self) -> str:
        """XRT flight model mirror material."""
        return self._mirror_data["MATERIAL"]

    @property
    @u.quantity_input
    def mirror_reflection(self) -> u.angstrom:
        """Reflection of a mirror."""
        return u.Quantity(self._mirror_data["REFL"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    @u.quantity_input
    def mirror_wavelength(self) -> u.angstrom:
        """Array of wavelengths for mirror reflectance."""
        return u.Quantity(self._mirror_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._mirror_data["LENGTH"]


class Filter:
    """
    Represents one of the focal plane filters of the X-Ray Telescope (XRT) on the Hinode spacecraft,
    which are mounted on two filter wheels.

    Parameters
    ----------
    index : int
        Index of the channel within the GENX data file.
    filter_number : int
        Specifies the filter wheel position (1 or 2).

    Attributes
    ----------
        Channel.Filter.filter_density : astropy.units.Quantity
            The density of the filter material in g/cm³.
        Channel.Filter.filter_material : str
            The material composition of the filter.
        Channel.Filter.filter_mesh_transmission : float
            The transmission efficiency of the filter's mesh.
        Channel.Filter.filter_name : str
            The descriptive name of the filter.
        Channel.Filter.filter_substrate : str
            The substrate material of the filter.
        Channel.Filter.filter_thickness : astropy.units.Quantity
            The thickness of the filter material, measured in Angstroms.
        Channel.Filter.filter_transmission : numpy.ndarray
            The transmission efficiency of the filter across different wavelengths.
        Channel.Filter.filter_wavelength : astropy.units.Quantity
            The wavelengths at which the filter's transmission data are measured, in Angstroms.
    """

    _genx_file = _genx_file

    def __init__(self, index, filter_number):
        self._channel_index = index
        self._fp_filter_data = self._genx_file[self._channel_index][
            f"FP_FILTER{filter_number}"
        ]

    @property
    @u.quantity_input
    def filter_density(self) -> u.g * u.cm**-3:
        """XRT filter density."""
        return u.Quantity(self._fp_filter_data["DENS"], u.g * u.cm**-3)

    @property
    def filter_material(self) -> str:
        """XRT filter material."""
        return self._fp_filter_data["MATERIAL"]

    @property
    def filter_mesh_transmission(self):
        """Mesh transmission for the focal plane filter."""
        return self._fp_filter_data["MESH_TRANS"]

    @property
    def filter_name(self) -> str:
        """XRT filter name."""
        return self._fp_filter_data["LONG_NAME"]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._fp_filter_data["LENGTH"]

    @property
    def filter_substrate(self) -> str:
        """XRT filter substrate."""
        return self._fp_filter_data["SUBSTRATE"]

    @property
    @u.quantity_input
    def filter_thickness(self) -> u.angstrom:
        """Filter thickness."""
        return u.Quantity(self._fp_filter_data["THICK"], u.angstrom)

    @property
    def filter_transmission(self):
        """Filter transmission."""
        return self._fp_filter_data["TRANS"][: self.number_of_wavelengths]

    @property
    @u.quantity_input
    def filter_wavelength(self) -> u.angstrom:
        """XRT filter wavelength in angstroms."""
        return u.Quantity(self._fp_filter_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]


class CCD:
    """
    Describes the Charge-Coupled Device (CCD) used in the X-Ray Telescope (XRT) on the Hinode spacecraft.

    The CCD is a crucial component of the XRT, responsible for capturing X-ray images. This class provides
    various properties and characteristics of the CCD, including its gain, quantum efficiency, and physical
    dimensions.

    Parameters
    ----------
    index : int
        Index of the channel within the GENX data file containing CCD-specific properties.

    Attributes
    ----------
        Channel.CCD.ccd_energy_per_electron : astropy.units.Quantity
            The energy required to dislodge a single electron, measured in electron-volts per electron.
        Channel.CCD.ccd_full_well : astropy.units.Quantity
            The full well capacity of the CCD in terms of the maximum number of electrons it can hold.
        Channel.CCD.ccd_gain_left : astropy.units.Quantity
            The gain of the CCD when reading from the left port, measured in electrons per digital number.
        Channel.CCD.ccd_gain_right : astropy.units.Quantity
            The gain of the CCD when reading from the right port, measured in electrons per digital number.
        Channel.CCD.ccd_name : str
            The name or identifier of the CCD.
        Channel.CCD.ccd_pixel_size : astropy.units.Quantity
            The size of individual pixels on the CCD, measured in micrometers.
        Channel.CCD.ccd_quantum_efficiency : numpy.ndarray
            The quantum efficiency of the CCD, representing the efficiency of photon-to-electron conversion.
        Channel.CCD.ccd_wavelength : astropy.units.Quantity
            The wavelengths at which the CCD's quantum efficiency is measured, in Angstroms.
    """

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._ccd_data = self._genx_file[self._channel_index]["CCD"]

    @property
    @u.quantity_input
    def ccd_energy_per_electron(self) -> u.eV / u.electron:
        """The energy necessary to dislodge one electron."""
        return u.Quantity(self._ccd_data["EV_PER_EL"], u.eV / u.electron)

    @property
    @u.quantity_input
    def ccd_full_well(self) -> u.electron:
        """Number of electrons for a CCD full well."""
        return u.Quantity(self._ccd_data["FULL_WELL"], u.electron)

    @property
    @u.quantity_input
    def ccd_gain_left(self) -> u.electron / u.DN:
        """Gain when reading the left port of the CCD."""
        return u.Quantity(self._ccd_data["GAIN_L"], u.electron / u.DN)

    @property
    @u.quantity_input
    def ccd_gain_right(self) -> u.electron / u.DN:
        """Gain when reading the right port of the CCD."""
        return u.Quantity(57.5, u.electron / u.DN)

    @property
    def ccd_name(self) -> str:
        """Hinode/XRT flight model CCD."""
        return self._ccd_data["LONG_NAME"]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._ccd_data["LENGTH"]

    @property
    @u.quantity_input
    def ccd_pixel_size(self) -> u.micron:
        """CCD pixel size in micrometers."""
        return u.Quantity(self._ccd_data["PIXEL_SIZE"], u.micron)

    @property
    def ccd_quantum_efficiency(self):
        """Quantum efficiency of the CCD."""
        return self._ccd_data["QE"][: self.number_of_wavelengths]

    @property
    @u.quantity_input
    def ccd_wavelength(self) -> u.angstrom:
        """Array of wavelengths for the CCD quantum efficiency in angstroms."""
        return u.Quantity(self._ccd_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]


class Channel:
    """
    Represents an XRT channel on the Hinode spacecraft.

    Available channels: "Al-mesh", "Al-poly", "C-poly", "Ti-poly", "Be-thin", "Be-med",
    "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-mesh", "Al-poly/Ti-poly", "Al-poly/Al-thick",
    "Al-poly/Be-thick", "C-poly/Ti-poly".

    Parameters
    ----------
    name : str
        The name of the filter for the XRT channel.

    Attributes
    ----------
    Geometry
        The geometric parameters of the XRT channel.
    EntranceFilter
        The entrance filter properties.
    Mirror 1
        Properties of the first mirror.
    Mirror 2
        Properties of the second mirror.
    Filter 1
        Properties of the first filter.
    Filter 2
        Properties of the second filter.
    CCD
        Properties of the CCD.
    Name
        Name of XRT X-Ray channel. Type: str
    Wavelength
        Array of wavelengths for every X-ray channel in angstroms. Type: astropy.units.Quantity
    Transmission
        Transmission of the channel.  Type: numpy.ndarray
    Number_of_wavelengths
        Length of the data.  Type: int
    Observatory
        Name of the spacecraft.  Type: str
    Instrument
        Name of the instrument (X-Ray Telescope -XRT).  Type: str

    """

    _genx_file = _genx_file

    def __init__(self, name):
        name = resolve_filter_name(name)
        if name in _channel_name_to_index_mapping:
            self._channel_index = _channel_name_to_index_mapping[name]
            self._channel_data = _genx_file[self._channel_index]
            self._geometry = Geometry(self._channel_index)
            self._entrancefilter = EntranceFilter(self._channel_index)
            self._mirror_1 = Mirror(self._channel_index, 1)
            self._mirror_2 = Mirror(self._channel_index, 2)
            self._filter_1 = Filter(self._channel_index, 1)
            self._filter_2 = Filter(self._channel_index, 2)
            self._ccd = CCD(self._channel_index)
        elif name.lower() == "open":  # Complete by adding remaining indices
            self._sample_channel_data = _genx_file[1]
            self._geometry = Geometry(1)
            self._channel_data = {
                "WAVE": self._sample_channel_data["WAVE"],
                "TRANS": np.ones_like(self._sample_channel_data["TRANS"]),
                "LENGTH": self._sample_channel_data["LENGTH"],
            }
        else:
            raise ValueError(
                f"{name} is not a valid channel. The available channels are: {list(_channel_name_to_index_mapping.keys())}"
            )

    @property
    def geometry(self) -> Geometry:
        """
        Geometric parameters of the XRT channel.
        """
        return self._geometry

    @property
    def entrancefilter(self) -> EntranceFilter:
        """
        Entrance filter properties.
        """
        return self._entrancefilter

    @property
    def mirror_1(self) -> Mirror:
        """
        Properties of the first mirror.
        """
        return self._mirror_1

    @property
    def mirror_2(self) -> Mirror:
        """
        Properties of the second mirror.
        """
        return self._mirror_2

    @property
    def filter_1(self) -> Filter:
        """
        Properties of the first filter.
        """
        return self._filter_1

    @property
    def filter_2(self) -> Filter:
        """
        Properties of the second filter.
        """
        return self._filter_2

    @property
    def ccd(self) -> CCD:
        """
        Properties of the CCD.
        """
        return self._ccd

    def __str__(self):
        """Readable printout."""
        return f"XRT Channel for {self.name}"

    def __repr__(self):
        """Code representation."""
        return f"Channel({self.name!r})"

    @property
    def name(self) -> str:
        """
        Name of XRT X-Ray channel.
        """
        return self._channel_data["NAME"]

    @property
    @u.quantity_input
    def wavelength(self) -> u.angstrom:
        """
        Array of wavelengths for every X-ray channel in angstroms.
        """
        return u.Quantity(self._channel_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def transmission(self):
        """
        Channel transmission.
        """
        return self._channel_data["TRANS"][: self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._channel_data["LENGTH"]

    @property
    def observatory(self) -> str:
        """
        The spacecraft name - Hinode.
        """
        return self._channel_data["OBSERVATORY"]

    @property
    def instrument(self) -> str:
        """
        Instrument - X-Ray Telescope -XRT.
        """
        return self._channel_data["INSTRUMENT"]
