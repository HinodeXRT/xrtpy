__all__ = [
    "Geometry",
    "EntranceFilter",
    "Mirror",
    "Filter",
    "CCD",
    "Channel",
    "resolve_filter_name",
]

import numpy as np
import pkg_resources
import sunpy.io.special

from astropy import units as u

filename = pkg_resources.resource_filename(
    "xrtpy", "data/channels/xrt_channels_v0016.genx"
)

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
    name = name.replace("_", "-")
    parts: list = name.split("/")
    new_parts: list = [part.capitalize() for part in name.split("/")]
    name: str = "/".join(new_parts)
    return name


class Geometry:
    """The physical geometric parameters of Hinode/XRT."""

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._geom_data = self._genx_file[self._channel_index]["GEOM"]

    @property
    def name(
        self,
    ) -> str:
        """Hinode/XRT flight model geometry."""
        return self._geom_data["LONG_NAME"]

    @property
    @u.quantity_input
    def focal_len(
        self,
    ) -> u.cm:
        """Hinode/XRT flight model geometry focal length."""
        return u.Quantity(self._geom_data["FOC_LEN"], u.cm)

    @property
    @u.quantity_input
    def aperture_area(self) -> u.cm ** 2:
        """Hinode/XRT flight model geometry aperture area."""
        return u.Quantity(self._geom_data["APERTURE_AREA"], u.cm ** 2)


class EntranceFilter:
    """
    Entrance filter properties.

    Thin prefilters cover the narrow annular entrance aperture of the XRT serving two main purposes:

    1. Reduce the visible light entering the instrument.
    2. Reduce the heat load in the instrument.
    """

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._en_filter_data = self._genx_file[self._channel_index]["EN_FILTER"]

    @property
    def entrancefilter_name(
        self,
    ) -> str:
        """Entrance filter name."""
        return self._en_filter_data["LONG_NAME"]

    @property
    @u.quantity_input
    def entrancefilter_thickness(self) -> u.angstrom:
        """XRT entrance filter material thickness."""
        return u.Quantity(self._en_filter_data["THICK"], u.angstrom)

    @property
    def entrancefilter_density(self) -> u.g * u.cm ** -3:
        """XRT entrance filter material density in g cm\ :math:`^{-3}`."""
        return u.Quantity(self._en_filter_data["DENS"], u.g * u.cm ** -3)

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._en_filter_data["LENGTH"]

    @property
    @u.quantity_input
    def entrancefilter_wavelength(self) -> u.angstrom:
        """Array of wavelengths for entrance filter transmission  in angstroms."""
        return u.Quantity(self._en_filter_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def entrancefilter_transmission(self):
        """Entrance filter transmission."""
        return self._en_filter_data["TRANS"][: self.number_of_wavelengths]

    @property
    def entrancefilter_mesh_transmission(self):
        """Transmission of mesh filter substrate."""
        return self._en_filter_data["MESH_TRANS"]

    @property
    def entrancefilter_substrate(self) -> str:
        """XRT entrance filter substrate."""
        return self._en_filter_data["SUBSTRATE"]

    @property
    def entrancefilter_material(self):
        """Filter material on entrance filter."""
        return self._en_filter_data["MATERIAL"]


class Mirror:
    """
    Grazing incidence mirror properties.

    Grazing-incidence optics used for soft X-ray imaging generally require a minimum of two surfaces.
    Since XRT is a two-bounce telescope, there are two mirror reflectivities.
    """

    _genx_file = _genx_file

    def __init__(self, index, mirror_number):
        self._channel_index = index
        self._mirror_data = self._genx_file[self._channel_index][
            f"MIRROR{mirror_number}"
        ]

    @property
    def name(
        self,
    ) -> str:
        """Hinode/XRT flight model mirror."""
        return self._mirror_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """XRT flight model mirror material."""
        return self._mirror_data["MATERIAL"]

    @property
    @u.quantity_input
    def density(self) -> u.g * u.cm ** -3:
        """Mirror mass density."""
        return u.Quantity(self._mirror_data["DENS"], u.g * u.cm ** -3)

    @property
    @u.quantity_input
    def graze_angle(
        self,
    ) -> u.deg:
        """Mirror graze angle in units of degrees."""
        return u.Quantity(self._mirror_data["GRAZE_ANGLE"], u.deg)

    @property
    @u.quantity_input
    def wavelength(
        self,
    ) -> u.angstrom:
        """Array of wavelengths for mirror reflectance."""
        return u.Quantity(self._mirror_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    @u.quantity_input
    def reflection(
        self,
    ) -> u.angstrom:
        """Reflection of a mirror."""
        return u.Quantity(self._mirror_data["REFL"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._mirror_data["LENGTH"]


class Filter:
    """
    X-ray filters using two filter wheels.

    The corresponding categories are used for both filter 1 and filter 2.
    """

    _genx_file = _genx_file

    def __init__(self, index, filter_number):
        self._channel_index = index
        self._fp_filter_data = self._genx_file[self._channel_index][
            f"FP_FILTER{filter_number}"
        ]

    @property
    def name(
        self,
    ) -> str:
        """XRT focal plane filter position."""
        return self._fp_filter_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """XRT focal plane filter material."""
        return self._fp_filter_data["MATERIAL"]

    @property
    @u.quantity_input
    def thickness(self) -> u.angstrom:
        """XRT  focal plane filter thickness."""
        return u.Quantity(self._fp_filter_data["THICK"], u.angstrom)

    @property
    @u.quantity_input
    def density(self) -> u.g * u.cm ** -3:
        """XRT  focal plane filter density."""
        return u.Quantity(self._fp_filter_data["DENS"], u.g * u.cm ** -3)

    @property
    @u.quantity_input
    def wavelength(self) -> u.angstrom:
        """Array of wavelength for every X-ray focal plane filter in angstroms."""
        return u.Quantity(self._fp_filter_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def transmission(self):
        """Get focal plane filter transmission."""
        return self._fp_filter_data["TRANS"][: self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._fp_filter_data["LENGTH"]

    @property
    def mesh_trans(self):
        """Mesh transmission for the focal plane filter."""
        return self._fp_filter_data["MESH_TRANS"]

    @property
    def substrate(self) -> str:
        """XRT substrate for the focal plane filter."""
        return self._fp_filter_data["SUBSTRATE"]

    @property
    def thickness(self):
        """Filter thickness."""
        return self._fp_filter_data["THICK"]


class CCD:
    """Charge-coupled device on board XRT."""

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._ccd_data = self._genx_file[self._channel_index]["CCD"]

    @property
    def ccd_name(
        self,
    ) -> str:
        """Hinode/XRT flight model CCD."""
        return self._ccd_data["LONG_NAME"]

    @property
    @u.quantity_input
    def ccd_ev_ore_electron(
        self,
    ) -> u.eV / u.electron:
        """The energy necessary to dislodge one electron."""
        return u.Quantity(self._ccd_data["EV_PER_EL"], u.eV / u.electron)

    @property
    @u.quantity_input
    def ccd_full_well(
        self,
    ) -> u.electron:
        """Number of electrons for a CCD full well."""
        return u.Quantity(self._ccd_data["FULL_WELL"], u.electron)

    @property
    @u.quantity_input
    def ccd_gain_left(
        self,
    ) -> u.electron:
        """Gain when reading the left port of the CCD."""
        return u.Quantity(self._ccd_data["GAIN_L"], u.electron)

    @property
    @u.quantity_input
    def ccd_gain_right(
        self,
    ) -> u.electron:
        """Gain when reading the right port of the CCD."""
        return u.Quantity(self._ccd_data["GAIN_R"], u.electron)

    @property
    def ccd_quantum_efficiency(self):
        """Quantum efficiency of the CCD."""
        return self._ccd_data["QE"][: self.number_of_wavelengths]

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
    @u.quantity_input
    def ccd_wavelength(self) -> u.angstrom:
        """Array of wavelengths for the CCD quantum efficiency in angstroms."""
        return u.Quantity(self._ccd_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]


class Channel:
    """
    XRTpy

    Available channels: ``"Al-mesh"``, ``"Al-poly"``,  ``"C-poly"``, ``"Ti-poly"``, ``"Be-thin"``, ``"Be-med"``, ``"Al-med"``, ``"Al-thick"``,  ``"Be-thick"`` ,
    ``"Al-poly/Al-mesh"``, ``"Al-poly/Ti-poly"``, ``"Al-poly/Al-thick"``, ``"Al-poly/Be-thick"`` , ``"C-poly/Ti-poly"``
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
        elif name.lower() == "open":  # Complete by adding remaining indexs
            self._sample_channel_data = _genx_file[1]
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
        return self._geometry

    @property
    def entrancefilter(self) -> EntranceFilter:
        return self._entrancefilter

    @property
    def mirror_1(self) -> Mirror:
        return self._mirror_1

    @property
    def mirror_2(self) -> Mirror:
        return self._mirror_2

    @property
    def filter_1(self) -> Filter:
        return self._filter_1

    @property
    def filter_2(self) -> Filter:
        return self._filter_2

    @property
    def ccd(self) -> CCD:
        return self._ccd

    def __str__(self):
        """Reable printout."""
        return f"XRT Channel for {self.name}"

    def __repr__(self):
        """Code representation."""
        return f"Channel({repr(self.name)})"

    @property
    def name(self) -> str:
        """Name of XRT X-Ray channel."""
        return self._channel_data["NAME"]

    @property
    @u.quantity_input
    def wavelength(self) -> u.angstrom:
        """Array of wavelengths for every X-ray channel in angstroms."""
        return u.Quantity(self._channel_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def transmission(self):
        """Get channel transmission."""
        return self._channel_data["TRANS"][: self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._channel_data["LENGTH"]

    @property
    def observatory(self) -> str:
        """Spacecraft: Hinode satellite."""
        return self._channel_data["OBSERVATORY"]

    @property
    def instrument(self) -> str:
        """X-Ray Telescope -XRT."""
        return self._channel_data["INSTRUMENT"]
