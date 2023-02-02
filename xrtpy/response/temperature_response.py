__all__ = [
    "TemperatureResponseFundamental",
]

import numpy as np
import scipy.io
import sunpy.time

from astropy import units as u
from astropy.constants import c, h
from datetime import datetime
from pathlib import Path
from scipy import integrate, interpolate

from xrtpy.response.channel import Channel, resolve_filter_name
from xrtpy.response.effective_area import effective_area
from xrtpy.util.time import epoch

_c_Å_per_s = c.to(u.angstrom / u.second).value
_h_eV_s = h.to(u.eV * u.s).value


_CHIANTI_filename = (
    Path(__file__).parent.absolute() / "data" / "XRT_emiss_model.default_CHIANTI.geny"
)
_CHIANTI_file = scipy.io.readsav(_CHIANTI_filename)
_XRT_emiss_model_file = _CHIANTI_file["p0"]


CHIANTI_file = {
    "abundance_model": _XRT_emiss_model_file["ABUND_MODEL"][0],
    "dens_model": _XRT_emiss_model_file["DENS_MODEL"][0],
    "ioneq_model": _XRT_emiss_model_file["IONEQ_MODEL"][0],
    "name": _XRT_emiss_model_file["NAME"][0],
    "spectra": _XRT_emiss_model_file["SPEC"],
    "spectra_units": _XRT_emiss_model_file["SPEC_UNITS"][0],
    "temperature": _XRT_emiss_model_file["TEMP"][0],
    "temp_units": _XRT_emiss_model_file["TEMP_UNITS"][0],
    "tlength": _XRT_emiss_model_file["TLENGTH"][0],
    "wlength": _XRT_emiss_model_file["WLENGTH"][0],
    "wavelength": _XRT_emiss_model_file["WAVE"][0],
    "wavelength_units": _XRT_emiss_model_file["WAVE_UNITS"][0],
}

_list_of_abundance = {
    "chianti": _XRT_emiss_model_file,
    "coronal": _XRT_coronal_chianti_emiss_model,
    "hybrid": 2,
    "photospheric": 3,
}

coronal_CHIANTI_file = {
    "logged_temperature": _XRT_coronal_chianti_emiss_model["LOGTE"],
    "wavelength": _XRT_coronal_chianti_emiss_model["LMBDA"],
    "corona_solar_spectra": _XRT_coronal_chianti_emiss_model["SOLSPEC"],
    "spectra_generated_information": _XRT_coronal_chianti_emiss_model["HEADER"]["TEXT"],
}

_corona_CHIANTI_filename = (
    Path(__file__).parent.absolute() / "data" / "solspec_ch1000_corona_chianti.genx"
)

_XRT_coronal_chianti_emiss_model = sunpy.io.special.genx.read_genx(
    _corona_CHIANTI_filename
)


def resolve_abundance_model_type(abundance_type):
    """Formats users abundance_model name."""
    if not isinstance(abundance_type, str):
        raise TypeError("Abundance model name must be a string")
    abundance_name = abundance_type.lower()
    list_of_abundance_name = ["coronal", "hybrid", "photospheric"]
    if abundance_name not in list_of_abundance_name:
        raise ValueError(
            f"\n{abundance_name} is not a current model for XRTpy.\n"
            "Available abundance models:\n"
            "Coronal, Hybrid and Photospheric.\n"
        )
    return abundance_name


class TemperatureResponseFundamental:
    """Produce the temperature response for each XRT x-ray channel, assuming a spectral emission model."""

    def __init__(self, filter_name, observation_date, abundance_type):
        self._name = resolve_filter_name(filter_name)
        self.observation_date = observation_date
        self._channel = Channel(self.name)
        self._abundance_type = resolve_abundance_model_type(abundance_type)

    @property
    def abundance_type(self):
        """Name of abundance model."""
        return self._abundance_type

    @property
    def name(self):
        """Name of searched filter."""
        return self._name

    @property
    def observation_date(self):
        """Users date of observation."""
        return self._observation_date

    @observation_date.setter
    def observation_date(self, date):
        """Validating users requested observation date."""
        astropy_time = sunpy.time.parse_time(date)  # Astropy time in utc
        observation_date = astropy_time.datetime
        if observation_date <= epoch:
            raise ValueError(
                rf"Invalid date: {observation_date}.\n Date must be after September 22nd, 2006 21:36:00."
            )
        self._observation_date = observation_date

    @property
    def CHIANTI_version(self):
        """Name of the emission model."""
        return CHIANTI_file["name"]

    @property
    def abundance_model(self):
        """A brief description of what abundance model was used in the creation of the emission spectra."""
        return CHIANTI_file["abundance_model"]

    @property
    def density_model(self):
        """A brief description of the plasma density, emission measure,or differential emission measure that was used in the creation of the emission spectra."""
        return CHIANTI_file["dens_model"]

    @property
    def ionization_model(self):
        """A brief description of that ionization equilibrium model was used in the creation of the emission spectra."""
        return CHIANTI_file["ioneq_model"]

    @property
    @u.quantity_input
    def CHIANTI_temperature(self):
        """CHIANTI temperatures in kelvin."""
        return u.Quantity(CHIANTI_file["temperature"] * u.K)

    @property
    def CHIANTI_file_spectra(self):
        """CHIANTI file spectra."""
        return CHIANTI_file["spectra"]

    @property
    @u.quantity_input
    def CHIANTI_wavelength(self):
        """CHIANTI file wavelength values in Å."""
        return u.Quantity(CHIANTI_file["wavelength"] * u.Angstrom)

    @property
    @u.quantity_input
    def channel_wavelength(self):
        """Array of wavelengths for every X-ray channel in Å."""
        return u.Quantity((Channel(self.name).wavelength[:3993]) * u.photon)

    @property
    def focal_len(self):
        """Focal length of the telescope in units of cm."""
        return Channel(self.name).geometry.geometry_focal_len

    @property
    def ev_per_electron(self):
        """Amount of energy it takes to dislodge 1 electron in the CCD."""
        return Channel(self.name).ccd.ccd_energy_per_electron

    @property
    @u.quantity_input
    def pixel_size(self) -> u.cm:
        """CCD pixel size. Units converted from μm to cm."""
        ccd_pixel_size = Channel(self.name).ccd.ccd_pixel_size
        return ccd_pixel_size.to(u.cm)

    @property
    @u.quantity_input
    def solid_angle_per_pixel(self) -> u.sr / u.pix:
        """This quantity represents the solid angle, which is given in units of steradians over pixel."""
        return (self.pixel_size / self.focal_len) ** 2 * (u.sr / u.pix)

    @u.quantity_input
    def spectra(self) -> u.photon * u.cm**3 / (u.sr * u.s * u.Angstrom):
        """Interpolation between the spectra wavelength onto the channel wavelength."""
        spectra_interpolate = []
        for i in range(61):
            interpolater = interpolate.interp1d(
                self.CHIANTI_wavelength, CHIANTI_file["spectra"][0][i], kind="linear"
            )
            spectra_interpolate.append(interpolater(self.channel_wavelength))
        return spectra_interpolate * (
            u.photon * u.cm**3 * (1 / u.sr) * (1 / u.s) * (1 / u.Angstrom)
        )

    @u.quantity_input
    def effective_area(self) -> u.cm**2:
        return effective_area(self.name, self.observation_date)

    @u.quantity_input
    def integration(self) -> u.electron * u.cm**5 / (u.s * u.pix):
        wavelength = (self.channel_wavelength).value
        constants = (_c_Å_per_s * _h_eV_s / self.channel_wavelength).value
        factors = (self.solid_angle_per_pixel / self.ev_per_electron).value
        effective_area = (self.effective_area()).value
        dwvl = wavelength[1:] - wavelength[:-1]
        dwvl = np.append(dwvl, dwvl[-1])
        # Simple summing like this is appropriate for binned data like in the current
        # spectrum file. More recent versions of Chianti include the line width,
        # which then makes the previous version that uses Simpson's method
        # to integrate more appropriate (10/05/2022)
        temp_resp_w_u_c = (
            self.spectra().value * effective_area * constants * factors * dwvl
        ).sum(axis=1)

        return temp_resp_w_u_c * (u.electron * u.cm**5 * (1 / u.s) * (1 / u.pix))

    @property
    @u.quantity_input
    def ccd_gain_right(self) -> u.electron / u.DN:
        """Provide the camera gain in electrons per data number."""
        return Channel(self.name).ccd.ccd_gain_right

    @u.quantity_input
    def temperature_response(self) -> u.DN * u.cm**5 / (u.s * u.pix):
        r"""Apply gain value to the Temperature Response in units of DN cm\ :sup:`5` s\ :sup:`-1` pix\ :sup:`-1`."""
        return self.integration() / self.ccd_gain_right
