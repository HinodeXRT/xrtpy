__all__ = [
    "TemperatureResponseFundamental",
]

import pkg_resources
import scipy.io
import sunpy.time

from astropy import units as u
from astropy.constants import c, h
from datetime import datetime
from scipy import integrate, interpolate

from xrtpy.response.channel import Channel, resolve_filter_name
from xrtpy.response.effective_area import effective_area
from xrtpy.util.time import epoch

_c_Å_per_s = c.to(u.angstrom / u.second).value
_h_eV_s = h.to(u.eV * u.s).value

_CHIANTI_filename = pkg_resources.resource_filename(
    "xrtpy", "response/data/XRT_emiss_model.default_CHIANTI.geny"
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


class TemperatureResponseFundamental:
    """Produce the temperature response for each XRT x-ray channel, assuming a spectral emission model."""

    def __init__(self, filter_name, observation_date):
        self._name = resolve_filter_name(filter_name)
        self.observation_date = observation_date
        self._channel = Channel(self.name)

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
                f"Invalid date: {observation_date}.\n Date must be after September 22nd, 2006 21:36:00."
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
        return Channel(self.name).ccd.ccd_ev_pre_electron

    @property
    @u.quantity_input
    def pixel_size(self) -> u.cm:
        """CCD pixel size. Units converted from microns to cm."""
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
        for i in range(0, 61):
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

        temp_resp_w_u_c = []
        for i in range(0, 61):
            temp_resp_w_u_c.append(
                integrate.simpson(
                    self.spectra()[i] * effective_area * constants * factors, wavelength
                )
            )
        return temp_resp_w_u_c * (u.electron * u.cm**5 * (1 / u.s) * (1 / u.pix))

    @property
    @u.quantity_input
    def ccd_gain_right(self) -> 1 / u.DN:
        """Provide the camera gain in electrons per data number."""
        return Channel(self.name).ccd.ccd_gain_right / u.DN

    @u.quantity_input
    def temperature_response(self) -> u.DN * u.cm**5 / (u.s * u.pix):
        """Apply gain value to the Temperature Response in units of DN cm\ :sup:`5` s\ :sup:`-1` pix\ :sup:`-1`."""
        return self.integration() / self.ccd_gain_right
