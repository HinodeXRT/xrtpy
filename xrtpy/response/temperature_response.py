__all__ = [
    "TemperatureResponseFundamental",
]

import numpy as np
import scipy.io
import sunpy.time

from astropy import units as u
from astropy.constants import c, h
from datetime import datetime
from numbers import Real
from pathlib import Path
from scipy import integrate, interpolate
from typing import Dict

from xrtpy.response.channel import Channel, resolve_filter_name
from xrtpy.response.effective_area import effective_area
from xrtpy.util.time import epoch

_c_Å_per_s = c.to(u.angstrom / u.second).value
_h_eV_s = h.to(u.eV * u.s).value


_abundance_model_file_path = {
    "coronal_abundance_path": Path(__file__).parent.absolute()
    / "data/chianti_emission_models"
    / "XRT_emiss_model.default_CHIANTI.geny",
    "hybrid_abundance_path": Path(__file__).parent.absolute()
    / "data/chianti_emission_models"
    / "XRT_emiss_model.default_CHIANTI_photospheric.geny",
    "photospheric_abundance_path": Path(__file__).parent.absolute()
    / "data/chianti_emission_models"
    / "XRT_emiss_model.default_CHIANTI_hybrid.geny",
}

_abundance_model_data = {
    "coronal": scipy.io.readsav(_abundance_model_file_path["coronal_abundance_path"])[
        "p0"
    ],
    "hybrid": scipy.io.readsav(_abundance_model_file_path["hybrid_abundance_path"])[
        "p0"
    ],
    "photospheric": scipy.io.readsav(
        _abundance_model_file_path["photospheric_abundance_path"]
    )["p0"],
}


def resolve_abundance_model_type(abundance_model):
    """Formats users abundance_model name."""
    if not isinstance(abundance_model, str):
        raise TypeError("Abundance model name must be a string")
    abundance_name = abundance_model.lower()
    list_of_abundance_name = ["coronal", "hybrid", "photospheric"]
    if abundance_name not in list_of_abundance_name:
        raise ValueError(
            f"\n{abundance_name} is not a current model for XRTpy.\n"
            "Available abundance models:\n"
            "Coronal, Hybrid, and Photospheric.\n"
        )
    return abundance_name


CHIANTI_file = {
    "abundance_model": _abundance_model_data["coronal"]["ABUND_MODEL"][0],
    "dens_model": _abundance_model_data["coronal"]["DENS_MODEL"][0],
    "ioneq_model": _abundance_model_data["coronal"]["IONEQ_MODEL"][0],
    "name": _abundance_model_data["coronal"]["NAME"][0],
    "spectra": _abundance_model_data["coronal"]["SPEC"],
    "spectra_units": _abundance_model_data["coronal"]["SPEC_UNITS"][0],
    "temperature": _abundance_model_data["coronal"]["TEMP"][0],
    "temp_units": _abundance_model_data["coronal"]["TEMP_UNITS"][0],
    "tlength": _abundance_model_data["coronal"]["TLENGTH"][0],
    "wlength": _abundance_model_data["coronal"]["WLENGTH"][0],
    "wavelength": _abundance_model_data["coronal"]["WAVE"][0],
    "wavelength_units": _abundance_model_data["coronal"]["WAVE_UNITS"][0],
}


class TemperatureResponseFundamental:
    """Produce the temperature response for each XRT x-ray channel, assuming a spectral emission model."""

    def __init__(self, filter_name, observation_date, abundance_model):
        self._name = resolve_filter_name(filter_name)
        self.observation_date = observation_date
        self._channel = Channel(self.filter_name)
        self._abundance_model = resolve_abundance_model_type(abundance_model)

    @property
    def filter_name(self):
        """Name of searched filter."""
        return self._name

    @property
    def abundances(self) -> Dict[str, Real]:
        return self._abundance_model

    @property
    def abundance_model(self):
        """A brief description of what abundance model was used in the creation of the emission spectra."""
        return self._abundance_model

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
    def get_abundance_data(self):
        abundance_type = self.abundances
        data = _abundance_model_data[abundance_type]
        if abundance_type == "coronal":
            # print(f'Requested data: {abundance_type}')
            return {
                "CHIANTI_abundance_model": data["ABUND_MODEL"][0],
                "dens_model": data["DENS_MODEL"][0],
                "ioneq_model": data["IONEQ_MODEL"][0],
                "name": data["NAME"][0],
                "spectra": data["SPEC"],
                "spectra_units": data["SPEC_UNITS"][0],
                "temperature": data["TEMP"][0],
                "temp_units": data["TEMP_UNITS"][0],
                "tlength": data["TLENGTH"][0],
                "wlength": data["WLENGTH"][0],
                "wavelength": data["WAVE"][0],
                "wavelength_units": data["WAVE_UNITS"][0],
            }
        else:
            ValueError("Unable to process data")
            # import pdb; pdb.set_trace()

    @property
    def CHIANTI_version(self):
        """Name of the emission model."""
        return CHIANTI_file["name"]

    @property
    def version(self):
        """Name of the emission model."""
        return self.get_abundance_data["name"]

    @property
    def abundance_model(self):
        """A brief description of what abundance model was used in the creation of the emission spectra."""
        return self.get_abundance_data["abundance_model"]

    @property
    def density_model(self):
        """A brief description of the plasma density, emission measure,or differential emission measure that was used in the creation of the emission spectra."""
        return self.get_abundance_data["dens_model"]

    @property
    def ionization_model(self):
        """A brief description of that ionization equilibrium model was used in the creation of the emission spectra."""
        return self.get_abundance_data["ioneq_model"]

    @property
    @u.quantity_input
    def temperature(self):
        """CHIANTI temperatures in kelvin."""
        return u.Quantity(self.get_abundance_data["temperature"] * u.K)

    @property
    def file_spectra(self):
        """CHIANTI file spectra."""
        return self.get_abundance_data["spectra"]

    @property
    @u.quantity_input
    def wavelength(self):
        """CHIANTI file wavelength values in Å."""
        return u.Quantity(self.get_abundance_data["wavelength"] * u.Angstrom)

    @property
    @u.quantity_input
    def channel_wavelength(self):
        """Array of wavelengths for every X-ray channel in Å."""
        return u.Quantity((Channel(self.filter_name).wavelength[:3993]) * u.photon)

    @property
    def focal_len(self):
        """Focal length of the telescope in units of cm."""
        return Channel(self.filter_name).geometry.geometry_focal_len

    @property
    def ev_per_electron(self):
        """Amount of energy it takes to dislodge 1 electron in the CCD."""
        return Channel(self.filter_name).ccd.ccd_energy_per_electron

    @property
    @u.quantity_input
    def pixel_size(self) -> u.cm:
        """CCD pixel size. Units converted from μm to cm."""
        ccd_pixel_size = Channel(self.filter_name).ccd.ccd_pixel_size
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
                self.wavelength,
                CHIANTI_file["spectra"][0][i],
                kind="linear",  ##########Correct##########
            )
            spectra_interpolate.append(interpolater(self.channel_wavelength))
        return spectra_interpolate * (
            u.photon * u.cm**3 * (1 / u.sr) * (1 / u.s) * (1 / u.Angstrom)
        )

    @u.quantity_input
    def effective_area(self) -> u.cm**2:
        return effective_area(self.filter_name, self.observation_date)

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
        return Channel(self.filter_name).ccd.ccd_gain_right

    @u.quantity_input
    def temperature_response(self) -> u.DN * u.cm**5 / (u.s * u.pix):
        r"""Apply gain value to the Temperature Response in units of DN cm\ :sup:`5` s\ :sup:`-1` pix\ :sup:`-1`."""
        return self.integration() / self.ccd_gain_right
