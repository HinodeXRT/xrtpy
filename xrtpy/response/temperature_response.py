__all__ = [
    "TemperatureResponseFundamental",
]

from pathlib import Path

import astropy.constants as const
import numpy as np
import scipy.io
from astropy import units as u
from scipy import interpolate
from xrtpy.response.channel import Channel
from xrtpy.response.effective_area import EffectiveAreaFundamental, parse_filter_input

_abundance_model_file_path = {
    "coronal_abundance_path": Path(__file__).parent.absolute()
    / "data/chianti_emission_models"
    / "XRT_emiss_model.default_CHIANTI.geny",
    "hybrid_abundance_path": Path(__file__).parent.absolute()
    / "data/chianti_emission_models"
    / "XRT_emiss_model.default_CHIANTI_hybrid.geny",
    "photospheric_abundance_path": Path(__file__).parent.absolute()
    / "data/chianti_emission_models"
    / "XRT_emiss_model.default_CHIANTI_photospheric.geny",
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

_list_of_abundance_name = ["coronal", "hybrid", "photospheric"]


def _resolve_abundance_model_type(abundance_model):
    """Formats and checks users abundance model input name."""
    if not isinstance(abundance_model, str):
        raise TypeError("Abundance model name must be a string")
    abundance_name = abundance_model.lower()
    if abundance_name not in _list_of_abundance_name:
        raise ValueError(
            f"\n{abundance_name} is not a current abundance model for XRTpy.\n"
            "Available abundance models:\n"
            "'coronal', 'hybrid', and 'photospheric'.\n"
        )
    return abundance_name


class TemperatureResponseFundamental:
    """Produce the temperature response for each XRT x-ray channel, assuming a spectral emission model."""

    def __init__(self, filter_name, observation_date, abundance_model="coronal"):
        """
        Initialize the TemperatureResponseFundamental class.

        Parameters
        ----------
        filter_name : str
            The name of the filter.
        observation_date : str or date-time object
            The date of the observation.
        abundance_model : str, optional
            The abundance model to use. Options are 'coronal' (default), 'hybrid', and 'photospheric'. Default abundance model is coronal.
        """
        parsed_filter = parse_filter_input(filter_name)
        self._name = filter_name  # keep original name (like "Al-poly/Ti-poly")
        self.filter1_name = parsed_filter.filter1
        self.filter2_name = parsed_filter.filter2
        self.is_combo = parsed_filter.is_combo
        self._channel = Channel(self.filter1_name)

        self._abundance_model = _resolve_abundance_model_type(abundance_model)
        self._effective_area_fundamental = EffectiveAreaFundamental(
            filter_name, observation_date
        )
        print(
            f"\nFilter1: {self.filter1_name}, Filter2: {self.filter2_name}, Combo: {self.is_combo}\n"
        )

    @property
    def filter_name(self):
        """Name of searched filter."""
        return self._name

    @property
    def abundances(self) -> str:
        """Defined the name of the requested abundance model as a string."""
        return self._abundance_model

    @property
    def observation_date(self):
        """Date of observation."""
        return self._effective_area_fundamental.observation_date

    @property
    def _get_abundance_data(self):
        """
        Return the requested abundance data used to calculate the temperature response.

        Returns
        -------
        dict
            A dictionary containing the abundance model data.
        """
        abundance_type_name = self.abundances
        data = _abundance_model_data[abundance_type_name]
        if abundance_type_name not in _list_of_abundance_name:
            raise ValueError("Unable to process data.")
        return {
            "abundance_model_info": data["ABUND_MODEL"][0],
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

    @property
    def chianti_abundance_version(self):
        """Version of the chianti abundance model."""
        return self._get_abundance_data["name"]

    @property
    def abundance_model_information(self):
        """A brief description of the abundance model used in the creation of the emission spectra."""
        return self._get_abundance_data["abundance_model_info"]

    @property
    def density_model(self):
        """A brief description of the plasma density, emission measure, or differential emission measure used in the creation of the emission spectra."""
        return self._get_abundance_data["dens_model"]

    @property
    def ionization_model(self):
        """A brief description of the ionization equilibrium model used in the creation of the emission spectra."""
        return self._get_abundance_data["ioneq_model"]

    @property
    @u.quantity_input
    def CHIANTI_temperature(self):
        """Emission model temperatures in kelvin."""
        return u.Quantity(self._get_abundance_data["temperature"] * u.K)

    @property
    def file_spectra(self):
        """Emission model file spectra."""
        return self._get_abundance_data["spectra"][0]

    @property
    @u.quantity_input
    def _wavelength_spectra(self):
        """Emission model file wavelength values in Å."""
        return u.Quantity(self._get_abundance_data["wavelength"] * u.Angstrom)

    @property
    @u.quantity_input
    def wavelength(self):
        """Array of wavelengths for every X-ray channel in Å."""
        return self._effective_area_fundamental.wavelength

    @property
    def focal_len(self):
        """Focal length of the telescope in units of cm."""
        return self._channel.geometry.geometry_focal_len

    @property
    def ev_per_electron(self):
        """Amount of energy it takes to dislodge 1 electron in the CCD."""
        return self._channel.ccd.ccd_energy_per_electron

    @property
    @u.quantity_input
    def pixel_size(self) -> u.cm:
        """CCD pixel size. Units converted from μm to cm."""
        ccd_pixel_size = self._channel.ccd.ccd_pixel_size
        return ccd_pixel_size.to(u.cm)

    @property
    @u.quantity_input
    def solid_angle_per_pixel(self) -> u.sr / u.pix:
        """This amount represents the solid angle, which is given in units of steradians over pixel."""
        return (self.pixel_size / self.focal_len) ** 2 * (u.sr / u.pix)

    @u.quantity_input
    def spectra(self) -> u.photon * u.cm**3 / (u.sr * u.s * u.Angstrom):
        """
        Interpolation between the spectra wavelength onto the channel wavelength.

        Returns
        -------
        numpy.ndarray
            Interpolated spectra values.
        """
        spectra_interpolate = []
        for i in range(61):
            interpolater = interpolate.interp1d(
                self._wavelength_spectra.to_value("AA"),
                self.file_spectra[i],
                kind="linear",
            )
            spectra_interpolate.append(interpolater(self.wavelength.to_value("AA")))
        return spectra_interpolate * u.Unit("photon cm3 sr-1 s-1 Angstrom-1")

    @u.quantity_input
    def effective_area(self) -> u.cm**2:
        """
        Calculate the effective area.

        Returns
        -------
        astropy.units.Quantity
            Effective area in cm^2.
        """
        return self._effective_area_fundamental.effective_area()

    @u.quantity_input
    def integration(self) -> u.electron * u.cm**5 / (u.s * u.pix):
        """
        Perform the integration of the temperature response.

        Returns
        -------
        astropy.units.Quantity
            Integrated temperature response in electron cm^5 / (s pix).
        """
        constants = const.h * const.c / self.wavelength / u.photon
        constants *= self.solid_angle_per_pixel / self.ev_per_electron
        # Simple summing like this is appropriate for binned data like in the current
        # spectrum file. More recent versions of Chianti include the line width,
        # which then makes the previous version that uses Simpson's method
        # to integrate more appropriate (10/05/2022)
        return (
            self.spectra()
            * self.effective_area()
            * constants
            * np.gradient(self.wavelength)
        ).sum(axis=1)

    @property
    @u.quantity_input
    def ccd_gain_right(self) -> u.electron / u.DN:
        """Provide the camera gain in electrons per data number."""
        return self._channel.ccd.ccd_gain_right

    @u.quantity_input
    def temperature_response(self) -> u.DN * u.cm**5 / (u.s * u.pix):
        """
        Apply gain value to the Temperature Response.

        Returns
        -------
        astropy.units.Quantity
            Temperature response in DN cm^5 / (s pix).
        """
        r"""Apply gain value to the Temperature Response in units of DN cm\ :sup:`5` s\ :sup:`-1` pix\ :sup:`-1`."""
        return self.integration() / self.ccd_gain_right
