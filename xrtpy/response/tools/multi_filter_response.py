from dataclasses import dataclass

from astropy import units as u

from xrtpy.response.temperature_response import TemperatureResponseFundamental


@dataclass
class _TemperatureResponseData:
    """
    Internal dataclass representing the temperature response for a single XRT filter channel.

    Attributes
    ----------
    filter_name : str
        Name of the XRT filter.
    temperature : astropy.units.Quantity
        Temperature grid in Kelvin.
    response : astropy.units.Quantity
        Instrument response in DN cm^5 / (pix s).
    """

    filter_name: str
    temperature: u.Quantity
    response: u.Quantity


def generate_temperature_responses(filters, obs_date, abundance="Coronal"):
    """
    Generate temperature response objects for one or more XRT filters.

    Parameters
    ----------
    filters : list of str
        Names of XRT filter channels (e.g., ["Al-mesh", "Be-thin", "Al-poly/Ti-poly"]).
    obs_date : str or astropy.time.Time
        Observation date to compute the instrument response (e.g., "2011-01-28T11:02:31").
    abundance : str, optional
        CHIANTI abundance model to use: "Coronal" (default), "Hybrid", or "Photospheric".

    Returns
    -------
    list of TempResponse
        List of response objects containing filter name, CHIANTI temperature grid, and the instrument response.

    Notes
    -----
    - Temperatures are returned in Kelvin.
    - Responses are in units of DN cm^5 s^-1 pix^-1.
    - Supports both single and combined filters (e.g., "Al-poly/Ti-poly").
    """
    responses = []
    for f in filters:
        obj = TemperatureResponseFundamental(f, obs_date, abundance)
        responses.append(
            _TemperatureResponseData(
                filter_name=obj.filter_name,
                temperature=obj.CHIANTI_temperature,
                response=obj.temperature_response(),
            )
        )
    return responses
