__all__ = [
    "epoch",
    "observation_date",
    "validating_data_observation_date",
]

import astropy.time
import os
import scipy.io
import sunpy.time

from datetime import datetime, timedelta
from pathlib import Path

# Hinode-XRT mission elapsed time "Epoch" is Sept 22, 2006 21:36:00.
epoch = astropy.time.Time("2006-09-22 21:36:00")


def xrt_data_time_to_dt(data_time: list, epoch: datetime) -> tuple:
    """
    Convert data times (float64) to datetime objects and seconds from epoch.

    Parameters
    ----------
    data_time : list of float
        A list of float values representing data times.
    epoch : datetime.datetime
        The reference datetime representing the mission epoch.

    Returns
    -------
    tuple
        A tuple containing two lists: data dates as datetime objects and
        data dates as seconds from the epoch.
    """
    data_dates_dt = []
    data_dates_seconds = []

    t0 = data_time[0]  # Initial time

    for time in data_time:
        dt = time - t0
        data_date_dt = epoch + timedelta(seconds=dt)
        data_dates_dt.append(data_date_dt)
        data_dates_seconds.append(float(data_date_dt.strftime("%s")))

    return data_dates_dt, data_dates_seconds


def observation_date(data_time: str) -> datetime:
    """
    Validate a user's requested observation date.

    Parameters
    ----------
    data_time : str
        A string representing the observation date in a recognized format.

    Returns
    -------
    datetime.datetime
        The validated observation date.

    Raises
    ------
    ValueError
        If the observation date is not after the defined epoch.

    Notes
    -----
    This function is used to validate a user's requested observation date before using
    it in functions related to effective area and temperature response functions.

    The `epoch` used for validation is defined as September 22, 2006 21:36:00.
    """

    observation_date = sunpy.time.parse_time(data_time)

    if observation_date <= epoch:
        raise ValueError(
            f"Invalid date: {observation_date.datetime}. "
            f"Date must be after {epoch}."
        )

    return observation_date


###############################################################################
############* Debating if this function belongs here *#########################
###############################################################################


# Absolute path  - Update method
_ccd_contam_filename = Path(
    "/Users/jvelasq/Projects/xrtpy/xrtpy/response/data/xrt_contam_on_ccd.geny"
)

# Read the contamination file
_ccd_contam_file = scipy.io.readsav(_ccd_contam_filename)

# CCD contam geny files keys for time and date.
_ccd_contamination_file_time = astropy.time.Time(
    _ccd_contam_file["p1"], format="utime", scale="utc"
)


def validating_data_observation_date(observation_date):
    """
    Validate the requested observation date against the available data.

    This function checks whether there is data available for the requested observation date.
    If the observation date is later than the last modification date of the contamination data file,
    it raises a ValueError indicating that no contamination data is available for the requested date.

    Parameters
    ----------
    observation_date : `datetime.datetime`
        The requested observation date to be validated.

    Returns
    -------
    observation_date : `datetime.datetime`
        The validated observation date.

    Raises
    ------
    ValueError
        If the observation date is later than the last modification date of the contamination data file,
        indicating that no contamination data is available for the requested date.

    """

    observation_date = observation_date()

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

    return observation_date


###############################################################################
###############################################################################
###############################################################################
##import pdb; pdb.set_trace()
