import astropy.time

__all__ = [
    "epoch",
]

from datetime import datetime, timedelta

# Hinode-XRT mission elapsed time "Epoch" is Sept 22, 2006 21:36:00.
epoch = astropy.time.Time("2006-09-22 21:36:00")


def xrt_data_time_to_dt(data_time, epoch):
    """
    Converting data time (float64) to a datetime object.

    Parameters
    ----------
    data_time : real number (?)
        Description...
    epoch : `datetime.datetime`
        This function will convert the requested date and time into a datetime
         object in seconds from the respected launched date to collect the correct date.
    """
    data_dates_dt = []
    data_dates_seconds = []

    for time in data_time:
        t0 = data_time[0]
        t1 = time
        dt = t1 - t0
        data_dates_dt.append(epoch + timedelta(0, dt))
        data_dates_seconds.append(float((epoch + timedelta(0, dt)).strftime("%s")))
    return (data_dates_dt, data_dates_seconds)
