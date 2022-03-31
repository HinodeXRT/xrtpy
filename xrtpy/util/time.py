__all__ = [
    "epoch",
    "xrt_data_time_to_dt",
]

from datetime import datetime 
from datetime import timedelta

#Hinode-XRT mission elapsed time "Epoch" is Sept 22, 2006 21:36:00.
epoch = datetime(year=2006, month=9, day=22, hour=21, minute=36, second=0)


def xrt_data_time_to_dt(data_time,epoch):
    '''Covering data time (float64) to a datetime object.'''
    data_dates_dt = []
    data_dates_seconds = []

    for time in data_time: 
        t0=data_time[0]
        t1=time
        dt = t1-t0
        data_dates_dt.append(epoch+timedelta(0,dt))
        data_dates_seconds.append(float((epoch+timedelta(0,dt)).strftime('%s')))
    return(data_dates_dt,data_dates_seconds)