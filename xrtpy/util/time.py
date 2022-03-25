__all__ = [
    "xrt_contam_time_to_dt",
]

import pkg_resources
import sunpy.io.special
import numpy as np
import smtplib , os, sys
from scipy import interpolate
import collections
import astropy.constants as const
import sunpy
import sunpy.map
from sunpy.data import manager
import scipy.io
import sunpy.io.special
from datetime import date,datetime 
from datetime import timedelta
from astropy.time import Time, TimeDelta

from astropy import units as u

#filename = pkg_resources.resource_filename("xrtpy", "data/channels/xrt_channels_v0016.genx")

#XRT Mission Elapsed Time "Epoch"
epoch = datetime(year=2006, month=9, day=22, hour=21, minute=36, second=0)


def xrt_contam_time_to_dt(year=2012, month=9, day=30, hour=21, minute=36, second=0):
    '''Converting inserted observation date to secounds with a proper reference date.'''
    observation_date_dt=[]
    observation_date_seconds=[]

    for time in data_time:
        t0=data_time[0]
        t1=time
        dt = t1-t0
        observation_date_dt.append((observation_date+timedelta(0,dt)))
        observation_date_seconds.append((observation_date +timedelta(0,dt)).strftime('%s'))
    return(observation_date_dt[0],observation_date_seconds[0])


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