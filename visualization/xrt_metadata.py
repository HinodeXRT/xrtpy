
from astropy.io import fits
import astropy.units as u

import numpy as np

from astropy.coordinates import SkyCoord

from sunpy.coordinates import frames
from astropy.time import Time
import scipy.io as sio
from scipy.io import readsav
import urllib.request
import sys
sys.path.append('/Users/ntrueba/SOLAR/code/xrtstuff/') #
import xrt_metadata_plot as xplt
import xrt_metadata_download as xfetch



from ipywidgets import Layout, interact, IntSlider,IntProgress, RadioButtons, FloatSlider,FloatRangeSlider,fixed

#This module is designed to download and manage XRT observation metadata for previewing purposes

# The first part contains the xrt_meta object, which is how we organize the metadata 
# After this, we have the download functions

class xrt_meta:
    # The xrt_meta object accepts a list of .fits headers and organizes them
    # Crucially, it separates metadata quantities by xrt filter  
 
    def __init__(self, xrt_downloaded_files):
        header_lis = xfetch.fetch_metadata(xrt_downloaded_files)
        self.head_lis = header_lis # stores the unfiltered header list 
        self.fkey_lis = ['EC_FW1_','EC_FW2_'] # xrt filter keywords in fits header, used by sort_xfilter function to do xrt filter sorting
        self.check_bool = True # Testing parameter, remove
        # upon initialization, the object filters observations and creates the .metadata dictonary 
        self.metadata = self.organize_metadata()


    def organize_metadata(self):
        # testing statement, remove
        if self.check_bool:
            print('Meta')

        # sort observations by filter.
        # filter_lis is a simple list containing strings for all filters in the data set (n_filter)
        # hbtest is a list of length n_filter, each containing a list of headers for all observations for each filter
        hbtest, filter_lis = self.sort_xfilters()


        # Get list of DATE_OBS for each observation
        # time_f_lis: obseration times for xrt obs, separated by xrt filter (same shape as hbtest)
        # time_lis is a flattened version of this, it is not sorted and should only be with the get_time_range() function
        time_lis, time_f_lis = self.get_time_grid(hbtest)

        # tmin and tmax are the time of first and last observations
        # dts: the time separation in seconds, the difference between first and last observations
        tmin, tmax, dts = self.get_time_range(time_lis)

        # filtered list of delta time in seconds - the difference between each observation time in seconds relative to the first observation
        # this is very useful for quickly sincronizing observations when plotting
        delta_t_lis = self.get_delta_T(time_f_lis, tmin)

        # The filtered metadata dictionary - more quantities can be easily added here 
        #asdasf
        xdata_set = {'HEADER_LIST': hbtest,
                     'FILTER_LIST': filter_lis,
                       'TIME_LIST': time_f_lis,
                       'MIN_TIME': tmin,
                       'MAX_TIME': tmax,
                       'TIME_RANGE':dts,
                       'DELTA_TIME_LIST': delta_t_lis}
        return xdata_set
        
    def sort_xfilters(self):
        #test statement
        if self.check_bool:
            print('X1')

        # this function separates an unfilted list of .fits headers from a set of xrt observations and separates them by filter
        hlen = len(self.head_lis) # length of input header list
        filter_lis = [] # this will store the list of xrt filters found in our data set
        header_f_lis = [] # this will store the xrt filter-sorted fits headers as separate lists
        
        for i in range(hlen):
            fkey0 = self.fkey_lis[0] # xrt header keyword
            fkey1 = self.fkey_lis[1] # xrt header keyword
            header_i = self.head_lis[i] # current header
            filter_n0 = header_i[fkey0] #extracting filter #1 from header
            filter_n1 = header_i[fkey1] #extracting filter #2 from header
            if False: # deprecated option to filter by FOV (partial vs FULL), can re-introduce
                fovx = header_i['FOVX'] 
                filter_n2 = ''
                if (fovx > 1000.0):
                    filter_n2 = '(FULL)'

            # the combined filter string 
            filter_n = filter_n0 +'/'+filter_n1#+filter_n2

            # is this a new filter?
            if filter_n not in filter_lis:
                #yes?
                filter_lis.append(filter_n) # add it to the list of filters we have in our data set
                header_f_lis.append([header_i]) # add a new list of headers for this filter, containing the current header
            else:
                # no?
                filter_i = filter_lis.index(filter_n) # what is the index of the filter?
                header_f_lis[filter_i].append(header_i) # append this header to the list corresponding to this filter
        # filter_lis are appended simultaneously, guaranteeing that filter_lis is organized the same as the first dimension of header_f_lis
        # this process only captures the filters in our observations with no need to hard code anythin
        return header_f_lis, filter_lis
    
    def get_time_grid(self, header_f_lis):
        #test statement
        if self.check_bool:
            print('X2')
        time_lis = []
        time_f_lis = []
        for i in range(len(header_f_lis)):
            tlis = []
            for j in range(len(header_f_lis[i])):
                htemp = header_f_lis[i][j]
                htime = htemp['DATE_OBS']
                tlis.append(htime)
                time_lis.append(htime)
            time_f_lis.append(tlis)
        return time_lis, time_f_lis
    
    def get_time_range(self, time_lis):
        # check statement
        if self.check_bool:
            print('X3')

        # this function gives you:
        # min_time : earliest obs
        # max_time : latest obs
        # time_range : difference between these two times in seconds 

        time_obj = Time(np.asarray(time_lis), format='isot', scale='utc')
        time_delt = time_obj - time_obj[0]
        time_delt = time_delt.value*24.0*3600.0
        
        min_i = np.argmin(time_delt)
        max_i = np.argmax(time_delt)
        min_time = time_lis[min_i]
        max_time = time_lis[max_i]

        time_range = Time(max_time, format='isot', scale='utc') - Time(min_time, format='isot', scale='utc')
        time_range = time_range.value*24.0*3600.0
        return min_time, max_time, time_range
    
    def get_delta_T(self, time_lis, min_time):

        # gives you an xrt-filter separated list:
        # For each filter, an array of observation times relative to the earliest observation in the data set, in seconds

        if self.check_bool:
            print('X4')

        dtlis = []
        for i in range(len(time_lis)):
            dtt = Time(time_lis[i], format='isot', scale='utc')-Time(min_time, format='isot', scale='utc')
            dtlis.append(dtt.value*24.0*3600.0)
        return dtlis 

    def get_frame_lis(self):

        # gives you an xrt-filter separated list:
        # For each filter, a list containing coordinates for the four corners of the FOV for each observation
        frame_f_lis = []
        bhl = self.metadata['HEADER_LIST']
        for i in range(len(bhl)):
            # i represents each filter_i
            frame_lis = []
            for j in range(len(bhl[i])):
                # j is each obs in filter_i
                frame_lis.append([self.get_xframet(bhl[i][j])]) # calling the function for each obs - currently using floats - see get_xframe for generalized version
            frame_f_lis.append(frame_lis)
        return frame_f_lis
    
    def get_xframet(self, head_i):

        # generates a square of coordines for the FOV for each obs
        # simple function, needs to be generalized to include spacecraft roll 
        # using float values instead of astropy coordinates in a specific frame - might be desirable to work in a sunpy ecosystem

        fovx, fovy = head_i['FOVX'], head_i['FOVY']
        xcen, ycen = head_i['XCEN'], head_i['YCEN']
        
        xc = (xcen + np.asarray([-1.0, 1.0, 1.0,-1.0,-1.0])*fovx*0.5) * u.arcsec
        yc = (ycen + np.asarray([-1.0,-1.0, 1.0,1.0,-1.0])*fovy*0.5) * u.arcsec
        return xc,yc
        
    def get_xframe(self, head_i, map_in = None):
        fovx, fovy = head_i['FOVX'], head_i['FOVY']
        xcen, ycen = head_i['XCEN'], head_i['YCEN']
        
        xc = (xcen + np.asarray([-1.0, 1.0, 1.0,-1.0,-1.0])*fovx*0.5) * u.arcsec
        yc = (ycen + np.asarray([-1.0,-1.0, 1.0,1.0,-1.0])*fovy*0.5) * u.arcsec
        obs_t = head_i['DATE_OBS']
        if map_in is None:
            coords_out = SkyCoord(xc, yc, frame = frames.Helioprojective, observer='earth', obstime=obs_t)
        else:
            coords_out = SkyCoord(xc, yc, frame = map_in.coordinate_frame)
            
        return coords_out
    

    def plot_preview(self):
        # plot preamble here, needs to be cleaned up
        # because we want to have an interactive plot, we don't want to recalculate some of these basic things every time the plot refreshes, so we calculate it here

        
        #self.head_lis[0]['RSUN_OBS']
        trange = self.metadata['TIME_RANGE']# Time from first frame to last in seconds
        filter_lis = self.metadata['FILTER_LIST'].copy() 
        flen = len(filter_lis)
        filter_lis_b = ['All']
        filter_lis_b = filter_lis_b + filter_lis #filter list for interactive plot, including the 'All' option to plot all filters

        fov_lis = self.get_frame_lis()# FOV list (n_filters, n_frames(filter))
        time_lis0 = self.metadata['DELTA_TIME_LIST'] # Frame time in seconds relative to first frame (n_filters, n_frames(filter))
        time_lis_abs = self.metadata['TIME_LIST'] # Absolute time for each frame (n_filters, n_frames(filter))
        col_vals = xplt.get_pcol(flen)

        # Set-up Solar limb and grid lines 
        ang_arr = np.linspace(0.0,6.28,100)
        lon_arr = np.sin(np.arange(0.0,0.01+np.pi/2.0,np.pi/12.0))
        rsun_p = 944.0 # this should not be hard-coded
        xca, yca = np.cos(ang_arr)*rsun_p, np.sin(ang_arr)*rsun_p
        yla, zla = [], []
        dt_bool = True #print the time distance between XRT frame and current frame?
        fill_bool = True #use fill between function? makes preview slower



        t0d = np.datetime64(self.metadata['MIN_TIME']) #T0 for the observations
        trp = [t0d,(np.timedelta64(int(np.round(trange)),'s')+t0d)] #Range for the full timeline

        for j in range(len(lon_arr)):
            zval = rsun_p*lon_arr[j]
            yline = rsun_p*((1.0 - lon_arr[j]**2.0)**0.5)
            zla.append(zval)
            yla.append(yline)
        inputs = {}
        inputs['xca']  = xca 
        inputs['yca']  = yca     
        inputs['rsun_p'] = rsun_p  
        inputs['zla'] = zla     
        inputs['lon_arr'] = lon_arr 
        inputs['yla'] = yla     
        inputs['xmeta'] = self
        inputs['flen'] = flen    
        inputs['t0d'] = t0d     
        inputs['trp'] = trp     
        inputs['dt_bool'] = dt_bool 
        inputs['fov_lis']= fov_lis 
        inputs['fill_bool']  =  fill_bool     
        inputs['filter_lis']  = filter_lis     
        inputs['time_lis0'] = time_lis0      
        inputs['col_vals'] = col_vals       
        inputs['time_lis_abs'] = time_lis_abs 

        # need to generalize when to use interactive vs animation etc
        interact(xplt.interactive_plot, 
                mini_tline = True,
                filter_selection = RadioButtons(options=filter_lis_b,index=0),
                zoom_v=FloatSlider(value=1.0,min=0.5,max=2.0, layout=Layout(width='70%')),
                tzoom = FloatRangeSlider(value=[0.0,1.0],min=0,max=1.0,step=0.01, layout=Layout(width='70%')),
                time_ind=IntSlider(value=50,min=0,max=99,step=1, layout=Layout(width='70%')),
                inputs = fixed(inputs))

        return