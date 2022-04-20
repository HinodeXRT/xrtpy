__all__ = [
    "EffectiveAreaFundamental",
    "effective_area",
]

import pkg_resources
import sunpy.io.special
import numpy as np
import sunpy.time
import math 
import scipy.io
import sunpy.io.special

from scipy import interpolate
from datetime import timedelta
from astropy import units as u
from functools import cached_property 

from xrtpy.util.time import epoch
from xrtpy.response.channel import Channel 
from xrtpy.response.channel import resolve_filter_name


index_mapping_to_fw1_name  = {
    "Open": 0,
    "Al-poly": 1,
    "C-poly": 2,
    "Be-thin": 3,
    "Be-med": 4,
    "Al-med": 5,
}

index_mapping_to_fw2_name = { 
    "Open": 0,
    "Al-mesh": 1,
    "Ti-poly": 2,
    "G-band": 3,
    "Al-thick": 4,
    "Be-thick": 5,
}


_ccd_contam_filename = pkg_resources.resource_filename( "xrtpy","data/channels/xrt_contam_on_ccd.geny")
_filter_contam_filename = pkg_resources.resource_filename( "xrtpy","data/channels/xrt_contam_on_filter.geny")

_ccd_contam_file = scipy.io.readsav(_ccd_contam_filename)
_filter_contam_file =scipy.io.readsav(_filter_contam_filename)

#CCD contam geny files keys for time and date.
_ccd_contamination_file_time = _ccd_contam_file['p1']
_ccd_contamination = _ccd_contam_file['p2']

#Filter contam geny files keys for time and date.
_filter_contamination_file_time =_filter_contam_file['p1']
_filter_contamination = _filter_contam_file['p2']

class EffectiveAreaFundamental:

    """
    Class for calculating the effective area.
    
    Parameters
    -----------
    filter_name : str
        The name of the filter.
        
    observation_date
        The date of the observation.  For valid date formats, look at the documentation for
        `sunpy.time.parse_time`.
    """
    def __init__(self,filter_name, observation_date):
        self._name = resolve_filter_name(filter_name)
        self.observation_date  = observation_date   
        self._channel = Channel(self.name)
    
    @property
    def name(self):
        """Name of Filter."""
        return self._name

    @property
    def observation_date(self):
        """Date of observation."""
        return self._observation_date 
    
    @observation_date.setter
    def observation_date(self,date):
        """Users observation date input."""
        astropy_time = sunpy.time.parse_time(date) #Astropy time in utc
        observation_date = astropy_time.datetime  

        if observation_date <= epoch:
            raise ValueError('Invalid date: {:}.\n Date must be after September 22nd, 2006 21:36:00.'.format(observation_date))
        self._observation_date = observation_date
    

    @property 
    def ccd_observation_date_to_seconds(self):
        """Converting users observation date into secounds. Used for interpolation."""

        ccd_observation_date_to_seconds = []

        for time in _ccd_contamination_file_time:
            t0=_ccd_contamination_file_time[0]
            t1=time
            dt = t1-t0
            ccd_observation_date_to_seconds.append(( self.observation_date + timedelta(0,dt) ).strftime('%s'))

        return ccd_observation_date_to_seconds[0]
    
    @property
    def ccd_data_dates_to_seconds(self):
        """Converting CCD data dates to datetimes."""
        
        ccd_data_dates_to_seconds = []

        for time in _ccd_contamination_file_time:
            t0=_ccd_contamination_file_time[0]
            t1=time
            dt = t1-t0
            ccd_data_dates_to_seconds.append(float((epoch+timedelta(0,dt)).strftime('%s')))

        return ccd_data_dates_to_seconds

    @property 
    def filter_observation_date_to_seconds(self): 
        """Converting users observation date into seconds. Used for interpolation."""

        filter_observation_date_to_seconds= []
        
        for time in _filter_contamination_file_time:
            t0=_filter_contamination_file_time[0]
            t1=time
            dt = t1-t0
            filter_observation_date_to_seconds.append(( self.observation_date + timedelta(0,dt)).strftime('%s'))
            
        return filter_observation_date_to_seconds[0]

    @property    
    def filter_data_dates_to_seconds(self):
        """Converting Filter contamination data dates to datetimes."""
    
        filter_data_dates_to_seconds = []

        for time in _filter_contamination_file_time:
            t0=_filter_contamination_file_time[0]
            t1=time
            dt = t1-t0
            filter_data_dates_to_seconds.append(float((epoch + timedelta(0,dt)).strftime('%s')))

        return filter_data_dates_to_seconds

    @property
    def contamination_on_CCD(self):
        """Calculation of contamination layer on the CCD, thickness given in Angstrom (Å)."""

        interpolater = scipy.interpolate.interp1d(self.ccd_data_dates_to_seconds, _ccd_contamination,kind='linear')
        ccd_contam_interpolated_date  =  interpolater(self.ccd_observation_date_to_seconds)
        
        return int(ccd_contam_interpolated_date)

    @property
    def filter_index_mapping_to_name(self):
        """Returns filter's corresponding number value."""
        if self.name in index_mapping_to_fw1_name:
            return( index_mapping_to_fw1_name.get(self.name))
        elif self.name in index_mapping_to_fw2_name:
            return(index_mapping_to_fw2_name.get(self.name) )
        
    @property
    def filter_wheel_number(self):
        """Defining choosen filter to its corresponding filter wheel."""
        return 0 if self.name in index_mapping_to_fw1_name else 1
        
    @property
    def filter_data(self):
        """Collecting filter data."""
        return (_filter_contamination[self.filter_index_mapping_to_name][self.filter_wheel_number])
    
    @property 
    def contamination_on_filter(self):
        """Calculation of contamination layer on a filter,thickness giving in Angstrom Å."""

        interpolater = scipy.interpolate.interp1d( self.filter_data_dates_to_seconds , self.filter_data ,kind='linear')
        filter_contam_interpolated_date = interpolater( self.filter_observation_date_to_seconds )

        return int(filter_contam_interpolated_date)

    @cached_property 
    def n_DEHP_attributes(self):
        """Diethylhexylphthalate: Wavelength (nm), Delta, Beta."""
        _n_DEHP_filename = pkg_resources.resource_filename( "xrtpy","response/data/n_DEHP.txt") 
        
        with open(_n_DEHP_filename, "r") as n_DEHP: 
            list_of_DEHP_attributes = []
            for line in n_DEHP:
                stripped_line = line.strip()
                line_list = stripped_line.split()
                list_of_DEHP_attributes.append(line_list)   

        return list_of_DEHP_attributes

    @cached_property
    def n_DEHP_wavelength(self):
        """Diethylhexylphthalate: Wavelength given in Angstrom (Å)."""
        
        wavelength_str = [] #nm

        for i in range(2,len(self.n_DEHP_attributes)):
            wavelength_str.append(self.n_DEHP_attributes[i][0])

        #Convert wavelength values from nanometers to Angstroms
        wavelength = np.array([float(i)*10 for i in wavelength_str])

        return(wavelength)

    @cached_property
    def n_DEHP_delta(self):
        """Diethylhexylphthalate: Delta."""

        delta_str = []
      
        for i in range(2,len(self.n_DEHP_attributes)):
            delta_str.append(self.n_DEHP_attributes[i][1])

        #Converting from str to float 
        delta_float = np.array([float(delta_str[i]) for i in range(0,len(self.n_DEHP_wavelength))])

        #Interpolate so ranges are the same
        delta = interpolate.interp1d(self.n_DEHP_wavelength, delta_float)(self.n_DEHP_wavelength)

        return(delta)

    @cached_property
    def n_DEHP_beta(self):
        """Diethylhexylphthalate: Beta."""

        beta_str =[]

        for i in range(2,len(self.n_DEHP_attributes)):
            beta_str.append(self.n_DEHP_attributes[i][2]) 

        #Converting from str to float 
        beta_float = np.array([float(beta_str[i]) for i in range(0,len(self.n_DEHP_wavelength))])
        
        #Interpolate so ranges are the same
        beta = interpolate.interp1d(self.n_DEHP_wavelength, beta_float)(self.n_DEHP_wavelength)
        
        return(beta)
   
    @cached_property
    def transmission_equation(self):
        """Defining equations that will be used to calculate the effective area."""
        
        n_o = 1.0 #index of medium at entrance of filter (assumed vacuum)
        n_t = 1.0 #index of medium at exit of filter (assumed vacuum)
        
        incidence_angle = 0 #Angle of incidence on Filter in radians
         
        wavelength_min =  1.   #Angstroms
        wavelength_max = 4000 #Max wavelength in Angstroms
        
        index = [( complex((1 - self.n_DEHP_delta[i]), (1.*self.n_DEHP_beta[i])) ) for i in range(0,4000)]
          
        #Snell's law    
        sin_a = [ ( (n_o* np.sin(incidence_angle)) /index[i]) for i in range(0,4000) ]
        
        #cos_a = [ (math.sqrt(1- (sin_a[i]**2) ) )  for i in range(0,wavelength_max) ]
        cos_a = 1
    
        return(index,sin_a,cos_a,wavelength_max,n_o,n_t,incidence_angle)
        
    @cached_property
    def angular_wavenumber_CCD(self):
        """Define angular wavenumber on CCD."""

        index,_,cos_a,wavelength_max,_,_,_ = self.transmission_equation
        
        #Define wavevector
        angular_wavenumber = np.array( [((2.*math.pi*index[i]* cos_a)/ self.n_DEHP_wavelength[i] )  for i in range(0,4000)] )
        
        #Multiply by thickness
        angular_wavenumber_thickness = angular_wavenumber* self.contamination_on_CCD 

        real_angular_wavenumber = [(float(angular_wavenumber_thickness[i].real))  for i in range(0,4000)] 
        imaginary_angular_wavenumber = [(angular_wavenumber_thickness[i].imag) for i in range(0,4000) ]                 
        
        kl = [ (complex(real_angular_wavenumber[i],imaginary_angular_wavenumber[i])) for i in range(0,4000) ]

        return kl
    
    @cached_property
    def filterwheel_angular_wavenumber(self): 
        """Define angular wavenumber for a filter."""
        index,_,cos_a,wavelength_max,_,_,_ = self.transmission_equation        
        
        #Define wavevector
        angular_wavenumber = np.array([((2.* math.pi*index[i]* cos_a)/ self.n_DEHP_wavelength[i] )  for i in range(0,4000)])
        
        #Multiply by thickness
        angular_wavenumber_thickness = angular_wavenumber* self.contamination_on_filter
        
        real_angular_wavenumber = [(float(angular_wavenumber_thickness[i].real))  for i in range(0,4000)]
        imaginary = [(angular_wavenumber_thickness[i].imag) for i in range(0,4000) ]
        
        kl = [ (complex(real_angular_wavenumber[i],imaginary[i])) for i in range(0,4000) ]

        return kl 
    
    @cached_property
    def CCD_contamination_transmission(self):
        """Calculate transmission matrix coefficient and transmittance on the CCD."""
               
        index,_,_,_,n_o,n_t,_= self.transmission_equation
    
        i_i=complex(0,1) #Define complex number
    
        #Define transfer matrix   
        M = [[[ np.cos(self.angular_wavenumber_CCD[i]), (-i_i*np.sin(self.angular_wavenumber_CCD[i]))/index[i]],[-i_i*np.sin(self.angular_wavenumber_CCD[i])*index[i], np.cos(self.angular_wavenumber_CCD[i])]] for i in range(0,4000)] 
        
        transmittance = [ (2*n_o/( (M[i][0][0]*n_o)+(M[i][0][1]*n_o*n_t)+(M[i][1][0])+(M[i][1][1]*n_t)) ) for i in range(0,4000)]

        transmission = [abs(transmittance[i]**2) for i in range(4000)]
        
        return transmission

    @property
    def channel_wavelength(self):
        """Array of wavelengths for every X-ray channel in Angstroms (Å)."""
        return Channel(self.name).wavelength
    
    @property
    def channel_geometry_aperture_area(self):
        """Hinode/XRT flight model geometry aperture area."""
        return Channel(self.name).geometry.aperture_area

    @property
    def channel_transmission(self):
        return Channel(self.name).transmission

    @property
    def interpolated_CCD_contamination_transmission(self):
        """Interpolate filter contam transmission to the wavelength."""
        CCD_contam_transmission = interpolate.interp1d( self.n_DEHP_wavelength,self.CCD_contamination_transmission)
        CCD_interpolate = CCD_contam_transmission( self.channel_wavelength)
        return CCD_interpolate
    
    @cached_property
    def filter_contamination_transmission(self):
        """Calculate transmission matrix coefficient and transmittance on a filter."""

        index,_,_,_,n_o,n_t,_= self.transmission_equation
        
        i_i = complex(0,1) #Define complex number
       
        #Define transfer matrix   
        M = [[[ np.cos(self.filterwheel_angular_wavenumber[i]), (-i_i*np.sin(self.filterwheel_angular_wavenumber[i]))/index[i]],[-i_i*np.sin(self.filterwheel_angular_wavenumber[i])*index[i], np.cos(self.filterwheel_angular_wavenumber[i])]] for i in range(0,4000)]
            
        transmittance = [ (2*n_o/((M[i][0][0]*n_o)+(M[i][0][1]*n_o*n_t)+(M[i][1][0])+(M[i][1][1]*n_t))) for i in range(0,4000)]

        transmission = [abs(transmittance[i]**2) for i in range(4000)]
        
        return transmission

    @property
    def interpolated_filter_contamination_transmission(self):
        """Interpolate filter contam transmission to the wavelength."""
        Filter_contam_transmission = interpolate.interp1d( self.n_DEHP_wavelength , self.filter_contamination_transmission)
        Filter_interpolate = Filter_contam_transmission( self.channel_wavelength)
        return Filter_interpolate

    @u.quantity_input
    def effective_area(self)-> u.cm**2:
        """Calculation of the Effective Area."""
        return (self.channel_geometry_aperture_area
                * self.channel_transmission
                * self.interpolated_CCD_contamination_transmission
                * self.interpolated_filter_contamination_transmission
               )


def effective_area(filter_name,observation_date):
    EAP = EffectiveAreaFundamental(filter_name,observation_date)
    return( EAP.effective_area() )