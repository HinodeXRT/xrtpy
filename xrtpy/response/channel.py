__all__ = [
    "channel",
    "Geometry",
    "EntranceFilter",
    "Mirror",
    "Filter",
    "CCD",
    "Channel",
    "resolve_filter_name",
    "EffectiveAreaPreparatory",
]

import pkg_resources
import sunpy.io.special
import numpy as np
import sunpy.time
import math 

import scipy.io
import sunpy.io.special

from scipy import interpolate
from datetime import date,datetime 
from datetime import timedelta

from astropy.time import Time, TimeDelta
from astropy import units as u

from functools import cached_property 

filename = pkg_resources.resource_filename(
    "xrtpy", "data/channels/xrt_channels_v0016.genx"
)

_channel_name_to_index_mapping = {
    "Al-mesh": 0,
    "Al-poly": 1,
    "C-poly": 2,
    "Ti-poly": 3,
    "Be-thin": 4,
    "Be-med": 5,
    "Al-med": 6,
    "Al-thick": 7,
    "Be-thick": 8,
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10,
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12,
    "C-poly/Ti-poly": 13,
}

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

filter_wheel_1 = ("Al-poly","C-poly","Be-thin","Be-med","Al-med")

filter_wheel_2 = ("Al-mesh", "Ti-poly", "G-band", "Al-thick", "Be-thick")

_genx_file = sunpy.io.special.genx.read_genx(filename)["SAVEGEN0"]


def resolve_filter_name(name):
    name = name.replace("_", "-")
    parts: list = name.split("/")
    new_parts: list = [part.capitalize() for part in name.split("/")]
    name: str = "/".join(new_parts)
    return name


class Geometry:
    """The physical geometric parameters of Hinode/XRT."""

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._geom_data = self._genx_file[self._channel_index]["GEOM"]

    @property
    def name(
        self,
    ) -> str:
        """Hinode/XRT flight model geometry."""
        return self._geom_data["LONG_NAME"]

    @property
    @u.quantity_input
    def focal_len(
        self,
    ) -> u.cm:
        """Hinode/XRT flight model geometry focal length."""
        return u.Quantity(self._geom_data["FOC_LEN"], u.cm)

    @property
    @u.quantity_input
    def aperture_area(self) -> u.cm ** 2:
        """Hinode/XRT flight model geometry aperture area."""
        return u.Quantity(self._geom_data["APERTURE_AREA"], u.cm ** 2)


class EntranceFilter:
    """
    Entrance filter properties.

    Thin prefilters cover the narrow annular entrance aperture of the XRT serving two main purposes:

    1. Reduce the visible light entering the instrument.
    2. Reduce the heat load in the instrument.
    """

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._en_filter_data = self._genx_file[self._channel_index]["EN_FILTER"]

    @property
    def entrancefilter_name(
        self,
    ) -> str:
        """Entrance filter name."""
        return self._en_filter_data["LONG_NAME"]

    @property
    def entrancefilter_material(self) -> str:
        """XRT entrance filter material."""
        return self._en_filter_data["MATERIAL"]

    @property
    @u.quantity_input
    def entrancefilter_thickness(self) -> u.angstrom:
        """XRT entrance filter material thickness."""
        return u.Quantity(self._en_filter_data["THICK"], u.angstrom)

    @property
    def entrancefilter_density(self) -> u.g * u.cm ** -3:
        """XRT entrance filter material density in g/cm\ :sup:`3`\ ."""  
        return u.Quantity(self._en_filter_data["DENS"], u.g * u.cm ** -3)

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._en_filter_data["LENGTH"]

    @property
    @u.quantity_input
    def entrancefilter_wavelength(self) -> u.angstrom:
        """Array of wavelengths for entrance filter transmission  in angstroms."""
        return u.Quantity(self._en_filter_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def entrancefilter_transmission(self):
        """Entrance filter transmission."""
        return self._en_filter_data["TRANS"][: self.number_of_wavelengths]

    @property
    def entrancefilter_mesh_transmission(self):
        """Transmission of mesh filter substrate."""
        return self._en_filter_data["MESH_TRANS"]

    @property
    def entrancefilter_substrate(self) -> str:
        """XRT entrance filter substrate."""
        return self._en_filter_data["SUBSTRATE"]

    @property
    def entrancefilter_material(self):
        """Filter material on entrance filter."""
        return self._en_filter_data["MATERIAL"]
    
    @property
    def entrancefilter_thickness(self):
        """Filter material thickness on entrance filter."""
        return self._en_filter_data["THICK"]

class Mirror:
    """
    Grazing incidence mirror properties.

    Grazing-incidence optics used for soft X-ray imaging generally require a minimum of two surfaces.
    Since XRT is a two-bounce telescope, there are two mirror reflectivities.
    """

    _genx_file = _genx_file

    def __init__(self, index, mirror_number):
        self._channel_index = index
        self._mirror_data = self._genx_file[self._channel_index][
            f"MIRROR{mirror_number}"
        ]

    @property
    def name(
        self,
    ) -> str:
        """Hinode/XRT flight model mirror."""
        return self._mirror_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """XRT flight model mirror material."""
        return self._mirror_data["MATERIAL"]

    @property
    @u.quantity_input
    def density(self) -> u.g * u.cm ** -3:
        """Mirror mass density."""
        return u.Quantity(self._mirror_data["DENS"], u.g * u.cm ** -3)

    @property
    @u.quantity_input
    def graze_angle(
        self,
    ) -> u.deg:
        """Mirror graze angle in units of degrees."""
        return u.Quantity(self._mirror_data["GRAZE_ANGLE"], u.deg)

    @property
    @u.quantity_input
    def wavelength(
        self,
    ) -> u.angstrom:
        """Array of wavelengths for mirror reflectance."""
        return u.Quantity(self._mirror_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    @u.quantity_input
    def reflection(
        self,
    ) -> u.angstrom:
        """Reflection of a mirror."""
        return u.Quantity(self._mirror_data["REFL"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._mirror_data["LENGTH"]


class Filter:
    """
    X-ray filters using two filter wheels.

    The corresponding categories are used for both filter 1 and filter 2.
    """

    _genx_file = _genx_file

    def __init__(self, index, filter_number):
        self._channel_index = index
        self._fp_filter_data = self._genx_file[self._channel_index][
            f"FP_FILTER{filter_number}"
        ]

    @property
    def name(
        self,
    ) -> str:
        """XRT focal plane filter position."""
        return self._fp_filter_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """XRT focal plane filter material."""
        return self._fp_filter_data["MATERIAL"]

    @property
    @u.quantity_input
    def thickness(self) -> u.angstrom:
        """XRT  focal plane filter thickness."""
        return u.Quantity(self._fp_filter_data["THICK"], u.angstrom)

    @property
    @u.quantity_input
    def density(self) -> u.g * u.cm ** -3:
        """XRT  focal plane filter density."""
        return u.Quantity(self._fp_filter_data["DENS"], u.g * u.cm ** -3)

    @property
    @u.quantity_input
    def wavelength(self) -> u.angstrom:
        """Array of wavelength for every X-ray focal plane filter in angstroms."""
        return u.Quantity(self._fp_filter_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def transmission(self):
        """Get focal plane filter transmission."""
        return self._fp_filter_data["TRANS"][: self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._fp_filter_data["LENGTH"]

    @property
    def mesh_trans(self):
        """Mesh transmission for the focal plane filter."""
        return self._fp_filter_data["MESH_TRANS"]

    @property
    def substrate(self) -> str:
        """XRT substrate for the focal plane filter."""
        return self._fp_filter_data["SUBSTRATE"]


    @property
    def thickness(self):
        """Filter thickness."""
        return self._fp_filter_data["THICK"]

class CCD:
    """Charge-coupled device on board XRT."""

    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._ccd_data = self._genx_file[self._channel_index]["CCD"]

    @property
    def ccd_name(
        self,
    ) -> str:
        """Hinode/XRT flight model CCD."""
        return self._ccd_data["LONG_NAME"]

    @property
    @u.quantity_input
    def ccd_ev_ore_electron(
        self,
    ) -> u.eV / u.electron:
        """The energy necessary to dislodge one electron."""
        return u.Quantity(self._ccd_data["EV_PER_EL"], u.eV / u.electron)

    @property
    @u.quantity_input
    def ccd_full_well(
        self,
    ) -> u.electron:
        """Number of electrons for a CCD full well."""
        return u.Quantity(self._ccd_data["FULL_WELL"], u.electron)

    @property
    @u.quantity_input
    def ccd_gain_left(
        self,
    ) -> u.electron:
        """Gain when reading the left port of the CCD."""
        return u.Quantity(self._ccd_data["GAIN_L"], u.electron)

    @property
    @u.quantity_input
    def ccd_gain_right(
        self,
    ) -> u.electron:
        """Gain when reading the right port of the CCD."""
        return u.Quantity(self._ccd_data["GAIN_R"], u.electron)

    @property
    def ccd_quantum_efficiency(self):
        """Quantum efficiency of the CCD."""
        return self._ccd_data["QE"][: self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._ccd_data["LENGTH"]

    @property
    @u.quantity_input
    def ccd_pixel_size(self) -> u.micron:
        """CCD pixel size in micrometers."""
        return u.Quantity(self._ccd_data["PIXEL_SIZE"], u.micron)

    @property
    @u.quantity_input
    def ccd_wavelength(self) -> u.angstrom:
        """Array of wavelengths for the CCD quantum efficiency in angstroms."""
        return u.Quantity(self._ccd_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

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

class EffectiveAreaPreparatory:

    """Data for a single filter at a given time"""
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
        """Users observation date imput."""
        astropy_time = sunpy.time.parse_time(date) #Astropy time in utc
        observation_date = astropy_time.datetime  

        mission_start_date = datetime(year=2006, month=9, day=22, hour=21, minute=36, second=0)
        if observation_date <= mission_start_date:
            raise ValueError('Invalid date: {:}.\n Date must be after September 22nd, 2006 21:36:00.'.format(observation_date))
        self._observation_date = observation_date
    

    @property 
    def ccd_observation_date_seconds(self):
        """Converting the observation date to seconds."""
        ccd_obs_date = []
        ccd_observation_date_seconds= []

        #Covering observation date into secounds to interpolate
        for time in _ccd_contamination_file_time:
            t0=_ccd_contamination_file_time[0]
            t1=time
            dt = t1-t0
            ccd_obs_date.append( ( self.observation_date + timedelta(0,dt) ) )
            ccd_observation_date_seconds.append(( self.observation_date + timedelta(0,dt) ).strftime('%s'))

        return ccd_observation_date_seconds[0]
    
    @property
    def ccd_data_dates_seconds(self):
        """Converting CCD data dates to datetimes."""
        ccd_data_dates_dt = []
        ccd_data_dates_seconds = []

        for time in _ccd_contamination_file_time:
            t0=_ccd_contamination_file_time[0]
            t1=time
            dt = t1-t0
            original_date = datetime(2006,9,22,21,36,0)
            ccd_data_dates_dt.append(original_date + timedelta(0,dt))
            ccd_data_dates_seconds.append(float((original_date+timedelta(0,dt)).strftime('%s')))

        return ccd_data_dates_seconds
    

    @property 
    def filter_observation_date_seconds(self): 
        """Converting the observation date to seconds."""

        filter_observation_date_seconds= []
        
        for time in _filter_contamination_file_time:
            t0=_filter_contamination_file_time[0]
            t1=time
            dt = t1-t0
            filter_observation_date_seconds.append(( self.observation_date + timedelta(0,dt)).strftime('%s'))
            
        return filter_observation_date_seconds[0]
        
    @property    
    def filter_data_dates_seconds(self):
        """Converting Filter contamination data dates to datetimes."""
        
        filter_data_dates_dt = []
        filter_data_dates_seconds = []

        for time in _filter_contamination_file_time:
            t0=_filter_contamination_file_time[0]
            t1=time
            dt = t1-t0
            original_date = datetime(2006,9,22,21,36,0)
            filter_data_dates_dt.append(original_date + timedelta(0,dt))
            filter_data_dates_seconds.append(float((original_date + timedelta(0,dt)).strftime('%s')))

        return filter_data_dates_seconds

    
    @property
    def contamination_on_CCD(self):
        """Calculation of contamination layer on the CCD, thickness given in Angstrom (Å)."""

        interpolater = scipy.interpolate.interp1d(self.ccd_data_dates_seconds, _ccd_contamination,kind='linear')
        ccd_contam_interpolated_date  =  interpolater(self.ccd_observation_date_seconds)
        
        return round(int(ccd_contam_interpolated_date))


    @property
    def filter_index_mapping_to_name(self):
        """Returns filter's corresponding number value."""
        if self.name in index_mapping_to_fw1_name.keys():
            return( index_mapping_to_fw1_name.get(self.name))
        elif self.name in index_mapping_to_fw2_name.keys():
            return(index_mapping_to_fw2_name.get(self.name) )
        
    @property
    def filter_wheel_number(self):
        """Defining choosen filter to its corresponding wheel."""
        return 0 if self.name in index_mapping_to_fw1_name.keys() else 1
        
    @property
    def filter_data(self):
        """Collecting filter data."""
        filter_data = _filter_contamination[self.filter_index_mapping_to_name][self.filter_wheel_number]
        return filter_data
    
    @property 
    def contamination_on_filter(self):
        """Calculation of contamination layer on a filter,thickness giving in Angstrom Å."""
        interpolater = scipy.interpolate.interp1d( self.filter_data_dates_seconds , self.filter_data ,kind='linear')
        filter_contam_interpolated_date = interpolater( self.filter_observation_date_seconds )
        return round(int(filter_contam_interpolated_date))

    @cached_property 
    def n_DEHP_attributes(self):
        """Diethylhexylphthalate: Wavelength (nm), Delta, Beta."""
        _n_DEHP_filename = pkg_resources.resource_filename( "xrtpy","response/data/n_DEHP.txt") 
        with open(_n_DEHP_filename, "r") as n_DEHP: 
            list_of_lists = []
            for line in n_DEHP:
                stripped_line = line.strip()
                line_list = stripped_line.split()
                list_of_lists.append(line_list)   
        return list_of_lists

    @cached_property
    def n_DEHP_wavelength(self):
        """Diethylhexylphthalate: Wavelength (nm)."""
        
        wavelength_str = [] #nm

        for i in range(2,len(self.n_DEHP_attributes)):
            wavelength_str.append(self.n_DEHP_attributes[i][0])

        #Convert nm to Angstroms
        wavelength = np.array([float(i)*10 for i in wavelength_str])

        return(wavelength)

    @cached_property
    def n_DEHP_delta(self):
        """Diethylhexylphthalate: Delta."""

        delta_str = []
      
        for i in range(2,len(self.n_DEHP_attributes)):
            delta_str.append(self.n_DEHP_attributes[i][1])

        #Covering from str to float 
        delta = np.array([float(delta_str[i]) for i in range(0,len(self.n_DEHP_wavelength))])

        #Interpolate so ranges are the same
        delta = interpolate.interp1d(self.n_DEHP_wavelength, delta)(self.n_DEHP_wavelength)

        return(delta)

    @cached_property
    def n_DEHP_beta(self):
        """Diethylhexylphthalate: Beta."""

        beta_str =[]

        for i in range(2,len(self.n_DEHP_attributes)):
            beta_str.append(self.n_DEHP_attributes[i][2]) 

        #Converting from str to float 
        beta = np.array([float(beta_str[i]) for i in range(0,len(self.n_DEHP_wavelength))])
        
        #Interpolate so ranges are the same
        beta = interpolate.interp1d(self.n_DEHP_wavelength, beta)(self.n_DEHP_wavelength)
        
        return(beta)

   
    @cached_property
    def transmission_equation(self):
        """Defining equations that will be used to calculate the effective area."""
        n_o = 1.0 #index of medium at entrance of filter (assumed vacuum)
        n_t = 1.0 #index of medium at exit of filter (assumed vacuum)
        
        incidence_angle = 0 #Angle of incidence on Filterin radians
         
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
        
        #Multiply be thickness
        angular_wavenumber_thickness = angular_wavenumber* self.contamination_on_filter
        
        real_angular_wavenumber = [(float(angular_wavenumber_thickness[i].real))  for i in range(0,4000)]
        imaginary = [(angular_wavenumber_thickness[i].imag) for i in range(0,4000) ]
        
        kl = [ (complex(real_angular_wavenumber[i],imaginary[i])) for i in range(0,4000) ]

        return kl 
    
    @cached_property
    def CCD_contamination_transmission(self):
        """Calculate transmission matrix coefficient and transmittance on the CCD."""
               
        index,_,_,_,n_o,n_t,_= self.transmission_equation
    
        i_i=complex(0,1) #Define complete number
    
        #Define transfer matrix   
        M = [[[ np.cos(self.angular_wavenumber_CCD[i]), (-i_i*np.sin(self.angular_wavenumber_CCD[i]))/index[i]],[-i_i*np.sin(self.angular_wavenumber_CCD[i])*index[i], np.cos(self.angular_wavenumber_CCD[i])]] for i in range(0,4000)] 
        
        transmittance = [ (2*n_o/( (M[i][0][0]*n_o)+(M[i][0][1]*n_o*n_t)+(M[i][1][0])+(M[i][1][1]*n_t)) ) for i in range(0,4000)]

        transmission = [abs(transmittance[i]**2) for i in range(4000)]
        
        return transmission

    @property
    def channel_wavelength(self):
        """Array of wavelengths for every X-ray channel in angstroms."""
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
        
        i_i = complex(0,1) #Define complete number
       
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
        """Calculation of the Effecive Area."""
        return (self.channel_geometry_aperture_area
                * self.channel_transmission
                * self.interpolated_CCD_contamination_transmission
                * self.interpolated_filter_contamination_transmission
               )

def effective_area(filter_name,observation_date):
    EAP = EffectiveAreaPreparatory(filter_name,observation_date)
    return( EAP.effective_area() )


class Channel:
    """
    XRTpy

    Available channels: ``"Al-mesh"``, ``"Al-poly"``,  ``"C-poly"``, ``"Ti-poly"``, ``"Be-thin"``, ``"Be-med"``, ``"Al-med"``, ``"Al-thick"``,  ``"Be-thick"`` ,
    ``"Al-poly/Al-mesh"``, ``"Al-poly/Ti-poly"``, ``"Al-poly/Al-thick"``, ``"Al-poly/Be-thick"`` , ``"C-poly/Ti-poly"``
    """

    _genx_file = _genx_file

    def __init__(self, name):
        name = resolve_filter_name(name)
        if name in _channel_name_to_index_mapping:
            self._channel_index = _channel_name_to_index_mapping[name]
            self._channel_data = _genx_file[self._channel_index]
            self._geometry = Geometry(self._channel_index)
            self._entrancefilter = EntranceFilter(self._channel_index)
            self._mirror_1 = Mirror(self._channel_index, 1)
            self._mirror_2 = Mirror(self._channel_index, 2)
            self._filter_1 = Filter(self._channel_index, 1)
            self._filter_2 = Filter(self._channel_index, 2)
            self._ccd = CCD(self._channel_index)
        elif name.lower() == "open": #Complete by adding remaining indexs
            self._sample_channel_data = _genx_file[1]
            self._geometry = Geometry(1)
            self._channel_data = {
                "WAVE": self._sample_channel_data["WAVE"],
                "TRANS":np.ones_like(self._sample_channel_data["TRANS"]),
                "LENGTH": self._sample_channel_data["LENGTH"],
                }
        else:
            raise ValueError(
                f"{name} is not a valid channel. The available channels are: {list(_channel_name_to_index_mapping.keys())}"
            )

    @property
    def geometry(self) -> Geometry:
        return self._geometry

    @property
    def entrancefilter(self) -> EntranceFilter:
        return self._entrancefilter

    @property
    def mirror_1(self) -> Mirror:
        return self._mirror_1

    @property
    def mirror_2(self) -> Mirror:
        return self._mirror_2

    @property
    def filter_1(self) -> Filter:
        return self._filter_1

    @property
    def filter_2(self) -> Filter:
        return self._filter_2

    @property
    def ccd(self) -> CCD:
        return self._ccd

    def __str__(self):
        """Reable printout."""
        return f"XRT Channel for {self.name}"

    def __repr__(self):
        """Code representation."""
        return f"Channel({repr(self.name)})"

    @property
    def name(self) -> str:
        """Name of XRT X-Ray channel."""
        return self._channel_data["NAME"]

    @property
    @u.quantity_input
    def wavelength(self) -> u.angstrom:
        """Array of wavelengths for every X-ray channel in angstroms."""
        return u.Quantity(self._channel_data["WAVE"], u.angstrom)[
            : self.number_of_wavelengths
        ]

    @property
    def transmission(self):
        """Get channel transmission."""
        return self._channel_data["TRANS"][: self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """Data number length."""
        return self._channel_data["LENGTH"]

    @property
    def observatory(self) -> str:
        """Spacecraft: Hinode satellite."""
        return self._channel_data["OBSERVATORY"]

    @property
    def instrument(self) -> str:
        """X-Ray Telescope -XRT."""
        return self._channel_data["INSTRUMENT"]