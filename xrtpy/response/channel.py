__all__ = ["Channel","channel","Geometry","EntranceFilter","Mirror_1","Mirror_2","Filter1","Filter2","CCD"]

import pkg_resources
import sunpy.io.special
import numpy as np

from astropy import units as u

filename = pkg_resources.resource_filename(
    "xrtpy", "data/channels/xrt_channels_v0016.genx"
)

_channel_name_to_index_mapping = {
    "Al-mesh": 0, 
    "Al-poly": 1, 
    "C-poly": 2 ,
    "Ti-poly": 3,
    "Be-thin": 4,
    "Be-med": 5,
    "Al-med": 6,
    "Al-thick": 7, 
    "Be-thick": 8 ,
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10, 
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12 , 
    "C-poly/Ti-poly": 13, 
}

_genx_file = sunpy.io.special.genx.read_genx(filename)['SAVEGEN0']


class Channel:
    """
    XRTpy
    """
    _genx_file = _genx_file

    def __init__(self,name):
        if name in _channel_name_to_index_mapping:
            self._channel_index = _channel_name_to_index_mapping[name]
            self._channel_data = _genx_file[self._channel_index]
        else:
            raise ValueError(f"{name} is not a valid channel.")

    def __str__(self):
        """Reable printout."""
        return f"XRT Channel for {self.name}"
    
    def __repr__(self):
        """Code representation."""
        return f"Channel({repr(self.name)})"

    @property
    def name(self) -> str:
        """
        Name of XRT x-ray channels.
        """
        return self._channel_data["NAME"]

    @property
    @u.quantity_input
    def wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._channel_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]
    
    @property
    def transmission(self):
        """
        Get channel transmission.
        """
        return self._channel_data['TRANS'][:self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._channel_data['LENGTH']
       
    @property
    def observatory(self) -> str:
        """
        Spacecraft: Hinode.
        """
        return self._channel_data['OBSERVATORY']
        
    @property
    def instrument(self) -> str:
        """
        X-Ray Telescope - XRT
        """
        return self._channel_data['INSTRUMENT']


class Geometry:
    _genx_file = _genx_file

    _geom_data = _genx_file[0]["GEOM"]

    def __init__(self, index=0):
        self._index = index

    @property
    def name(self,) -> str:
        """
        Hinode/XRT flight model Geometry.
        """
        return self._geom_data["LONG_NAME"]

    @property
    @u.quantity_input
    def focal_len(self,)-> u.cm:  
        """
        Hinode/XRT flight model geometry focual length.
        """
        return u.Quantity(self._geom_data['FOC_LEN'], u.cm)

    @property
    @u.quantity_input
    def aperture_area(self)-> u.cm**2:
        """
        Hinode/XRT flight model geometry aperture area.
        """
        return u.Quantity(self._geom_data['APERTURE_AREA'], u.cm**2)

class EntranceFilter:
    _genx_file = _genx_file

    _en_filter_data = _genx_file[0]["EN_FILTER"]

    def __init__(self, index=0):
        self._index = index

    @property
    def name(self,) -> str:
        """
        Entrance filters name.
        """
        return self._en_filter_data["LONG_NAME"]

    @property
    def filter_material(self) -> str:
        """
        XRT channel filter material.
        """
        return self._en_filter_data[0]['MATERIAL']

    @property
    def filter_thickness(self):
        """
        XRT channel filter material thickness.
        """
        return self._en_filter_data[0]['THICK']
    
    @property
    def filter_density(self):
        """
        XRT channel filter material density.
        """
        return self._en_filter_data[0]['DENS']

    @property
    @u.quantity_input
    def wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._en_filter_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

    @property
    def transmission(self):
        """
        Entrance filter channel transmission.
        """
        return self._en_filter_data['TRANS'][:self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._en_filter_data[0]['LENGTH']

    @property
    def mesh_trans(self):
        """
        Mech Trans
        """
        return self._en_filter_data['MESH_TRANS']

    @property
    def substrate(self) -> str:
        """
        XRT subtrate.
        """
        return self._en_filter_data['SUBSTRATE']

class Mirror_1:
    _genx_file = _genx_file
    
    _mirror1_data = _genx_file[0]["MIRROR1"]

    def __init__(self, index=0):
        self._index = index

    @property
    def name(self,) -> str:
        """
        Hinode/XRT flight model mirror-1.
        """
        return self._mirror1_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """
        XRT flight model mirror-1 material.
        """
        return self._mirror1_data['MATERIAL']

    @property
    @u.quantity_input
    def density(self)-> u.g * u.cm**3: 
        """
        Mirror-1 density.
        """
        return u.Quantity(self._mirror1_data['DENS'], u.g * u.cm**3)
        
    @property
    @u.quantity_input
    def graze_angle(self)-> u.deg: 
        """
        Mirror-1 graze angle in units of degrees.
        """
        return u.Quantity(self._mirror1_data['GRAZE_ANGLE'], u.deg)

    @property
    @u.quantity_input
    def wavelength(self) -> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._mirror1_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

    @property
    def reflection(self):
        """
        Miiror-1 reflection.
        """
        return self._mirror1_data['REFL']
    
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._genx_file[0]['LENGTH']


class Mirror_2:
    _genx_file = _genx_file

    _mirror2_data = _genx_file[0]["MIRROR2"]

    def __init__(self, index=0):
        self._index = index

    @property
    def name(self,) -> str:
        """
        Hinode/XRT Flight Model Mirror-1
        """
        return self._mirror2_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """
        XRT flight model mirror-1 material
        """
        return self._mirror2_data['MATERIAL']

    @property
    @u.quantity_input
    def density(self)-> u.g * u.cm**3:
        """
        Mirror-1 density.
        """
        return u.Quantity(self._mirror2_data['DENS'], u.g * u.cm**3) 

    @property
    @u.quantity_input
    def graze_angle(self)-> u.deg:
        """
        Mirror-1 graze angle-degrees.
        """
        return u.Quantity(self._mirror2_data['GRAZE_ANGLE'], u.deg)

    @property
    @u.quantity_input
    def wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._mirror2_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

    @property
    def reflection(self):
        """
        Mirror-2 reflection.
        """
        return self._mirror2_data['REFL']
    
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._genx_file[0]['LENGTH']

class Filter1:
    _genx_file = _genx_file

    _fp_filter1_data = _genx_file[0]["FP_FILTER1"]

    def __init__(self, index=0):
        self._index = index

    @property
    def name(self,) -> str:
        """
        XRT channel Filter-1 Position
        """
        return self._fp_filter1_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """
        XRT channel Filter-1 Material
        """
        return self._fp_filter1_data['MATERIAL']

    @property  
    def filter1_thickness(self):
        """
        XRT channel Filter-1 Thickness: Unts? A?
        """
        return self._fp_filter1_data[0]['THICK']

    @property
    @u.quantity_input
    def filter_density(self)-> u.g * u.cm**3:
        """
        XRT channel Filter1 Density:g cm^-3
        """
        return u.Quantity(self._fp_filter1_data[0]['DENS'], u.g * u.cm**3)

    @property
    @u.quantity_input
    def wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel-Filter1 in Angstroms.
        """
        return u.Quantity(self._fp_filter1_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]
    
    @property
    def transmission(self):
        """
        Get channel transmission filter-1.
        """
        return self._fp_filter1_data['TRANS'][:self.number_of_wavelengths] 
      
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._fp_filter1_data[0]['LENGTH']
    
    @property
    def mesh_trans(self):
        """
        Mech Trans channel filter-1.
        """
        return self._fp_filter1_data['MESH_TRANS']

    @property
    def substrate(self) -> str:
        """
        XRT substrate channel filter-1.
        """
        return self._fp_filter1_data['SUBSTRATE']


class Filter2:

    _genx_file = _genx_file

    _fp_filter2_data = _genx_file[0]["FP_FILTER2"]

    def __init__(self, index=0):
        self._index = index

    @property
    def name(self,) -> str:
        """
        XRT channel Filter-2 position.
        """
        return self._fp_filter2_data["LONG_NAME"]

    @property
    def material(self) -> str:
        """
        XRT channel filter-2 material.
        """
        return self._fp_filter2_data['MATERIAL']

    @property  
    def filter1_thickness(self): 
        """
        XRT channel filter-2 thickness.
        """
        return self._fp_filter2_data[0]['THICK']

    @property
    def filter_density(self):
        """
        XRT channel filter-2 density.
        """
        return self._fp_filter2_data[0]['DENS']

    @property
    @u.quantity_input
    def wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel Filter-2 in Angstroms.
        """
        return u.Quantity(self._fp_filter2_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]
    
    @property
    def transmission(self):
        """
        Get channel transmission filter-2.
        """
        return self._fp_filter2_data['TRANS'][:self.number_of_wavelengths] 
      
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._fp_filter2_data[0]['LENGTH']
    
    @property
    def mesh_trans(self):
        """
        Mech Trans channel filter-2.
        """
        return self._fp_filter2_data['MESH_TRANS']

    @property
    def substrate(self) -> str:
        """
        XRT substrate channel filter-2.
        """
        return self._fp_filter2_data['SUBSTRATE']


class CCD:
    _genx_file = _genx_file
    
    _ccd_data = _genx_file[0]["CCD"]

    def __init__(self, index=0):
        self._index = index

    @property
    def name(self,) -> str:
        """
        Hinode/XRT flight model CCD.
        """
        return self._ccd_data["LONG_NAME"]

    @property
    @u.quantity_input
    def ev_ore_el(self,)->u.eV:
        """
        XRT EV_PER_EL.
        """
        return u.Quantity(self._ccd_data['EV_PER_EL'],u.eV)

    @property
    @u.quantity_input
    def full_well(self,)->u.electron: 
        """
        Units electrons.
        """
        return u.Quantity(self._ccd_data['FULL_WELL'],u.electron)

    @property
    def gain_l(self):
        """
        UNITS el per 12 bit ADU.
        """
        return self._ccd_data['GAIN_L']

    @property
    def gain_r(self):
        """
        UNITS el per 12 bit ADU.
        """
        return self._ccd_data['GAIN_R']

    @property
    def quantum_efficiency(self):
        """
        Quantum efficiency.
        """
        return self._ccd_data['QE'][:self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._genx_file[0]['LENGTH']
    
    @property
    @u.quantity_input
    def ev_pre_el(self)-> u.eV:
        """
        Units eV
        """
        return u.Quantity(self._ccd_data['EV_PRE_EL'], u.eV)

    @property
    @u.quantity_input
    def pixel_size(self)->u.micron:
        """
        CCD pixel size.
        """
        return u.Quantity(self._ccd_data['PIXEL_SIZE'],u.micron)

    @property
    @u.quantity_input
    def wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._ccd_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]
