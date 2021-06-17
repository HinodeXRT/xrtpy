__all__ = ["channel","Geometry","EntranceFilter","Mirror","Filter","CCD","Channel","resolve_filter_name"]

import pkg_resources
import sunpy.io.special
import numpy as np

from astropy import units as u

filename = pkg_resources.resource_filename( "xrtpy", "data/channels/xrt_channels_v0016.genx")

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

def resolve_filter_name(name):
    name = name.replace("_","-")
    parts: list = name.split("/")
    new_parts: list = [part.capitalize() for part in name.split("/")]
    name: str = "/".join(new_parts)
    return(name)

class Geometry:
    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._geom_data = self._genx_file[self._channel_index]["GEOM"]

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

    def __init__(self, index):
        self._channel_index = index
        self._en_filter_data = self._genx_file[self._channel_index]["EN_FILTER"]

    @property
    def entrancefilter_name(self,) -> str:
        """
        Entrance filter name.
        """
        return self._en_filter_data["LONG_NAME"]

    @property
    def entrancefilter_material(self) -> str:
        """
        XRT channel filter material.
        """
        return self._en_filter_data['MATERIAL']

    @property 
    @u.quantity_input
    def entrancefilter_thickness(self) -> u.angstrom:
        """
        XRT channel filter material thickness.
        """
        return u.Quantity(self._en_filter_data['THICK'],u.angstrom)
    
    @property
    def entrancefilter_density(self)-> u.g * u.cm**-3 :
        """
        XRT channel filter material density. g cm^-3
        """
        return u.Quantity(self._en_filter_data['DENS'],u.g * u.cm**-3)

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._en_filter_data['LENGTH']

    @property
    @u.quantity_input
    def entrancefilter_wavelength(self) -> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._en_filter_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

    @property 
    def entrancefilter_transmission(self):
        """
        Entrance filter channel transmission.
        """
        return self._en_filter_data['TRANS'][:self.number_of_wavelengths]

    @property
    def entrancefilter_mesh_transmission(self):
        """
        Transmission of mech filter substrate.
        """
        return self._en_filter_data['MESH_TRANS']

    @property
    def entrancefilter_substrate(self) -> str:
        """
        XRT channel filter substrate
        """
        return self._en_filter_data['SUBSTRATE']

class Mirror:
    _genx_file = _genx_file

    def __init__(self, index,mirror_number):
        self._channel_index = index
        self._mirror_data = self._genx_file[self._channel_index][f"MIRROR{mirror_number}"]

    @property
    def mirror_name(self,) -> str:
        """
        Hinode/XRT flight model mirror.
        """
        return self._mirror_data["LONG_NAME"]

    @property
    def mirror_material(self) -> str:
        """
        XRT flight model mirror material.
        """
        return self._mirror_data['MATERIAL']

    @property
    @u.quantity_input
    def mirror_density(self)-> u.g * u.cm**-3: 
        """
        Mirror density.
        """
        return u.Quantity(self._mirror_data['DENS'], u.g * u.cm**-3)
        
    @property
    @u.quantity_input
    def mirror_graze_angle(self,)-> u.deg: 
        """
        Mirror graze angle in units of degrees.
        """
        return u.Quantity(self._mirror_data['GRAZE_ANGLE'], u.deg)

    @property
    @u.quantity_input
    def mirror_wavelength(self,) -> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._mirror_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

    @property
    @u.quantity_input
    def mirror_reflection1(self,) -> u.angstrom:
        """
        Reflection of miiror.
        """
        return u.Quantity(self._mirror_data['REFL'], u.angstrom)[:self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._mirror_data['LENGTH']

class Filter:
    _genx_file = _genx_file

    def __init__(self, index, filter_number):
        self._channel_index = index
        self._fp_filter_data = self._genx_file[self._channel_index][f"FP_FILTER{filter_number}"]

    @property
    def filter_name(self,) -> str:
        """
        XRT channel Filter Position
        """
        return self._fp_filter_data["LONG_NAME"]

    @property
    def filter_material(self) -> str:
        """
        XRT channel filter material.
        """
        return self._fp_filter_data['MATERIAL']

    @property  
    @u.quantity_input
    def filter_thickness(self)-> u.angstrom:
        """
        XRT channel filter thickness.
        """
        return u.Quantity(self._fp_filter_data['THICK'], u.angstrom)

    @property
    @u.quantity_input
    def filter_density(self)-> u.g * u.cm**-3:
        """
        XRT channel filter fensity.
        """
        return u.Quantity(self._fp_filter_data['DENS'], u.g * u.cm**-3)

    @property
    @u.quantity_input
    def filter_wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel-Filter in Angstroms.
        """
        return u.Quantity(self._fp_filter_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]
    
    @property
    def filter_transmission(self):
        """
        Get channel transmission filter.
        """
        return self._fp_filter_data['TRANS'][:self.number_of_wavelengths] 
      
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._fp_filter_data['LENGTH']
    
    @property
    def filter_mesh_trans(self):
        """
        Mech Trans channel filter.
        """
        return self._fp_filter_data['MESH_TRANS']

    @property
    def filter_substrate(self) -> str:
        """
        XRT substrate channel filter.
        """
        return self._fp_filter_data['SUBSTRATE']

class CCD:
    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._ccd_data =self._genx_file[self._channel_index]["CCD"]

    @property
    def ccd_name(self,) -> str:
        """
        Hinode/XRT flight model CCD.
        """
        return self._ccd_data["LONG_NAME"]

    @property
    @u.quantity_input
    def ccd_ev_ore_electron(self,)->u.eV/u.electron:
        """
        The number of eVs necessary to dislodge one electron.
        """
        return u.Quantity(self._ccd_data['EV_PER_EL'],u.eV/u.electron)

    @property
    @u.quantity_input
    def ccd_full_well(self,)->u.electron: 
        """
        Units electrons.
        """
        return u.Quantity(self._ccd_data['FULL_WELL'],u.electron)

    @property
    @u.quantity_input
    def ccd_gain_left(self,)->u.electron: 
        """
        Gain when reading the left port of the CCD
        """
        return u.Quantity(self._ccd_data['GAIN_L'],u.electron)

    @property
    @u.quantity_input
    def ccd_gain_right(self,)->u.electron: 
        """
        Gain when reading the right port of the CCD
        """
        return u.Quantity(self._ccd_data['GAIN_R'],u.electron)

    @property
    def ccd_quantum_efficiency(self):
        """
        Quantum efficiency.
        """
        return self._ccd_data['QE'][:self.number_of_wavelengths] 

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._ccd_data['LENGTH']

    @property
    @u.quantity_input
    def ccd_pixel_size(self)->u.micron:
        """
        CCD pixel size.
        """
        return u.Quantity(self._ccd_data['PIXEL_SIZE'],u.micron)

    @property
    @u.quantity_input
    def ccd_wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._ccd_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

class Channel:
    """
    XRTpy

    Available channels: ``"Al-mesh"``, ``"Al-poly"``,  ``"C-poly"``, ``"Ti-poly"``, ``"Be-thin"``, ``"Be-med"``, ``"Al-med"``, ``"Al-thick"``,  ``"Be-thick"`` ,
    ``"Al-poly/Al-mesh"``, ``"Al-poly/Ti-poly"``, ``"Al-poly/Al-thick"``, ``"Al-poly/Be-thick"`` , ``"C-poly/Ti-poly"`` 
    """
    _genx_file = _genx_file

    def __init__(self,name):
        name = resolve_filter_name(name)
        if name in _channel_name_to_index_mapping:
            self._channel_index = _channel_name_to_index_mapping[name]
            self._channel_data = _genx_file[self._channel_index] 
            self._geometry = Geometry(self._channel_index)
            self._entrancefilter = EntranceFilter(self._channel_index)
            self._mirror_1 = Mirror(self._channel_index,1)
            self._mirror_2 = Mirror(self._channel_index,2)
            self._filter_1 = Filter(self._channel_index,1)
            self._filter_2 = Filter(self._channel_index,2)
            self._ccd = CCD(self._channel_index)
        else:
            raise ValueError(f"{name} is not a valid channel. The available channels are: {list(_channel_name_to_index_mapping.keys())}")
    
    @property
    def geometry(self)-> Geometry:
        return self._geometry

    @property
    def entrancefilter(self)-> EntranceFilter:
        return self._entrancefilter

    @property
    def mirror_1(self)-> Mirror:
        return self._mirror_1

    @property
    def mirror_2(self)-> Mirror:
        return self._mirror_2

    @property
    def filter_1(self)-> Filter:
        return self._filter_1
    
    @property
    def filter_2(self)-> Filter:
        return self._filter_2

    @property
    def ccd(self)-> CCD:
        return self._ccd

    def __str__(self):
        """Reable printout."""
        return f"XRT Channel for {self.name}"
    
    def __repr__(self):
        """Code representation."""
        return f"Channel({repr(self.name)})"

    @property
    def name(self) -> str:
        """
        Name of XRT X-Ray channels.
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
