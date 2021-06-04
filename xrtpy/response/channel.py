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
        Entrance filters name.
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
    def entrancefilter_density(self)-> u.g * u.cm**3 :
        """
        XRT channel filter material density. g cm^-3
        """
        return u.Quantity(self._en_filter_data['DENS'],u.g * u.cm *u.cm * u.cm)

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

class Mirror_1:
    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._mirror1_data = self._genx_file[self._channel_index]["MIRROR1"]

    @property
    def mirror1_name(self,) -> str:
        """
        Hinode/XRT flight model mirror-1.
        """
        return self._mirror1_data["LONG_NAME"]

    @property
    def mirror1_material(self) -> str:
        """
        XRT flight model mirror-1 material.
        """
        return self._mirror1_data['MATERIAL']

    @property
    @u.quantity_input
    def mirror1_density(self)-> u.g * u.cm**3: 
        """
        Mirror-1 density.
        """
        return u.Quantity(self._mirror1_data['DENS'], u.g * u.cm**3)
        
    @property
    @u.quantity_input
    def mirror1_graze_angle(self,)-> u.deg: 
        """
        Mirror-1 graze angle in units of degrees.
        """
        return u.Quantity(self._mirror1_data['GRAZE_ANGLE'], u.deg)

    @property
    @u.quantity_input
    def mirror1_wavelength(self,) -> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._mirror1_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

    @property
    @u.quantity_input
    def mirror1_reflection1(self,) -> u.angstrom:
        """
        Reflection of miiror-1.
        """
        return u.Quantity(self._mirror1_data['REFL'], u.angstrom)[:self.number_of_wavelengths]

    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._mirror1_data['LENGTH']

class Mirror_2:
    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._mirror2_data = self._genx_file[self._channel_index]["MIRROR2"]

    @property
    def mirror2_name(self,) -> str:
        """
        Hinode/XRT Flight Model Mirror-1
        """
        return self._mirror2_data["LONG_NAME"]

    @property
    def mirror2_material(self) -> str:
        """
        XRT flight model mirror-1 material
        """
        return self._mirror2_data['MATERIAL']

    @property
    @u.quantity_input
    def mirror2_density(self)-> u.g * u.cm**3:
        """
        Mirror-1 density.
        """
        return u.Quantity(self._mirror2_data['DENS'], u.g * u.cm**3) 

    @property
    @u.quantity_input
    def mirror2_graze_angle(self)-> u.deg:
        """
        Mirror-1 graze angle-degrees.
        """
        return u.Quantity(self._mirror2_data['GRAZE_ANGLE'], u.deg)

    @property
    @u.quantity_input
    def mirror2_wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        return u.Quantity(self._mirror2_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]

    @property
    def mirror2_reflection(self):
        """
        Mirror-2 reflection.
        """
        return self._mirror2_data['REFL'][:self.number_of_wavelengths]
    
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._mirror2_data['LENGTH']

class Filter1:
    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._fp_filter1_data = self._genx_file[self._channel_index]["FP_FILTER1"]

    @property
    def filter1_name(self,) -> str:
        """
        XRT channel Filter-1 Position
        """
        return self._fp_filter1_data["LONG_NAME"]

    @property
    def filter1_material(self) -> str:
        """
        XRT channel filter-1 material.
        """
        return self._fp_filter1_data['MATERIAL']

    @property  
    @u.quantity_input
    def filter1_thickness(self)-> u.angstrom:
        """
        XRT channel filter-1 thickness.
        """
        return u.Quantity(self._fp_filter1_data['THICK'], u.angstrom)

    @property
    @u.quantity_input
    def filter1_density(self)-> u.g * u.cm**3:
        """
        XRT channel filter-1 fensity.
        """
        return u.Quantity(self._fp_filter1_data['DENS'], u.g * u.cm**3)

    @property
    @u.quantity_input
    def filter1_wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel-Filter1 in Angstroms.
        """
        return u.Quantity(self._fp_filter1_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]
    
    @property
    def filter1_transmission(self):
        """
        Get channel transmission filter-1.
        """
        return self._fp_filter1_data['TRANS'][:self.number_of_wavelengths] 
      
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._fp_filter1_data['LENGTH']
    
    @property
    def filter1_mesh_trans(self):
        """
        Mech Trans channel filter-1.
        """
        return self._fp_filter1_data['MESH_TRANS']

    @property
    def filter1_substrate(self) -> str:
        """
        XRT substrate channel filter-1.
        """
        return self._fp_filter1_data['SUBSTRATE']

class Filter2:
    _genx_file = _genx_file

    def __init__(self, index):
        self._channel_index = index
        self._fp_filter2_data = self._genx_file[self._channel_index]["FP_FILTER2"]

    @property
    def filter2_name(self,) -> str:
        """
        XRT channel Filter-2 position.
        """
        return self._fp_filter2_data["LONG_NAME"]

    @property
    def filter2_material(self) -> str:
        """
        XRT channel filter-2 material.
        """
        return self._fp_filter2_data['MATERIAL']

    @property 
    @u.quantity_input 
    def filter2_thickness(self)-> u.angstrom:
        """
        XRT channel filter-2 thickness.
        """
        return u.Quantity(self._fp_filter2_data['THICK'], u.angstrom)

    @property 
    @u.quantity_input
    def filter2_density(self)-> u.g * u.cm**3:
        """
        XRT channel filter-2 density.
        """
        return u.Quantity(self._fp_filter2_data['DENS'], u.g * u.cm**3)

    @property
    @u.quantity_input
    def filter2_wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel Filter-2 in Angstroms.
        """
        return u.Quantity(self._fp_filter2_data['WAVE'], u.angstrom)[:self.number_of_wavelengths]
    
    @property
    def filter2_transmission(self):
        """
        Get channel transmission filter-2.
        """
        return self._fp_filter2_data['TRANS'][:self.number_of_wavelengths] 
      
    @property
    def number_of_wavelengths(self):
        """
        Data number length.
        """
        return self._fp_filter2_data['LENGTH']
    
    @property
    def filter2_mesh_trans(self):
        """
        Mech Trans channel filter-2.
        """
        return self._fp_filter2_data['MESH_TRANS']

    @property
    def filter2_substrate(self) -> str:
        """
        XRT substrate channel filter-2.
        """
        return self._fp_filter2_data['SUBSTRATE']

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
    """
    _genx_file = _genx_file

    def __init__(self,name):
        if name in _channel_name_to_index_mapping:
            self._channel_index = _channel_name_to_index_mapping[name]
            self._channel_data = _genx_file[self._channel_index] 
            self._geometry = Geometry(self._channel_index)
            self._entrancefilter = EntranceFilter(self._channel_index)
            self._mirror_1 = Mirror_1(self._channel_index)
            self._mirror_2 = Mirror_2(self._channel_index)
            self._filter1 = Filter1(self._channel_index)
            self._filter2 = Filter2(self._channel_index)
            self._ccd = CCD(self._channel_index)

        else:
            raise ValueError(f"{name} is not a valid channel.")
    
    @property
    def geometry(self)-> Geometry:
        return self._geometry

    @property
    def entrancefilter(self)-> EntranceFilter:
        return self._entrancefilter

    @property
    def mirror_1(self)-> Mirror_1:
        return self._mirror_1

    @property
    def mirror_2(self)-> Mirror_2:
        return self._mirror_2

    @property
    def filter1(self)-> Filter1:
        return self._filter1
    
    @property
    def filter2(self)-> Filter2:
        return self._filter2

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
