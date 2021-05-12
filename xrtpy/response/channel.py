__all__ = ["Channel"]

import pkg_resources
import sunpy.io.special

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

    def __init__(self,name):
        if name in _channel_name_to_index_mapping:
            self._channel_index = _channel_name_to_index_mapping[name]
            self._channel_data = _genx_file[self._channel_index]
        else:
            raise ValueError(f"{name} is not a valid channel.")

    @property
    def name(self):
        """
        Name of XRT x-ray channels
        """
        return self._channel_data["NAME"]

    @property
    @u.quantity_input
    def wavelength(self)-> u.angstrom:
        """
        Array of wavelengths for every x-ray channel in Angstroms.
        """
        #for i range(array[:length]):
            #wave.append(v6_genx_s[0]["WAVE"][i])
        self.wave = v6_genx_s['WAVE']
        return u.Quantity(self.v6_genx_s['WAVE'], u.angstrom)
        #return self.v6_genx_s['WAVE']

    @property
    def trans(self):
        """
        Get channel transmissions
        """
        self.trans = v6_genx_s['TRANS']
        #for i range(0,v6_genx_s[0]["LENGTH"]:
            #TRANS.append(v6_genx_s[0]["TRANS"][i])
        return self.v6_genx_s['TRANS']

    @property
    def length(self):
        self.trans = v6_genx_s['LENGTH']
        """
        Data number length.
        """
        return self.v6_genx_s['LENGTH']