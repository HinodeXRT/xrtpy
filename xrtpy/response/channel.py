__all__ = ["Channel"]

import pkg_resources
import sunpy.io.special

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

_genx_file = sunpy.io.special.genx.read_genx(filename)["SAVEGEN0"]


class Channel:
    """
    XRTpy
    """

    def __init__(self, name):
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
