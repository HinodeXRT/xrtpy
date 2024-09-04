from xrtpy.util import filename2repo_path, make_exposure_map, time
from xrtpy.util.time import epoch

__all__ = [
    "epoch",
    "time",
    "filename2repo_path",
    "make_exposure_map",
    "SSW_MIRRORS",
]

SSW_MIRRORS = (
    "https://sohoftp.nascom.nasa.gov/solarsoft/",
    "https://hesperia.gsfc.nasa.gov/ssw/",
)
