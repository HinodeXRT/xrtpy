from xrtpy.util.filename2repo_path import filename2repo_path
from xrtpy.util.make_exposure_map import make_exposure_map
from xrtpy.util.time import epoch

__all__ = [
    "epoch",
    "filename2repo_path",
    "make_exposure_map",
    "SSW_MIRRORS",
]

SSW_MIRRORS = (
    "https://sohoftp.nascom.nasa.gov/solarsoft/",
    "https://hesperia.gsfc.nasa.gov/ssw/",
)
