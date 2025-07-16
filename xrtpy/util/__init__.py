from xrtpy.util.filename2repo_path import filename2repo_path
from xrtpy.util.make_exposure_map import make_exposure_map
from xrtpy.util.time import epoch

from .filters import solve_filter_name, validate_and_format_filters

__all__ = [
    "epoch",
    "filename2repo_path",
    "make_exposure_map",
    "SSW_MIRRORS",
    "solve_filter_name",
    "validate_and_format_filters"
]

SSW_MIRRORS = (
    "https://sohoftp.nascom.nasa.gov/solarsoft/",
    "https://hesperia.gsfc.nasa.gov/ssw/",
)
