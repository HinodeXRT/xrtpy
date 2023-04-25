"""
Functionality for subtracting light leak (visible stray light) image from XRT synoptic composite images.
"""
__all__ = ["lightleak_correction"]

from pathlib import Path

_data_file_path = {
    "lightleak_calibration_data": Path(__file__).parent.absolute()
    / "data/lightleak_calibration_data"
}


class lightleak_correction:
    """Functionality for subtracting light leak (visible stray light) image from XRT synoptic composite images."""

    def __init__(self, data_file, in_idx, in_da, abundance_model="coronal"):
        self._data_file = data_file
        self._in_idx = in_idx  # Fact check with new function
        self._in_da = in_da

    @property
    def index_data(self):
        """This is the index for a synoptic composite image."""
        return self.self._in_idx

    @property
    def index_data_array(self):
        """This is the 2D data array for the synoptic composite"""
        return self._in_idx
