"""
Shared utilities for loading IDL DEM .sav files and parsing their filenames.

Filename convention expected:
    xrt_IDL_dem_<DATE>_<Filter1><Intensity1>_<Filter2><Intensity2>_..._.sav

Example:
    xrt_IDL_dem_20071213T0401_Bemed603.875886_Bethin150.921435_Alpoly2412.34_.sav

Compact filter name → XRTpy filter name mapping is handled automatically.
Will be updated once I start testing with MC. 
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import readsav

# Filter name mapping — compact filename token → XRTpy filter name
# Keys are sorted longest-first at parse time to handle compound filters
# (e.g. "AlpolyTipoly" must be matched before "Alpoly")

_FILTER_MAP: dict[str, str] = {
    "AlpolyAlmesh":  "Al-poly/Al-mesh",
    "AlpolyTipoly":  "Al-poly/Ti-poly",
    "AlpolyAlthick": "Al-poly/Al-thick",
    "AlpolyBethick": "Al-poly/Be-thick",
    "CpolyTipoly":   "C-poly/Ti-poly",
    "Bemed":         "Be-med",
    "Bethin":        "Be-thin",
    "Bethick":       "Be-thick",
    "Alpoly":        "Al-poly",
    "Tipoly":        "Ti-poly",
    "Almesh":        "Al-mesh",
    "Almed":         "Al-med",
    "Althick":       "Al-thick",
    "Cpoly":         "C-poly",
}

# Pre-sorted (longest first) for unambiguous prefix matching
_FILTER_KEYS_SORTED = sorted(_FILTER_MAP, key=len, reverse=True)



# Data classes
@dataclass(frozen=True)
class SavCase:
    """
    All information needed to run one IDL vs XRTpy comparison.

    Attributes
    ----------
    sav_path : Path
        Path to the IDL .sav file.
    observation_date : str
        ISO-8601 observation date, e.g. ``"2008-01-04T11:04:26"``.
    filters : list[str]
        XRTpy filter names, in the same order as intensities.
    intensities : list[float]
        Observed intensities in DN/s/pix, matching filters.
    label : str
        Short human-readable label derived from the filename.
    """
    sav_path: Path
    observation_date: str
    filters: list[str]
    intensities: list[float]
    label: str

    @property
    def intensities_array(self) -> np.ndarray:
        return np.array(self.intensities, dtype=float)


@dataclass(frozen=True)
class IDLResult:
    """Loaded IDL DEM solution."""
    logT: np.ndarray   # (nT,) log10 K
    dem:  np.ndarray   # (nT,) cm^-5 K^-1


# Filename parser
def _parse_filter_token(token: str) -> tuple[str, float] | tuple[None, None]:
    """
    Parse a single token like ``'Bemed603.875886'`` into
    ``('Be-med', 603.875886)``.  Returns ``(None, None)`` on failure.
    """
    for key in _FILTER_KEYS_SORTED:
        if token.startswith(key):
            remainder = token[len(key):]
            try:
                return _FILTER_MAP[key], float(remainder)
            except ValueError:
                continue
    return None, None


def parse_sav_filename(sav_path: str | Path) -> SavCase:
    """
    Parse an IDL DEM .sav filename and return a :class:`SavCase`.

    The expected filename format is::

        xrt_IDL_dem_<YYYYMMDDTHHMI>_<Token1>_<Token2>_..._.sav

    where each token is a compact filter name immediately followed by its
    intensity value (e.g. ``Bemed603.875886``).

    Parameters
    ----------
    sav_path : str or Path
        Full path to the .sav file.

    Returns
    -------
    SavCase
    """
    sav_path = Path(sav_path)
    stem = sav_path.stem  # filename without .sav

    # Strip the fixed prefix
    base = stem.removeprefix("xrt_IDL_dem_").strip("_")

    parts = [p for p in base.split("_") if p]
    if not parts:
        raise ValueError(f"Cannot parse filename: {sav_path.name}")

    # --- Date ---
    date_token = parts[0]
    m = re.match(r"(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})", date_token)
    if m:
        obs_date = (
            f"{m.group(1)}-{m.group(2)}-{m.group(3)}"
            f"T{m.group(4)}:{m.group(5)}"
        )
    else:
        obs_date = date_token  # fall back to raw string

    # --- Filter / intensity pairs ---
    filters: list[str] = []
    intensities: list[float] = []

    for token in parts[1:]:
        fname, val = _parse_filter_token(token)
        if fname is not None:
            filters.append(fname)
            intensities.append(val)

    if not filters:
        raise ValueError(
            f"No filter/intensity pairs found in filename: {sav_path.name}"
        )

    label = date_token  # e.g. "20071213T0401"

    return SavCase(
        sav_path=sav_path,
        observation_date=obs_date,
        filters=filters,
        intensities=intensities,
        label=label,
    )


# SAV loade
def load_idl_sav(sav_path: str | Path) -> IDLResult:
    """
    Load an IDL DEM .sav file and return ``(logT, dem_base)`` as 1-D arrays.

    Tries common key name variants used across different IDL save scripts.

    Parameters
    ----------
    sav_path : str or Path

    Returns
    -------
    IDLResult
    """
    sav_path = Path(sav_path)
    data = readsav(str(sav_path), python_dict=True)

    # logT
    logT = None
    for key in ("logt", "logT", "logT_out", "logt_out", "logT_mc"):
        if key in data:
            logT = np.array(data[key]).ravel().astype(float)
            break
    if logT is None:
        raise KeyError(
            f"logT key not found in {sav_path.name}. "
            f"Available keys: {sorted(data.keys())}"
        )

    # DEM base
    dem = None
    for key in ("dem", "dem_out", "dem0", "dem_mc"):
        if key in data:
            arr = np.array(data[key])
            if arr.ndim == 2:
                arr = arr[:, 0]   # first column = base DEM
            dem = arr.ravel().astype(float)
            break
    if dem is None:
        raise KeyError(
            f"DEM key not found in {sav_path.name}. "
            f"Available keys: {sorted(data.keys())}"
        )

    return IDLResult(logT=logT, dem=dem)


# Convenience: discover all .sav files in a directory
def discover_cases(data_dir: str | Path) -> list[SavCase]:
    """
    Find all ``xrt_IDL_dem_*.sav`` files in *data_dir* and parse each one.

    Parameters
    ----------
    data_dir : str or Path
        Directory to search (non-recursive).

    Returns
    -------
    list[SavCase]
        One entry per valid .sav file found.
    """
    data_dir = Path(data_dir)
    cases: list[SavCase] = []
    for sav_file in sorted(data_dir.glob("xrt_IDL_dem_*.sav")):
        try:
            cases.append(parse_sav_filename(sav_file))
        except ValueError as exc:
            print(f"  Warning: skipping {sav_file.name} — {exc}")
    return cases
