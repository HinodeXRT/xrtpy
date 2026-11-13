"""
Extract the ">N% contribution" temperature list from each synthetic intensity
file and write a compact CSV the plotting scripts can read.

Each XRT_intensities_DEM_{id}.txt contains a line (which may WRAP onto
continuation lines) like:
    Temperatures with > 2% contribution      6.80000      6.90000      7.00000
          7.40000      7.50000      7.60000
The number of temperatures varies and may have gaps (multi-component DEMs).

Reads by CONTENT (finds the line starting "Temperature"), then absorbs any
following numbers-only continuation lines until the "Filter" header. Does
NOT depend on line position, so the 4-header-line runners are unaffected.

Output: contribution_bands.csv, one row per DEM:
    dem_id, n_temps, temps
`temps` = space-separated qualifying logT values (empty if none/absent).

Usage (from synthetic_dem_validation/):
    python extract_contribution_bands.py
"""

import csv
import re
from pathlib import Path

BASE_DIR = Path(__file__).parent
INTENSITY_DIR = BASE_DIR / "data" / "synthetic_DEM_condition_data"
OUT_CSV = BASE_DIR / "data" / "contribution_bands.csv"

# Matches a float like 5.50000, 6.1, 1.2e3, etc.
FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def extract_temps(txt_path: Path):
    """Return (dem_id, [logT floats]); reads label line + continuation lines."""
    lines = txt_path.read_text().splitlines()

    dem_id = None
    temps = []
    in_contrib = False
    for line in lines:
        stripped = line.strip()
        low = stripped.lower()

        # DEM index line
        if dem_id is None and low.startswith("dem index"):
            m = re.search(r"index\s*:\s*(\d+)", stripped, re.IGNORECASE)
            if m:
                dem_id = int(m.group(1))
            continue

        # Contribution label line (wording may vary)
        if low.startswith("temperature"):
            key = "contribution"
            if key in low:
                idx = low.rfind(key) + len(key)
                tail = stripped[idx:]
            else:
                tail = re.sub(r"^[^0-9]*\d+\s*%?\s*[^0-9.+-]*", "", stripped)
            temps = [float(x) for x in FLOAT_RE.findall(tail)]
            in_contrib = True          # continuation lines may follow
            continue

        # Absorb numbers-only continuation lines; stop at Filter/any letters
        if in_contrib:
            if low.startswith("filter") or re.search(r"[a-zA-Z]", stripped):
                in_contrib = False
            elif stripped:
                temps.extend(float(x) for x in FLOAT_RE.findall(stripped))

    if dem_id is None:
        m = re.search(r"DEM_(\d+)\.txt$", txt_path.name)
        dem_id = int(m.group(1)) if m else -1

    return dem_id, temps


def _dem_id_from_path(p: Path) -> int:
    m = re.search(r"DEM_(\d+)\.txt$", p.name)
    return int(m.group(1)) if m else -1


def main():
    files = sorted(
        INTENSITY_DIR.glob("XRT_intensities_DEM_*.txt"),
        key=_dem_id_from_path,
    )
    if not files:
        raise FileNotFoundError(f"No intensity files in {INTENSITY_DIR}")

    rows = []
    n_with = 0
    for p in files:
        dem_id, temps = extract_temps(p)
        rows.append((dem_id, temps))
        if temps:
            n_with += 1

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dem_id", "n_temps", "temps"])
        for dem_id, temps in sorted(rows):
            temps_str = " ".join(f"{t:.5f}" for t in temps)
            w.writerow([dem_id, len(temps), temps_str])

    print(f"Wrote {OUT_CSV}")
    print(f"{len(rows)} DEMs, {n_with} with a contribution line, "
          f"{len(rows) - n_with} without.")

    print("\nExamples (incl. known multi-line cases):")
    show = {r[0]: r[1] for r in rows}
    for dem_id in [0, 99, 101, 137, 143]:
        if dem_id in show:
            t = show[dem_id]
            print(f"  DEM {dem_id}: {len(t)} temps -> "
                  f"{t[:8]}{' ...' if len(t) > 8 else ''}")


if __name__ == "__main__":
    main()
