
""" 
Three-way synthetic DEM validation.

Compare Python, IDL, and known synthetic DEMs using a common
emission-measure-weighted discrepancy metric.
"""

import numpy as np


def temperature_bin_widths(log_temperature):
    """Return linear-temperature integration weights, in kelvin."""
    log_temperature = np.asarray(log_temperature, dtype=float)

    if log_temperature.ndim != 1 or len(log_temperature) < 2:
        raise ValueError("Temperature grid must be one-dimensional with >=2 bins.")

    spacing = np.diff(log_temperature)
    if np.any(spacing <= 0):
        raise ValueError("Temperature grid must be strictly increasing.")

    if not np.allclose(spacing, spacing[0], rtol=1e-5, atol=1e-7):
        raise ValueError("Expected a uniformly spaced log-temperature grid.")

    delta_log_temperature = float(np.mean(spacing))

    return (
        10.0**log_temperature
        * np.log(10.0)
        * delta_log_temperature
    )


def dem_discrepancy(dem_a, dem_b, weights, mask=None):
    """Calculate symmetric, EM-weighted DEM discrepancy (range 0 to 2)."""
    a = np.asarray(dem_a, dtype=float)
    b = np.asarray(dem_b, dtype=float)
    weights = np.asarray(weights, dtype=float)

    if a.ndim != 1 or a.shape != b.shape or a.shape != weights.shape:
        raise ValueError("DEM arrays and weights must have identical 1D shapes.")

    if not (
        np.all(np.isfinite(a))
        and np.all(np.isfinite(b))
        and np.all(np.isfinite(weights))
    ):
        raise ValueError("DEM arrays and weights must contain finite values.")

    if np.any(a < 0) or np.any(b < 0) or np.any(weights <= 0):
        raise ValueError("DEM values must be nonnegative and weights positive.")

    if mask is None:
        mask = np.ones(a.shape, dtype=bool)
    else:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != a.shape:
            raise ValueError("Mask must match the DEM shape.")

    if not np.any(mask):
        return np.nan

    numerator = 2.0 * np.sum(weights[mask] * np.abs(a[mask] - b[mask]))
    denominator = np.sum(weights[mask] * (a[mask] + b[mask]))

    if denominator == 0.0:
        return np.nan

    return float(numerator / denominator)




def contribution_band_mask(log_temperature, band_temperatures):
    """Select temperature-grid bins included in the contribution band."""
    log_temperature = np.asarray(log_temperature, dtype=float)
    band_temperatures = np.asarray(band_temperatures, dtype=float)

    if log_temperature.ndim != 1:
        raise ValueError("Temperature grid must be one-dimensional.")

    if band_temperatures.ndim != 1:
        raise ValueError("Band temperatures must be one-dimensional.")

    mask = np.zeros(log_temperature.shape, dtype=bool)

    for temperature in band_temperatures:
        matches = np.isclose(
            log_temperature,
            temperature,
            rtol=0.0,
            atol=1e-5,
        )

        if np.count_nonzero(matches) != 1:
            raise ValueError(
                f"Contribution-band temperature {temperature} "
                "does not match exactly one grid bin."
            )

        mask |= matches

    return mask


def integrated_em_bias(dem_recovered, dem_reference, weights, mask=None):
    """Calculate signed relative integrated emission-measure bias."""
    recovered = np.asarray(dem_recovered, dtype=float)
    reference = np.asarray(dem_reference, dtype=float)
    weights = np.asarray(weights, dtype=float)

    if (
        recovered.ndim != 1
        or recovered.shape != reference.shape
        or recovered.shape != weights.shape
    ):
        raise ValueError("DEM arrays and weights must have identical 1D shapes.")

    if not (
        np.all(np.isfinite(recovered))
        and np.all(np.isfinite(reference))
        and np.all(np.isfinite(weights))
    ):
        raise ValueError("DEM arrays and weights must contain finite values.")

    if (
        np.any(recovered < 0)
        or np.any(reference < 0)
        or np.any(weights <= 0)
    ):
        raise ValueError("DEM values must be nonnegative and weights positive.")

    if mask is None:
        mask = np.ones(recovered.shape, dtype=bool)
    else:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != recovered.shape:
            raise ValueError("Mask must match the DEM shape.")

    if not np.any(mask):
        return np.nan

    recovered_em = np.sum(recovered[mask] * weights[mask])
    reference_em = np.sum(reference[mask] * weights[mask])

    if reference_em == 0.0:
        return np.nan

    return float((recovered_em - reference_em) / reference_em)



def load_contribution_bands(path):
    """Read the existing per-DEM contribution-band temperatures."""
    import csv

    bands = {}
    with path.open(newline="") as file:
        for row in csv.DictReader(file):
            dem_id = int(row["dem_id"])
            if dem_id in bands:
                raise ValueError(f"Duplicate contribution band for DEM {dem_id}")
            bands[dem_id] = [
                float(value) for value in row["temps"].split()
            ]
    return bands


def validate_one_dem(dem_id, data_dir, bands):
    """Compare Python, IDL, and synthetic truth for one DEM."""
    from scipy.io import readsav

    truth_path = data_dir / "synthetic_dems_data" / f"DEM_{dem_id}.txt"
    idl_path = (
        data_dir / "all_IDL_synthetic_dem_data"
        / f"idl_synthetic_dem_idx{dem_id}_allfilters_base.sav"
    )
    python_path = (
        data_dir / "output_data" / "xrtpy_base_dem_results"
        / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"
    )

    for path in (truth_path, idl_path, python_path):
        if not path.is_file():
            raise FileNotFoundError(f"DEM {dem_id}: missing {path}")

    if dem_id not in bands:
        raise ValueError(f"DEM {dem_id}: contribution band missing")

    truth = np.loadtxt(truth_path)
    if truth.ndim != 2 or truth.shape[1] != 2:
        raise ValueError(f"DEM {dem_id}: invalid truth file")

    truth_logt = truth[:, 0]
    truth_dem = truth[:, 1]

    idl = readsav(str(idl_path), python_dict=True)
    idl_logt = np.asarray(idl["logt_out"], dtype=float).ravel()
    idl_dem = np.asarray(idl["dem_out"], dtype=float).ravel()
    idl_chisq = float(np.asarray(idl["chisq"]).ravel()[0])

    with np.load(python_path) as result:
        python_logt = np.asarray(result["logT"], dtype=float).ravel()
        python_dem = np.asarray(result["dem"], dtype=float).ravel()
        python_chisq = float(result["chisq"])

    # Never compare DEMs on mismatched grids.
    for name, grid in [
        ("IDL", idl_logt),
        ("Python", python_logt),
    ]:
        if grid.shape != truth_logt.shape or not np.allclose(
            grid, truth_logt, rtol=0.0, atol=1e-5
        ):
            raise ValueError(
                f"DEM {dem_id}: {name} temperature grid mismatch"
            )

    if not np.all(np.isfinite([idl_chisq, python_chisq])):
        raise ValueError(f"DEM {dem_id}: nonfinite chi-square")

    weights = temperature_bin_widths(truth_logt)
    band_mask = contribution_band_mask(
        truth_logt, bands[dem_id]
    )

    if not np.any(band_mask):
        raise ValueError(f"DEM {dem_id}: empty contribution band")

    pairs = {
        "python_vs_truth": (python_dem, truth_dem),
        "idl_vs_truth": (idl_dem, truth_dem),
        "python_vs_idl": (python_dem, idl_dem),
    }

    row = {
        "dem_id": dem_id,
        "n_temperature_bins": len(truth_logt),
        "n_gray_band_bins": int(band_mask.sum()),
        "python_raw_chisq": python_chisq,
        "idl_raw_chisq": idl_chisq,
    }

    readable_rows = []

    for region_name, mask, label in [
        ("gray_band", band_mask, "Gray band"),
        ("full_grid", None, "Full grid"),
    ]:
        for pair_name, (a, b) in pairs.items():
            distance = dem_discrepancy(a, b, weights, mask)
            bias = integrated_em_bias(a, b, weights, mask)

            if not np.isfinite(distance) or not np.isfinite(bias):
                raise ValueError(
                    f"DEM {dem_id}: undefined {pair_name} "
                    f"metric in {region_name}"
                )

            prefix = f"{pair_name}_{region_name}"
            row[f"{prefix}_distance"] = distance
            row[f"{prefix}_em_bias_percent"] = 100.0 * bias

            readable_rows.append({
                "dem_id": dem_id,
                "comparison": pair_name.replace("_", " "),
                "region": label,
                "dem_distance_0_to_2": distance,
                "em_bias_percent": 100.0 * bias,
                "python_raw_chisq": python_chisq,
                "idl_raw_chisq": idl_chisq,
            })

    return row, readable_rows


def main():
    """Validate all 162 existing synthetic DEM cases."""
    import csv
    from pathlib import Path

    base_dir = Path(__file__).resolve().parent
    data_dir = base_dir / "data"
    output_dir = data_dir / "output_data" / "scientific_validation"
    bands_path = data_dir / "contribution_bands.csv"

    bands = load_contribution_bands(bands_path)

    all_metrics = []
    all_readable = []

    for dem_id in range(162):
        row, readable_rows = validate_one_dem(
            dem_id, data_dir, bands
        )
        all_metrics.append(row)
        all_readable.extend(readable_rows)

        print(
            f"[{dem_id + 1:3d}/162] DEM {dem_id:3d} "
            f"Python chi2={row['python_raw_chisq']:.5g}  "
            f"IDL chi2={row['idl_raw_chisq']:.5g}"
        )

    # Write only after all cases have passed validation.
    output_dir.mkdir(parents=True, exist_ok=True)

    wide_path = output_dir / "synthetic_triplet_metrics.csv"
    readable_path = output_dir / "synthetic_triplet_readable.csv"

    with wide_path.open("w", newline="") as file:
        writer = csv.DictWriter(
            file, fieldnames=list(all_metrics[0].keys())
        )
        writer.writeheader()
        writer.writerows(all_metrics)

    with readable_path.open("w", newline="") as file:
        writer = csv.DictWriter(
            file, fieldnames=list(all_readable[0].keys())
        )
        writer.writeheader()
        writer.writerows(all_readable)

    print("\nScientific validation completed.")
    print(f"DEM cases: {len(all_metrics)}")
    print(f"Comparison rows: {len(all_readable)}")
    print(f"Saved: {wide_path}")
    print(f"Saved: {readable_path}")


if __name__ == "__main__":
    main()
