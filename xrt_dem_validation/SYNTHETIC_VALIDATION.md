# Synthetic DEM Validation: Python, IDL, and Known Truth

## Purpose

This study compares the existing XRTpy and IDL base-solution differential emission measure (DEM) inversions against 162 synthetic DEMs with known truth.

The analysis addresses three distinct questions:

1. **Python recovery:** How closely does the XRTpy base result reproduce the synthetic truth?
2. **IDL recovery:** How closely does the IDL base result reproduce the same truth?
3. **Implementation agreement:** How closely do the XRTpy and IDL base results agree with each other?

The analysis uses existing standalone base-solution outputs. No inversion is rerun, and the numerical solver is not modified.

## Inputs and temperature regions

All comparisons use 26 temperature bins on a common logarithmic grid, log₁₀(T/K) = 5.5 to 8.0, with a spacing of 0.1 dex.

Inputs are obtained from:

- `data/synthetic_dems_data/DEM_{id}.txt` — synthetic truth.
- `data/output_data/xrtpy_base_dem_results/xrtpy_dem_idx{id}_allfilters_base.npz` — Python standalone base solution.
- `data/all_IDL_synthetic_dem_data/idl_synthetic_dem_idx{id}_allfilters_base.sav` — IDL standalone base solution.
- `data/contribution_bands.csv` — temperature-bin selections for the contribution bands.

Two regions are evaluated:

- **Gray contribution band (primary):** Only the log-temperature bin centers specified for each DEM in `contribution_bands.csv`.
- **Full grid (secondary):** All 26 temperature bins.

Disconnected contribution-band intervals are supported. The supplied band selection is accepted as a fixed input to this validation; its underlying physical derivation is not independently assessed here.

Temperature grids are checked before DEM arrays are compared. The implementation does not silently interpolate or reorder mismatched arrays.

## Quantitative metrics

### Symmetric DEM discrepancy

For two nonnegative DEMs, A and B, the emission-measure-weighted discrepancy over region R is:

\[
D_R(A,B)=\frac{2\sum_{i\in R}\Delta T_i|A_i-B_i|}{\sum_{i\in R}\Delta T_i(A_i+B_i)}.
\]

The integration weights follow the inversion's logarithmic-grid convention:

\[
\Delta T_i=T_i\ln(10)\Delta\log_{10}T.
\]

The discrepancy lies between 0 (identical) and 2 (maximally different) when the denominator is positive. Smaller values indicate closer DEM agreement. The metric is symmetric and is not interpreted as an accuracy percentage.

An empty selected region or zero denominator yields an undefined (`NaN`) metric. The production analysis checks for undefined metrics instead of silently treating them as valid measurements.

### Signed integrated emission-measure bias

For recovered DEM A relative to reference DEM B:

\[
B_R(A\mid B)=\frac{EM_R(A)-EM_R(B)}{EM_R(B)},
\qquad
EM_R(A)=\sum_{i\in R}\Delta T_i A_i.
\]

The CSV reports bias as a percentage. A positive bias indicates excess integrated emission measure; a negative bias indicates deficient emission measure.

Synthetic truth is the reference for Python–truth and IDL–truth comparisons. IDL is the reference for Python–IDL comparisons.

### Raw χ²

The original raw χ² values from the Python and IDL standalone base results are retained for intensity-fit diagnostics.

These are **raw χ² values, not reduced χ² values**. A low χ² indicates that a solution fits the input intensities under its residual definition; it does not independently establish accurate recovery of the synthetic DEM.

## Results: 162 synthetic DEM cases

All 162 cases were processed, generating 972 individual comparison records (three pairs × two temperature regions). Processing completion establishes that the stored results could be compared successfully; it does **not** establish that all underlying inversions converged.

### DEM discrepancy distributions

| Region | Comparison | Median | 25th percentile | 75th percentile | 90th percentile |
|---|---|---:|---:|---:|---:|
| Gray band | Python vs truth | 0.6241 | 0.3031 | 1.0752 | 2.0000 |
| Gray band | IDL vs truth | 0.4765 | 0.2996 | 0.8154 | 1.9283 |
| Gray band | Python vs IDL | 0.5706 | 0.3703 | 1.9888 | 2.0000 |
| Full grid | Python vs truth | 0.8907 | 0.3433 | 1.5628 | 2.0000 |
| Full grid | IDL vs truth | 0.6837 | 0.4174 | 1.2381 | 1.9064 |
| Full grid | Python vs IDL | 0.7821 | 0.4719 | 1.9742 | 2.0000 |

In this dataset, IDL has a lower median truth-recovery discrepancy than Python in both regions. The median discrepancy is also lower inside the contribution band than across the full grid for all three comparisons.

These are descriptive sample-level findings; they do not establish statistical significance or general performance on real XRT observations.

### Paired truth-recovery comparisons

| Diagnostic | Gray band | Full grid |
|---|---:|---:|
| Python closer to truth | 74 | 80 |
| IDL closer to truth | 88 | 82 |
| Numerical ties (absolute tolerance 1e-10) | 0 | 0 |
| Median paired discrepancy difference (Python − IDL) | +0.0754 | +0.0109 |
| Python distances ≥ 1.999999 | 39 | 39 |
| IDL distances ≥ 1.999999 | 0 | 0 |

The positive median paired differences indicate larger Python discrepancies in the corresponding regions. The near-maximal-distance count is a **numerically defined diagnostic**, not a verified count of optimizer nonconvergence.

The 39 Python cases with near-maximal discrepancies constitute an important failure tail and should not be excluded from the primary analysis.

### Illustrative cases

**DEM 14:** Both implementations yield relatively small raw χ² values. The gray-band DEM discrepancies are 0.1957 for Python–truth, 0.3087 for IDL–truth, and 0.3887 for Python–IDL. This demonstrates that both inversions can fit the intensities while producing distinct DEM estimates.

**DEM 0:** Python's raw χ² is approximately 95.3832 versus 0.00317872 for IDL. Python's gray-band discrepancy relative to truth is approximately 2, compared with 0.2431 for IDL. The stored Python result has negligible recovered emission measure relative to truth. Earlier numerical diagnostics identified optimizer difficulties for this example, but this study does not independently establish their cause.

Extreme raw χ² values also occur in other Python and IDL cases. The logarithmic χ² visualization retains those cases rather than concealing them through linear-axis scaling.

## Reproducibility

From the `synthetic_dem_validation/` directory, with the existing input data and Python dependencies available:

```bash
python validate_synthetic_triplet.py
python plot_synthetic_triplet_validation.py
python -m pytest -q test_synthetic_triplet_metrics.py
```

The mathematical unit-test suite currently contains eight tests and has passed.

### Generated tables

Under `data/output_data/scientific_validation/`:

- `synthetic_triplet_metrics.csv` — one wide record per synthetic DEM, including all pairwise regional metrics and both raw χ² values.
- `synthetic_triplet_readable.csv` — six descriptive records per DEM, suitable for spreadsheet inspection and plotting.
- `synthetic_triplet_summary.csv` — median, quartiles, 90th percentile, and integrated EM-bias summaries for each comparison and region.

### Generated figures

Under `plots/scientific_validation/`:

- `gray_band_dem_discrepancy.png` — distribution of gray-band discrepancies.
- `gray_band_vs_full_grid.png` — paired comparison between the selected contribution band and the entire temperature grid.
- `python_vs_idl_chisq.png` — comparison of stored raw χ² on logarithmic axes.

## Scope and limitations

- This is a retrospective validation of **existing standalone base inversions**, not a new inversion experiment.
- The 162 synthetic configurations do not establish performance across all operational observing conditions.
- The contribution bands are supplied inputs and are not independently derived or physically validated in this study.
- The analysis does not test Monte Carlo uncertainty distributions or estimate confidence intervals on reconstruction accuracy.
- The stored χ² values are used as reported and are not independently recalculated from forward-modeled intensities.
- Solver convergence status is not comprehensively established by the current result files or summary tables; near-maximal discrepancies are therefore not automatically classified as optimizer failures.
- The study describes the specified 162 cases without imposing an arbitrary pass/fail accuracy threshold.
- Accurate intensity fitting does not uniquely determine the underlying DEM; comparisons with known synthetic truth are necessary to assess recovery.

## Interpretation for release

The analysis provides a reproducible, quantitative baseline for the existing XRTpy synthetic DEM results. It demonstrates substantial variation in truth-recovery accuracy and implementation agreement, including a pronounced Python discrepancy tail.

These findings support transparent characterization of current solver behavior. They do not establish equivalence to IDL, universal recovery accuracy, or readiness for every scientific application.

The generated tables and plots provide a baseline against which future solver changes can be evaluated without redefining the validation metrics.