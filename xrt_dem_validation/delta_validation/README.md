# delta_validation

Applies the DEM goodness-of-fit parameter Δ from Reeves & Weber (draft,
Eq. 3) to the 162 synthetic DEMs, for both XRTpy and IDL.

    Δ = (1 / n_temp) * sqrt( Σ (DEM_true − DEM_calc)² / σ² )

We use the paper's method, not its setup: our filter set, responses, and
MC count stay the same as the rest of the synthetic validation campaign.

## Choices made

- Linear DEM space, as the equation is written.
- `n_temp` is outside the square root, as written. Each bin off by 1σ gives
  Δ = 1/√n_temp, not 1.
- σ = plain median absolute deviation of the MC runs in each bin (no 1.4826 factor).
- The band is computed **per solver** from that solver's own DEM: bins
  contributing ≥2% of the modeled signal in at least one filter.

## Steps

1. [x] `delta_metric.py`: `compute_delta()` plus tests
2. [ ] σ per bin from the MC runs (`sigma_from_mc()`)
3. [ ] 2% band per solver (`contribution_band()`)
4. [ ] Run script: Δ for all 162 DEMs × {XRTpy, IDL} → CSV
5. [ ] Check Δ visually on a few DEMs
6. [ ] Figure 4-style plot: log χ²_obs vs log Δ, with R and slope

## Run tests

    cd delta_validation
    pytest test_delta_metric.py -v
    ruff check .
