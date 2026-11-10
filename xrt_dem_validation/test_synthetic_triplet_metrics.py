"""Unit tests for synthetic DEM validation metrics."""

import numpy as np
import pytest

from validate_synthetic_triplet import (
    temperature_bin_widths,
    dem_discrepancy,
    contribution_band_mask,
    integrated_em_bias,
)


@pytest.fixture
def grid_and_weights():
    logt = np.round(np.arange(5.5, 8.01, 0.1), 1)
    return logt, temperature_bin_widths(logt)


def test_temperature_integration_weights(grid_and_weights):
    logt, weights = grid_and_weights
    expected = 10.0**logt * np.log(10.0) * 0.1

    assert len(weights) == 26
    np.testing.assert_allclose(weights, expected, rtol=1e-12)


def test_identical_dems_have_zero_distance(grid_and_weights):
    _, weights = grid_and_weights
    a = np.arange(1.0, 27.0)

    assert dem_discrepancy(a, a, weights) == pytest.approx(0.0)


def test_symmetric_distance_and_amplitude_scaling(grid_and_weights):
    _, weights = grid_and_weights
    a = np.ones(26)
    b = 2.0 * a

    assert dem_discrepancy(a, b, weights) == pytest.approx(2.0 / 3.0)
    assert dem_discrepancy(b, a, weights) == pytest.approx(2.0 / 3.0)


def test_maximum_distance_and_zero_denominator(grid_and_weights):
    _, weights = grid_and_weights
    a = np.ones(26)
    zero = np.zeros(26)

    assert dem_discrepancy(a, zero, weights) == pytest.approx(2.0)
    assert np.isnan(dem_discrepancy(zero, zero, weights))


def test_disconnected_contribution_band(grid_and_weights):
    logt, _ = grid_and_weights

    # Includes two separated temperature intervals.
    mask = contribution_band_mask(logt, [5.5, 5.6, 6.0, 6.1])

    assert mask.sum() == 4
    np.testing.assert_allclose(logt[mask], [5.5, 5.6, 6.0, 6.1])


def test_gray_band_changes_selected_region(grid_and_weights):
    logt, weights = grid_and_weights
    band = contribution_band_mask(logt, [6.4, 6.5])

    a = np.ones(26)
    b = a.copy()
    b[~band] = 10.0

    assert dem_discrepancy(a, b, weights, band) == pytest.approx(0.0)
    assert dem_discrepancy(a, b, weights) > 0


def test_integrated_em_bias(grid_and_weights):
    logt, weights = grid_and_weights
    band = contribution_band_mask(logt, [6.4, 6.5])

    a = np.ones(26)

    assert integrated_em_bias(a, a, weights) == pytest.approx(0.0)
    assert integrated_em_bias(2 * a, a, weights) == pytest.approx(1.0)
    assert integrated_em_bias(a, 2 * a, weights) == pytest.approx(-0.5)
    assert integrated_em_bias(2 * a, a, weights, band) == pytest.approx(1.0)
    assert np.isnan(integrated_em_bias(a, np.zeros(26), weights))


def test_reject_invalid_inputs(grid_and_weights):
    logt, weights = grid_and_weights

    with pytest.raises(ValueError):
        temperature_bin_widths(logt[::-1])

    with pytest.raises(ValueError):
        contribution_band_mask(logt, [9.0])

    with pytest.raises(ValueError):
        dem_discrepancy(np.ones(25), np.ones(26), weights)

    with pytest.raises(ValueError):
        dem_discrepancy(-np.ones(26), np.ones(26), weights)

    with pytest.raises(ValueError):
        integrated_em_bias(np.ones(26), np.ones(25), weights)
