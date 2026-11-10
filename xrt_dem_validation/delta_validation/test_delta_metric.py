import numpy as np
import pytest

from delta_metric import compute_delta


def test_every_bin_off_by_one_sigma():
    # 4 bins, each off by exactly 1 sigma -> sqrt(4) / 4 = 0.5
    true = np.array([10.0, 20.0, 30.0, 40.0])
    sigma = np.array([1.0, 2.0, 3.0, 4.0])
    calc = true + sigma
    mask = np.ones(4, dtype=bool)

    delta, n = compute_delta(true, calc, sigma, mask)

    assert n == 4
    assert delta == pytest.approx(0.5)


def test_perfect_match_is_zero():
    true = np.array([1e20, 5e21, 3e20])
    sigma = np.array([1e19, 1e20, 1e19])
    mask = np.ones(3, dtype=bool)

    delta, n = compute_delta(true, true, sigma, mask)

    assert n == 3
    assert delta == 0.0


def test_only_band_bins_are_used():
    # Bins outside the band are wildly wrong and have sigma=0;
    # they must be ignored completely.
    true = np.array([5.0, 10.0, 10.0, 5.0])
    calc = np.array([999.0, 12.0, 8.0, -999.0])
    sigma = np.array([0.0, 2.0, 2.0, 0.0])
    mask = np.array([False, True, True, False])

    delta, n = compute_delta(true, calc, sigma, mask)

    # Each band bin is off by 1 sigma -> sqrt(2) / 2
    assert n == 2
    assert delta == pytest.approx(np.sqrt(2) / 2)


def test_empty_band_returns_nan():
    x = np.ones(3)
    delta, n = compute_delta(x, x, x, np.zeros(3, dtype=bool))

    assert n == 0
    assert np.isnan(delta)


def test_zero_sigma_in_band_raises():
    x = np.ones(3)
    sigma = np.array([1.0, 0.0, 1.0])

    with pytest.raises(ValueError, match="sigma <= 0"):
        compute_delta(x, x, sigma, np.ones(3, dtype=bool))


def test_shape_mismatch_raises():
    with pytest.raises(ValueError, match="Shape mismatch"):
        compute_delta(np.ones(3), np.ones(4), np.ones(3), np.ones(3, dtype=bool))
