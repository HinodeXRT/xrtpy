from pathlib import Path

import numpy as np
import pytest
import sunpy.map
from scipy.io import readsav

from xrtpy.response.temperature_from_filter_ratio import temperature_from_filter_ratio


def get_observed_data():
    directory = (
        Path(__file__).parent.absolute()
        / "data"
        / "temperature_from_filter_ratio_testing_files"
    )
    data_files = sorted(Path(directory).glob("L1_XRT20110128_*.*.fits"))

    return data_files


def get_IDL_results_data():
    directory = (
        Path(__file__).parent.absolute()
        / "data"
        / "temperature_from_filter_ratio_testing_files"
    )
    results_files = sorted(Path(directory).glob("IDL_results_*.sav"))

    return results_files


def rebin_image(data, binfac=1):
    """
    Given a data array and a binning factor return the data array rebinned by
    the binning factor.
    """

    s = data.shape
    ns = (s[0] // binfac, s[1] // binfac)
    rbs = (ns[0], binfac, ns[1], binfac)
    # sums the data in binfac x binfac sized regions
    drbin = data.reshape(rbs).mean(-1).mean(1)
    # for a boolean mask, this makes a pixel masked if any of the summed
    # pixels is masked. If we want to mask only if all the pixels are masked
    # then we could use prod in place of sum above

    if data.dtype == bool:
        # need to convert back to bool after summing
        drbin = drbin.astype(bool)
    return drbin


def test_standard_case():
    """
    Test case with all default values:
    no binning, no masking, no temperature range, standard thresholds
    """

    data_files = get_observed_data()
    # it turns out that the IDL test data was generated with the inverse order
    # of the files
    file1 = data_files[1]
    file2 = data_files[0]
    map1 = sunpy.map.Map(file1)
    map2 = sunpy.map.Map(file2)

    T_e, EM, Terr, EMerr = temperature_from_filter_ratio(map1, map2)
    T_EM = temperature_from_filter_ratio(map1, map2)
    T_e = T_EM.Tmap
    Terr = T_EM.Terrmap
    EM = T_EM.EMmap
    EMerr = T_EM.EMerrmap

    testdata = get_IDL_results_data()

    # This is needed because there are multiple test data sets, though so far
    # only a test written for the standard case
    fnames = [td.name for td in testdata]
    idata1 = fnames.index("IDL_results_bin1.sav")
    testdata1 = testdata[idata1]

    idldata = readsav(testdata1)
    goodT = (T_e.data > 0.0) & (idldata.te > 0.0)
    goodE = (EM.data > 0.0) & (idldata.em > 0.0)
    assert np.allclose(
        10.0 ** T_e.data[goodT], 10.0 ** idldata.te[goodT], atol=2.0e5, rtol=0.02
    )
    assert np.allclose(
        10.0 ** EM.data[goodE], 10.0 ** idldata.em[goodE], atol=4.0e44, rtol=0.03
    )
    assert np.allclose(
        10.0 ** Terr.data[goodT], 10.0 ** idldata.et[goodT], atol=1.0e4, rtol=0.08
    )
    assert np.allclose(
        10.0 ** EMerr.data[goodE], 10.0 ** idldata.ee[goodE], atol=4.0e43, rtol=0.02
    )


def test_binning_case():
    """
    Test case with following parameters:
    binning by a factor of 2
    no masking
    no temperature range
    standard thresholds
    """

    data_files = get_observed_data()
    # it turns out that the IDL test data was generated with the inverse order
    # of the files
    file1 = data_files[1]
    file2 = data_files[0]
    map1 = sunpy.map.Map(file1)
    map2 = sunpy.map.Map(file2)

    T_EM = temperature_from_filter_ratio(map1, map2, binfac=2)
    T_e = T_EM.Tmap
    Terr = T_EM.Terrmap
    EM = T_EM.EMmap
    EMerr = T_EM.EMerrmap

    testdata = get_IDL_results_data()

    # This is needed because there are multiple test data sets
    fnames = [td.name for td in testdata]
    idata1 = fnames.index("IDL_results_bin2.sav")
    testdata1 = testdata[idata1]

    idldata = readsav(testdata1)
    idlTe = rebin_image(idldata.te, 2)
    idlEM = rebin_image(idldata.em, 2)
    idlTerr = rebin_image(idldata.et, 2)
    idlEMerr = rebin_image(idldata.ee, 2)
    goodT = (T_e.data > 0.0) & (idlTe > 0.0)
    goodE = (EM.data > 0.0) & (idlEM > 0.0)

    assert np.allclose(
        10.0 ** T_e.data[goodT], 10.0 ** idlTe[goodT], atol=2.0e5, rtol=0.02
    )
    assert np.allclose(
        10.0 ** EM.data[goodE], 10.0 ** idlEM[goodE], atol=1.0e44, rtol=0.05
    )
    assert np.allclose(
        10.0 ** Terr.data[goodT], 10.0 ** idlTerr[goodT], atol=1.0e4, rtol=0.1
    )
    assert np.allclose(
        10.0 ** EMerr.data[goodE], 10.0 ** idlEMerr[goodE], atol=2.0e43, rtol=0.03
    )


def test_expmap_shape_mismatch_raises():
    """
    Test that an exposure map whose shape does not match its image raises a
    clear ValueError up front, rather than failing later with a confusing
    broadcasting error.
    """

    data_files = get_observed_data()
    file1 = data_files[1]
    file2 = data_files[0]
    map1 = sunpy.map.Map(file1)
    map2 = sunpy.map.Map(file2)

    good_expmap1 = np.full(map1.data.shape, map1.meta["EXPTIME"])
    good_expmap2 = np.full(map2.data.shape, map2.meta["EXPTIME"])
    bad_expmap = np.full((10, 10), 1.0)

    with pytest.raises(ValueError, match="expmap1 must match map1 shape"):
        temperature_from_filter_ratio(
            map1, map2, expmap1=bad_expmap, expmap2=good_expmap2
        )

    with pytest.raises(ValueError, match="expmap2 must match map2 shape"):
        temperature_from_filter_ratio(
            map1, map2, expmap1=good_expmap1, expmap2=bad_expmap
        )


def test_expmap_with_binfac_matches_scalar_exptime():
    """
    Test that a uniform exposure map combined with binning gives the same
    result as passing no exposure map at all (which uses the scalar EXPTIME).

    This guards the fix for combining expmap1/expmap2 with binfac > 1, and
    also catches a change from mean to sum when binning the exposure maps,
    since that would shift EM by log10(binfac**2).
    """

    data_files = get_observed_data()
    file1 = data_files[1]
    file2 = data_files[0]
    map1 = sunpy.map.Map(file1)
    map2 = sunpy.map.Map(file2)

    expmap1 = np.full(map1.data.shape, map1.meta["EXPTIME"])
    expmap2 = np.full(map2.data.shape, map2.meta["EXPTIME"])

    with_exp = temperature_from_filter_ratio(
        map1,
        map2,
        no_threshold=True,
        binfac=2,
        expmap1=expmap1,
        expmap2=expmap2,
    )
    without_exp = temperature_from_filter_ratio(
        map1,
        map2,
        no_threshold=True,
        binfac=2,
    )

    expected_shape = (map1.data.shape[0] // 2, map1.data.shape[1] // 2)
    assert with_exp.Tmap.data.shape == expected_shape

    assert np.allclose(with_exp.Tmap.data, without_exp.Tmap.data)
    assert np.allclose(with_exp.EMmap.data, without_exp.EMmap.data)
