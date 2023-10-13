import numpy as np
import sunpy.map

from pathlib import Path
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

    delta = 10.0 ** T_e.data[goodT] - 10.0 ** idlTe[goodT]  # noqa
    x = 10.0 ** idlTe[goodT]  # noqa

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
