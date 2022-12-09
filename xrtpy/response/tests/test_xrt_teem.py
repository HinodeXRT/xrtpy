from scipy.io import readsav
import numpy as np
from astropy.io import fits
import pkg_resources
from pathlib import Path
import pytest

from xrtpy.response.xrt_teem import xrt_teem


def get_observed_data():

    directory = pkg_resources.resource_filename(
        "xrtpy", "response/tests/data/xrt_teem_testing_files"
    )
    data_files = sorted(Path(directory).glob('L1_XRT20110128_*.*.fits'))

    return data_files


def get_IDL_results_data():

    directory = pkg_resources.resource_filename(
        "xrtpy", "response/tests/data/xrt_teem_testing_files"
    )
    results_files = sorted(Path(directory).glob('IDL_results_*.sav'))
    
    return results_files


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
    hdu1 = fits.open(file1)
    data1 = hdu1[0].data
    hdr1 = hdu1[0].header
    hdu2 = fits.open(file2)
    hdr2 = hdu2[0].header
    data2 = hdu2[0].data
    hdu1.close()
    hdu2.close()

    Te,EM,Terr,EMerr = xrt_teem(hdr1, data1, hdr2, data2)

    testdata = get_IDL_results_data()

    # This is needed because there are multiple test data sets, though so far
    # only a test written for the standard case
    fnames = [td.name for td in testdata]
    idata1 = fnames.index('IDL_results_bin1.sav')
    testdata1 = fnames[idata1]

    idldata = readsav(testdata1)
    goodT = (Te > 0.) & (idldata.te > 0.)
    goodE = (EM > 0.) & (idldata.em > 0.)
    assert np.allclose(10.**Te[goodT], 10.**idldata.te[goodT], atol=1.E3,
            rtol=0.02)
    assert np.allclose(10.**EM[goodE], 10.**idldata.em[goodE], atol=8.E43,
            rtol=0.1)
    assert np.allclose(10.**Terr[goodT], 10.**idldata.et[goodT], atol=2.E4,
            rtol=0.05)
    assert np.allclose(10.**EMerr[goodT], 10.**idldata.ee[goodT], atol=2.E43,
            rtol=0.15)

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
    hdu1 = fits.open(file1)
    data1 = hdu1[0].data
    hdr1 = hdu1[0].header
    hdu2 = fits.open(file2)
    hdr2 = hdu2[0].header
    data2 = hdu2[0].data
    hdu1.close()
    hdu2.close()

    Te,EM,Terr,EMerr = xrt_teem(hdr1, data1, hdr2, data2, binfac=2)

    testdata = get_IDL_results_data()

    # This is needed because there are multiple test data sets, though so far
    # only a test written for the standard case
    fnames = [td.name for td in testdata]
    idata1 = fnames.index('IDL_results_bin2.sav')
    testdata1 = fnames[idata1]

    idldata = readsav(testdata1)
    goodT = (Te > 0.) & (idldata.te > 0.)
    goodE = (EM > 0.) & (idldata.em > 0.)
    assert np.allclose(10.**Te[goodT], 10.**idldata.te[goodT], atol=1.E3,
            rtol=0.02)
    assert np.allclose(10.**EM[goodE], 10.**idldata.em[goodE], atol=8.E43,
            rtol=0.1)
    assert np.allclose(10.**Terr[goodT], 10.**idldata.et[goodT], atol=2.E4,
            rtol=0.055)
    assert np.allclose(10.**EMerr[goodT], 10.**idldata.ee[goodT], atol=2.E43,
            rtol=0.15)
