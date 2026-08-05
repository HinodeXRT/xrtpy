import sys

import pytest
from sunpy.net import Fido
from sunpy.net import attrs as a

from xrtpy.visualization_GUI.fov import metadata_downloader as xfetch

pytestmark = pytest.mark.skipif(sys.version_info >= (3, 13), reason="zeep/lxml is not installable on Python 3.13 due to missing prebuilt wheels"
)

time_range = a.Time("2011-06-07 06:00:00", "2011-06-07 07:30:54")
instrument = a.Instrument("xrt")

xrt_downloaded_files = Fido.search(time_range, instrument)
n_files = len(xrt_downloaded_files[0])

html_str = xfetch._get_html_str()


@pytest.mark.xfail(xfail_strict=True)
def test_download():
    assert n_files == 0


@pytest.mark.xfail(xfail_strict=True)
def test_html():
    assert len(html_str) == 0


@pytest.mark.parametrize(
    ("argument", "expected_result"),
    [
        (
            "20061017",
            "xrt20061017_0000_NE7.geny",
        ),  # pair the argument with the expected result
        ("20090405", "xrt20090405_0000_NE394.geny"),
    ],
)
def test_get_urls(argument, expected_result):
    geny_lis = xfetch._get_urls([argument], html_str)
    genystr = geny_lis[0]
    assert genystr == expected_result


file_lis, cat_fi = xfetch._date_to_metafile(xrt_downloaded_files)


def test_date_to_metafile_a():
    assert file_lis[0] == "20110607"


def test_date_to_metafile_b():
    assert len(file_lis) == 1


def test_date_to_metafile_c():
    assert len(cat_fi) == n_files


@pytest.mark.xfail(xfail_strict=True)
def test_date_to_metafile_d():
    assert len(cat_fi) == 0


@pytest.mark.parametrize(
    ("argument", "expected_result"),
    [
        (True, n_files),  # pair the argument with the expected result
        (False, n_files),
    ],
)
def test_fetch_metadata(argument, expected_result):
    hdul = xfetch.fetch_metadata(xrt_downloaded_files, fast_bool=argument)
    assert len(hdul) == expected_result
