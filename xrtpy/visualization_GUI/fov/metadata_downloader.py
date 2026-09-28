import urllib.request

import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.io import readsav


def fetch_metadata(xrt_downloaded_files, fast_bool=True):
    """
    Query Hinode/XRT data repositories to retrieve either level 0 (fast) or level 1 (slow) metadata

    Parameters:
    -----------
    xrt_downloaded_files : ~sunpy.net.fido_factory.UnifiedResponse
        fido query result for XRT data

    fast_bool : boolean, optional
        if True, do a fast retrieval of level 0 metadata

    Returns:
    --------
    hdul_lis : ~astropy.io.fits.PrimaryHDU
        HDU list containing the fits header for all XRT observations in your fido search
    """

    if fast_bool:
        print("Fast Metadata (Level 0)")
        return _download_metadata_fast(xrt_downloaded_files)
    else:
        print("Slow Metadata (Level 1)")
        hdul = _download_metadata_slow(xrt_downloaded_files)
        # hlis = []
        # for i in range(len(xrt_downloaded_files[0])):
        #    hlis.append(hdul[i + 2].header)
        hlis = [hdul[i + 2].header for i in range(len(xrt_downloaded_files[0]))]
        return hlis


def _download_metadata_slow(xrt_downloaded_files):
    """
    Query Hinode/XRT data repositories to retrieve level 1 metadata

    Parameters:
    -----------
    xrt_downloaded_files : ~sunpy.net.fido_factory.UnifiedResponse
        fido query result for XRT data

    Returns:
    --------
    hdul_lis : ~astropy.io.fits.PrimaryHDU
        HDU list containing the fits header for all XRT observations in your fido search
    """

    url_lis = xrt_downloaded_files[0][:]["fileid"]
    url_str_lis = list(url_lis)
    primary_hdu = fits.PrimaryHDU(data=np.ones((3, 3)))
    c1 = fits.Column(name="URL", array=url_str_lis, format="100A")
    c2 = fits.Column(
        name="header_int", array=np.asarray(range(len(url_str_lis))) + 2, format="J"
    )
    table_hdu = fits.BinTableHDU.from_columns([c1, c2])

    hdu_lis = fits.HDUList([primary_hdu, table_hdu])

    for i in range(len(url_str_lis)):
        if i % 10 == 0:
            print(int(1000.0 * i / len(url_str_lis)) / 10.0, "%")

        fsspec_kwargs = {"block_size": 100_000, "cache_type": "bytes"}
        with fits.open(
            url_lis[i], use_fsspec=True, fsspec_kwargs=fsspec_kwargs
        ) as hdul:
            # Download a single header
            t_header = hdul[0].header
            image_hdu = fits.ImageHDU(
                data=np.ones((100, 100)), header=t_header, name="header" + str(i)
            )
        hdu_lis.append(image_hdu)
    return hdu_lis


def _download_metadata_fast(xrt_downloaded_files):
    """
    Query Hinode/XRT data repositories to retrieve level 0 metadata

    Parameters:
    -----------
    xrt_downloaded_files : ~sunpy.net.fido_factory.UnifiedResponse
        fido query result for XRT data

    Returns:
    --------
    hdul_lis : ~astropy.io.fits.PrimaryHDU
        HDU list containing the fits header for all XRT observations in your fido search
    """

    html_str = _get_html_str()
    file_lis, cat_fi = _date_to_metafile(xrt_downloaded_files)
    genyl = _get_urls(file_lis, html_str)
    print("downloading")
    data_dict_lis = _get_metafile(genyl)
    print("done")
    hdu_lis = _match_vso_to_cat(data_dict_lis, cat_fi, xrt_downloaded_files)
    return hdu_lis


def _get_html_str():
    """
    Retrieve list of file names for .geny xrtcat files.

    Parameters:
    -----------

    Returns:
    --------
    html_lis : list
        list of strings for daily xrtcat .geny files
    """
    url_start = "https://sot.lmsal.com/data/sot/metadata/sswdb/hinode/xrt/xrt_genxcat/"
    with urllib.request.urlopen(url_start) as response:  # noqa : S310
        html = response.read()
        html_str = html.decode("utf-8")
    return html_str


def _date_to_metafile(xrt_download_list):
    """
    (Part of fast metadata retrieval)
    Get a partial list of daily xrtcat file names cat contain all observations requested

    Parameters:
    -----------
    xrt_downloaded_files : ~sunpy.net.fido_factory.UnifiedResponse
        fido query result for XRT data

    Returns:
    --------
    file_lis : list
        list containing the unique date string portion for each daily .geny file that contain the requested xrt observatio0ns
    cat_fi: integer list
        list containing indexes for each observation pointing at the correct daily .geny file
    """
    time_lis = xrt_download_list[0]["Start Time"]
    year_lis = xrt_download_list[0]["Start Time"].ymdhms.year
    month_lis = xrt_download_list[0]["Start Time"].ymdhms.month
    day_lis = xrt_download_list[0]["Start Time"].ymdhms.day
    cat_fi = []
    file_lis = []
    for i in range(len(time_lis)):
        year_str = str(year_lis[i])
        month_str = str(month_lis[i])
        day_str = str(day_lis[i])
        if len(day_str) < 2:
            day_str = "0" + day_str
        if len(month_str) < 2:
            month_str = "0" + month_str
        ndatei = year_str + month_str + day_str
        if ndatei in file_lis:
            cat_fi.append(file_lis.index(ndatei))
        else:
            file_lis.append(ndatei)
            cat_fi.append(file_lis.index(ndatei))

    return file_lis, cat_fi


def _get_urls(file_n_lis, html_str):
    """
    (Part of fast metadata retrieval)
    Retrieve the complete file name for each xrtcat file use the unique daily string for each file to search the complete file list

    Parameters:
    -----------
    file_n_lis : string list
        list with strings containing the date portion of the file name for each xrtcat file
    html_str :
        a long string containing the unformatted filenames for all xrtcat files

    Returns:
    --------
    geny_lis: string list
        list containing the complete xrtcat file name for all xrtcat files that contain the requested observations
    """

    nfile = len(file_n_lis)
    geny_lis = []
    for i in range(nfile):
        find_url = "xrt" + file_n_lis[i]
        findex = html_str.find(find_url)
        gen_fn = html_str[findex : findex + 35]
        findex2 = gen_fn.find("geny")
        gen_fn = gen_fn[: findex2 + 4]
        geny_lis.append(gen_fn)
    return geny_lis


def _get_metafile(geny_lis):
    """
    (Part of fast metadata retrieval)
    Retrieve metadata for all observations contained within each .geny files listed in geny_lis

    Parameters:
    -----------
    geny_lis : string list
        list containing the complete xrtcat file name for all xrtcat files that contain the requested observations

    Returns:
    --------
    data_dict_lis: dictionary list
        list containing a dictionary for each .geny file, each dictionary containing metadata for all observations listed in each .geny file

    """

    url_start = "https://sot.lmsal.com/data/sot/metadata/sswdb/hinode/xrt/xrt_genxcat/"
    ngeny = len(geny_lis)
    data_dict_lis = []
    for i in range(ngeny):
        gen_fn = geny_lis[i]
        f, h = urllib.request.urlretrieve(url_start + gen_fn)  # noqa : S310
        data2 = readsav(f)["p0"]
        data_dict2 = {k: data2[k] for k in data2.dtype.names}
        data_dict_lis.append(data_dict2)
    return data_dict_lis


def _match_vso_to_cat(data_dict_lis, cat_fi, xrt_download):
    """
    (Part of fast metadata retrieval)
    Match list of requested observations to the metadata retrieved from xrtcat

    Parameters:
    -----------
    data_dict_lis : list of dictionaries
        list containing a dictionary for each .geny file, each dictionary containing metadata for all observations listed in each .geny file
    cat_fi: integer list
        list containing indexes for each observation pointing at the correct daily .geny file (or at the correct dictionary in data_dict_lis)
    xrt_download : ~sunpy.net.fido_factory.UnifiedResponse
        fido query result for XRT data
    Returns:
    --------
    header_lis: list of dictionaries
        HDU list containing the fits header for all XRT observations in your fido search
    """
    n_dict = len(data_dict_lis)
    cat_time_lis = []
    for i in range(n_dict):
        data_dict = data_dict_lis[i]
        date_obs_cat = data_dict["DATE_OBS"]
        cat_str = [cat_bin.decode("ascii") for cat_bin in date_obs_cat]
        cat_time = Time(np.asarray(cat_str), format="isot", scale="utc")
        cat_time_lis.append(cat_time)

    min_ti_lis = []
    delt_lis = []
    delt_lisp = []
    delt_lism = []
    header_lis = []
    for i in range(len(xrt_download[0]["Start Time"])):
        cat_time = cat_time_lis[cat_fi[i]]
        stime = xrt_download[0]["Start Time"][i]
        dealt = cat_time - stime
        dealt = dealt.value * 24.0 * 3600.0  # convert from days to seconds
        min_ti = np.argmin(np.abs(dealt))
        min_ti_lis.append(min_ti)
        delt_lis.append(dealt[min_ti])
        try:
            delt_lisp.append(dealt[min_ti + 1])
            delt_lism.append(dealt[min_ti - 1])
        except IndexError:
            print()
        header_lis.append(_meta_to_dict(data_dict_lis[cat_fi[i]], min_ti))
    return header_lis


def _meta_to_dict(data_dict, di):
    """
    (Part of fast metadata retrieval)
    Convert ascii metadata for each obseration to a python dictionary

    Parameters:
    -----------
    data_dict : dictionary
        dictionary containing all of the metadata for all observations in a given .geny file
    di : integer
        index that matches the observation entry from fido search to the file entry in data_dict

    Returns:
    --------
    hdict : dictionary
        dictionary of all metadata for a single obseration corresponding to index di

    """
    dkeys = data_dict.keys()
    hdict = {}
    for dki in dkeys:
        try:
            hdict[dki] = data_dict[dki][di].decode("ascii")
        except IndexError:  # noqa : PERF203
            hdict[dki] = data_dict[dki][di]
    return hdict
