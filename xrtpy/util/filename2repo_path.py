import logging
from pathlib import Path


def filename2repo_path(  # noqa: C901
    filename, urlroot="https://xrt.cfa.harvard.edu/", join=False, verbose=False
):
    """
    Given a filename with standard Level 1 or synoptics naming convention for
    files return its url in the CfA repository. Currently works for Level 1
    data, Level 1 quality data, Level 1 JPEG200 images, Level 2 synoptics
    composite fits files and Level 2 grade maps.
    Files in the archive that are not currently supported include: Level 0
    data, standard jpeg images (neither Level 1 nor synoptics), Level 2
    synoptics that are not composites (synop_XRT*.fits files)

    Parameters:
    -----------
    filename    string, required
                filename comforming to standard naming conventions

    urlroot     string
                root of the directory tree on the repository; by default the
                root on the CfA XRT repository

    join        boolean, optional
                if True, the joined path and filename are returned rather than
                just the url path

    verbose     boolean, optional
                if True, prints out inferred observation date and time for the data file

    Returns:
    --------
    path        string
                full url path for file (not including the file name)

    """

    if verbose:
        logging.basicConfig(format="%(funcName)s: %(message)s", level=logging.INFO)

    # deal with issue of Path converting // to /
    if "//" in urlroot:
        trans, root = urlroot.split("//")
        root = Path(root)
    else:
        trans = "https:"
        root = Path(urlroot)

    filename = Path(filename)
    if filename.suffix == ".jp2":
        # need special handling for JPEG2000 files
        parts = filename.split("_")
        yr = parts[0]
        mo = parts[1]
        dy = parts[2]
        hr = "H" + parts[4] + "00"
    else:
        # remove extentsion and then take everything after XRT as the date+time
        try:
            dt = filename.stem.split("XRT")[1]
        except ValueError:
            print("filename does not correspond to XRT data naming conventions")
        date, time = dt.split("_")
        yr = date[0:4]
        mo = date[4:6]
        dy = date[6:]
        hr = "H" + time[0:2] + "00"
    logging.info(f"yr, mo, dy, hr = {yr}, {mo}, {dy}, {hr}")
    if filename.name[:6] == "L1_XRT":
        if filename.name[-9:] == "qual.fits":
            path = root / "data_products" / "Level1_Qual" / yr / mo / dy / hr
        else:
            path = root / "level1" / yr / mo / dy / hr
    elif filename.name[:3] == "XRT":
        path = root / "level0" / yr / mo / dy / hr
        raise NotImplementedError(
            "Level 0 data is not currently available on" " the CfA website"
        )
    elif filename.name[:8] == "comp_XRT":
        # Note "synop_XRT", i.e. non-composite synoptics, are not currently
        # available on the web site
        if "gmap" in filename.name:
            path = root / "data_products" / "Level2_gmap" / yr / mo / dy / hr
        else:
            path = root / "data_products" / "Level2_Synoptics" / yr / mo / dy / hr
    elif filename.name.split("_")[-1] == "XRT.jp2":
        path = root / "data_products" / "jp2" / yr / mo / dy / hr
    elif (
        filename.name.split("_")[-1] == "COMP.jp2"
        or filename.name.split("_")[-1] == "SYNOP.jp2"
    ):
        path = root / "level2" / "synoptics" / yr / mo / dy / hr
        raise NotImplementedError(
            "Synoptics jpeg2000 data is not currently" " available on the CfA website"
        )
    else:
        path = None
        raise ValueError(
            "files with this naming convention are not available on the CfA website"
        )
    urlpath = trans + "//" + str(path) + "/"
    if join:
        urlpath = urlpath + str(filename)
    return urlpath
