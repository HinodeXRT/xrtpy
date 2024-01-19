import logging
import numpy as np
import requests

from astropy.io import fits
from pathlib import Path

from xrtpy.util.filename2repo_path import filename2repo_path


def make_exposure_map(comp_image_file, qualfiles=None, retsatpix=False, verbose=False):
    """
    Make an exposure map for a composite image. That is, using data from the
    header of a composite image, create an image where the pixel values are
    the exposure times used in the composite image

    Parameters:
    -----------
    comp_image_file : string or pathlib.Path object
        input composite image. Can be a double or a triple composite

    qualfiles : list of strings or Path objects, optional
        data quality filenames or paths that correspond to the images that
        make up the composite image. If the composite is a triple, should be
        list of two values, medium and long. If the composite is a double,
        only a single values should be given, the long exposure qualfile.
        Must be in order, [medium, long] if composite is a triple.

    retsatpix : boolean, optional
        if True, return the saturated pixel (boolean) images for the long
        image, and if the composite is a triple, also the medium image. If
        those images are returned then it is in addtiion to the exposure map,
        as a tuple: expmap, lng_sat or expmap, med_sat, lng_sat

    verbose : boolean, optional
        if True, print diagnostic information

    Returns:
    --------
    exp_map : numpy array
        image in the shape of the input composite image where each pixel
        has the exposure time that corresponds to that pixel in the composite
        image

    med_sat : boolean array, optional
        if retsatpix is True and the input image is a triple composite, then
        this array is returned with the pixels for which the medium exposure
        is saturated set to True and all others False

    lng_sat : boolean array, optional
        if retsatpix is True, then this array is returned with the pixels for
        which the long exposure is saturated set to True and all others False
    """

    if verbose:
        logging.basicConfig(format="%(funcName)s: %(message)s", level=logging.INFO)
    comp_header = fits.getheader(comp_image_file)
    short_exp_filepath = Path(comp_header["SRTFNAME"])
    logging.info(f"short exp. file: {short_exp_filepath}")
    short_exp_urlpath = filename2repo_path(short_exp_filepath, join=True)
    logging.info(f"short exp. url: {short_exp_urlpath}")
    response = requests.get(short_exp_urlpath)
    if response.status_code != requests.codes.ok:
        raise requests.exceptions.HTTPError(f"Error downloading {short_exp_urlpath}")
    with open(short_exp_filepath, "wb") as outfile:
        outfile.write(response.content)

    if "MEDFNAME" in comp_header:
        triple = True
        logging.info("Composite image is a triple")
        if qualfiles is None:
            medium_exp_filename = Path(comp_header["MEDFNAME"])
            medium_exp_qualpath = Path(
                medium_exp_filename.stem + ".qual" + medium_exp_filename.suffix
            )
            medium_exp_urlpath = filename2repo_path(medium_exp_qualpath, join=True)
            logging.info(f"medium exp. url: {medium_exp_urlpath}")
            response = requests.get(medium_exp_urlpath)
            if response.status_code != requests.codes.ok:
                raise requests.exceptions.HTTPError(
                    f"Error downloading {medium_exp_urlpath}"
                )
            with open(medium_exp_qualpath, "wb") as outfile:
                outfile.write(response.content)
        else:
            medium_exp_qualpath = Path(qualfiles[0])
    else:
        medium_exp_filename = None
        triple = False
        logging.info("Composite image is a double")

    if qualfiles is None:
        long_exp_filename = Path(comp_header["LNGFNAME"])
        long_exp_qualpath = Path(
            long_exp_filename.stem + ".qual" + long_exp_filename.suffix
        )
        long_exp_urlpath = filename2repo_path(long_exp_qualpath, join=True)
        response = requests.get(long_exp_urlpath)
        if response.status_code != requests.codes.ok:
            raise requests.exceptions.HTTPError(f"Error downloading {long_exp_urlpath}")
        with open(long_exp_qualpath, "wb") as outfile:
            outfile.write(response.content)
    else:
        if triple:
            long_exp_qualpath = Path(qualfiles[1])
        else:
            long_exp_qualpath = Path(qualfiles)

    naxis1 = comp_header["NAXIS1"]
    naxis2 = comp_header["NAXIS2"]

    srt_exp = fits.getval(short_exp_filepath, "EXPTIME")
    logging.info(f"Short image exposure time = {srt_exp}")
    lng_hdu = fits.open(long_exp_qualpath)
    lng_exp = lng_hdu[0].header["EXPTIME"]
    lng_gmap = lng_hdu[0].data
    lng_hdu.close()

    # do bit arithmetic to get saturated pixels
    lng_sat = (lng_gmap & 1).astype(bool)
    logging.info(f"No. of saturated pixels in long exposure = {np.sum(lng_sat)}")
    exp_map = np.ones((naxis1, naxis2)) * lng_exp
    exp_map[lng_sat] = srt_exp

    if triple:
        med_hdu = fits.open(medium_exp_qualpath)
        med_gmap = med_hdu[0].data
        med_exp = med_hdu[0].header["EXPTIME"]
        med_hdu.close()
        med_sat = (med_gmap & 1).astype(bool)
        exp_map[~med_sat & lng_sat] = med_exp
        logging.info(
            "No. of saturated pixels in medium exposure = " f"{np.sum(med_sat)}"
        )
        logging.info(f"Medium image exposure time = {med_exp}")
    logging.info(f"Long image exposure time = {lng_exp}")

    if retsatpix:
        if triple:
            return exp_map, med_sat, lng_sat
        else:
            return exp_map, lng_sat
    else:
        return exp_map
