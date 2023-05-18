"""
Functionality for removing the visible light leak from XRT composite image data.
"""
__all__ = ["xrt_remove_lightleak"]

import numpy as np

from datetime import datetime
from pathlib import Path
from sunpy.map import Map


def xrt_remove_lightleak(in_map, kfact=1.0, leak_image=None, verbose=False):
    r"""
    Subtract light leak (visible stray light) image from XRT synoptic
    composite images.

    Parameters:
    -----------
    in_map : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the synoptic composite image which the visible stray light
        will be subtracted from

    kfact : float, default: 1.0
        k-factor to apply when subtracting the light leak image:
        ``out_data = in_data - k * [leak_img]``


    leak_image : float array, dimensions 1024x1024,  optional
        A leak image to subtract, if a non-default image is desired.
        It's assumed that the image is 1024x1024, prepped and exposure
        normalized.

    verbose : boolean, optional
        If True, print out extra messages.

    Returns:
    --------
    out_map : ~sunpy.map.sources.hinode.XRTMap
        |Map| of input image with the light leak removed. The metadata HISTORY
        is also updated to reflect the fact that the light leak was removed.

    Example:
    --------
    >>> file = ``"comp_XRT20200220_061539.6.fits"`` # doctest: +SKIP
    >>> in_map = Map(file) # doctest: +SKIP
    >>> out_map =xrt_remove_lightleak(in_map) # doctest: +SKIP

    Notes:
    ------
    (Taken from the IDL routine xrt_synleaksub.pro)
    1. XRT images obtained after 9-May-2021 suffer visible stray light
    contamination (light leak) due to the pre-filter failure, i.e.  tiny
    rupture development occurred multiple times as follows:

    phase 1 :  9-May-2012 12:00
    phase 2 : 14-Jun-2015 12:30
    phase 3 : 27-May-2017 11:00
    phase 4 : 29-May-2018 00:00
    phase 5 :  8-Jun-2022 12:40

    2. The light leak correction is done by simply subtracting the light leak
    image (visible stray light component included in each X-ray filter pair)
    obtained during the Hinode satellite's eclipse season, which occurs
    roughly May to August each year.  The light leak image has the following
    characteristics:

    * Leak pattern and intensity differ with each filter and also with the
      satellite pointing.
    * Leak pattern and intensity of the same filter and pointing are roughly
      constant during each stray light phase (above), but vary by ~10%
      depending on the growth of the contamination layer on the CCD (that
      repeats gradual increase and jump down after CCD bakeouts).

    3. The image for the light leak correction should therefore be selected
    for the same filter, same satellite pointing, same stray light phase then
    ideally be adjusted for the level of CCD contamination at the time of
    observation.  However, as a practical matter, preparing the leak image for
    every possible pointing is hard to achieve, while we have good amount of
    light leak measurements at the disk center pointing. The full-disk
    composite images are therefore corrected most reliably for the light leak.

    4. Intensity variation of the leak image due to the growth of CCD
    contamination layer is well determined for Ti_poly at the stray light
    phase 1 by using the intensity correlation between Ti_poly and Al_mesh
    images (cg. Takeda et al. 2016, SolPhys.  291, p.317). The resulting
    k-factor is obtained with the function GET_SLCORFACT_RAW.PRO, and the leak
    image subtraction has been already performed only for the Ti_poly SCIA
    images at the phase 1 (as of Feb-2022).

    """

    # Check to see if the lightleak has already been subtracted
    history = in_map.meta["HISTORY"]
    if "Light leak subtraction: DONE" in history:
        print(
            "HISTORY indicates light leak subtraction already done on image"
            ", returning input map"
        )
        return in_map
    # ********* select leak image from the archive *********
    if leak_image is None:
        dir_leak = Path(__file__).parent.absolute() / "data" / "leak_fits"
        filt1 = in_map.meta["EC_FW1"]
        filt2 = in_map.meta["EC_FW2"]
        fpair = str(filt1) + str(filt2)
        sfilt1 = in_map.meta["EC_FW1_"]
        sfilt2 = in_map.meta["EC_FW2_"]
        sfpair = f"{sfilt1}/{sfilt2}"
        sl_phase = check_sl_phase(in_map.meta["date_obs"])
        if sl_phase == 6:
            print(
                "Warning: light leak images for this period are not yet"
                " available. Defaulting to previous phase."
            )
            sl_phase = 5

    if fpair == "01":  # Open/Al_mesh
        sl_phase_dict = {
            2: "term_p2am_20150718_160913.fits",
            3: "term_p3am_20170808_180126.fits",
            4: "term_p4am_20180712_171919.fits",
            5: "term_p5am_20220709_180901.fits",
        }
    elif fpair == "10":  # Al_poly/Open
        sl_phase_dict = {
            2: "term_p2ap_20150620_172818.fits",
            3: "term_p3ap_20170809_183821.fits",
            4: "term_p4ap_20180712_171928.fits",
            5: "term_p5ap_20220709_180910.fits",
        }
    elif fpair == "20":  # C_poly/Open
        sl_phase_dict = {2: "term_p2cp_20150620_190645.fits"}
    elif fpair == "02":  # Open/Ti_poly
        sl_phase_dict = {
            1: "term_p1tp_20140515_182503.fits",
            2: "term_p2tp_20150718_160921.fits",
        }
    else:
        sl_phase_dict = {}
        print(f"No leak image available for this filter pair {sfpair}.")
        if verbose:
            info = (
                f"{in_map.meta['date_obs']}   {in_map.measurement}   "
                + f"{in_map.meta['naxis1']}x{in_map.meta['naxis2']}"
            )
            print(info)
        print("Output data are identical with the input")
        return in_map

    try:
        leak_fn = sl_phase_dict[sl_phase]
    except KeyError:
        print(f"No leak image available for phase {sl_phase}")
        if verbose:
            print(f"Filter pair used: {sfpair}")
        print("Output data are identical with the input")
        return in_map

    leak_map = Map(dir_leak / leak_fn)

    leak_image = leak_map.data
    if in_map.meta["naxis1"] == 2048:
        # assumes leak_image is 1024 x 1024 (I think)
        leak_image2 = rebin(leak_image, (2048, 2048)) * 0.25 * kfact
    else:
        leak_image2 = leak_image * kfact

    out_image = in_map.data - leak_image2

    hist_entry = f"Light leak subtraction: DONE with {leak_fn} and kfactor:{kfact}"
    out_meta = in_map.meta
    out_meta["history"] += f"\n{hist_entry}"

    if verbose:
        print(f"appended to the history : {hist_entry}")

    return Map((out_image, out_meta))


def check_sl_phase(date_obs):
    """
    Get the stray light phase number

    Parameters:
    -----------
    obs_date : string
        Observation date string formatted in ISO format: YYYY-MM-DDTHH:mm
        where YYYY is the year, MM is the month number (zero-padded), DD is
        the day of the month (zero-padded), HH is the hour and mm is the
        minute (note T separator between date and time; may also include
        seconds). This is the standard format for the value of date_obs in the
        FITS headers of XRT images.

    Returns:
    --------
    phase no. : integer
        The phase of the light leak as follows:

    phase 0 : 23-Oct-2006 10:00
    phase 1 :  9-May-2012 12:00
    phase 2 : 14-Jun-2015 12:30
    phase 3 : 27-May-2017 11:00
    phase 4 : 29-May-2018 00:00
    phase 5 :  8-Jun-2022 12:40
    phase 6 :  5-May-2023  4:39
    """

    time_p1 = datetime.strptime("09-May-2012 12:00", "%d-%b-%Y %H:%M")
    time_p2 = datetime.strptime("14-Jun-2015 12:30", "%d-%b-%Y %H:%M")
    time_p3 = datetime.strptime("27-May-2017 11:00", "%d-%b-%Y %H:%M")
    time_p4 = datetime.strptime("29-May-2018 00:00", "%d-%b-%Y %H:%M")
    time_p5 = datetime.strptime("08-Jun-2022 12:40", "%d-%b-%Y %H:%M")
    time_p6 = datetime.strptime("05-May-2023 04:39", "%d-%b-%Y %H:%M")

    intime = datetime.fromisoformat(date_obs)
    if intime >= time_p6:
        phase = 6
    elif intime >= time_p5:
        phase = 5
    elif intime >= time_p4:
        phase = 4
    elif intime >= time_p3:
        phase = 3
    elif intime >= time_p2:
        phase = 2
    elif intime >= time_p1:
        phase = 1
    else:
        phase = 0

    return phase


def my_rebin(image, newdims):
    """Rebin an array to a new shape"""

    newimage = np.zeros(newdims)
    origdims = image.shape
    ifac = newdims[0] // origdims[0]
    jfac = newdims[1] // origdims[1]
    for i in range(ifac):
        for j in range(jfac):
            newimage[i::ifac, j::jfac] = image[:]
    return newimage


def rebin(image, newshape):
    """
    Rebin an array to a new shape. Copied from:
    https://scipy-cookbook.readthedocs.io/items/Rebinning.html
    Should behave similarly to IDL rebin

    Parameters:
    -----------
    image : 2D numpy array
        image to be rebinned

    newshape : tuple
        dimensions of rebinned image. Each dimension should be an integer
        multiple of the corresponding dimension in the input image, though
        they need not be the same multiple

    Returns:
    --------
    newimage : 2D numpy array
        rebinned image
    """

    assert len(image.shape) == len(newshape)

    slices = [
        slice(0, old, float(old) / new) for old, new in zip(image.shape, newshape)
    ]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype("i")
    return image[tuple(indices)]
