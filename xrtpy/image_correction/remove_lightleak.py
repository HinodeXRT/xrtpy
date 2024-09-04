"""
Functionality for removing the visible light leak from XRT composite image data.
"""

__all__ = ["remove_lightleak"]

import warnings
from urllib.parse import urljoin

import astropy.time
import astropy.units as u
from sunpy.data import manager
from sunpy.map import Map
from sunpy.time import parse_time

from xrtpy.util import SSW_MIRRORS

LL_FILE_HASHES = {
    "term_p1cp_20140527_204601.fits": "bdb924a6ae62292980b266cec0fb96c3626d921efe8148a13b26f964513ea533",
    "term_p1tp_20140515_182503.fits": "0f5211633069cfb6a2d9d95f42d42bf52c4e16e1c47b89d32c567b430cd794b0",
    "term_p2am_20150718_160913.fits": "b2ba35e7b59d6785b5290900405b82814ebe3ca25c321d55d3853b452844df1f",
    "term_p2ap_20150620_172818.fits": "f08c96b86a5c2ed6e980d4a565a536476fe941b48877d2b22709b9694c2ce1c7",
    "term_p2cp_20150620_190645.fits": "0e7f72da2d408894980b8a14ead08e0fcd8706cb94c4aa4e256c60576c0c7678",
    "term_p2tp_20150718_160921.fits": "9fae8a831c52ecbe3d53d7e4d25caf78702a1aa5a117693d1bf75f0f200b30e6",
    "term_p3am_20170808_180126.fits": "ab8a92f728bfa7f68837c349923247a744ff5f86f2d996efc2c54b833455e3f3",
    "term_p3ap_20170809_183821.fits": "069d54e292d070322da0095efb903e79a1605371996a791e16e595420f557734",
    "term_p4am_20180712_171919.fits": "fd73856d8ea6eee83285ef91f0feb25f8568e6d0b88f77a2c35f4381c888abb3",
    "term_p4ap_20180712_171928.fits": "d57a92975e708999fabe22dbd4d491a9a10d9491ee4068a62c73ed4555f7e5dd",
    "term_p5am_20220709_180901.fits": "a052bc352cddfe2e7da69f90de83dd8549a457b8ef52acc18932ad41e255f410",
    "term_p5ap_20220709_180910.fits": "ab078b4d9573c4b737e7fb41311f824eb540c44c1fd8e3a588cd4e6e0ac2cb39",
}


def _get_stray_light_phase(date_obs):
    """
    Get the stray light phase number.

    The stray light phase is computed according to the date as follows:

    phase 0 : 23-Oct-2006 10:00
    phase 1 :  9-May-2012 12:00
    phase 2 : 14-Jun-2015 12:30
    phase 3 : 27-May-2017 11:00
    phase 4 : 29-May-2018 00:00
    phase 5 :  8-Jun-2022 12:40
    phase 6 :  5-May-2023  4:39

    Parameters:
    -----------
    obs_date : Any parse-time compatible format
        Observation date. This can be specified in any format understood by
        `sunpy.time.parse_time`.

    Returns:
    --------
    phase : `int`
        The phase of the light leak
    """
    phase_periods = astropy.time.Time(
        [
            "2012-05-09 12:00",
            "2015-06-14 12:30",
            "2017-05-27 11:00",
            "2018-05-29 00:00",
            "2022-06-08 12:40",
            "2023-05-05 04:39",
        ],
        format="iso",
        scale="utc",
    )

    date_obs = parse_time(date_obs)
    if date_obs >= phase_periods[5]:
        phase = 6
    elif date_obs >= phase_periods[4]:
        phase = 5
    elif date_obs >= phase_periods[3]:
        phase = 4
    elif date_obs >= phase_periods[2]:
        phase = 3
    elif date_obs >= phase_periods[1]:
        phase = 2
    elif date_obs >= phase_periods[0]:
        phase = 1
    else:
        phase = 0

    if phase == 6:
        warnings.warn(  # noqa: B028
            "light leak images for this period are not yet"
            " available. Defaulting to previous phase."
        )
        phase = 5

    return phase


def _select_lightleak_file(filter_wheel_1, filter_wheel_2, date):
    phase = _get_stray_light_phase(date)
    file_dict = {
        ("open", "al_mesh"): {
            2: "term_p2am_20150718_160913.fits",
            3: "term_p3am_20170808_180126.fits",
            4: "term_p4am_20180712_171919.fits",
            5: "term_p5am_20220709_180901.fits",
        },
        ("al_poly", "open"): {
            2: "term_p2ap_20150620_172818.fits",
            3: "term_p3ap_20170809_183821.fits",
            4: "term_p4ap_20180712_171928.fits",
            5: "term_p5ap_20220709_180910.fits",
        },
        ("c_poly", "open"): {
            2: "term_p2cp_20150620_190645.fits",
        },
        ("open", "ti_poly"): {
            1: "term_p1tp_20140515_182503.fits",
            2: "term_p2tp_20150718_160921.fits",
        },
    }
    fw_tuple = (filter_wheel_1.lower(), filter_wheel_2.lower())
    if fw_tuple not in file_dict:
        raise ValueError(
            f"No leak image available for {filter_wheel_1}/{filter_wheel_2}."
        )
    if phase not in file_dict[fw_tuple]:
        raise ValueError(
            f"No leak image available for {date} and {filter_wheel_1}/{filter_wheel_2}"
        )

    return file_dict[fw_tuple][phase]


def remove_lightleak(in_map, scale=1.0, leak_map=None):
    r"""
    Subtract visible stray light image from XRT synoptic composite images.

    Parameters:
    -----------
    in_map : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the synoptic composite image which the visible stray light
        will be subtracted from
    scale : float, default: 1.0
        Scaling factor to apply when subtracting the light leak image:
        ``out_data = in_data - scale * leak_img``
    leak_map : `~sunpy.map.Map`,  optional
        A leak image to subtract, if a non-default image is desired.
        It's assumed that the image is 1024x1024, prepped and exposure
        normalized.

    Returns:
    --------
    out_map : `~sunpy.map.sources.hinode.XRTMap`
        |Map| of input image with the light leak removed. The metadata HISTORY
        is also updated to reflect the fact that the light leak was removed.

    Example:
    --------
    >>> file = ``"comp_XRT20200220_061539.6.fits"`` # doctest: +SKIP
    >>> in_map = Map(file) # doctest: +SKIP
    >>> out_map = remove_lightleak(in_map) # doctest: +SKIP

    Notes:
    ------
    (Taken from the IDL routine :file:`xrt_synleaksub.pro`)
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
    images (cg. :cite:p:`Takeda:2016`, p.317). The resulting
    k-factor is obtained with the function :file:`GET_SLCORFACT_RAW.PRO`, and the leak
    image subtraction has been already performed only for the Ti_poly SCIA
    images at the phase 1 (as of Feb-2022).
    """
    if "Light leak subtraction: DONE" in in_map.meta["HISTORY"]:
        warnings.warn("HISTORY indicates light leak subtraction already done on image.")  # noqa: B028

    if leak_map is None:
        fw1 = in_map.meta['EC_FW1_']
        fw2 = in_map.meta['EC_FW2_']
        leak_filename = _select_lightleak_file(
            fw1, fw2, in_map.date
        )

        # NOTE: This function is being defined inline because the filename is only known once
        # the date and filter wheel combination are known.
        @manager.require(
            "ll_file",
            [
                urljoin(mirror, f"hinode/xrt/idl/util/leak_fits/{leak_filename}")
                for mirror in SSW_MIRRORS
            ],
            LL_FILE_HASHES[leak_filename],
        )
        def get_ll_file():
            return manager.get("ll_file")

        leak_map = Map(get_ll_file())
    else:
        # NOTE: This is to fill in the HISTORY in the resulting file
        leak_filename = "unknown"

    # case of 2048 x 2048 input image - for full resolution image, the light leak image flux
    # in each pixel must be split into four pixels each, thus the factor of 0.25 below
    if in_map.dimensions[0] == 2 * leak_map.dimensions[0]:
        leak_map = leak_map.resample(u.Quantity(in_map.dimensions))
        leak_map *= 0.25
    leak_map *= scale

    out_map = in_map - leak_map.quantity

    hist_entry = (
        f"Light leak subtraction: DONE with {leak_filename} and kfactor:{scale}"
    )
    out_map.meta["history"] += f"\n{hist_entry}"

    return out_map
