"""
Functionality for diagnosing plasma temperature through the filter
ratio technique.
"""

from datetime import datetime
from typing import Any, NamedTuple

import numpy as np
from astropy import units as u
from astropy.constants import c, h
from sunpy.coordinates.sun import B0, angular_radius
from sunpy.image.resample import reshape_image_to_4d_superpixel
from sunpy.map import Map

from xrtpy.response.temperature_response import TemperatureResponseFundamental

__all__ = ["temperature_from_filter_ratio"]


class TempEMdata(NamedTuple):
    Tmap: Any
    EMmap: Any
    Terrmap: Any
    EMerrmap: Any


def temperature_from_filter_ratio(
    map1,
    map2,
    abundance_model="coronal",
    binfac=1,
    Trange=None,
    no_threshold=False,
    Te_err_threshold=0.5,
    photon_noise_threshold=0.2,
    mask=None,
    expmap1=None,
    expmap2=None,
):
    r"""
    Get coronal temperatures and emission measures from a pair of images using
    the filter ratio method.

    .. note::

        The program uses solar spectra calculated with CHIANTI database ver.
        10.0 (density: :math:`10^9` cm\ :sup:`-3`\ , ionization
        equilibrium: ``chianti.ioneq``)

    Parameters:
    -----------
    map1 : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the first XRT level 1 data image.  If the image is
        normalized, then it is assumed that the un-normalized image can be
        recovered by multiplying by the exposure time (exposure time must be
        available in the metadata). It is also assumed that the metadata
        history will contain the string ``"XRT_RENORMALIZE"`` if the image has been
        normalized.

    map2 : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the second image (must use different filters from the first
        image). The image shape should match that in ``map1``. The same
        considerations apply as for ``map1``.

    abundance_model : string, optional
        Chianti abundance model for the spectrum assumed when deriving the
        plasma temperatures and emission measures. The default is coronal
        abundances. Other currently available choices are ``"photospheric"`` and
        ``"hybrid"``.

    binfac : integer, optional (default = 1)
        spatial binning factor

    Trange : 2 element sequence containing floats [Optional]
        Range of log10(temperature) values to examine. Must be in order from
        lower to higher.

    Te_err_threshold : float [Optional]
        Threshold value for the temperature error. (default = 0.5 [50%])
        Values of temperature with error exceeding the threshold are set to 0.

    photon_noise_threshold : float [Optional]
        Threshold value for ratio of photon noise to signal. (default = 0.2
        [20%]). If photon noise/signal exceeds this value, the corresponding
        temperatures are set to 0.

    mask : Boolean array of the same shape as the images [Optional]
        If provided, masks out parts of the images from the analysis. Note:
        pixels to be masked out should be `True`, unmasked should be `False`.

    no_threshold : Boolean [Optional]
        If True, no thresholds are set. (default = False)

    expmap1 : numpy array [Optional]
        if provided, gives exposure time (s) for each pixel in image 1. This
        is useful for composite images in which different parts of the image
        have different exposure times

    expmap2 : numpy array [Optional]
        if provided, gives exposure time (s) for each pixel in image 2.

    Returns:
    --------
    TempEMdata : namedtuple of ~sunpy.map.sources.hinode.XRTMap objects
        namedtuple containing the attributes Tmap, EMmap, Terrmap and EMerrmap
        where the maps correspond to images and metadata with:
        log10 of the derived electron temperature [K], log10 of the derived
        volume emission measure [cm^-3], uncertainty in log10 temperature [K],
        and uncertainty in log10 volume emission measure [cm^-3].

    Examples:
    ---------
    Using this function, you can derive the coronal temperature using
    filter ratio method.

    >>> T_EM = temperature_from_filter_ratio(map1, map2) # doctest: +SKIP

    If you want to bin the image data in space to reduce photon noise, set
    binfac to the factor by which you want to bin.  For example to bin the
    data by a factor of 3 do:

    >>> T_EM = temperature_from_filter_ratio(map1, map2, binfac=3) # doctest: +SKIP

    The data is binned first and then the temperature is derived. Note that
    the image size is reduced by the factor binfac in each dimension, which
    contrasts with the behavior of the IDL routine xrt_teem.pro.

    Notes
    -----
    The returned values of pixels where the temperature cannot be derived
    or the error is greater than the threshold or the photon noise is
    greater than the threshold are set to 0. The EM for those pixels is
    also set to 0. The masks included in the output maps mask out pixels for
    which the data for either image is 0 or for which temperature cannot be
    determined, as well as for those masked out by an input mask (if
    provided).

    The details of the coronal-temperature-diagnostic capability of
    Hinode/XRT is described in
    Narukage et al. 2011, Solar Phys., 269, 169.
    http://adsabs.harvard.edu/doi/10.1007/s11207-010-9685-2
    and
    Narukage et al. 2013, Solar Phys.,
    http://adsabs.harvard.edu/doi/10.1007/s11207-013-0368-7
    These two papers are the reference papers of this program.

    Modification History:

    IDL routine written by N.Narukage (NAOJ). See original IDL code for more details.
    """

    hdr1 = map1.meta
    hdr2 = map2.meta
    data1 = map1.data
    data2 = map2.data
    if data1.shape != data2.shape:
        raise ValueError("The input images must be the same size")
    data1 = data1.astype(float)
    data2 = data2.astype(float)

    n1 = "XRT_RENORMALIZE" in hdr1["HISTORY"]
    n2 = "XRT_RENORMALIZE" in hdr2["HISTORY"]
    # This allows use of normalized data (contrary to original IDL code):
    if n1:
        if expmap1 is None:
            data1 = data1 * hdr1["EXPTIME"]
        else:
            data1 = data1 * expmap1
    if n2:
        if expmap2 is None:
            data2 = data2 * hdr2["EXPTIME"]
        else:
            data2 = data2 * expmap2

    if mask is None:
        mask = np.zeros_like(data1, dtype=bool)
    # if the input data have already been masked then this preserves those masks
    if map1.mask is not None:
        mask1 = mask | map1.mask
    else:
        mask1 = mask
    if map2.mask is not None:
        mask2 = mask | map2.mask
    else:
        mask2 = mask
    mask = mask1 | mask2
    map1 = Map(data1, map1.meta, mask=mask)
    map2 = Map(data2, map2.meta, mask=mask)

    if binfac > 1:
        map1 = map1.superpixel([binfac, binfac] * u.pix)
        map2 = map2.superpixel([binfac, binfac] * u.pix)
        data1 = map1.data
        data2 = map2.data
        # fix for issue with binning in superpixel - mask not summed
        mask = (
            reshape_image_to_4d_superpixel(mask, [binfac, binfac], [0, 0]).sum(3).sum(1)
        )
        mask = mask.astype(bool)
        map1.mask = mask
        map2.mask = mask
        # binning updates header keywords
        hdr1 = map1.meta
        hdr2 = map2.meta
    # input mask for data should be False in parts of the images to be used and
    # True in places we want to mask out
    # Here we additionally mask out pixels in which the data in either image
    # is <= 0
    dmask = (data1 <= 0.0) | (data2 <= 0.0)
    mask = mask | dmask
    map1 = Map(data1, map1.meta, mask=mask)
    map2 = Map(data2, map2.meta, mask=mask)

    filt1 = measurement_to_filtername(map1.measurement)
    date_obs1 = hdr1["DATE_OBS"]
    tresp1 = TemperatureResponseFundamental(filt1, date_obs1, abundance_model)

    filt2 = measurement_to_filtername(map2.measurement)
    date_obs2 = hdr2["DATE_OBS"]
    tresp2 = TemperatureResponseFundamental(filt2, date_obs2, abundance_model)

    if filt1 == filt2:
        raise ValueError("Filters for the two images cannot be the same")

    T_e, EM, model_ratio, ok_pixel = _derive_temperature(
        map1,
        map2,
        tresp1,
        tresp2,
        expmap1=expmap1,
        expmap2=expmap2,
        binfac=binfac,
        Trange=Trange,
    )

    T_error, EMerror, K1, K2 = calculate_TE_errors(
        map1, map2, T_e, EM, model_ratio, tresp1, tresp2, Trange=Trange
    )

    T_e = T_e.filled(0.0)
    EM = EM.filled(0.0)
    T_error = T_error.filled(0.0)
    EMerror = EMerror.filled(0.0)
    if not no_threshold:
        ok_wothr = ok_pixel.copy()
        Kd1 = np.sqrt(K1 / data1)
        Kd1 = Kd1.filled(0.0)
        Kd2 = np.sqrt(K2 / data2)
        Kd2 = Kd2.filled(0.0)
        tthr = (T_error - T_e) <= np.log10(Te_err_threshold)
        k1thr = Kd1 <= photon_noise_threshold
        k2thr = Kd2 <= photon_noise_threshold
        ok_pixel = (
            ok_pixel
            & ((T_error - T_e) <= np.log10(Te_err_threshold))
            & (Kd1 <= photon_noise_threshold)
            & (Kd2 <= photon_noise_threshold)
        )
        npix_thresh = np.sum(~ok_pixel)
        npix_Terr = np.sum(~tthr)
        npix_d1_noise = np.sum(~k1thr)
        npix_d2_noise = np.sum(~k2thr)
        npix_badpix = np.sum(~ok_wothr)
        extra_metadata = {
            "NPTHRESH": (npix_thresh, "number of pixels ruled out by threshold"),
            "NPTERR": (npix_Terr, "number of pixels ruled out by T_e errors"),
            "NPD1NOIS": (npix_d1_noise, "number of pixels ruled out by d1 noise"),
            "NPD2NOIS": (npix_d2_noise, "number of pixels ruled out by d2 noise"),
            "NPBPIX": (npix_badpix, "number of bad pixels before threshold"),
        }
        mask = mask | ~ok_pixel
        T_e[mask] = 0.0
        EM[mask] = 0.0
        T_error[mask] = 0.0
        EMerror[mask] = 0.0

        Tmodel = tresp1.CHIANTI_temperature.value
        if Trange:
            Tmodel = Tmodel[
                (np.log10(Tmodel) >= Trange[0]) & (np.log10(Tmodel) <= Trange[1])
            ]
        extra_metadata["TEMIN"] = (Tmodel.min(), "lowest T_e examined")
        extra_metadata["TEMAX"] = (Tmodel.max(), "highest T_e examined")
        extra_metadata["TETHRESH"] = (Te_err_threshold, "T_e error threshold")
        extra_metadata["PHNOISE"] = (photon_noise_threshold, "Photon noise threshold")
    else:
        Tmodel = tresp1.CHIANTI_temperature.value
        extra_metadata = {
            "TEMIN": (Tmodel.min(), "lowest T_e examined"),
            "TEMAX": (Tmodel.max(), "highest T_e examined"),
        }
    Tmap, EMmap, Terrmap, EMerrmap = make_results_maps(
        hdr1, hdr2, T_e, EM, T_error, EMerror, extra_metadata
    )
    return TempEMdata(Tmap, EMmap, Terrmap, EMerrmap)


def deriv(x, y):
    """
    Use three-point Lagrangian interpolation to compute the first derivative
    of an array of data

    Method taken from IDL routine of the same name.

    Parameters:
    -----------
    x : float array
        abscissa of the data
    y : float array
        ordinate of the data

    Returns:
    --------
    float array :
        first derivative at the data points
    """
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    y0 = y[:-2]
    y1 = y[1:-1]
    y2 = y[2:]
    x01 = x0 - x1
    x02 = x0 - x2
    x12 = x1 - x2
    dydx1 = (
        y0 * x12 / (x01 * x02) + y1 * (1.0 / x12 - 1.0 / x01) - y2 * x01 / (x02 * x12)
    )
    dydx0 = (
        y0[0] * (x01[0] + x02[0]) / (x01[0] * x02[0])
        - y1[0] * x02[0] / (x01[0] * x12[0])
        + y2[0] * x01[0] / (x02[0] * x12[0])
    )
    dydxN = (
        -y0[-1] * x12[-1] / (x01[-1] * x02[-1])
        + y1[-1] * x02[-1] / (x01[-1] * x12[-1])
        - y2[-1] * (x02[-1] + x12[-1]) / (x02[-1] * x12[-1])
    )
    return np.append(np.insert(dydx1, 0, dydx0), dydxN)


def _derive_temperature(
    map1, map2, tresp1, tresp2, expmap1=None, expmap2=None, binfac=1, Trange=None
):
    """
    Given two XRT Level 1 images, their associated metadata and the
    TemperatureResponseFundamental objects associated with them, derive the
    temperatures and emission measures in the image using the filter ratio
    method.

    Parameters:
    -----------
    map1 : ~sunpy.map.sources.hinode.XRTMap
        map for the first XRT level 1 data image

    map2 : ~sunpy.map.sources.hinode.XRTMap
        map for the second XRT level 1 data image

    tresp1 : ~xrtpy.response.temperature_response.TemperatureResponseFundamental
        temperature response for first image

    tresp2 : ~xrtpy.response.temperature_response.TemperatureResponseFundamental
        temperature response for second image

    expmap1 : numpy array, Optional
        exposure time map for data image 1

    expmap2 : numpy array, Optional
        exposure time map for data image 2

    binfac : integer, Optional (default = 1)
        spatial binning factor

    Trange : 2 element sequence containing floats, Optional
        Range of log10(temperature) values to examine. Must be in order from
        lower to higher. (Passed from temperature_from_filter_ratio.)


    Returns:
    --------
    T_e : 2D float array
        derived temperatures for the images

    EM : 2D float array
        derived emission measures for the images

    model_ratio : 1D float array
        emission ratios for the filters for the images as a function of
        temperature

    ok_pixel : 2D boolean array
        indicates which pixels are okay (True) and which aren't (False) based
        on whether the derived T_e corresponds to one (and not more than one)
        model temperature based on the flux ratio
    """

    Tmodel = tresp1.CHIANTI_temperature.value
    logTmodel = np.log10(Tmodel)

    flux1 = tresp1.temperature_response().value
    flux2 = tresp2.temperature_response().value
    ratio = flux1 / flux2
    t2mk = np.argmin(np.abs(Tmodel - 2.0e6))
    rev_ratio = ratio[t2mk] > 1.0

    # Need to convert column emission measure to volume emission measure - to
    # do that we multiply the EM by the area (in cm^2) of a sky pixel
    plate_scale = map1.meta["PLATESCL"]
    lsun = 6.95700e10 / angular_radius(map1.meta["DATE_OBS"]).value
    em2vem = (plate_scale * lsun) ** 2
    flux1 /= em2vem
    flux2 /= em2vem

    # Note: this requires Trange to be a 2 element iterable ordered low to
    # high (Trange is log10 of temperature limits).
    if Trange:
        in_trange = (logTmodel >= Trange[0]) & (logTmodel <= Trange[1])
        if np.any(in_trange):
            Tmodel = Tmodel[in_trange]
            logTmodel = np.log10(Tmodel)
            flux1 = flux1[in_trange]
            flux2 = flux2[in_trange]
        else:
            raise ValueError(
                "The temperature response does not include"
                " any of the input temperatures in Trange"
            )

    data1 = np.ma.masked_where(map1.mask, map1.data)
    data2 = np.ma.masked_where(map2.mask, map2.data)
    mask = map1.mask

    if expmap1 is None:
        exptime1 = map1.meta["EXPTIME"]
    else:
        exptime1 = expmap1
    if expmap2 is None:
        exptime2 = map2.meta["EXPTIME"]
    else:
        exptime2 = expmap2

    if rev_ratio:
        data_ratio = (data2 / exptime2) / (data1 / exptime1)
        model_ratio = flux2 / flux1
    else:
        data_ratio = (data1 / exptime1) / (data2 / exptime2)
        model_ratio = flux1 / flux2
    ok_num = np.zeros(data_ratio.shape, dtype=int)
    ok_cnt = np.zeros(data_ratio.shape, dtype=int)
    for i, m in enumerate(model_ratio[1:]):
        # n is True for all pixels that have flux ratios that are between the
        # model ratio for i and i + 1
        n = (data_ratio >= min(model_ratio[i], m)) & (
            data_ratio < max(model_ratio[i], m)
        )
        if np.sum(n) > 0:
            ok_num[n] = i + 1
            # ok_cnt keeps track of the number of model ratios that match the
            # data ratio for each pixel - if it ends up being more than 1 then
            # the model ratio is non-monotonic with temperature and those
            # pixels are assigned a temperature of 0
            ok_cnt[n] += 1

    # Performs linear interpolation between two adjacent temperatures in Tmodel,
    # based on the observed filter ratio (`data_ratio`) and the model ratios (`model_ratio`).
    # `a` and `b` represent differences between the observed and model ratios.
    # The equation weighs the temperatures proportionally, with more weight given
    # to the temperature whose model ratio is closer to the observed filter ratio.
    # This ensures a smooth transition between temperature values, providing a more
    # accurate estimate when the observed filter ratio lies between two model ratios.
    ok_num[ok_cnt != 1] = 0
    a = np.abs(model_ratio[ok_num] - data_ratio)
    b = np.abs(model_ratio[np.maximum((ok_num - 1), 0)] - data_ratio)
    T_e = (
        np.log10(Tmodel[ok_num]) * b + np.log10(Tmodel[np.maximum((ok_num - 1), 0)]) * a
    ) / (a + b)
    T_e = np.ma.masked_where(T_e <= 0.0, T_e)

    # This masks out pixels where more than one temperature is consistent with
    # the observed flux ratio - so temperatures for which the model filter
    # ratio is non-monotonic are effectively removed
    ok_pixel = (ok_cnt == 1) & (~mask)
    T_e[~ok_pixel] = 0.0

    # Note that T_e is the log10 of the electron temperature here
    DN = np.interp(T_e, logTmodel, flux1, left=0.0, right=0.0)
    DN = np.ma.masked_where(((DN <= 0.0) | (T_e <= 0.0)), DN)
    EM = np.ma.log10(data1 / (DN * exptime1)) - np.log10(binfac**2)
    EM[~ok_pixel] = 0.0
    return T_e, EM, model_ratio, ok_pixel


def calculate_TE_errors(map1, map2, T_e, EM, model_ratio, tresp1, tresp2, Trange=None):
    """
    Given values for T_e and EM derived via the filter ratio method and the
    values of the temperature and model flux ratio, return the errors (i.e.
    uncertainties) in T_e and EM.

    Parameters:
    -----------
    map1 : ~sunpy.map.sources.hinode.XRTMap
        map for the first XRT level 1 data image

    map2 : ~sunpy.map.sources.hinode.XRTMap
        map for the second XRT level 1 data image

    T_e : 2D float array
        Previously derived temperatures for an image pair

    EM : 2D float array
        Previously derived volume emission measures for an image pair

    model_ratio : 1D float array
        emission ratios for the filters for the images as a function of
        temperature

    tresp1: ~xrtpy.response.temperature_response.TemperatureResponseFundamental
        container for model data for the filters for image 1

    tresp2: ~xrtpy.response.temperature_response.TemperatureResponseFundamental
        container for model data for the filters for image 2

    Trange : 2 element sequence containing floats [Optional]
        Range of log10(temperature) values to examine. Must be in order from
        lower to higher. (Passed from temperature_from_filter_ratio.)

    Returns:
    --------
    T_error : 2D float array
        Estimate of the statistical errors in the temperature determination

    EMerror : 2D float array
        Estimate of the statistical errors in the volume emission measure
        determination

    K1 : 2D float array
        Narukage's K factor (K(2)) for image 1

    K2 : 2D float array
        Narukage's K factor for image 2
    """

    wvl = tresp1.wavelength
    eVe = tresp1.ev_per_electron
    gain = tresp1.ccd_gain_right
    # (h*c/lambda) * 1/(eV per electron) * 1/gain
    e2dn = (h.to(u.eV * u.s) * c.to(u.angstrom / u.s) / (wvl * eVe * gain)).value
    dwvl = wvl[1:] - wvl[:-1]
    dwvl = np.append(dwvl, dwvl[-1]).value

    Tmodel = tresp1.CHIANTI_temperature.value
    logTmodel = np.log10(Tmodel)
    effarea1 = tresp1.effective_area().value
    effarea2 = tresp2.effective_area().value
    spect1 = tresp1.spectra().value
    spect2 = tresp2.spectra().value
    flux1 = tresp1.temperature_response().value
    flux2 = tresp2.temperature_response().value
    if Trange:
        in_trange = (logTmodel >= Trange[0]) & (logTmodel <= Trange[1])
        if np.any(in_trange):
            Tmodel = Tmodel[in_trange]
            logTmodel = np.log10(Tmodel)
            spect1 = spect1[in_trange, :]
            spect2 = spect2[in_trange, :]
            flux1 = flux1[in_trange]
            flux2 = flux2[in_trange]
        else:
            raise ValueError(
                "The temperature response does not include"
                " any of the input temperatures in Trange"
            )

    dlnR_dlnT_mod = np.abs(deriv(np.log(Tmodel), np.log(model_ratio)))
    dlnR_dlnT = np.interp(T_e, logTmodel, dlnR_dlnT_mod, left=0.0, right=0.0)
    dlnR_dlnT = np.ma.masked_where(((dlnR_dlnT == 0.0) | (T_e <= 0.0)), dlnR_dlnT)
    # These sums should be converted to integrals if the spectrum includes
    # line broadening (as is the default in ChiantiPy).
    K1_mod = np.array(
        [
            (s1 * effarea1 * e2dn**2 * dwvl).sum() / (s1 * effarea1 * e2dn * dwvl).sum()
            for s1 in spect1
        ]
    )
    K2_mod = np.array(
        [
            (s2 * effarea2 * e2dn**2 * dwvl).sum() / (s2 * effarea2 * e2dn * dwvl).sum()
            for s2 in spect2
        ]
    )
    K1 = np.interp(T_e, logTmodel, K1_mod, left=0.0, right=0.0)
    K2 = np.interp(T_e, logTmodel, K2_mod, left=0.0, right=0.0)
    K1 = np.ma.masked_where(((K1 == 0.0) | (T_e <= 0.0)), K1)
    K2 = np.ma.masked_where(((K2 == 0.0) | (T_e <= 0.0)), K2)

    dlnf1_dlnT_mod = deriv(np.log(Tmodel), np.log(flux1))
    dlnf1_dlnT = np.interp(T_e, logTmodel, dlnf1_dlnT_mod, left=0.0, right=0.0)
    dlnf1_dlnT = np.ma.masked_where(((dlnf1_dlnT == 0.0) | (T_e <= 0.0)), dlnf1_dlnT)

    dlnf2_dlnT_mod = deriv(np.log(Tmodel), np.log(flux2))
    dlnf2_dlnT = np.interp(T_e, logTmodel, dlnf2_dlnT_mod, left=0.0, right=0.0)
    dlnf2_dlnT = np.ma.masked_where(((dlnf2_dlnT == 0.0) | (T_e <= 0.0)), dlnf2_dlnT)

    data1 = np.ma.masked_where(map1.mask, map1.data)
    data2 = np.ma.masked_where(map2.mask, map2.data)
    T_error = np.ma.log10(np.sqrt(K1 / data1 + K2 / data2) / dlnR_dlnT) + T_e

    EMerror = (
        np.ma.log10(
            np.ma.sqrt(dlnf2_dlnT**2 * K1 / data1 + dlnf1_dlnT**2 * K2 / data2)
            / dlnR_dlnT
        )
        + EM
    )
    return T_error, EMerror, K1, K2


def make_results_maps(hdr1, hdr2, T_e, EM, T_error, EMerror, extra_metadata):
    """
    Create SunPy Map objects from the image metadata and temperature, volume
    emission measure, temperature uncertainty and emission measure uncertainty
    data derived for a pair of XRT Level 1 images.

    Parameters:
    -----------
    hdr1 : metadata dictionary
        metadata associated with image 1

    hdr2 : metadata dictionary
        metadata associated with image 2

    T_e : 2D float array
        image containing the temperatures (log10(T(K))) derived for images 1 & 2

    EM : 2D float array
        image containing the volume emission measures (log10(e.m.(cm^-3)))
        derived for images 1 & 2

    T_error : 2D float array
        image containing the uncertainties in T_e derived for the images

    EMerror : 2D float array
        image containing the uncertainties in EM derived for the images

    Returns:
    --------
    Tmap : ~sunpy.map.sources.hinode.XRTMap
        Map containing T_e and metadata

    EMmap : ~sunpy.map.sources.hinode.XRTMap
        Map containing EM and metadata

    Terrmap : ~sunpy.map.sources.hinode.XRTMap
        Map containing T_error and metadata

    EMerrmap : ~sunpy.map.sources.hinode.XRTMap
        Map containing EMerror and metadata
    """

    date_obs1 = hdr1["date_obs"]
    datename = "".join(date_obs1.split("T")[0].split("-"))
    timename = "".join(date_obs1.split("T")[1].split(":"))[:8]
    # This is the data file name associated with the give date_obs
    filename1 = f"L1_XRT{datename}_{timename}.fits"
    date_obs2 = hdr2["date_obs"]
    datename = "".join(date_obs2.split("T")[0].split("-"))
    timename = "".join(date_obs2.split("T")[1].split(":"))[:8]
    filename2 = f"L1_XRT{datename}_{timename}.fits"

    rsun_ref = 6.95700e08
    rsun_obs = angular_radius(hdr1["DATE_OBS"]).value
    dsun = rsun_ref / np.sin(rsun_obs * np.pi / 6.48e5)
    solarb0 = B0(hdr1["DATE_OBS"]).value
    hdr1["RSUN_REF"] = hdr1.get("RSUN_REF", rsun_ref)
    hdr1["RSUN_OBS"] = hdr1.get("RSUN_OBS", rsun_obs)
    hdr1["DSUN_OBS"] = hdr1.get("DSUN_OBS", dsun)
    hdr1["SOLAR_B0"] = hdr1.get("SOLAR_B0", solarb0)
    hdr1["HGLN_OBS"] = 0.0
    hdr1["HGLT_OBS"] = solarb0
    new_hdr = {}
    kw_to_copy = [
        "naxis",
        "naxis1",
        "naxis2",
        "date_obs",
        "time-obs",
        "ctime",
        "date_end",
        "crpix1",
        "crpix2",
        "crval1",
        "crval2",
        "cdelt1",
        "cdelt2",
        "cunit1",
        "cunit2",
        "ctype1",
        "ctype2",
        "dsun_obs",
        "rsun_ref",
        "rsun_obs",
        "solar_b0",
        "hgln_obs",
        "hglt_obs",
        "crota1",
        "crota2",
        "platescl",
        "keycomments",
    ]
    for kw in kw_to_copy:
        new_hdr[kw] = hdr1[kw]
    # Ensure CTYPE1 and CTYPE2 use the updated coordinate system
    new_hdr["L1_file1"] = filename1
    new_hdr["L1_file2"] = filename2
    
    extra_values, extra_comments = split_values_comments(extra_metadata)
    new_hdr.update(extra_values)
    create_date = datetime.now().ctime()
    new_hdr["history"] = f"Created by temperature_from_filter_ratio {create_date}\n"
    Thdr = new_hdr.copy()
    Thdr["BUNIT"] = "log10(K)"
    Thdr["history"] = (
        new_hdr["history"] + "Temperature derived using filter ratio method"
    )
    Tmap = Map(T_e, Thdr)
    Tmap.meta["keycomments"] = extra_comments
    Tmap.nickname = "Log Derived Temperature (K)"
    EMhdr = new_hdr.copy()
    EMhdr["BUNIT"] = r"log10(cm$^{{-3}}$)"
    EMhdr["history"] = (
        new_hdr["history"] + "Volume emission measure derived using filter ratio method"
    )
    EMmap = Map(EM, EMhdr)
    EMmap.meta["keycomments"] = extra_comments
    EMmap.nickname = r"Log Derived Volume E.M. (cm$^{{-3}}$)"
    Terrhdr = new_hdr.copy()
    Terrhdr["BUNIT"] = "log10(K)"
    Terrhdr["history"] = (
        new_hdr["history"] + "Temperature uncertainty derived using filter ratio method"
    )
    Terrmap = Map(T_error, Terrhdr)
    Terrmap.meta["keycomments"] = extra_comments
    Terrmap.nickname = "Log Derived Temperature Errors (K)"
    EMerrhdr = new_hdr.copy()
    EMerrhdr["BUNIT"] = r"log10(cm$^{{-3}}$)"
    EMerrhdr["history"] = (
        new_hdr["history"]
        + "Volume emission measure uncertainty derived using filter ratio method"
    )
    EMerrmap = Map(EMerror, EMerrhdr)
    EMerrmap.meta["keycomments"] = extra_comments
    EMerrmap.nickname = r"Log Derived V.E.M. Errors (cm$^{{-3}}$)"
    return Tmap, EMmap, Terrmap, EMerrmap


def measurement_to_filtername(measurement):
    """
    Convert filter combination as formatted in the measurement attribute for a
    SunPy |Map| to the filtername format as required for input to
    TemperatureResponseFundamental
    """
    fw1, fw2 = measurement.replace(" ", "_").split("-")
    if fw1 != "Open":
        filtername = fw1
        if fw2 != "Open":
            filtername = f"{fw1}/{fw2}"
    elif fw2 != "Open":
        filtername = fw2
    else:
        raise ValueError("Invalid filter values, both fw1 and fw2 are Open")
    return filtername


def split_values_comments(metadict):
    """
    Given a dictionary of keys and values where some of the values are tuples
    of the actual value and a string comment, return a dictionary with only
    the actual values and another one with the comments.
    """
    comment_dict = {}
    value_dict = {}
    for k, v in metadict.items():
        if isinstance(v, tuple):
            value_dict[k] = v[0]
            comment_dict[k] = v[1]
        else:
            value_dict[k] = v
    return value_dict, comment_dict
