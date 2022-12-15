import numpy as np
import sys

from astropy import units as u
from astropy.constants import c, h
from sunpy.coordinates.sun import angular_radius

from xrtpy.response.temperature_response import TemperatureResponseFundamental


def xrt_teem(
    hdr1,
    data1,
    hdr2,
    data2,
    binfac=1,
    Trange=None,
    no_threshold=False,
    Te_err_threshold=0.5,
    photon_noise_threshold=0.2,
    mask=None,
    verbose=False,
):
    r"""
    Get coronal temperatures and emission measures from a pair of images using
    the filter ratio method.

    Currently this program uses the solar spectrum calculated with CHIANTI
    database ver. 6.0.1 (density: :math:`10^9` cm\ :sup:`-3`\ , ionization
    equilibrium: chianti.ioneq, abundance: sun_coronal_ext), because this
    is the only spectrum available in xrtpy. We expect this to change.

    Parameters:
    -----------
    hdr1 : fits header
        header for the first image

    data1 : array
        XRT level1 data for first image. If the image is normalized, then it
        is assumed that the un-normalized image can be recovered by
        multiplying by the exposure time (exposure time must be available). It
        is also assumed that the header history will contain the string
        XRT_RENORMALIZE if the image has been normalized.

    hdr2 : fits header
        header for the second image (must use different filters from the first
        image)

    data2 : array
        XRT level1 data for second image. Shape must match data1.

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

    mask : Boolean array of shape of images [Optional]
        If provided, masks out parts of the images from the analysis. Note:
        pixels to be masked out should be True, unmasked should be False

    no_threshold : Boolean [Optional]
        If True, no thresholds are set. (default = False)

    verbose : Boolean [Optional]
        If True, information is printed


    Returns:
    --------
    T_e : 2-dimensional float array
        log10 of the derived electron temperature [K].
    EM : 2-dimensional float array
        log10 of the derived volume emission measure [cm^-3].
    T_error : 2-dimensional float array
        error of log10 temperature [K].
    EM_error : 2-dimensional float array
        error of log10 volume emission measure [cm^-3].

    Examples
    --------
    Using this function, you can derive the coronal temperature using
    filter ratio method.

    >>> T_e, EM, Terror, EMerror = xrt_teem(hdr1, data1, hdr2, data2) # doctest: +SKIP

    If you want to bin the image data in space to reduce photon noise, set
    binfac to the factor by which you want to bin.  For example to bin the
    data by a factor of 3 do:

    >>> T_e, EM, Terror, EMerror = xrt_teem(hdr1, data1, hdr2, data2, binfac=3) # doctest: +SKIP

    The data is binned first and then the temperature is derived. Note that
    the image size is not reduced, but pixels within 3×3 squares are set to
    the same value, which results from averaging over those pixels.

    Notes
    -----
    The returned values of pixels where the temperature cannot be derived
    or the error is greater than the threshold or the photon noise is
    greater than the threshold are set to 0. The EM for those pixels is
    also set to 0.

    The details of the coronal-temperature-diagnostic capability of
    Hinode/XRT is described in
    Narukage et al. 2011, Solar Phys., 269, 169.
    http://adsabs.harvard.edu/doi/10.1007/s11207-010-9685-2
    and
    Narukage et al. 2013, Solar Phys.,
    http://adsabs.harvard.edu/doi/10.1007/s11207-013-0368-7
    These two papers are the reference papers of this program.

    Modification History:

    IDL routine written by N.Narukage (NAOJ). Converted to python by Jonathan
    D. Slavin (SAO). See original IDL code for more details.
    """

    n1 = "XRT_RENORMALIZE" in hdr1["HISTORY"]
    n2 = "XRT_RENORMALIZE" in hdr2["HISTORY"]
    # This allows use of normalized data (contrary to original IDL code):
    if n1 or n2:
        if n1 and n2:
            data1 = data1 * hdr1["EXPTIME"]
            data2 = data2 * hdr2["EXPTIME"]
        elif n1:
            data1 = data1 * hdr1["EXPTIME"]
        else:
            data2 = data2 * hdr2["EXPTIME"]
    if data1.shape != data2.shape:
        raise ValueError("The input images must be the same size")

    if mask is None:
        mask = np.zeros_like(data1, dtype=bool)
    data1 = data1.astype(float)
    data2 = data2.astype(float)
    data1 = np.ma.masked_where(mask, data1)
    data2 = np.ma.masked_where(mask, data2)
    if binfac > 1:
        s = data1.shape
        ns = (s[0] // binfac, s[1] // binfac)
        rbs = (ns[0], binfac, ns[1], binfac)
        # sums the data in binfac x binfac sized regions
        d1 = data1.reshape(rbs).sum(-1).sum(1)
        d2 = data2.reshape(rbs).sum(-1).sum(1)
        # this makes a pixel masked if any of the summed pixels is masked
        # if we want to mask only if all the pixels are masked then we could
        # use prod in place of sum here
        msk = mask.reshape(rbs).sum(-1).sum(1)
        # need to convert back to bool after summing
        msk = msk.astype(bool)
        # This restores the image to the size of the original images as in the
        # IDL code:
        d1tmp = np.zeros_like(data1)
        d2tmp = np.zeros_like(data2)
        mtmp = np.zeros_like(mask, dtype=bool)
        for i in range(binfac):
            for j in range(binfac):
                d1tmp[i::binfac, j::binfac] = d1[:]
                d2tmp[i::binfac, j::binfac] = d2[:]
                mtmp[i::binfac, j::binfac] = msk[:]
        data1 = d1tmp
        data2 = d2tmp
        mask = mtmp

    # input mask for data should be False in parts of the images to be used and
    # True in places we want to mask out
    # Here we additionally mask out pixels in which the data in either image
    # is <= 0
    dmask = (data1 <= 0.0) | (data2 <= 0.0)
    mask = mask | dmask

    fw1 = hdr1["EC_FW1_"]
    fw2 = hdr1["EC_FW2_"]
    if fw1 != "Open":
        filt1 = fw1
        if fw2 != "Open":
            filt1 = f"{fw1}/{fw2}"
    elif fw2 != "Open":
        filt1 = fw2
    else:
        raise ValueError("Invalid filter values, both fw1 and fw2 are Open")
    date_obs1 = hdr1["DATE_OBS"]
    tresp1 = TemperatureResponseFundamental(filt1, date_obs1)
    fw1 = hdr2["EC_FW1_"]
    fw2 = hdr2["EC_FW2_"]
    if fw1 != "Open":
        filt2 = fw1
        if fw2 != "Open":
            filt2 = f"{fw1}/{fw2}"
    elif fw2 != "Open":
        filt2 = fw2
    else:
        raise ValueError("Invalid filter values, both fw1 and fw2 are Open")
    date_obs2 = hdr2["DATE_OBS"]
    tresp2 = TemperatureResponseFundamental(filt2, date_obs2)

    if filt1 == filt2:
        raise ValueError("Filters for the two images cannot be the same")

    flux1 = tresp1.temperature_response().value
    flux2 = tresp2.temperature_response().value

    # Need to convert column emission measure to volume emission measure - to
    # do that we multiply the EM by the area (in cm^2) of a sky pixel
    plate_scale = hdr1["PLATESCL"]
    lsun = 6.95700e10 / angular_radius(hdr1["DATE_OBS"]).value
    em2vem = (plate_scale * lsun) ** 2
    flux1 /= em2vem
    flux2 /= em2vem
    ratio = flux1 / flux2
    Tmodel = tresp1.CHIANTI_temperature.value
    t2mk = np.argmin(np.abs(Tmodel - 2.0e6))
    logTmodel = np.log10(Tmodel)
    rev_ratio = ratio[t2mk] > 1.0
    spect1 = tresp1.spectra().value
    spect2 = tresp2.spectra().value
    effarea1 = tresp1.effective_area().value  # in cm^2
    effarea2 = tresp2.effective_area().value  # in cm^2

    # Note: this requires Trange to be a 2 element iterable ordered low to
    # high (Trange is log10 of temperature limits).
    if Trange:
        in_trange = (logTmodel >= Trange[0]) & (logTmodel <= Trange[1])
        if np.any(in_trange):
            Tmodel = Tmodel[in_trange]
            logTmodel = np.log10(Tmodel)
            flux1 = flux1[in_trange]
            flux2 = flux2[in_trange]
            spect1 = spect1[in_trange, :]
            spect2 = spect2[in_trange, :]
        else:
            raise ValueError(
                "The temperature response does not include"
                " any of the input temperatures in Trange"
            )

    if rev_ratio:
        data_ratio = (data2 / hdr2["EXPTIME"]) / (data1 / hdr1["EXPTIME"])
        model_ratio = flux2 / flux1
    else:
        data_ratio = (data1 / hdr1["EXPTIME"]) / (data2 / hdr2["EXPTIME"])
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
    EM = np.ma.log10(data1 / (DN * hdr1["EXPTIME"])) - np.log10(binfac**2)
    EM[~ok_pixel] = 0.0

    # derive T_e error, EM error:
    dlnR_dlnT_mod = np.abs(deriv(np.log(Tmodel), np.log(model_ratio)))
    dlnR_dlnT = np.interp(T_e, logTmodel, dlnR_dlnT_mod, left=0.0, right=0.0)
    dlnR_dlnT = np.ma.masked_where(((dlnR_dlnT == 0.0) | (T_e <= 0.0)), dlnR_dlnT)

    wvl = tresp1.channel_wavelength
    eVe = tresp1.ev_per_electron
    gain = tresp1.ccd_gain_right
    # (h*c/lambda) * 1/(eV per electron) * 1/gain
    e2dn = (h.to(u.eV * u.s) * c.to(u.angstrom / u.s) / (wvl * eVe * gain)).value
    dwvl = wvl[1:] - wvl[:-1]
    dwvl = np.append(dwvl, dwvl[-1]).value

    # These sums should be converted to integrals if the spectrum includes
    # line broadening (as is the default in ChiantiPy).
    K1_mod = np.array(
        [
            (s1 * effarea1 * e2dn**2 * dwvl).sum()
            / (s1 * effarea1 * e2dn * dwvl).sum()
            for s1 in spect1
        ]
    )
    K2_mod = np.array(
        [
            (s2 * effarea2 * e2dn**2 * dwvl).sum()
            / (s2 * effarea2 * e2dn * dwvl).sum()
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

    T_error = np.ma.log10(np.sqrt(K1 / data1 + K2 / data2) / dlnR_dlnT) + T_e

    EMerror = (
        np.ma.log10(
            np.ma.sqrt(dlnf2_dlnT**2 * K1 / data1 + dlnf1_dlnT**2 * K2 / data2)
            / dlnR_dlnT
        )
        + EM
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
        if verbose:
            print(f"number of pixels ruled out by threshold = " f"{np.sum(~ok_pixel)}")
            print(f"number of pixels ruled out by T_e errors = {np.sum(~tthr)}")
            print(f"number of pixels ruled out by d1 noise  = {np.sum(~k1thr)}")
            print(f"number of pixels ruled out by d2 noise  = {np.sum(~k2thr)}")
            print(f"number of bad pixels before threshold   = " f"{np.sum(~ok_wothr)}")
        mask = mask | ~ok_pixel
        T_e[mask] = 0.0
        EM[mask] = 0.0
        T_error[mask] = 0.0
        EMerror[mask] = 0.0

        if verbose:
            print("from xrt_teem:")
            print(f"Examined T_e range: {Tmodel.min()} - {Tmodel.max()} K")
            print(f"Applied thresholds: - T_e error < {Te_err_threshold*100.} %")
            print(
                f"                    - Photon noise < "
                f"{photon_noise_threshold*100.} %"
            )
    else:
        if verbose:
            print("from xrt_teem:")
            print(f"Examined T_e range: {Tmodel.min()} - {Tmodel.max()} K")
            print("No thresholds applied")
    return T_e, EM, T_error, EMerror


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
