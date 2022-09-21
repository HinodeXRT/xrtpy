import sys
import numpy as np
from scipy.integrate import simpson
from astropy.io import fits
from astropy.constants import c, h
from astropy import units as u
from xrtpy.response.temperature_response import TemperatureResponseFundamental
from sunpy.coordinates.sun import angular_radius


def xrt_teem(hdr1, data1, hdr2, data2, binfac=1, trange=None,
             no_threshold=False, te_err_threshold=0.5,
             photon_noise_threshold=0.2, mask=None, verbose=False):
    """
    Get coronal temperatures in an image using the filter ratio method.

    By default, this program uses the solar spectrum calculated with CHIANTI
    database ver. 6.0.1 (density: 10^9 [cm^-3], ionization equilibrium:
    chianti.ioneq, abundance: sun_coronal_ext).
    NOTE: Currently this is the only spectrum that can be used because of
    limitations in xrtpy. We expect this to change.

    Parameters:
    -----------
    hdr1 : fits header
        header for the first image

    data1 : array
        XRT *un-normalized* level1 data for first image (now normalized images
        are also accepted)

    hdr2 : fits header
        header for the second image (should use different filter than
        the first image)

    data2 : array
        XRT *un-normalized* level1 data for second image. Shape must match
        data1.

    binfac : integer, optional (default = 1)
        spatial binning factor

    trange : 2 element sequence containing floats [Optional]
        Range of log10(temperature) values to examine. Must be in order from
        lower to higher.

    te_err_threshold : float [Optional]
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
    te : 2-dimensional float array
        log10 of the derived temperature [K].
    em : 2-dimensional float array
        log10 of the derived volume emission measure [cm^-3].
    et : 2-dimensional float array
        error of log10 temperature [K].
    ee : 2-dimensional float array
        error of log10 volume emission measure [cm^-3].

    EXAMPLES:

       Using this function, you can derive the coronal temperature using
       filter ratio method.
       >> te, em, et, ee = xrt_teem(hdr1, data1, hdr2, data2)

       If you want to bin data in space to collect photons (to reduce photon
       noise), set bin to True as follows:
       In this case, data is binned as 3x3 in pixels first. After this,
       temperature is derived with binned data.
       >> te, em, et, ee = xrt_teem(hdr1, data1, hdr2, data2, binfac=3)

    NOTES:

       The returned value of pixels where temperature cannot be derived or
       error is greater than the threshold is set to 0.

       The details of the coronal-temperature-diagnostic capability of
       Hinode/XRT is described in
         Narukage et al. 2011, Solar Phys., 269, 169.
         http://adsabs.harvard.edu/doi/10.1007/s11207-010-9685-2
       and
         Narukage et al. 2013, Solar Phys., in press
         http://adsabs.harvard.edu/doi/10.1007/s11207-013-0368-7
       These two papers are the reference papers of this program.

    Modification History:

    IDL routine written by N.Narukage (NAOJ). Converted to python by Jonathan
    D. Slavin (SAO). See original IDL code for more details.
    """
    n1 = ('XRT_RENORMALIZE' in hdr1['HISTORY'])
    n2 = ('XRT_RENORMALIZE' in hdr2['HISTORY'])
    # This allows use of normalized data (contrary to original IDL code):
    if n1 or n2:
        if n1 and n2:
            data1 = data1*hdr1['EXPTIME']
            data2 = data2*hdr2['EXPTIME']
        elif n1:
            data1 = data1*hdr1['EXPTIME']
        else:
            data2 = data2*hdr2['EXPTIME']
    if data1.shape != data2.shape:
        raise ValueError('The input images must be the same size')

    if mask is None:
        mask = np.zeros_like(data1, dtype=bool)
    data1 = data1.astype(float)
    data2 = data2.astype(float)
    data1 = np.ma.masked_where(mask, data1)
    data2 = np.ma.masked_where(mask, data2)
    if binfac > 1:
        s = data1.shape
        ns = (s[0]//binfac, s[1]//binfac)
        rbs = (ns[0], binfac, ns[1], binfac)
        d1 = data1.reshape(rbs).sum(-1).sum(1)
        d2 = data2.reshape(rbs).sum(-1).sum(1)
        msk = mask.reshape(rbs).sum(-1).sum(1)
        # need to convert back to bool after summing - this makes a pixel
        # masked if any of the summed pixels is masked
        # if we want to mask only if all the pixels are masked then we could
        # use prod in place of sum above
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
    dmask = ((data1 <= 0.) | (data2 <= 0.))
    mask = (mask | dmask)

    fw1 = hdr1['EC_FW1_']
    fw2 = hdr1['EC_FW2_']
    if fw1 != 'Open':
        filt1 = fw1
        if fw2 != 'Open':
            filt1 = f'{fw1}/{fw2}'
    elif fw2 != 'Open':
        filt1 = fw2
    else:
        raise ValueError('Invalid filter values, both fw1 and fw2 are Open')
    date_obs1 = hdr1['DATE_OBS']
    tresp1 = TemperatureResponseFundamental(filt1, date_obs1)
    fw1 = hdr2['EC_FW1_']
    fw2 = hdr2['EC_FW2_']
    if fw1 != 'Open':
        filt2 = fw1
        if fw2 != 'Open':
            filt2 = f'{fw1}/{fw2}'
    elif fw2 != 'Open':
        filt2 = fw2
    else:
        raise ValueError('Invalid filter values, both fw1 and fw2 are Open')
    date_obs2 = hdr2['DATE_OBS']
    tresp2 = TemperatureResponseFundamental(filt2, date_obs2)

    spect1 = tresp1.spectra().value
    ea1 = tresp1.effective_area().value    # in cm^2
    spect2 = tresp2.spectra().value
    ea2 = tresp2.effective_area().value    # in cm^2
    wl = tresp1.channel_wavelength.value   # in Angstoms
    eVe = tresp1.ev_per_electron
    gain = tresp1.ccd_gain_right
    # old version:
    # flux1_alt = (spect1*ea1*e2dn).sum(axis=1)
    # flux2_alt = (spect2*ea2*e2dn).sum(axis=1)
    flux1 = tresp1.temperature_response().value
    flux2 = tresp2.temperature_response().value
    # Need to convert column emission measure to volume emission measure - to
    # do that we multiply the em by the area (in cm^2) of a sky pixel
    plate_scale = hdr1['PLATESCL']
    lsun = 6.95700E10/angular_radius(hdr1['DATE_OBS']).value
    em2vem = (plate_scale*lsun)**2
    flux1 /= em2vem
    flux2 /= em2vem
    ratio = flux1/flux2
    t = tresp1.CHIANTI_temperature.value
    t2mk = np.argmin(np.abs(t - 2.E6))
    rev_ratio = (ratio[t2mk] > 1.)

    # Note: this requires trange to be a 2 element iterable ordered low to
    # high (trange is log10 of temperature limits).
    if trange:
        in_trange = ((np.log10(t) >= trange[0]) &
                     (np.log10(t) <= trange[1]))
        if np.any(in_trange):
            t = t[in_trange]
            flux1 = flux1[in_trange]
            flux2 = flux2[in_trange]
            spect1 = spect1[in_trange, :]
            spect2 = spect2[in_trange, :]
        else:
            raise ValueError('The temperature response does not include'
                             ' any of the input temperatures in trange')

    if rev_ratio:
        data_ratio = (data2/hdr2['EXPTIME'])/(data1/hdr1['EXPTIME'])
        model_ratio = flux2/flux1
    else:
        data_ratio = (data1/hdr1['EXPTIME'])/(data2/hdr2['EXPTIME'])
        model_ratio = flux1/flux2
    ok_num = np.zeros(data_ratio.shape, dtype=int)
    ok_cnt = np.zeros(data_ratio.shape, dtype=int)
    for i, m in enumerate(model_ratio[1:]):
        n = ((data_ratio >= min(model_ratio[i], m)) &
             (data_ratio <= max(model_ratio[i], m)))
        if np.sum(n) > 0:
            ok_num[n] = i + 1
            ok_cnt[n] += 1

    ok_num[ok_cnt != 1] = 0
    a = np.abs(model_ratio[ok_num] - data_ratio)
    b = np.abs(model_ratio[np.maximum((ok_num - 1), 0)] - data_ratio)
    te = ((np.log10(t[ok_num])*b + np.log10(t[np.maximum((ok_num - 1), 0)])*a)
          / (a + b))
    te = np.ma.masked_where(te <= 0., te)

    # This masks out pixels where more than one temperature is consistent with
    # the observed flux ratio - so temperatures for which the filter ratio is
    # non-monotonic are effectively removed
    OK_pixel = ((ok_cnt == 1) & (~mask))
    te[~OK_pixel] = 0.

    # Note that te is the log10 of the electron temperature here
    DN = np.interp(te, np.log10(t), flux1, left=0., right=0.)
    DN = np.ma.masked_where(((DN <= 0.) | (te <= 0.)), DN)
    em = np.ma.log10(data1/(DN*hdr1['EXPTIME'])) - np.log10(binfac**2)
    em[~OK_pixel] = 0.

    # derive Te error, em error:
    dlnR_dlnT = np.abs(deriv(np.log(t), np.log(model_ratio)))
    dlnR_dlnT = np.interp(te, np.log10(t), dlnR_dlnT, left=0., right=0.)
    dlnR_dlnT = np.ma.masked_where(((dlnR_dlnT == 0.) | (te <= 0.)), dlnR_dlnT)

    e2dn = (h.to(u.eV*u.s)*c.to(u.angstrom/u.s)/(wl*eVe*gain)).value
    K1 = np.array([(s1*ea1*e2dn**2).sum()/(s1*ea1*e2dn).sum() for s1 in
                  spect1])
    K2 = np.array([(s2*ea2*e2dn**2).sum()/(s2*ea2*e2dn).sum() for s2 in
                  spect2])
    K1 = np.interp(te, np.log10(t), K1, left=0., right=0.)
    K2 = np.interp(te, np.log10(t), K2, left=0., right=0.)
    K1 = np.ma.masked_where(((K1 == 0.) | (te <= 0.)), K1)
    K2 = np.ma.masked_where(((K2 == 0.) | (te <= 0.)), K2)

    dlnf1_dlnT = deriv(np.log(t), np.log(flux1))
    dlnf1_dlnT = np.interp(te, np.log10(t), dlnf1_dlnT, left=0., right=0.)
    dlnf1_dlnT = np.ma.masked_where(((dlnf1_dlnT == 0.) | (te <= 0.)),
                                    dlnf1_dlnT)

    dlnf2_dlnT = deriv(np.log(t), np.log(flux2))
    dlnf2_dlnT = np.interp(te, np.log10(t), dlnf2_dlnT, left=0., right=0.)
    dlnf2_dlnT = np.ma.masked_where(((dlnf2_dlnT == 0.) | (te <= 0.)),
                                    dlnf2_dlnT)

    et = np.ma.log10(np.sqrt(K1/data1 + K2/data2)/dlnR_dlnT) + te

    ee = np.ma.log10(np.ma.sqrt(dlnf2_dlnT**2 * K1/data1
                     + dlnf1_dlnT**2*K2/data2)/dlnR_dlnT) + em
    te = te.filled(0.)
    em = em.filled(0.)
    et = et.filled(0.)
    ee = ee.filled(0.)

    if not no_threshold:
        OK_wothr = OK_pixel.copy()
        Kd1 = np.sqrt(K1/data1)
        Kd1 = Kd1.filled(0.)
        Kd2 = np.sqrt(K2/data2)
        Kd2 = Kd2.filled(0.)
        tthr = ((et - te) <= np.log10(te_err_threshold))
        k1thr = (Kd1 <= photon_noise_threshold)
        k2thr = (Kd2 <= photon_noise_threshold)
        OK_pixel = (OK_pixel & ((et - te) <= np.log10(te_err_threshold))
                    & (Kd1 <= photon_noise_threshold)
                    & (Kd2 <= photon_noise_threshold))
        if verbose:
            print(f'number of pixels ruled out by threshold = '
                  f'{np.sum(~OK_pixel)}')
            print(f'number of pixels ruled out by te errors = {np.sum(~tthr)}')
            print(f'number of pixels ruled out by d1 noise  = {np.sum(~k1thr)}')
            print(f'number of pixels ruled out by d2 noise  = {np.sum(~k2thr)}')
            print(f'number of bad pixels before threshold   = '
                  f'{np.sum(~OK_wothr)}')
        mask = (mask | ~OK_pixel)
        te[mask] = 0.
        em[mask] = 0.
        et[mask] = 0.
        ee[mask] = 0.

        if verbose:
            print('from xrt_teem:')
            print(f'Examined Te range: {t.min()} - {t.max()} K')
            print(f'Applied thresholds: - Te error < {te_err_threshold*100.} %')
            print(f'                    - Photon noise < '
                  f'{photon_noise_threshold*100.} %')
    else:
        if verbose:
            print('from xrt_teem:')
            print(f'Examined Te range: {t.min()} - {t.max()} K')
            print('No thresholds applied')
    return te, em, et, ee


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
    dydx1 = y0*x12/(x01*x02) + y1*(1./x12 - 1./x01) - y2*x01/(x02*x12)
    dydx0 = (y0[0]*(x01[0]
             + x02[0])/(x01[0]*x02[0])
             - y1[0]*x02[0]/(x01[0]*x12[0])
             + y2[0]*x01[0]/(x02[0]*x12[0]))
    dydxN = (-y0[-1]*x12[-1]/(x01[-1]*x02[-1])
             + y1[-1]*x02[-1]/(x01[-1]*x12[-1])
             - y2[-1]*(x02[-1]
             + x12[-1])/(x02[-1]*x12[-1]))
    return np.append(np.insert(dydx1, 0, dydx0), dydxN)
