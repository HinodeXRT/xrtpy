"""
Functionality for deconvolving XRT image data with the point spread function.
"""
__all__ = ["xrt_deconvolve"]

import numpy as np

from datetime import datetime
from numpy.fft import fft2, fftshift, ifft2
from sunpy.data import manager
from sunpy.image.resample import resample
from sunpy.image.transform import affine_transform
from sunpy.map import Map


@manager.require(
    "PSF560.fits",
    [
        "https://sohoftp.nascom.nasa.gov/solarsoft/hinode/xrt/idl/util/XRT20170324_151721.0.PSF560.fits"
    ],
    "0eaa5da6fb69661e7f46d1f0c463e4b3b1745426a399a4fbc53fc0c0ae87dd0d",
)
@manager.require(
    "PSF1000.fits",
    [
        "https://sohoftp.nascom.nasa.gov/solarsoft/hinode/xrt/idl/util/XRT20170324_161721.0.PSF1000.fits"
    ],
    "95590a7174692977a2f111b932811c9c7ae105a59b93bfe6c96fba862cefacf1",
)
def xrt_deconvolve(image_map, niter=5, verbose=False, psf1keV=False):
    """
    Use the XRT mirror model point spread function (PSF) to deconvolve an XRT
    image

    Parameters:
    -----------
    image_map : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the input XRT image

    niter : integer, optional (default = 5)
        number of iterations to perform

    verbose : boolean, optional
        if True, extra information is printed

    psf1keV : boolean, optional
        if True, use the 1.0 keV PSF instead of the default 560 eV PSF.

    Returns:
    --------
    deconv_map : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the output deconvolved image

    """

    if psf1keV:
        used_psf = manager.get("PSF560.fits")
    else:
        used_psf = manager.get("PSF1000.fits")

    psf_map = Map(used_psf)
    psf_meta = psf_map.meta

    if verbose:
        print(f"XRT_DECONVOLVE: Using PSF in\n{used_psf}")

    image_meta = image_map.meta
    data_chip_sum = image_meta["chip_sum"]
    # place PSF & data on same chip sum.

    if psf_meta["chip_sum"] != data_chip_sum:
        if verbose:
            print(
                "XRT_DECONVOLVE: Input data and PSF have different"
                " chip sums. Binning PSF..."
            )
        psf_map = rebin_psf(psf_map, image_map.meta)

    data = np.clip(image_map.data, 0.0, None)

    deconvolve_hist = (
        f"(XRT_DECONVOLVE) Deconvolved with {used_psf}, "
        f"bin={data_chip_sum}, niter={niter}."
    )

    # if input image size is smaller than the psf (2048x2048), put the
    # image in the center of zero array the same size as the psf
    extract_data = False
    naxis1 = image_map.meta["naxis1"]
    naxis2 = image_map.meta["naxis2"]
    psf_naxis1 = psf_map.meta["naxis1"]
    psf_naxis2 = psf_map.meta["naxis2"]
    if (naxis1 * data_chip_sum != psf_naxis1 * psf_meta["chip_sum"]) or (
        naxis2 * data_chip_sum != psf_naxis1 * psf_meta["chip_sum"]
    ):
        extract_data = True
        if verbose:
            print("Input data not same size as PSF. Dropping image" " in zero array.")
        tmp_data = np.zeros((psf_naxis1, psf_naxis2))
        ddx = naxis1 // 2
        ddy = naxis2 // 2
        xcen = psf_naxis1 // 2 - 1
        ycen = psf_naxis2 // 2 - 1
        tmp_data[xcen - ddx : xcen + ddx, ycen - ddy : ycen + ddy] = data
    else:
        tmp_data = data
    tmp_deconv = richardson_lucy_deconvolution(tmp_data, psf_map.data, num_iter=5)

    if extract_data:
        tmp_deconv = tmp_deconv[xcen - ddx : xcen + ddx, ycen - ddy : ycen + ddy]
    deconv_data = np.minimum(tmp_deconv, 2500.0)

    date = datetime.now().ctime()
    added_hist = f"{__name__}: ({date}) " + deconvolve_hist
    deconv_meta = image_meta
    deconv_meta["history"] = image_meta["history"] + added_hist
    return Map(deconv_data, deconv_meta)


def xrt_fft_2dim_convolution(image1, image2, correlation=False):
    """
    Convolve (or optionally correlate) two images
    """
    if correlation:
        fftres = ifft2(fft2(image1) * np.conj(fft2(image2)))
    else:
        fftres = ifft2(fft2(image1) * fft2(image2))
    return fftshift(fftres)


def richardson_lucy_deconvolution(image, psf, num_iter=5):
    """
    Use the Richardson-Lucy algorithm to deconvolve an image.
    """
    psfnorm = xrt_fft_2dim_convolution(psf, np.ones_like(psf))
    ohat = np.cdouble(image)
    for i in range(num_iter):
        ihat = xrt_fft_2dim_convolution(psf, ohat)
        ohat *= xrt_fft_2dim_convolution(image / ihat, psf, correlation=True) / psfnorm
    return np.abs(ohat)


def rebin_psf(psf_map, image_meta):
    """
    Rebin the point spread function (psf) to match the dimensions of an image.
    It's assumed that the image is smaller than the psf array

    Parameters:
    -----------
    psf_map : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the point spread function

    image_meta : ~sunpy.util.metadata.MetaDict
        meta data for an image

    Returns:
    --------
    downsampled map : ~sunpy.map.sources.hinode.XRTMap
        |Map| for the psf that has been binned to match the pixel size of the
        data
    """
    psf_image = psf_map.data
    psf_meta = psf_map.meta
    center = (np.array(psf_image.shape)[::-1] - 1) / 2.0
    new_center = (center[0] - 0.5, center[1] - 0.5)
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])
    psf_shifted = affine_transform(
        psf_image,
        rmatrix,
        recenter=True,
        image_center=new_center,
        order=3,
        method="scikit-image",
        missing=0.0,
    )
    psf_shifted = psf_shifted / psf_shifted.sum()
    dims = (
        psf_meta["naxis1"] // image_meta["chip_sum"],
        psf_meta["naxis2"] // image_meta["chip_sum"],
    )
    downsampled = resample(psf_shifted, dims, method="spline", center=True)
    psf = downsampled / downsampled.sum()
    psf_meta["chip_sum"] = image_meta["chip_sum"]
    psf_meta["naxis1"] = psf.shape[0]
    psf_meta["naxis2"] = psf.shape[1]
    return Map(psf, psf_meta)
