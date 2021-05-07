import logging

import numpy as np
import numpy.ma as ma
from scipy import ndimage, signal
from multiprocessing import Pool, cpu_count

import topo_descriptors.helpers as hlp
from topo_descriptors import CFG

logger = logging.getLogger(__name__)


def compute_tpi(dem_da, scales, smth_factors=None, ind_nans=[], crop=None):
    """Wrapper to 'tpi' function to launch computations for all scales and save
    outputs as netCDF files.

    Parameters
    ----------
    dem_da : xarray DataArray representing the DEM and its grid coordinates.
    scales : scalar or list of scalars
        Scale(s) in meters on which we want to compute the TPI.
        Corresponds to the diameter of the kernel used to compute it.
    smth_factors (optional) : scalar or None or list with a combination of both.
        Fraction(s) of the scale(s) at which the DEM is smoothed first (with a
        gaussian filter). If None (default), no prior smoothing is performed.
        If a scalar, the same fraction is used to determine the smoothing scale
        of all specified scales. If a list, must match the length of arg 'scales'.
    ind_nans (optional) : tuple of two 1D arrays
        Contains the (row, column) indices of the NaNs in the original DEM to be
        reassigned after computations. NaNs in the original DEM should be
        interpolated prior computations as they propagate in convolutions with
        the fast fourrier transform method (scipy.signal.convolve).
    crop (optional) : dict
        If specified the outputs are cropped to the given extend. Keys should be
        the coordinates labels of dem_da and values should be slices of [min,max]
        extend. Default is None.

    See also
    --------
    tpi, _tpi_kernel
    """

    hlp.check_dem(dem_da)
    logger.info(f"***Starting TPI computation for scales {scales} meters***")
    if not hasattr(scales, "__iter__"):
        scales = [scales]
    if not hasattr(smth_factors, "__iter__"):
        smth_factors = [smth_factors] * len(scales)

    scales_pxl, _ = hlp.scale_to_pixel(scales, dem_da)
    sigmas = hlp.get_sigmas(smth_factors, scales_pxl)

    for idx, scale_pxl in enumerate(scales_pxl):
        logger.info(
            f"Computing scale {scales[idx]} meters with smoothing factor"
            f" {smth_factors[idx]} ..."
        )
        name = _tpi_name(scales[idx], smth_factors[idx])
        array = tpi(dem=dem_da.values, size=scale_pxl, sigma=sigmas[idx])

        array[ind_nans] = np.nan
        hlp.to_netcdf(array, dem_da.coords, name, crop)
        del array


@hlp.timer
def tpi(dem, size, sigma=None):
    """Compute the TPI over a digital elevation model. The TPI represents
    the elevation difference of a pixel relative to its neighbors.

    Parameters
    ----------
    dem : array representing the DEM.
    size : int
        Size of the kernel for the convolution. Represents the diameter (i.e. scale)
        in pixels at which the TPI is computed.
    sigma (optional) : scalar
        If provided, the DEM is first smoothed with a gaussian filter of standard
        deviation sigma (in pixel size).

    Returns
    -------
    array with TPI values

    See also
    --------
    scipy.signal.convolve, scipy.ndimage.gaussian_filter
    """

    kernel = _tpi_kernel(size)
    if sigma:
        dem = ndimage.gaussian_filter(dem, sigma)

    conv = signal.convolve(
        dem, kernel, mode="same"
    )  # ndimage.convolve(dem, kernel, mode='reflect')
    return dem - conv / np.sum(kernel)


def _tpi_name(scale, smth_factor):
    """Return name for the array in output of the tpi function"""

    add = f"_SMTHFACT{smth_factor:.3g}" if smth_factor else ""
    return f"TPI_{scale}M{add}"


def _tpi_kernel(size):
    """Generate a circular kernel to compute TPI.

    Parameters
    ----------
    size : int
        Size of the kernel.

    Returns
    -------
    2-D array with kernel values
    """

    middle = int(size / 2)
    if size < 5:
        kernel = np.ones((size, size), dtype=np.float32)
    else:
        xx, yy = np.mgrid[:size, :size]
        circle = (xx - middle) ** 2 + (yy - middle) ** 2
        kernel = np.asarray(circle <= (middle ** 2), dtype=np.float32)

    kernel[middle, middle] = 0
    return kernel


def compute_valley_ridge(
    dem_da,
    scales,
    mode,
    flat_list=[0, 0.15, 0.3],
    smth_factors=None,
    ind_nans=[],
    crop=None,
):
    """Wrapper to 'valley_ridge' function to launch computations for all scales
    and save outputs as netCDF files.

    Parameters
    ----------
    dem_da : xarray DataArray representing the DEM and its grid coordinates.
    scales : scalar or list of scalars
        Scale(s) in meters over which we want to compute the valley or the ridge
        index. Corresponds to the size of the squared kernel used to compute it.
    mode : {valley, ridge}
        Whether to compute the valley or ridge index.
    flat_list (optional) : list of floats in [0,1[
        Fractions of flat along the center ligne of the V-shape kernels. A certain
        amount of flat is use to approximate the shape of glacial valleys.
        Default is [0, 0.15, 0.3].
    smth_factors (optional) : scalar or None or list with a combination of both.
        Fraction(s) of the scale(s) at which the DEM is smoothed first (with a
        gaussian filter). If None (default), no prior smoothing is performed.
        If a scalar, the same fraction is used to determine the smoothing scale
        of all specified scales. If a list, must match the length of arg 'scales'.
    ind_nans (optional) : tuple of two 1D arrays
        Contains the (row, column) indices of the NaNs in the original DEM to be
        reassigned after computations. NaNs in the original DEM should be
        interpolated prior computations as they propagate in convolutions with
        the fast fourrier transform method (scipy.signal.convolve).
    crop (optional) : dict
        If specified the outputs are cropped to the given extend. Keys should be
        the coordinates labels of dem_da and values should be slices of [min,max]
        extend. Default is None.

    See also
    --------
    valley_ridge, _valley_kernels, _valley_ridge_names
    """

    hlp.check_dem(dem_da)
    logger.info(f"***Starting {mode} index computation for scales {scales} meters***")
    if not hasattr(scales, "__iter__"):
        scales = [scales]
    if not hasattr(smth_factors, "__iter__"):
        smth_factors = [smth_factors] * len(scales)

    scales_pxl, _ = hlp.scale_to_pixel(scales, dem_da)
    sigmas = hlp.get_sigmas(smth_factors, scales_pxl)

    pool = Pool(processes=min(len(scales_pxl), cpu_count()))
    for idx, scale_pxl in enumerate(scales_pxl):
        logger.info(
            f"Computing scale {scales[idx]} meters with smoothing factor"
            f" {smth_factors[idx]} ..."
        )
        names = _valley_ridge_names(scales[idx], mode, smth_factors[idx])
        pool.apply_async(
            _valley_ridge_wrap,
            args=(
                dem_da,
                scale_pxl,
                mode,
                flat_list,
                sigmas[idx],
                names,
                ind_nans,
                crop,
            ),
        )

    pool.close()
    pool.join()


def _valley_ridge_wrap(dem_da, size, mode, flat_list, sigma, names, ind_nans, crop):
    """Wrapper to valley_ridge and hlp.to_netcdf functions to ease the parallelization
    of the different scales"""

    arrays = valley_ridge(dem_da.values, size, mode, flat_list, sigma)
    for array, name in zip(arrays, names):
        array[ind_nans] = np.nan
        hlp.to_netcdf(array, dem_da.coords, name, crop)


@hlp.timer
def valley_ridge(dem, size, mode, flat_list=[0, 0.15, 0.3], sigma=None):
    """Compute the valley or ridge index over a digital elevation model.
    The valley/ridge index highlights valley or ridges  at various scales.

    Parameters
    ----------
    dem : array representing the DEM.
    size : int
        Size of the kernel for the convolution. Represents the width (i.e. scale)
        in pixels of the valleys we are trying to highlight.
    mode : {valley, ridge}
        Whether to compute the valley or ridge index.
    flat_list (optional) : list of floats in [0,1[
        Fractions of flat along the center ligne of the V-shape kernels. A certain
        amount of flat is use to approximate the shape of glacial valleys.
        Default is [0, 0.15, 0.3].
    sigma (optional) : scalar
        If provided, the DEM is first smoothed with a gaussian filter of standard
        deviation sigma (in pixel size).

    Returns
    -------
    list of two arrays :
        First element is the norm and second the direction of the valley or ridge
        index. The direction in degrees is defined from 0° to 179°, increasing
        clockwise. W-E oriented valleys have a direction close to 0° or 180°,
        while S-N oriented valleys have a direction close to 90°.

    See also
    --------
    scipy.signal.convolve, scipy.ndimage.gaussian_filter
    """

    if mode not in ("valley", "ridge"):
        raise ValueError(f"Unknown mode {mode!r}")

    if sigma:
        dem = ndimage.gaussian_filter(dem, sigma)

    dem = (dem - dem.mean()) / dem.std()
    n_y, n_x = dem.shape
    dem = np.broadcast_to(dem, (len(flat_list), n_y, n_x))
    angles = np.arange(0, 180, dtype=np.float32)
    index_norm = np.zeros((n_y, n_x), dtype=np.float32) - np.inf
    index_dir = np.empty((n_y, n_x), dtype=np.float32)

    if mode == "ridge":
        kernels = _ridge_kernels(size, flat_list)
    else:
        kernels = _valley_kernels(size, flat_list)

    for angle in angles:  # 0° = E-W valleys, 90° = S-N valleys

        kernels_rot = _rotate_kernels(kernels, angle)
        conv = signal.convolve(dem, kernels_rot, mode="same")
        conv = np.max(conv, axis=0)
        bool_greater = conv > index_norm
        index_norm[bool_greater] = conv[bool_greater]
        index_dir[bool_greater] = angle
        del bool_greater
        if angle % 20 == 0:
            logger.info(f"angle {int(angle)}/180 finished")

    index_norm = np.ndarray.clip(index_norm, min=0)
    return [index_norm, index_dir]


def _valley_ridge_names(scale, mode, smth_factor):
    """Return names for the arrays in output of the valley_ridge function"""

    add = f"_SMTHFACT{smth_factor:.3g}" if smth_factor else ""
    name_norm = f"{mode}_NORM_{scale}M{add}"
    name_dir = f"{mode}_DIR_{scale}M{add}"

    return [name_norm, name_dir]


def _valley_kernels(size, flat_list):
    """Generate normalized V-shape and U-shape kernels to compute valley index.

    Parameters
    ----------
    size : int
        Size of the kernel.
    flat_list : list of floats in [0,1[
        Fractions of flat along the center ligne of the V-shape kernels. A certain
        amount of flat is use to approximate the shape of glacial valleys.

    Returns
    -------
    3-D array with 2-D kernels for each specified flat fraction.
    """

    middle = int(np.floor(size / 2))
    kernel_tmp = np.broadcast_to(np.arange(0, middle + 1), (size, middle + 1)).T
    kernel_tmp = np.concatenate(
        (np.flip(kernel_tmp[1:, :], axis=0), kernel_tmp), axis=0
    )
    kernel_tmp = np.asarray(kernel_tmp, dtype=np.float32)
    kernels = np.broadcast_to(kernel_tmp, (len(flat_list), size, size)).copy()

    for ind, flat in enumerate(flat_list):
        halfwidth = int(np.floor(np.floor(size * flat / 2) + 0.5))
        kernels[ind, middle - halfwidth : middle + halfwidth + 1, :] = kernels[
            ind, middle - halfwidth, 0
        ]
        kernels = (kernels - np.mean(kernels, axis=(1, 2), keepdims=True)) / np.std(
            kernels, axis=(1, 2), keepdims=True
        )

    return kernels


def _ridge_kernels(size, flat_list):
    """Generate normalized flipped V-shape and U-shape kernels to compute ridge index.

    Parameters
    ----------
    size : int
        Size of the kernel.
    flat_list : list of floats in [0,1[
        Fractions of flat along the center ligne of the V-shape kernels. A certain
        amount of flat is use to approximate the shape of glacial valleys.

    Returns
    -------
    3-D array with 2-D kernels for each specified flat fraction.
    """

    return _valley_kernels(size, flat_list) * -1


def _rotate_kernels(kernel, angle):
    """Rotate a 3-D kernel in the plane given by the last two axes"""

    kernels_rot = ndimage.rotate(
        kernel, angle, axes=(1, 2), reshape=True, order=2, mode="constant", cval=-9999
    )
    kernels_rot = ma.masked_array(kernels_rot, mask=kernels_rot == -9999)
    kernels_rot = (
        kernels_rot - np.mean(kernels_rot, axis=(1, 2), keepdims=True)
    ) / np.std(kernels_rot, axis=(1, 2), keepdims=True)
    return ma.MaskedArray.filled(kernels_rot, 0).astype(np.float32)


def compute_gradient(dem_da, scales, sig_ratios=1, ind_nans=[], crop=None):
    """Wrapper to 'gradient' function to launch computations for all scales
    and save outputs as netCDF files.

    Parameters
    ----------
    dem_da : xarray DataArray representing the DEM and its grid coordinates.
    scales : scalar or list of scalars
        Scale(s) in meters on which we want to compute the the valley or ridge
        index. Corresponds to the size of the squared kernel used to compute it.
    sig_ratios (optional) : scalar or list of scalars.
        Ratios w.r.t scales to define the smoothing scale in the perpendicular
        direction of the directional derivatives. If a list, must match the length
        of arg 'scales'. Default is 1 (i.e. same smoothing on both directions for
        all scales)
    ind_nans (optional) : tuple of two 1D arrays
        Contains the (row, column) indices of the NaNs in the original DEM to be
        reassigned after computations. NaNs in the original DEM should be
        interpolated prior computations as they propagate in convolutions.
    crop (optional) : dict
        If specified the outputs are cropped to the given extend. Keys should be
        the coordinates labels of dem_da and values should be slices of [min,max]
        extend. Default is None.

    See also
    --------
    gradient, sobel, _gradient_names
    """

    hlp.check_dem(dem_da)
    logger.info(f"***Starting gradients computation for scales {scales} meters***")
    if not hasattr(scales, "__iter__"):
        scales = [scales]
    if not hasattr(sig_ratios, "__iter__"):
        sig_ratios = [sig_ratios] * len(scales)

    scales_pxl, res_meters = hlp.scale_to_pixel(scales, dem_da)
    sigmas = scales_pxl / CFG.scale_std

    for idx, sigma in enumerate(sigmas):
        logger.info(
            f"Computing scale {scales[idx]} meters with sigma ratio "
            f"{sig_ratios[idx]} ..."
        )
        names = _gradient_names(scales[idx], sig_ratios[idx])
        arrays = gradient(
            dem=dem_da.values,
            sigma=sigma,
            res_meters=res_meters,
            sig_ratio=sig_ratios[idx],
        )

        for array, name in zip(arrays, names):
            array[ind_nans] = np.nan
            hlp.to_netcdf(array, dem_da.coords, name, crop)

        del arrays


@hlp.timer
def gradient(dem, sigma, res_meters, sig_ratio=1):
    """Compute directional derivatives, slope and aspect over a digital elevation
    model.

    Parameters
    ----------
    dem : array representing the DEM.
    sigma : scalar
        Standard deviation for the gaussian filters. This is set at 1/4 of the
        scales to which topo descriptors are computed (i.e. scale = 4*sigma)
    res_meters : dict with two 1-D or 2-D arrays
        Resolution in meters of each DEM grid points in the x and y directions.
        This is the second element returned by the function hlp.scale_to_pixel.
    sig_ratio (optional) : scalar
        Ratio w.r.t sigma to define the standard deviation of the gaussian in
        the direction perpendicular to the derivative. Default is 1 (i.e. same
        smoothing on both directions).

    Returns
    -------
    list of four arrays :
        First element is the W-E derivative, second the S-N derivative, third
        the slope (i.e. magnitude of the gradient) and fourth the aspect (i.e.
        direction of the gradient).

    See also
    --------
    scipy.ndimage.gaussian_filter, np.gradient, topo_helpers.scale_to_pixel
    """

    if sigma <= 1:
        dx, dy = sobel(dem)  # for lowest scale, use sobel filter instead
    elif sig_ratio == 1:  # faster solution when sig_ratio is 1
        dy, dx = np.gradient(ndimage.gaussian_filter(dem, sigma))
    else:
        sigma_perp = sigma * sig_ratio
        dx = np.gradient(ndimage.gaussian_filter(dem, (sigma_perp, sigma)), axis=1)
        dy = np.gradient(ndimage.gaussian_filter(dem, (sigma, sigma_perp)), axis=0)

    _normalize_dxy(dx, dy, res_meters)

    slope = np.sqrt(dx ** 2, dy ** 2)
    aspect = (
        180 + np.degrees(np.arctan2(dx, dy))
    ) % 360  # north faces = 0°, east faces = 90°

    return [dx, dy, slope, aspect]


def _gradient_names(scale, sig_ratio):
    """Return names for the arrays in output of the gradient function"""

    name_dx = f"WE_DERIVATIVE_{scale}M_SIGRATIO{sig_ratio:.3g}"
    name_dy = f"SN_DERIVATIVE_{scale}M_SIGRATIO{sig_ratio:.3g}"
    name_slope = f"SLOPE_{scale}M_SIGRATIO{sig_ratio:.3g}"
    name_aspect = f"ASPECT_{scale}M_SIGRATIO{sig_ratio:.3g}"

    return [name_dx, name_dy, name_slope, name_aspect]


def sobel(dem):
    """Compute directional derivatives, based on the sobel filter. The sobel is
    a combination of a derivative and a smoothing filter. It is defined over a
    3x3 kernel.

    Parameters
    ----------
    dem : array representing the DEM.

    Returns
    -------
    dx : 2-D array
        Derivative in the x direction (i.e. along second axis)
    dy : 2-D array
        Derivative in the y direction (i.e. along first axis)

    See also
    --------
    ndimage.convolve
    """

    sobel_kernel = np.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]], dtype=np.float32)

    sobel_kernel /= np.sum(np.abs(sobel_kernel))
    dx = ndimage.convolve(dem, sobel_kernel)
    dy = ndimage.convolve(dem, sobel_kernel.T)

    return dx, dy


def _normalize_dxy(dx, dy, res_meters):
    """Normalize directional derivatives, based on (projected) grid resolution.
    This is useful when the DEM does not come on an equidistant projected grid.
    Normalization occurs 'in place'.

    Parameters
    ----------
    dx : 2-D array
        Derivative in the x direction (i.e. along second axis)
    dy : 2-D array
        Derivative in the y direction (i.e. along first axis)
    res_meters : dict with two 1-D or 2-D arrays
        Resolution in meters of each DEM grid points in the x and y directions.
        This is the second element returned by the function hlp.scale_to_pixel.

    See also
    --------
    topo_helpers.scale_to_pixel
    """

    mean_res = np.mean(np.abs([res_meters["x"].mean(), res_meters["y"].mean()]))
    x_res = res_meters["x"] / mean_res
    y_res = res_meters["y"] / mean_res
    if len(y_res.shape) == 1:
        y_res = y_res[:, np.newaxis]

    dx /= x_res
    dy /= y_res


# TODO: Sx and Relief
