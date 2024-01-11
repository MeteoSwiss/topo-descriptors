import datetime as dt
import functools
import logging
import time
from pathlib import Path

import numpy as np
import xarray as xr
import utm

from topo_descriptors import CFG


logger = logging.getLogger(__name__)


def get_dem_netcdf(path_dem):
    """Load the DEM into a xarray Dataset and filter NaNs

    Parameters
    ----------
    path_dem: string
        absolute or relative path to the DEM netCDF file.

    Returns
    -------
    xarray Dataset with the DEM values.
    """

    dem_ds = xr.open_dataset(path_dem).astype(np.float32).squeeze(drop=True)
    return dem_ds.where(dem_ds > CFG.min_elevation)


def to_netcdf(array, dem_ds, name, crop=None, outdir=".", units=None):
    """Save an array of topographic descriptors in NetCDF. It has the same coordinates
    and attributes as the input DEM and a specified name.

    Parameters
    ----------
    array : array to be saved as netCDF
    dem_ds : xarray Dataset
        Original DEM dataset from which coordinates and attributes are copied.
    name : string
        Name for the array
    crop (optional) : dict
        The array is cropped to the given extend before being saved. Keys should
        be coordinates labels as in coords and values should be slices of [min,max]
        extend. Default is None.
    outdir (optional) : string
        The path to the output directory. Save to working directory by default.
    units (optional) : string
        Units of the topographic descriptor array added as variable attribute.
    """

    name = str.upper(name)
    outdir = Path(outdir)
    ds = xr.Dataset(
        {name: (get_da(dem_ds).dims, array)}, coords=dem_ds.coords, attrs=dem_ds.attrs
    ).sel(crop)
    if units is not None:
        ds[name].attrs.update(units=units)

    filename = f"topo_{name}.nc"
    ds.to_netcdf(outdir / filename)
    logger.info(f"saved: {outdir / filename}")


def scale_to_pixel(scales, dem_ds):
    """Convert distances in meters to the closest odd number of pixels based on
    the DEM resolution.

    Parameters
    ----------
    scales : list of scalars
        Scales in meters on which we want to compute the topographic descriptor.
        Corresponds to the size of the squared kernel used to compute it.
    dem_ds : xarray Dataset representing the DEM and its grid coordinates.
    Coordinates must be projected or in WGS84 and named 'x', 'y'. In the latter case,
    they are reprojected to UTM to derive the average resolution in meters.

    Returns
    -------
    list of int :
        Contain the corresponding scales in pixel size.
    dict with two 1-D or 2-D arrays :
        Resolution in meters of each DEM grid points in the x and y directions.
    """
    check_dem(dem_ds)
    x_coords, y_coords = dem_ds["x"].values, dem_ds["y"].values
    epsg_code = dem_ds.attrs["crs"].lower()
    if "epsg:4326" in epsg_code:
        logger.debug(
            "Reprojecting coordinates from WGS84 to UTM to obtain units of meters"
        )
        x_coords, y_coords = np.meshgrid(x_coords, y_coords)
        x_coords, y_coords, _, _ = utm.from_latlon(y_coords, x_coords)
        x_coords, y_coords = x_coords.astype(np.float32), y_coords.astype(np.float32)

    n_dims = len(x_coords.shape)
    x_res = np.gradient(x_coords, axis=n_dims - 1)
    y_res = np.gradient(y_coords, axis=0)
    mean_res = np.mean(np.abs([x_res.mean(), y_res.mean()]))
    logger.debug(f"Estimated resolution: {mean_res:.0f} meters.")

    return round_up_to_odd(np.array(scales) / mean_res), {"x": x_res, "y": y_res}


def round_up_to_odd(f):
    """round float to the nearest odd integer"""

    return np.asarray(np.round((f - 1) / 2) * 2 + 1, dtype=np.int64)


def get_sigmas(smth_factors, scales_pxl):
    """Return scales expressed in standard deviations for gaussian filters.

    Parameters
    ----------
    smth_factors : list of scalars or None elements or a combination of both.
        Factors by which the scales in pixel must be multiplied. None or zeros
        results in None in the output.
    scales_pxl : list of int
        Scales expressed in pixels.

    Returns
    -------
    list of scalars/None elements representing scales in standard deviations.
    """

    sigmas = (
        [fact if fact else np.nan for fact in smth_factors] * scales_pxl / CFG.scale_std
    )

    return [None if np.isnan(sigma) else sigma for sigma in sigmas]


def fill_na(dem_ds):
    """get indices of NaNs and interpolates them.

    Parameters
    ----------
    dem_ds : xarray Dataset containing the elevation data.

    Returns
    -------
    ind_nans : tuple of two 1D arrays
        Contains the row / column indices of the NaNs in the original dem.
    Xarray Dataset with interpolated NaNs in x direction using "nearest" method.
    """

    ind_nans = np.where(np.isnan(get_da(dem_ds)))
    return ind_nans, dem_ds.interpolate_na(
        dim="x", method="nearest", fill_value="extrapolate"
    )


def timer(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        t_start = time.monotonic()
        value = func(*args, **kwargs)
        t_elapsed = str(dt.timedelta(seconds=time.monotonic() - t_start)).split(".", 2)[
            0
        ]
        logger.info(f"Computed in {t_elapsed} (HH:mm:ss)")
        return value

    return wrapper_timer


def check_dem(dem):
    """
    Check that the input dem conforms to the data model, namely:
      - instance of xarray.Dataset
      - 2D field
      - y and x dimensions
      - crs attribute specifying an EPSG code.
    """
    if not isinstance(dem, xr.Dataset):
        raise ValueError("dem must be a xr.Dataset")
    if dem[list(dem)[0]].dims != ("y", "x"):
        raise ValueError("dem dimensions must be ('y', 'x')")
    if not "crs" in dem.attrs:
        raise KeyError("missing 'crs' (case sensitive) attribute in dem")
    if not "epsg:" in dem.attrs["crs"].lower():
        raise ValueError(
            "missing 'epsg:' (case insensitive) key in the 'crs' attribute"
        )


def get_da(dem_ds):
    """
    Return the xarray.DataArray with DEM values without knowing its name.
    """

    return dem_ds[list(dem_ds)[0]]
