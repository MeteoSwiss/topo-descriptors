import numpy as np
import xarray as xr
import pyproj
import functools
import time
import datetime as dt
from env_setting import DATAPATH
import logging
logger = logging.getLogger(__name__)


def get_dem_netcdf(path_dem):
    """Load the DEM into a xarray DataArray and filter NaNs
    
    Parameters
    ----------
    path_dem: string
        absolute or relative path to the DEM netCDF file.
        
    Returns
    -------
    xarray DataArray with the DEM values.
    """

    dem_ds = xr.open_dataset(path_dem, decode_times=False)
    dem_da = dem_ds.to_array().isel(
                                variable=0, 
                                drop=True).reset_coords(
                                                drop=True).astype(np.float32)
    
    return dem_da.where(dem_da > -100)


def to_netcdf(array, coords, name, crop=None):
    """Save an array of topographic descriptors in NetCDF. It is first converted
    into a xarray DataArray with the same coordinates as the input DEM DataArray 
    and a specified name.
    
    Parameters
    ----------
    array : array to be saved as netCDF
    coords : dict
        Coordinates for the array (i.e. those of the DEM).
    name : string
        Name for the array
    crop (optional) : dict
        The array is cropped to the given extend before being saved. Keys should 
        be coordinates labels as in coords and values should be slices of [min,max] 
        extend. Default is None.
    """
    
    name = str.upper(name)
    da = xr.DataArray(array, coords=coords, name=name).sel(crop)
    filename = f'topo_{name}.nc'
    da.to_dataset().to_netcdf(DATAPATH / filename)


def scale_to_pixel(scales, dem_da):
    """Convert distances in meters to the closest odd number of pixels based on 
    the DEM resolution.
    
    Parameters
    ----------
    scales : list of scalars
        Scales in meters on which we want to compute the topographic descriptor. 
        Corresponds to the size of the squared kernel used to compute it.
    dem_da : xarray DataArray representing the DEM and its grid coordinates. 
    Coordinates must be in LV03/LV95 and named 'chx', 'chy' or in WGS84 and named
    'lon', 'lat'. In the latter case, they are reprojected in LV03 to derive the
    average resolution in meters.
        
    Returns
    -------
    list of int :
        Contain the corresponding scales in pixel size.
    dict with two 1-D or 2-D arrays :
        Resolution in meters of each DEM grid points in the x and y directions.
    """
    
    dim_x, dim_y = dem_da.dims[1], dem_da.dims[0]
    x_coords, y_coords = dem_da[dim_x].values, dem_da[dim_y].values
    if dim_y == 'chy':
        pass
    elif dim_y == 'lat':
        logger.warning(
        f'Reprojecting coordinates from WGS84 to LV03 to have units of meters')
        x_coords, y_coords = np.meshgrid(x_coords, y_coords)
        x_coords, y_coords = pyproj.transform(pyproj.Proj('+init=EPSG:4326'), 
                                              pyproj.Proj('+init=EPSG:21781'), 
                                              x_coords, y_coords)
        x_coords, y_coords = x_coords.astype(np.float32), y_coords.astype(np.float32)
    else:
        raise ValueError(f'Unknown coordinates reference system. Cannot convert' 
                         f' to meters to find the average pixel resolution')
      
    n_dims = len(x_coords.shape)
    x_res = np.gradient(x_coords, axis=n_dims - 1)
    y_res = np.gradient(y_coords, axis=0)
    mean_res = np.mean(np.abs([x_res.mean(), y_res.mean()]))
    
    return round_up_to_odd(np.array(scales) / mean_res), {'x' : x_res, 'y' : y_res}


def round_up_to_odd(f):
    """round float to the nearest odd integer"""
    
    return np.asarray(np.ceil(f) // 2 * 2 + 1, dtype = np.int64)


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
    
    sigmas = [fact if fact else np.nan 
                      for fact in smth_factors] * scales_pxl / 4 # scale = 4std

    return [None if np.isnan(sigma) else sigma for sigma in sigmas]


def fill_na(dem_da):
    """get indices of NaNs and interpolates them.
    
    Parameters
    ----------
    dem_da : xarray DataArray containing the elevation data.
    
    Returns
    -------
    ind_nans : tuple of two 1D arrays
        Contains the row / column indices of the NaNs in the original dem.
    Xarray DataArray with interpolated NaNs in x direction using "nearest" method.
    """
    
    ind_nans = np.where(np.isnan(dem_da))
    return ind_nans, dem_da.interpolate_na(dim=dem_da.dims[1], 
                               method='nearest', 
                               fill_value='extrapolate')


def timer(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        t_start = time.monotonic()
        value = func(*args, **kwargs)
        t_elapsed = str(
                    dt.timedelta(seconds=time.monotonic() - t_start)
                    ).split('.', 2)[0]
        logger.info(f'Computed in {t_elapsed} (HH:mm:ss)') 
        return value
    return wrapper_timer
