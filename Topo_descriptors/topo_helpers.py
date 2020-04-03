import numpy as np
import xarray as xr
from env_setting import DATAPATH, RESOURCES

def get_dem_netcdf(filename):
    """Load the DEM into a xarray DataArray and filter NaNs"""

    dem_ds = xr.open_dataset(RESOURCES/filename, decode_times=False)
    dem_da = dem_ds.to_array().isel(variable=0, drop =True)
    
    return dem_da.where(dem_da != -9999.)

def save_topo_netcdf(da):
    """Save a DataArray of topographic descriptors in NetCDF"""
    
    filename = 'topo_descriptor_' + da.name + '.nc'
    da.to_dataset().to_netcdf(DATAPATH / filename)

def _to_da(array, dem_da, name):
    """Convert a numpy array representing a topographic descriptor into a xarray 
    DataArray with the same coordinates as the input DEM DataArray and a specified
    name
    """
    return xr.DataArray(array, coords=dem_da.coords, name=str.upper(name))
 
def dist_to_pixel(dist_list, dem_da):
    """Convert distances in KM to the closest odd number of pixels based on the DEM 
    resolution.
    
    Parameters
    ----------
    dist_list : list
        Scales in KM on which we want to compute the topographic descriptor. 
        Corresponds to the size of the squared kernel used to compute it.
    dem_da : xarray DataArray representing the DEM and its grid coordinates
    in meters (CH coordinates).
        
    Returns
    -------
    list with the corresponding scales in pixel size.
    """
    
    res_0 = dem_da['chx'].diff('chx').mean()
    res_1 = dem_da['chy'].diff('chy').mean()
    res = np.mean(np.abs([res_0, res_1]))
    return round_up_to_odd(dist_list * 1000 / res)

def round_up_to_odd(f):
    """round float to the nearest odd integer"""
    return np.array(np.ceil(f) // 2 * 2 + 1, dtype = np.int64)