import numpy as np
import xarray as xr
    
def get_dem_netcdf(filepath):

    dem_ds = xr.open_dataset(filepath, decode_times=False)
    return dem_ds.to_array().isel(variable=0, drop =True)

def round_up_to_odd(f):
    return np.array(np.ceil(f) // 2 * 2 + 1, dtype = np.int64)