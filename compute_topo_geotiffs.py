import numpy as np
from env_setting import RESOURCES
import Topo_descriptors.tpi_valley as tpivalley
import Topo_descriptors.sx_relief as sxrelief

#%%
# Script to compute spatial descriptors from DEM

#%% Define parameters

# Wich DEM to use
dem_filename = 'dhm25_clean.nc'
path_dem = RESOURCES /dem_filename

# Sx parameters
theta = 5
alpha = 5
height = 10
d_min_sx = np.array([0,100,100,100]) # dist in meters
d_max_sx = np.array([500,1000,2000,4000]) # dist in meters

# Tpi parameters
dist_tpi = [500,2000,4000]  # diameter in meter

# Valley parameters
dist_valley = [2000,6000,10000,20000] # convolution size in meters
flat_valley = [0.15,0.3] # fraction of zeros in the middle -> version with one 0 automatically included.

#%% Launch computation
       
#if __name__ == '__main__':
#    sxrelief.get_sxrelief_geotiff(theta, alpha, d_min_sx, d_max_sx, height, path_dem)

tpivalley.tpi_netcdf(path_dem, dist_tpi)
tpivalley.valley_netcdf(path_dem, dist_valley, flat_valley)

