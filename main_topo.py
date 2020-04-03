import numpy as np
from env_setting import RESOURCES
import Topo_descriptors.topo as tp
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
dist_tpi = [0.1, 0.5,2,4]  # convolution size (diameter) in km

# Valley parameters
dist_valley = [2,6,10,20] # convolution size in km
flat_valley = [0.15,0.3] # fraction of zeros in the middle -> version with one 0 automatically included.

#%% Launch computation
       
#if __name__ == '__main__':
#    sxrelief.get_sxrelief_geotiff(theta, alpha, d_min_sx, d_max_sx, height, path_dem)

tp.compute_tpi(path_dem, dist_tpi)
tp.compute_valley(path_dem, dist_valley)

