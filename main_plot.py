'''Script to plot topographic descriptors maps'''

from env_setting import RESOURCES
from env_setting import DATAPATH
import Topo_descriptors.topo_helpers as hlp
from matplotlib import pyplot as plt
import xarray as xr
import numpy as np

#%% get the DEM

path_dem = RESOURCES / 'DEM_CCS4_50m_buffer_v-0-3.nc'
dem_da = hlp.get_dem_netcdf(path_dem)
inca_domain = {'chx':slice(255000,965000), 'chy':slice(-160000,480000)}

dem_da= dem_da.sel(inca_domain)

#%% make plots

test1 = xr.open_dataset(str(DATAPATH / 'topo_TPI_60000M.nc'))
test2 = xr.open_dataset(str(DATAPATH / 'topo_TPI_60000M_SMTHFACT0.5.nc'))

area = {'chx':slice(615000,675000), 'chy':slice(130000,170000)}

plt.imshow(np.flipud(dem_da))
plt.figure()
plt.imshow(np.flipud(test1['TPI_60000M']))
plt.figure()
plt.imshow(np.flipud(test2['TPI_60000M_SMTHFACT0.5']))

