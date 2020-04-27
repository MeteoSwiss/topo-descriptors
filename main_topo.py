#%%
# Script to compute spatial descriptors from DEM

from env_setting import RESOURCES
import Topo_descriptors.topo as tp
import Topo_descriptors.topo_helpers as hlp
import logging
logging.basicConfig(level=logging.INFO)
logging.captureWarnings(True)
logger = logging.getLogger(__name__)

#%% get the DEM

path_dem = RESOURCES / 'DEM_CCS4_50m_buffer_v-0-3.nc'
dem_da = hlp.get_dem_netcdf(path_dem)
inca_domain = {'chx':slice(255000,965000), 'chy':slice(-160000,480000)}

#dem_da= dem_da.sel(chx=slice(300000,800000), chy=slice(00000,300000))
ind_nans, dem_da = hlp.fill_na(dem_da)

#%% Define parameters

topo_scales = [100, 300, 500, 1000, 2000, 4000, 6000, 10000, 20000, 30000, 
               60000, 100000] # convolution scale in meters

#%% Launch computation
       
# TPI
tp.compute_tpi(dem_da, 
               topo_scales, 
               smth_factors=None, 
               ind_nans=ind_nans, 
               crop=inca_domain)

# Gradients
tp.compute_gradient(dem_da, 
                    topo_scales, 
                    sig_ratios=1, 
                    ind_nans=ind_nans, 
                    crop=inca_domain)

# Valley Index
tp.compute_valley_ridge(dem_da, 
                        topo_scales,
                        mode='valley',
                        flat_list=[0, 0.2, 0.4], 
                        smth_factors=None, 
                        ind_nans=ind_nans,
                        crop=inca_domain)

# Ridge Index
tp.compute_valley_ridge(dem_da, 
                        topo_scales,
                        mode='ridge',
                        flat_list=[0, 0.1, 0.2], 
                        smth_factors=None, 
                        ind_nans=ind_nans,
                        crop=inca_domain)

