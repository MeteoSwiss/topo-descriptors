'''
This scripts creates a single netCDF file out of a list of SRTM zipped tiles.
'''

import numpy as np
import xarray as xr

from Topo_descriptors.preprocessing import srtm_mosaic

# path to SRTM archive (containing zipped SRTM90 v4.1 tiles)
# >>> ls -l /repos/repos_data/climate/cmsaf/topo/SRTM_CGIAR/GeoTiff
# -rw-rw-r-- 1 res gsx-apk   713075 Jun 18  2009 srtm_01_02.zip
# -rw-rw-r-- 1 res gsx-apk   130923 Jun 18  2009 srtm_01_07.zip
# ...
# -rw-rw-r-- 1 res gsx-apk  3348530 Jun 19  2009 srtm_72_21.zip
# -rw-rw-r-- 1 res gsx-apk   136698 Jun 19  2009 srtm_72_22.zip
dir = '/repos/repos_data/climate/cmsaf/topo/SRTM_CGIAR/GeoTiff/'

# compute the SRTM mosaic from a list of relevant tiles
tiles = ['srtm_37_03','srtm_38_03','srtm_39_03', 'srtm_37_04','srtm_38_04','srtm_39_04']
raster, extent, res = srtm_mosaic(dir, tiles)

# extract all coordinates of the pixel centers from the grid extent, shape and
# resolution
shape = raster.shape
coord_lat = np.linspace(extent[2] + 0.5 * res[0], extent[3] - 0.5 * res[0], shape[0])
coord_lon = np.linspace(extent[0] + 0.5 * res[1], extent[1] - 0.5 * res[1], shape[1])

# round up to a given precision
precision = 1e-7
coord_lon = coord_lon // precision * precision
coord_lat = coord_lat // precision * precision

# change value assigned to missing data
raster[raster == -32768] = -9999

# as xarray
srtm = xr.Dataset(
        {'SRTM90': (['lat', 'lon'], raster[::-1, :].astype(np.int16))},
        coords={
            'lat': ('lat', coord_lat.astype(np.float32)),
            'lon': ('lon', coord_lon.astype(np.float32)),
            }
        # TODO: add more attributes
    )
print(srtm)

# save into working directory
print('saving to netcdf...', end=' ', flush=True)
srtm.to_netcdf('srtm_mosaic_WSG84.nc')
print('done!')
